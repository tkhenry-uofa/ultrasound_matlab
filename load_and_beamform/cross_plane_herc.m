clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end
            
vrs_path   = "vrs_data/melanoma_mouse/";

col_vrs_path = vrs_path + "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxColumn";
col_vrs_name = "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxColumn_Intensity_00";

row_vrs_path = vrs_path + "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxRow";
row_vrs_name = "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxRow_Intensity_00";

row_vrs_name = row_vrs_name + ".vrs";
col_vrs_name = col_vrs_name + ".vrs";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

vrs_col  = VRSFile(fullfile(col_vrs_path, col_vrs_name));
col_data = vrs_col.GetData();

vrs_row  = VRSFile(fullfile(row_vrs_path, row_vrs_name));
row_data = vrs_row.GetData();


load(fullfile(col_vrs_path, "postEnv.mat"), "Trans", "Receive", "Resource", "TW");
load(fullfile(col_vrs_path, "preEnv.mat"), "sparseElements", "scan");

receive             = Receive([Receive.bufnum] == 1);
receive_orientation = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation;
recieve_elements    = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation.GetElementCount(scan.Die);



bp.rf_raw_dim     = struct('x', size(col_data, 1), 'y', size(col_data, 2));
bp.dec_data_dim.x = max(1 + [receive.endSample] - [receive.startSample], [], "all");
bp.dec_data_dim.y = recieve_elements;
bp.dec_data_dim.z = max([receive.acqNum]);
bp.dec_data_dim.w = 0; % Averaging

tx_count = bp.dec_data_dim.z;

bp.sampling_frequency = receive(1).samplesPerWave * Trans.frequency(1) * 1e6;
bp.center_frequency   = Trans.frequency * 1e6 *2;
bp.speed_of_sound     = Resource.Parameters.speedOfSound;

bp.time_offset = TW(1).Parameters(3) / bp.center_frequency;


if (exist('sparseElements'))
	bp.uforces_channels = sparseElements - 1;
else
	bp.uforces_channels = 0:(recieve_elements - 1);
end

full_channel_mapping  = uint16(Trans.ConnectorES - 1);
channel_offset = 1 + recieve_elements * uint16(receive_orientation.Contains(tobe.Orientation.Column));
col_channel_mapping = full_channel_mapping(1:recieve_elements);
row_channel_mapping = full_channel_mapping(recieve_elements+1:end);



%%
dieSize       = scan.Die.GetSize();

bp.xdc_transform = single([
                    1,0,0,dieSize(1)/2;
                    0,1,0,dieSize(2)/2;
                    0,0,1,0;
                    0,0,0,1;
                    ]);
bp.xdc_element_pitch = scan.Die.Pitch;


bp.focal_depths  = ones(1,tx_count) .* scan.TransmitEvents(1, 1).FocalDepth;
bp.transmit_angles  = single(zeros(1,tx_count));

% NOTE: plane and position along plane normal for beamforming 2D HERCULES
bp.beamform_plane = 0;
bp.off_axis_pos = 0;
fc = bp.center_frequency;
fs = bp.sampling_frequency;

bp.f_number = 1;

% The end of every channel has padded zeros not from any transmit
expected_length = bp.dec_data_dim.x * 128;

shuffled_data = col_data(1:expected_length, :);

shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(col_data, 2)];
shuffled_data = reshape(shuffled_data, shuffled_size);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-30, 30]/1000;
y_range = [-30, 30]/1000;
z_range = [10, 40]/1000;

bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.00015;
axial_resolution = 0.00015;

bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution );
bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution );
bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution );
bp.output_points.w = 1; % Number of frames for averaging


bp.readi_group_size = tx_count;
bp.readi_group_id = 0;


if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end
loadlibrary('cuda_transfer')

fprintf("Sending data\n")

%% Column beamform

bp.channel_mapping = col_channel_mapping;

output_counts_xyz.x = bp.output_points.x;
output_counts_xyz.y = bp.output_points.y;
output_counts_xyz.z = bp.output_points.z;


% Complex volumes aren't supported so they're interleaved
interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];

calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

fprintf("Beamforming volume 1 \n");

col_volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

[~,~,~,col_volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
    pipe_name, smem_name, col_data, bp.rf_raw_dim, output_counts_xyz, col_volume_ptr);
fprintf("Received Response\n")

volume = reshape(col_volume_ptr,interleaved_volume_size );

real_vol = volume(1:2:end,:,:);
im_vol = volume(2:2:end,:,:);

col_volume = complex(real_vol, im_vol);

%% Row Beamform
if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end
bp.channel_mapping = row_channel_mapping;
calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);
fprintf("Beamforming volume 2 \n");

row_volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

[~,~,~,row_volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
    pipe_name, smem_name, row_data, bp.rf_raw_dim, output_counts_xyz, row_volume_ptr);
fprintf("Received Response\n")

volume = reshape(row_volume_ptr,interleaved_volume_size );

real_vol = volume(1:2:end,:,:);
im_vol = volume(2:2:end,:,:);

row_volume = complex(real_vol, im_vol);

    
unloadlibrary('cuda_transfer')


%% Post processing 

dynamic_range = 55;
threshold = 160;

total_volume = col_volume + permute(row_volume, [2, 1, 3]);

processed_total = process_volume(total_volume, dynamic_range, threshold);
processed_col = process_volume(col_volume, dynamic_range, threshold);
processed_row = process_volume(row_volume, dynamic_range, threshold);


processed_col = flip(processed_col,3);
processed_row = flip(processed_row,3);
processed_total = flip(processed_total,3);



%% 
 volumeViewer(processed_total)




%%

function processed_volume = process_volume(volume,dynamic_range, threshold)

    mag_volume = abs(volume); 

    threshold_value = 10^(threshold/20);

    % threshold_min = 10^(-120/20);

     mag_volume = min(mag_volume, threshold_value);
    
    % mag_volume = mag_volume/max(mag_volume,[],"all");
    
    processed_volume = 20*log10(mag_volume);

    processed_volume = processed_volume - max(processed_volume,[],"all");

    processed_volume = max(processed_volume, -dynamic_range);
end