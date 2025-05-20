clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end
            
vrs_path   = "vrs_data/";

vrs_path = vrs_path + "Resolution_HERCULES-Diverging-TxColumn";
vrs_name = "240905_ATS539_Resolution_HERCULES-Diverging-TxColumn_Intensity_49";

% vrs_path = vrs_path + "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxColumn";
% vrs_name = "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxColumn_Intensity_00";

% vrs_path = vrs_path + "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxRow";
% vrs_name = "250204_MN32-1_M18_Side_herc_diverging_HERCULES-Diverging-TxRow_Intensity_00";

vrs_name   = vrs_name + ".vrs";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

vrs  = VRSFile(fullfile(vrs_path, vrs_name));
data = vrs.GetData();

load(fullfile(vrs_path, "postEnv.mat"), "Trans", "Receive", "Resource", "TW");
load(fullfile(vrs_path, "preEnv.mat"), "sparseElements", "scan");

receive             = Receive([Receive.bufnum] == 1);
receive_orientation = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation;
recieve_elements    = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation.GetElementCount(scan.Die);



bp.rf_raw_dim     = struct('x', size(data, 1), 'y', size(data, 2));
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

bp.channel_mapping = full_channel_mapping(channel_offset:channel_offset+recieve_elements-1);

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

% bp.das_shader_id = 0; % Forces
bp.das_shader_id = 1; % Hercules

cs_stages = uint8([
	OGLShaderStage.HADAMARD, ...
	OGLShaderStage.DAS, ...
	OGLShaderStage.MIN_MAX, ...
]);

group_cells = cell(1,16);

% The end of every channel has padded zeros not from any transmit
expected_length = bp.dec_data_dim.x * 128;

shuffled_data = data(1:expected_length, :);

shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(data, 2)];
shuffled_data = reshape(shuffled_data, shuffled_size);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-10, 10]/1000;
y_range = [-30, 30]/1000;
z_range = [5, 100]/1000;

bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.0002;
axial_resolution = 0.0002;

bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution );
bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution );
bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution );
bp.output_points.w = 1; % Number of frames for averaging


%% Data prep
readi_group_count = 1;

bp.readi_group_size = tx_count/readi_group_count;
bp.readi_group_id = 0;

% Break up transmits into Readi groups
for i=1:readi_group_count
    end_tx = i * bp.readi_group_size;
    start_tx = end_tx - (bp.readi_group_size - 1);

    slice = shuffled_data(:,start_tx:end_tx,:);
    slice = reshape(slice, [], bp.rf_raw_dim.y);

    group_cells{i} = slice;
end

bp.rf_raw_dim.x = length(slice);


if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end

loadlibrary('cuda_transfer')
calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

fprintf("Sending data\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Readi loop

output_counts_xyz.x = bp.output_points.x;
output_counts_xyz.y = bp.output_points.y;
output_counts_xyz.z = bp.output_points.z;

low_res_cells = cell(1,readi_group_count);

% Complex volumes aren't supported so they're interleaved
interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];


for g = 1:readi_group_count

    calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

    fprintf("Readi group %i. \n", g);

    volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

    [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
        pipe_name, smem_name, group_cells{g}, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
    fprintf("Received Response\n")
    
    volume = reshape(volume_ptr,interleaved_volume_size );

    real_vol = volume(1:2:end,:,:);
    im_vol = volume(2:2:end,:,:);
    
    low_res_cells{g} = complex(real_vol, im_vol);

    pause(0.5);

    bp.readi_group_id = bp.readi_group_id + 1;
   
end

unloadlibrary('cuda_transfer')


%% Post processing 
volume_size = interleaved_volume_size;
volume_size(1) = interleaved_volume_size(1)/2;
high_res_volume = zeros(volume_size);

dynamic_range = 40;

threshold = 210;

processed_low_res = low_res_cells;
low_res_images = low_res_cells;
for g = 1:readi_group_count
    high_res_volume = high_res_volume + low_res_cells{g};
    processed_low_res{g} = process_volume(low_res_cells{g}, dynamic_range, threshold); 

    low_res_images{g} = squeeze(processed_low_res{g}(:,round(output_counts_xyz.y/2),:));
end
processed_volume = process_volume(high_res_volume, dynamic_range, threshold);

x_image = squeeze(processed_volume(round(output_counts_xyz.x/2),:,:)).';
y_image = squeeze(processed_volume(:, round(output_counts_xyz.y/2),:)).';

processed_volume = flip(processed_volume,3);

figure();
imagesc(x_image);
colormap("gray");
axis('image');
clim([-dynamic_range, 0]); 
colorbar;
title("Hercules X Axis")

figure();
imagesc(y_image);
colormap("gray");
axis('image');
clim([-dynamic_range, 0]); 
colorbar;
title("Hercules Y Axis")

%% 
 volumeViewer(processed_volume)
%%
% figure(2);
% 
% % Main title
% % sgtitle(speed_str + ' Readi Low Resolution Images','FontSize',16);
% sgtitle('Stationary Readi Low Resolution Images','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(2, 2, i); 
%     colormap("gray");
%     imagesc((low_res_images{i}).');
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end



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