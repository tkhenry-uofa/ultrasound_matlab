clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path   = "../real_data/vrs_data/match_filter_test/";

data_path = data_path + "250514_MN32-5_reso_FORCES-TxColumn/";
data_file_name = "250514_MN32-5_reso_FORCES-TxColumn_00";

params_path = data_path + '250514_MN32-5_reso_FORCES-TxColumn.bp';

data_file_name   = data_file_name + ".zst";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

data_file = fopen(data_path + data_file_name, "r");
raw_data = fread(data_file, '*uint8');
fclose(data_file); 

data = ornot_zstd_decompress_mex(raw_data);

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));
data = reshape(data, raw_bp.raw_data_dim(1),raw_bp.raw_data_dim(2));


% data(1:2,:) = 0;

%%
bp.decode          = raw_bp.decode_mode;
bp.beamform_plane  = raw_bp.beamform_mode;

% bp.rf_raw_dim = raw_bp.rf_raw_dim;
bp.rf_raw_dim.x      = raw_bp.raw_data_dim(1);
bp.rf_raw_dim.y      = raw_bp.raw_data_dim(2);
% bp.dec_data_dim    = raw_bp.dec_data_dim;
bp.dec_data_dim.x    = raw_bp.decoded_data_dim(1);
bp.dec_data_dim.y    = raw_bp.decoded_data_dim(2);
bp.dec_data_dim.z    = raw_bp.decoded_data_dim(3);
bp.dec_data_dim.w    = raw_bp.decoded_data_dim(4);

% Map transducer properties
bp.xdc_element_pitch = raw_bp.transducer_element_pitch;
bp.xdc_transform     = raw_bp.transducer_transform_matrix;

% Map channel and angle related fields
bp.channel_mapping  = raw_bp.channel_mapping;
bp.transmit_angles  = raw_bp.steering_angles;
bp.focal_depths     = raw_bp.focal_depths;

% Map acoustic parameters
bp.speed_of_sound    = raw_bp.speed_of_sound;
bp.center_frequency  = raw_bp.center_frequency;
bp.sampling_frequency= raw_bp.sampling_frequency;
bp.time_offset       = raw_bp.time_offset;

% Map transmit mode
bp.transmit_mode     = raw_bp.transmit_mode;


%%
% 0 = Forces, 1 = Uforces, 2 = Hercules
bp.das_shader_id = 0;
% bp.das_shader_id = 2;

% NOTE: plane and position along plane normal for beamforming 2D HERCULES
bp.beamform_plane = 0;
bp.off_axis_pos = 0;
fc = bp.center_frequency;
fs = bp.sampling_frequency;

bp.f_number = 0.5;

group_cells = cell(1,16);

% The end of every channel has padded zeros not from any transmit
expected_length = bp.dec_data_dim.x * bp.dec_data_dim.z;
shuffled_data = data(1:expected_length, :);

shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(data, 2)];
shuffled_data = reshape(shuffled_data, shuffled_size);

start_depth = 5/1000;

start_sample = round(2 * start_depth * bp.sampling_frequency / bp.speed_of_sound);
tx_end_sample = start_sample -1;

shuffled_data(1:tx_end_sample,:,:) = 0;

% shuffled_data = shuffled_data(start_sample:end,:,:);
% 
% bp.dec_data_dim.x = bp.dec_data_dim.x - tx_end_sample;
% bp.rf_raw_dim.x = bp.rf_raw_dim.x - tx_end_sample * 256;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-20, 20]/1000;
y_range = [-20, 20]/1000;
z_range = [5, 100]/1000;

if bp.das_shader_id == 0
    y_range = [0,0]; % Forces
end

bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.0002;
axial_resolution = 0.0002;

bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution );
bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution );
bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution );
bp.output_points.w = 1; % Number of frames for averaging

if bp.output_points.y <= 0
    bp.output_points.y = 1;
end

x_range = linspace(x_range(1), x_range(2), bp.output_points.x);
y_range = linspace(y_range(1), y_range(2), bp.output_points.y);
z_range = linspace(z_range(1), z_range(2), bp.output_points.z);

%% Data prep
readi_group_count = 1;

tx_count = bp.dec_data_dim.z;
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

    pause(0.005);

    bp.readi_group_id = bp.readi_group_id + 1;
   
end

unloadlibrary('cuda_transfer')


%% Post processing 
volume_size = interleaved_volume_size;
volume_size(1) = interleaved_volume_size(1)/2;
high_res_volume = zeros(volume_size);

dynamic_range = 50;

threshold = 180;

processed_low_res = low_res_cells;
low_res_images = low_res_cells;
for g = 1:readi_group_count
    high_res_volume = high_res_volume + low_res_cells{g};
    processed_low_res{g} = process_volume(low_res_cells{g}, dynamic_range, threshold); 

    low_res_images{g} = squeeze(processed_low_res{g}(:,round(output_counts_xyz.y/2),:));
end
processed_volume = process_volume(high_res_volume, dynamic_range, threshold);

x_image = squeeze(processed_volume(:, round(output_counts_xyz.y/2),:)).';
y_image = squeeze(processed_volume(round(output_counts_xyz.x/2),:,:)).';

processed_volume = flip(processed_volume,3);



%% 
% volumeViewer(processed_volume)

%%
figure();
% sgtitle("FORCES Image");

subplot 121
imagesc(x_range * 1000, z_range * 1000, x_image);
title("No Filter")
colormap("gray");
axis('image');
% clim([-dynamic_range, 0]); 
colorbar;

subplot 122
imagesc(x_range * 1000, z_range * 1000, x_image);
title("With Filter")
colormap("gray");
axis('image');
% clim([-dynamic_range, 0]); 
colorbar;


% %
% figure();
% 
% % Main title
% % sgtitle(speed_str + ' Readi Low Resolution Images','FontSize',16);
% sgtitle('Stationary Readi Low Resolution Images','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(4, 4, i); 
%     colormap("gray");
%     imagesc((low_res_images{i}).');
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
%     axis image;
% 
%     % clim([-dynamic_range, 0]); 
% 
% 
% end
