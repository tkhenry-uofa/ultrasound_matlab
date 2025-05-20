clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');
path(path, "..\simulations")

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path   = "vrs_data/readi/stage_motion/250423_MN32-5_reso_motion_FORCES-TxColumn/";
data_file_name = "250423_MN32-5_reso_motion_FORCES-TxColumn_06.zst";
params_path = data_path + '250423_MN32-5_reso_motion_FORCES-TxColumn.bp';

static_data_path   = "vrs_data/readi/stage_motion/250423_MN32-5_reso_motion_FORCES-TxColumn/";
static_data_file_name = "250423_MN32-5_reso_motion_FORCES-TxColumn_00.zst";
static_params_path = static_data_path + '250423_MN32-5_reso_motion_FORCES-TxColumn.bp';

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

data_file = fopen(data_path + data_file_name, "r");
raw_data = fread(data_file, '*uint8');
fclose(data_file); 
data = ornot_zstd_decompress_mex(raw_data);

static_data_file = fopen(static_data_path + static_data_file_name, "r");
raw_data = fread(static_data_file, '*uint8');
fclose(static_data_file); 
static_data = ornot_zstd_decompress_mex(raw_data);

flow_raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));
static_raw_bp = ornot_bp_load_mex(convertStringsToChars(static_params_path));


data = reshape(data, flow_raw_bp.raw_data_dim(1),flow_raw_bp.raw_data_dim(2));
static_data = reshape(static_data, flow_raw_bp.raw_data_dim(1),flow_raw_bp.raw_data_dim(2));

clear("raw_data");

%%
bp.decode          = flow_raw_bp.decode_mode;
bp.beamform_plane  = flow_raw_bp.beamform_mode;

% bp.rf_raw_dim = raw_bp.rf_raw_dim;
bp.rf_raw_dim.x      = flow_raw_bp.raw_data_dim(1);
bp.rf_raw_dim.y      = flow_raw_bp.raw_data_dim(2);
% bp.dec_data_dim    = raw_bp.dec_data_dim;
bp.dec_data_dim.x    = flow_raw_bp.decoded_data_dim(1);
bp.dec_data_dim.y    = flow_raw_bp.decoded_data_dim(2);
bp.dec_data_dim.z    = flow_raw_bp.decoded_data_dim(3);
bp.dec_data_dim.w    = flow_raw_bp.decoded_data_dim(4);

% Map transducer properties
bp.xdc_element_pitch = flow_raw_bp.transducer_element_pitch;
bp.xdc_transform     = flow_raw_bp.transducer_transform_matrix;

% Map channel and angle related fields
bp.channel_mapping  = flow_raw_bp.channel_mapping;
bp.transmit_angles  = flow_raw_bp.steering_angles;
bp.focal_depths     = flow_raw_bp.focal_depths;

% Map acoustic parameters
bp.speed_of_sound    = flow_raw_bp.speed_of_sound;
bp.center_frequency  = flow_raw_bp.center_frequency;
bp.sampling_frequency= flow_raw_bp.sampling_frequency;
bp.time_offset       = flow_raw_bp.time_offset;

% Map transmit mode
bp.transmit_mode     = flow_raw_bp.transmit_mode;


%%
% 0 = Forces, 1 = Uforces, 2 = Hercules
bp.das_shader_id = 0;

fc = bp.center_frequency;
fs = bp.sampling_frequency;

% The end of every channel has padded zeros not from any transmit
expected_signal_length = bp.dec_data_dim.x * 128;

shuffled_data = data(1:expected_signal_length, :);
shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(data, 2)];
shuffled_data = reshape(shuffled_data, shuffled_size);

static_shuffled_data = static_data(1:expected_signal_length, :);
static_shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(data, 2)];
static_shuffled_data = reshape(static_shuffled_data, static_shuffled_size);

start_depth = 3/1000;

start_sample = round(2 * start_depth * bp.sampling_frequency / bp.speed_of_sound);
tx_end_sample = start_sample -1;

shuffled_data(1:tx_end_sample,:,:) = 0;
% shuffled_data = shuffled_data(start_sample:end,:,:);

static_shuffled_data(1:tx_end_sample,:,:) = 0;
% static_shuffled_data = static_shuffled_data(start_sample:end,:,:);

% bp.dec_data_dim.x = bp.dec_data_dim.x - tx_end_sample;
% bp.rf_raw_dim.x = bp.rf_raw_dim.x - tx_end_sample * 256;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-25, 25]/1000;
y_range = [0, 0]/1000;
z_range = [5, 100]/1000;

bp.f_number = 1.0;

bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.0001;
axial_resolution = 0.0001;

bp.output_points.x = max(1,floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution ));
bp.output_points.y = max(1,floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution ));
bp.output_points.z = max(1,floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution ));
bp.output_points.w = 1; % Not used but needs to be in the struct

x_range = linspace(x_range(1), x_range(2), bp.output_points.x);
y_range = linspace(y_range(1), y_range(2), bp.output_points.y);
z_range = linspace(z_range(1), z_range(2), bp.output_points.z);

%% Data prep
readi_group_count = 1;
tx_count = bp.dec_data_dim.z;
readi_group_size = tx_count/readi_group_count;

bp.readi_group_size = readi_group_size;
bp.readi_group_id = 0;
group_data_cells = cell(1,readi_group_count);

% Break up transmits into Readi groups
for i=1:readi_group_count
    end_tx = i * bp.readi_group_size;
    start_tx = end_tx - (bp.readi_group_size - 1);

    slice = shuffled_data(:,start_tx:end_tx,:);
    slice = reshape(slice, [], bp.rf_raw_dim.y);

    static_slice = shuffled_data(:,start_tx:end_tx,:);
    static_slice = reshape(static_slice, [], bp.rf_raw_dim.y);

    group_data_cells{i} = slice;
end

readi_expected_signal_length = length(slice);
bp.rf_raw_dim.x = readi_expected_signal_length;


if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end

loadlibrary('cuda_transfer')


fprintf("Sending data\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Starting Beamforming\n");
output_counts_xyz.x = bp.output_points.x;
output_counts_xyz.y = bp.output_points.y;
output_counts_xyz.z = bp.output_points.z;

low_res_array = cell(1,readi_group_count);

% Complex volumes aren't supported so they're interleaved
interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];

bp.f_number = 0.5;

bp.readi_group_size = 128;
bp.focal_depths = static_raw_bp.focal_depths;
bp.time_offset = static_raw_bp.time_offset;
bp.rf_raw_dim.x = expected_signal_length;

calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);
fprintf("Static Forces\n");
volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

[~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
    pipe_name, smem_name, static_data, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
fprintf("Received Response\n")

volume = reshape(volume_ptr,interleaved_volume_size );

real_vol = volume(1:2:end,:,:);
im_vol = volume(2:2:end,:,:);

static_image_raw = squeeze(complex(real_vol, im_vol)).';

bp.readi_group_size = readi_group_size;
bp.focal_depths = flow_raw_bp.focal_depths;
bp.time_offset = flow_raw_bp.time_offset;
bp.rf_raw_dim.x = readi_expected_signal_length;

% Readi loop
for g = 1:readi_group_count
    calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

    fprintf("Readi group %i. \n", g);
    volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

    [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
        pipe_name, smem_name, group_data_cells{g}, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
    fprintf("Received Response\n")
    
    volume = reshape(volume_ptr,interleaved_volume_size );

    real_vol = volume(1:2:end,:,:);
    im_vol = volume(2:2:end,:,:);
    
    low_res_array{g} = squeeze(complex(real_vol, im_vol)).';

    pause(0.05);
    bp.readi_group_id = bp.readi_group_id + 1;
   
end


unloadlibrary('cuda_transfer')


%% Motion compensation

shifted_images = low_res_array;
ncc_block_matching;

% phase_match_test;
fprintf("motion compensation done\n");



%% Processing
dynamic_range = 40;
threshold = 200;

shifted_image_raw = complex(zeros(size(low_res_array{1})));
forces_image_raw = shifted_image_raw;
processed_low_res = low_res_array;
processed_shifted_array = low_res_array;
for i = 1:readi_group_count
    processed_low_res{i} = process_volume(low_res_array{i}, dynamic_range,threshold); 
    processed_shifted_array{i} = process_volume(shifted_images{i}, dynamic_range,threshold); 
    shifted_image_raw = shifted_image_raw + shifted_images{i};
    forces_image_raw = forces_image_raw + low_res_array{i};
end

static_image = process_volume(static_image_raw,dynamic_range,threshold);
forces_image = process_volume(forces_image_raw,dynamic_range,threshold);

processed_shifted_image = process_volume(shifted_image_raw,dynamic_range,threshold);
processed_shifted_image = imgaussfilt(processed_shifted_image,2);


%% Plotting

figure()

% sgtitle(speed_str + " Point Spread Function", 'FontSize',16)
% sgtitle("Stationary Point Spread Function", 'FontSize',16)

subplot 211
colormap("gray");
imagesc(x_range*1000, z_range*1000, static_image);
subtitle('Static', 'FontSize',14);
axis("image")
clim([-dynamic_range, 0]); 
xlabel("Lateral Distance (mm)");
ylabel("Axial Distance (mm)");

subplot 212
colormap("gray");
imagesc(x_range*1000, z_range*1000, forces_image);
subtitle('Motion', 'FontSize',14)
axis("image")
clim([-dynamic_range, 0]); 
xlabel("Lateral Distance (mm)");
ylabel("Axial Distance (mm)");



colorbar('Position', [0.92 0.15 0.02 0.7]);

figure();

% Main title
% sgtitle(speed_str + ' Uncompensated Low Resolution Images','FontSize',16);
sgtitle('Uncompensated','FontSize',16);

% Loop through the images and display each in the 2x8 grid
for i = 1:readi_group_count
    subplot(4, 4, i); 
    colormap("gray");
    imagesc(x_range*1000, z_range*1000,(processed_low_res{i}));
    title(['Image ' num2str(i)]);
    xlabel("mm");
    ylabel("mm");

    clim([-dynamic_range, 0]); 


end
colorbar('Position', [0.92 0.15 0.02 0.7]);

% figure();
% 
% % Main title
% % sgtitle(speed_str + ' Compensated Low Resolution Images','FontSize',16);
% sgtitle('Phase Match Warped','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(2, 4, i); 
%     colormap("gray");
%     imagesc(x_range*1000, z_range*1000,(processed_shifted_array{i}));
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end
% colorbar('Position', [0.92 0.15 0.02 0.7]);

