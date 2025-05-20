clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');
addpath("C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\simulations")
path(path, "..\simulations\motion_comp")

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

% data_path   = "vrs_data/readi/stage_motion/250410_MN32-5_reso_axial_motion_FORCES-TxRow/";
% data_file_name = "250410_MN32-5_reso_axial_motion_FORCES-TxRow_01.zst";
% params_path = data_path + '250410_MN32-5_reso_axial_motion_FORCES-TxRow.bp';
% 
% data_path   = "vrs_data/readi/stage_motion/250411_C1_reso_hand_lateral_FORCES-TxColumn/";
% data_file_name = "250411_C1_reso_hand_lateral_FORCES-TxColumn_02.zst";
% params_path = data_path + '250411_C1_reso_hand_lateral_FORCES-TxColumn.bp';

% data_path   = "vrs_data/readi/flow/250423_MN32-5_reso_motion_FORCES-TxColumn/";
% data_file_name = "250423_MN32-5_reso_motion_FORCES-TxColumn_06.zst";
% params_path = data_path + '250423_MN32-5_reso_motion_FORCES-TxColumn.bp';

data_path   = "vrs_data/readi/flow/250423_MN32-5_flow_4mm_360_stopping_FORCES-TxColumn/";
data_file_name = "250423_MN32-5_flow_4mm_360_stopping_FORCES-TxColumn_05.zst";
params_path = data_path + '250423_MN32-5_flow_4mm_360_stopping_FORCES-TxColumn.bp';

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

data_file = fopen(data_path + data_file_name, "r");
raw_data = fread(data_file, '*uint8');
fclose(data_file); 
data = ornot_zstd_decompress_mex(raw_data);

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));

data = reshape(data, raw_bp.raw_data_dim(1),raw_bp.raw_data_dim(2));

clear("raw_data");

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

fc = bp.center_frequency;
fs = bp.sampling_frequency;




% The end of every channel has padded zeros not from any transmit
expected_signal_length = bp.dec_data_dim.x * bp.dec_data_dim.z;

shuffled_data = data(1:expected_signal_length, :);
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

% 250411_C1_reso_hand_lateral_FORCES-TxColumn_02
% x_range = [-30, 30]/1000;
% y_range = [0, 0]/1000;
% z_range = [15, 110]/1000;
% lateral_resolution = 0.00015;
% axial_resolution = 0.00015;
% bp.f_number = 0.0;

% 250410_MN32-5_reso_axial_motion_FORCES-TxRow_01
% x_range = [-20, 20]/1000;
% y_range = [0, 0]/1000;
% z_range = [15, 80]/1000;
% lateral_resolution = 0.00005;
% axial_resolution = 0.00005;
% bp.f_number = 0.5;

% % 250423_MN32-5_reso_motion_2_FORCES-TxColumn
% x_range = [-30, 30]/1000;
% y_range = [0, 0]/1000;
% z_range = [10, 90]/1000;
% lateral_resolution = 0.00005;
% axial_resolution = 0.00005;
% bp.f_number = 0.5;

% 250423_MN32-5_reso_motion_2_FORCES-TxColumn
x_range = [-15, 15]/1000;
y_range = [0, 0]/1000;
z_range = [10, 30]/1000;
lateral_resolution = 0.00005;
axial_resolution = 0.00005;
bp.f_number = 0.5;


bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

bp.output_points.x = max(1,floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution ));
bp.output_points.y = max(1,floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution ));
bp.output_points.z = max(1,floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution ));
bp.output_points.w = 1; % Not used but needs to be in the struct

x_range = linspace(x_range(1), x_range(2), bp.output_points.x);
y_range = linspace(y_range(1), y_range(2), bp.output_points.y);
z_range = linspace(z_range(1), z_range(2), bp.output_points.z);

%% Data prep
readi_group_count = 8;
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
    if i == 1
        slice(:,1,:) = zeros(size(slice(:,1,:)));
        % slice(:,1,:) = slice(:,1,:)./sqrt(64);
    end
    slice = reshape(slice, [], bp.rf_raw_dim.y);

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

    pause(0.5);
    bp.readi_group_id = bp.readi_group_id + 1;
   
end


unloadlibrary('cuda_transfer')


%% Motion compensation
tic;
shifted_images = low_res_array;
ncc_block_matching;
% reverse_block_matching;

elapsed = toc;
% phase_match_test;
fprintf("Compensation Time: %.3f\n",elapsed);

%% Processing
dynamic_range = 40;
threshold = 400;


shifted_image_raw = complex(zeros(size(low_res_array{1})));
forces_image_raw = shifted_image_raw;

processed_low_res = cell(size(low_res_array));
processed_shifted_array = cell(size(low_res_array));
for i = 1:readi_group_count
    processed_low_res{i} = process_volume(low_res_array{i}, dynamic_range,threshold); 
    processed_shifted_array{i} = process_volume(shifted_images{i}, dynamic_range,threshold + 40); 
    
    forces_image_raw = forces_image_raw + low_res_array{i};
    shifted_image_raw = shifted_image_raw + shifted_images{i};
end

forces_image = process_volume(forces_image_raw,dynamic_range,threshold);
processed_shifted_image = process_volume(shifted_image_raw,dynamic_range,threshold + 40);

% processed_shifted_image = imgaussfilt(processed_shifted_image,1);


%% Plotting

figure();

% Main title
% sgtitle(speed_str + ' Uncompensated Low Resolution Images','FontSize',16);
% sgtitle('Readi Low Resolution Images','FontSize',16);

% Loop through the images and display each in the 2x8 grid
for i = 1:readi_group_count
    subplot(2, 4, i); 
    colormap("gray");
    imagesc(x_range*1000, z_range*1000,(processed_low_res{i}));
    title(['Image ' num2str(i)]);
    % xlabel("mm");
    % ylabel("mm");

    clim([-dynamic_range, 0]); 


end
colorbar('Position', [0.92 0.15 0.02 0.7]);

figure();

% Main title
% sgtitle(speed_str + ' Compensated Low Resolution Images','FontSize',16);
sgtitle('Warped','FontSize',16);

% Loop through the images and display each in the 2x8 grid
for i = 1:readi_group_count
    subplot(2, 4, i); 
    colormap("gray");
    imagesc(x_range*1000, z_range*1000,(processed_shifted_array{i}));
    title(['Image ' num2str(i)]);
    xlabel("mm");
    ylabel("mm");

    clim([-dynamic_range, 0]); 


end
colorbar('Position', [0.92 0.15 0.02 0.7]);

figure();

% sgtitle(speed_str + " Point Spread Function", 'FontSize',16)
% sgtitle("Stationary Point Spread Function", 'FontSize',16)

subplot 121
colormap("gray");
imagesc(x_range*1000, z_range*1000, forces_image);
% subtitle('Uncompensated', 'FontSize',14)
axis("image")
% clim([-dynamic_range, 0]); 
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");

subplot 122
colormap("gray");
imagesc(x_range*1000, z_range*1000, processed_shifted_image);
% subtitle('Compenstated', 'FontSize',14)
axis("image")
% clim([-dynamic_range, 0]); 
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");

colorbar('Position', [0.92 0.15 0.02 0.7]);

