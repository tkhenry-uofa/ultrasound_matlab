clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end


data_path_root   = "vrs_data/readi/stage_motion/";
dataset_name = "250423_MN32-5_reso_motion_FORCES-TxColumn";
data_file_range = 0:0;

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

bp = load_and_parse_bp(params_path);
frame_data = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);
frame_count = length(frame_data);

sample_count = single(bp.dec_data_dim(1));
rx_channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;

tx_region = 2.5/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data, bp.rf_raw_dim] = crop_and_blank_tx(frame_data, bp, tx_region);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_ranges = [-20, 20;  % X
                  0,  0;   % Y
                 10,  80] / 1000;
resolution = 0.0001;

bp.f_number = 0.5;
[bp, x_range, y_range, z_range] = configure_output_points(volume_ranges, resolution, resolution, bp);

image_size = [length(z_range), length(x_range)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readi data prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

readi_group_count = 16;
[readi_group_data, readi_bp] = readi_data_breakup(frame_data, bp, readi_group_count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end

raw_images = cell(size(readi_group_data));

for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d, Group %d\n", f, g);
        current_raw_image = cuda_beamform(readi_group_data{f,g}, readi_bp);
        fprintf("Done\n");
        raw_images{f,g} = current_raw_image.';
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
end

if true, unloadlibrary('cuda_transfer'); end



%% Motion compensation
tic;
reference_image_id = 5;
patch_size = 32;   
v_margin = 16;
h_margin = 16;
motion_cell_size = 4; 
var_min = 0.005;
corr_threshold = 1.05;
corr_dr = 35;
img_threshold = 3000;

search_area = [ 1, 1; image_size];

shifted_images = forward_block_matching( ...
    raw_images, reference_image_id, patch_size, motion_cell_size, ...
    search_area,v_margin, h_margin, corr_threshold, var_min, corr_dr, img_threshold);

elapsed = toc;
% phase_match_test;
fprintf("Compensation Time: %.3f\n",elapsed);

%% Processing
plot_frame = 1;
dynamic_range = 50;
threshold = 300;
power_image = false;

image_dims = size(raw_images{plot_frame,1});

shifted_image_raw = complex(zeros(size(raw_images{plot_frame,1})));
forces_image_raw = shifted_image_raw;

processed_low_res = cell(1,readi_group_count);
processed_shifted_array = processed_low_res;
for i = 1:readi_group_count

    current_raw_image = raw_images{plot_frame,i};
    current_shifted_image = shifted_images{i};
    processed_low_res{i} = process_volume(current_raw_image, dynamic_range,threshold); 
    processed_shifted_array{i} = process_volume(current_shifted_image, dynamic_range,threshold + 40); 
    
    forces_image_raw = forces_image_raw + current_raw_image;
    shifted_image_raw = shifted_image_raw + current_shifted_image;
end

forces_image = process_volume(forces_image_raw,dynamic_range,threshold);
processed_shifted_image = process_volume(shifted_image_raw,dynamic_range,threshold + 40);

% processed_shifted_image = imgaussfilt(processed_shifted_image,1);


%% Plotting

plot_image_grid(processed_low_res, [2,4],x_range, z_range, "Readi Low Resolution Images");
plot_image_grid(processed_shifted_array, [2,4],x_range, z_range, "Motion Compensated Low Reolution Images");

figure();
subplot 121
plot_bmode(forces_image, x_range, z_range, "Uncompensated FORCES")

subplot 122
plot_bmode(processed_shifted_image, x_range, z_range, "Compensated FORCES")

