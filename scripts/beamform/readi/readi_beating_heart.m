clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\client_lib\output")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');


if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

%% Data loading
data_path_root = "vrs_data/readi/beating_heart/";

% dataset_name = "250521_MN32-5_beating_heart_FORCES-TxRow";
dataset_name = "250521_MN32-5_beating_heart_FORCES-TxColumn";

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

data_file_range = 0:15;
data_file_paths = data_path + dataset_name + compose('_%02i.zst', data_file_range).';

bp = load_and_parse_bp(params_path);
frame_data = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);
frame_count = length(frame_data);

sample_count = single(bp.dec_data_dim(1));
rx_channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;

tx_region = 0/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data, bp.rf_raw_dim] = crop_and_blank_tx(frame_data, bp, tx_region);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volume_ranges = [-40, 40; % X
                  0,  0; % Y
                  40, 120] / 1000; % Z

lateral_resolution = 0.00015;
axial_resolution = lateral_resolution;

bp.f_number = 0.5;

[bp, x_range, y_range, z_range] = configure_output_points(volume_ranges, lateral_resolution, axial_resolution, bp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readi data prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readi_group_count = 16;

[readi_group_data, readi_bp] = readi_data_breakup(frame_data, bp, readi_group_count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_images = cell(size(readi_group_data));

total_duration = 0;
for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d, Group %d\n", f, g);
        tic;
        raw_image = cuda_beamform(readi_group_data{f,g}, readi_bp);
        duration = toc;
        fprintf("Elapsed time: %4.3f seconds.\n", duration);
        total_duration = total_duration + duration;
        raw_images{f,g} = raw_image.';
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
end

fprintf("Beamforming complete, total duration: %4.3f seconds.\n", total_duration);

clear("readi_group_data");

%% SVD Filtering
tissue_sv_cutoff = 30;
noise_sv_cutoff = 50;
filtered_images = svd_filter(raw_images,tissue_sv_cutoff,noise_sv_cutoff);


%% Processing
plot_frame = 5;

shifted_images = raw_images(plot_frame,:);

dr = 45;
tr = 220;

p_forces = cell(1,numel(raw_images));
p_low_res = cell(size(raw_images));
for f=1:frame_count
    [p_forces{f}, p_low_res(f,:)] = process_readi(raw_images(f,:),dr,tr,false);
end

p_filt_forces = cell(1,numel(raw_images));
p_filt_low_res = cell(size(raw_images));
for f=1:frame_count
    [p_filt_forces{f}, p_filt_low_res(f,:)] = process_readi(filtered_images(f,:),dr,tr,false);
end


%% Plot Image Components 

figure();
subplot 121
plot_bmode(p_forces{plot_frame}, x_range, z_range, "Forces");
subplot 122
plot_bmode(p_filt_forces{plot_frame}, x_range, z_range, "Filtered Forces");

plot_image_grid(p_low_res(plot_frame,:), [4,4], x_range, z_range, "Unfiltered Readi Images");

plot_image_grid(p_filt_low_res(plot_frame,:), [4,4], x_range, z_range, "Filtered Readi Images");


%%
animate_frames(reshape(p_filt_low_res.',[],1),64);
% animate_frames(reshape(p_low_res.',[],1),128);

