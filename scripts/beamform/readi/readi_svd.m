clear all;

addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "vrs_data/readi/better_flow/";

dataset_name = "250425_MN32-5_flow_6mm_30_half_f_FORCES-TxColumn";
% dataset_name = "250425_MN32-5_flow_6mm_60_FORCES-TxColumn";
data_file_range = 0:15;
frame_count = length(data_file_range);

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

bp = load_and_parse_bp(params_path);

frame_data = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);

sample_count = single(bp.dec_data_dim(1));
rx_channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;

tx_region = 5/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data, bp.rf_raw_dim] = crop_and_blank_tx(frame_data, bp, tx_region);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_ranges = [-15, 15; % X
                  0,  0; % Y
                  10, 25] / 1000; % Z

lateral_resolution = 0.00005;
axial_resolution = lateral_resolution;

bp.f_number = 1.0;

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

for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        tic;
        fprintf("Beamforming Image %d, Group %d\n", f, g);
        raw_image = cuda_beamform(readi_group_data{f,g}, readi_bp);
        toc
        raw_images{f,g} = raw_image.';
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
end

if true, unloadlibrary('cuda_transfer_matlab'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SVD
tissue_sv_cutoff = 30;
noise_sv_cutoff = 80;
tissue_start = 2;
filtered_images = svd_filter(raw_images,tissue_sv_cutoff,noise_sv_cutoff, tissue_start);

%% Plotting
plot_frame = 1;

shifted_images = raw_images(plot_frame,:);

dr = 40;
tr = 220;

p_forces = cell(1,size(raw_images,1));
p_low_res = cell(size(raw_images));
for f=1:frame_count
    [p_forces{f}, p_low_res(f,:)] = process_readi(raw_images(f,:),dr,tr,false);
end

p_filt_forces = cell(1,size(raw_images,1));
p_filt_low_res = cell(size(raw_images));
for f=1:frame_count
    [p_filt_forces{f}, p_filt_low_res(f,:)] = process_readi(filtered_images(f,:),dr,tr,false);
end


%% Plot Image Components 

% figure();
% subplot 121
% plot_bmode(p_forces{plot_frame}, x_range, z_range, "Forces");
% subplot 122
% plot_bmode(p_filt_forces{plot_frame}, x_range, z_range, "Filtered Forces");
% plot_image_grid(p_low_res(plot_frame,:), [4,4], x_range, z_range, "Unfiltered Readi Images");
% plot_image_grid(p_filt_low_res(plot_frame,:), [4,4], x_range, z_range, "Filtered Readi Images");


%%
% animate_frames(p_forces,8);
% animate_frames(reshape(p_filt_low_res.',[],1),100);
% imagej_animate_frames(reshape(p_filt_low_res.',[],1));

frame_cell = reshape(p_filt_low_res.',[],1);
frame_count = numel(frame_cell);
image_size = size(frame_cell{1}');
dr = -1 * min(frame_cell{1},[],"all");

frame_array_size = [image_size, frame_count];
frame_array = uint8(zeros(frame_array_size,"uint8"));

for i = 1:frame_count
    frame = frame_cell{i};
    frame_array(:,:,i) =  uint8((frame + dr) * 255/dr)';
end

%%

if ~exist("IJM", "var")
    addpath("C:\Users\tkhen\Fiji\scripts");
    ImageJ;
end

IJM.show('frame_array')
