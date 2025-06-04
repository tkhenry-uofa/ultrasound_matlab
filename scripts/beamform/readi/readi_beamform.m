clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "vrs_data/match_filter_test/";

dataset_name = "250514_MN32-5_reso_FORCES-TxColumn";
data_file_range = 0:0;
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
volume_ranges = [-30, 30; % X
                  0,  0; % Y
                  5, 80] / 1000; % Z

lateral_resolution = 0.0002;
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

if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end

raw_images = cell(size(readi_group_data));

for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d, Group %d\n", f, g);
        raw_image = cuda_beamform(readi_group_data{f,g}, readi_bp);
        fprintf("Done\n");
        raw_images{f,g} = raw_image;
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
end

if true, unloadlibrary('cuda_transfer'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_frame = 1;
dynamic_range = 50;
threshold = 300;
power_image = false;

image_dims = size(raw_images{plot_frame,1});

processed_image_array = cell(1,readi_group_count);
forces_image_raw = zeros(image_dims);

for g = 1:readi_group_count
    raw_image = raw_images{plot_frame,g};
    processed_image = process_volume(raw_image,dynamic_range,threshold,power_image).';
    processed_image_array{g} = processed_image;
    forces_image_raw = forces_image_raw + raw_image;
end

processed_forces_image = process_volume(forces_image_raw,dynamic_range,threshold,power_image).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();
plot_bmode(processed_forces_image, x_range, z_range, "Forces Image");
colorbar;
plot_image_grid(processed_image_array, [4, 4], x_range, z_range, "Readi Groups");