clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "vrs_data/readi/beating_heart/";

% dataset_name = "250521_MN32-5_beating_heart_static_FORCES-TxColumn";
% dataset_name = "250521_MN32-5_beating_heart_static_FORCES-TxRow";
% dataset_name = "250521_MN32-5_beating_heart_FORCES-TxColumn";
dataset_name = "250521_MN32-5_beating_heart_FORCES-TxRow";

data_file_range = 0:15;

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

bp = load_and_parse_bp(params_path);
frame_data = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);
frame_count = length(frame_data);

sample_count = single(bp.dec_data_dim(1));
channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;

tx_region = 3/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data, bp.rf_raw_dim] = crop_and_blank_tx(frame_data, bp, tx_region);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_ranges = [-40, 40; % X
                  0,  0; % Y
                  30, 120] / 1000; % Z

lateral_resolution = 0.0002;
axial_resolution = lateral_resolution;

bp.f_number = 0.5;

[bp, x_range, y_range, z_range] = configure_output_points(volume_ranges, lateral_resolution, axial_resolution, bp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_images = cell(1,frame_count);

for i = 1:frame_count
    fprintf("Beamforming Image %d\n", i);
    raw_image = cuda_beamform(frame_data{i}, bp);
    fprintf("Done\n");
    raw_images{i} = raw_image;
end

if true, unloadlibrary('cuda_transfer_matlab'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamic_range = 50;
threshold = 300;
power_image = false;

processed_image = process_volume(raw_images{1},dynamic_range,threshold,power_image).';

processed_image_array = cell(1,frame_count);
for i = 1:frame_count
    processed_image_array{i} = process_volume(raw_images{i},dynamic_range,threshold,power_image).';
end

figure();
plot_bmode(processed_image, x_range, z_range);
colorbar;

plot_image_grid(processed_image_array, [4, 4], x_range, z_range, "Sequence");


