clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "vrs_data/match_filter_test/";

dataset_name = "250514_MN32-5_reso_FORCES-TxColumn";
data_file_range = 0:3;

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

data_file_paths = data_path + dataset_name + compose('_%02i.zst', data_file_range).';

[bp, arrays] = load_and_parse_bp(params_path);
bp.channel_mapping = arrays.channel_mapping;
bp.focal_depths = arrays.focal_depths;


total_frames = length(data_file_paths);
frame_data = cell(1,total_frames);
for i = 1:total_frames
    data_file = fopen(data_file_paths(i), "r");
    raw_data = fread(data_file, '*uint8');
    data = ornot_zstd_decompress_mex(raw_data);
    frame_data{i} = reshape(data, bp.rf_raw_dim(1),bp.rf_raw_dim(2));
    fclose(data_file); 
end



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
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end

raw_images = cell(1,total_frames);

for i = 1:total_frames
    fprintf("Beamforming Image %d\n", i);
    raw_image = cuda_beamform_real(frame_data{i}, bp);
    fprintf("Done\n");
    raw_images{i} = raw_image;
end

%% Post processing 
dynamic_range = 50;
threshold = 300;
power_image = false;

processed_image = process_volume(raw_images{1},dynamic_range,threshold,power_image).';

processed_image_array = cell(1,total_frames);
for i = 1:total_frames
    processed_image_array{i} = process_volume(raw_images{i},dynamic_range,threshold,power_image).';
end

figure();
plot_bmode(processed_image, x_range, z_range);
colorbar;

plot_image_grid(processed_image_array, x_range, z_range, [2, 2],"Sequence");

if true, unloadlibrary('cuda_transfer'); end

