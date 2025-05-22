clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');


if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "../vrs_data/match_filter_test/";

dataset_name = "250514_MN32-5_reso_FORCES-TxColumn";
% dataset_name = "250514_MN32-5_reso_FORCES-TxColumn-Chirp-2e-05";
% dataset_name = "250514_MN32-5_cyst_FORCES-TxColumn";
% dataset_name = "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-5e-06";

% dataset_name = "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-2e-05";

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

data_file_num = 0;
data_file_name = dataset_name + sprintf('_%02i.zst', data_file_num);

%% Load Data and BP
data_file_h = fopen(data_path + data_file_name, "r");
raw_data = fread(data_file_h, '*uint8');
fclose(data_file_h); 
data = ornot_zstd_decompress_mex(raw_data);

bp = load_and_parse_bp(params_path);
data = reshape(data, bp.rf_raw_dim(1),bp.rf_raw_dim(2));

sample_count = single(bp.dec_data_dim(1));
rx_channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;


%% Match Filter setup

chirp_length = 2e-5;
% chirp_length = 5e-6;

load(data_path + "postVsx.mat","Scan");

match_filter = generate_chirp_filter(fc,fs,Scan.Die.Bandwidth,chirp_length);

%% Remove trailing zeros and blank transmit

% The end of every channel has padded zeros not from any transmit
expected_length = sample_count * transmit_count;

cropped_data = data(1:expected_length, :);

bp.rf_raw_dim(1) = expected_length;

% Breaks the first dim up into sample_count and transmit_count
% Doesn't reorder anything in memory so we can pass this right to the lib
shuffled_size = [sample_count, transmit_count, bp.rf_raw_dim(2)];
cropped_data = reshape(cropped_data, shuffled_size);

start_depth = 5/1000;
start_sample = round(2 * start_depth * bp.sampling_frequency / bp.speed_of_sound);
tx_end_sample = start_sample-1;

cropped_data(1:tx_end_sample,:,:) = 0;


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Beamforming

if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end
loadlibrary('cuda_transfer')

%% Unfiltered
bp.filter_length = 0;
bp.match_filter = single(zeros(1024,1));

fprintf("Beamforming Unfiltered\n");
unfiltered_raw_img = cuda_beamform_real(cropped_data, bp);
fprintf("Done\n");

%% Filtered
filtered_bp = bp;
filtered_bp.filter_length = length(match_filter);
filtered_bp.match_filter(1:length(match_filter)) = match_filter;
filtered_bp.time_offset = bp.time_offset + (length(match_filter)-1) * 1.0 / bp.sampling_frequency;

fprintf("Beamforming Filtered\n");
filtered_raw_img = cuda_beamform_real(cropped_data,filtered_bp);
fprintf("Done\n");

unloadlibrary('cuda_transfer')


%% Post processing 
volume_size = size(filtered_raw_img);

dynamic_range = 50;

threshold = 800000;

processed_unfiltered = process_volume(unfiltered_raw_img,dynamic_range,threshold,false).';
processed_filtered = process_volume(filtered_raw_img,dynamic_range,threshold,false).';


%% 
% volumeViewer(processed_volume)

%%
figure();
% sgtitle("FFT, No Pulse Delay");

subplot 121
imagesc(x_range * 1000, z_range * 1000, processed_unfiltered);
title("No Filter")
colormap("gray");
axis('image');
% clim([-dynamic_range, 0]); 
colorbar;

subplot 122
imagesc(x_range * 1000, z_range * 1000, processed_filtered);
title("Filtered")
colormap("gray");
axis('image');
% clim([-dynamic_range, 0]); 
colorbar;

