clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\client_lib\output")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

data_path_root   = "vrs_data/readi/beating_heart/";

% dataset_name = "250521_MN32-5_beating_heart_static_FORCES-TxColumn";
% dataset_name = "250521_MN32-5_beating_heart_static_FORCES-TxRow";
% dataset_name = "250521_MN32-5_beating_heart_FORCES-TxColumn";
dataset_name = "250521_MN32-5_beating_heart_FORCES-TxRow";
data_file_range = 0:0;

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

tx_region = 0/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data, bp.rf_raw_dim] = crop_and_blank_tx(frame_data, bp, tx_region);

bp.data_type = int32(BeamformerDataType.I16); % int16


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_ranges = [-20, 20; % X
                  0,  0; % Y
                  50, 100] / 1000; % Z

lateral_resolution = 0.0005;
axial_resolution = lateral_resolution;

bp.f_number = 1.0;

[bp, x_range, y_range, z_range] = configure_output_points(volume_ranges, lateral_resolution, axial_resolution, bp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readi data prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readi_group_count = 16;


readi_bp = bp;

readi_bp.readi_group_count = readi_group_count;
readi_group_data = cell(frame_count,readi_group_count);
readi_group_size = bp.dec_data_dim(3) / readi_group_count;

readi_bp.readi_group_id = 0;

% Break up transmits into Readi groups
for f=1:frame_count
    full_frame_data = frame_data{f};
    for i=1:readi_group_count
        end_tx = i * readi_group_size;
        start_tx = end_tx - (readi_group_size - 1);
        readi_group_data{f,i} = full_frame_data(:,start_tx:end_tx,:);
    end
end

sample_count = bp.dec_data_dim(1);
readi_bp.rf_raw_dim(1) = readi_group_size * sample_count;
readi_bp.dec_data_dim(3) = readi_group_size;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if true && libisloaded('cuda_transfer_matlab'), unloadlibrary('cuda_transfer_matlab'); end

raw_images = cell(size(readi_group_data));

for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d\n", g );
        raw_image = cuda_beamform(readi_group_data{f,g}, readi_bp);
        fprintf("Done\n");
        raw_images{f,g} = raw_image;
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
    
end

if true && libisloaded('cuda_transfer_matlab'), unloadlibrary('cuda_transfer_matlab'); end

fprintf("Beamforming Complete\n");

% clear("readi_group_data");


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

