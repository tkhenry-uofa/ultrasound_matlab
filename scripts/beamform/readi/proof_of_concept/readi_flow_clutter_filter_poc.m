clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');


if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

%% Data loading
data_path_root = "vrs_data/readi/better_flow/";

% dataset_name = "250424_MN32-5_flow_6mm_static_FORCES-TxColumn";
% dataset_name = "250424_MN32-5_flow_6mm_10_FORCES-TxColumn";
% dataset_name = "250424_MN32-5_flow_6mm_60_FORCES-TxColumn";
% dataset_name = "250424_MN32-5_flow_6mm_180_FORCES-TxColumn";
dataset_name = "250425_MN32-5_flow_6mm_10_FORCES-TxColumn";

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

data_file_range = 0:7;
data_file_paths = data_path + dataset_name + compose('_%02i.zst', data_file_range).';

bp = load_and_parse_bp(params_path);
frame_data = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);
frame_count = length(frame_data);

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

volume_ranges = [-10, 10; % X
                  0,  0; % Y
                  10, 25] / 1000; % Z

lateral_resolution = 0.00005;
axial_resolution = lateral_resolution;

bp.f_number = 0.5;

[bp, x_range, y_range, z_range] = configure_output_points(volume_ranges, lateral_resolution, axial_resolution, bp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readi data prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

readi_bp = bp;

readi_group_count = 16;
readi_group_data = cell(frame_count,readi_group_count);

readi_bp.readi_group_size = transmit_count/readi_group_count;
readi_bp.readi_group_id = 0;

% Break up transmits into Readi groups
for f=1:frame_count
    full_frame_data = frame_data{f};
    for i=1:readi_group_count
        end_tx = i * readi_bp.readi_group_size;
        start_tx = end_tx - (readi_bp.readi_group_size - 1);
        readi_group_data{f,i} = full_frame_data(:,start_tx:end_tx,:);
    end
end
readi_bp.rf_raw_dim(1) = readi_bp.readi_group_size * sample_count;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end

raw_images = cell(size(readi_group_data));

for f = 1:frame_count
    readi_bp.readi_group_id = 0;
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d, Group %d\n", f, g);
        raw_image = cuda_beamform_real(readi_group_data{f,g}, readi_bp);
        fprintf("Done\n");
        raw_images{f,g} = raw_image.';
        readi_bp.readi_group_id = readi_bp.readi_group_id + 1;
    end
end

if true, unloadlibrary('cuda_transfer'); end


fprintf("Beamforming Complete\n");

clear("readi_group_data");

%% Motion compensation

low_res_array = raw_images;

vessel_row_start = 102;
vessel_row_end = 210;
vessel_row_range = vessel_row_start:vessel_row_end;
vessel_col_range = 1:size(raw_image.',2);

svd_test;

%%
tic;
shifted_images = low_res_array(1,:);

forces_frame_id =7;
image_array = filtered_frame_cell(forces_frame_id,:);
ncc_block_matching;

elapsed = toc;
% phase_match_test;
fprintf("Compensation Time: %.3f\n",elapsed);

%% Processing
dynamic_range = 35;
threshold = 220;


shifted_image_raw = complex(zeros(size(low_res_array{1,1})));
forces_image_raw = shifted_image_raw;

processed_low_res = cell(1,readi_group_count);
processed_shifted_array = cell(1,readi_group_count);
for i = 1:readi_group_count
    processed_low_res{i} = process_volume(filtered_frame_cell{forces_frame_id,i}, dynamic_range,threshold);  
    forces_image_raw = forces_image_raw + low_res_array{forces_frame_id,i};

    processed_shifted_array{i} = process_volume(shifted_images{i},dynamic_range,threshold);
    shifted_image_raw = shifted_image_raw + shifted_images{i};
end



processed_forces_image = process_volume(forces_image_raw,dynamic_range,threshold);
processed_shifted_image = process_volume(shifted_image_raw,dynamic_range,threshold);

shifted_image_raw(1:vessel_row_start,:) = 0;
shifted_image_raw(vessel_row_end:end,:) = 0;
shifted_image_raw = shifted_image_raw ./ max(shifted_image_raw,[],"all");

% forces_image_raw(91:230,:) = 0;
forces_image_raw = forces_image_raw ./ max(forces_image_raw,[],"all");
combined_image_raw = forces_image_raw + shifted_image_raw;
processed_combined_image = process_volume(combined_image_raw,dynamic_range,threshold);

%% Average motion

motion_range = [1:5,11:16];
average_motion_array = zeros(size(processed_combined_image));
for i=motion_range
    
    motion_array = motion_cell{i};
    motion_array = imgaussfilt(motion_array,1.0);
    full_motion_array = imresize(motion_array, image_size, 'bilinear');
    
    full_motion_array(1:vessel_row_start,:,:) = 0;
    full_motion_array(vessel_row_end:end,:,:) = 0;

    if(i > reference_readi_group)
        full_motion_array(:,:,2) = -full_motion_array(:,:,2);
    end 
    average_motion_array = average_motion_array + full_motion_array;
end

average_motion_array = average_motion_array ./ length(motion_range);

full_x_motion = squeeze(average_motion_array(:,:,2));
full_z_motion = squeeze(average_motion_array(:,:,1));
full_motion = sqrt(full_x_motion.^2 + full_z_motion.^2);

x_motion_max = max(max(full_x_motion));

cross_section_col = 325;
average_flow_curve = sum(full_x_motion,2)./ length(x_range);
flow_velocity = abs(full_x_motion(:,cross_section_col));

profile_row_range = vessel_row_start:vessel_row_end+20;
average_flow_curve = average_flow_curve(profile_row_range);
flow_velocity = flow_velocity(profile_row_range);

%% Plot Flow Profiles
% figure();
% sgtitle("Axial Flow Profiles");
% subplot 121;
% plot(z_range(profile_row_range).*1000,average_flow_curve./10);
% title("Laterally-Averaged Velocity");
% ylim([0 1.2]);
% xlabel("Depth (mm)");
% ylabel("Velocity (cm/s)");
% 
% subplot 122
% plot(z_range(profile_row_range).*1000, flow_velocity./10);
% title("Velocity Profile at X=6.2mm");
% xlabel("Depth (mm)");
% ylabel("Velocity (cm/s)");

%% Plot Image Components 

% figure();
% 
% sgtitle('1.2 cm/s Peak Flow', 'FontSize',14)
% 
% subplot 221
% colormap("gray");
% imagesc(x_range*1000, z_range*1000, processed_forces_image);
% title("Forces Image");
% axis("image")
% ylabel("Depth (mm)");
% xlabel("Lateral Position (mm)");
% 
% subplot 222
% colormap("gray");
% imagesc(x_range*1000, z_range*1000, processed_shifted_image);
% title('Motion Compensated Speckle')
% axis("image")
% ylabel("Depth (mm)");
% xlabel("Lateral Position (mm)");
% 
% subplot 223
% colormap("gray");
% imagesc(x_range*1000, z_range*1000, processed_combined_image);
% title('Combined Image')
% axis("image")
% ylabel("Depth (mm)");
% xlabel("Lateral Position (mm)");
% 
% ax = subplot(2,2,4);
% imagesc(x_range*1000, z_range*1000, abs(full_x_motion)./10);
% title('Lateral Motion Map')
% colormap(ax,"hot");
% axis("image")
% colorbar('Position', [0.92 0.12 0.02 0.3]);
% ylabel("Depth (mm)");
% xlabel("Lateral Position (mm)");

%% Make Layered Image
base_image = processed_combined_image;
color_image = abs(full_x_motion) ./ 10;
alpha = 0.6;  % Opacity for the top image

composite_rgb = layer_heatmap(processed_combined_image, color_image, vessel_row_range, vessel_col_range, alpha );

fg = figure();
imagesc(x_range.*1000,z_range.*1000,composite_rgb);
fg.Colormap = hot(256);

cb = colorbar('Position', [0.92 0.15 0.02 0.7]);
cb.Ticks = linspace(0,1,7);
cb.TickLabels = linspace(0,1.2,7);

title("1.2 cm/s Lateral Flow", "FontSize",14)
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

