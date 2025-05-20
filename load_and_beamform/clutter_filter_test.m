clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');
addpath("C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\simulations")
path(path, "..\simulations\motion_comp")

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

%% Data loading

data_range = 0:15;
frame_count = length(data_range);

data_root = "vrs_data/readi/better_flow/";

% data_file_name = "250424_MN32-5_flow_6mm_static_FORCES-TxColumn";
% data_file_name = "250424_MN32-5_flow_6mm_10_FORCES-TxColumn";
% data_file_name = "250424_MN32-5_flow_6mm_60_FORCES-TxColumn";
% data_file_name = "250424_MN32-5_flow_6mm_180_FORCES-TxColumn";

data_file_name = "250425_MN32-5_flow_6mm_10_FORCES-TxColumn";

data_folder   =  data_root + data_file_name + "/";
data_paths = data_folder + data_file_name + compose('_%02i.zst', data_range).';
params_path = data_folder + data_file_name + '.bp';
% params_path = data_folder + 'parameters_fixed.bp';

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));

frame_data_cell = cell(1,frame_count);
for i = 1:frame_count
    data_file = fopen(data_paths(i), "r");
    raw_data = fread(data_file, '*uint8');
    data = ornot_zstd_decompress_mex(raw_data);
    frame_data_cell{i} = reshape(data, raw_bp.raw_data_dim(1),raw_bp.raw_data_dim(2));
    fclose(data_file); 
end

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

prf = 1000;

expected_signal_length = bp.dec_data_dim.x * bp.dec_data_dim.z;

cropped_data_cell = cell(1,frame_count);
% The end of every channel has padded zeros not from any transmit, remove
% them and optionally blank out the start
for f = 1:frame_count

    cropped_data = frame_data_cell{f};
    cropped_data = cropped_data(1:expected_signal_length, :);
    
    % Blank out the start of the signal if the transmit is too loud
    cropped_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, size(cropped_data, 2)];
    cropped_data = reshape(cropped_data, cropped_size);
    start_depth = 5/1000;
    start_sample = round(2 * start_depth * bp.sampling_frequency / bp.speed_of_sound);
    tx_end_sample = start_sample -1;
    cropped_data(1:tx_end_sample,:,:) = 0;

    cropped_data_cell{f} = cropped_data;
end



% shuffled_data = shuffled_data(start_sample:end,:,:);
% 
% bp.dec_data_dim.x = bp.dec_data_dim.x - tx_end_sample;
% bp.rf_raw_dim.x = bp.rf_raw_dim.x - tx_end_sample * 256;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_range = [-10, 10]/1000;
y_range = [0, 0]/1000;
z_range = [10, 25]/1000;
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
readi_group_count = 16;
tx_count = bp.dec_data_dim.z;
readi_group_size = tx_count/readi_group_count;

bp.readi_group_size = readi_group_size;
bp.readi_group_id = 0;
group_data_cells = cell(frame_count,readi_group_count);

% Break up transmits into Readi groups
for f=1:frame_count
    
    frame_data = cropped_data_cell{f};
    for i=1:readi_group_count
        end_tx = i * bp.readi_group_size;
        start_tx = end_tx - (bp.readi_group_size - 1);
    
        slice = frame_data(:,start_tx:end_tx,:);
        if i == 1
            slice(:,1,:) = zeros(size(slice(:,1,:)));
        %     slice(:,1,:) = slice(:,1,:)./sqrt(64);
        end
        slice = reshape(slice, [], bp.rf_raw_dim.y);
    
        group_data_cells{f,i} = slice;
    end
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

low_res_array = cell(frame_count,readi_group_count);

% Complex volumes aren't supported so they're interleaved
interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];

% Readi loop
for f = 1:frame_count
    fprintf("\nFrame %i/%i: \n", f,frame_count);
    for g = 1:readi_group_count
        calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);
    
        fprintf("\tReadi group %02i", g);
        volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));
    
        [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
            pipe_name, smem_name, group_data_cells{f,g}, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
        fprintf("........Received Response\n")
        
        volume = reshape(volume_ptr,interleaved_volume_size );
    
        real_vol = volume(1:2:end,:,:);
        im_vol = volume(2:2:end,:,:);
        
        low_res_array{f,g} = squeeze(complex(real_vol, im_vol)).';
    
        pause(0.05);
        bp.readi_group_id = bp.readi_group_id + 1;
       
    end
    bp.readi_group_id = 0;
end

clear('data');
clear('cropped_data_cell')
clear('frame_data_cell')

unloadlibrary('cuda_transfer')

%% Motion compensation

vessel_row_start = 102;
vessel_row_end = 210;
vessel_row_range = vessel_row_start:vessel_row_end;

svd_test;

%%
tic;
shifted_images = low_res_array(1,:);

forces_frame_id =7;
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

%%
figure();
sgtitle("Axial Flow Profiles");
subplot 121;
plot(z_range(profile_row_range).*1000,average_flow_curve./10);
title("Laterally-Averaged Velocity");
ylim([0 1.2]);
xlabel("Depth (mm)");
ylabel("Velocity (cm/s)");

subplot 122
plot(z_range(profile_row_range).*1000, flow_velocity./10);
title("Velocity Profile at X=6.2mm");
xlabel("Depth (mm)");
ylabel("Velocity (cm/s)");

%% Plotting

figure();

sgtitle('1.2 cm/s Peak Flow', 'FontSize',14)

subplot 221
colormap("gray");
imagesc(x_range*1000, z_range*1000, processed_forces_image);
title("Forces Image");
axis("image")
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

subplot 222
colormap("gray");
imagesc(x_range*1000, z_range*1000, processed_shifted_image);
title('Motion Compensated Speckle')
axis("image")
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

subplot 223
colormap("gray");
imagesc(x_range*1000, z_range*1000, processed_combined_image);
title('Combined Image')
axis("image")
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

ax = subplot(2,2,4);
imagesc(x_range*1000, z_range*1000, abs(full_x_motion)./10);
title('Lateral Motion Map')
colormap(ax,"hot");
axis("image")
colorbar('Position', [0.92 0.12 0.02 0.3]);
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

%% Layered Image
base_image = processed_combined_image;
color_image = abs(full_x_motion) ./ 10;
alpha = 0.6;  % Opacity for the top image

alpha_map = zeros(size(processed_combined_image));
alpha_map(vessel_row_start:vessel_row_end,:) = alpha;

base_norm = mat2gray(base_image);
color_norm  = mat2gray(color_image);

% 2. Map each to their respective colormaps
% 256 is default size of MATLAB colormaps
gray_map = gray(256);
color_map = hot(256);

base_mapped = ind2rgb(round(base_norm * 255) + 1, gray_map);
top_mapped  = ind2rgb(round(color_norm  * 255) + 1, color_map);

% 3. Blend the two using per-pixel alpha
% alpha_map = mat2gray(alpha_map); % Ensure alpha between 0 and 1
composite_rgb = (1 - alpha_map) .* base_mapped + alpha_map .* top_mapped;

fg = figure();
imagesc(x_range.*1000,z_range.*1000,composite_rgb);
fg.Colormap = color_map;

cb = colorbar('Position', [0.92 0.15 0.02 0.7]);
cb.Ticks = linspace(0,1,7);
cb.TickLabels = linspace(0,1.2,7);

title("1.2 cm/s Lateral Flow", "FontSize",14)
ylabel("Depth (mm)");
xlabel("Lateral Position (mm)");

combined_images(:,:,:,forces_frame_id) = composite_rgb; 

%%

figure();

sgtitle('Filtered READI Flow','FontSize',16);

% Loop through the images and display each in the 2x8 grid
for i = 1:16


    subplot(4, 4, i); 
    colormap("gray");
    imagesc(x_range*1000, z_range*1000,(squeeze(processed_frame_array(:,:,frame_range(i)))));
    title(['Image ' num2str(i)]);
    xlabel("mm");
    ylabel("mm");

    % clim([-dynamic_range, 0]); 


end
% colorbar('Position', [0.92 0.15 0.02 0.7]);
