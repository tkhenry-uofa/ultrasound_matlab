clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path   = "vrs_data/match_filter_test/";

% data_path = data_path + "250514_MN32-5_reso_FORCES-TxColumn/";
% data_file_name = "250514_MN32-5_reso_FORCES-TxColumn_00";
% params_path = data_path + '250514_MN32-5_reso_FORCES-TxColumn.bp';

% data_path = data_path + "250514_MN32-5_reso_FORCES-TxColumn-Chirp-2e-05/";
% data_file_name = "250514_MN32-5_reso_FORCES-TxColumn-Chirp-2e-05_00";
% params_path = data_path + '250514_MN32-5_reso_FORCES-TxColumn-Chirp-2e-05.bp';

% data_path = data_path + "250514_MN32-5_cyst_FORCES-TxColumn/";
% data_file_name = "250514_MN32-5_cyst_FORCES-TxColumn_00";
% params_path = data_path + '250514_MN32-5_cyst_FORCES-TxColumn.bp';
% 
data_path = data_path + "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-2e-05/";
data_file_name = "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-2e-05_00";
params_path = data_path + '250514_MN32-5_cyst_FORCES-TxColumn-Chirp-2e-05.bp';
% 
% data_path = data_path + "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-5e-06/";
% data_file_name = "250514_MN32-5_cyst_FORCES-TxColumn-Chirp-5e-06_00";
% params_path = data_path + '250514_MN32-5_cyst_FORCES-TxColumn-Chirp-5e-06.bp';


chirp_length = 2e-5;
% chirp_length = 5e-6;


load(data_path + "postVsx.mat","Scan");
bandwidth = Scan.Die.Bandwidth;

data_file_name   = data_file_name + ".zst";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

data_file = fopen(data_path + data_file_name, "r");
raw_data = fread(data_file, '*uint8');
fclose(data_file); 

data = ornot_zstd_decompress_mex(raw_data);

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));
data = reshape(data, raw_bp.raw_data_dim(1),raw_bp.raw_data_dim(2));


% data(1:2,:) = 0;
fc = raw_bp.center_frequency;
fs = raw_bp.sampling_frequency;


%% Match Filter setup
% Impulse response, 1 cyc hamming modulated sin
cycles = 1;
impulse_response = sin(2*pi*fc*(0:1/fs:cycles/fc)).';
impulse_response = impulse_response .* hamming(length(impulse_response));


% Sin 8 excitation
% cycles = 8;
% excitation = sin(2*pi*fc*(0:1/fs:cycles/fc)).';

% Chirp excitation
f1 = bandwidth(1); 
f2 = bandwidth(2);  
BW = f2 - f1;                       
tapering = 0.20;         
% Create Chirp excitation
t = 0:1/fs:chirp_length-1/fs;
F = f1 + BW/(2*chirp_length) * t;
excitation = sin(2*pi*F.*t);
excitation = excitation.*tukeywin(length(excitation),tapering)';

match_filter = excitation;
% match_filter = conv(excitation, impulse_response);

match_filter = single(fliplr(conv(match_filter, impulse_response)));
% match_filter = match_filter./sum(abs(match_filter));

% figure();plot(abs(fftshift(fft(match_filter))));

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

sample_count = single(bp.dec_data_dim.x);
rx_channel_count = single(bp.dec_data_dim.y);
transmit_count = single(bp.dec_data_dim.z);

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
bp.readi_group_size = bp.dec_data_dim.z;
bp.readi_group_id = 0;


%%
% 0 = Forces, 1 = Uforces, 2 = Hercules
bp.das_shader_id = 0;
% bp.das_shader_id = 2;

% NOTE: plane and position along plane normal for beamforming 2D HERCULES
bp.beamform_plane = 0;
bp.off_axis_pos = 0;
fc = bp.center_frequency;
fs = bp.sampling_frequency;

bp.f_number = 0.5;


% The end of every channel has padded zeros not from any transmit
expected_length = bp.dec_data_dim.x * bp.dec_data_dim.z;

cropped_data = data(1:expected_length, :);

bp.rf_raw_dim.x = expected_length;

shuffled_size = [bp.dec_data_dim.x, bp.dec_data_dim.z, bp.rf_raw_dim.y];
cropped_data = reshape(cropped_data, shuffled_size);

start_depth = 15/1000;
start_sample = round(2 * start_depth * bp.sampling_frequency / bp.speed_of_sound);
tx_end_sample = start_sample -1;

cropped_data(1:tx_end_sample,:,:) = 0;

transmit_matrix = permute(cropped_data, [1 3 2]);
tx_mat_size = size(transmit_matrix);
transmit_matrix = reshape(transmit_matrix, [], bp.dec_data_dim.z);

decoded_data = int16(single(transmit_matrix) * hadamard(128));

% Keep everything in sample-transmit-channel order for consistency
decoded_data = permute(reshape(decoded_data, tx_mat_size),[1 3 2]);


% shuffled_data = shuffled_data(start_sample:end,:,:);
% 
% bp.dec_data_dim.x = bp.dec_data_dim.x - tx_end_sample;
% bp.rf_raw_dim.x = bp.rf_raw_dim.x - tx_end_sample * 256;

filtered_data = match_filter_data(cropped_data,match_filter,fs);
filtered_data = int16(filtered_data);

% Match filtering can smear noise into the transmit region
filtered_data(1:tx_end_sample,:,:) = 0; 

% filtered_dec_data = match_filter_data(decoded_data,match_filter,fs);
% filtered_dec_data = int16(filtered_dec_data);
% 
% filtered_dec_data(1:tx_end_sample,:,:) = 0; 

%%

% data_start = 1;
% 
% c = raw_bp.speed_of_sound;
% 
% meters_per_sample = c / fs;
% depth_per_sample = meters_per_sample/2;
% max_depth = depth_per_sample * sample_count;
% 
% time_offset_depth = raw_bp.time_offset * c / 2;
% depth_range = linspace(0,max_depth,sample_count) - time_offset_depth;
% 
% 
% center_channel = bp.channel_mapping(81)+1;
% 
% figure(); 
% subplot 311
% plot(match_filter);
% title("Match Filter");
% 
% subplot 312
% plot(depth_range,cropped_data(data_start:end,1,center_channel));
% title("Unfiltered")
% 
% subplot 313
% plot(depth_range,filtered_data(data_start:end,1,center_channel));
% title("Filtered");



% 
% for i=1:128
%     channel = bp.channel_mapping(i) +1;
% 
%     cropped_signal = (cropped_data(data_start:end,1,channel));
%     filtered_signal = (filtered_data(data_start:end,1,channel));
% 
%     % % decoded_signal = 20*log10(decoded_signal);
%     % % decoded_signal = decoded_signal - max(decoded_signal);
%     % % decoded_signal = max(decoded_signal, -60);
%     % % 
%     % % filtered_dec_signal = 20*log10(filtered_dec_signal);
% 
%     subplot 211
%     plot(depth_range,cropped_signal);
%     title_str = sprintf("Channel %i",i);
%     title(title_str);
%     % pause(0.1)
% 
%     subplot 212
%     % channel = bp.channel_mapping(i);
%     plot(depth_range, filtered_signal);
%     title_str = sprintf("Channel %i",i);
%     title(title_str);
%     pause(0.05)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-30, 30]/1000;
y_range = [-20, 20]/1000;
z_range = [5, 80]/1000;

if bp.das_shader_id == 0
    y_range = [0,0]; % Forces
end

bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.0002;
axial_resolution = 0.0002;

bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution );
bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution );
bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution );
bp.output_points.w = 1; % Number of frames for averaging

if bp.output_points.y <= 0
    bp.output_points.y = 1;
end

x_range = linspace(x_range(1), x_range(2), bp.output_points.x);
y_range = linspace(y_range(1), y_range(2), bp.output_points.y);
z_range = linspace(z_range(1), z_range(2), bp.output_points.z);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Beamforming


output_counts_xyz.x = bp.output_points.x;
output_counts_xyz.y = bp.output_points.y;
output_counts_xyz.z = bp.output_points.z;


% Complex volumes aren't supported so they're interleaved
interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];

if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end

loadlibrary('cuda_transfer')

%% Unfiltered

bp.filter_length = 0;
bp.match_filter = single(zeros(1024,1));

calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

fprintf("Beamforming Unfiltered\n");

volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

test = reshape(cropped_data,[],256);

[~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
    pipe_name, smem_name, cropped_data, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
fprintf("Received Response\n")

volume = reshape(volume_ptr,interleaved_volume_size );

real_vol = volume(1:2:end,:,:);
im_vol = volume(2:2:end,:,:);

unfiltered_raw_img = squeeze(complex(real_vol, im_vol));

%% Filtered
fprintf("Beamforming Filtered\n");

bp.filter_length = length(match_filter);
bp.match_filter(1:length(match_filter)) = match_filter;

bp.time_offset = bp.time_offset + (length(match_filter)-1) * 1.0 / bp.sampling_frequency;

calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));

[~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
    pipe_name, smem_name, cropped_data, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
fprintf("Received Response\n")

volume = reshape(volume_ptr,interleaved_volume_size );

real_vol = volume(1:2:end,:,:);
im_vol = volume(2:2:end,:,:);

filtered_raw_img = squeeze(complex(real_vol, im_vol));

   

unloadlibrary('cuda_transfer')


%% Post processing 
volume_size = interleaved_volume_size;
volume_size(1) = interleaved_volume_size(1)/2;

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


% %
% figure();
% 
% % Main title
% % sgtitle(speed_str + ' Readi Low Resolution Images','FontSize',16);
% sgtitle('Stationary Readi Low Resolution Images','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(4, 4, i); 
%     colormap("gray");
%     imagesc((low_res_images{i}).');
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
%     axis image;
% 
%     % clim([-dynamic_range, 0]); 
% 
% 
% end

