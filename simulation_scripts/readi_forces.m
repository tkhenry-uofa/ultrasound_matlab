%% General Configuration
clear;
clc;
% close("all")

addpath(genpath(pwd));

path(path, "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\field_ii\Field_II_ver_3_30_windows")

field_init(-1)

fs = 50e6; % Sampling frequency (Hz)

% f0 = 2.5e6;
c = 1454;   % Speed of sound (m/s)

%%% A#02
array_name = 'a02';
f0 = 7.8e6;  % Transducer center frequency (Hz)
lambda = c/f0; % Wavelength (m)
pitch = 1.5e-4;
width = 1.2e-4;

% %%% QL#02
% array_name = "ql02";
% f0 = 3.85e6;  % Transducer center frequency (Hz)
% pitch = 2.5e-4;
% width = 2.2e-4;

lambda = c/f0; % Wavelength (m)
kerf = pitch - width;
height = width; % Height of the elements (m)
row_count = 128; % Total from two side by side arrays
column_count = row_count;
samples_per_meter = fs/c/2;

set_field('fs', fs);
set_field('c',c);

%% Tx waveform generation

impulse_t = 0:1/fs:2/f0;
impulse_response=sin(2*pi*f0*(impulse_t));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

cyc =8 ;
tx_t = 0:1/fs:cyc/f0;
excitation=sin(2*pi*f0*(tx_t));

f1 = 6e6;
f2 = 10e6;
B = f2 - f1;                       % Bandwidth [MHz]
T = 30e-6;                                  % exc time [sec]
tapering = 0.20;                        % amplitude tappering [%]
% Create Chirp excitation
t = 0:1/fs:T-1/fs;
F = f1 + B/(2*T) * t;
exc = sin(2*pi*F.*t);
% excitation to create match filter later
excitation = exc.*tukeywin(length(exc),tapering)';

% figure();plot(abs(fft(exc)));

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

% plot(emission)

returned_wave = conv(emission, impulse_response);
returned_time = length(returned_wave)/fs;

match_filter = single(fliplr(emission));

pulse_delay = length(emission)/(2*fs);

% Delay to add to our calculated samples to hit the middle of the return
filter_pulse_delay = 0;%pulse_delay + length(match_filter) * 1 / fs;

%% Tx delays and apodizations

transmit_type = TransmitType.Elevation_Focus;
sequence_type = SequenceType.Readi;

src_loc = [0 0 113.6]/1000; % m
forces_sources = 1:column_count;
discard_transmits = [];

prf = c / (src_loc(3) * 2);

% tx_apo = hamming(column_count)*ones(1,row_count);
% tx_apo = ones(column_count,1) * hamming(row_count)';
tx_apo = hamming(column_count)*hamming(row_count)';
% tx_apo = ones(column_count, row_count);
no_transmits = 128;

tx_config = struct(...
    'f0', f0,...
    'fs',fs,...
    'cols',column_count,...
    'rows',row_count,...
    'width',width,...
    'kerf',kerf,...
    'pitch',pitch,...
    'imp',impulse_response,...
    'ex',excitation,...
    'apo',tx_apo,...
    'transmit',transmit_type,...
    'sequence',sequence_type,...
    'src',src_loc,...
    'forces_srcs', forces_sources,...
    'discard_txs', discard_transmits,...
    'c',c,...
    'pulse_delay', pulse_delay, ...
    'l_angle',0,...
    'r_angle',0,...
    'x_min',0,...
    'x_max',0,...
    'y_min',0,...
    'y_max',0,...
    'no_transmits',no_transmits,...
    'cross_offset',0,...
    'curve_radius',0,...
    'print',true);

[tx_array, tx_config] = create_tx_array(tx_config, 2);
% These angles describe the region hit by the diverging wave
if transmit_type == TransmitType.Ydiv
    tx_config.l_angle = atan( abs(src_loc(2) - tx_config.y_min)/abs(src_loc(3)) );
    tx_config.r_angle = atan( abs(src_loc(2) - tx_config.y_max)/abs(src_loc(3)) );
elseif transmit_type == TransmitType.Xdiv
    tx_config.l_angle = atan( abs(src_loc(1) - tx_config.x_min)/abs(src_loc(3)) );
    tx_config.r_angle = atan( abs(src_loc(1) - tx_config.x_max)/abs(src_loc(3)) );
end

% show_xdc_mod(tx_array,'delay');
tx_config.print = false;

% See the array with apo or delays (slow)
% show_xdc_mod(tx_array,"apo")

%% Point scatter generation
% [points, amps] = point_grid( [-20:5:20]/1000, [0]/1000, [5:5:100]/1000);

% [points, amps] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[70, 80, 90]/1000,15000,2.5/1000);
% [points, amps] = single_cyst(150000, [-25, 25]/1000,[-1, 1]/1000,[60, 90]/1000, 2/1000);

% scatter3(points(:,1),points(:,2),points(:,3), 10, 'filled')

% points = point_grid( [-2.5]/1000, [0]/1000, [50]/1000);

point1 = [0, 0, 80]/1000;
% point2 = [5, 0, 47.5]/1000;
% point3 = [0, 0, 45]/1000;
% point4 = [-1.8, 0, 53.2]/1000;

% amps = [1;1;1;1];
amps = 1;

point_cell = cell(1,128);

speed = (0/1000)/128; % m/transmit
true_speed = speed * prf; % m/s

for i=1:tx_config.no_transmits
    point_cell{i} = [point1 + [speed, 0, 0] * (i-1)]; 
                     % point2 + [0, 0, speed] * (i-1);
                     % point3 + [0, 0, 0] * (i-1);
                     % point4 + [speed/sqrt(2), 0, speed/sqrt(2)] * (i-1)]; 
end

%% Aquisition Simulation
parallel = true;

% Put results into cells to play nicely with parfor
all_scans = cell(tx_config.no_transmits,1);
rx_element_locs = cell(tx_config.no_transmits,1);
signal_lengths = cell(tx_config.no_transmits,1);
max_values = cell(tx_config.no_transmits,1);
tic;
for T = 1:tx_config.no_transmits
    fprintf("Transmission: %d \n", T)
    [all_scans{T}, rx_element_locs{T}, signal_lengths{T}] = aquisition_simulation(tx_config,ones(row_count),point_cell{T},amps,parallel,T);
    max_values{T} = max(all_scans{T},[],"all");
end
fprintf('Simulation complete, elapsed time: %.2f seconds\n\n', toc);

max_value = max(cell2mat(max_values),[],"all");

for T = 1:tx_config.no_transmits

    all_scans{T} = all_scans{T} ./ max_value;

end


%% Processing

signal_lengths = cell2mat(signal_lengths);

max_length = max(signal_lengths)+200;
data_array = zeros(max_length,row_count,no_transmits);

% Pad each transmit with 0s so they are the same length.
for T = 1:no_transmits
    padded_data = zeros(max_length, column_count);
    padded_data(1:signal_lengths(T), :) = all_scans{T};
    
    data_array(:,:,T) = padded_data;
end

% filtered_data = match_filter_data(data_array, match_filter);
filtered_data = data_array;

%%
% data_start = 5300;
% 
% figure(); 
% subplot 311
% plot(emission);
% title("Emission");
% 
% subplot 312
% plot(data_array(data_start:end,64,1));
% title("Unfiltered PSF")
% 
% subplot 313
% plot(filtered_data(data_start:end,64,1));
% title("Filtered PSF");



%% Volume configuration
% x_range = [-10, 10]/1000;
% y_range = [0, 0]/1000;
% z_range = [40, 60]/1000;

x_range = [-10, 10]/1000;
y_range = [0, 0]/1000;
z_range = [60, 100]/1000;


resolution = 0.0002; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

%% Readi Beamforming
tic;
f_number = 1.0;
readi_group_count = 1;
readi_group_size = row_count/readi_group_count;

% bp.das_shader_id = 0; % Forces
% bp.das_shader_id = 2; % Hercules
sequence_type = 0;
cuda_beamform = true;

if cuda_beamform == true
    low_res_array = cuda_processing_f2(data_array,tx_config,readi_group_count,vol_config, f_number, sequence_type);
else
    low_res_array = readi_forces_beamform(tx_config, vol_config, data_array, rx_element_locs,...
    1:vol_config.x_count, 1:vol_config.z_count, readi_group_count, f_number);
end

fprintf('Readi beamforming complete, elapsed time: %.2f seconds\n\n', toc);

% Forces images are already 2D, but X is the leading dim (across rows) 
% So rotate it so that X is acreoss the cols 
for i = 1:readi_group_count
    low_res_array{i} = low_res_array{i}.';
end

rows = vol_config.z_count;
cols = vol_config.x_count;

readi_image_raw = complex(zeros(size(low_res_array{1})));
for i = 1:readi_group_count
    readi_image_raw = readi_image_raw + low_res_array{i};
end

%% Forces Beamforming
tic;
fprintf("Volume Building\n")

tx_config.pulse_delay = filter_pulse_delay;

% all_scans = lpf_rf_data(all_scans,f0,fs);
cuda_beamform = true;
if cuda_beamform == true
    cuda_cell = cuda_processing_f2(filtered_data,tx_config,1,vol_config, f_number, sequence_type);
    forces_image_raw = cuda_cell{1}.';
else
    % Forces Decoding
    data_array_force = reshape(data_array, [], no_transmits);
    H = hadamard(no_transmits);

    We need a matrix with size (data length*channel count) X transmit count

    data_array_force = data_array_force * H;
    data_array_force = reshape(data_array_force, max_length, column_count, no_transmits);
    data_array_force = hilbert(data_array_force);

    Split into cells to work nicely with parfor
    data_cell = mat2cell(data_array_force, max_length, column_count, ones(1, no_transmits));

    forces_image_raw = beamform_volume(tx_config, vol_config, data_cell, rx_element_locs,...
    1:vol_config.x_count, vol_config.y_mid, 1:vol_config.z_count,parallel, f_number);
end

% forces_image_raw = squeeze(forces_image_raw(:,vol_config.y_mid,:)).';



%% Motion compensation

shifted_images = low_res_array;
% ncc_block_matching;

% phase_match_test;
fprintf("motion compensation done\n");



%% Processing
dynamic_range = 60;

shifted_image_raw = complex(zeros(size(shifted_images{1})));
processed_low_res = low_res_array;
processed_shifted_array = low_res_array;
for i = 1:readi_group_count
    processed_low_res{i} = process_volume(low_res_array{i}, dynamic_range); 
    processed_shifted_array{i} = process_volume(shifted_images{i}, dynamic_range); 
    shifted_image_raw = shifted_image_raw + shifted_images{i};
end

readi_image = process_volume(readi_image_raw,dynamic_range);
forces_image = process_volume(forces_image_raw,dynamic_range);
% processed_shifted_image = process_volume(shifted_image_raw,dynamic_range);

% processed_shifted_image = imgaussfilt(processed_shifted_image,2);




%% Plotting

% speed_str = sprintf("%0.2f m/s",true_speed);
figure()

% sgtitle(speed_str + " Point Spread Function", 'FontSize',16)
sgtitle("Convolution, No Pulse Delay", 'FontSize',16)

tx_edge = rx_element_locs{1}(1,128)*1000;


subplot 121
colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, readi_image);
subtitle('No Filter', 'FontSize',14);
axis("image")
clim([-dynamic_range, 0]); 
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");

subplot 122
% colormap("gray");
% imagesc(vol_config.x_range*1000, vol_config.z_range*1000, processed_shifted_image);
% % subtitle('Shifted', 'FontSize',14)
% axis("image")
% clim([-dynamic_range, 0]); 
% % xlabel("Lateral Distance (mm)");
% % ylabel("Axial Distance (mm)");

% subplot 313
colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, forces_image);
subtitle('Filtered', 'FontSize',14);
axis("image")
clim([-dynamic_range, 0]); 
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");



colorbar('Position', [0.92 0.15 0.02 0.7]);



% 
% figure();
% 
% % Main title
% % sgtitle(speed_str + ' Uncompensated Low Resolution Images','FontSize',16);
% sgtitle('Uncompensated','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(4, 4, i); 
%     colormap("gray");
%     imagesc(vol_config.x_range*1000, vol_config.z_range*1000,(processed_low_res{i}));
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end
% % colorbar('Position', [0.92 0.15 0.02 0.7]);
% 
% figure();
% 
% % Main title
% % sgtitle(speed_str + ' Compensated Low Resolution Images','FontSize',16);
% sgtitle('Phase Match Warped','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(4, 4, i); 
%     colormap("gray");
%     imagesc(vol_config.x_range*1000, vol_config.z_range*1000,(processed_shifted_array{i}));
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end
% colorbar('Position', [0.92 0.15 0.02 0.7]);



%% Export data
% filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\real_data\field_ii\readi\cyst_sim_move.mat";
% overwrite = false;
% export_scan_data(tx_config,all_scans,filepath,overwrite)




