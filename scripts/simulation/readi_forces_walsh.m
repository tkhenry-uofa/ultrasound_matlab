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

cyc = 2;
tx_t = 0:1/fs:cyc/f0;
excitation=sin(2*pi*f0*(tx_t));

% figure();plot(abs(fft(exc)));

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

% plot(emission)

returned_wave = conv(emission, impulse_response);
returned_time = length(returned_wave)/fs;

pulse_delay = length(returned_wave)/(2*fs);

%% Tx delays and apodizations

transmit_type = TransmitType.Elevation_Focus;
sequence_type = SequenceType.FORCES;

src_loc = [0 0 60]/1000; % m
forces_sources = 1:column_count;
discard_transmits = [];

prf = c / (src_loc(3) * 2);

half_hamming = zeros(1,column_count);
half_hamming(33:96) = hamming(64);

tx_apo = ones(row_count,column_count);
% tx_apo = hamming(column_count)*hamming(row_count).';
% tx_apo = hamming(column_count)*half_hamming;
no_transmits = 128;
walsh_ordering = true;

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
    'walsh_ordering',walsh_ordering);

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

point1 = [0, 0, 55]/1000;


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

%%
% sample = all_scans{1};
% sample = sample(:,64);
% figure();
% plot(sample(3750:end));
% title("10 mm rx aperture")


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


%% Volume configuration
x_range = [-10, 10]/1000;
y_range = [0, 0]/1000;
z_range = [45, 65]/1000;


resolution = 0.0001; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

vol_config.f_number = 1;

%% Readi Beamforming
tic;
readi_group_count = 4;

vol_config.readi_group_count = readi_group_count;
readi_group_size = row_count/readi_group_count;

% bp.das_shader_id = 0; % Forces
% bp.das_shader_id = 2; % Hercules
sequence_type = 0;
cuda_beamform = true;

if cuda_beamform == true
    low_res_array = cuda_processing_f2(data_array,tx_config,vol_config, true);
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




%% Processing
dynamic_range = 50;

readi_image_raw = complex(zeros(size(low_res_array{1})));
processed_low_res = low_res_array;
for i = 1:readi_group_count
    processed_low_res{i} = process_volume(low_res_array{i}, dynamic_range); 
    readi_image_raw = readi_image_raw + low_res_array{i};
end

readi_image = process_volume(readi_image_raw,dynamic_range);


%% Plotting

figure()
plot_bmode(readi_image, x_range, z_range, "Walsh FORCES");
colorbar;

plot_image_grid(processed_low_res, [2 2], x_range, z_range, "Walsh Low Res")

%% Export data
% filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\real_data\field_ii\readi\cyst_sim_move.mat";
% overwrite = false;
% export_scan_data(tx_config,all_scans,filepath,overwrite)




