%% General Configuration
clear;
clc;
close("all")

path(path, "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\Field_II_ver_3_30_windows")
field_init(-1)

fs = 50e6; % Sampling frequency (Hz)
f0 = 2.5e6;  % Transducer center frequency (Hz)

c = 1540;   % Speed of sound (m/s)
lambda = c/f0; % Wavelength (m)
x_pitch = lambda/2;
y_pitch = lambda/2;
kerf = 28e-6; % Kerf (m)
width = x_pitch - kerf; % Width of the element (m)
height = y_pitch - kerf; % Height of the elements (m)
row_count = 48; % Total from two side by side arrays
column_count = 64;



x_subs = 2;
y_subs = 4;

sub_array_count = x_subs*y_subs;

row_width = 2;
row_length = column_count/x_subs;


set_field('fs', fs);
set_field('c',c);

%% Tx waveform generation
cyc = 2;
tx_t = 0:1/fs:cyc/f0;

imp_t = 0:1/fs:cyc/f0;

impulse_response=sin(2*pi*f0*(imp_t));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

excitation=sin(2*pi*f0*(tx_t));

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

returned_wave = conv(emission, impulse_response);

returned_time = (1:length(returned_wave))/fs;

% Delay to add to our calculated samples to hit the middle of the return
pulse_delay = returned_time( round(length(returned_time)/2));

% tx_apro = [zeros(column_count/4,1); hann(column_count/2); zeros(column_count/4,1)]*[zeros(row_count/4,1); hann(row_count/2); zeros(row_count/4,1)]';

%% Tx delays and aprodizations
tx_apro = hamming(column_count)*ones(1,row_count);
% tx_apro = ones(column_count,1) * hamming(row_count)';
% tx_apro = hamming(column_count)*hamming(row_count)';
% tx_apro = ones(column_count, row_count);

% 'plane', 'xLine', 'yLine'
% MUST USE SINGLE QUOTES OR THIS CANT BE PARSED IN C/C++
transmit_type = 'yLine';
src_loc = [0 0 -40]/1000; % m

tx_config = struct(...
    'f0', f0,...
    'fs',fs,...
    'cols',column_count,...
    'rows',row_count,...
    'width',width,...
    'kerf',kerf,...
    'x_pitch',x_pitch,...
    'y_pitch', y_pitch,...
    'imp',impulse_response,...
    'ex',excitation,...
    'apro',tx_apro,...
    'transmit',transmit_type,...
    'src',src_loc,...
    'c',c,...
    'pulse_delay', pulse_delay, ...
    'l_angle',0,...
    'r_angle',0,...
    'x_min',0,...
    'x_max',0,...
    'y_min',0,...
    'y_max',0,...
    'x_range',0,...
    'y_range',0,...
    'no_transmits',0,...
    'cross_offset',0,...
    'print',true);

[tx_array, tx_config] = create_tx_array(tx_config);
% These angles describe the region hit by the diverging wave
if strcmp(transmit_type, 'xLine')
    tx_config.l_angle = atan( abs(src_loc(2) - tx_config.y_min)/abs(src_loc(3)) );
    tx_config.r_angle = atan( abs(src_loc(2) - tx_config.y_max)/abs(src_loc(3)) );
elseif strcmp(transmit_type, 'yLine')
    tx_config.l_angle = atan( abs(src_loc(1) - tx_config.x_min)/abs(src_loc(3)) );
    tx_config.r_angle = atan( abs(src_loc(1) - tx_config.x_max)/abs(src_loc(3)) );
end

show_xdc_mod(tx_array,'delay');
tx_config.print = false;

% See the array with apro or delays (slow)
% show_xdc_mod(tx_array,"delays")

%% Rx sequence generation
tx_config.no_transmits = 3;
rx_sequence = octopus_sequence();

% tx_config.no_transmits = 12;
% rx_sequence = cardiac_full_sample();

print_sequence(rx_sequence);

%% Point scatter generation
% [points, amps] = point_grid( [0]/1000, [0]/1000, [100, 110, 120]/1000);

% [points, amps] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[100, 110, 120]/1000,15000,3/1000);

% [points2, amps2] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[70, 80, 90]/1000,15000,3/1000);
% [points3, amps3] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[90, 100, 110]/1000,15000,3/1000);
% 
% [points4, amps4] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[110, 120, 130]/1000,15000,3/1000);
% [points5, amps5] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[130, 140, 150]/1000,15000,3/1000);
% 
% points2(:, [1, 2, 3]) = points2(:,[2, 1, 3]);
% points4(:, [1, 2, 3]) = points4(:,[2, 1, 3]);
%
% points = [points1;points2;points3;points4;points5];
% amps = [amps1;amps2;amps3;amps4;amps5];

% points = [points4; points5];
% amps = [amps4; amps5];

[points, amps] = point_grid( (-5:5:40)/1000, (-5:5:40)/1000, (5:5:195)/1000);

%% Aquisition Simulation
parallel = true;
all_scans = cell(tx_config.no_transmits,1);
rx_element_locs = cell(tx_config.no_transmits,1);

tic;
parfor T = 1:tx_config.no_transmits
    fprintf("Transmission: %d \n", T)
    [rx_scans, rx_element_locs_tx] = aquisition_simulation(tx_config,rx_sequence{T},points,amps,parallel);

    all_scans{T} = single(hilbert(rx_scans));
    % all_scans{T} = rx_scans;
    rx_element_locs{T} = single(rx_element_locs_tx);
    % Simulate the electronic summing of elements in the same column

    % scans_temp = zeros( size(rx_scans,1), size(rx_scans,2)/2);
    % elements_temp = zeros(3, size(rx_element_locs_tx,2)/2);
    % for s = 1:sub_array_count
    %     for e = 1:row_length
    %         first_index = (s-1) * row_length * row_width + e;
    %         second_index = first_index + row_length;
    % 
    %         output_index = (s-1) * row_length + e;
    % 
    %         scans_temp(:, output_index) = rx_scans(:, first_index) + rx_scans(:,second_index);
    % 
    %         element1 = rx_element_locs_tx(:,first_index);
    %         element2 = rx_element_locs_tx(:,second_index);
    %         element_loc = element1 + (element2 - element1)/2;
    %         elements_temp(:, output_index) = element_loc;
    % 
    %     end
    % 
    % 
    % end
    % 
    % all_scans{T} = single(hilbert(rx_scans));
    % rx_element_locs{T} = single(elements_temp);
end

% Update the y locations
% tx_config.y_range = tx_config.y_range(1:2:end-1) + tx_config.y_pitch/2;

fprintf('Elapsed Time: %.2f seconds\n\n', toc);

%% 
filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\mixes\data\oct\large_psf_grid_full.mat";
overwrite = false;
export_scan_data(tx_config,rx_element_locs,all_scans,filepath,overwrite)

%% Volume configuration
x_range = [-5, 40]/1000;
y_range = [-1, 1]/1000;
z_range = [0, 200]/1000;

resolution = 0.0001; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

%% Volume building
tic;
fprintf("Volume Building\n")
aprodize = false;

% X slice
volume = beamform_volume(tx_config, vol_config, all_scans, rx_element_locs,...
    1:vol_config.x_count, vol_config.y_mid, 1:vol_config.z_count,parallel,aprodize);

% volume(vol_config.x_mid, vol_config.y_mid, :) = 0; % Don't wanna beamform this twice
% 
% % Y slice
% volume = volume + beamform_volume(tx_config, vol_config, all_scans, rx_element_locs,...
%     vol_config.x_mid, 1:vol_config.y_count, 1:vol_config.z_count,parallel,aprodize);
% fprintf('Elapsed Time: %.2f seconds\n\n', toc);

%% Processing
dynamic_range = 40;
processed_volume = process_volume(volume,vol_config,dynamic_range);

%% Plotting
% figure_num = 1;
x_slice = squeeze(processed_volume(:,vol_config.y_mid,:))';
y_slice = squeeze(processed_volume(vol_config.x_mid,:,:))';
% figure_num = plot_slice(x_slice, figure_num,vol_config.x_range, vol_config.z_range);
% title("X Axis")
% figure_num = plot_slice(y_slice, figure_num,vol_config.y_range, vol_config.z_range);
% title("Y Axis")

%%
figure(1)

% subplot 121

colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, x_slice);
title('X axis')
axis("image")
clim([-dynamic_range, 0]); colorbar;

% Get the nearest line to the array edges
max_i = round(interp1(vol_config.x_range, 1:vol_config.x_count, tx_config.x_max,"linear",vol_config.x_count));
min_i = round(interp1(vol_config.x_range, 1:vol_config.x_count, tx_config.x_min,"linear",1));



% xline(vol_config.x_range(max_i)*1000,Color="yellow");
% xline(vol_config.x_range(min_i)*1000,Color="yellow");
% 
% subplot 122
% 
% colormap("gray");
% imagesc(vol_config.y_range*1000, vol_config.z_range*1000, y_slice);
% axis("image")
% title('Y axis')
% clim([-dynamic_range, 0]); colorbar;
% sgtitle('Octopus focus at 110mm')

% xline(vol_config.x_range(max_i)*1000,Color="yellow");
% xline(vol_config.x_range(min_i)*1000,Color="yellow");



%% Extra

% x_y_compare;


