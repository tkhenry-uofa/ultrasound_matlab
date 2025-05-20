%% General Configuration
clear;
clc;
close("all")

path(path, "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\Field_II_ver_3_30_windows")
field_init(-1)

fs = 50e6; % Sampling frequency (Hz)
% f0 = 2.5e6;
c = 1454;   % Speed of sound (m/s)

% %%% A#02
% array_name = 'a02';
% f0 = 7.8e6;  % Transducer center frequency (Hz)
% lambda = c/f0; % Wavelength (m)
% pitch = 1.5e-4;
% pitch = lambda;
% width = pitch - 3e-5;
% width = 1.2e-4;

%%% QL#02
array_name = "ql02";
f0 = 3.85e6;  % Transducer center frequency (Hz)
pitch = 2.5e-4;
width = 2.2e-4;

lambda = c/f0; % Wavelength (m)
kerf = pitch - width;
height = width; % Height of the elements (m)
row_count = 128; % Total from two side by side arrays
column_count = row_count;

sub_array_size = 64;

set_field('fs', fs);
set_field('c',c);

%% Tx waveform generation
cyc = 2;
tx_t = 0:1/fs:cyc/f0;

impulse_response=sin(2*pi*f0*(tx_t));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

excitation=sin(2*pi*f0*(tx_t));

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

returned_wave = conv(emission, impulse_response);
returned_time = (1:length(returned_wave))/fs;

% Delay to add to our calculated samples to hit the middle of the return
pulse_delay = returned_time( round(length(returned_time)/2));

% tx_apo = [zeros(column_count/4,1); hann(column_count/2); zeros(column_count/4,1)]*[zeros(row_count/4,1); hann(row_count/2); zeros(row_count/4,1)]';

%% Tx delays and apodizations
% tx_apo = hamming(column_count)*ones(1,row_count);
% tx_apo = ones(column_count,1) * hamming(row_count)';
tx_apo = hamming(column_count)*hamming(row_count)';
% tx_apo = ones(column_count, row_count);

%% Tx delays and apodizations

transmit_type = TransmitType.Plane;
sequence_type = SequenceType.Mixes;

src_loc = [0 0 0]/1000; % m

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
    'src',src_loc,...
    'c',c,...
    'pulse_delay', pulse_delay, ...
    'transmit',transmit_type,...
    'sequence',sequence_type,...
    'l_angle',0,...
    'r_angle',0,...
    'x_min',0,...
    'x_max',0,...
    'y_min',0,...
    'y_max',0,...
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

% show_xdc_mod(tx_array,'delay');
tx_config.print = false;

% See the array with apo or delays (slow)
% show_xdc_mod(tx_array,"delays")

%% Rx sequence generation
tx_config.no_transmits = 16;
% If a line crosses two sub arrays offset it by this many elements
% (Helps simulate higher spatial sampling and remove grating lobes)
tx_config.cross_offset = 0; 

active_rows = 1:(sub_array_size/tx_config.no_transmits):(sub_array_size-tx_config.cross_offset);
rx_sequence = diagonal_cross_sequence( tx_config.no_transmits, sub_array_size, active_rows, tx_config.cross_offset );

% rx_sequence = tobe_full_sample(128);

print_sequence(rx_sequence);

% full_rx = zeros(128);
% for i=1:tx_config.no_transmits
%     full_rx = full_rx + rx_sequence{i};
% end
% imagesc(full_rx);
% 
% % Compute 2D FFT
% fftResult = fft2(full_rx);
% 
% % Shift the zero-frequency component to the center
% fftShifted = fftshift(fftResult);
% 
% % Compute the magnitude of the FFT
% fftMagnitude = abs(fftShifted);
% 
% % Optionally take the logarithm for better visibility
% fftLogMagnitude = 20*log(fftMagnitude); % Adding 1 to avoid log(0)
% 
% fftLogMagnitude = max(fftLogMagnitude, -50);
% 
% % Display the FFT result
% figure;
% imagesc(fftMagnitude);
% colormap('jet'); % Optional: Choose a colormap
% colorbar;
% title('8 Transmit Offset 4');


%% Volume configuration
x_range = [-20, 20]/1000;
y_range = [-5, 5]/1000;
z_range = [40, 70]/1000;

resolution = 0.0001; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

%% Point scatter generation
[points, amps] = point_grid( [0]/1000, [0]/1000, [55]/1000);
% [points, amps] = point_grid( [0]/1000, [0]/1000, [10 40 70 100 130 160 190 ]/1000);
% [points, amps] = cyst_pht(120000);
% [points1, amps1] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[70, 80, 90]/1000,15000,3/1000);
% 
% [points2, amps2] = cyst2x2([5, 15, 25]/1000, [-0.5, 0.5]/1000,[90, 100, 110]/1000,15000,3/1000);
% [points3, amps3] = cyst2x2([5, 15, 25]/1000, [-0.5, 0.5]/1000,[90, 100, 110]/1000,15000,3/1000);
% 
% [points4, amps4] = cyst2x2([10, 20, 30]/1000, [-0.5, 0.5]/1000,[110, 120, 130]/1000,15000,3/1000);
% [points5, amps5] = cyst2x2([10, 20, 30]/1000, [-0.5, 0.5]/1000,[110, 120, 130]/1000,15000,3/1000);
% 
% points3(:, [1, 2, 3]) = points3(:,[2, 1, 3]);
% points5(:, [1, 2, 3]) = points5(:,[2, 1, 3]);
% 
% 
% points = [points1;points2;points3;points4;points5];
% amps = [amps1;amps2;amps3;amps4;amps5];

% [points, amps] = point_grid( (-5:5:40)/1000, (-5:5:40)/1000, (5:5:195)/1000);
% [points, amps] = single_cyst(200000, [-5, 5]/1000,[0, 30]/1000,[100, 110]/1000, 3/1000);

% points = [points1;points2;points3];
% amps = [amps1;amps2;amps3];

%% Aquisition Simulation
parallel = true;
all_scans = cell(tx_config.no_transmits,1);

all_scans_normal = cell(tx_config.no_transmits,1);
rx_element_locs = cell(tx_config.no_transmits,1);

if parallel == true
    num_workers = 16;
else
    num_workers = 0;
end
tic;
for T = 1:tx_config.no_transmits
    fprintf("Transmission: %d \n", T)
    [rx_scans, rx_element_locs{T}] = aquisition_simulation(tx_config,rx_sequence{T},points,amps,parallel,T);
    all_scans{T} = hilbert(rx_scans.* 1e28);
    all_scans_normal{T} = rx_scans;
end

fprintf('Elapsed Time: %.2f seconds\n\n', toc);

%% 
filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\mixes\data\64_tobe_grid\cyst_tube_16_offset.mat";
overwrite = false;
% export_scan_data(tx_config,rx_element_locs,all_scans,filepath,overwrite)


%% Volume building
tic;
fprintf("Volume Building\n")

% all_scans = lpf_rf_data(all_scans,f0,fs);
apodize = false;

volume = beamform_volume(tx_config, vol_config, all_scans, rx_element_locs,...
    1:vol_config.x_count, vol_config.y_mid, 1:vol_config.z_count,parallel, apodize);



%% Processing
dynamic_range = 40;
processed_volume = process_volume(volume,dynamic_range);

%% Plotting
% figure_num = 1;
% p_slice = squeeze(processed_volume(:,vol_config.y_mid,:))';
% figure_num = plot_slice(p_slice, figure_num, x_range, z_range);




%%
x_slice = squeeze(processed_volume(:,vol_config.y_mid,:))';
y_slice = squeeze(processed_volume(vol_config.x_mid,:,:))';
figure(1)

colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, x_slice);
title('X axis')
axis("image")
clim([-dynamic_range, 0]); colorbar;



