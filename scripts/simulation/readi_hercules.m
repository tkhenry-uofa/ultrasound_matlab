%% General Configuration
clear;
clc;
close("all")

addpath(genpath(pwd));

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

% %%% QL#02
% array_name = "ql02";
% f0 = 3.85e6;  % Transducer center frequency (Hz)
% pitch = 2.5e-4;
% width = 2.2e-4;

%%% M32-1
array_name = 'm32-1';
f0 = 6.25e6;
pitch = 2.5e-4;
width = 2.2e-4;

lambda = c/f0; % Wavelength (m)
kerf = pitch - width;
height = width; % Height of the elements (m)
row_count = 128; % Total from two side by side arrays
column_count = row_count;

set_field('fs', fs);
set_field('c',c);

%% Tx waveform generation
cyc =2 ;
tx_t = 0:1/fs:cyc/f0;

impulse_response=sin(2*pi*f0*(tx_t));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';

excitation=sin(2*pi*f0*(tx_t));

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

% plot(emission)

returned_wave = conv(emission, impulse_response);
returned_time = (1:length(returned_wave))/fs;

% Delay to add to our calculated samples to hit the middle of the return
pulse_delay = returned_time( round(length(returned_time)/2));

%% Tx delays and apodizations

transmit_type = TransmitType.Plane;
sequence_type = SequenceType.Hercules;

src_loc = [0 0 -30]/1000; % m
% src_loc = [0 0 inf]/1000;

prf = c / (113.6/1000 * 2);

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

[tx_array, tx_config] = create_tx_array(tx_config, 4);
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
[points, amps] = point_grid( [0]/1000, [0]/1000, [60]/1000);

point_cell = cell(1,128);

speed = ([0,0,0]/1000)/128; % m/transmit

true_speed = norm(speed) * prf; % m/s

for i=1:tx_config.no_transmits
    point_cell{i} = points + speed * (i-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aquisition Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Staggered transmits
readi_group_count = 1;
readi_group_size = no_transmits / readi_group_count;

staggered_order = zeros(no_transmits,1);
for i = 0:(no_transmits-1)
    
    group_number = mod(i,readi_group_count);
    pos_in_group = floor(i/readi_group_count);

    staggered_order(group_number * readi_group_size + pos_in_group + 1) = i + 1;
end

parallel = true;

% Put results into cells to play nicely with parfor
all_scans = cell(tx_config.no_transmits,1);
rx_element_locs = cell(tx_config.no_transmits,1);
signal_lengths = cell(tx_config.no_transmits,1);

tic;
parfor T = 1:tx_config.no_transmits
    fprintf("Transmission: %d \n", T)
    [data, rx_element_locs{T}, signal_lengths{T}] = aquisition_simulation(tx_config,ones(row_count),point_cell{T},amps,parallel,staggered_order(T));
    % [data, rx_element_locs{T}, signal_lengths{T}] = aquisition_simulation(tx_config,ones(row_count),point_cell{T},amps,parallel,T);
    
    all_scans{T} = data * 1e23;
end


fprintf('Simulation complete, elapsed time: %.2f seconds\n\n', toc);


%% Ensuring each transmit is the same length

signal_lengths = cell2mat(signal_lengths);

max_length = max(signal_lengths);
data_array = single(zeros(max_length,row_count,no_transmits));

data_cell = cell(1,no_transmits);

% Pad each transmit with 0s so they are the same length.
for T = 1:no_transmits
    padded_data = zeros(max_length, column_count);
    padded_data(1:signal_lengths(T), :) = all_scans{T};
    
    data_array(:,:,T) = padded_data;
    data_cell{T} = single(padded_data);
end

%% Export data
% filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\simulations\data\readi\psf_17_cmps.mat";
% overwrite = false;
% export_scan_data(tx_config,data_cell,filepath,overwrite)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volume configuration
x_range = [-50, 50]/1000; 
y_range = [-10, 10]/1000;
z_range = [45, 65]/1000;

resolution = 0.0001; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

%% Readi Beamforming
tic;

cuda_beamform = true;

f_number = 0.0;

% bp.das_shader_id = 0; % Forces
% bp.das_shader_id = 2; % Hercules

sequence_type = 2;

if cuda_beamform == true
    low_res_array = cuda_processing_f2(data_array,tx_config,readi_group_count,vol_config, f_number, sequence_type);
else
    low_res_array = readi_hercules_beamform(tx_config, vol_config, data_array, rx_element_locs,...
        vol_config.x_mid, 1:vol_config.y_count, 1:vol_config.z_count,readi_group_count, f_number);
end

fprintf('Readi beamforming complete, elapsed time: %.2f seconds\n\n', toc);


%% Normal hercules decoding

decoded_data = reshape(data_array, [], no_transmits);
H = hadamard(no_transmits);

% We need a matrix with size (data length*channel count) X transmit count

decoded_data = decoded_data * H;
decoded_data = reshape(decoded_data, max_length, column_count, no_transmits);
decoded_data = hilbert(decoded_data);

% Split into cells to work nicely with parfor
data_cell = mat2cell(decoded_data, max_length, column_count, ones(1, no_transmits));

%% Standard volume building
tic;
fprintf("Volume Building\n")

% all_scans = lpf_rf_data(all_scans,f0,fs);
apodize = true;

unstaggered_data = zeros(size(data_array));

% Put it back in normal order for hercules
for i = 1:length(staggered_order)
    unstaggered_data(:,:,staggered_order(i)) = data_array(:,:,i);
end

hercules_volume_raw = single(zeros(size(low_res_array{1})));

if cuda_beamform == true
    hercules_cell = cuda_processing_f2(unstaggered_data,tx_config,1,vol_config, f_number, sequence_type);
    hercules_volume_raw = hercules_cell{1};
else
    hercules_volume_raw = beamform_volume(tx_config, vol_config, data_cell, rx_element_locs,...
        1:vol_config.x_count, vol_config.y_mid, 1:vol_config.z_count,parallel, apodize);
end



fprintf('Elapsed Time: %.2f seconds\n\n', toc);

%% Motion compensation

shifted_image_raw = readi_peak_motion_comp(low_res_array,readi_group_count,speed,resolution);


%% Processing

readi_h = hadamard(readi_group_count);

high_res_image = real(zeros(size(low_res_array{1})));
for i = 1:readi_group_count
    for j = 1
        high_res_image = high_res_image + low_res_array{i} .* readi_h(j,i);
    end
end

dynamic_range = 50;

processed_readi = process_volume(high_res_image,dynamic_range);
processed_hercules = process_volume(hercules_volume_raw,dynamic_range);
processed_shifted = process_volume(shifted_image_raw,dynamic_range);

% volumeViewer(processed_readi)
% volumeViewer(processed_shifted)

readi_slice = squeeze(processed_readi(:,vol_config.y_mid,:)).';
hercules_slice = squeeze(processed_hercules(:,vol_config.y_mid,:)).';
shifted_slice = squeeze(processed_shifted(:,vol_config.y_mid,:)).';

% readi_slice = squeeze(processed_readi(vol_config.x_mid,:,:)).';
% hercules_slice = squeeze(processed_hercules(vol_config.x_mid,:,:)).';
% shifted_slice = squeeze(processed_shifted(vol_config.x_mid,:,:)).';

processed_low_res = low_res_array;
y_low_res_images = low_res_array;
for i = 1:readi_group_count
    low_res_volume = process_volume(low_res_array{i}, dynamic_range); 

    % processed_low_res{i} = squeeze(low_res_image(:,vol_config.y_mid,:)).';
    processed_low_res{i} = low_res_volume;
    y_low_res_images{i} = squeeze(low_res_volume(vol_config.x_mid,:,:)).';
    x_low_res_images{i} = squeeze(low_res_volume(:,vol_config.y_mid,:)).';
end

% volumeViewer(processed_low_res{1})
% volumeViewer(processed_low_res{2})

%% Plotting

speed_str = sprintf("%0.2f m/s",true_speed);
% figure()

% sgtitle(speed_str + " Point Spread Function", 'FontSize',16)
% sgtitle("Stationary Point Spread Function", 'FontSize',16)

tx_edge = rx_element_locs{1}(1,128)*1000;

% subplot 211
% colormap("gray");
% imagesc(vol_config.x_range*1000, vol_config.z_range*1000, hercules_slice);
% subtitle('Hercules Y Axis', 'FontSize',14);
% axis("image")
% clim([-dynamic_range, 0]); 
% % xline(-tx_edge,'Color','red');
% % xline(tx_edge,'Color','red');
% % legend('Array edge');
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");
figure();
% subplot 212
colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, readi_slice);
subtitle('Readi Stationary Y Axis', 'FontSize',14);
axis("image")
clim([-dynamic_range, 0]); 
% xline(-tx_edge,'Color','red');
% xline(tx_edge,'Color','red');
% legend('Array edge');
xlabel("Lateral Distance (mm)");
ylabel("Axial Distance (mm)");

% figure();
% % subplot 313
% colormap("gray");
% imagesc(vol_config.x_range*1000, vol_config.z_range*1000, shifted_slice);
% subtitle('Readi Motion Compensated Y Axis', 'FontSize',14)
% axis("image")
% clim([-dynamic_range, 0]); 
% % xline(-tx_edge,'Color','red');
% % xline(tx_edge,'Color','red');
% % legend('Array edge');
% xlabel("Lateral Distance (mm)");
% ylabel("Axial Distance (mm)");


colorbar('Position', [0.92 0.15 0.02 0.7]);
% 
% 
% figure();
% 
% % Main title
% sgtitle(speed_str + ' Readi Low Resolution Images Y Axis','FontSize',16);
% % sgtitle('Stationary Readi Low Resolution Images','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(2, 4, i); 
%     colormap("gray");
%     imagesc(vol_config.y_range*1000, vol_config.z_range*1000,(y_low_res_images{i}));
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
%     axis('image')
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end
% colorbar('Position', [0.92 0.15 0.02 0.7]);
