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

returned_wave = conv(emission, impulse_response);
returned_time = (1:length(returned_wave))/fs;

% Delay to add to our calculated samples to hit the middle of the return
pulse_delay = returned_time( round(length(returned_time)/2));

% tx_apro = [zeros(column_count/4,1); hann(column_count/2); zeros(column_count/4,1)]*[zeros(row_count/4,1); hann(row_count/2); zeros(row_count/4,1)]';

%% Tx delays and aprodizations
% tx_apro = hamming(column_count)*ones(1,row_count);
% tx_apro = ones(column_count,1) * hamming(row_count)';
% tx_apro = hamming(column_count)*hamming(row_count)';
tx_apro = ones(column_count, row_count);

% 'plane', 'xLine', 'yLine'
% MUST USE SINGLE QUOTES OR THIS CANT BE PARSED IN C/C++
transmit_type = 'plane';
src_loc = [0 0 -25]/1000; % m
curve_radius = -src_loc(3);

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
    'no_transmits',0,...
    'cross_offset',0,...
    'curve_radius', curve_radius,...
    'sub_element_count',1,....
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

% See the array with apro or delays (slow)
% show_xdc_mod(tx_array,"delays")

%% Rx sequence generation
tx_config.no_transmits = 16;
% If a line crosses two sub arrays offset it by this many elements
% (Helps simulate higher spatial sampling and remove grating lobes)
tx_config.cross_offset = 0; 

% active_rows = 2:(row_count/tx_config.no_transmits):row_count;
% rx_sequence = single_array_diag( tx_config.no_transmits, row_count, active_rows );

rx_sequence = tobe_full_sample(row_count);

% print_sequence(rx_sequence)



%% Point scatter generation
% [points, amps] = point_grid( [0]/1000, [0]/1000, [80]/1000);
[points, amps] = point_grid( [-0]/1000, [0]/1000, [50:10:150]/1000);
% [points, amps] = point_grid(  (-16:5:16)/1000, (-200:5:200)/1000, (5:5:200)/1000);
% [points, amps] = cyst_pht(120000);
% [points, amps] = cyst2x2([-10, 0, 10]/1000, [-0.5, 0.5]/1000,[70, 80, 90]/1000,15000,3/1000);

% [points, amps] = single_cyst(50000, [0, 20]/1000,[-5, 5]/1000,[80, 90]/1000, 2/1000);


%% Aquisition Simulation
parallel = true;
all_scans = cell(tx_config.no_transmits,1);
rx_element_locs = cell(tx_config.no_transmits,1);
rx_starts = cell(tx_config.no_transmits,1);

tic;


parfor T = 1:tx_config.no_transmits
    fprintf("Transmission: %d \n", T)
    config = tx_config;
    if parallel == true
        field_init(-1)
        set_field('fs', config.fs);
        set_field('c',config.c);
    end

    

    element_count = config.rows * config.cols;
    subs = config.sub_element_count;

    tx_array = create_tx_array(config);
    array_details = xdc_get(tx_array);

    physical_pos = array_details(24:26,1:subs:end);

    active_elements = 1:element_count;
    active_elements = active_elements(logical(reshape(rx_sequence{T},[],1)));

    physical_pos = physical_pos(:,active_elements);

    rect_get = array_details(:, ismember(array_details(1,:) + 1,active_elements) );

    rectangles = zeros(19, length(active_elements)*subs);

    rectangles(1, :) = rect_get(1,:) + 1 - (length(active_elements)*(T-1)) ; % 1 indexed physical element
    rectangles(2:13, :) = rect_get(11:22, :); % 4 corners coordinates 
    rectangles(14, :) = 1; % Apodize
    rectangles(15, :) = rect_get(3,:); % Width
    rectangles(16, :) = rect_get(4,:); % Height
    rectangles(17:19, :) = rect_get(8:10, :); % Center


    
    rx_array = xdc_rectangles(rectangles.',physical_pos.',[0,0,10000000000]);

    
    % show_xdc_mod(rx_array,'apro')

    xdc_impulse(rx_array,config.imp);

    [scans, rx_start] = calc_scat_multi( tx_array, rx_array, points, amps );

    rx_starts{T} = rx_start;
    % Add zeros for the time before the first reception
    zero_padding = zeros(length(0:1/config.fs:rx_start), length(physical_pos));

    scan_lines = single([zero_padding; scans]);
    
    rx_element_locs{T} = physical_pos;
    all_scans{T} = hilbert(scan_lines) .* 1e21;
end

for T = 1:tx_config.no_transmits

end

fprintf('Elapsed Time: %.2f seconds\n\n', toc);

%% 
filepath = "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\vrs_transfers\vrs_data\curved\psf.mat";
overwrite = false;
export_scan_data(tx_config,all_scans,filepath,overwrite)


%% Volume configuration
x_range = [-50, 50]/1000;
y_range = [-50, 50]/1000;
z_range = [40, 150]/1000;

resolution = 0.0003; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);


%% Volume building
tic;
fprintf("Volume Building\n")

% all_scans = lpf_rf_data(all_scans,f0,fs);
aprodize = true;



%%


%% Processing
dynamic_range = 50;

volume = beamform_volume(tx_config, vol_config, all_scans, rx_element_locs,...
    vol_config.x_mid, 1:vol_config.y_count, 1:vol_config.z_count,parallel, aprodize);
fprintf('Elapsed Time: %.2f seconds\n\n', toc);
processed_volume = process_volume(volume,vol_config,dynamic_range);


%% Plotting

x_slice = squeeze(processed_volume(vol_config.x_mid,:,:))';

figure(1)

colormap("gray");
imagesc(vol_config.y_range*1000, vol_config.z_range*1000, x_slice);
title('y')
axis("image")
clim([-dynamic_range, 0]); colorbar;

%% Processing


apro_volume = beamform_volume(tx_config, vol_config, all_scans, rx_element_locs,...
    1:vol_config.x_count, vol_config.y_mid, 1:vol_config.z_count,parallel, aprodize);
fprintf('Elapsed Time: %.2f seconds\n\n', toc);


processed_apro = process_volume(apro_volume,vol_config,dynamic_range);
y_slice = squeeze(processed_apro(:,vol_config.y_mid,:))';

figure(2)
colormap("gray");
imagesc(vol_config.x_range*1000, vol_config.z_range*1000, y_slice);
title('x')
axis("image")
clim([-dynamic_range, 0]); colorbar;



% 
% subplot 122
% 
% colormap("gray");
% imagesc(vol_config.x_range*1000, vol_config.z_range*1000, apro_slice);
% axis("image")
% title('With rx apro')
% clim([-dynamic_range, 0]); colorbar;
% sgtitle('128x128 mixes')
% 

