%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script simulates a a 2x2 array of 64x64 TOBE arrays side by side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close("all")

path(path, "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\Field_II_ver_3_30_windows")
field_init(-1)

fs = 50e6; % Sampling frequency (Hz)
f0 = 2.5e6;  % Transducer center frequency (Hz)
c = 1540;   % Speed of sound (m/s)
lambda = c/f0; % Wavelength (m)


row_count = 128; % Total from two side by side arrays
column_count = row_count;
array_size = row_count/2;

pitch = lambda;
kerf_x = 10e-6; % Kerf in x-direction (m)
kerf_y = 10e-6; % Kerf in y-direction
width = pitch - kerf_x; % Width of the element (m)
height = width; % Height of the elements (m)

% Plane = straight plane wave
% xDiv = diverging wave in x direction (line source in y)
% yDiv = diverging wave in y direction
transmit_type = "plane";

src_loc = [0 0 -30]/1000; % m

set_field('fs', fs);
set_field('c',c);

%% Volume configuration
% Volume range (m)
x_min = -15/1000;
x_max = 15/1000; 

y_min = -15/1000;
y_max = 15/1000;

z_min = 35/1000;
z_max = 65/1000;

resolution = 0.00015;

lateral_resolution = resolution; % One cell per 100 microns

lateral_fs = 1/lateral_resolution;

axial_resolution = resolution;

x_range = x_min:lateral_resolution:x_max;
y_range = y_min:lateral_resolution:y_max;
z_range = z_min:axial_resolution:z_max;

x_count = length(x_range);
x_mid = round(x_count/2);
y_count = length(y_range);
y_mid = round(y_count/2);
z_count = length(z_range);
z_mid = round(z_count/2);

volume = single(zeros(x_count, y_count, z_count));

% [-25 -12.5 0 12.5 25]
%% Point Generation
[points, amps] = point_grid( [0]/1000, [0]/1000, [50]/1000);
% [points, amps] = cyst_pht(600000);



%% Transmit array generation
tx_focus = [0 0 100000000000];
tx_apro = ones(row_count, row_count);
tx_array = xdc_2d_array(row_count,row_count,width,height,kerf_x,kerf_y,tx_apro,1,1,tx_focus);
xdc_center_focus(tx_array, [0,0,0]);
xdc_focus(tx_array,0,tx_focus);



% Impulse, a basic cos wave with an exponential envelope 
cyc = 1;
tx_t = 0:1/fs:cyc/f0;

impulse_response=sin(2*pi*f0*(tx_t));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (tx_array, impulse_response);

excitation=sin(2*pi*f0*(tx_t));
xdc_excitation (tx_array, excitation);

emission = conv(impulse_response, excitation);
emission_time = (1:length(emission))/fs;

returned_wave = conv(emission, impulse_response);
returned_time = (1:length(returned_wave))/fs;

pulse_delay = returned_time( round(length(returned_time)/2));

% Apodization
hamming_window_x = hamming(row_count)*ones(1,row_count);
hamming_window_y = ones(row_count,1)*hamming(row_count)';
hamming_window_both = hamming(row_count)*hamming(row_count)';
ele_apodization(tx_array,(1:(row_count)*(row_count))',hamming_window_x(:));

% tx_delays(tx_array,src_loc,transmit_type,row_count,column_count,c,true);


%% Receive array generation
% If true each cross on the 4 arrays is offset from its two neighbours
% so there is never a straight 128 element line
offset_crosses = 2; 
no_transmits = 16;


active_rows = round(linspace(1,array_size-offset_crosses, no_transmits));

active_rows = 2:(array_size/no_transmits):(array_size-offset_crosses);

cross_sequence = diagonal_cross_sequence( no_transmits, array_size, active_rows, offset_crosses );

% print_sequence(cross_sequence)

transmit_weighting = hamming(no_transmits);
rx_apro = hanning(256);
rx_apro = rx_apro(129:256);



%% Aquisition Simulation
all_scans = cell(no_transmits,1);
rx_element_locs= cell(no_transmits,1);
tic;
parfor T = 1:no_transmits

    field_init(-1)
    set_field('fs', fs);
    set_field('c',c);
    fprintf("Transmission: %d \n", T)
    
    tx_array = xdc_2d_array(row_count,row_count,width,height,kerf_x,kerf_y,tx_apro,1,1,tx_focus);
    xdc_center_focus(tx_array, [0,0,0]);
    xdc_focus(tx_array,0,tx_focus)
    % ele_apodization(tx_array,(1:(row_count)*(row_count))',hamming_window_both(:));

    % tx_delays(tx_array,src_loc,transmit_type,row_count,column_count,c,false);


    xdc_impulse (tx_array, impulse_response);
    xdc_excitation (tx_array, excitation);



    rx_cross = cross_sequence{T};
    rx_focus = tx_focus;
    rx_array = xdc_2d_array(row_count,row_count,width,height,kerf_x,kerf_y,ceil(rx_cross),1,1,rx_focus);

    xdc_impulse(rx_array,impulse_response);
    element_locs = xdc_get(rx_array);
    element_locs = element_locs(24:26,:);
    element_count = length(element_locs);

    rx_element_locs{T} = single(element_locs);

    [scans, rx_start] = calc_scat_multi( tx_array, rx_array, points, amps );

    % Add zeros for the time before the first reception
    zero_scans = zeros(length(0:1/fs:rx_start), element_count);

    transformed = single(hilbert([zero_scans; scans]));
    all_scans{T} = transformed;

end




%% Image building
fprintf("Volume Building\n")
loop_average = 0;
variance_volume = single(volume);
max_dist = sqrt(x_max^2 + z_max^2 );
parfor i = 1:x_count
    tic;
    % remaining = loop_average * (x_count + 1 - i );
    % min_remaining = floor(remaining/60);
    % sec_remaining = mod(remaining, 60);
    % percent_remaining = (i-1)*100/x_count;
    % fprintf("Slice %d/%d (%.2f%%)    %dm %.2fs remaining (%.2f s/slice)\n", i, x_count, percent_remaining, min_remaining, sec_remaining, loop_average);
    x_loc = x_range(i);

    for j = y_mid
        % fprintf("Slice: %d \n", j)
        y_loc = y_range(j);

        for k = 1:z_count
            z_loc = z_range(k);

            cell_values = zeros(no_transmits,1);
            for T = 1:no_transmits
                element_locs = rx_element_locs{T};
                scans = all_scans{T};

                [full_dist, lateral_dist, tx_dist] = calc_dist(x_loc,y_loc,z_loc,element_locs,src_loc,transmit_type); 
                row_indices = round((full_dist/c + pulse_delay)*fs);

                row_indices(row_indices > length(scans)) = 1;
                indices = sub2ind(size(scans), row_indices, 1:length(element_locs))';

                lateral_ratio = lateral_dist./max_dist;
                lateral_ratio(lateral_ratio>1)=1;
                apros = rx_apro(round(lateral_ratio*(length(rx_apro)-1))+1);
                values = scans(indices);
                % values = values.*apros;
                % values = values * (tx_dist/max_dist);

                cell_values(T) = sum(values);

            end
            
            variance_volume(i,j,k) = var(cell_values);
            volume(i,j,k) = sum(cell_values);

        end
    end

    % loop_average = (2*loop_average + toc)/3;
end
elapsedTime = toc; 
fprintf('Elapsed Time: %.2f seconds\n', elapsedTime);

%%

dr = 60;

combined_volume = abs(volume);

combined_volume = combined_volume/max(max(max(combined_volume)));
variance_volume = variance_volume/max(max(max(variance_volume)));


processed_volume = 20*log10(combined_volume);
processed_volume = max(processed_volume, -dr);

slice = squeeze(processed_volume(:,y_mid,:))';

[~, ind ]= max(slice,[],"all");
[peak_row, peak_column] = ind2sub(size(slice),ind);
%
% % For cyc = 2
% target_row = peak_row + 4;
% target_column = peak_column - 9;

target_row = peak_row;
target_column = peak_column;


figure(2)
colormap("gray");
x_axis = [y_min x_max];
z_axis = [z_min z_max];
imagesc(x_axis, z_axis, slice);
axis("image")
clim([-dr 0]); colorbar;
xlabel("Horizontal Distance (cm)"); ylabel("Vertical Distance (cm)");
% yline(z_range(target_row),"Color","yellow")
% xline(y_range(target_column), "Color", "yellow")

%%
% plot the yz-plane from the volume
figure(3)
subplot 411

plot(emission_time, emission);
title("Transmitted Signal");

subplot 412
colormap("gray");
x_axis = [x_min x_max];
z_axis = [z_min z_max];
imagesc(x_axis, z_axis, slice);
clim([-dr 0]); colorbar;
xlabel("Horizontal Distance (cm)"); ylabel("Vertical Distance (cm)");


% yline(z_range(target_row),"Color","yellow")
subplot 413
trace = slice(target_row,:);
plot(x_range, trace);

subplot 414
trace_c = slice(:,target_column);
plot(z_range, trace_c);








