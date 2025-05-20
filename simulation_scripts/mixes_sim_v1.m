%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script simulates a a 2x2 array of 64x64 TOBE arrays side by side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close("all")

tic;

path(path, "C:\Users\tkhen\OneDrive\Documents\MATLAB\lab\Field_II_ver_3_30_windows")
field_init(-1)

fs = 100e6; % Sampling frequency (Hz)
f0 = 2.5e6;  % Transducer center frequency (Hz)
c = 1540;   % Speed of sound (m/s)
lambda = c/f0; % Wavelength (m)
row_count = 128; % Total from two side by side arrays

pitch = lambda;
kerf_x = 10e-6; % Kerf in x-direction (m)
kerf_y = 10e-6; % Kerf in y-direction
width = pitch - kerf_x; % Width of the element (m)
height = width; % Height of the elements (m)


%% Volume configuration
% Volume range (m)
x_min = -20/1000;
x_max = 20/1000; 

y_min = -20/1000;
y_max = 20/1000;

z_min = 20/1000;
z_max = 80/1000;

resolution = 0.0001; % One cell per 100 microns

x_range = x_min:resolution:x_max;
y_range = y_min:resolution:y_max;
z_range = z_min:resolution:z_max;

x_count = length(x_range);
y_count = length(y_range);
z_count = length(z_range);

volume = zeros(x_count, y_count, z_count);


%% Point Generation
% Creates a grid of lateralPoints in the xy plane
% at each level in zPoints.
xPoints = [0]/1000;
yPoints = [-10, 0, 10]/1000;
zPoints = [ 40, 55, 70 ]/1000;

p_x_len = length(xPoints);
p_y_len = length(yPoints);
p_z_len = length(zPoints);

pointArray = zeros(p_x_len,p_y_len,p_z_len,3);
for i = 1:p_x_len
    for j = 1:p_y_len
        for k = 1:p_z_len
            pointArray(i,j,k,:) = [xPoints(i) yPoints(j) zPoints(k)];
        end
    end
end

points = reshape(pointArray, [], 3);
amps = ones(p_z_len*p_y_len*p_x_len,1);

% Add an additional scatter far away to increase simulation length.
points = [points; [0 0 250]/1000];
amps = [amps; 0.001];


%% Transmit array generation
tx_focus = [0 0 1000];
tx_apro = ones(row_count, row_count);
tx_array = xdc_2d_array(row_count,row_count,width,height,kerf_x,kerf_y,tx_apro,1,1,tx_focus);

% Impulse, a basic cos wave with an exponential envelope 
cyc = 1;
tx_t = -cyc/f0:1/fs:cyc/f0;
impulse_cos = 2*cos(2*pi*f0*tx_t).*exp(-tx_t.^2/(cyc*0.5*(1/f0)^2));
xdc_impulse(tx_array,impulse_cos);

% Apodization
hamming_window = hamming(row_count)*hamming(row_count)';
ele_apodization(tx_array,(1:(row_count)*(row_count))',hamming_window(:));

% calc_int(tx_array,fs)

%% Receive array generation
% If true each cross on the 4 arrays is offset from its two neighbours
% so there is never a straight 128 element line
offset_crosses = 2; 
no_transmits = 8;

array_size = round(row_count/2);
active_rows = round(linspace(1,array_size-offset_crosses, no_transmits));

cross_sequence = diagonal_cross_sequence( no_transmits, array_size, active_rows, offset_crosses );

% print_sequence(cross_sequence)

%% Simulation
for T = 1:no_transmits
    fprintf("Transmission: %d \n", T)
    
    rx_cross = cross_sequence{T};
    rx_focus = tx_focus;
    rx_array = xdc_2d_array(row_count,row_count,width,height,kerf_x,kerf_y,rx_cross,1,1,rx_focus);

    xdc_impulse(rx_array,impulse_cos);
    element_locs = xdc_get(rx_array);
    element_locs = element_locs(24:26,:);
    element_count = length(element_locs);

    [scans, rx_start] = calc_scat_multi( tx_array, rx_array, points, amps );

    % Normalize
    scans = scans / max(max(scans));

    % Add zeros for the time before the first reception
    zero_scans = zeros(length(0:1/fs:rx_start), element_count);
    scans = [zero_scans; scans];

    % Image building
    fprintf("Volume Building\n")

    for i = 201
        fprintf("T: %d, I: %d\n", T,i);
        x_loc = x_range(i);
        for j = 1:y_count
            y_loc = y_range(j);
            for k = 1:z_count
                z_loc = z_range(k);

                [dist, ratio] = calc_dist(x_loc,y_loc,z_loc,element_locs); 

                index = round(dist/c*fs);
                index = (index) + (0:element_count-1) * size(scans,1);
                slices = scans(index);
                summed = sum(slices);

                volume(i,j,k) = volume(i,j,k) + summed;
            end
        end
    end
    elapsedTime = toc;     
end

fprintf('Elapsed Time: %.2f seconds\n', elapsedTime);

% Post processing
volume = abs(hilbert(volume));
volume = volume/max(max(max(volume)));
volume = 20*log10(volume);

%plot the yz-plane from the volume
figure(1)
colormap("gray");
x_axis = [x_min x_max];
z_axis = [z_min z_max];
imagesc(x_axis, z_axis, flip(squeeze(volume(201,:,:))'));
caxis([-40 0]); colorbar;
xlabel("Horizontal Distance (cm)"); ylabel("Vertical Distance (cm)");



function [full_path, ratio] = calc_dist(x_loc,y_loc,z_loc,Th_loc)

    delta_x = x_loc-Th_loc(1,:);
    delta_y = y_loc-Th_loc(2,:);
    
    %find location from each element to point
    dist = ((delta_x).^2 + (delta_y).^2 + (z_loc).^2).^(1/2);
    ratio = dist.^(-1) .* z_loc;
    
    full_path = dist + z_loc;

end

function sequence = diagonal_cross_sequence( no_transmits, array_size, active_rows, offset_crosses )
    
sequence = cell(no_transmits,1);
full = zeros(array_size*2);
    
    for i = 1:no_transmits
        % Initialize matrices for each cross
        cross1 = zeros(array_size); % No offset
        cross2 = zeros(array_size); % Offset in y-direction
        cross3 = zeros(array_size); % Offset in x-direction
        cross4 = zeros(array_size); % Offset in both directions
    
        % Define cross1 with no offsets
        cross1(active_rows(i), :) = 1;
        cross1(:, active_rows(i)) = 1;
    
        % Define cross2 with offset in y-direction only
        cross2((active_rows(i) + offset_crosses), :) = 1;
        cross2(:, active_rows(i)) = 1;
    
        % Define cross3 with offset in x-direction only
        cross3(active_rows(i), :) = 1;
        cross3(:, (active_rows(i) + offset_crosses)) = 1;
    
        % Define cross4 with offsets in both x and y directions
        cross4((active_rows(i) + offset_crosses), :) = 1;
        cross4(:, (active_rows(i) + offset_crosses)) = 1;
    
        % Arrange crosses into the final pattern
        cross = [cross1, cross2; cross3, cross4];
    
        % Store each cross sequence
        sequence{i} = cross;
        full = full + cross;
    end
end

function print_sequence(sequence)

    for i = 1:length(sequence)

        imagesc(sequence{i})
        pause(1)
    end

end