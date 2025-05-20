% Finds brightest point on a slice and plots all the RF data at that delay
% Run after running the main simulation

close all;


if(exist('vol_config','var'))
    y_mid = vol_config.y_mid;
    num_transmits = tx_config.no_transmits;
else
    num_transmits = no_transmits;
end
slice = squeeze(processed_volume(:,y_mid,:))';

    
[~, ind ]= max(slice,[],"all");
[peak_row, peak_column] = ind2sub(size(slice),ind);
element_locs = rx_element_locs{1};

i = peak_column;
k = peak_row;
% 
i = 118;
k = 99;
% % k= 62;


if(exist('cross_sequence','var'))
    rx_sequence = cross_sequence;
    x_loc = x_range(i);
    y_loc = y_range(y_mid);
    z_loc = z_range(k);
    
else
    x_loc = vol_config.x_range(i);
    y_loc = vol_config.y_range(vol_config.y_mid);
    z_loc = vol_config.z_range(k);
end

vox_loc = [x_loc, y_loc, z_loc];
element_count = size(element_locs,2);

window_length = 2* length(returned_time);
window_mid = round(window_length/2);

cell_values = zeros(num_transmits,1);
traces = zeros(num_transmits, element_count, 4*window_length + 1);
for T = 1:num_transmits
    element_locs = rx_element_locs{T};

    scans = data_array{T};
    sample_count = size(scans,1);

    [full_dist, lateral_dist, tx_dist] = calc_dist(vox_loc,element_locs,src_loc,transmit_type);

    row_indices = round((full_dist/c + pulse_delay)*fs);

    row_indices(row_indices > length(scans)) = 1;

    for e = 1:element_count
        range = (row_indices(e) - window_length):(row_indices(e) + 3*window_length);
        
        range(range<1) = 1;
        range(range>sample_count) = 1;
        traces(T,e,:) = scans(range,e);
    end

end


% test = squeeze(traces(1,200,:));
% L = length(test);
% mid = round(L/2)+1;
% Y = fftshift(fft(test));
% Ya = abs(Y);
% f_range = fs/L*(-L/2:L/2-1);
% plot(f_range,Ya)
% xlabel("f (Hz)")
% ylabel("|fft(X)|")
% 
% [~,index ]=max(Ya);
% f = f_range(index)

%%

plot_slice(slice,1);
xline(i,Color="red")
yline(k,Color="red")

title_string = sprintf('X location: %0.2g', vol_config.x_range(i));

% Get the nearest line to the array edges
max_i = interp1(vol_config.x_range, 1:vol_config.x_count, tx_config.x_max,"linear",vol_config.x_count);
min_i = interp1(vol_config.x_range, 1:vol_config.x_count, tx_config.x_min,"linear",1);

xline(max_i,Color="yellow");
xline(min_i,Color="yellow")

title(title_string)

% figure(2);
% plot(returned_time, returned_wave);
% title("Expected Psf");


trace_fig = figure(3);
array_fig = figure(4);
total_fig = figure(5);

[AX,AY] = ndgrid(tx_config.x_range, tx_config.y_range);
grid_points = [AX(:), AY(:)];
array_pic = zeros(size(AX));
total = zeros(1,size(traces,3));
for i = 1:element_count
    figure(trace_fig);
    clf(trace_fig)
    hold on;
    offset = 0;  % Initialize offset
    sample_locs = zeros(num_transmits,2);
    for T = 1:4:num_transmits
        sample_locs(T,:) = rx_element_locs{T}(1:2,i);

        trace = squeeze(traces(T,i,:));
        sample = real(trace);
        abs_sample = abs(trace);
        
        total = total + trace;

        plot(sample + offset, Color="blue");
        plot(abs_sample + offset, Color="red",LineStyle="--")
        [~, peak] = max(abs_sample);
        scatter(peak, abs_sample(peak)+offset)
        offset = offset + max(abs_sample) - min(sample);  % Adjust offset

    end
    xline(window_length+1)
    title("Element: ", i)
    pause(0.1)
    hold off

    figure(array_fig);
    [~, loc_index] = ismember(round(grid_points,4), round(sample_locs,4),"rows");
    array_pic = array_pic + reshape(loc_index>0, size(AX));

    imagesc(tx_config.y_range, tx_config.x_range, array_pic)
    axis("image")

    figure(total_fig)
    clf(total_fig)
    title("Running Total")
    hold on 
    total_slice = total(1:2*window_length);
    plot(real(total_slice),Color="blue");
    plot(abs(total_slice),Color="red",LineStyle="--")
    xline(window_length+1)
    
    hold off

end

