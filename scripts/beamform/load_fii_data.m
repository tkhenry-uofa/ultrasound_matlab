%% Setup

clear all;

load("data\readi\hercules_psf.mat");
% rx_scans = rx_scans * 1e27; % Scale up fieldii values
addpath(genpath(pwd));


do_hercules = false;
readi_group_count = 4;

x_range = [-20, 20]/1000;
y_range = [-20, 20]/1000;
z_range = [55, 65]/1000;

resolution = 0.0001; % Spatial voxel size

vol_config = volume_config(x_range,y_range,z_range,resolution);

%% Beamforming
fprintf("Starting readi beamform\n")
low_res_cells = cuda_processing_f2(rx_scans,tx_config,readi_group_count,vol_config); 


if do_hercules
    fprintf("Starting hercules beamform\n")
    hercules_cell = cuda_processing_f2(rx_scans,tx_config,1,vol_config); 
    hercules_volume = hercules_cell{1};
end


%% Post processing 
volume_size = vol_config.size;
high_res_volume = zeros(volume_size);

dynamic_range = 50;

processed_low_res = low_res_cells;
low_res_images = low_res_cells;
for g = 1:readi_group_count
    high_res_volume = high_res_volume + low_res_cells{g};
    processed_low_res{g} = process_volume(low_res_cells{g}, dynamic_range); 

    low_res_images{g} = squeeze(processed_low_res{g}(:,vol_config.y_mid,:));
    % low_res_images{g} = squeeze(processed_low_res{g}(vol_config.x_mid, :,:));
end


processed_readi = process_volume(high_res_volume, dynamic_range);

if do_hercules
    processed_hercules = process_volume(hercules_volume, dynamic_range);
end




image = squeeze(processed_readi(:, vol_config.y_mid,:)).';
% image = squeeze(processed_readi(vol_config.x_mid, :,:)).';

figure();

imagesc(image);
title("Readi High Res Image",'FontSize',16);
colormap("gray");
axis('image');

if do_hercules
    hercules_image = squeeze(processed_hercules(:, vol_config.y_mid,:)).';
    % hercules_image = squeeze(processed_hercules(vol_config.x_mid,:,:)).';

    figure();
    imagesc(hercules_image);
    title("Hercules Image",'FontSize',16);
    colormap("gray");
    axis('image');
end



%% 

figure();

% Main title
% sgtitle(speed_str + ' Readi Low Resolution Images','FontSize',16);
sgtitle('Stationary Readi Low Resolution Images','FontSize',16);

% Loop through the images and display each in the 2x8 grid
for i = 1:readi_group_count
    subplot(2, 2, i); 
    colormap("gray");
    imagesc(vol_config.x_range*1000, vol_config.z_range*1000, (low_res_images{i}).');
    title(['Image ' num2str(i)]);
    xlabel("mm");
    ylabel("mm");
    
    clim([-dynamic_range, 0]); 
    

end


%% Extra 


combined_volume = low_res_cells{1} + low_res_cells{2};
combined_volume = process_volume(combined_volume, dynamic_range);
