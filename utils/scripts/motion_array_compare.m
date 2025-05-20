
image_size = size(processed_shifted_image);
cell_index = 14;

motion_array = motion_cell{cell_index};

motion_array = imgaussfilt(motion_array,1.0);

pulse_difference = abs(cell_index - reference_readi_group) * double(readi_group_size);
dt = pulse_difference/prf;
dx_pixel = lateral_resolution;

motion_array = motion_array .* lateral_resolution ./ dt;


motion_size = size(motion_array);
full_motion_array = imresize(motion_array, image_size, 'bilinear');

% full_motion_array = imgaussfilt(full_motion_array,2);

x_motion = squeeze(motion_array(:,:,2));
z_motion = squeeze(motion_array(:,:,1));
[mX, mZ] = meshgrid(1:motion_size(2),1:motion_size(1));

full_x_motion = squeeze(full_motion_array(:,:,2));
full_z_motion = squeeze(full_motion_array(:,:,1));

full_motion = sqrt(full_z_motion.^2 + full_x_motion.^2);

% cross_section_row = 26;
cross_section_row = 360;
% flow_velocity = sqrt(x_motion(:,cross_section_row).^2 + z_motion(:,cross_section_row).^2);
flow_velocity = abs(full_x_motion(:,cross_section_row));


figure();
% 
% % contour(X,Y,full_x_motion);
% 
% % subplot 121
% quiver(mX,-1.*mZ,x_motion,z_motion,0.9);
% 
% % figure();
% imagesc(abs(full_x_motion));
% % quiver(X,-1.*Y,full_x_motion,full_z_motion);
% axis('image');
% title("Motion Map")
% % 
% figure();
% % subplot 122
% plot(z_range * 1000, flow_velocity * 100);
% title("1.2 cm/s Flow Profile")
% xlabel("Depth (mm)");
% ylabel("Speed (cm/s)");



subplot 221
imagesc(processed_low_res{cell_index});
% imagesc(x_range*1000, z_range*1000,(processed_low_res{cell_index}));
axis('image');
title("Original");

subplot 222
imagesc(processed_low_res{reference_readi_group});
% imagesc(x_range*1000, z_range*1000,(processed_low_res{reference_frame}));
axis('image');
title("Reference");

subplot 223
quiver(mX,-1.*mZ,x_motion,z_motion);
% quiver(x_downsampled,-1.*z_downsampled,x_motion,z_motion);
% quiver(X,-1.*Y,x_motion,z_motion);
axis('image');
title("Motion Map")

subplot 224
% imagesc(x_range*1000, z_range*1000,(processed_shifted_array{cell_index}));
% imagesc(processed_shifted_array{cell_index});
% axis('image');
% title("Shifted");
plot(z_range * 1000, flow_velocity * 100);