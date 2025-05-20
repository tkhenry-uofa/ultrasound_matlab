
motion_array = motion_cell{1};

U = squeeze(motion_array(:,:,1));
V = squeeze(motion_array(:,:,2));


figure();
quiver(U,V);