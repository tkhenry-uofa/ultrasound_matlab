
dynamic_range = 80;

patch_size = 32;
search_area = 256;
motion_grid_size = 8;

reference_frame = low_res_array{4};
compressed_reference_frame = process_volume(reference_frame,dynamic_range);
image_count = length(low_res_array);

image_size = size(reference_frame);
row_count = image_size(1);
col_count = image_size(2);

patch_row_corners = 1:motion_grid_size:(row_count - patch_size + 1);
patch_col_corners = 1:motion_grid_size:(col_count - patch_size + 1);

patch_x_count = length(patch_col_corners);


for G = 1:image_count
    
    target_image = low_res_array{G};
    compressed_target_image = process_volume(target_image, dynamic_range);






end