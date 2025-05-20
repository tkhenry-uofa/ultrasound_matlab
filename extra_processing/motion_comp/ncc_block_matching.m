%% NCC Block Matching
% Parameters (adjust as needed)

% 250411_C1_reso_hand_lateral_FORCES-TxColumn_02
% patch_size = 64;         
% search_area = 64;
% motion_grid_size = 8; 
% var_min = 0.016;
% threshold = 1.05;
% test_dr = 60;
% reference_frame = 3;

% 250410_MN32-5_reso_axial_motion_FORCES-TxRow_01
% patch_size = 16;         
% search_area = 64;
% motion_grid_size = 8; 
% var_min = 0.016;
% threshold = 1.1;
% test_dr = 60;
% reference_frame = 3;

% 3 m/s simulation
% patch_size = 64;         
% search_area = 256;
% motion_grid_size = 16; 
% var_min = 0.0;
% threshold = 1.1;
% test_dr = 40;
% reference_frame = 8;

% 250424_MN32-5_flow_6mm_10_FORCES-TxColumn
% frame 10
patch_size = 32;         
search_area = 16;
vertical_margin = 16;
horizontal_margin = 16;
motion_grid_size = 4; 
var_min = 0.005;
threshold = 1.05;
test_dr = 35;
% forces_frame_id = 8;
power_compress = false;
power_threshold = 3000;


% 250424_MN32-5_flow_6mm_60_FORCES-TxColumn
% frame 10
% forces_frame_id = 9;
% patch_size = 32;         
% vertical_margin = 16;
% horizontal_margin = 78;
% motion_grid_size = 8; 
% var_min = 0.01;
% threshold = 1.1;
% test_dr = 35;
% power_compress = true;
% power_threshold = 130;




reference_readi_group = 8;


image_array = filtered_frame_cell(forces_frame_id,:);


raw_reference_image = image_array{reference_readi_group};
reference_image = process_volume(raw_reference_image,test_dr,power_threshold,power_compress);


[row_count, col_count] = size(reference_image);
image_size = [row_count,col_count];
[X,Y] = meshgrid(1:col_count, 1:row_count);
num_patches_y = ceil((row_count - patch_size + 1) / motion_grid_size);
num_patches_x = ceil((col_count - patch_size + 1) / motion_grid_size);

motion_row_start = floor(vessel_row_start / motion_grid_size);
motion_row_end = ceil(vessel_row_end / motion_grid_size);

motion_array = zeros(num_patches_y, num_patches_x, 2);

motion_cell = cell(1, numel(image_array));
motion_cell{1} = motion_array;

fprintf("Starting motion compenstation\n");

shifted_images = cell(size(motion_cell));

% Process each non-reference image

parfor G = 1:numel(image_array)

    current_image = image_array{G};
    motion_array = zeros(num_patches_y, num_patches_x,2);

    if G == reference_readi_group
        motion_cell{G} = motion_array;
        shifted_images{G} = current_image;
        continue;
    end
    current_image_processed = process_volume(current_image,test_dr,power_threshold,power_compress);
    
    current_shifted_image = complex(zeros(image_size));

    fprintf("Image: %i\n", G);

    peak_index = [0 0];

    current_h_margin = horizontal_margin;% * (abs(G - reference_readi_group));

    % Loop over each patch in the current image
    for i = motion_row_start:motion_row_end

        for j = 1:num_patches_x
            
            patch_row_start = (i-1)*motion_grid_size + 1;
            patch_col_start = (j-1)*motion_grid_size + 1;

            patch_row_end = min(row_count, patch_row_start + patch_size - 1);
            patch_col_end = min(col_count, patch_col_start + patch_size - 1);

            patch_row_range = patch_row_start:patch_row_end;
            patch_col_range = patch_col_start:patch_col_end;

            patch_row_count = length(patch_row_range);
            patch_col_count = length(patch_col_range);

            search_row_start = max(1, patch_row_start - vertical_margin);
            search_row_end = min(row_count, patch_row_start + patch_size + vertical_margin - 1);
            search_row_range = search_row_start:search_row_end;

            search_col_start = max(1, patch_col_start - current_h_margin);
            search_col_end = min(col_count, patch_col_start + patch_size + current_h_margin - 1);
            search_col_range = search_col_start:search_col_end;

            patch_corner = [patch_row_start, patch_col_start];
            search_corner = [search_row_start, search_col_start];

            reference_patch = current_image_processed(patch_row_range, patch_col_range);
            search_patch = reference_image(search_row_range, search_col_range);
            
            reference_var = var(reshape(reference_patch,[],1));
            search_var = var(reshape(search_patch,[],1));

            if(reference_var ~= 0)
                % norm_corr_mat = corr_mat/ sqrt(patch_var * ref_var);
                % no_shift_point = patch_corner - search_corner + size(reference_patch);
                no_shift_point = patch_corner - search_corner + [1,1];
                norm_corr_mat = normxcorr2(reference_patch,search_patch);
                 
                norm_corr_mat = norm_corr_mat(patch_row_count:(end-patch_row_count+1),patch_col_count:(end-patch_col_count+1));

                % corr_mat = xcorr2(reference_patch,search_patch);
    
                no_shift_c = norm_corr_mat(no_shift_point(1), no_shift_point(2));
                corr_var = var( reshape(norm_corr_mat,[],1));
                [peak_c, index] = max(norm_corr_mat, [],"all");

                % if(i == 34 && j == 30 )
                % 
                %     figure();
                %     imagesc(norm_corr_mat);
                %     fprintf("\n");
                % end
    
                if (peak_c > no_shift_c * threshold) && (corr_var > var_min)
                    [peak_index(1),peak_index(2)] = ind2sub(size(norm_corr_mat),index);
                    motion_vector = peak_index - no_shift_point;
                    motion_array(i,j,:) = motion_vector;
                else
                    motion_vector = [0,0];
                end
            else
                motion_vector = [0,0];
            end

            dy = motion_vector(1);
            dx = motion_vector(2);

            % New top-left corner after shift
            new_row = patch_row_start + dy;
            new_col = patch_col_start + dx;
            new_row_end = patch_row_end + dy;
            new_col_end = patch_col_end + dx;

            % Compute valid destination bounds in image
            dst_row1 = max(1, new_row);
            dst_col1 = max(1, new_col);
            dst_row2 = min(row_count, new_row_end);
            dst_col2 = min(col_count, new_col_end);

            src_row1 = dst_row1 - dy;
            src_col1 = dst_col1 - dx;
            src_row2 = dst_row2 - dy;
            src_col2 = dst_col2 - dx;

            % Add patch to image with bounds checking
            current_shifted_image(dst_row1:dst_row2, dst_col1:dst_col2) = ...
                current_shifted_image(dst_row1:dst_row2, dst_col1:dst_col2) + ...
                current_image(src_row1:src_row2, src_col1:src_col2);


        end
    end
 
    % Save results
    motion_cell{G} = motion_array;
    shifted_images{G} = current_shifted_image;
end
