function [shifted_images, motion_cell] = forward_block_matching( ...
    image_array, reference_image_id, patch_size, m_grid_size, search_area, v_margin, h_margin, ...
    corr_threshold, min_var, dr, img_threshold)
    

    image_count = numel(image_array);
    raw_reference_image = image_array{reference_image_id};
    reference_image = process_volume(raw_reference_image,dr,img_threshold,false);
    shifted_images = cell(size(image_array));
    motion_cell = cell(size(image_array));
    
    
    [row_count, col_count] = size(reference_image);
    image_size = [row_count,col_count];
    total_patches = round(image_size/m_grid_size);
    total_patches_y = round(row_count / m_grid_size);
    total_patches_x = round(col_count / m_grid_size);

    search_patches = round(search_area/m_grid_size);

    % This many rows/columns overlap the edge of the image
    edge_patch_margin = floor((patch_size + m_grid_size) / (m_grid_size * 2));

    search_patches = max(search_patches, edge_patch_margin);
    search_patches = min(search_patches, total_patches - edge_patch_margin);

    patch_range_z = search_patches(1,1):search_patches(2,1);
    patch_range_x = search_patches(1,2):search_patches(2,2);
    
    parfor G = 1:image_count
        fprintf("Image: %i\n", G);

        current_image = image_array{G};
        motion_array = zeros(total_patches_y, total_patches_x,2);
    
        if G == reference_image_id
            motion_cell{G} = motion_array;
            shifted_images{G} = current_image;
            continue;
        end
        current_image_processed = process_volume(current_image,dr,img_threshold,false);
        current_shifted_image = complex(zeros(image_size));
    
        current_h_margin = h_margin;
        current_v_margin = v_margin;
    
        peak_index = [0 0];
        % TODO: Break this into an estimation and a warp pass
        for i = patch_range_z
    
            for j = patch_range_x
                
                patch_row_start = (i-1)*m_grid_size + 1;
                patch_col_start = (j-1)*m_grid_size + 1;
    
                patch_row_end = min(row_count, patch_row_start + patch_size - 1);
                patch_col_end = min(col_count, patch_col_start + patch_size - 1);
    
                patch_row_range = patch_row_start:patch_row_end;
                patch_col_range = patch_col_start:patch_col_end;
    
                patch_row_count = length(patch_row_range);
                patch_col_count = length(patch_col_range);
    
                search_row_start = max(1, patch_row_start - current_v_margin);
                search_row_end = min(row_count, patch_row_start + patch_size + current_v_margin - 1);
                search_row_range = search_row_start:search_row_end;
    
                search_col_start = max(1, patch_col_start - current_h_margin);
                search_col_end = min(col_count, patch_col_start + patch_size + current_h_margin - 1);
                search_col_range = search_col_start:search_col_end;
    
                patch_corner = [patch_row_start, patch_col_start];
                search_corner = [search_row_start, search_col_start];
    
                reference_patch = current_image_processed(patch_row_range, patch_col_range);
                search_patch = reference_image(search_row_range, search_col_range);
    
                reference_var = var(reshape(reference_patch,[],1));
    
                if(reference_var ~= 0)
                    no_shift_point = patch_corner - search_corner + [1,1];
                    corr_map = normxcorr2(reference_patch,search_patch);
                     
                    corr_map = corr_map(patch_row_count:(end-patch_row_count+1),patch_col_count:(end-patch_col_count+1));
        
                    no_shift_c = corr_map(no_shift_point(1), no_shift_point(2));
                    corr_var = var( reshape(corr_map,[],1));
                    [peak_c, index] = max(corr_map, [],"all");
        
                    if (peak_c > no_shift_c * corr_threshold) && (corr_var > min_var)
                        [peak_index(1),peak_index(2)] = ind2sub(size(corr_map),index);
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

end