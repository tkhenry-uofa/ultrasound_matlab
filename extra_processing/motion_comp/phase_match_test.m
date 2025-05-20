%% Phase Correlation Block Matching and Reverse Warping
% Parameters

patch_size = 32;         % Patch size (square)
search_area = 128;
motion_grid_size = 16; 
margin = 0; % Stay this far away from the image edge


threshold = 1.1; 
sharpness_thresh = 0.01;
test_dr = 50;

reference_frame = 5;
raw_reference_image = low_res_array{reference_frame};
reference_image = abs(raw_reference_image); % use magnitude only

[row_count, col_count] = size(reference_image);
image_size = [row_count, col_count];
cropped_image_size = [row_count - 2*margin, col_count - 2*margin];
[X,Y] = meshgrid(1:col_count, 1:row_count);
num_patches_y = floor((row_count - patch_size - 2*margin + 1 ) / motion_grid_size);
num_patches_x = floor((col_count - patch_size - 2*margin + 1 ) / motion_grid_size);

motion_cell = cell(1, numel(low_res_array));
shifted_images = cell(size(motion_cell));

fprintf("Starting motion compensation\n");

% Process each non-reference image
parfor G = 1:numel(low_res_array)
    current_image = low_res_array{G};
    motion_array = zeros(num_patches_y, num_patches_x, 2);
    fprintf("Image: %i\n", G);

    % Loop over each patch in the current image
    for i = 1:num_patches_y
        for j = 1:num_patches_x

            patch_row_start = margin + (i-1)*motion_grid_size + 1;
            patch_col_start = margin + (j-1)*motion_grid_size + 1;
            patch_row_range = patch_row_start:patch_row_start+patch_size-1;
            patch_col_range = patch_col_start:patch_col_start+patch_size-1;

            search_row_start = max(1, patch_row_start - search_area);
            search_row_end = min(row_count, patch_row_start + patch_size + search_area - 1);

            search_col_start = max(1, patch_col_start - search_area);
            search_col_end = min(col_count, patch_col_start + patch_size + search_area - 1);

            search_row_range = search_row_start:search_row_end;
            search_col_range = search_col_start:search_col_end;

            ref_patch = abs(raw_reference_image(patch_row_range, patch_col_range));
            search_patch = abs(current_image(search_row_range, search_col_range));

            if var(ref_patch(:)) == 0
                continue;
            end

            % Phase correlation
            fft_ref = fft2(ref_patch, size(search_patch,1), size(search_patch,2));
            fft_search = fft2(search_patch);
            R = fft_ref .* conj(fft_search);
            R = R ./ max(abs(R), eps);
            pcorr = real(ifft2(R));

            [peak_val, peak_idx] = max(pcorr(:));
            [p_row, p_col] = ind2sub(size(pcorr), peak_idx);

            [Hcorr, Wcorr] = size(pcorr);
            r_win = max(1, p_row-1):min(Hcorr, p_row+1);
            c_win = max(1, p_col-1):min(Wcorr, p_col+1);
            local_vals = pcorr(r_win, c_win);
            sharpness = peak_val - mean(local_vals(:));

            no_shift_point = [patch_row_start, patch_col_start] - [search_row_start, search_col_start] + size(ref_patch);
            no_shift_value = pcorr(no_shift_point(1),no_shift_point(2));

            if (sharpness > sharpness_thresh) 
                
                motion_vector = [p_row, p_col] - no_shift_point;
                motion_array(i,j,:) = motion_vector;
            end
        end
    end

    % Upsample motion to full resolution
    full_motion_array = imresize(motion_array, cropped_image_size, 'bicubic');
    full_motion_array(:,:,1) = imgaussfilt(full_motion_array(:,:,1), 2); % dy
    full_motion_array(:,:,2) = imgaussfilt(full_motion_array(:,:,2), 2); % dx

    full_motion_array = padarray(full_motion_array,[margin, margin],0,"both");

    % Reverse warp
    [X, Y] = meshgrid(1:col_count, 1:row_count);
    sample_x = X - full_motion_array(:,:,2);
    sample_y = Y - full_motion_array(:,:,1);
    sample_x = min(max(sample_x, 1), col_count);
    sample_y = min(max(sample_y, 1), row_count);

    real_interp = interp2(X, Y, real(current_image), sample_x, sample_y, 'cubic', 0);
    imag_interp = interp2(X, Y, imag(current_image), sample_x, sample_y, 'cubic', 0);
    current_shifted_image = complex(real_interp, imag_interp);

    motion_cell{G} = motion_array;
    shifted_images{G} = current_shifted_image;
end