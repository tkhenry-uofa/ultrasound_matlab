%% NCC Block Matching
% Parameters (adjust as needed)

% 250411_C1_reso_hand_lateral_FORCES-TxColumn_02
patch_size = 64;         
search_margin = 64;
motion_grid_size = 8; 
var_min = 0.016;
threshold = 1.00;
test_dr = 60;
reference_frame = 3;

% 250410_MN32-5_reso_axial_motion_FORCES-TxRow_01
% patch_size = 32;         
% search_area = 32;
% motion_grid_size = 8; 
% var_min = 0.01;
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


search_size = search_margin * 2 + patch_size;
corr_size = search_size + patch_size - 1;

cropped_corr_size = search_size - patch_size + 1;

corr_range = linspace(-1,1,cropped_corr_size);

bias_variance = 0.25;

bias_envelope = exp(-corr_range.^2./(2*bias_variance));
bias_env_2d = bias_envelope .* bias_envelope.';
% figure(); plot(bias_envelope)


raw_reference_image = low_res_array{reference_frame};
reference_image = process_volume(raw_reference_image,test_dr);


[row_count, col_count] = size(reference_image);
image_size = [row_count,col_count];
[X,Y] = meshgrid(1:col_count, 1:row_count);
num_patches_y = ceil((row_count - patch_size + 1) / motion_grid_size);

num_patches_x = ceil((col_count - patch_size + 1) / motion_grid_size);

motion_array = zeros(num_patches_y, num_patches_x, 2);

motion_cell = cell(1, numel(low_res_array));
motion_cell{1} = motion_array;

fprintf("Starting motion compenstation\n");

shifted_images = cell(size(motion_cell));

% Process each non-reference image
parfor G = 1:numel(low_res_array)
    current_image = low_res_array{G};

    current_image_processed = process_volume(current_image,test_dr);
    motion_array = zeros(num_patches_y, num_patches_x,2);

    current_shifted_image = complex(zeros(image_size));

    fprintf("Image: %i\n", G);

    peak_index = [0 0];

    % Loop over each patch in the current image
    for i = 1:num_patches_y

        for j = 1:num_patches_x

            patch_row_start = (i-1)*motion_grid_size + 1;
            patch_col_start = (j-1)*motion_grid_size + 1;

            patch_row_end = min(row_count, patch_row_start + patch_size - 1);
            patch_col_end = min(col_count, patch_col_start + patch_size - 1);

            patch_row_range = patch_row_start:patch_row_end;
            patch_col_range = patch_col_start:patch_col_end;

            patch_row_count = length(patch_row_range);
            patch_col_count = length(patch_col_range);

            search_row_start = max(1, patch_row_start - search_margin);
            search_row_end = min(row_count, patch_row_start + patch_size + search_margin - 1);
            search_row_range = search_row_start:search_row_end;

            search_col_start = max(1, patch_col_start - search_margin);
            search_col_end = min(col_count, patch_col_start + patch_size + search_margin - 1);
            search_col_range = search_col_start:search_col_end;

            patch_corner = [patch_row_start, patch_col_start];
            search_corner = [search_row_start, search_col_start];

            reference_patch = reference_image(patch_row_range, patch_col_range);
            search_patch = current_image_processed(search_row_range, search_col_range);
%%
            if i == 3000 && j == 32
                inspect_patch
            end
            
            reference_var = var(reshape(reference_patch,[],1));
            search_var = var(reshape(search_patch,[],1));

            if(reference_var ~= 0)

                no_shift_point = patch_corner - search_corner + [1,1];
                norm_corr_mat = normxcorr2(reference_patch,search_patch);
               
                % Discarding the edges where the patches don't completely
                % overlap
                norm_corr_mat_cropped = norm_corr_mat(patch_row_count:(end-patch_row_count+1),patch_col_count:(end-patch_col_count+1));
                
                if(size(norm_corr_mat_cropped) == size(bias_env_2d))
                    norm_corr_mat_cropped = norm_corr_mat_cropped .* bias_env_2d;
                end
                

                no_shift_c = norm_corr_mat_cropped(no_shift_point(1), no_shift_point(2));
                corr_var = var( reshape(norm_corr_mat_cropped,[],1));
                [peak_c, index] = max(norm_corr_mat_cropped, [],"all");
    
                if (peak_c > no_shift_c * threshold) && (corr_var > var_min)
                    [peak_index(1),peak_index(2)] = ind2sub(size(norm_corr_mat_cropped),index);
                    motion_vector = peak_index - no_shift_point;
                    motion_array(i,j,:) = motion_vector;
                else
                    motion_vector = [0,0];
                end
            else
                motion_vector = [0,0];
            end
        end
    end

    % Upsample motion to full resolution (nearest keeps patch edges sharp)
    full_motion_array = imresize(motion_array, image_size, 'nearest');

    % Optional: Smooth the flow a bit (edge-aware or Gaussian)
    full_motion_array(:,:,1) = imgaussfilt(full_motion_array(:,:,1), 2);  % dy
    full_motion_array(:,:,2) = imgaussfilt(full_motion_array(:,:,2), 2);  % dx

    % Build reverse sampling grid
    [X, Y] = meshgrid(1:col_count, 1:row_count);

    % Reverse warping: pull from source location
    sample_x = X - full_motion_array(:,:,2);  % subtract dx
    sample_y = Y - full_motion_array(:,:,1);  % subtract dy

    % Clamp sampling points to valid range (optional but safe)
    sample_x = min(max(sample_x, 1), col_count);
    sample_y = min(max(sample_y, 1), row_count);

    % Interpolate real & imaginary parts separately (interp2 does not support complex)
    real_interp = interp2(X, Y, real(current_image), sample_x, sample_y, 'nearest', 0);
    imag_interp = interp2(X, Y, imag(current_image), sample_x, sample_y, 'nearest', 0);
    current_shifted_image = complex(real_interp, imag_interp);


    % Save results
    motion_cell{G} = motion_array;
    shifted_images{G} = current_shifted_image;
end
