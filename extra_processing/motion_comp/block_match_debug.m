%% Setup: Load and preprocess image
G = 2;  % Index of image you want to warp
test_dr = 80;

patch_size = 32;
search_area = 200;
motion_grid_size = 8;
threshold = 1.05;

reference_frame = 5;
raw_ref = low_res_array{reference_frame};
raw_current = low_res_array{G};

ref_mag = abs(raw_ref);     % Use raw magnitude for motion estimation
current_mag = abs(raw_current);

image_size = size(raw_ref);

[M, N] = size(ref_mag);
[X, Y] = meshgrid(1:N, 1:M);

num_y = floor((M - patch_size + 1) / motion_grid_size);
num_x = floor((N - patch_size + 1) / motion_grid_size);
motion_array = zeros(num_y, num_x, 2);

%% Parameters
threshold = 1.05;          % Peak must exceed no-shift point
sharpness_thresh = 0.1;    % Peak must rise above neighborhood

%% Step 1: Block-based motion estimation on magnitude with confidence
for i = 1:num_y
    for j = 1:num_x
        r0 = (i-1)*motion_grid_size + 1;
        c0 = (j-1)*motion_grid_size + 1;

        r_patch = r0:r0+patch_size-1;
        c_patch = c0:c0+patch_size-1;

        r_search = max(1, r0 - search_area):min(M, r0 + patch_size + search_area - 1);
        c_search = max(1, c0 - search_area):min(N, c0 + patch_size + search_area - 1);

        ref_patch = ref_mag(r_patch, c_patch);
        search_block = current_mag(r_search, c_search);

        if var(ref_patch(:)) == 0
            continue;
        end

        % Normalize
        ref_patch = ref_patch - mean(ref_patch(:));
        search_block = search_block - mean(search_block(:));

        corr = normxcorr2(ref_patch, search_block);
        no_shift_pt = [r0, c0] - [r_search(1), c_search(1)] + size(ref_patch);
        no_shift_val = corr(no_shift_pt(1), no_shift_pt(2));

        [peak_val, peak_idx] = max(corr(:));
        [p_row, p_col] = ind2sub(size(corr), peak_idx);

        % Get correlation matrix size
        [Hcorr, Wcorr] = size(corr);
        
        % Sharpness check (local 3x3 window)
        r_win = max(1, p_row-1):min(Hcorr, p_row+1);
        c_win = max(1, p_col-1):min(Wcorr, p_col+1);
        local_vals = corr(r_win, c_win);
        local_mean = mean(local_vals(:));
        sharpness = peak_val - local_mean;

        if peak_val > no_shift_val * threshold && sharpness > sharpness_thresh
            motion_array(i,j,:) = [p_row, p_col] - no_shift_pt;
        end
    end
end


%% Step 2: Upsample motion field and smooth
full_motion = imresize(motion_array, [M, N], 'bilinear');

% Optional: smooth motion
full_motion(:,:,1) = imgaussfilt(full_motion(:,:,1), 2);
full_motion(:,:,2) = imgaussfilt(full_motion(:,:,2), 2);

% Visualize motion field
figure;
quiver(X(1:8:end,1:8:end), Y(1:8:end,1:8:end), ...
       full_motion(1:8:end,1:8:end,2), ...
      -full_motion(1:8:end,1:8:end,1));
axis ij; title("Motion Field");

%% Step 3: Reverse warp complex-valued source image
sample_x = X - full_motion(:,:,2);  % dx
sample_y = Y - full_motion(:,:,1);  % dy

sample_x = min(max(sample_x, 1), N);
sample_y = min(max(sample_y, 1), M);

% Interpolate complex-valued image
real_interp = interp2(X, Y, real(raw_current), sample_x, sample_y, 'linear', 0);
imag_interp = interp2(X, Y, imag(raw_current), sample_x, sample_y, 'linear', 0);
warped_image = complex(real_interp, imag_interp);

%% Step 4: Visualize results
figure;
subplot(1,3,1); imagesc(process_volume(raw_current,test_dr)); title("Original"); axis image;
subplot(1,3,2); imagesc(process_volume(raw_ref,test_dr)); title("Reference"); axis image;
subplot(1,3,3); imagesc(process_volume(warped_image,test_dr)); title("Warped to Ref"); axis image;
colormap gray;
