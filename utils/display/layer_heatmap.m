function composite_rgb = layer_heatmap(base_image,heatmap, row_range, col_range, alpha)

    
    alpha_map = zeros(size(base_image));
    alpha_map(row_range,col_range) = alpha;
    
    base_norm = mat2gray(base_image);
    color_norm  = mat2gray(heatmap);
    
    % 2. Map each to their respective colormaps
    % 256 is default size of MATLAB colormaps
    gray_map = gray(256);
    color_map = hot(256);
    
    base_mapped = ind2rgb(round(base_norm * 255) + 1, gray_map);
    top_mapped  = ind2rgb(round(color_norm  * 255) + 1, color_map);
    
    % 3. Blend the two using per-pixel alpha
    % alpha_map = mat2gray(alpha_map); % Ensure alpha between 0 and 1
    composite_rgb = (1 - alpha_map) .* base_mapped + alpha_map .* top_mapped;

end