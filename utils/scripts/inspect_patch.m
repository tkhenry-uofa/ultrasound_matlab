figure();

                % Reference and search patches
                subplot(2,2,1);
                imagesc(reference_patch);
                axis image;
                title("Reference Patch");
            
                subplot(2,2,2);
                imagesc(search_patch);
                axis image;
                title("Search Area");
                hold on;
                % Green box at the no-shift location
                rectangle('Position', [patch_col_range(1)-search_col_range(1)+1, ...
                                       patch_row_range(1)-search_row_range(1)+1, ...
                                       length(patch_col_range)-1, ...
                                       length(patch_row_range)-1], ...
                          'EdgeColor', 'g', 'LineWidth', 1.5);
            
                % Full reference image with red box
                subplot(2,2,3);
                imagesc(reference_image);
                axis image;
                title("Full Reference Image");
                hold on;
                rectangle('Position', [patch_col_range(1), patch_row_range(1), ...
                                       length(patch_col_range)-1, length(patch_row_range)-1], ...
                          'EdgeColor', 'r', 'LineWidth', 1.5);
                hold off;
            
                % Full current image with red (current) and green (no-shift) boxes
                subplot(2,2,4);
                imagesc(current_image_processed);
                axis image;
                title("Full Search Image");
                hold on;
                % Green box = no-shift
                rectangle('Position', [patch_col_range(1), patch_row_range(1), ...
                                       length(patch_col_range)-1, length(patch_row_range)-1], ...
                          'EdgeColor', 'g', 'LineWidth', 1.5);
                % Red box = actual search area
                rectangle('Position', [search_col_range(1), search_row_range(1), ...
                                       length(search_col_range)-1, length(search_row_range)-1], ...
                          'EdgeColor', 'r', 'LineWidth', 1.5);
                hold off;