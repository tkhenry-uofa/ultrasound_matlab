% Parameters
patch_size = 1;                             % Grid spacing in pixels

image_id = 1;

figure();

% subplot 121
img = generate_walsh(32);
imagesc(img); axis("image"); hold on;


% Get image dimensions
[rows, cols, ~] = size(img);

% Overlay grid and index labels
for row = 1:patch_size:rows
    for col = 1:patch_size:cols
        % Draw grid lines
        if col == 1
            line([1 cols], [row row], 'Color', 'r', 'LineWidth', 0.5);  % horizontal
        end
        if row == 1
            line([col col], [1 rows], 'Color', 'r', 'LineWidth', 0.5);  % vertical
        end
        
        % % Compute patch indices
        % row_idx = ceil(row / patch_size);
        % col_idx = ceil(col / patch_size);
        % 
        % % Overlay index text
        % text(col + 2, row + 12, ...
        %      sprintf('(%d,%d)', row_idx, col_idx), ...
        %      'Color', 'y', 'FontSize', 8, 'FontWeight', 'bold');
    end
end

title(['Grid overlay with spacing = ' num2str(patch_size) ' px and patch indices']);
% 
% subplot 122
% imagesc(img); axis("image"); 
