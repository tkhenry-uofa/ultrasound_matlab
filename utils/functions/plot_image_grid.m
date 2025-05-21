function plot_image_grid(images,x, z, subplot_dims, title_str)

    figure();
    image_count = length(images);

    if nargin == 5
        % Main title
        sgtitle(title_str,'FontSize',16);
    end

    % Loop through the images and display each in the 2x8 grid
    for i = 1:image_count
        subplot(subplot_dims(1), subplot_dims(2), i); 
    
        plot_bmode(images{i}, x, z, sprintf('Frame %d', i));
        xlabel('mm');
        ylabel('mm');
    end

    colorbar('Position', [0.92 0.15 0.02 0.7]);
end