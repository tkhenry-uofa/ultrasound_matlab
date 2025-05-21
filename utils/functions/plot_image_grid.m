function plot_image_grid(images, subplot_dims, x, z, title_str)

    figure();
    image_count = length(images);

    
    t = tiledlayout(subplot_dims(1), subplot_dims(2), 'TileSpacing', 'compact');
    for i = 1:image_count
        ax = nexttile;
        plot_bmode(images{i}, x, z);
        % axis off;
        xlabel("")
        ylabel("")
        title(sprintf("Image %d",i))
    end

    if nargin == 5
        % Main title
        title(t,title_str,'FontSize',16);
    end


    colorbar('Position', [0.93 0.15 0.02 0.7]);
end