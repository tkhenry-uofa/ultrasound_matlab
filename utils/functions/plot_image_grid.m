function plot_image_grid(images,title, x_range, z_range, subplot_dims)

image_count = length(images);

figure();

sgtitle(title,'FontSize',16);

    % Loop through the images and display each in the 2x8 grid
    for i = 1:image_count
        subplot(subplot_dims(1), subplot_dims(2), i); 
        colormap("gray");
        imagesc(x_range, z_range, images{i});
        title(['Image ' num2str(i)]);
        xlabel("Lateral (mm)");
        ylabel("Depth (mm)");
        axis image;
    
    
    end

end