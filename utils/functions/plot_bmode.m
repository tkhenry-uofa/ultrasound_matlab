function plot_bmode(image, x, z, title_str)
    % Assumes a figure is already open
    imagesc(x * 1000, z * 1000, image);
    colormap(gray);
    axis image;
    xlabel('Lateral Position (mm)');
    ylabel('Depth (mm)');

    if nargin == 4
        title(title_str);
    end
end