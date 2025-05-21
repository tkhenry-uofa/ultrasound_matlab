function plot_bmode(image, x, z, title_str)

    % Assumes a figure is already open
    h = imagesc(image);
    colormap(gray);
    axis image;

    if nargin >= 3
        set(h, 'XData', x * 1000, 'YData', z * 1000); 
        xlabel('Lateral Position (mm)');
        ylabel('Depth (mm)');
    end

    if nargin == 4
        title(title_str);
    end
end