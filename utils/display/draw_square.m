function draw_square(image,offsets,patch_size, dr)

    processed_image = process_volume(image,dr);
    figure();
    imagesc(processed_image);
    drawrectangle(gca,'Position',[offsets(2),offsets(1),patch_size(2),patch_size(1)], ...
    'FaceAlpha',0);

end