function [compound_image,readi_images] = process_readi(raw_images,dr,tr,power)

    image_count = numel(raw_images);
    readi_images = cell(size(raw_images));

    compound_image_raw = zeros(size(raw_images{1}));
    
    for i=1:image_count
        raw_image = raw_images{i};
        compound_image_raw = compound_image_raw + raw_image;
        readi_images{i} = process_volume(raw_image,dr,tr,power);
    end
    
    compound_image = process_volume(compound_image_raw,dr,tr,power);
end

