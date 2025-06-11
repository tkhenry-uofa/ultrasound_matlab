
function filtered_images_cell = svd_filter(images_cell, tissue_sv_cutoff, noise_sv_cutoff, tissue_start)


    if ~exist("tissue_start", "var")
        tissue_start = 1;
    end

    cell_size = size(images_cell);

    images_cell = reshape(images_cell,1,[]);
    image_count = numel(images_cell);
    image_size = size(images_cell{1});
    image_array_size = [image_size, image_count];
    image_array = zeros(image_array_size);
    
    for i=1:image_count
        image_array(:,:,i) = images_cell{i};
    end
    
    flat_images = reshape(image_array,[],image_count);
    
    [U,S,V] = svd(flat_images,"econ");

    S_filter = S;

    if(tissue_sv_cutoff > 0)
        S_filter(tissue_start:tissue_sv_cutoff,:) = 0;
    end

    if(noise_sv_cutoff > 0)
        S_filter(noise_sv_cutoff:end,:) = 0;
    end
    
    filtered_images_array = single(U*S_filter*V');
    filtered_images_array = reshape(filtered_images_array,image_array_size);
    
    filtered_images_cell = mat2cell(filtered_images_array, image_size(1), image_size(2), ones(1,image_count));
    filtered_images_cell = reshape(filtered_images_cell, cell_size);
end