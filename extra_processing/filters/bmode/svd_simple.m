
function filtered_images_cell = svd_simple(images_cell)

image_count = numel(images_cell);
image_size = size(images_cell{1});
image_array_size = [image_size, image_count];
image_array = zeros(image_array_size);

for i=1:image_count
    image_array(:,:,frame) = images_cell{i};
end

flat_images = reshape(image_array,[],image_count);

[U,S,V] = svd(flat_images,"econ");
%%
S_filter = S;

% 10 ml
tissue_svd_cutoff = 17;
noise_svd_cutoff = 100;

% 60 ml
% tissue_svd_cutoff = 20;
% noise_svd_cutoff = 80;

S_filter(1:tissue_svd_cutoff,:) = 0;
S_filter(noise_svd_cutoff:end,:) = 0;

filtered_images_array = U*S_filter*V';
filtered_images_array = reshape(filtered_images_array,frame_size);

filtered_images_cell = mat2cell(filtered_images_array, image_size(1), image_size(2), ones(1,image_count));

end