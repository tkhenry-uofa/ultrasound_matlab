
frame_count_real = frame_count;
% frame_count = 12;
image_size = size(low_res_array{1,1});

image_count = frame_count * readi_group_count;
frame_size = [image_size, image_count];

image_array = zeros(frame_size);

frame = 1;
for f=1:frame_count
    for i=1:readi_group_count
        image_array(:,:,frame) = low_res_array{f,i};
        frame = frame + 1;
    end
end

flat_frame = reshape(image_array,[],image_count);

[U,S,V] = svd(flat_frame,"econ");
%%
S_filter = S;

% 10 ml
% tissue_svd_cutoff = 17;
% noise_svd_cutoff = 120;

% 60 ml
% tissue_svd_cutoff = 20;
% noise_svd_cutoff = 80;
% 

tissue_svd_cutoff = 50;
noise_svd_cutoff = 4;
% 
% S_filter(1:tissue_svd_cutoff,:) = 0;
S_filter(noise_svd_cutoff:end,:) = 0;

filtered_flat_frame = U*S_filter*V';
% filtered_flat_frame = flat_frame;

filtered_frame_array = reshape(filtered_flat_frame,frame_size);

filtered_frame_cell = mat2cell(filtered_frame_array, image_size(1), image_size(2), ones(1,image_count));
filtered_frame_cell = reshape(filtered_frame_cell, readi_group_count, frame_count)';

processed_frame_array = uint8(zeros(size(filtered_frame_array)));
filtered_image = zeros(image_size);
filtered_forces_images = uint8(zeros([image_size, frame_count]));

dr = 45; power_thr = 2000;

for f=1:frame_count
    filtered_forces_image = zeros(image_size);
    for g = 1:readi_group_count
            
        i = (f-1) * readi_group_count + g;
        current_image = squeeze(filtered_frame_array(:,:,i));
        filtered_image = filtered_image + current_image;
        filtered_forces_image = filtered_forces_image + current_image; 

        current_image = uint8((process_volume(current_image,dr,power_thr,false) + dr) *255/dr);
        processed_frame_array(:,:,i) = current_image;
    end
    
    filtered_forces_images(:,:,f) = uint8((process_volume(filtered_forces_image,dr) + dr) *255/dr);
end
%%

% processed_frame_array = uint8(process_volume(filtered_frame_array,dr,power_thr,true)*255);
% processed_image = process_volume(filtered_image,40);

% figure();imagesc(processed_image);colormap gray; axis image;

% view_frame = 10;
% frame_range = (1:readi_group_count) + (view_frame - 1) * readi_group_count; 


% implay(processed_frame_array(:,:,frame_range),2);

% implay(readi_compare_array,4);
% 


implay(processed_frame_array,4*8);

% implay(filtered_forces_images,2);

%%
frame_count = frame_count_real;