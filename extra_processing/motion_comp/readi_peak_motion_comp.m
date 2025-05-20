function shifted_image_raw = readi_peak_motion_comp(low_res_array,readi_group_count,speed,resolution)

    if(readi_group_count <= 1)
        shifted_image_raw = low_res_array{1};
        return;
    end


    middle_image = floor(readi_group_count/2);
    volume = low_res_array{middle_image};

    image_2d = false;
    if(length(size(volume)) == 2)
        image_2d = true;

        % Make it 3d so we don't need two version of the code
        volume = reshape(volume,size(volume,1),size(volume,2),1);
    end
    
    [~, max_idx] = max(volume(:));
    left_max = [0,0,0];
    [left_max(1),left_max(2),left_max(3)] = ind2sub(size(volume), max_idx);
    
    volume = low_res_array{middle_image+1};

    if image_2d
        volume = reshape(volume,size(volume,1),size(volume,2),1);
    end
    
    [~, max_idx] = max(volume(:));
    right_max = left_max;
    [right_max(1),right_max(2),right_max(3)] = ind2sub(size(volume), max_idx);
    
    avg_max = round((right_max + left_max)/2);
    
    img_size = size(volume);
    
    shifted_image_raw = zeros(img_size);
    frame_max = zeros(size(avg_max));
    for i = 1:readi_group_count
    
        volume = low_res_array{i};

        if image_2d
            volume = reshape(volume,size(volume,1),size(volume,2),1);
        end
    
        [~, max_idx] = max(volume(:));
        [frame_max(1),frame_max(2),frame_max(3)] = ind2sub(size(volume), max_idx);
    
        diff = avg_max - frame_max;
    
        % diff = speed * (i-1)/resolution;
    
    
        x_src = round(max(1, 1-diff(1)):min(img_size(1), img_size(1)-diff(1)));
        y_src = round(max(1, 1-diff(2)):min(img_size(2), img_size(2)-diff(2)));
        z_src = round(max(1, 1-diff(3)):min(img_size(3), img_size(3)-diff(3)));
        x_dest = round(x_src + diff(1));
        y_dest = round(y_src + diff(2));
        z_dest = round(z_src + diff(3));
    
        shifted_image_raw(x_dest,y_dest,z_dest) = shifted_image_raw(x_dest,y_dest,z_dest) + volume(x_src,y_src,z_src);
    
    end

    if image_2d
        shifted_image_raw = reshape(shifted_image_raw,size(shifted_image_raw,1),size(shifted_image_raw,2));
    end

end

