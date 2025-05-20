function speed = pixels_to_speed(pixels,prf,group_size,resolution)
    speed = pixels .* prf .* resolution ./ group_size;
end

