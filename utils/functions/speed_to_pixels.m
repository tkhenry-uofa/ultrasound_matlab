function pixels = speed_to_pixels(speed,prf,group_size,resolution)
pixels = speed .* group_size ./ (prf .* resolution);
end

