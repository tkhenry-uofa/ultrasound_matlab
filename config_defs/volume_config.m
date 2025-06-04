function config = volume_config(x_range,y_range,z_range,resolution)

    config = struct();

    config.x_min = x_range(1);
    config.x_max = x_range(2); 
    
    config.y_min = y_range(1);
    config.y_max = y_range(2); 
    
    config.z_min = z_range(1);
    config.z_max = z_range(2); 

    config.lateral_resolution = resolution;
    config.axial_resolution = resolution;

    config.x_count = floor((x_range(2)-x_range(1))/config.lateral_resolution);
    config.y_count = floor((y_range(2)-y_range(1))/config.lateral_resolution);
    config.z_count = floor((z_range(2)-z_range(1))/config.axial_resolution);

    if config.x_count == 0
        config.x_count = 1;
    end
    
    if config.y_count == 0
        config.y_count = 1;
    end
    
    if config.z_count == 0
        config.z_count = 1;
    end

    config.x_range = config.x_min:config.lateral_resolution:config.x_max;
    config.y_range = config.y_min:config.lateral_resolution:config.y_max;
    config.z_range = config.z_min:config.axial_resolution:config.z_max;

    % GPU beamforming will count like this even if it means there is no
    % cell at x_max
    if(length(config.x_range) > config.x_count)
        config.x_range = config.x_range(1:config.x_count);
    end

    if(length(config.y_range) > config.y_count)
        config.y_range = config.y_range(1:config.y_count);
    end

    if(length(config.z_range) > config.z_count)
        config.z_range = config.z_range(1:config.z_count);
    end

    config.x_mid = round(config.x_count/2);
    config.y_mid = round(config.y_count/2);
    config.z_mid = round(config.z_count/2);
    
    config.size = [config.x_count,config.y_count,config.z_count];

    config.f_number = 0;
    config.readi_group_count = 0;

end