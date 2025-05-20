function aprodizations = apro_no_fnum(vox_loc,element_locs,tx_config)

    element_count = length(element_locs);
    l = tx_config.l_angle;
    r = tx_config.r_angle;
    type = tx_config.transmit;
    x_diam = 2*tx_config.x_max; % All coordinates assume the center of the array is 0,0,0
    y_diam = 2*tx_config.y_max; % so these should be symmetric with the minimums
    angle = 0;
    % Find the widest angle that the wave leaves the apature at
    if tx_config.l_angle > tx_config.r_angle
       angle = tx_config.l_angle;
    else
       angle = tx_config.r_angle;
    end

    depth = vox_loc(3);
    
    N = 128;
    hann_window = hann(2*N);
    hann_window = hann_window(N+1:end);

    if type == "xLine"
        max_dist = y_diam + depth*tan(angle);
        y_dists = abs(element_locs(2,:) - vox_loc(2));
        indecies = round(y_dists/max_dist*N);
    elseif type == "yLine"
        max_dist = x_diam + depth*tan(angle);
        x_dists = abs(element_locs(1,:) - vox_loc(1));
        indecies = round(x_dists/max_dist*N);
    else
        max_dist = sqrt((x_diam/2)^2 + (y_diam/2)^2);
        dists = sqrt( (element_locs(1,:) - vox_loc(1)).^2 + (element_locs(2,:) - vox_loc(2)).^2);
        indecies = round(dists/max_dist*N);
    end
    
    indecies = indecies + 1; % One based indexing is evil
    indecies( indecies > N ) = N;
    aprodizations = hann_window(indecies);
    
end