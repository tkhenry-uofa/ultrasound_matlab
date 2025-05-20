function [full_path,tx_dist,rx_dist] = calc_dist(vox_loc,th_loc,src_loc,type)

    rx_vec = th_loc - vox_loc.';
    
    x = vox_loc(1);
    y = vox_loc(2);
    z = vox_loc(3);

    xs = src_loc(1);
    ys = src_loc(2);
    zs = src_loc(3);

    %find location from each element to point
    rx_dist = vecnorm(rx_vec);

    src_sign = 1;

    if z < zs
        src_sign = -1;
    end

    if( type == TransmitType.Xdiv || type == TransmitType.Elevation_Focus || type == TransmitType.Focused )

        % Source is a line on the y axis, we only need its x and z pos
        tx_dist = src_sign * sqrt(( zs - z).^2 + ( xs - x ).^2 ) + zs;

    elseif type == TransmitType.Ydiv
        % Source is a line on the x axis, we only need its y and z pos
        tx_dist = src_sign * sqrt( (zs - z).^2 + ( ys - y).^2 ) + zs;
    else % type == "plane"
        tx_dist = vox_loc(3);
    end
    
    full_path = rx_dist + tx_dist;

end