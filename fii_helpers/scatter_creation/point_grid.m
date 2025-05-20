function [points, amps] = point_grid( x_points, y_points, z_points )


    p_x_len = length(x_points);
    p_y_len = length(y_points);
    p_z_len = length(z_points);
    
    pointArray = zeros(p_x_len,p_y_len,p_z_len,3);
    for i = 1:p_x_len
        for j = 1:p_y_len
            for k = 1:p_z_len
                pointArray(i,j,k,:) = [x_points(i) y_points(j) z_points(k)];
            end
        end
    end
    
    points = reshape(pointArray, [], 3);
    amps = ones(p_z_len*p_y_len*p_x_len,1);

end