function [points, amps] = cyst2x2(x_range, y_range, z_range, samples, radius)

    [points1,amps1] = single_cyst(samples,x_range(1:2),y_range,z_range(1:2),radius);
    [points2,amps2] = single_cyst(samples,x_range(2:3),y_range,z_range(1:2),radius);
    [points3,amps3] = single_cyst(samples,x_range(1:2),y_range,z_range(2:3),radius);
    [points4,amps4] = single_cyst(samples,x_range(2:3),y_range,z_range(2:3),radius);
    
    points = [points1;points2;points3;points4];
    amps = [amps1;amps2;amps3;amps4];
        
end