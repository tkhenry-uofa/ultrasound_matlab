function [bp,x_range, y_range, z_range] = configure_output_points(ranges,lateral_resolution, axial_resolution, bp)

    bp.output_min_coordinate(1:3) = ranges(:,1);
    bp.output_max_coordinate(1:3) = ranges(:,2);
    
    
    resolutions = [lateral_resolution, lateral_resolution, axial_resolution, 1];
    
    bp.output_points = floor( (bp.output_max_coordinate - bp.output_min_coordinate) ./ resolutions);
    bp.output_points( bp.output_points == 0 ) = 1;
    
    x_range = linspace(ranges(1,1), ranges(1,2), bp.output_points(1));
    y_range = linspace(ranges(2,1), ranges(2,2), bp.output_points(2));
    z_range = linspace(ranges(3,1), ranges(3,2), bp.output_points(3));
    

end