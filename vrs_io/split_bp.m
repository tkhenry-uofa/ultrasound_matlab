function [head, ui, tail, arrays] = split_bp(bp)

    % --- Head ---
    head.rf_raw_dim          = bp.rf_raw_dim;
    head.dec_data_dim        = bp.dec_data_dim;
    head.xdc_element_pitch   = bp.xdc_element_pitch;
    head.xdc_transform       = bp.xdc_transform;
    head.transmit_mode       = bp.transmit_mode;
    head.decode              = bp.decode;
    head.das_shader_id       = bp.das_shader_id;
    head.time_offset         = bp.time_offset;

    % --- UI ---
    ui.beamform_plane        = bp.beamform_plane;
    ui.speed_of_sound        = bp.speed_of_sound;
    ui.center_frequency      = bp.center_frequency;
    ui.sampling_frequency    = bp.sampling_frequency;
    ui.interpolate           = bp.interpolate;
    ui.output_min_coordinate = bp.output_min_coordinate;
    ui.output_max_coordinate = bp.output_max_coordinate;
    ui.output_points         = bp.output_points;
    ui.off_axis_pos          = bp.off_axis_pos;
    ui.f_number              = bp.f_number;

    % --- Tail ---
    tail.readi_group_size   = 0;
    tail.readi_group_id     = bp.readi_group_id;

    % --- Arrays ---
    sparse_elements  = bp.sparse_elements;
    transmit_count = head.dec_data_dim(3);
    if(sparse_elements(1) == -1)
        sparse_elements(1:transmit_count) = 1:transmit_count;
    end

    focal_depths     = bp.focal_depths;
    transmit_angles  = bp.transmit_angles;

    focal_vectors = zeros(1, length(focal_depths) * 2);
    focal_vectors(1:2:end) = transmit_angles;
    focal_vectors(2:2:end) = focal_depths;
    
    arrays.focal_vectors    = focal_vectors;
    arrays.sparse_elements  = sparse_elements;
    arrays.channel_mapping  = bp.channel_mapping;
end