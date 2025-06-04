function bp = load_and_parse_bp(path)

    raw_bp = ornot_bp_load_mex(convertStringsToChars(path));


    %% Head Properties
                            
    bp.rf_raw_dim               = raw_bp.raw_data_dim(1:2);
    bp.dec_data_dim             = raw_bp.decoded_data_dim;

    bp.xdc_element_pitch        = raw_bp.transducer_element_pitch;
    bp.xdc_transform            = raw_bp.transducer_transform_matrix;
    
    bp.transmit_mode            = raw_bp.transmit_mode;
    bp.decode                   = raw_bp.decode_mode;
    bp.das_shader_id            = raw_bp.beamform_mode;
    bp.time_offset              = raw_bp.time_offset;

    %% UI Properties

    bp.speed_of_sound           = raw_bp.speed_of_sound;
    bp.center_frequency         = raw_bp.center_frequency;
    bp.sampling_frequency       = raw_bp.sampling_frequency;

    % bp.interpolate              = false;
    % bp.coherency_weighting      = false;

    bp.output_min_coordinate    = [0, 0, 0, 0];
    bp.output_max_coordinate    = [0, 0, 0, 0];
    bp.output_points            = [0, 0, 0, 0];
    bp.off_axis_pos             = 0;
    bp.beamform_plane           = 0;
    bp.f_number                 = 0;
    
    

    %% Tail Properties    
    bp.readi_group_count = 1;
    % bp.readi_group_size = 128;
    bp.readi_group_id = 0;

    bp.data_type = int32(BeamformerDataType.I16);

    %% Arrays
    bp.channel_mapping  = raw_bp.channel_mapping;
    bp.transmit_angles  = raw_bp.steering_angles;
    bp.focal_depths     = raw_bp.focal_depths;
    bp.sparse_elements  = raw_bp.sparse_elements;

end