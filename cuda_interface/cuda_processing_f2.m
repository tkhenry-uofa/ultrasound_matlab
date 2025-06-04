function low_res_array = cuda_processing_f2(rf_data,tx_config,vol_config,plane)
    
    data_size = size(rf_data);
    sample_count = data_size(1);
    channel_count = data_size(2);
    tx_count = data_size(3);

    % Beamformer expects time x transmits x channels
    rf_data_shuffled = permute(rf_data, [1 3 2]);

    %% Head properties
    bp.rf_raw_dim = [sample_count * tx_count, channel_count];
    bp.dec_data_dim = [data_size,1];

    bp.xdc_element_pitch = [tx_config.pitch, tx_config.pitch];
    
    bp.xdc_transform = eye(4);
    bp.xdc_transform(1,4) = tx_config.x_max;
    bp.xdc_transform(2,4) = tx_config.y_max;

    bp.time_offset = tx_config.pulse_delay;
    bp.das_shader_id = int32(tx_config.sequence);

    bp.transmit_mode = int32(acquisition.TransmitModes.RowTxColRx); % 1

    hadamard_sequences = [  SequenceType.FORCES, ...
                            SequenceType.UFORCES, ...
                            SequenceType.HERCULES];

    if(ismember(tx_config.sequence, hadamard_sequences))
        bp.decode = true;
    else
        bp.decode = false;
    end

    bp.data_type = int32(BeamformerDataType.F32);

    %% UI Properties

    bp.sampling_frequency = tx_config.fs;
    bp.center_frequency   = tx_config.f0;
    bp.speed_of_sound     = tx_config.c;
    bp.f_number           = vol_config.f_number;

    bp.output_min_coordinate = [vol_config.x_min, vol_config.y_min, vol_config.z_min, 0];
    bp.output_max_coordinate = [vol_config.x_max, vol_config.y_max, vol_config.z_max, 0];
    
    bp.output_points = [floor((vol_config.x_max - vol_config.x_min)/vol_config.lateral_resolution ), ...
                        floor((vol_config.y_max - vol_config.y_min)/vol_config.lateral_resolution ), ...
                        floor((vol_config.z_max - vol_config.z_min)/vol_config.axial_resolution ), ...
                        1]; % Number of frames for averaging

    if(plane)
        bp.output_points(2) = 1;
        bp.output_min_coordinate(2) = 0;
        bp.output_max_coordinate(2) = 0;
    end
   
    %% Readi
    
    if(isfield(vol_config, "readi_group_count"))
        readi_group_count = vol_config.readi_group_count;
    else
        readi_group_count = 1;
    end

    if readi_group_count < 1
        readi_group_count = 1;
    end

    data_group_cells = cell(1,readi_group_count);
    readi_group_size = tx_count/readi_group_count;

    % Break up transmits into Readi groups
    for i=1:readi_group_count
        end_tx = i * readi_group_size;
        start_tx = end_tx - (readi_group_size - 1);
        data_group_cells{i} = rf_data_shuffled(:,start_tx:end_tx,:);
    end
    bp.rf_raw_dim(1) = bp.rf_raw_dim(1) / readi_group_count;
    bp.dec_data_dim(3) = readi_group_size;
    tx_count = readi_group_size;

    bp.readi_group_count = readi_group_count;
    bp.readi_group_id = 0;


    %% Arrays
    bp.channel_mapping  = 0:(channel_count-1);
    bp.focal_depths  = ones(1,tx_count).*tx_config.src(3);
    bp.transmit_angles = zeros(1,channel_count);
    bp.sparse_elements = ones(1,tx_count) * -1;

    %% Beamforming
    low_res_array = cell(1,readi_group_count);
    for g = 1:readi_group_count
        fprintf("Beamforming Image %d\n", g );
        raw_image = cuda_beamform(data_group_cells{g}, bp);
        fprintf("Done\n");
        low_res_array{g} = raw_image;
        bp.readi_group_id = bp.readi_group_id + 1;
    end
    
    unloadlibrary('cuda_transfer_matlab')

end

