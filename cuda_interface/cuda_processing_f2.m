function low_res_array = cuda_processing_f2(rf_data,tx_config,readi_group_count,vol_config, f_number, sequence_type)

    addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")

    pipe_name = '\\.\pipe\beamformer_data_fifo';
    smem_name = 'Local\ogl_beamformer_parameters';
    pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn
    
    
    bp.output_min_coordinate = struct('x', vol_config.x_min, 'y', vol_config.y_min, 'z', vol_config.z_min,   'w', 0);
    bp.output_max_coordinate = struct('x',  vol_config.x_max, 'y', vol_config.y_max, 'z', vol_config.z_max, 'w', 0);
    
    bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/vol_config.lateral_resolution );
    bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/vol_config.lateral_resolution );
    bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/vol_config.axial_resolution );
    bp.output_points.w = 1; % Number of frames for averaging
    
    bp.das_shader_id = sequence_type;

    % bp.das_shader_id = 0; % Forces
    % bp.das_shader_id = 2; % Hercules
    if(sequence_type == 0)
        bp.output_min_coordinate.y = 0;
        bp.output_max_coordinate.y = 0;
        bp.output_points.y = 1;
    end

    bp.rf_raw_dim = struct('x', size(rf_data, 1) * size(rf_data, 3), 'y', size(rf_data, 2));
    
    bp.dec_data_dim.x = size(rf_data,1);
    bp.dec_data_dim.y = size(rf_data,2);
    bp.dec_data_dim.z = size(rf_data,3);
    bp.dec_data_dim.w = 0; % Averaging
    
    bp.sampling_frequency = tx_config.fs;
    bp.center_frequency   = tx_config.f0;
    bp.speed_of_sound     = tx_config.c;
    
    bp.time_offset = tx_config.pulse_delay;
    
    bp.channel_mapping  = ones(1,256);

    bp.xdc_element_pitch = [tx_config.pitch, tx_config.pitch];
    
    bp.xdc_transform = eye(4);
    bp.xdc_transform(1,4) = tx_config.x_max;
    bp.xdc_transform(2,4) = tx_config.y_max;
    
    bp.focal_depths  = ones(1,size(rf_data,3)).*tx_config.src(3);
    
    % NOTE: plane and position along plane normal for beamforming 2D HERCULES
    bp.beamform_plane = 0;
    bp.off_axis_pos = 0;
    fc = bp.center_frequency;
    fs = bp.sampling_frequency;
    
    bp.f_number = 1;
    
    data_group_cells = cell(1,16);
   
    bp.readi_group_size = bp.dec_data_dim.z/readi_group_count;
    bp.readi_group_id = 0;

    bp.data_type = 1; % 0 = int16, 1 = float32

    bp.f_number = f_number;
    
    % Break up transmits into Readi groups
    for i=1:readi_group_count
        end_tx = i * bp.readi_group_size;
        start_tx = end_tx - (bp.readi_group_size - 1);
    
        data = rf_data(:,:,start_tx:end_tx);
        data_group_cells{i} = reshape(data, [], bp.readi_group_size);
    end
    
    bp.rf_raw_dim.x = bp.rf_raw_dim.x / readi_group_count;
    
    if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end
    
    loadlibrary('cuda_transfer')
    
    fprintf("Setting params\n")
    calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

    fprintf("Sending data\n")
    
    
    %% Readi loop
    
    output_counts_xyz.x = bp.output_points.x;
    output_counts_xyz.y = bp.output_points.y;
    output_counts_xyz.z = bp.output_points.z;
    
    low_res_array = cell(1,readi_group_count);
    
    % Complex volumes aren't supported so they're interleaved
    interleaved_volume_size = [bp.output_points.x*2, bp.output_points.y, bp.output_points.z];
    
    
    for g = 1:readi_group_count
    
        calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);
        fprintf("Readi group %i. \n", g);
    
        volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));
    
        try
            [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_f32', ...
            pipe_name, smem_name, data_group_cells{g}, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
            
        catch Error
            unloadlibrary('cuda_transfer');
            rethrow(Error);
        end
        
        
        fprintf("Received Response\n")
        volume = reshape(volume_ptr,interleaved_volume_size );
    
        real_vol = volume(1:2:end,:,:);
        im_vol = volume(2:2:end,:,:);
        
        low_res_array{g} = squeeze(complex(real_vol, im_vol));
    
        % pause(0.5);
    
        bp.readi_group_id = bp.readi_group_id + 1;
       
    end
    
    unloadlibrary('cuda_transfer')

end

