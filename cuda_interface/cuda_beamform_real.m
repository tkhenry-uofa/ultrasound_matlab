function [volume] = cuda_beamform_real(data_array,bp)


    pipe_name = '\\.\pipe\beamformer_data_fifo';
    smem_name = 'Local\ogl_beamformer_parameters';


    output_counts_xyz.x = bp.output_points(1);
    output_counts_xyz.y = bp.output_points(2);
    output_counts_xyz.z = bp.output_points(3);
    output_counts_xyz.w = 1;
    
    
    % Complex volumes aren't supported so they're interleaved
    interleaved_volume_size = [bp.output_points(1)*2, bp.output_points(2), bp.output_points(3)];
    
    if ~libisloaded('cuda_transfer'), loadlibrary('cuda_transfer'); end
      
    calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

    
    volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));
    
    % TODO: Fix in API
    rf_raw_dim_struct.x = bp.rf_raw_dim(1);
    rf_raw_dim_struct.y = bp.rf_raw_dim(2);

    [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
        pipe_name, smem_name, data_array, rf_raw_dim_struct, output_counts_xyz, volume_ptr);
    
    volume = reshape(volume_ptr,interleaved_volume_size );
    
    real_vol = volume(1:2:end,:,:);
    im_vol = volume(2:2:end,:,:);
    
    volume = squeeze(complex(real_vol, im_vol));

end