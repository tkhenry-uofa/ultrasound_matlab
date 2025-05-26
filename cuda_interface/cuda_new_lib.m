function [volume] = cuda_new_lib(data_array,bp)
   
    % Complex volumes aren't supported so they're interleaved
    interleaved_volume_size = [bp.output_points(1)*2, bp.output_points(2), bp.output_points(3)];
    
    if ~libisloaded('cuda_transfer_new_matlab'), loadlibrary('cuda_transfer_new_matlab'); end
   
    volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));
    
    [~,volume_ptr] = calllib('cuda_transfer_new_matlab', 'beamform_i16', data_array, bp, volume_ptr);
    

    volume = reshape(volume_ptr,interleaved_volume_size );
    
    real_vol = volume(1:2:end,:,:);
    im_vol = volume(2:2:end,:,:);
    
    volume = squeeze(complex(real_vol, im_vol));

end