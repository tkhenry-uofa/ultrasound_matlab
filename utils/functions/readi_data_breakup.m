function [readi_group_data,readi_bp] = readi_data_breakup(frame_data,bp,readi_group_count)
    frame_count = numel(frame_data);
    readi_bp = bp;
    
    readi_bp.readi_group_count = readi_group_count;
    readi_group_data = cell(frame_count,readi_group_count);
    readi_group_size = bp.dec_data_dim(3) / readi_group_count;
    
    readi_bp.readi_group_id = 0;
    
    % Break up transmits into Readi groups
    for f=1:frame_count
        full_frame_data = frame_data{f};
        for i=1:readi_group_count
            end_tx = i * readi_group_size;
            start_tx = end_tx - (readi_group_size - 1);
            readi_group_data{f,i} = full_frame_data(:,start_tx:end_tx,:);
        end
    end
    readi_bp.rf_raw_dim(1) = readi_group_size * bp.dec_data_dim(1);
    readi_bp.dec_data_dim(3) = readi_group_size;
end

