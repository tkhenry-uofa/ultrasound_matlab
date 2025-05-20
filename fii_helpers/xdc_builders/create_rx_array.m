function [rx_array, rx_subs] = create_rx_array(config, transmit_no,active_rx_els)

    encoded_sequences = [SequenceType.Hercules, SequenceType.Hercforce];
    linear_sequences = [SequenceType.Forces, SequenceType.Readi];

    rx_subs = config.sub_element_count;
    if(ismember(config.sequence, linear_sequences))
    
        rx_subs = 32;
        % Just make a long 1D array because all the rows are active
        % anyway
        rx_array = xdc_linear_array(config.cols, config.width, config.pitch * config.rows, ...
            config.kerf, 1, rx_subs, [0 0 10000000000000]);
    else
        rx_array = xdc_2d_array(config.cols,config.rows,config.width,config.width,...
            config.kerf,config.kerf,active_rx_els,1,1,[0 0 100000000000]);
    end

    % Tx on rows (Y), rx on columns (X), bias pattern is across rows 
    if (ismember(config.sequence, encoded_sequences))
        H = hadamard(config.rows);

 
        apo_pattern = repmat(squeeze(H(transmit_no,:)),config.rows,1);
        apo_pattern = reshape(apo_pattern, [], 1);

        ele_apodization(rx_array,(1:(config.rows*config.cols)).',apo_pattern);
        
    end

    xdc_impulse(rx_array,config.imp);
end