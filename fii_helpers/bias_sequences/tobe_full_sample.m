function sequence = tobe_full_sample(row_count)
    
    max_rx_per_tx = 1024; % Current cuda beamformer limitations

    total_count = row_count^2;

    num_transmits = ceil(total_count / max_rx_per_tx);

    sequence = cell(num_transmits,1);
    for T = 1:num_transmits
        
        full_array = zeros(row_count);

        start_i = (T-1) * max_rx_per_tx + 1;
        end_i = T * max_rx_per_tx;
        full_array(start_i:end_i) = 1;
    
        % Store each cross sequence
        sequence{T} = full_array;
    end

end

