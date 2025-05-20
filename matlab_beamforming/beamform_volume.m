function volume = beamform_volume(tx_config,vol_config,rx_scans,rx_element_locs,xr_range,yr_range,zr_range, parallel, f_number)
    no_transmits = tx_config.no_transmits;
    loop_average = 0;
    volume_array = cell(no_transmits,1);

    apodize = f_number > 0;

    parfor T = 1:no_transmits
        transmit_volume = single(zeros(vol_config.size));
        element_locs = rx_element_locs{T};
        element_count = length(element_locs);
        scans = rx_scans{T};
        fprintf("Transmit: %d \n",T);
        
        source = tx_config.src;

        if tx_config.sequence == SequenceType.Forces || tx_config.sequence == SequenceType.Readi
            if ismember(tx_config.discard_txs, T)
                continue;
            end
            source = element_locs(:,tx_config.forces_srcs(T));

        end

        for i = xr_range
            tic;
            if parallel == false
                remaining = loop_average * (length(xr_range) + 1 - i );
                min_remaining = floor(remaining/60);
                sec_remaining = mod(remaining, 60);
                percent_remaining = (i-1)*100/length(xr_range);
                fprintf("Slice %d/%d (%.2f%%)    %dm %.2fs remaining (%.2f s/slice)\n", i, length(xr_range), percent_remaining, min_remaining, sec_remaining, loop_average);
            end
            x_loc = vol_config.x_range(i);
        
            for j = yr_range
                % fprintf("Slice: %d \n", j)
                y_loc = vol_config.y_range(j);
        
                for k = zr_range
                    z_loc = vol_config.z_range(k);

                    vox_pos = [x_loc, y_loc, z_loc];

                    [full_dist, ~, ~] = calc_dist(vox_pos,element_locs,source,tx_config.transmit); 
                    row_floats = (full_dist/tx_config.c + tx_config.pulse_delay)*tx_config.fs;
                    row_floats(row_floats > size(scans,1)) = 1; % The first few samples will always be 0 because of the padding
                    row_indices = round(row_floats);

                    indices = sub2ind(size(scans), row_indices, 1:length(element_locs))';

                    values = scans(indices);
                    
                    if apodize == true
                        apos = hann_apo([x_loc, y_loc, z_loc],element_locs,tx_config,f_number);
                        values = values.*apos;
                    end
                    
                    transmit_volume(i,j,k) = single(sum(values));
        
                end
            end
            if parallel == false
                loop_average = (2*loop_average + toc)/3;
            end
        end

        volume_array{T} = transmit_volume;
    end
    
    volume = single(zeros(vol_config.size));
    for i = 1:no_transmits
        volume = volume + volume_array{i};
    end

end