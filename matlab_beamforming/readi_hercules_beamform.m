function low_res_array= readi_hercules_beamform(tx_config,vol_config,rx_scans,rx_element_locs,xr_range,yr_range,zr_range, readi_group_count, f_number)
    no_transmits = tx_config.no_transmits;
    loop_average = 0;
    
    H_128 = hadamard(128);

    group_count = 16;
    group_size = 8;

    low_res_array = cell(group_count,1);

    src_loc = tx_config.src;

    data_cell = cell(1,group_count);

    if (f_number == 0)
        apodize = false;
    else
        apodize = true;
    end

    for i = 1:group_count
        data_cell{i} = rx_scans(:,:,(i-1)*8+1 : i*8);
    end

    % Data from the 128 virtual tx elements is spread evenly across 
    % the 128 events, we are breaking it into groups of 8 events and
    % pulling data from all 128 virtual sources out of those 8 events to
    % make a volume
    parfor event_group = 1:group_count
        fprintf("Group %i\n",event_group);
        low_res_image = single(zeros(vol_config.x_count, vol_config.y_count, vol_config.z_count));

        group_data = data_cell{event_group};
        data_size = size(group_data);

        % Stack each rx element in the first dimension so we can
        % 2d matrix multiply 
        group_data = reshape(group_data, [], group_size);
        
        for source_group = 1:group_count
            fprintf("Event group: %i, Source group: %i\n", event_group, source_group);
            
            hadamard_group = H_128((event_group-1)*8+1:(event_group*8),(source_group-1)*8+1:(source_group*8));
            
            decoded_data = group_data * hadamard_group;
            data_size(3) = 8;
            decoded_data = reshape(decoded_data, data_size);
            decoded_data = hilbert(decoded_data);

            for v = 1:group_size
                src_idx = (source_group -1) * group_size + v;
                src_data = decoded_data(:,:,v);
                rx_locs = rx_element_locs{src_idx};

                for i = xr_range
                    x_loc = vol_config.x_range(i);

                    for j = yr_range
                        y_loc = vol_config.y_range(j);
                            
                        for k = zr_range
                            z_loc = vol_config.z_range(k);
        
                            vox_pos = [x_loc, y_loc, z_loc];
        
                            [full_dist, ~, ~] = calc_dist(vox_pos,rx_locs,src_loc,tx_config.transmit); 
                            row_floats = (full_dist/tx_config.c + tx_config.pulse_delay)*tx_config.fs;
                            row_floats(row_floats > size(src_data,1)) = 1; % The first few samples will always be 0 because of the padding
                            row_indices = round(row_floats);
    
                            indices = sub2ind(size(src_data), row_indices, 1:length(rx_locs))';
                            values = src_data(indices);
                            
                            if apodize == true
        
                                apos = hann_apo([x_loc, y_loc, z_loc],rx_locs,tx_config, f_number);
                                values = values.*apos;
                            end
                            
                            low_res_image(i,j,k) = low_res_image(i,j,k) + single(sum(values));
                    
                        end
                    end
                end
            end
        end

        low_res_array{event_group} = low_res_image;
        % high_res_image = high_res_image + low_res_image;
    end
end