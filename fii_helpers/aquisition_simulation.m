function [scan_lines, rx_element_locs, signal_length] = aquisition_simulation(config, active_rx_els, points, amps, parallel, transmit_no)

    %% Field II need to be reinitialized on each thread
    if parallel == true
        field_init(-1)
        set_field('fs', config.fs);
        set_field('c',config.c);
    end
    
    element_count = config.rows * config.cols;
    tx_subs = config.sub_element_count;

    tx_array = create_tx_array(config, transmit_no);
    [rx_array, rx_subs] = create_rx_array(config, transmit_no,active_rx_els);
    
    % show_xdc_mod(rx_array,"apo")
    array_info = xdc_get(rx_array);
    rx_element_locs_full = array_info(24:26,1:rx_subs:end);
    element_count = length(rx_element_locs_full);

    [scans, rx_start] = calc_scat_multi( tx_array, rx_array, points, amps );


    % Add zeros for the time before the first reception
    zero_padding = zeros(length(0:1/config.fs:rx_start), element_count);

    all_scans = [zero_padding; scans];

    % Scale up to avoid float underflow
    % all_scans = all_scans * 1e25;

    signal_length = size(all_scans,1);

    % We need to sum up each column, (rows are biased)
    % calc_scat output is sequential rows back to front
    if config.sequence == SequenceType.HERCULES
        
        scan_lines = zeros(signal_length,config.rows);
        for i=1:config.cols
            start_channel = i;
            end_channel = config.cols * (config.rows - 1) + i;
            scan_lines(:,i) = sum(all_scans(:, start_channel:config.cols:end_channel), 2);
        end

        % When decoded these will be the elements of the T'th row (back to front)

        starting_loc = (transmit_no-1) * config.cols + 1;
        ending_loc = transmit_no * config.cols;

        rx_element_locs = rx_element_locs_full(:, starting_loc:ending_loc);
    else
        scan_lines = all_scans;
        rx_element_locs = rx_element_locs_full;
    end

    scan_lines = single(scan_lines);

end