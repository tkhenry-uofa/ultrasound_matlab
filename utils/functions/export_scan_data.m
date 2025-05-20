function export_scan_data(tx_config,all_scans_cell,filepath,force)
    num_elements = size(all_scans_cell{1},2);
    lengths = zeros(tx_config.no_transmits,1);
    for i = 1:tx_config.no_transmits
        lengths(i) = size(all_scans_cell{i},1);
    end
    max_length = max(lengths);

    rx_scans = complex(single(zeros(max_length, num_elements,tx_config.no_transmits)),0);
    for i = 1:tx_config.no_transmits
        scans = all_scans_cell{i};
    
        difference = max_length - length(scans);
        padded_zeros = complex(zeros(difference,num_elements),0);
        
        padded_scans = [scans; padded_zeros];
        rx_scans(:,:,i) = padded_scans;
    end
    
    if(isfile(filepath) && force == false)
        fprintf("File exists: %s\n",filepath);
    else
        fprintf("Saving file: '%s'\n",filepath);
        save(filepath, "rx_scans","tx_config",'-v7.3')
    end
end
 