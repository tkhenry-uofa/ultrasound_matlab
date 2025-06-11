function [tx_array, config] = create_tx_array(config, transmit_no)

    element_count = config.rows*config.cols;
    sub_element_count = 1;
    if(isfield(config, 'curve_radius') && config.curve_radius ~= 0 )
        
        tx_array = xdc_convex_focused_multirow(config.cols, config.width, config.rows, config.width * ones(config.rows,1).',...
        config.kerf,config.kerf, config.curve_radius, 10000000, 3,1,[0 0 0]);

        sub_element_count = 3;
        
    else
        tx_array = xdc_2d_array(config.cols,config.rows,config.width,config.width,...
        config.kerf,config.kerf,ones(config.cols,config.rows),1,1,[0 0 0]);
    end

    
    xdc_center_focus(tx_array, [0,0,0]);

    % Force the automatic focusing to be flat
    xdc_focus(tx_array,0,[0,0,100000000000])
    xdc_impulse (tx_array, config.imp);
    xdc_excitation (tx_array, config.ex);

    apo = config.apo;

    encoded_sequences = [SequenceType.FORCES, SequenceType.UFORCES];
    if(ismember(config.sequence, encoded_sequences))
        if isfield(config,"walsh_ordering")&& config.walsh_ordering
            H = generate_walsh(config.no_transmits);
        else
            H = hadamard(config.no_transmits);
        end
        
        tx_pattern = H(:,transmit_no);

        apo = apo .* tx_pattern;
    end

    apo_array = zeros(element_count,sub_element_count);
    for i = 1:sub_element_count
        apo_array(:,i) = reshape(apo,[],1);
    end

    ele_apodization(tx_array,(1:element_count)',apo_array);
    set_tx_delays(tx_array,config.src,config.transmit,config.rows,config.cols,config.c,sub_element_count);

    element_locs = xdc_get(tx_array);
    element_locs = element_locs(24:25, 1:sub_element_count:end);

    min = element_locs(1:2,1);
    max = element_locs(1:2,end);

    config.x_range = element_locs(1,1:config.cols);
    config.y_range = element_locs(2,1:config.cols:(length(element_locs)-config.cols+1));

    config.x_min = min(1);
    config.y_min = min(2);

    config.x_max = max(1);
    config.y_max = max(2);

    config.sub_element_count = sub_element_count;
end