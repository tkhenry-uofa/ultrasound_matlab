function [tx_array, config] = create_curved_array(config)

    tx_array = xdc_2d_array(config.cols,config.rows,config.width,config.width,...
        config.kerf,config.kerf,ones(config.cols,config.rows),1,1,[0 0 0]);

    tx_array = xdc_convex_focused_multirow(config.cols,config.width, config.rows, config.width * ones(config.rows,1),...
        config.kerf,config.kerf, config.curve_radius, 0, ones(config.cols,config.rows),1,1,[0 0 0]);

    xdc_center_focus(tx_array, [0,0,0]);

    % Force the automatic focusing to be flat
    xdc_focus(tx_array,0,[0,0,100000000000])
    xdc_impulse (tx_array, config.imp);
    xdc_excitation (tx_array, config.ex);
    ele_apodization(tx_array,(1:(config.cols)*(config.rows))',config.apro(:));
    set_tx_delays(tx_array,config.src,config.transmit,config.rows,config.cols,config.c,config.print);

    element_locs = xdc_get(tx_array);
    element_locs = element_locs(24:25,:);

    min = element_locs(1:2,1);
    max = element_locs(1:2,end);

    config.x_range = element_locs(1,1:config.cols);
    config.y_range = element_locs(2,1:config.cols:(length(element_locs)-config.cols+1));

    config.x_min = min(1);
    config.y_min = min(2);

    config.x_max = max(1);
    config.y_max = max(2);
end