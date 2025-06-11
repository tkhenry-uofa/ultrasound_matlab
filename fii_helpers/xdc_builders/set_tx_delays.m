function t = set_tx_delays(tx_array,src_loc,transmit_type,row_count,column_count,c,sub_element_count)

    element_locs = xdc_get(tx_array);
    element_locs = element_locs(:,1:sub_element_count:end);
    element_count = length(element_locs);

    element_x_locs = element_locs(24,1:column_count);
    element_y_locs = element_locs(25,1:column_count:element_count);
    element_z_locs = element_locs(26,1:column_count);

    [~,Y] = ndgrid(element_z_locs, element_y_locs);

    delays = zeros(column_count,row_count);
    if transmit_type == TransmitType.Xdiv
        dist = sqrt( (element_x_locs-src_loc(1)).^2 + (element_z_locs-src_loc(3)).^2);
        element_delays = dist/c;
        element_delays = element_delays - min(element_delays, [],"all");

        for i = 1:row_count
            delays(:,i) = element_delays;
        end

    elseif transmit_type == TransmitType.Ydiv || transmit_type == TransmitType.Elevation_Focus
        dist = sqrt( (Y-src_loc(2)).^2 + (-src_loc(3)).^2 );
        element_delays = dist/c;

        if( src_loc(3) < 0 )
            element_delays = element_delays - min(element_delays, [],"all");
        else
            element_delays = max(element_delays, [],"all") - element_delays;
        end


        delays = element_delays;

        % for i = 1:column_count
        %     delays(i,:) = element_delays;
        % end
    
    elseif transmit_type == TransmitType.Focused
        dist = sqrt( (X-src_loc(1)).^2 + src_loc(3).^2);
        element_delays = -dist/c;
        element_delays = element_delays - max(element_delays, [],"all");

        for i = 1:row_count
            delays(:,i) = element_delays;
        end
    else
        t = 0;
        return
    end
   
    t = min(element_delays);
    delays = reshape(delays, [],1);
    delay_array = zeros(element_count,sub_element_count);
    for i = 1:sub_element_count
        delay_array(:,i) = delays;
    end

    ele_delay(tx_array, (1:element_count).', delay_array);
end

