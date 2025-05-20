function processed_volume = process_volume(volume,dynamic_range,threshold,power)

    if(~exist("power","var"))
        power = false;
    end

    gamma = 1;
    mag_volume = abs(volume); 

    rf_dynamic_range = 10^(-dynamic_range/20);
    if(exist("threshold","var"))
        threshold_value = 10^(threshold/20);
        mag_volume = min(mag_volume, threshold_value);   
    end
    max_value = max(mag_volume,[],"all");
    power_volume = (mag_volume/max_value).^(2 * gamma);

    if power
        % processed_volume = max(power_volume, rf_dynamic_range);
        processed_volume = power_volume;
    else
        processed_volume = 10*log10(power_volume);
        processed_volume = max(processed_volume, -dynamic_range);
    end
    
end