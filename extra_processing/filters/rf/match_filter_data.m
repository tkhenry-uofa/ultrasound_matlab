function filtered_data_array = match_filter_data(data_array,filter,fs)

    data_dim = size(data_array);
    signal_length = data_dim(1);

    use_fft = true;
    
    if(use_fft)
        padded_filter = zeros(1,signal_length);
        padded_filter(1:length(filter)) = filter;
        
        filter_spectrum = fft(padded_filter);
        data_spectrum = fft(data_array,signal_length,1);
        
        % The reshape lets MATLAB broadcast the filter across the other two dims
        filtered_data_spectrum = data_spectrum .* reshape((filter_spectrum), [], 1, 1);
        
        filtered_data_array = single(real(ifft(filtered_data_spectrum,signal_length,1)));

        f_x = linspace(-fs/2,fs/2,signal_length);

        % figure();
        % subplot 511
        % plot(data_array(:,1,235));
        % title("Raw Data")
        % subplot 512
        % plot(f_x,abs(fftshift(data_spectrum(:,1,235))));
        % title("Data Spectrum");
        % 
        % subplot 513
        % plot(f_x,abs(fftshift(filter_spectrum)));
        % title("Filter Spectrum");
        % subplot 514
        % plot(f_x,abs(fftshift(filtered_data_spectrum(:,1,235))));
        % title("Filtered Spectrum");
        % 
        % subplot 515
        % plot(filtered_data_array(:,1,64));
        % title("Filtered Data")
    else
        
        filtered_data_array = zeros(data_dim,'single');
        
        for i=1:data_dim(2)
            for j=1:data_dim(3)
                filtered_data_array(:,i,j) = conv(data_array(:,i,j),filter,'same');
            end
        end

    end

end