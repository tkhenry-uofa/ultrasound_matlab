function sequence = single_array_diag( no_transmits, array_size, active_rows )
    sequence = cell(no_transmits,1);
    full = zeros(array_size);
    
    for i = 1:no_transmits
        % Initialize matrices for each cross
        cross1 = zeros(array_size); % No offset
    
        % Define cross1 with no offsets
        cross1(active_rows(i), :) = 1;
        cross1(:, active_rows(i)) = 1;
        cross1(active_rows(i), active_rows(i)) = 0; % Middle element can't be read
    
    
        % Store each cross sequence
        sequence{i} = cross1;
        full = full + cross1;
    end
end