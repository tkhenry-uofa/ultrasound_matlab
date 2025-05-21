function sequence = alt_lines( no_transmits, array_size, active_rows, offset_crosses )
    
    sequence = cell(no_transmits,1);
    full = zeros(array_size*2);
    
    for i = 1:no_transmits
        % Initialize matrices for each cross
        cross1 = zeros(array_size); % No offset
        cross2 = zeros(array_size); % Offset in y-direction
        cross3 = zeros(array_size); % Offset in x-direction
        cross4 = zeros(array_size); % Offset in both directions

        if mod(i,2) == 0

            % Define cross1 with no offsets
            cross1(active_rows(i), :) = 1;
        
            % Define cross2 with offset in y-direction only
            cross2(active_rows(i)+offset_crosses, :) = 1;
    
            % Define cross1 with no offsets
            cross3(active_rows(i),:) = 1;
        
            % Define cross2 with offset in y-direction only
            cross4(active_rows(i)+offset_crosses, :) = 1;
        
        else
            % Define cross1 with no offsets
            cross2(:,active_rows(i)) = 1;
        
            % Define cross2 with offset in y-direction only
            cross1(:,active_rows(i)) = 1;
    
            % Define cross1 with no offsets
            cross4(:,active_rows(i)+offset_crosses) = 1;
        
            % Define cross2 with offset in y-direction only
            cross3(:,active_rows(i)+offset_crosses) = h1;


        end
    
        
    
        % Arrange crosses into the final pattern
        cross = [cross1, cross2; cross3, cross4];

        imagesc(cross)
    
        % Store each cross sequence
        sequence{i} = cross;
        full = full + cross;
    end
end
