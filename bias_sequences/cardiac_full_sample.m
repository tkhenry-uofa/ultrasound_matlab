function sequence = cardiac_full_sample()
    
    num_transmits = 12;
    row_count = 48;
    column_count = 64;

    % Sub array dimensions
    sub_row_count = 12;
    sub_col_count = 32;

    % Number of rows grouped together for one rx
    row_groups = 1;

    % Offset between the left and right rows
    offset = 0;

    sequence = cell(num_transmits,1);
    for T = 1:num_transmits
        % Initialize matrices for each cross
        left_subs = zeros(sub_col_count, sub_row_count);
        right_subs = zeros(sub_col_count, sub_row_count);

        left_i = (T-1) * row_groups * offset + 1;
        right_i = left_i + offset;

        left_subs(:, T) = 1;
        right_subs(:, T) = 1;

        % left_subs(:, left_i:left_i+1) = 1;
        % right_subs(:, right_i:right_i+1) = 1;
    
        % Arrange crosses into the final pattern
        full_array = [left_subs, left_subs, left_subs, left_subs; right_subs, right_subs, right_subs, right_subs];
    
        % Store each cross sequence
        sequence{T} = full_array;
    end

end

