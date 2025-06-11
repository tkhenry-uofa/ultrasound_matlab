function walsh_matrix = generate_walsh(mat_size)
    
    H = hadamard(mat_size);

    walsh_matrix = zeros(mat_size);

    % Reorder the rows by zero crossings
    for i=1:mat_size
        H_row = H(i,:);
        zero_crossings = 0;
        for j = 1:(mat_size-1)
            zero_crossings = zero_crossings + (sign(H_row(j)) ~= sign(H_row(j+1)));
        end
        walsh_matrix(zero_crossings+1,:) = H(i,:); 
    end
end