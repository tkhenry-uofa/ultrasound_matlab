function [cropped_data_cell, final_rf_dims] = crop_and_blank_tx(data_cell, bp, tx_depth)

    % Remove the extra zeros at the end of each channel and blank the transmit region
    frame_count = length(data_cell);
    cropped_data_cell = cell(1, frame_count);

    sample_count = bp.dec_data_dim(1);
    transmit_count = bp.dec_data_dim(3);
    three_D_dims = [sample_count, transmit_count, bp.rf_raw_dim(2)];

    expected_length = sample_count * transmit_count;
    final_rf_dims = bp.rf_raw_dim;
    final_rf_dims(1) = expected_length;

    for i = 1:frame_count
        cropped_data = data_cell{i};
        cropped_data = cropped_data(1:expected_length, :);

        cropped_data = reshape(cropped_data, three_D_dims);
        tx_end_sample = round(2 * tx_depth * bp.sampling_frequency / bp.speed_of_sound);
        
        cropped_data(1:tx_end_sample,:,:) = 0;

        cropped_data_cell{i} = cropped_data;
    end

end