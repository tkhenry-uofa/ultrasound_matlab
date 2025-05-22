function rf_data_cell = load_vrs_data(path_stub,frame_range,rf_data_dim)

    frame_count = length(frame_range);
    rf_data_cell = cell(1,frame_count);

    file_paths = path_stub + compose('_%02i.zst', frame_range).';

    for i = 1:frame_count
        data_file = fopen(file_paths(i), "r");
        raw_data = fread(data_file, '*uint8');
        data = ornot_zstd_decompress_mex(raw_data);
        rf_data_cell{i} = reshape(data, rf_data_dim(1),rf_data_dim(2));
        fclose(data_file); 
    end

end