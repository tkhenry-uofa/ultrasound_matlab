function animate_frames(frame_cell, frame_rate)

    if(~exist("frame_rate", "var"))
        frame_rate = 4;
    end
    
    frame_count = numel(frame_cell);
    image_size = size(frame_cell{1});
    
    dr = -1 * min(frame_cell{1},[],"all");

    frame_array_size = [image_size, frame_count];
    frame_array = zeros(frame_array_size, "uint8");
    
    for i = 1:frame_count
        frame = frame_cell{i};
        frame_array(:,:,i) = uint8((frame + dr) * 255/dr);
    end

    implay(frame_array,frame_rate);

end