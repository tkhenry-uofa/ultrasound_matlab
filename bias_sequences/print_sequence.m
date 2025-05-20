function print_sequence(sequence)

    pause_time = 0.5;
    total = zeros(size(sequence{1}));

    transmit = ones(128);
    figure()

    tx_title = sprintf("Bias Pattern \n Transmit");
    rx_title = sprintf("Bias Pattern \n Receive: ");

    for i = 1:length(sequence)
        clims = [0 1];

        
        imagesc(sequence{i}', clims)
        axis("image")
        % title(rx_title + string(i));
        title("Receive Bias Pattern");
        pause(pause_time);

        
        total = total + sequence{i};
    end
    imagesc(total')
    title("All Receive Elements");
    axis("image")
end
