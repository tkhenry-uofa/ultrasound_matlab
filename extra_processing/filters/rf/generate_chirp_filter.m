function filter = generate_chirp_filter(fc,fs,bandwidth,chirp_length)

    % Impulse response, 1 cyc hamming modulated sin
    cycles = 1;
    impulse_response = sin(2*pi*fc*(0:1/fs:cycles/fc)).';
    impulse_response = impulse_response .* hamming(length(impulse_response));
    
    % Chirp excitation
    f1 = bandwidth(1); 
    f2 = bandwidth(2);  
    BW = f2 - f1;                       
    tapering = 0.20;         
    % Create Chirp excitation
    t = 0:1/fs:chirp_length-1/fs;
    F = f1 + BW/(2*chirp_length) * t;
    excitation = sin(2*pi*F.*t);
    excitation = excitation.*tukeywin(length(excitation),tapering)';
    
    filter = single(fliplr(conv(excitation, impulse_response)));

end