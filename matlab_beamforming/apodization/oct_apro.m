% Rx distances are mapped to a hann window s.t. the max distance is 0 and
% the min is 1 for this voxel. This will give a different sized apature
% depending on voxel position.
function aprodizations = oct_apro(vox_loc,element_locs,tx_config)


    % rx_dists = [abs(vox_loc(1) - element_locs(1,:)); abs(vox_loc(2) - element_locs(2,:))]'; 
    rx_dists = sqrt((vox_loc(1) - element_locs(1,:)).^2 + (vox_loc(2) - element_locs(2,:)).^2)'; 

    max_dist = max(rx_dists);
    min_dist = min(rx_dists);

    % Normalize to 0-1
    rx_dists = ( rx_dists - min_dist )./(max_dist-min_dist);
    
    N = 128;
    hann_window = hann(2*N);
    pos_window = hann_window(N+1:end);
    neg_window = hann_window(1:N);
    
    indecies = round(rx_dists * N);
    
    indecies = indecies + 1; % One based indexing is evil
    indecies( indecies > N ) = N;
    pos_apros = pos_window(indecies);
    neg_apros = neg_window(indecies);

    aprodizations = pos_apros .* neg_apros;
    aprodizations = aprodizations./max(aprodizations);


    
end