% Rx distances are mapped to a hann window s.t. the max distance is 0 and
% the min is 1 for this voxel. This will give a different sized apature
% depending on voxel position.
function aprodizations = hann_apo(vox_loc,element_locs,tx_config, f_number)

    lat_dists = sqrt((vox_loc(1) - element_locs(1,:)).^2 + (vox_loc(2) - element_locs(2,:)).^2)'; 

    % Normalize to 0-1
    ratio = ( f_number * lat_dists )./vox_loc(3);
    
    aprodizations = cos(pi.*ratio).^2;
    aprodizations(ratio > 0.5) = 0;
    
end