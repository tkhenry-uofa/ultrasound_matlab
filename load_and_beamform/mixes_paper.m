clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\test_app\matlab_lib")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

average_range = 0:7;
total_frames = length(average_range);

data_folder   = "vrs_data/mixes_paper/cyst_4mhz/hypo_div_tx_col/";
data_file_name = "250220_MN32-1_ATS539_Cyst_HERCULES-Diverging-45-TxColumn";

data_paths = data_folder + data_file_name + compose('_%02i.zst', average_range).';
% params_path = data_folder + data_file_name + '.bp';
params_path = data_folder + 'parameters_fixed.bp';

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));

frame_data = cell(1,total_frames);
for i = 1:total_frames
    data_file = fopen(data_paths(i), "r");
    raw_data = fread(data_file, '*uint8');
    data = ornot_zstd_decompress_mex(raw_data);
    frame_data{i} = reshape(data, raw_bp.raw_data_dim(1),raw_bp.raw_data_dim(2));
    fclose(data_file); 
end

% data(1:2,:) = 0;

%%
bp.decode          = raw_bp.decode_mode;
bp.beamform_plane  = raw_bp.beamform_mode;

% bp.rf_raw_dim = raw_bp.rf_raw_dim;
bp.rf_raw_dim.x      = raw_bp.raw_data_dim(1);
bp.rf_raw_dim.y      = raw_bp.raw_data_dim(2);
% bp.dec_data_dim    = raw_bp.dec_data_dim;
bp.dec_data_dim.x    = raw_bp.decoded_data_dim(1);
bp.dec_data_dim.y    = raw_bp.decoded_data_dim(2);
bp.dec_data_dim.z    = raw_bp.decoded_data_dim(3);
bp.dec_data_dim.w    = raw_bp.decoded_data_dim(4);

% Map transducer properties
bp.xdc_element_pitch = raw_bp.transducer_element_pitch;
bp.xdc_transform     = raw_bp.transducer_transform_matrix;

% Map channel and angle related fields
bp.channel_mapping  = raw_bp.channel_mapping;
bp.transmit_angles  = raw_bp.steering_angles;
bp.focal_depths     = raw_bp.focal_depths;

% Map acoustic parameters
bp.speed_of_sound    = raw_bp.speed_of_sound;
bp.center_frequency  = raw_bp.center_frequency;
bp.sampling_frequency= raw_bp.sampling_frequency;
bp.time_offset       = raw_bp.time_offset;

% Map transmit mode
bp.transmit_mode     = raw_bp.transmit_mode;

% NOTE: plane and position along plane normal for beamforming 2D HERCULES
bp.beamform_plane = 0;
bp.off_axis_pos = 0;
fc = bp.center_frequency;
fs = bp.sampling_frequency;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = [-10, 10]/1000;
y_range = [-40, 40]/1000;
z_range = [5, 80]/1000;


bp.output_min_coordinate = struct('x', x_range(1), 'y', y_range(1), 'z', z_range(1),   'w', 0);
bp.output_max_coordinate = struct('x', x_range(2), 'y', y_range(2), 'z', z_range(2), 'w', 0);

lateral_resolution = 0.0003;
axial_resolution = 0.0003;

bp.output_points.x = floor((bp.output_max_coordinate.x - bp.output_min_coordinate.x)/lateral_resolution );
bp.output_points.y = floor((bp.output_max_coordinate.y - bp.output_min_coordinate.y)/lateral_resolution );
bp.output_points.z = floor((bp.output_max_coordinate.z - bp.output_min_coordinate.z)/axial_resolution );
bp.output_points.w = 1; % Number of frames for averaging

if bp.output_points.y <= 0
    bp.output_points.y = 1;
end

bp.readi_group_size = 128;
bp.readi_group_id = 0;

x_range = linspace(x_range(1), x_range(2), bp.output_points.x);
y_range = linspace(y_range(1), y_range(2), bp.output_points.y);
z_range = linspace(z_range(1), z_range(2), bp.output_points.z);

%% Data prep
if libisloaded('cuda_transfer'), unloadlibrary('cuda_transfer'); end

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

loadlibrary('cuda_transfer')

fprintf("Sending data\n")

output_counts_xyz.x = bp.output_points.x;
output_counts_xyz.y = bp.output_points.y;
output_counts_xyz.z = bp.output_points.z;

%%
% 0 = Forces, 1 = Uforces, 2 = Hercules, 100 = Mixes
%bp.das_shader_id = 0;
bp.das_shader_id = 100;

bp.f_number = 1;

mixes_versions = ["Fully Sampled", "MiXES 64", "MiXES 32", "MiXES 16", "MiXES 8", "MiXES 32 Staggered", "MiXES 16 Staggered", "MiXES 8 Staggered"];
mixes_counts = [128, 64, 32, 16, 8, 64, 32, 16]; % The staggered mixes have half the transmits but the same amount of crosses                                                % As it would require 4 64x64 arrays
mixes_offsets = [0, 0, 0, 0, 0, 2, 4, 8];
mixes_spacings = round(128./mixes_counts);

mixes_rows = cell(1, length(mixes_versions));

for i=1:length(mixes_versions)
    mixes_rows{i} = 0:mixes_spacings(i):127;
end

volumes_cell = cell(1,total_frames);

% Complex volumes aren't supported so they're interleaved
volume_size = [bp.output_points.x, bp.output_points.y, bp.output_points.z];
interleaved_volume_size = volume_size;
interleaved_volume_size(1) = interleaved_volume_size(1) * 2;


calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);

for j = 1:length(mixes_versions)

    fprintf("Beamforming %s.\n", mixes_versions(j));

    average_volume = zeros(volume_size);

    bp.mixes_count = mixes_counts(j);
    bp.mixes_offset = mixes_offsets(j);
    bp.mixes_rows = mixes_rows{j};

    for i = 1:total_frames

        calllib('cuda_transfer', 'set_beamformer_parameters', smem_name, bp);
    
        volume_ptr = libpointer('singlePtr', single(zeros(interleaved_volume_size)));
    
        [~,~,~,volume_ptr] = calllib('cuda_transfer', 'beamform_i16', ...
            pipe_name, smem_name, frame_data{i}, bp.rf_raw_dim, output_counts_xyz, volume_ptr);
        fprintf("Received Response\n")
        
        volume = reshape(volume_ptr,interleaved_volume_size );
        
        real_vol = volume(1:2:end,:,:);
        im_vol = volume(2:2:end,:,:);
        
        average_volume = average_volume + complex(real_vol, im_vol);    
        pause(0.3);
    
    end

    volumes_cell{j} = average_volume;
end

unloadlibrary('cuda_transfer')

%% Post processing 
average_volume = zeros(volume_size);

dynamic_range = 40;

threshold = 700;

processed_volumes = volumes_cell;
processed_images = cell(2,length(volumes_cell));

for g = 1:length(volumes_cell)
    processed_volume = process_volume(volumes_cell{g}, dynamic_range, threshold);
    processed_images{1,g} = squeeze(processed_volume(:, round(output_counts_xyz.y/2) + 4,:)).';
    processed_images{2,g} = squeeze(processed_volume(round(output_counts_xyz.x/2),:,:)).';
    processed_volumes{g} = flip(processed_volume,3);
end

save_folder = data_folder + "/figures_no_av";
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

for idx = 1:size(processed_images, 2)
    fig = figure;
    t = tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'none');
    sgtitle(mixes_versions(idx));  % overall title
    
    % First tile: y_image
    ax1 = nexttile;
    imagesc(y_range * 1000, z_range * 1000, processed_images{2, idx});
    colormap(ax1, 'gray');
    axis(ax1, 'image');
    caxis(ax1, [-dynamic_range, 0]);
    colorbar(ax1);
    
    % Second tile: x_image
    ax2 = nexttile;
    imagesc(x_range * 1000, z_range * 1000, processed_images{1, idx});
    colormap(ax2, 'gray');
    axis(ax2, 'image');
    caxis(ax2, [-dynamic_range, 0]);
    colorbar(ax2);
    
    % Save the figure to the provided folder using the version as filename
    filename = fullfile(save_folder, sprintf('%s.png', mixes_versions(idx)));
    saveas(fig, filename);
    filename = fullfile(save_folder, sprintf('%s.fig', mixes_versions(idx)));
    saveas(fig, filename);
    close(fig);


    volume_path = sprintf(save_folder + "/%s_volume.mat", mixes_versions(idx));

    current_volume = processed_volumes{idx};
    save(volume_path, 'current_volume');

    % fig = uifigure();
    % view = viewer3d(fig);
    % view.CameraPosition = [volume_size(2)/2 volume_size(1)*2 volume_size(3)/2];
    % view.CameraUpVector = [0 0 1];
    % view.CameraTarget = [volume_size(2)/2 0 volume_size(3)/2];
    % view.CameraZoom = 0.9119;
    % vol_view = volshow(processed_volumes{idx}, Parent=view);
    % 
    % 
    % 
    % vol_view.RenderingStyle = "MaximumIntensityProjection";
    % obj = ancestor(vol_view,'figure','toplevel');
    % I = getframe(obj);
    % imshow(I.cdata);
    % 
    % filename = fullfile(save_folder, sprintf('%s_volume.png', mixes_versions(idx)));
    % saveas(fig, filename);
    % filename = fullfile(save_folder, sprintf('%s_volume.fig', mixes_versions(idx)));
    % saveas(fig, filename);
    % close(fig);
    % 
      % Close the figure to avoid cluttering the workspace
end


%%
% figure(2);
% 
% % Main title
% % sgtitle(speed_str + ' Readi Low Resolution Images','FontSize',16);
% sgtitle('Stationary Readi Low Resolution Images','FontSize',16);
% 
% % Loop through the images and display each in the 2x8 grid
% for i = 1:readi_group_count
%     subplot(2, 2, i); 
%     colormap("gray");
%     imagesc((low_res_images{i}).');
%     title(['Image ' num2str(i)]);
%     xlabel("mm");
%     ylabel("mm");
% 
%     clim([-dynamic_range, 0]); 
% 
% 
% end



%%

function processed_volume = process_volume(volume,dynamic_range, threshold)

    mag_volume = abs(volume); 
    threshold_value = 10^(threshold/20);
    mag_volume = min(mag_volume, threshold_value);
    
    processed_volume = 20*log10(mag_volume/threshold_value);
    processed_volume = processed_volume - max(processed_volume,[],"all");
    processed_volume = max(processed_volume, -dynamic_range);
end