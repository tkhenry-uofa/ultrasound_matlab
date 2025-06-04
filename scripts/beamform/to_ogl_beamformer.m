clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\ogl_beamforming_cuda\out")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path_root   = "vrs_data/readi/stage_motion/";
dataset_name = "250423_MN32-5_reso_motion_FORCES-TxColumn";
data_file_range = 0:0;

data_path = data_path_root + dataset_name + "/";
params_path = data_path + dataset_name + ".bp";

bp = load_and_parse_bp(params_path);
frame_data_cell = load_vrs_data(data_path + dataset_name, data_file_range, bp.rf_raw_dim);
frame_count = length(frame_data_cell);

sample_count = single(bp.dec_data_dim(1));
rx_channel_count = single(bp.dec_data_dim(2));
transmit_count = single(bp.dec_data_dim(3));

fc = bp.center_frequency;
fs = bp.sampling_frequency;

tx_region = 2.5/1000; % How far down to crop to avoid hearing the transmit pulse
[frame_data_cell, bp.rf_raw_dim] = crop_and_blank_tx(frame_data_cell, bp, tx_region);


% Render settings
x_range = [-20, 20]/1000;
z_range = [5, 70]/1000;
y_range = [0, 0]/1000;


bp.output_min_coordinate = [x_range(1), y_range(1), z_range(1), 0];
bp.output_max_coordinate = [x_range(2), y_range(2), z_range(2), 0];

bp.output_points = [ 512 1 1024 0 ];

bp.interpolate = false; % bool

bp.off_axis_pos = 0; 
bp.beamform_plane = 0;
bp.f_number = 0;

[bp_head, bp_ui, bp_tail, arrays] = split_bp(bp);

transmit_mode = acquisition.TransmitModes(bp.transmit_mode);
if(transmit_mode == acquisition.TransmitModes.RowTxRowRx || ...
   transmit_mode == acquisition.TransmitModes.ColTxRowRx )
    image_plane = 1;
else
    image_plane = 0;
end


cs_stages = uint8([
	OGLShaderStage.CUDA_DECODE, ...
    OGLShaderStage.CUDA_HILBERT, ...
	OGLShaderStage.DAS
]);
% 
% cs_stages = uint8([
% 	OGLShaderStage.CUDA_DECODE, ...
% 	OGLShaderStage.DAS
% ]);

% cs_stages = uint8([
% 	OGLShaderStage.DECODE, ...
%     OGLShaderStage.CUDA_HILBERT, ...
% 	OGLShaderStage.DAS
% ]);

% cs_stages = uint8([
% 	OGLShaderStage.DECODE, ...
% 	OGLShaderStage.DAS
% ]);



if ~libisloaded('ogl_beamformer_lib'), loadlibrary('ogl_beamformer_lib'); end

fprintf("Setting params\n")

timeout_ms = 0;

%% Sending Params
calllib('ogl_beamformer_lib', 'beamformer_push_parameters_ui', bp_ui, timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_parameters_head', bp_head, timeout_ms);
calllib('ogl_beamformer_lib', 'set_beamformer_pipeline', cs_stages, numel(cs_stages));

calllib('ogl_beamformer_lib', 'beamformer_push_channel_mapping', arrays.channel_mapping, numel(arrays.channel_mapping), timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_sparse_elements', arrays.sparse_elements, numel(arrays.sparse_elements), timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_focal_vectors', arrays.focal_vectors, transmit_count, timeout_ms);

%% Sending data


frame_data = frame_data_cell{1};
info = whos('frame_data');
data_size = info.bytes;

for i = 1:frame_count

    frame_data = frame_data_cell{i};
    if calllib('ogl_beamformer_lib', 'beamformer_push_data', frame_data, data_size, timeout_ms)
        calllib('ogl_beamformer_lib', 'beamformer_start_compute', image_plane);
    else
        warning("OGL Beamfroming Push Data Failed");
    end

    pause(0.5);
end

if true, unloadlibrary('ogl_beamformer_lib'); end



