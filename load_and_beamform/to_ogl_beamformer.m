clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\ogl_beamforming_cuda\out")
addpath('C:\Users\tkhen\source\repos\ornot\core\lib');

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

data_path   = "vrs_data/match_filter_test/";

data_path = data_path + "250514_MN32-5_cyst_FORCES-TxColumn/";
data_file_name = "250514_MN32-5_cyst_FORCES-TxColumn_00";
params_path = data_path + '250514_MN32-5_cyst_FORCES-TxColumn.bp';

data_file_name   = data_file_name + ".zst";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

data_file = fopen(data_path + data_file_name, "r");
raw_data = uint8(fread(data_file));
fclose(data_file);

data = ornot_zstd_decompress_mex(raw_data);

raw_bp = ornot_bp_load_mex(convertStringsToChars(params_path));
data = reshape(data, raw_bp.raw_data_dim(1:2));
data_info = whos('data');
data_size = data_info.bytes;

%% BF parameter construction
% Sequence-specific data
bp_head.decode          = raw_bp.decode_mode;
bp_head.rf_raw_dim      = raw_bp.raw_data_dim(1:2);

bp_head.dec_data_dim    = raw_bp.decoded_data_dim;

transmit_count = bp_head.dec_data_dim(3);

% Map transducer properties
bp_head.xdc_element_pitch = raw_bp.transducer_element_pitch;
bp_head.xdc_transform     = raw_bp.transducer_transform_matrix;

bp_head.time_offset       = raw_bp.time_offset;

bp_head.transmit_mode     = raw_bp.transmit_mode;
bp_head.das_shader_id     = raw_bp.beamform_mode;

transmit_mode = acquisition.TransmitModes(bp_head.transmit_mode);

if(transmit_mode == acquisition.TransmitModes.RowTxRowRx || ...
   transmit_mode == acquisition.TransmitModes.ColTxRowRx )
    image_plane = 1;
else
    image_plane = 0;
end

% Readi stuff
bp_tail.readi_group_id = 0;
bp_tail.readi_group_size = transmit_count;

% Render settings
x_range = [-20, 20]/1000;
z_range = [5, 70]/1000;
y_range = [20, 20]/1000;


bp_ui.output_min_coordinate = [x_range(1), y_range(1), z_range(1), 0];
bp_ui.output_max_coordinate = [x_range(2), y_range(2), z_range(2), 0];

bp_ui.output_points = [ 512 1 1024 0 ];

bp_ui.speed_of_sound    = raw_bp.speed_of_sound;
bp_ui.center_frequency  = raw_bp.center_frequency;
bp_ui.sampling_frequency= raw_bp.sampling_frequency;

bp_ui.interpolate = false; % bool

bp_ui.off_axis_pos = 0; 
bp_ui.beamform_plane = 0;

bp_ui.f_number = 0;


channel_mapping = raw_bp.channel_mapping;
sparse_elements = raw_bp.sparse_elements;
if(sparse_elements(1) == -1)
    sparse_elements(1:transmit_count) = 1:transmit_count;
end
focal_depths = raw_bp.focal_depths;
% focal_depths(:) = focal_depths(1);

focal_angles = raw_bp.steering_angles;

focal_vectors = zeros(1, length(focal_depths) * 2);
focal_vectors(1:2:end) = focal_angles;
focal_vectors(2:2:end) = focal_depths;


cs_stages = uint8([
	OGLShaderStage.CUDA_DECODE, ...
    OGLShaderStage.CUDA_HILBERT, ...
	OGLShaderStage.DAS
]);
% 
cs_stages = uint8([
	OGLShaderStage.CUDA_DECODE, ...
	OGLShaderStage.DAS
]);



if libisloaded('ogl_beamformer_lib'), unloadlibrary('ogl_beamformer_lib'); end

loadlibrary('ogl_beamformer_lib')

fprintf("Setting params\n")

timeout_ms = 0;

%% Sending Params
calllib('ogl_beamformer_lib', 'beamformer_push_channel_mapping', channel_mapping, numel(channel_mapping), timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_sparse_elements', sparse_elements, numel(sparse_elements), timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_focal_vectors', focal_vectors, transmit_count, timeout_ms);

calllib('ogl_beamformer_lib', 'beamformer_push_parameters_ui', bp_ui, timeout_ms);
calllib('ogl_beamformer_lib', 'beamformer_push_parameters_head', bp_head, timeout_ms);
calllib('ogl_beamformer_lib', 'set_beamformer_pipeline', cs_stages, numel(cs_stages));

%% Sending data

if calllib('ogl_beamformer_lib', 'beamformer_push_data', data, data_size, timeout_ms)
    calllib('ogl_beamformer_lib', 'beamformer_start_compute', image_plane);
else
    warning("OGL Beamfroming Push Data Failed");
end

unloadlibrary('ogl_beamformer_lib')



