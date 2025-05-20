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

%%
bp.decode          = raw_bp.decode_mode;
bp.beamform_plane  = raw_bp.beamform_mode;

bp.rf_raw_dim.x      = raw_bp.raw_data_dim(1);
bp.rf_raw_dim.y      = raw_bp.raw_data_dim(2);

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

bp.center_frequency = 6.25e6;

% Map transmit mode
bp.transmit_mode     = raw_bp.transmit_mode;


bp.f_number = 1;

bp.das_shader_id = 0; % Forces
% bp.das_shader_id = 2; % Hercules

bp.output_min_coordinate = struct('x', -20e-3, 'y', 0, 'z', 5e-3,   'w', 0);
bp.output_max_coordinate = struct('x',  20e-3, 'y', 0, 'z', 60e-3, 'w', 0);

bp.output_points.x = 1024;
bp.output_points.y = 1;
bp.output_points.z = 1024;
bp.output_points.w = 1; % Number of frames for averaging


cs_stages = uint8([
	OGLShaderStage.CUDA_DECODE, ...
    OGLShaderStage.CUDA_HILBERT, ...
	OGLShaderStage.DAS, ...
	OGLShaderStage.MIN_MAX, ...
]);

% cs_stages = uint8([
% 	OGLShaderStage.CUDA_DECODE, ...
% 	OGLShaderStage.DAS, ...
% 	OGLShaderStage.MIN_MAX, ...
% ]);


if libisloaded('ogl_beamformer_lib'), unloadlibrary('ogl_beamformer_lib'); end

loadlibrary('ogl_beamformer_lib')

fprintf("Setting params\n")
calllib('ogl_beamformer_lib', 'set_beamformer_parameters', smem_name, bp);
calllib('ogl_beamformer_lib', 'set_beamformer_pipeline', smem_name, cs_stages, numel(cs_stages));
fprintf("Sending data\n")

sync = false;

if sync

    output_counts_xyz.x = bp.output_points.x;%#ok<UNRCH>
    output_counts_xyz.y = bp.output_points.y;
    output_counts_xyz.z = bp.output_points.z;
    volume = single(zeros(bp.output_points.x, bp.output_points.y, bp.output_points.z));
    volumePtr = libpointer('singlePtr', volume);

    [~,~,~,volumePtr] = calllib('ogl_beamformer_lib', 'beamform_data_synchronized', ...
        pipe_name, smem_name, data, bp.rf_raw_dim, output_counts_xyz, volumePtr);
    fprintf("Received Response\n")

else
    calllib('ogl_beamformer_lib', 'send_data', pipe_name, smem_name, data, bp.rf_raw_dim); %#ok<UNRCH>
end


unloadlibrary('ogl_beamformer_lib')

% %%
% volume = reshape(volumePtr, size(volume));
% 
% processed_volume = process_volume(volume, 50);
% 
% volumeViewer(processed_volume);

