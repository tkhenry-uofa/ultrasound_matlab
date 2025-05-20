clear all;

addpath("C:\Users\tkhen\source\repos\cuda_toolkit\ogl_beamforming_cuda\helpers")

% addpath("C:\Users\tkhen\source\repos\ogl_beamforming\helpers")

if isempty(matlab.project.currentProject)
    proj = openProject("../Ultrasound-Beamforming/Ultrasound-Beamforming.prj");
end

vrs_path   = "vrs_data/mixes_paper/";

vrs_path = vrs_path + "reso_herc_div_col";
vrs_name = "250204_MN32-1_8MHz_ATS539_Resolution_HERCULES-Diverging-TxColumn_Intensity_00";


% bp.das_shader_id = 0; % Forces
bp.das_shader_id = 2; % Hercules

vrs_name   = vrs_name + ".vrs";

pipe_name = '\\.\pipe\beamformer_data_fifo';
smem_name = 'Local\ogl_beamformer_parameters';
pipe_output = '\\.\pipe\beamformer_output_fifo'; % hardcoded in the lib rn

vrs  = VRSFile(fullfile(vrs_path, vrs_name));
data = vrs.GetData();

load(fullfile(vrs_path, "postEnv.mat"), "Trans", "Receive", "Resource", "TW");
load(fullfile(vrs_path, "preEnv.mat"), "sparseElements", "scan");

receive             = Receive([Receive.bufnum] == 1);
receive_orientation = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation;
recieve_elements    = scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation.GetElementCount(scan.Die);

transmit_count = max([receive.acqNum]);

bp.output_min_coordinate = struct('x', -20e-3, 'y', 0, 'z', 5e-3,   'w', 0);
bp.output_max_coordinate = struct('x',  20e-3, 'y', 0, 'z', 60e-3, 'w', 0);

bp.output_points.x = 1024;
bp.output_points.y = 1;
bp.output_points.z = 1024;
bp.output_points.w = 1; % Number of frames for averaging



bp.rf_raw_dim     = struct('x', size(data, 1), 'y', size(data, 2));
bp.dec_data_dim.x = max(1 + [receive.endSample] - [receive.startSample], [], "all");
bp.dec_data_dim.y = recieve_elements;
bp.dec_data_dim.z = transmit_count;
bp.dec_data_dim.w = 0; % Averaging

bp.sampling_frequency = receive(1).samplesPerWave * Trans.frequency(1) * 1e6;
bp.center_frequency   = Trans.frequency * 1e6;
bp.speed_of_sound     = Resource.Parameters.speedOfSound;

bp.time_offset = TW(1).Parameters(3) / bp.center_frequency;

bp.channel_mapping  = Trans.ConnectorES(scan.TransmitEvents(1).ImagingPattern.ReceiveOrientation.GetElements(scan.Die)) - 1;
if (exist('sparseElements'))
	bp.uforces_channels = sparseElements - 1;
else
	bp.uforces_channels = 0:(recieve_elements - 1);
end

% bp.channel_offset = 0 + recieve_elements * uint16(receive_orientation.Contains(tobe.Orientation.Column));

dieSize       = scan.Die.GetSize();

bp.xdc_transform = single([
                    1,0,0,dieSize(1)/2;
                    0,1,0,dieSize(2)/2;
                    0,0,1,0;
                    0,0,0,1;
                    ]);
bp.xdc_element_pitch = scan.Die.Pitch;


bp.focal_depths  = ones(1,transmit_count) .* scan.TransmitEvents(1, 1).FocalDepth;
bp.transmit_angles  = single(zeros(1,transmit_count));

% NOTE: plane and position along plane normal for beamforming 2D HERCULES
bp.beamform_plane = 0;
bp.off_axis_pos = 0;
fc = bp.center_frequency;
fs = bp.sampling_frequency;

bp.f_number = 1;


% cs_stages = uint8([
% 	OGLShaderStage.CUDA_DECODE, ...
%     OGLShaderStage.CUDA_HILBERT, ...
% 	OGLShaderStage.DAS, ...
% 	OGLShaderStage.MIN_MAX, ...
% ]);

% cs_stages = uint8([
% 	OGLShaderStage.CUDA_DECODE, ...
% 	OGLShaderStage.DAS, ...
% 	OGLShaderStage.MIN_MAX, ...
% ]);



cs_stages = uint8([
	OGLShaderStage.HADAMARD, ...
	OGLShaderStage.DAS, ...
	OGLShaderStage.MIN_MAX, ...
]);

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

