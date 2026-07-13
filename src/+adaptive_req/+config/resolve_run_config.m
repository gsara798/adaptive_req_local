function RUN = resolve_run_config(CFG, varargin)
%RESOLVE_RUN_CONFIG Resolve a scalar configuration into a RUN struct.
%
% This function assumes that CFG contains scalar values for the fields used
% by one aperture sweep condition.
%
% It builds:
%   RUN.cfg
%   RUN.feat_cfg
%   RUN.req_options
%   RUN.req_preview
%   RUN.EXP
%   RUN.REQ
%   RUN.PLOT
%   RUN.SAVE
%   RUN.OUTPUT

p = inputParser;
p.FunctionName = 'adaptive_req.config.resolve_run_config';

addRequired(p, 'CFG', @isstruct);

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'MakeDirs', true, @(x) islogical(x) || isnumeric(x));

parse(p, CFG, varargin{:});

root_dir = char(p.Results.RootDir);
output_dir = char(p.Results.OutputDir);
make_dirs = logical(p.Results.MakeDirs);

if ~isfield(CFG.EXP, 'timestamp') || isempty(CFG.EXP.timestamp)
    CFG.EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
end

cfg = adaptive_req.config.default_sim_config();
cfg = apply_struct_overrides(cfg, CFG.SIM);

REQ = CFG.REQ;

validate_scalar_req_field(REQ, 'M');
validate_scalar_req_field(REQ, 'cs_guess');
validate_scalar_req_field(REQ, 'gamma_win');
validate_scalar_req_field(REQ, 'pad_factor');

feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', REQ.M, ...
    'cs_guess', REQ.cs_guess, ...
    'gamma_win', REQ.gamma_win, ...
    'pad_factor', REQ.pad_factor);

req_options = { ...
    'Nbins', REQ.Nbins, ...
    'Nbins_auto_oversample', REQ.Nbins_auto_oversample, ...
    'Nbins_min', REQ.Nbins_min, ...
    'smooth_sigma', REQ.smooth_sigma, ...
    'use_donut', REQ.use_donut, ...
    'donut_cs_min', REQ.donut_cs_min, ...
    'donut_cs_max', REQ.donut_cs_max, ...
    'donut_taper_rel', REQ.donut_taper_rel, ...
    'apply_donut_to_final_map', REQ.apply_donut_to_final_map};

[req_preview, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

PLOT = CFG.PLOT;

PLOT.store_wavefields = ...
    PLOT.show_selected_step_wavefield || ...
    PLOT.show_all_step_wavefields;

if isfield(CFG, 'OUTPUT')
    OUTPUT = CFG.OUTPUT;
else
    OUTPUT = struct();
end

if ~isfield(OUTPUT, 'store_req_curve')
    OUTPUT.store_req_curve = false;
end

if ~isfield(OUTPUT, 'store_req_mapping')
    OUTPUT.store_req_mapping = true;
end

if ~isfield(OUTPUT, 'compute_global_req')
    OUTPUT.compute_global_req = false;
end

if ~isfield(OUTPUT, 'store_global_req_mapping')
    OUTPUT.store_global_req_mapping = true;
end

if ~isfield(OUTPUT, 'store_req_metadata')
    OUTPUT.store_req_metadata = true;
end

if ~isfield(OUTPUT, 'store_feature_struct')
    OUTPUT.store_feature_struct = false;
end

if ~isfield(OUTPUT, 'verbose')
    OUTPUT.verbose = true;
end

SAVE = CFG.SAVE;

if isempty(output_dir)

    SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
    SAVE.run_name = sprintf('%s_%s_%s', ...
        CFG.EXP.name, CFG.EXP.sampling_mode, CFG.EXP.timestamp);

    SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);

else

    SAVE.output_dir = output_dir;
    [SAVE.root_dir, SAVE.run_name] = fileparts(SAVE.output_dir);

end

SAVE.figure_dir = fullfile(SAVE.output_dir, 'figures');
SAVE.step_diag_dir = fullfile(SAVE.figure_dir, 'step_diagnostics');
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

if make_dirs
    make_dir_if_needed(SAVE.output_dir);
    make_dir_if_needed(SAVE.figure_dir);
    make_dir_if_needed(SAVE.step_diag_dir);
    make_dir_if_needed(SAVE.table_dir);
    make_dir_if_needed(SAVE.data_dir);
end

RUN = struct();

RUN.CFG = CFG;
RUN.cfg = cfg;
RUN.feat_cfg = feat_cfg;
RUN.req_options = req_options;
RUN.req_preview = req_preview;

RUN.EXP = CFG.EXP;
RUN.REQ = CFG.REQ;
RUN.PLOT = PLOT;
RUN.SAVE = SAVE;
RUN.OUTPUT = OUTPUT;

end

function S = apply_struct_overrides(S, overrides)

names = fieldnames(overrides);

for i = 1:numel(names)

    name_i = names{i};
    value_i = overrides.(name_i);

    if isempty(value_i)
        continue;
    end

    S.(name_i) = value_i;

end

end

function validate_scalar_req_field(REQ, field_name)

value = REQ.(field_name);

if ~(isnumeric(value) && isscalar(value))
    error('REQ.%s must be scalar when resolving one RUN condition.', field_name);
end

end

function make_dir_if_needed(folder_path)

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

end
