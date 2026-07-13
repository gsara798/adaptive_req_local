function RUN = load_run_config(profile_name, varargin)
%LOAD_RUN_CONFIG Load and resolve a full adaptive_req run configuration.
%
% Usage:
%   RUN = adaptive_req.config.load_run_config('default', 'RootDir', root_dir);
%
%   RUN = adaptive_req.config.load_run_config( ...
%       'test_03_baseline', ...
%       'RootDir', root_dir);
%
% This function returns a resolved RUN struct with:
%   RUN.CFG
%   RUN.cfg
%   RUN.feat_cfg
%   RUN.req_options
%   RUN.req_preview
%   RUN.EXP
%   RUN.REQ
%   RUN.PLOT
%   RUN.SAVE

p = inputParser;
p.FunctionName = 'adaptive_req.config.load_run_config';

addRequired(p, 'profile_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'ConfigDir', '', @(x) ischar(x) || isstring(x));

parse(p, profile_name, varargin{:});

profile_name = char(p.Results.profile_name);
root_dir = char(p.Results.RootDir);
config_dir = char(p.Results.ConfigDir);

if isempty(config_dir)
    config_dir = fullfile(root_dir, 'configs');
end

if exist(config_dir, 'dir')
    addpath(config_dir);
end

CFG = adaptive_req.config.default_run_config();

if ~strcmpi(profile_name, 'default')

    if exist(profile_name, 'file') ~= 2
        error('Configuration profile not found: %s', profile_name);
    end

    profile_fun = str2func(profile_name);
    CFG = profile_fun(CFG);

end

CFG.EXP.profile_name = profile_name;
CFG.EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));

cfg = adaptive_req.config.default_sim_config();
cfg = apply_struct_overrides(cfg, CFG.SIM);

REQ = CFG.REQ;

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

SAVE = CFG.SAVE;

SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
SAVE.run_name = sprintf('%s_%s_%s', ...
    CFG.EXP.name, CFG.EXP.sampling_mode, CFG.EXP.timestamp);

SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.figure_dir = fullfile(SAVE.output_dir, 'figures');
SAVE.step_diag_dir = fullfile(SAVE.figure_dir, 'step_diagnostics');
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

make_dir_if_needed(SAVE.output_dir);
make_dir_if_needed(SAVE.figure_dir);
make_dir_if_needed(SAVE.step_diag_dir);
make_dir_if_needed(SAVE.table_dir);
make_dir_if_needed(SAVE.data_dir);

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
RUN.OUTPUT = CFG.OUTPUT;

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

function make_dir_if_needed(folder_path)

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

end
