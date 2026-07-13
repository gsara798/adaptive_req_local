%% test_04_mc_config_preview.m
% Preview the MC baseline configuration.
%
% This test:
%   1. Loads the test_04_mc_baseline profile.
%   2. Reads CFG.SWEEP.paths.
%   3. Builds the Cartesian parameter matrix.
%   4. Prints the expected number of MC output rows.
%   5. Builds one scalar baseline preview using the first value of each
%      swept parameter.
%
% This script does not run simulations.

clear; clc; close all;
format compact;

%% ========================================================================
% 1. Project setup
% ========================================================================

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% ========================================================================
% 2. Load raw MC baseline profile
% ========================================================================

CFG = adaptive_req.config.load_profile_config( ...
    'test_04_mc_baseline', ...
    'RootDir', root_dir);

fprintf('\nLoaded MC configuration profile.\n');
fprintf('Profile name = %s\n', CFG.EXP.profile_name);
fprintf('Experiment name = %s\n', CFG.EXP.name);
fprintf('Run mode = %s\n', CFG.EXP.run_mode);
fprintf('Sampling mode = %s\n', CFG.EXP.sampling_mode);

%% ========================================================================
% 3. Build sweep preview
% ========================================================================

[T_axes, T_conditions, CFG_base] = ...
    adaptive_req.config.build_sweep_preview(CFG);

fprintf('\nActive sweep axes.\n');
disp(T_axes);

fprintf('\nFirst parameter conditions.\n');
n_preview = min(20, height(T_conditions));
disp(T_conditions(1:n_preview, :));

%% ========================================================================
% 4. Print MC size estimate
% ========================================================================

n_conditions = height(T_conditions);

rows_per_condition = ...
    CFG.EXP.num_steps * ...
    CFG.EXP.num_realizations * ...
    CFG.EXP.num_patches;

expected_total_rows = n_conditions * rows_per_condition;

fprintf('\nMC size estimate.\n');
fprintf('Number of sweep conditions = %d\n', n_conditions);
fprintf('Rows per condition = %d\n', rows_per_condition);
fprintf('Expected total rows = %d\n', expected_total_rows);

fprintf('\nRows per condition breakdown.\n');
fprintf('num_steps = %d\n', CFG.EXP.num_steps);
fprintf('num_realizations = %d\n', CFG.EXP.num_realizations);
fprintf('num_patches = %d\n', CFG.EXP.num_patches);

%% ========================================================================
% 5. Build scalar baseline preview
% ========================================================================
% CFG_base is obtained from CFG by taking the first value of each swept
% parameter. This is useful for checking derived quantities such as
% lambda_guess, win_size, and k0_true.

cfg_preview = adaptive_req.config.default_sim_config();
cfg_preview = apply_struct_overrides(cfg_preview, CFG_base.SIM);

REQ = CFG_base.REQ;

feat_cfg_preview = adaptive_req.config.default_feature_config( ...
    'M', REQ.M, ...
    'cs_guess', REQ.cs_guess, ...
    'gamma_win', REQ.gamma_win, ...
    'pad_factor', REQ.pad_factor);

req_options_preview = { ...
    'Nbins', REQ.Nbins, ...
    'Nbins_auto_oversample', REQ.Nbins_auto_oversample, ...
    'Nbins_min', REQ.Nbins_min, ...
    'smooth_sigma', REQ.smooth_sigma, ...
    'use_donut', REQ.use_donut, ...
    'donut_cs_min', REQ.donut_cs_min, ...
    'donut_cs_max', REQ.donut_cs_max, ...
    'donut_taper_rel', REQ.donut_taper_rel, ...
    'apply_donut_to_final_map', REQ.apply_donut_to_final_map};

[req_preview, feat_cfg_preview] = adaptive_req.config.default_req_config( ...
    cfg_preview, feat_cfg_preview, req_options_preview{:});

fprintf('\nScalar baseline preview.\n');
fprintf('f0 = %.1f Hz\n', cfg_preview.f0);
fprintf('cs_bg = %.3f m/s\n', cfg_preview.cs_bg);
fprintf('M = %.2f\n', feat_cfg_preview.M);
fprintf('cs_guess_used = %.3f m/s\n', feat_cfg_preview.cs_guess_used);
fprintf('lambda_guess = %.3f cm\n', feat_cfg_preview.lambda_guess_used * 100);
fprintf('win_size = %d pixels\n', feat_cfg_preview.win_size);
fprintf('patch size = %.3f cm\n', feat_cfg_preview.win_size * cfg_preview.dx * 100);
fprintf('k0_true = %.3f rad/m\n', req_preview.k0_true);
fprintf('Nbins = %s\n', string(req_preview.Nbins));
fprintf('smooth_sigma = %.4f\n', req_preview.smooth_sigma);

%% ========================================================================
% 6. Save configuration preview
% ========================================================================

SAVE = struct();

SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
SAVE.run_name = sprintf('%s_%s', CFG.EXP.name, CFG.EXP.timestamp);
SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

make_dir_if_needed(SAVE.output_dir);
make_dir_if_needed(SAVE.table_dir);
make_dir_if_needed(SAVE.data_dir);

save(fullfile(SAVE.data_dir, 'test_04_mc_config_preview.mat'), ...
    'CFG', ...
    'CFG_base', ...
    'T_axes', ...
    'T_conditions', ...
    'cfg_preview', ...
    'feat_cfg_preview', ...
    'req_preview', ...
    'expected_total_rows', ...
    '-v7.3');

writetable(T_axes, fullfile(SAVE.table_dir, ...
    'test_04_sweep_axes.csv'));

writetable(T_conditions, fullfile(SAVE.table_dir, ...
    'test_04_sweep_conditions.csv'));

fprintf('\nSaved config preview to:\n%s\n', ...
    fullfile(SAVE.data_dir, 'test_04_mc_config_preview.mat'));

fprintf('\nTest 04 config preview completed successfully.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

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