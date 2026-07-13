%% test_03_config_loader.m
% Config-driven aperture sweep using adaptive_req.config.load_run_config.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load resolved run configuration

RUN = adaptive_req.config.load_run_config( ...
    'test_03_baseline', ...
    'RootDir', root_dir);

print_run_summary(RUN);

%% Run aperture sweep

[T_raw, sweep] = adaptive_req.studies.run_aperture_sweep(RUN);

%% Basic output summary

fprintf('\nOutput table summary.\n');
fprintf('Number of rows = %d\n', height(T_raw));
fprintf('Number of columns = %d\n', width(T_raw));

expected_rows = ...
    RUN.EXP.num_steps * ...
    RUN.EXP.num_realizations * ...
    RUN.EXP.num_patches;

fprintf('Expected rows = %d\n', expected_rows);

%% Summary plots

FIGS = adaptive_req.figures.plot_run_summary(T_raw, sweep, RUN);

%% Save outputs

save(fullfile(RUN.SAVE.data_dir, 'test_03_results.mat'), ...
    'T_raw', ...
    'sweep', ...
    'RUN', ...
    'FIGS', ...
    '-v7.3');

T_csv = remove_heavy_table_columns(T_raw);

writetable(T_csv, fullfile(RUN.SAVE.table_dir, 'test_03_table.csv'));

fprintf('\nSaved MAT output to:\n%s\n', ...
    fullfile(RUN.SAVE.data_dir, 'test_03_results.mat'));

fprintf('\nSaved CSV table to:\n%s\n', ...
    fullfile(RUN.SAVE.table_dir, 'test_03_table.csv'));

fprintf('\nTest 03 completed successfully.\n');

%% Local helper functions

function print_run_summary(RUN)

fprintf('\nLoaded run configuration.\n');
fprintf('Profile = %s\n', RUN.EXP.profile_name);
fprintf('Experiment name = %s\n', RUN.EXP.name);
fprintf('Sampling mode = %s\n', RUN.EXP.sampling_mode);
fprintf('f0 = %.1f Hz\n', RUN.cfg.f0);
fprintf('cs_bg = %.3f m/s\n', RUN.cfg.cs_bg);
fprintf('M = %.2f\n', RUN.feat_cfg.M);
fprintf('win_size = %d pixels\n', RUN.feat_cfg.win_size);
fprintf('patch size = %.3f cm\n', RUN.feat_cfg.win_size * RUN.cfg.dx * 100);
fprintf('Nbins = %s\n', string(RUN.req_preview.Nbins));
fprintf('smooth_sigma = %.4f\n', RUN.req_preview.smooth_sigma);

end

function T = remove_heavy_table_columns(T)

vars_to_remove = { ...
    'req_curve', ...
    'feat', ...
    'feature_struct'};

vars_to_remove = vars_to_remove(ismember(vars_to_remove, ...
    T.Properties.VariableNames));

if ~isempty(vars_to_remove)
    T(:, vars_to_remove) = [];
end

end