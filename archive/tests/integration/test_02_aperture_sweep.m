%% test_02_aperture_sweep.m
% Test the full aperture sweep workflow for adaptive_req.
%
% This script verifies that:
%   1. The project setup works.
%   2. Simulation, feature, and REQ settings can be defined cleanly.
%   3. The aperture sweep runs through adaptive_req.studies.run_aperture_sweep.
%   4. The output table contains q_theory, features, and speed estimates.
%   5. Summary plots and numerical outputs can be generated.

clear; clc; close all;
format compact;

%% ========================================================================
% 1. Locate project root and set up package path
% ========================================================================
% This script is assumed to live in:
%
%   adaptive_req/tests/integration/
%
% Therefore:
%   this_dir -> adaptive_req/tests/integration
%   root_dir -> adaptive_req
%
% setup_adaptive_req adds:
%
%   adaptive_req/src/
%
% to the MATLAB path.

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% ========================================================================
% 2. Experiment settings
% ========================================================================

EXP = struct();

EXP.name = 'test_02_aperture_sweep';
EXP.sampling_mode = 'cone';

EXP.num_steps = 10;
EXP.num_realizations = 1;
EXP.num_patches = 3;

EXP.step_indices = [];

EXP.seed_base = 1000;

EXP.selected_step = 1;
EXP.selected_patch = 1;

EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));

%% ========================================================================
% 3. Output settings
% ========================================================================

SAVE = struct();

SAVE.root_dir = fullfile(root_dir, 'outputs', EXP.name);
SAVE.run_name = sprintf('%s_%s_%s', ...
    EXP.name, EXP.sampling_mode, EXP.timestamp);

SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.figure_dir = fullfile(SAVE.output_dir, 'figures');
SAVE.step_diag_dir = fullfile(SAVE.figure_dir, 'step_diagnostics');
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

SAVE.save_png = true;
SAVE.save_pdf = true;
SAVE.save_fig = false;
SAVE.png_resolution = 300;
SAVE.close_after_save = false;

make_dir_if_needed(SAVE.output_dir);
make_dir_if_needed(SAVE.figure_dir);
make_dir_if_needed(SAVE.step_diag_dir);
make_dir_if_needed(SAVE.table_dir);
make_dir_if_needed(SAVE.data_dir);

%% ========================================================================
% 4. Plot switches
% ========================================================================
% Step diagnostics are generated inside the aperture sweep.
% These do not require StoreWavefields = true because the wavefield exists
% in memory during each step.

PLOT = struct();

PLOT.show_step_diagnostics = false;
PLOT.save_step_diagnostics = false;
PLOT.step_diagnostic_visible = false;
PLOT.close_step_diagnostics_after_save = true;

% Summary plots generated after the sweep.
PLOT.show_q_vs_aperture = true;
PLOT.show_feature_space = true;
PLOT.show_feature_vs_q = true;
PLOT.show_feature_grid = true;

% These require storing all wavefields in sweep.sims.
% Keep them false for a cleaner and lighter integration test.
PLOT.show_selected_step_wavefield = false;
PLOT.show_all_step_wavefields = false;

PLOT.store_wavefields = ...
    PLOT.show_selected_step_wavefield || PLOT.show_all_step_wavefields;

PLOT.save_summary_figures = false;

%% ========================================================================
% 5. Simulation configuration
% ========================================================================

cfg = adaptive_req.config.default_sim_config();

%% ========================================================================
% 6. Feature and REQ configuration
% ========================================================================
% M controls the local window size:
%
%   L_patch approximately equals M * lambda_guess
%
% Nbins controls the radial discretization used to compute the cumulative
% radial energy curve.

REQ = struct();

REQ.M = 3;
REQ.cs_guess = 3.0;
REQ.gamma_win = 1.0;
REQ.pad_factor = 2.0;

REQ.Nbins = 'auto';
REQ.Nbins_auto_oversample = 1;
REQ.Nbins_min = 16;

REQ.smooth_sigma = 0.001;

REQ.use_donut = false;
REQ.donut_cs_min = 1.0;
REQ.donut_cs_max = 5.0;
REQ.donut_taper_rel = 0.06;
REQ.apply_donut_to_final_map = false;

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

fprintf('\nREQ configuration preview.\n');
fprintf('M = %.2f\n', feat_cfg.M);
fprintf('cs_guess_used = %.3f m/s\n', feat_cfg.cs_guess_used);
fprintf('lambda_guess = %.3f cm\n', feat_cfg.lambda_guess_used * 100);
fprintf('win_size = %d pixels\n', feat_cfg.win_size);
fprintf('patch size = %.3f cm\n', feat_cfg.win_size * cfg.dx * 100);
fprintf('Nbins requested = %s\n', string(req_preview.Nbins));
fprintf('Nbins auto oversample = %.2f\n', req_preview.Nbins_auto_oversample);
fprintf('Nbins min = %d\n', req_preview.Nbins_min);
fprintf('smooth_sigma = %.4f\n', req_preview.smooth_sigma);

%% ========================================================================
% 7. Run aperture sweep
% ========================================================================

[T_raw, sweep] = adaptive_req.studies.run_aperture_sweep( ...
    cfg, feat_cfg, ...
    'SamplingMode', EXP.sampling_mode, ...
    'NumSteps', EXP.num_steps, ...
    'StepIndices', EXP.step_indices, ...
    'NumRealizations', EXP.num_realizations, ...
    'NumPatches', EXP.num_patches, ...
    'SeedBase', EXP.seed_base, ...
    'ReqOptions', req_options, ...
    'StoreWavefields', PLOT.store_wavefields, ...
    'StoreReqCurve', true, ...
    'StoreFeatureStruct', false, ...
    'PlotStepDiagnostics', PLOT.show_step_diagnostics || PLOT.save_step_diagnostics, ...
    'SaveStepDiagnostics', PLOT.save_step_diagnostics, ...
    'StepDiagnosticPatchIndex', EXP.selected_patch, ...
    'StepDiagnosticDir', SAVE.step_diag_dir, ...
    'StepDiagnosticVisible', PLOT.step_diagnostic_visible, ...
    'CloseStepDiagnosticsAfterSave', PLOT.close_step_diagnostics_after_save, ...
    'SaveDiagnosticPNG', SAVE.save_png, ...
    'SaveDiagnosticPDF', SAVE.save_pdf, ...
    'SaveDiagnosticFIG', SAVE.save_fig, ...
    'DiagnosticResolution', SAVE.png_resolution, ...
    'Verbose', true);

%% ========================================================================
% 8. Display output table summary
% ========================================================================

fprintf('\nOutput table summary.\n');
fprintf('Number of rows = %d\n', height(T_raw));
fprintf('Number of columns = %d\n', width(T_raw));

expected_rows = EXP.num_steps * EXP.num_realizations * EXP.num_patches;

fprintf('Expected rows = %d\n', expected_rows);

if height(T_raw) ~= expected_rows
    warning('Unexpected number of rows in T_raw.');
end

display_vars = { ...
    'step_idx', ...
    'realization', ...
    'patch_idx', ...
    'Omega_sr', ...
    'q_theory', ...
    'cs_true', ...
    'cs_req_q_theory', ...
    'radial_entropy', ...
    'ang_entropy', ...
    'Nbins_effective'};

display_vars = display_vars(ismember(display_vars, T_raw.Properties.VariableNames));

disp(T_raw(:, display_vars));

%% ========================================================================
% 9. Summary plots
% ========================================================================

if PLOT.show_q_vs_aperture

    figQ = adaptive_req.figures.plot_q_vs_aperture( ...
        T_raw, sweep, ...
        'UseOmega', true, ...
        'ShowRaw', true, ...
        'ShowErrorbar', true, ...
        'Title', 'Reference quantile vs aperture');

    save_summary_figure_if_requested( ...
        figQ, PLOT, SAVE, 'q_vs_aperture');
end

if PLOT.show_feature_space

    figFeatureSpace = adaptive_req.figures.plot_feature_space( ...
        T_raw, ...
        'XFeature', 'ang_entropy', ...
        'YFeature', 'radial_entropy', ...
        'ColorVariable', 'q_theory', ...
        'Title', 'Feature space colored by q_{theory}');

    save_summary_figure_if_requested( ...
        figFeatureSpace, PLOT, SAVE, 'feature_space_q_theory');
end

if PLOT.show_feature_vs_q

    figAng = adaptive_req.figures.plot_feature_vs_q( ...
        T_raw, ...
        'ang_entropy', ...
        'ColorBy', 'Omega_sr');

    figRad = adaptive_req.figures.plot_feature_vs_q( ...
        T_raw, ...
        'radial_entropy', ...
        'ColorBy', 'Omega_sr');

    save_summary_figure_if_requested( ...
        figAng, PLOT, SAVE, 'ang_entropy_vs_q');

    save_summary_figure_if_requested( ...
        figRad, PLOT, SAVE, 'radial_entropy_vs_q');
end

if PLOT.show_feature_grid

    [figFeatures, fit_table] = adaptive_req.figures.plot_features_vs_q_grid( ...
        T_raw, ...
        'ColorBy', 'Omega_sr');

    fprintf('\nFeature vs q fit summary.\n');
    disp(fit_table);

    save_summary_figure_if_requested( ...
        figFeatures, PLOT, SAVE, 'features_vs_q_grid');
end

%% ========================================================================
% 10. Optional post-sweep wavefield plots
% ========================================================================

if PLOT.show_selected_step_wavefield

    if PLOT.store_wavefields

        selected_step = min(EXP.selected_step, sweep.num_steps);

        adaptive_req.figures.plot_step_wavefield( ...
            sweep, selected_step, ...
            'Realization', 1, ...
            'ShowPatches', true, ...
            'NumPatches', EXP.num_patches);

    else

        warning(['Selected-step wavefield plot requested, but ', ...
                 'PLOT.store_wavefields is false.']);

    end
end

if PLOT.show_all_step_wavefields

    if PLOT.store_wavefields

        adaptive_req.figures.plot_sweep_wavefields( ...
            sweep, ...
            'Realization', 1);

    else

        warning(['All-step wavefield plot requested, but ', ...
                 'PLOT.store_wavefields is false.']);

    end
end

%% ========================================================================
% 11. Save numerical results
% ========================================================================

save(fullfile(SAVE.data_dir, 'test_02_aperture_sweep_results.mat'), ...
    'T_raw', ...
    'sweep', ...
    'EXP', ...
    'REQ', ...
    'PLOT', ...
    'SAVE', ...
    'cfg', ...
    'feat_cfg', ...
    'req_preview', ...
    'req_options', ...
    '-v7.3');

T_csv = T_raw;

vars_to_remove = { ...
    'req_curve', ...
    'feat', ...
    'feature_struct'};

vars_to_remove = vars_to_remove(ismember(vars_to_remove, ...
    T_csv.Properties.VariableNames));

if ~isempty(vars_to_remove)
    T_csv(:, vars_to_remove) = [];
end

writetable(T_csv, fullfile(SAVE.table_dir, ...
    'test_02_aperture_sweep_table.csv'));

fprintf('\nSaved MAT output to:\n%s\n', ...
    fullfile(SAVE.data_dir, 'test_02_aperture_sweep_results.mat'));

fprintf('\nSaved CSV table to:\n%s\n', ...
    fullfile(SAVE.table_dir, 'test_02_aperture_sweep_table.csv'));

fprintf('\nTest 02 completed successfully.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function make_dir_if_needed(folder_path)

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

end

function save_summary_figure_if_requested(fig_handle, PLOT, SAVE, base_name)

if ~PLOT.save_summary_figures
    return;
end

if isempty(fig_handle) || ~isvalid(fig_handle)
    return;
end

adaptive_req.figures.save_figure_bundle( ...
    fig_handle, ...
    SAVE.figure_dir, ...
    base_name, ...
    'SavePNG', SAVE.save_png, ...
    'SavePDF', SAVE.save_pdf, ...
    'SaveFIG', SAVE.save_fig, ...
    'Resolution', SAVE.png_resolution, ...
    'CloseAfterSave', SAVE.close_after_save);

end