%% analyze_test_12_level04_map_prediction_sweep.m
% Test 12 Level 04: map-level adaptive-q REQ validation.
%
% This script runs new homogeneous simulations over frequency, REQ window
% size, true background SWS, and aperture. It applies trained Test 12
% adaptive-q model to produce SWS maps, error maps, q maps, center-ROI MAPE,
% and center-ROI CoV.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLegendFontSize', 11);

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = adaptive_req.config.load_profile_config( ...
    'test_12_level04_map_prediction_sweep', ...
    'RootDir', root_dir);

%% Locate Test 12 run

PATHS12 = locate_latest_test12_paths(root_dir);
fprintf('\nUsing Test 12 run folder without loading full MC table:\n%s\n', ...
    PATHS12.run_dir);

%% Output folders

analysis_tag = sprintf('level_04_map_prediction_sweep_%s_%s', ...
    sanitize_filename(CFG.MAP04.MODEL.ModelName), ...
    sanitize_filename(CFG.MAP04.MODEL.FeatureSet));
analysis_dir = fullfile(PATHS12.analysis_dir, ...
    analysis_tag);
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
data_dir = fullfile(analysis_dir, 'data');
condition_dir = fullfile(data_dir, 'conditions');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(data_dir);
make_dir_if_needed(condition_dir);

%% Load trained Level 01 models

MODEL_SPECS = load_level01_deployment_models(PATHS12, CFG);

checkpoint_file = fullfile(data_dir, ...
    'level12_level04_map_prediction_checkpoint.mat');

%% Build condition table

T_conditions = build_condition_table(CFG);

T_all_pred = table();
T_all_roi = table();
T_all_overall = table();
completed_conditions = [];

if exist(checkpoint_file, 'file')
    fprintf('\nLoading Level 04 checkpoint:\n%s\n', checkpoint_file);
    S = load(checkpoint_file, 'T_all_pred', 'T_all_roi', ...
        'T_all_overall', 'completed_conditions');
    T_all_pred = S.T_all_pred;
    T_all_roi = S.T_all_roi;
    T_all_overall = S.T_all_overall;
    completed_conditions = S.completed_conditions;
    fprintf('Completed conditions in checkpoint: %d\n', ...
        numel(completed_conditions));
end

%% Run map sweep

for ci = 1:height(T_conditions)
    cond = T_conditions(ci, :);

    if any(completed_conditions == cond.condition_id)
        fprintf('Skipping completed condition %03d.\n', cond.condition_id);
        if CFG.MAP04.FIG.SavePerConditionMaps
            condition_file = fullfile(condition_dir, sprintf( ...
                'level12_level04_condition_%03d.mat', cond.condition_id));
            if exist(condition_file, 'file')
                S_cond = load(condition_file, 'T_pred_cond', 'cond');
                plot_condition_maps(S_cond.T_pred_cond, CFG, fig_dir, ...
                    S_cond.cond);
                fprintf('Replotted completed condition %03d.\n', ...
                    cond.condition_id);
            end
        end
        continue;
    end

    fprintf('\n=== Level 04 map condition %03d / %03d ===\n', ...
        ci, height(T_conditions));
    disp(cond(:, {'SIM_f0', 'REQ_M', 'SIM_cs_bg', ...
        'aperture_label', 'ConeHalfAngleDeg'}));

    cfg = build_sim_cfg(CFG, cond);
    [feat_cfg, req_options] = build_req_settings(CFG, cond);

    fprintf('Running homogeneous simulation with %d waves...\n', cfg.Nwaves);
    t_sim = tic;
    sim = adaptive_req.simulate.run_single_simulation(cfg);
    fprintf('Simulation completed in %.2f s.\n', toc(t_sim));

    fprintf('Extracting REQ map features with StepX=%d, StepZ=%d...\n', ...
        CFG.MAP04.REQ.StepX, CFG.MAP04.REQ.StepZ);
    t_feat = tic;
    [T_feat, req_out, global_req] = extract_map_feature_table( ...
        sim, cfg, feat_cfg, req_options, cond, CFG);
    fprintf('REQ/features completed in %.2f s for %d windows.\n', ...
        toc(t_feat), height(T_feat));

    fprintf('Applying adaptive-q models to REQ mappings...\n');
    t_pred = tic;
    T_pred_cond = table();
    for mi = 1:numel(MODEL_SPECS)
        T_pred_i = predict_model_map(MODEL_SPECS(mi), T_feat, cfg);
        T_pred_cond = concat_tables(T_pred_cond, T_pred_i);
    end
    fprintf('Adaptive-q prediction completed in %.2f s.\n', toc(t_pred));

    T_pred_cond = add_error_metrics(T_pred_cond);
    T_roi_cond = summarize_roi_metrics(T_pred_cond);
    T_overall_cond = summarize_roi_metrics(add_all_roi(T_pred_cond));

    T_all_pred = concat_tables(T_all_pred, T_pred_cond);
    T_all_roi = concat_tables(T_all_roi, T_roi_cond);
    T_all_overall = concat_tables(T_all_overall, T_overall_cond);

    save(fullfile(condition_dir, sprintf( ...
        'level12_level04_condition_%03d.mat', cond.condition_id)), ...
        'CFG', 'cfg', 'cond', 'T_feat', 'T_pred_cond', ...
        'T_roi_cond', 'T_overall_cond', 'req_out', 'global_req', '-v7.3');

    if CFG.MAP04.FIG.SavePerConditionMaps
        plot_condition_maps(T_pred_cond, CFG, fig_dir, cond);
    end

    completed_conditions(end + 1, 1) = cond.condition_id; %#ok<SAGROW>
    save(checkpoint_file, 'T_all_pred', 'T_all_roi', ...
        'T_all_overall', 'completed_conditions', 'T_conditions', ...
        'MODEL_SPECS', 'CFG', '-v7.3');
end

%% Save combined outputs

writetable(remove_cell_columns(T_all_pred), fullfile(table_dir, ...
    'level12_level04_map_predictions.csv'));
writetable(T_all_roi, fullfile(table_dir, ...
    'level12_level04_center_roi_metrics.csv'));
writetable(T_all_overall, fullfile(table_dir, ...
    'level12_level04_overall_map_metrics.csv'));
writetable(T_conditions, fullfile(table_dir, ...
    'level12_level04_condition_table.csv'));

save(fullfile(data_dir, 'level12_level04_map_prediction_sweep.mat'), ...
    'CFG', 'MODEL_SPECS', 'T_conditions', 'T_all_pred', ...
    'T_all_roi', 'T_all_overall', '-v7.3');

%% Summary figures

if CFG.MAP04.FIG.SaveSummaryFigures
    plot_summary_figures(T_all_roi, T_all_overall, CFG, fig_dir);
end

fprintf('\nTest 12 Level 04 map prediction sweep complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
fprintf('\nCenter ROI metrics preview.\n');
disp(sortrows(T_all_roi(:, {'model_name', 'SIM_f0', 'REQ_M', ...
    'SIM_cs_bg', 'aperture_label', 'MAPE_pct', 'CoV_pct'}), ...
    {'model_name', 'MAPE_pct'}));

%% Local functions

function MODEL_SPECS = load_level01_deployment_models(PATHS12, CFG)

deployment_dir = fullfile(PATHS12.analysis_dir, ...
    'level_01_model_comparison', 'models', 'deployment');

requested_names = string(CFG.MAP04.MODEL.ModelName);
if CFG.MAP04.MODEL.IncludeLocalOnly
    requested_names = [requested_names; "LocalOnly"];
end

MODEL_SPECS = struct([]);
for i = 1:numel(requested_names)
    try
        [MODEL_i, INFO_i, file_i] = adaptive_req.analysis.load_q_model_deployment( ...
            deployment_dir, ...
            'ModelName', requested_names(i), ...
            'FeatureSet', CFG.MAP04.MODEL.FeatureSet, ...
            'ModelType', CFG.MAP04.MODEL.ModelType);
    catch ME
        error(['Could not load the deployment model for %s | %s | %s.\n', ...
            'Expected deployment directory:\n%s\n\n', ...
            'Run the updated analyze_test_12_level01_model_comparison.m ', ...
            'to export lightweight deployment model files.\n\nOriginal error:\n%s'], ...
            requested_names(i), CFG.MAP04.MODEL.FeatureSet, ...
            CFG.MAP04.MODEL.ModelType, deployment_dir, ME.message);
    end

    MODEL_SPECS(i).model_name = string(INFO_i.model_name);
    MODEL_SPECS(i).feature_set = string(INFO_i.feature_set);
    MODEL_SPECS(i).model_type = string(INFO_i.model_type);
    MODEL_SPECS(i).model_file = string(file_i);
    MODEL_SPECS(i).model = MODEL_i;
end

fprintf('\nLoaded deployment models:\n');
disp(struct2table(rmfield(MODEL_SPECS, 'model')));

end

function PATHS = locate_latest_test12_paths(root_dir)

experiment_name = 'test_12_cs_guess_window_sweep';
output_root = fullfile(root_dir, 'outputs', experiment_name);

if ~exist(output_root, 'dir')
    error('Output root not found: %s', output_root);
end

D = dir(output_root);
D = D([D.isdir]);
names = string({D.name});
D = D(names ~= "." & names ~= "..");

if isempty(D)
    error('No Test 12 run folders found in: %s', output_root);
end

[~, ord] = sort([D.datenum], 'descend');
run_dir = fullfile(output_root, D(ord(1)).name);

PATHS = struct();
PATHS.root_dir = root_dir;
PATHS.experiment_name = experiment_name;
PATHS.output_root = output_root;
PATHS.run_dir = run_dir;
PATHS.data_dir = fullfile(run_dir, 'data');
PATHS.table_dir = fullfile(run_dir, 'tables');
PATHS.figure_dir = fullfile(run_dir, 'figures');
PATHS.analysis_dir = fullfile(run_dir, 'analysis');

if ~exist(PATHS.analysis_dir, 'dir')
    mkdir(PATHS.analysis_dir);
end

end

function T = build_condition_table(CFG)

rows = struct([]);
idx = 0;

for ai = 1:numel(CFG.MAP04.ConeHalfAngleDeg)
    for ci = 1:numel(CFG.MAP04.cs_bg_list)
        for mi = 1:numel(CFG.MAP04.M_list)
            for fi = 1:numel(CFG.MAP04.f_list)
                idx = idx + 1;
                rows(idx).condition_id = idx;
                rows(idx).SIM_f0 = CFG.MAP04.f_list(fi);
                rows(idx).REQ_M = CFG.MAP04.M_list(mi);
                rows(idx).SIM_cs_bg = CFG.MAP04.cs_bg_list(ci);
                rows(idx).aperture_idx = ai;
                rows(idx).aperture_label = CFG.MAP04.aperture_labels(ai);
                rows(idx).ConeHalfAngleDeg = CFG.MAP04.ConeHalfAngleDeg(ai);
            end
        end
    end
end

T = struct2table(rows);

end

function cfg = build_sim_cfg(CFG, cond)

S = CFG.MAP04.SIM;
cfg = adaptive_req.config.default_sim_config( ...
    'WaveModel', S.WaveModel, ...
    'f0', cond.SIM_f0, ...
    'cs_bg', cond.SIM_cs_bg, ...
    'Nwaves', S.Nwaves, ...
    'SNR', S.SNR, ...
    'AmpJitter', S.AmpJitter, ...
    'Lx', S.Lx, ...
    'Lz', S.Lz, ...
    'dx', S.dx, ...
    'dz', S.dz, ...
    'SourceSampling', S.SourceSampling, ...
    'AngularSamplingMethod', S.AngularSamplingMethod, ...
    'ConeAxis', S.ConeAxis, ...
    'ConeHalfAngleDeg', cond.ConeHalfAngleDeg, ...
    'ForceInPlaneWave', S.ForceInPlaneWave, ...
    'UseParfor', S.UseParfor, ...
    'Seed', S.SeedBase + cond.condition_id - 1);

end

function [feat_cfg, req_options] = build_req_settings(CFG, cond)

R = CFG.MAP04.REQ;
feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', cond.REQ_M, ...
    'cs_guess', R.cs_guess, ...
    'gamma_win', R.gamma_win, ...
    'pad_factor', R.pad_factor);

req_options = { ...
    'Nbins', R.Nbins, ...
    'Nbins_auto_oversample', R.Nbins_auto_oversample, ...
    'Nbins_min', R.Nbins_min, ...
    'smooth_sigma', R.smooth_sigma};

end

function [T, OUT, global_req] = extract_map_feature_table( ...
    sim, cfg, feat_cfg, req_options, cond, CFG)

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

% Avoid nested tiny parfor calls in the per-window radial average. The
% simulation can still use parfor, but local REQ windows are faster serially.
cfg_req = cfg;
cfg_req.UseParfor = false;

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
        sim.Uxz, cfg_req, req_cfg, feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features(global_curve);

global_req = struct();
global_req.q = q_global;
global_req.mapping = adaptive_req.quantile.make_req_mapping(global_curve);
global_req.features = global_features;
global_req.shape_features = global_shape;

OUT = adaptive_req.estimators.req_estimator_map( ...
    sim.Uxz, cfg_req, feat_cfg, ...
    'StepX', CFG.MAP04.REQ.StepX, ...
    'StepZ', CFG.MAP04.REQ.StepZ, ...
    'EdgeMode', CFG.MAP04.REQ.EdgeMode, ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', true, ...
    'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, ...
    'UseWindowParfor', true, ...
    'StoreReqCurves', false, ...
    'Verbose', false);

T = OUT.feature_table;
T.condition_id = cond.condition_id * ones(height(T), 1);
T.condition_label = repmat(string(sprintf( ...
    'f%d_M%g_cs%g_%s', cond.SIM_f0, cond.REQ_M, ...
    cond.SIM_cs_bg, cond.aperture_label)), height(T), 1);
T.step_idx = cond.aperture_idx * ones(height(T), 1);
T.realization_idx = ones(height(T), 1);
T.frequency_hz = cond.SIM_f0 * ones(height(T), 1);
T.aperture_label = repmat(string(cond.aperture_label), height(T), 1);
T.ConeHalfAngleDeg = cond.ConeHalfAngleDeg * ones(height(T), 1);
T.Omega_sr = cone_omega_sr(cond.ConeHalfAngleDeg) * ones(height(T), 1);
T.omega_mean = T.Omega_sr;
T.SIM_Nwaves = cfg.Nwaves * ones(height(T), 1);
T.REQ_cs_guess = feat_cfg.cs_guess_used * ones(height(T), 1);
T.M_eff_guess = feat_cfg.M * ones(height(T), 1);
T.M_eff_true_diag = feat_cfg.M * feat_cfg.cs_guess_used ./ ...
    cond.SIM_cs_bg * ones(height(T), 1);
T.lambda_guess = feat_cfg.cs_guess_used / cond.SIM_f0 * ones(height(T), 1);
T.lambda_true = cond.SIM_cs_bg / cond.SIM_f0 * ones(height(T), 1);
T.cs_true = cond.SIM_cs_bg * ones(height(T), 1);
T.roi_name = classify_center_roi(T.x_center_m, T.z_center_m, CFG.MAP04.ROI);
T.global_req_mapping = repmat({global_req.mapping}, height(T), 1);
T.q_global_req = q_global * ones(height(T), 1);
T.global_REQ_Nbins_effective = global_curve.Nbins_effective * ones(height(T), 1);

T = assign_global_feature_columns(T, global_features, "global_");
T = assign_global_feature_columns(T, global_shape, "global_");

end

function T_pred = predict_model_map(MODEL_SPEC, T_feat, cfg)

T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_SPEC.model, T_feat, ...
    'ModelType', MODEL_SPEC.model_type, ...
    'ModelName', MODEL_SPEC.model_name);

T_pred = make_prediction_base(T_feat);
T_pred.model_name = repmat(MODEL_SPEC.model_name, height(T_feat), 1);
T_pred.feature_set = repmat(MODEL_SPEC.feature_set, height(T_feat), 1);
T_pred.model_type = repmat(MODEL_SPEC.model_type, height(T_feat), 1);
T_pred.model_role = repmat("operational_map", height(T_feat), 1);
T_pred.q_pred = min(max(T_q.q_pred, 0.001), 0.999);
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, cfg.f0);

end

function T = make_prediction_base(T_feat)

vars = ["condition_id", "condition_label", "step_idx", ...
    "realization_idx", "frequency_hz", "SIM_f0", "SIM_cs_bg", ...
    "SIM_WaveModel", "SIM_Nwaves", "REQ_M", "REQ_cs_guess", ...
    "M_eff_guess", "M_eff_true_diag", "aperture_label", ...
    "ConeHalfAngleDeg", "Omega_sr", "patch_idx", "map_iz", ...
    "map_ix", "x_center_m", "z_center_m", "roi_name", "cs_true"];
vars = vars(ismember(vars, string(T_feat.Properties.VariableNames)));
T = T_feat(:, cellstr(vars));

end

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)

cs_pred = nan(numel(q_pred), 1);
for i = 1:numel(q_pred)
    if isempty(mappings{i}) || ~isfinite(q_pred(i))
        continue;
    end
    cs_pred(i) = adaptive_req.quantile.quantile_to_cs( ...
        mappings{i}, q_pred(i), f0);
end

end

function T = add_error_metrics(T)

T.cs_error = T.cs_pred - T.cs_true;
T.cs_abs_error = abs(T.cs_error);
T.cs_error_pct = 100 * T.cs_error ./ T.cs_true;
T.cs_abs_error_pct = abs(T.cs_error_pct);

end

function T_roi = summarize_roi_metrics(T)

T = T(T.roi_name ~= "outside_roi", :);
if isempty(T)
    T_roi = table();
    return;
end

[G, T_roi] = findgroups(T(:, {'condition_id', 'model_name', ...
    'feature_set', 'model_type', 'roi_name', 'SIM_f0', 'REQ_M', ...
    'SIM_cs_bg', 'REQ_cs_guess', 'M_eff_guess', 'M_eff_true_diag', ...
    'aperture_label', 'ConeHalfAngleDeg'}));
T_roi.n = splitapply(@numel, T.cs_error_pct, G);
T_roi.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), ...
    T.cs_error_pct, G);
T_roi.bias_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.CoV_pct = splitapply(@cov_pct, T.cs_pred, G);
T_roi.cs_true = splitapply(@(x) median(x, 'omitnan'), T.cs_true, G);
T_roi.cs_pred_mean = splitapply(@(x) mean(x, 'omitnan'), T.cs_pred, G);
T_roi.cs_pred_median = splitapply(@(x) median(x, 'omitnan'), T.cs_pred, G);
T_roi.cs_pred_std = splitapply(@(x) std(x, 'omitnan'), T.cs_pred, G);

end

function T = add_all_roi(T)

T.roi_name = repmat("all_valid_map", height(T), 1);

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function roi = classify_center_roi(x, z, roi_spec)

roi = repmat("outside_roi", numel(x), 1);
in_roi = x >= roi_spec.xlim_m(1) & x <= roi_spec.xlim_m(2) & ...
    z >= roi_spec.zlim_m(1) & z <= roi_spec.zlim_m(2);
roi(in_roi) = string(roi_spec.name);

end

function T = assign_global_feature_columns(T, values, prefix)

if ~isstruct(values)
    return;
end

names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if isnumeric(value_i) && isscalar(value_i)
        T.(char(string(prefix) + string(names{i}))) = ...
            double(value_i) * ones(height(T), 1);
    elseif islogical(value_i) && isscalar(value_i)
        T.(char(string(prefix) + string(names{i}))) = ...
            double(value_i) * ones(height(T), 1);
    end
end

end

function omega = cone_omega_sr(half_angle_deg)

theta = deg2rad(half_angle_deg);
omega = 2*pi*(1 - cos(theta));

end

function plot_condition_maps(T_pred, CFG, fig_dir, cond)

models = unique(T_pred.model_name, 'stable');
for mi = 1:numel(models)
    Tm = T_pred(T_pred.model_name == models(mi), :);

    cs_map = map_from_table(Tm, 'cs_pred');
    err_map = map_from_table(Tm, 'cs_error_pct');
    q_map = map_from_table(Tm, 'q_pred');
    x_cm = unique(Tm.x_center_m, 'stable') * 100;
    z_cm = unique(Tm.z_center_m, 'stable') * 100;

    base = sprintf('level12_level04_cond%03d_%s_%s_%s', ...
        cond.condition_id, sanitize_filename(models(mi)), ...
        sanitize_filename(Tm.feature_set(1)), ...
        sanitize_filename(cond.aperture_label));

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 30 11]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl);
    imagesc(ax, x_cm, z_cm, cs_map);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal', 'FontSize', 11);
    colorbar(ax);
    clim(ax, sws_color_limits(CFG.MAP04.cs_bg_list));
    title(ax, 'Predicted SWS', 'FontSize', 12, 'FontWeight', 'normal');
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    draw_roi(ax, CFG.MAP04.ROI);

    ax = nexttile(tl);
    imagesc(ax, x_cm, z_cm, err_map);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal', 'FontSize', 11);
    colorbar(ax);
    clim(ax, CFG.MAP04.FIG.ErrorLimitsPct);
    title(ax, 'SWS error (%)', 'FontSize', 12, 'FontWeight', 'normal');
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    draw_roi(ax, CFG.MAP04.ROI);

    ax = nexttile(tl);
    imagesc(ax, x_cm, z_cm, q_map);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal', 'FontSize', 11);
    colorbar(ax);
    clim(ax, [0 1]);
    title(ax, 'Predicted q', 'FontSize', 12, 'FontWeight', 'normal');
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    draw_roi(ax, CFG.MAP04.ROI);

    title(tl, sprintf('%s | %s | f=%d Hz | M=%g | c_s=%g | %s', ...
        models(mi), Tm.feature_set(1), cond.SIM_f0, cond.REQ_M, ...
        cond.SIM_cs_bg, cond.aperture_label), ...
        'Interpreter', 'none', 'FontSize', 13, 'FontWeight', 'normal');

    exportgraphics(fig, fullfile(fig_dir, base + "_maps.png"), ...
        'Resolution', 220, 'BackgroundColor', 'white');
    close(fig);
end

end

function A = map_from_table(T, value_var)

nz = max(T.map_iz);
nx = max(T.map_ix);
A = nan(nz, nx);
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
A(idx) = T.(value_var);

end

function draw_roi(ax, roi)

hold(ax, 'on');
x0 = roi.xlim_m(1) * 100;
z0 = roi.zlim_m(1) * 100;
w = diff(roi.xlim_m) * 100;
h = diff(roi.zlim_m) * 100;
rectangle(ax, 'Position', [x0 z0 w h], ...
    'EdgeColor', 'w', 'LineWidth', 1.5);
hold(ax, 'off');

end

function lim = sws_color_limits(cs_values)

cs_values = double(cs_values(:));
lo = min(cs_values);
hi = max(cs_values);

if ~isfinite(lo) || ~isfinite(hi)
    lim = [0 1];
elseif hi > lo
    lim = [lo hi];
else
    pad = max(0.25, 0.15 * abs(lo));
    lim = [lo - pad, hi + pad];
end

end

function plot_summary_figures(T_roi, T_all, CFG, fig_dir)

plot_metric_heatmaps(T_roi, 'MAPE_pct', CFG, fig_dir, ...
    'level12_level04_center_roi_mape_heatmaps.png', ...
    'Center ROI MAPE (%)');
plot_metric_heatmaps(T_roi, 'CoV_pct', CFG, fig_dir, ...
    'level12_level04_center_roi_cov_heatmaps.png', ...
    'Center ROI CoV (%)');
plot_metric_heatmaps(T_all, 'MAPE_pct', CFG, fig_dir, ...
    'level12_level04_overall_mape_heatmaps.png', ...
    'Overall map MAPE (%)');

plot_metric_box(T_roi, 'MAPE_pct', fig_dir, ...
    'level12_level04_center_roi_mape_box_by_params.png');
plot_metric_box(T_roi, 'CoV_pct', fig_dir, ...
    'level12_level04_center_roi_cov_box_by_params.png');
plot_frequency_curves(T_roi, fig_dir);

end

function plot_metric_heatmaps(T, metric_var, ~, fig_dir, file_name, title_text)

models = unique(T.model_name, 'stable');
apertures = unique(T.aperture_label, 'stable');
cs_values = unique(T.SIM_cs_bg, 'stable');

for mi = 1:numel(models)
    for ai = 1:numel(apertures)
        fig = figure('Color', 'w', 'Units', 'centimeters', ...
            'Position', [2 2 24 18]);
        tl = tiledlayout(2, 2, 'TileSpacing', 'compact', ...
            'Padding', 'compact');
        for ci = 1:numel(cs_values)
            ax = nexttile(tl);
            Ti = T(T.model_name == models(mi) & ...
                T.aperture_label == apertures(ai) & ...
                T.SIM_cs_bg == cs_values(ci), :);
            M_values = unique(T.REQ_M, 'stable');
            f_values = unique(T.SIM_f0, 'stable');
            Z = nan(numel(M_values), numel(f_values));
            for ii = 1:numel(M_values)
                for jj = 1:numel(f_values)
                    idx = Ti.REQ_M == M_values(ii) & Ti.SIM_f0 == f_values(jj);
                    if any(idx)
                        Z(ii, jj) = mean(Ti.(metric_var)(idx), 'omitnan');
                    end
                end
            end
            imagesc(ax, f_values, M_values, Z);
            set(ax, 'YDir', 'normal');
            colorbar(ax);
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, 'REQ M');
            title(ax, sprintf('c_s = %g m/s', cs_values(ci)));
        end
        title(tl, sprintf('%s | %s | aperture=%s', ...
            title_text, models(mi), apertures(ai)), 'Interpreter', 'none');
        exportgraphics(fig, fullfile(fig_dir, ...
            replace(file_name, ".png", "_" + sanitize_filename(models(mi)) + ...
            "_" + sanitize_filename(apertures(ai)) + ".png")), ...
            'Resolution', 250, 'BackgroundColor', 'white');
        close(fig);
    end
end

end

function plot_metric_box(T, metric_var, fig_dir, file_name)

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 28 12]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
boxchart(categorical(string(T.REQ_M)), T.(metric_var), ...
    'GroupByColor', categorical(T.model_name));
xlabel('REQ M');
ylabel(strrep(metric_var, '_', '\_'));
title('By window size');
grid on;

nexttile;
boxchart(categorical(string(T.SIM_cs_bg)), T.(metric_var), ...
    'GroupByColor', categorical(T.model_name));
xlabel('True c_s (m/s)');
ylabel(strrep(metric_var, '_', '\_'));
title('By true SWS');
grid on;

nexttile;
boxchart(categorical(T.aperture_label), T.(metric_var), ...
    'GroupByColor', categorical(T.model_name));
xlabel('Aperture');
ylabel(strrep(metric_var, '_', '\_'));
title('By aperture');
grid on;
legend('Location', 'best', 'Interpreter', 'none');

exportgraphics(fig, fullfile(fig_dir, file_name), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_frequency_curves(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 28 16]);
tl = tiledlayout(numel(models), 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for mi = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(mi), :);
    groups = findgroups(Ti.REQ_M, Ti.SIM_cs_bg, Ti.aperture_label);
    hold(ax, 'on');
    for gi = 1:max(groups)
        Tg = Ti(groups == gi, :);
        Tg = sortrows(Tg, 'SIM_f0');
        plot(ax, Tg.SIM_f0, Tg.MAPE_pct, '-o', ...
            'DisplayName', sprintf('M=%g c_s=%g %s', ...
            Tg.REQ_M(1), Tg.SIM_cs_bg(1), Tg.aperture_label(1)));
    end
    title(ax, sprintf('Center ROI MAPE vs frequency | %s', models(mi)), ...
        'Interpreter', 'none');
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'MAPE (%)');
    grid(ax, 'on');
end

exportgraphics(fig, fullfile(fig_dir, ...
    'level12_level04_center_roi_mape_frequency_curves.png'), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

end

function T = concat_tables(A, B)

if isempty(A)
    T = B;
    return;
end
if isempty(B)
    T = A;
    return;
end

vars_all = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');
A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);
T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];

end

function T = add_missing_columns(T, vars_all)

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    if ismember(vars_all(i), vars)
        continue;
    end
    name_i = char(vars_all(i));
    string_like = any(endsWith(vars_all(i), ...
        ["name", "label", "type", "role", "set"]));
    if string_like
        T.(name_i) = strings(height(T), 1);
    else
        T.(name_i) = nan(height(T), 1);
    end
end

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end

function name = sanitize_filename(x)

name = regexprep(string(x), '[^A-Za-z0-9_]+', '_');
name = matlab.lang.makeValidName(char(name));
name = string(name);

end
