%% analyze_test_11_level02_bilayer_transfer.m
% Test 11 Level 02: apply Test 11 models to a simple bilayer phantom.
%
% Three transfer tests are compared:
%   1. GlobalQSingle: one global quantile applied to each local REQ mapping.
%   2. LocalOnly:     Test 11 local-only model.
%   3. Hybrid:        Test 11 local + global context model.
%
% Patch extraction uses valid edge handling: patch centers are only placed
% where the full REQ window stays inside the simulated field.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

%% Locate Test 11 Level 01 outputs

[~, ~, PATHS11] = adaptive_req.analysis.load_mc_results( ...
    'test_11_global_req_features', ...
    'RootDir', root_dir, ...
    'Verbose', true);

level01_dir = fullfile(PATHS11.analysis_dir, ...
    'level_01_global_vs_local');
model_file = fullfile(level01_dir, 'models', ...
    'level11_global_vs_local_models.mat');
metrics_file = fullfile(level01_dir, 'tables', ...
    'level11_sws_metrics.csv');

if ~exist(model_file, 'file')
    error('Test 11 Level 01 model file not found:\n%s', model_file);
end

S = load(model_file, 'MODELS');
MODELS = S.MODELS;

MODEL_local = get_named_model(MODELS, "LocalOnly");
MODEL_global = get_named_model(MODELS, "GlobalOnly");
MODEL_hybrid = get_named_model(MODELS, "HybridLocalGlobal");

%% Output folders

analysis_dir = fullfile(PATHS11.analysis_dir, ...
    'level_02_bilayer_transfer');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
data_dir = fullfile(analysis_dir, 'data');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(data_dir);

%% Improvement summary from homogeneous Test 11

T_improvement = compute_improvement_summary(metrics_file);
writetable(T_improvement, fullfile(table_dir, ...
    'level11_level02_homogeneous_improvement_summary.csv'));

fprintf('\nHomogeneous Test 11 improvement summary.\n');
disp(T_improvement);

%% Bilayer simulation

cfg = adaptive_req.config.default_sim_config( ...
    'WaveModel', 'spherical', ...
    'f0', 400, ...
    'cs_bg', 2.0, ...
    'Nwaves', 2000, ...
    'SNR', Inf, ...
    'AmpJitter', 0, ...
    'Lx', 0.05, ...
    'Lz', 0.05, ...
    'dx', 1e-4, ...
    'dz', 1e-4, ...
    'SourceSampling', 'ranges', ...
    'AngularSamplingMethod', 'random', ...
    'ConeAxis', [-1 0 0], ...
    'ConeHalfAngleDeg', 180, ...
    'ForceInPlaneWave', false, ...
    'UseParfor', true, ...
    'Seed', 100);

cfg.MaskType = 'bilayer';
cfg.cs_inc = 3.0;
cfg.MaskParams = struct( ...
    'Bi_Angle', deg2rad(0), ...
    'Bi_Offset', 0.025, ...
    'SigmaEdge', 1e-6);

feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', 3, ...
    'cs_guess', 3.0, ...
    'gamma_win', 1.0, ...
    'pad_factor', 1.0);

req_options = { ...
    'Nbins', 'auto', ...
    'Nbins_auto_oversample', 1, ...
    'Nbins_min', 16, ...
    'smooth_sigma', 1};

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

fprintf('\nRunning bilayer simulation...\n');
sim = adaptive_req.simulate.run_single_simulation(cfg);

patch_pack = build_dense_valid_patch_pack(cfg, feat_cfg, ...
    'StepX', 7, ...
    'StepZ', 7);

roi_specs = build_roi_specs();

fprintf('\nExtracting local/global REQ features with valid edge mode...\n');
[T_feat, global_req] = extract_bilayer_feature_table( ...
    sim, cfg, feat_cfg, req_options, patch_pack, roi_specs);

%% Prediction tests

theory_field_type = 'Diffuse3D';
q_theory_discrete = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    cfg.dx, cfg.dz, cfg.f0, feat_cfg.cs_guess_used, ...
    'M', feat_cfg.M, ...
    'Gamma', feat_cfg.gamma_win, ...
    'PadFactor', feat_cfg.pad_factor, ...
    'Nbins', req_cfg.Nbins, ...
    'SmoothSigma', req_cfg.smooth_sigma, ...
    'TheoryMode', 'S2D', ...
    'FieldType', theory_field_type, ...
    'Plot', false);

T_pred_theory = predict_fixed_q_to_sws( ...
    T_feat, "TheoryDiffuse3D", q_theory_discrete.q_th, cfg.f0);
T_pred_global = predict_global_model_single_q( ...
    MODEL_global, T_feat, cfg);
T_pred_local = predict_model_to_sws( ...
    MODEL_local, T_feat, "LocalOnly", cfg.f0);
T_pred_hybrid = predict_model_to_sws( ...
    MODEL_hybrid, T_feat, "HybridLocalGlobal", cfg.f0);

T_pred = [T_pred_theory; T_pred_global; T_pred_local; T_pred_hybrid];
T_pred = add_error_metrics(T_pred);

T_roi = summarize_roi_metrics(T_pred);
T_overall = summarize_roi_metrics(add_all_roi(T_pred));

writetable(remove_cell_columns(T_feat), fullfile(table_dir, ...
    'level11_level02_bilayer_features.csv'));
writetable(T_pred, fullfile(table_dir, ...
    'level11_level02_bilayer_predictions.csv'));
writetable(T_roi, fullfile(table_dir, ...
    'level11_level02_bilayer_roi_metrics.csv'));
writetable(T_overall, fullfile(table_dir, ...
    'level11_level02_bilayer_overall_metrics.csv'));

save(fullfile(data_dir, 'level11_level02_bilayer_transfer.mat'), ...
    'cfg', 'feat_cfg', 'req_cfg', 'sim', 'patch_pack', 'roi_specs', ...
    'global_req', 'T_feat', 'T_pred', 'T_roi', 'T_overall', ...
    'T_improvement', 'q_theory_discrete', 'theory_field_type', '-v7.3');

%% Figures

plot_bilayer_overview(sim, patch_pack, roi_specs, fig_dir);
plot_model_maps(sim, T_pred, roi_specs, fig_dir);
plot_model_error_maps(sim, T_pred, roi_specs, fig_dir);
plot_roi_bars(T_roi, fig_dir);
plot_roi_distributions(T_pred, fig_dir);
plot_true_vs_pred_by_model(T_pred, fig_dir);
plot_roi_frequency_style_by_model(T_pred, cfg, fig_dir);

fprintf('\nBilayer transfer analysis complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
fprintf('\nROI metrics.\n');
disp(T_roi(:, {'model_name', 'roi_name', 'n', 'MAPE_pct', ...
    'RMSE_pct', 'bias_pct', 'CoV_pct', 'cs_true_median', ...
    'cs_pred_median'}));

%% Local functions

function MODEL = get_named_model(MODELS, name)

names = string({MODELS.name});
idx = find(names == string(name), 1);
if isempty(idx)
    error('Model not found in MODELS: %s', name);
end
MODEL = MODELS(idx).model;

end

function T = compute_improvement_summary(metrics_file)

Tmetrics = readtable(metrics_file, 'TextType', 'string');
Tbag = Tmetrics(Tmetrics.model_type == "bagged_trees" & ...
    Tmetrics.model_role == "operational", :);

local = Tbag(Tbag.model_name == "LocalOnly", :);
global_row = Tbag(Tbag.model_name == "GlobalOnly", :);
hybrid = Tbag(Tbag.model_name == "HybridLocalGlobal", :);

T = table();
T.comparison = [
    "Hybrid vs Local"
    "Hybrid vs Global"];
T.MAPE_from_pct = [
    local.MAPE_pct
    global_row.MAPE_pct];
T.MAPE_to_pct = [
    hybrid.MAPE_pct
    hybrid.MAPE_pct];
T.MAPE_absolute_improvement_pct_points = ...
    T.MAPE_from_pct - T.MAPE_to_pct;
T.MAPE_relative_improvement_pct = ...
    100 * T.MAPE_absolute_improvement_pct_points ./ T.MAPE_from_pct;
T.HighError_from_pct = [
    local.HighError_gt20_pct
    global_row.HighError_gt20_pct];
T.HighError_to_pct = [
    hybrid.HighError_gt20_pct
    hybrid.HighError_gt20_pct];
T.HighError_absolute_improvement_pct_points = ...
    T.HighError_from_pct - T.HighError_to_pct;
T.HighError_relative_reduction_pct = ...
    100 * T.HighError_absolute_improvement_pct_points ./ ...
    T.HighError_from_pct;

end

function roi_specs = build_roi_specs()

roi_specs = struct([]);

roi_specs(1).name = "soft_roi";
roi_specs(1).label = "Soft ROI";
roi_specs(1).xlim_m = [0.0100 0.0170];
roi_specs(1).zlim_m = [0.0215 0.0285];
roi_specs(1).expected_layer = "soft";

roi_specs(2).name = "hard_roi";
roi_specs(2).label = "Hard ROI";
roi_specs(2).xlim_m = [0.0330 0.0400];
roi_specs(2).zlim_m = [0.0215 0.0285];
roi_specs(2).expected_layer = "hard";

end

function patch_pack = build_dense_valid_patch_pack(cfg, feat_cfg, varargin)

p = inputParser;
addParameter(p, 'StepX', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'StepZ', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parse(p, varargin{:});

step_x = round(p.Results.StepX);
step_z = round(p.Results.StepZ);

Nx = round(cfg.Lx / cfg.dx) + 1;
Nz = round(cfg.Lz / cfg.dz) + 1;

win_size = round(feat_cfg.win_size);
if mod(win_size, 2) == 0
    win_size = win_size + 1;
end

half_win = floor(win_size / 2);

cx_list = (1 + half_win):step_x:(Nx - half_win);
cz_list = (1 + half_win):step_z:(Nz - half_win);

if isempty(cx_list) || isempty(cz_list)
    error('Dense valid patch grid is empty. Reduce M or stride.');
end

[CX, CZ] = meshgrid(cx_list, cz_list);
cx_all = reshape(CX.', 1, []);
cz_all = reshape(CZ.', 1, []);
n = numel(cx_all);

x_idx_list = cell(1, n);
z_idx_list = cell(1, n);
patch_labels = strings(1, n);

for i = 1:n
    x_idx_list{i} = (cx_all(i) - half_win):(cx_all(i) + half_win);
    z_idx_list{i} = (cz_all(i) - half_win):(cz_all(i) + half_win);
    patch_labels(i) = "valid_" + sprintf('%04d', i);
end

patch_pack = struct();
patch_pack.pattern = 'dense_valid';
patch_pack.edge_mode = 'valid';
patch_pack.step_x = step_x;
patch_pack.step_z = step_z;
patch_pack.n_patches = n;
patch_pack.cx_list = cx_all;
patch_pack.cz_list = cz_all;
patch_pack.x_idx_list = x_idx_list;
patch_pack.z_idx_list = z_idx_list;
patch_pack.patch_labels = patch_labels;
patch_pack.win_size = win_size;
patch_pack.half_win = half_win;
patch_pack.n_x = numel(cx_list);
patch_pack.n_z = numel(cz_list);
patch_pack.x_cent_m = (cx_list - 1) * cfg.dx;
patch_pack.z_cent_m = (cz_list - 1) * cfg.dz;

end

function [T, global_req] = extract_bilayer_feature_table( ...
    sim, cfg, feat_cfg, req_options, patch_pack, roi_specs)

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
        sim.Uxz, cfg, adaptive_req.config.default_req_config( ...
            cfg, feat_cfg, req_options{:}), feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features( ...
    global_curve);

global_req = struct();
global_req.q = q_global;
global_req.curve = global_curve;
global_req.mapping = adaptive_req.quantile.make_req_mapping(global_curve);
global_req.features = global_features;
global_req.shape_features = global_shape;

OUT = adaptive_req.estimators.req_estimator_map( ...
    sim.Uxz, cfg, feat_cfg, ...
    'StepX', patch_pack.step_x, ...
    'StepZ', patch_pack.step_z, ...
    'EdgeMode', 'valid', ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', true, ...
    'ReturnFeatureTable', true, ...
    'Verbose', false);

T = OUT.feature_table;
T.condition_id = ones(height(T), 1);
T.condition_position = ones(height(T), 1);
T.condition_label = repmat("bilayer_transfer", height(T), 1);
T.step_idx = ones(height(T), 1);
T.realization_idx = ones(height(T), 1);
T.patch_label = "valid_" + compose('%04d', T.patch_idx);
T.roi_name = strings(height(T), 1);
T.cs_true_patch_mean = nan(height(T), 1);
T.cs_true_patch_median = nan(height(T), 1);
T.cs_true_patch_std = nan(height(T), 1);
T.cs_true = nan(height(T), 1);
T.SIM_Nwaves = cfg.Nwaves * ones(height(T), 1);
T.SIM_SNR = cfg.SNR * ones(height(T), 1);
T.global_REQ_Nbins_effective = global_curve.Nbins_effective * ...
    ones(height(T), 1);
T.Omega_sr = cone_omega_sr(cfg.ConeHalfAngleDeg) * ones(height(T), 1);
T.omega_mean = T.Omega_sr;
T.q_global_theory = q_global * ones(height(T), 1);
T.global_req_mapping = repmat({global_req.mapping}, height(T), 1);

for i = 1:height(T)
    xi = patch_pack.x_idx_list{i};
    zi = patch_pack.z_idx_list{i};
    cs_patch = sim.cs_map(zi, xi);

    T.roi_name(i) = classify_roi_by_center( ...
        T.x_center_m(i), T.z_center_m(i), roi_specs);
    T.cs_true_patch_mean(i) = mean(cs_patch(:), 'omitnan');
    T.cs_true_patch_median(i) = median(cs_patch(:), 'omitnan');
    T.cs_true_patch_std(i) = std(cs_patch(:), 'omitnan');
    T.cs_true(i) = T.cs_true_patch_median(i);
end

T = assign_global_feature_columns(T, global_features, "global_");
T = assign_global_feature_columns(T, global_shape, "global_");

end

function rows = assign_scalar_fields_to_row(rows, row_idx, values, prefix)

if ~isstruct(values)
    return;
end

names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if isnumeric(value_i) && isscalar(value_i)
        rows(row_idx).(char(string(prefix) + string(names{i}))) = ...
            double(value_i);
    elseif islogical(value_i) && isscalar(value_i)
        rows(row_idx).(char(string(prefix) + string(names{i}))) = ...
            double(value_i);
    end
end

end

function roi = classify_roi_by_center(x, z, roi_specs)

roi = "outside_roi";
for i = 1:numel(roi_specs)
    in_x = x >= roi_specs(i).xlim_m(1) && x <= roi_specs(i).xlim_m(2);
    in_z = z >= roi_specs(i).zlim_m(1) && z <= roi_specs(i).zlim_m(2);
    if in_x && in_z
        roi = string(roi_specs(i).name);
        return;
    end
end

end

function T_pred = predict_global_model_single_q(MODEL_global, T_feat, cfg)

T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_global, T_feat, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "GlobalQSingleModel");

q_single = median(T_q.q_pred, 'omitnan');
q_map = q_single * ones(height(T_feat), 1);
cs_pred = q_to_cs_for_table(q_map, T_feat.req_mapping, cfg.f0);

T_pred = make_prediction_table_base(T_feat);
T_pred.model_name = repmat("GlobalQSingleModel", height(T_feat), 1);
T_pred.model_type = repmat("global_model_single_q", height(T_feat), 1);
T_pred.model_role = repmat("operational_global_q", height(T_feat), 1);
T_pred.q_pred = q_map;
T_pred.cs_pred = cs_pred;
T_pred.cs_global_constant = nan(height(T_feat), 1);

end

function T_pred = predict_fixed_q_to_sws(T_feat, model_name, q_fixed, f0)

q_map = q_fixed * ones(height(T_feat), 1);

T_pred = make_prediction_table_base(T_feat);
T_pred.model_name = repmat(string(model_name), height(T_feat), 1);
T_pred.model_type = repmat("fixed_theory_q", height(T_feat), 1);
T_pred.model_role = repmat("theory_baseline", height(T_feat), 1);
T_pred.q_pred = q_map;
T_pred.cs_pred = q_to_cs_for_table(q_map, T_feat.req_mapping, f0);
T_pred.cs_global_constant = nan(height(T_feat), 1);

end

function T_pred = predict_model_to_sws(MODEL, T_feat, model_name, f0)

T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL, T_feat, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', model_name);

T_pred = make_prediction_table_base(T_feat);
T_pred.model_name = repmat(string(model_name), height(T_feat), 1);
T_pred.model_type = repmat("bagged_trees", height(T_feat), 1);
T_pred.model_role = repmat("operational_transfer", height(T_feat), 1);
T_pred.q_pred = T_q.q_pred;
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, f0);
T_pred.cs_global_constant = nan(height(T_feat), 1);

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

function T = make_prediction_table_base(T_feat)

vars = {'patch_idx', 'patch_label', 'cx', 'cz', ...
    'x_center_m', 'z_center_m', 'roi_name', ...
    'cs_true_patch_mean', 'cs_true_patch_median', ...
    'cs_true_patch_std', 'cs_true', 'SIM_f0', 'SIM_cs_bg', ...
    'SIM_WaveModel', 'REQ_M', 'Omega_sr'};
T = T_feat(:, vars);

end

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)

n = numel(q_pred);
cs_pred = nan(n, 1);

for i = 1:n
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
[G, T_roi] = findgroups(T(:, {'model_name', 'model_type', ...
    'model_role', 'roi_name'}));
T_roi.n = splitapply(@numel, T.cs_error_pct, G);
T_roi.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), ...
    T.cs_error_pct, G);
T_roi.bias_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.CoV_pct = splitapply(@cov_pct, T.cs_pred, G);
T_roi.cs_true_median = splitapply(@(x) median(x, 'omitnan'), ...
    T.cs_true, G);
T_roi.cs_pred_median = splitapply(@(x) median(x, 'omitnan'), ...
    T.cs_pred, G);
T_roi.cs_pred_mean = splitapply(@(x) mean(x, 'omitnan'), ...
    T.cs_pred, G);
T_roi.cs_pred_std = splitapply(@(x) std(x, 'omitnan'), ...
    T.cs_pred, G);

end

function T = add_all_roi(T)

T.roi_name = repmat("all", height(T), 1);

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function plot_bilayer_overview(sim, patch_pack, roi_specs, fig_dir)

fig = figure('Color', 'w', 'Position', [100 100 950 760]);
imagesc(sim.x * 100, sim.z * 100, sim.cs_map);
axis image;
set(gca, 'YDir', 'reverse');
colormap parula;
cb = colorbar;
cb.Label.String = 'c_s (m/s)';
title('Bilayer true c_s map, valid centers, ROIs', ...
    'Interpreter', 'tex');
xlabel('x (cm)');
ylabel('z (cm)');
hold on;
scatter(sim.x(patch_pack.cx_list) * 100, ...
    sim.z(patch_pack.cz_list) * 100, 6, 'w', 'filled', ...
    'MarkerEdgeColor', 'k');
draw_rois(roi_specs);
save_fig(fig, fig_dir, 'level11_level02_bilayer_true_map_rois');

end

function plot_model_maps(sim, T_pred, roi_specs, fig_dir)

models = unique(T_pred.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [80 80 1550 520]);
tl = tiledlayout(fig, 1, numel(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T_pred(T_pred.model_name == models(i), :);
    [map_i, xc, zc] = table_to_map(Ti, 'cs_pred');
    imagesc(ax, xc * 100, zc * 100, map_i);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    xlim(ax, [min(sim.x) max(sim.x)] * 100);
    ylim(ax, [min(sim.z) max(sim.z)] * 100);
    colormap(ax, parula);
    clim(ax, [min(sim.cs_map(:)) max(sim.cs_map(:))]);
    cb = colorbar(ax);
    cb.Label.String = 'c_s (m/s)';
    title(ax, models(i), 'Interpreter', 'none');
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    hold(ax, 'on');
    draw_rois(roi_specs, ax);
end

title(tl, 'Bilayer predicted c_s maps', 'Interpreter', 'tex');
save_fig(fig, fig_dir, 'level11_level02_bilayer_predicted_sws_maps');

end

function plot_model_error_maps(sim, T_pred, roi_specs, fig_dir)

models = unique(T_pred.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [80 80 1550 520]);
tl = tiledlayout(fig, 1, numel(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T_pred(T_pred.model_name == models(i), :);
    [map_i, xc, zc] = table_to_map(Ti, 'cs_error_pct');
    imagesc(ax, xc * 100, zc * 100, map_i);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    xlim(ax, [min(sim.x) max(sim.x)] * 100);
    ylim(ax, [min(sim.z) max(sim.z)] * 100);
    colormap(ax, parula);
    clim(ax, [-50 50]);
    cb = colorbar(ax);
    cb.Label.String = 'c_s error (%)';
    title(ax, models(i), 'Interpreter', 'none');
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    hold(ax, 'on');
    draw_rois(roi_specs, ax);
end

title(tl, 'Bilayer signed c_s error maps (%)', 'Interpreter', 'tex');
save_fig(fig, fig_dir, 'level11_level02_bilayer_error_maps');

end

function plot_roi_bars(T_roi, fig_dir)

fig = figure('Color', 'w', 'Position', [100 100 1350 540]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

labels = categorical(T_roi.model_name + " / " + T_roi.roi_name);

ax = nexttile(tl);
bar(ax, labels, T_roi.MAPE_pct);
ylabel(ax, 'MAPE (%)');
title(ax, 'ROI accuracy');
grid(ax, 'on');
xtickangle(ax, 30);

ax = nexttile(tl);
bar(ax, labels, T_roi.CoV_pct);
ylabel(ax, 'CoV (%)');
title(ax, 'ROI precision');
grid(ax, 'on');
xtickangle(ax, 30);

save_fig(fig, fig_dir, 'level11_level02_bilayer_roi_mape_cov');

end

function [map, xc, zc] = table_to_map(T, value_var)

xc = unique(T.x_center_m, 'stable');
zc = unique(T.z_center_m, 'stable');

xc = sort(xc(:).');
zc = sort(zc(:));

map = nan(numel(zc), numel(xc));

[~, ix] = ismember(T.x_center_m, xc);
[~, iz] = ismember(T.z_center_m, zc);

for i = 1:height(T)
    if ix(i) > 0 && iz(i) > 0
        map(iz(i), ix(i)) = T.(value_var)(i);
    end
end

end

function plot_roi_distributions(T_pred, fig_dir)

T = T_pred(T_pred.roi_name ~= "outside_roi", :);
fig = figure('Color', 'w', 'Position', [100 100 1250 540]);
boxchart(categorical(T.model_name + " / " + T.roi_name), ...
    T.cs_abs_error_pct);
ylabel('|c_s error| (%)', 'Interpreter', 'tex');
title('Bilayer ROI patch error distributions', 'Interpreter', 'tex');
grid on;
xtickangle(30);
save_fig(fig, fig_dir, ...
    'level11_level02_bilayer_roi_error_distributions');

end

function plot_true_vs_pred_by_model(T_pred, fig_dir)

fig = figure('Color', 'w', 'Position', [100 100 1250 460]);
models = unique(T_pred.model_name, 'stable');
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T_pred(T_pred.model_name == models(i), :);
    scatter(ax, Ti.cs_true, Ti.cs_pred, 22, ...
        categorical(Ti.roi_name), 'filled');
    hold(ax, 'on');
    lo = min([Ti.cs_true; Ti.cs_pred], [], 'omitnan');
    hi = max([Ti.cs_true; Ti.cs_pred], [], 'omitnan');
    plot(ax, [lo hi], [lo hi], 'k--');
    axis(ax, 'equal');
    xlabel(ax, 'True c_s (m/s)', 'Interpreter', 'tex');
    ylabel(ax, 'Predicted c_s (m/s)', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    grid(ax, 'on');
end

save_fig(fig, fig_dir, 'level11_level02_bilayer_true_vs_pred');

end

function plot_roi_frequency_style_by_model(T_pred, cfg, fig_dir)

T = T_pred(T_pred.roi_name ~= "outside_roi", :);
models = unique(T.model_name, 'stable');
f0 = cfg.f0;

for i = 1:numel(models)
    Ti = T(T.model_name == models(i), :);
    soft = Ti(Ti.roi_name == "soft_roi", :);
    hard = Ti(Ti.roi_name == "hard_roi", :);

    fig = figure('Color', 'w', 'Position', [100 100 760 620]);
    ax = axes(fig);
    hold(ax, 'on');

    h_soft = errorbar(ax, f0, mean(soft.cs_pred, 'omitnan'), ...
        std(soft.cs_pred, 0, 'omitnan'), 'o-', ...
        'LineWidth', 2.2, ...
        'MarkerSize', 9, ...
        'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
        'Color', [0.8500 0.3250 0.0980], ...
        'DisplayName', 'Soft');

    h_hard = errorbar(ax, f0, mean(hard.cs_pred, 'omitnan'), ...
        std(hard.cs_pred, 0, 'omitnan'), 'o-', ...
        'LineWidth', 2.2, ...
        'MarkerSize', 9, ...
        'MarkerFaceColor', [0 0.4470 0.7410], ...
        'Color', [0 0.4470 0.7410], ...
        'DisplayName', 'Hard');

    yline(ax, cfg.cs_bg, '--', 'Color', [0.8500 0.3250 0.0980], ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
    yline(ax, cfg.cs_inc, '--', 'Color', [0 0.4470 0.7410], ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'c_s (m/s)', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    legend(ax, [h_soft h_hard], {'Soft', 'Hard'}, 'Location', 'best');
    xlim(ax, [f0 - 75, f0 + 75]);
    ylim(ax, [min(cfg.cs_bg, cfg.cs_inc) - 0.4, ...
        max(cfg.cs_bg, cfg.cs_inc) + 0.4]);

    adaptive_req.templates.apply_paper_style(fig, ax, ...
        'Times New Roman', 22);
    save_fig(fig, fig_dir, "level11_level02_frequency_style_" + ...
        sanitize_filename(models(i)));
end

end

function draw_rois(roi_specs, ax)

if nargin < 2
    ax = gca;
end

for i = 1:numel(roi_specs)
    x_cm = roi_specs(i).xlim_m * 100;
    z_cm = roi_specs(i).zlim_m * 100;
    rectangle(ax, 'Position', [x_cm(1), z_cm(1), ...
        diff(x_cm), diff(z_cm)], ...
        'EdgeColor', 'w', 'LineWidth', 2.2, 'LineStyle', '-');
    text(ax, x_cm(1), z_cm(2) + 0.12, roi_specs(i).label, ...
        'Color', 'w', 'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

end

function omega = cone_omega_sr(half_angle_deg)

theta = deg2rad(half_angle_deg);
omega = 2*pi*(1 - cos(theta));

end

function name = sanitize_filename(x)

name = regexprep(string(x), '[^A-Za-z0-9_]+', '_');
name = matlab.lang.makeValidName(char(name));
name = string(name);

end

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end

function save_fig(fig, fig_dir, base_name)

adaptive_req.templates.apply_paper_style(fig, [], 'Times New Roman', 18);
adaptive_req.templates.paper_export(fig, fig_dir, base_name, ...
    'DPI', 300, ...
    'VectorPDF', true, ...
    'WidthCm', 28, ...
    'HeightCm', 11);
close(fig);

end
