%% analyze_test_11_level03_bilayer_frequency_sweep.m
% Test 11 Level 03: bilayer transfer over frequency.
%
% This analysis uses a dedicated config profile:
%   configs/test_11_level03_bilayer_frequency_sweep.m
%
% It compares:
%   1. TheoryDiffuse3D: fixed theoretical discrete q.
%   2. GlobalQSingleModel: one global ML q applied to all windows.
%   3. LocalOnly: local patch features only.
%   4. HybridLocalGlobal: local patch features plus global context.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = adaptive_req.config.load_profile_config( ...
    'test_11_level03_bilayer_frequency_sweep', ...
    'RootDir', root_dir);

%% Locate Test 11 Level 01 models

[~, ~, PATHS11] = adaptive_req.analysis.load_mc_results( ...
    'test_11_global_req_features', ...
    'RootDir', root_dir, ...
    'Verbose', true);

level01_dir = fullfile(PATHS11.analysis_dir, ...
    'level_01_global_vs_local');
model_file = fullfile(level01_dir, 'models', ...
    'level11_global_vs_local_models.mat');

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
    'level_03_bilayer_frequency_sweep');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
data_dir = fullfile(analysis_dir, 'data');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(data_dir);

%% Run frequency sweep

roi_specs = CFG.BILAYER.ROIS;
model_order = CFG.BILAYER.MODELS;

T_all_feat = table();
T_all_pred = table();
T_all_roi = table();
T_all_overall = table();
sim_summary = struct([]);

for fi = 1:numel(CFG.BILAYER.f_list)
    f0 = CFG.BILAYER.f_list(fi);
    fprintf('\n=== Bilayer frequency sweep: f0 = %d Hz ===\n', f0);

    cfg = build_bilayer_sim_cfg(CFG, f0, fi);
    [feat_cfg, req_options] = build_req_settings(CFG);
    [req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
        cfg, feat_cfg, req_options{:});

    fprintf('Running bilayer simulation...\n');
    sim = adaptive_req.simulate.run_single_simulation(cfg);

    fprintf('Extracting local and global REQ features...\n');
    [T_feat, global_req, req_out] = extract_frequency_feature_table( ...
        sim, cfg, feat_cfg, req_cfg, req_options, roi_specs, CFG, fi);

    fprintf('Applying theory and ML q models...\n');
    q_theory_discrete = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, feat_cfg.cs_guess_used, ...
        'M', feat_cfg.M, ...
        'Gamma', feat_cfg.gamma_win, ...
        'PadFactor', feat_cfg.pad_factor, ...
        'Nbins', req_cfg.Nbins, ...
        'SmoothSigma', req_cfg.smooth_sigma, ...
        'TheoryMode', CFG.BILAYER.REQ.TheoryMode, ...
        'FieldType', CFG.BILAYER.REQ.TheoryFieldType, ...
        'Plot', false);

    T_pred_theory = predict_fixed_q_to_sws( ...
        T_feat, "TheoryDiffuse3D", q_theory_discrete.q_th, cfg.f0);
    T_pred_global = predict_global_model_single_q( ...
        MODEL_global, T_feat, cfg.f0);
    T_pred_local = predict_model_to_sws( ...
        MODEL_local, T_feat, "LocalOnly", cfg.f0);
    T_pred_hybrid = predict_model_to_sws( ...
        MODEL_hybrid, T_feat, "HybridLocalGlobal", cfg.f0);

    T_pred_f = [T_pred_theory; T_pred_global; ...
        T_pred_local; T_pred_hybrid];
    T_pred_f = add_error_metrics(T_pred_f);
    T_pred_f = order_prediction_models(T_pred_f, model_order);

    T_roi_f = summarize_roi_metrics(T_pred_f);
    T_overall_f = summarize_roi_metrics(add_all_roi(T_pred_f));

    T_all_feat = [T_all_feat; remove_cell_columns(T_feat)]; %#ok<AGROW>
    T_all_pred = [T_all_pred; T_pred_f]; %#ok<AGROW>
    T_all_roi = [T_all_roi; T_roi_f]; %#ok<AGROW>
    T_all_overall = [T_all_overall; T_overall_f]; %#ok<AGROW>

    sim_summary(fi).f0 = f0;
    sim_summary(fi).cfg = cfg;
    sim_summary(fi).global_req = global_req;
    sim_summary(fi).q_theory_discrete = q_theory_discrete;
    sim_summary(fi).x_centers_m = req_out.x_m;
    sim_summary(fi).z_centers_m = req_out.z_m;

    save(fullfile(data_dir, sprintf( ...
        'level11_level03_bilayer_frequency_%03dHz.mat', f0)), ...
        'CFG', 'cfg', 'feat_cfg', 'req_cfg', 'sim', 'req_out', ...
        'global_req', 'q_theory_discrete', 'T_feat', 'T_pred_f', ...
        'T_roi_f', 'T_overall_f', '-v7.3');

    plot_frequency_maps(sim, T_pred_f, roi_specs, model_order, CFG, ...
        fig_dir, f0, req_out);
end

%% Save combined tables

writetable(T_all_feat, fullfile(table_dir, ...
    'level11_level03_bilayer_frequency_features.csv'));
writetable(T_all_pred, fullfile(table_dir, ...
    'level11_level03_bilayer_frequency_predictions.csv'));
writetable(T_all_roi, fullfile(table_dir, ...
    'level11_level03_bilayer_frequency_roi_metrics.csv'));
writetable(T_all_overall, fullfile(table_dir, ...
    'level11_level03_bilayer_frequency_overall_metrics.csv'));

save(fullfile(data_dir, 'level11_level03_bilayer_frequency_sweep.mat'), ...
    'CFG', 'T_all_feat', 'T_all_pred', 'T_all_roi', 'T_all_overall', ...
    'sim_summary', '-v7.3');

%% Cross-frequency figures

plot_cross_frequency_figures(T_all_pred, T_all_roi, roi_specs, ...
    model_order, CFG, fig_dir);

fprintf('\nLevel 03 bilayer frequency sweep complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
fprintf('\nROI metrics by frequency.\n');
disp(T_all_roi(:, {'frequency_hz', 'model_name', 'roi_name', ...
    'n', 'MAPE_pct', 'CoV_pct', 'cs_true_median', ...
    'cs_pred_mean'}));

%% Local functions

function MODEL = get_named_model(MODELS, name)

names = string({MODELS.name});
idx = find(names == string(name), 1);
if isempty(idx)
    error('Model not found in MODELS: %s', name);
end
MODEL = MODELS(idx).model;

end

function cfg = build_bilayer_sim_cfg(CFG, f0, freq_idx)

S = CFG.BILAYER.SIM;
cfg = adaptive_req.config.default_sim_config( ...
    'WaveModel', S.WaveModel, ...
    'f0', f0, ...
    'cs_bg', S.cs_bg, ...
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
    'ConeHalfAngleDeg', S.ConeHalfAngleDeg, ...
    'ForceInPlaneWave', S.ForceInPlaneWave, ...
    'UseParfor', S.UseParfor, ...
    'Seed', S.SeedBase + freq_idx - 1);

cfg.MaskType = CFG.BILAYER.MaskType;
cfg.cs_inc = S.cs_inc;
cfg.MaskParams = CFG.BILAYER.MaskParams;

end

function [feat_cfg, req_options] = build_req_settings(CFG)

R = CFG.BILAYER.REQ;
feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', R.M, ...
    'cs_guess', R.cs_guess, ...
    'gamma_win', R.gamma_win, ...
    'pad_factor', R.pad_factor);

req_options = { ...
    'Nbins', R.Nbins, ...
    'Nbins_auto_oversample', R.Nbins_auto_oversample, ...
    'Nbins_min', R.Nbins_min, ...
    'smooth_sigma', R.smooth_sigma};

end

function [T, global_req, OUT] = extract_frequency_feature_table( ...
    sim, cfg, feat_cfg, req_cfg, req_options, roi_specs, CFG, frequency_idx)

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
        sim.Uxz, cfg, req_cfg, feat_cfg);
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
    'StepX', CFG.BILAYER.REQ.StepX, ...
    'StepZ', CFG.BILAYER.REQ.StepZ, ...
    'EdgeMode', CFG.BILAYER.REQ.EdgeMode, ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', true, ...
    'ReturnFeatureTable', true, ...
    'Verbose', false);

T = OUT.feature_table;
T.condition_id = frequency_idx * ones(height(T), 1);
T.condition_position = frequency_idx * ones(height(T), 1);
T.condition_label = repmat(string(sprintf( ...
    'bilayer_frequency_%03dHz', cfg.f0)), height(T), 1);
T.step_idx = frequency_idx * ones(height(T), 1);
T.realization_idx = ones(height(T), 1);
T.frequency_hz = cfg.f0 * ones(height(T), 1);
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
T.q_global_req = q_global * ones(height(T), 1);
T.global_req_mapping = repmat({global_req.mapping}, height(T), 1);

for i = 1:height(T)
    xi = (T.cx(i) - OUT.half_win):(T.cx(i) + OUT.half_win);
    zi = (T.cz(i) - OUT.half_win):(T.cz(i) + OUT.half_win);
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

function T_pred = predict_global_model_single_q(MODEL_global, T_feat, f0)

T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_global, T_feat, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "GlobalQSingleModel");

q_single = median(T_q.q_pred, 'omitnan');
q_map = q_single * ones(height(T_feat), 1);
cs_pred = q_to_cs_for_table(q_map, T_feat.req_mapping, f0);

T_pred = make_prediction_table_base(T_feat);
T_pred.model_name = repmat("GlobalQSingleModel", height(T_feat), 1);
T_pred.model_type = repmat("global_model_single_q", height(T_feat), 1);
T_pred.model_role = repmat("operational_global_q", height(T_feat), 1);
T_pred.q_pred = q_map;
T_pred.cs_pred = cs_pred;

end

function T_pred = predict_fixed_q_to_sws(T_feat, model_name, q_fixed, f0)

q_map = q_fixed * ones(height(T_feat), 1);

T_pred = make_prediction_table_base(T_feat);
T_pred.model_name = repmat(string(model_name), height(T_feat), 1);
T_pred.model_type = repmat("fixed_theory_q", height(T_feat), 1);
T_pred.model_role = repmat("theory_baseline", height(T_feat), 1);
T_pred.q_pred = q_map;
T_pred.cs_pred = q_to_cs_for_table(q_map, T_feat.req_mapping, f0);

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

end

function T = make_prediction_table_base(T_feat)

vars = {'condition_id', 'condition_position', 'condition_label', ...
    'step_idx', 'realization_idx', 'frequency_hz', ...
    'patch_idx', 'patch_label', 'cx', 'cz', ...
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
[G, T_roi] = findgroups(T(:, {'frequency_hz', 'model_name', ...
    'model_type', 'model_role', 'roi_name'}));
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

function T = order_prediction_models(T, model_order)

model_order = string(model_order);
T.model_name = categorical(string(T.model_name), model_order, ...
    'Ordinal', true);
T = sortrows(T, {'model_name', 'patch_idx'});
T.model_name = string(T.model_name);

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

function plot_frequency_maps(sim, T_pred_f, roi_specs, model_order, CFG, ...
    fig_dir, f0, req_out)

F = CFG.BILAYER.FIG;
cs_limits = [CFG.BILAYER.SIM.cs_bg, CFG.BILAYER.SIM.cs_inc];

fig = adaptive_req.figures.plot_bilayer_true_map( ...
    sim, roi_specs, ...
    'XCentersM', req_out.x_m, ...
    'ZCentersM', req_out.z_m, ...
    'Title', sprintf('Bilayer true c_s map, %d Hz', f0), ...
    'ColorLimits', cs_limits, ...
    'FigureSize', [16 13], ...
    'FontSize', F.FontSize, ...
    'ShowPatchCenters', F.ShowPatchCenters);
save_fig(fig, fig_dir, sprintf( ...
    'level11_level03_bilayer_true_map_rois_%03dHz', f0), ...
    16, 13);

fig = adaptive_req.figures.plot_bilayer_model_maps( ...
    T_pred_f, roi_specs, ...
    'ValueVar', 'cs_pred', ...
    'Models', model_order, ...
    'Title', sprintf('Bilayer predicted c_s maps, %d Hz', f0), ...
    'ColorLimits', cs_limits, ...
    'ColorbarLabel', 'c_s (m/s)', ...
    'FigureSize', [F.MapWidthCm, F.MapHeightCm], ...
    'FontSize', F.FontSize, ...
    'ShowRoiLabels', false);
save_fig(fig, fig_dir, sprintf( ...
    'level11_level03_bilayer_predicted_sws_maps_%03dHz', f0), ...
    F.MapWidthCm, F.MapHeightCm);

fig = adaptive_req.figures.plot_bilayer_model_maps( ...
    T_pred_f, roi_specs, ...
    'ValueVar', 'cs_error_pct', ...
    'Models', model_order, ...
    'Title', sprintf('Bilayer signed c_s error maps, %d Hz', f0), ...
    'ColorLimits', [-30 30], ...
    'ColorbarLabel', 'c_s error (%)', ...
    'FigureSize', [F.MapWidthCm, F.MapHeightCm], ...
    'FontSize', F.FontSize, ...
    'ShowRoiLabels', false);
save_fig(fig, fig_dir, sprintf( ...
    'level11_level03_bilayer_error_maps_%03dHz', f0), ...
    F.MapWidthCm, F.MapHeightCm);

end

function plot_cross_frequency_figures(T_all_pred, T_all_roi, roi_specs, ...
    model_order, CFG, fig_dir)

F = CFG.BILAYER.FIG;
fig = adaptive_req.figures.plot_bilayer_roi_error_distribution( ...
    T_all_pred, ...
    'Title', 'Bilayer ROI patch error distributions, all frequencies', ...
    'FigureSize', [32 11], ...
    'FontSize', F.FontSize);
save_fig(fig, fig_dir, ...
    'level11_level03_bilayer_roi_error_distributions_all_freq', ...
    32, 11);

fig = adaptive_req.figures.plot_bilayer_roi_bars( ...
    T_all_roi, ...
    'Title', 'Bilayer ROI accuracy and precision, all frequencies', ...
    'FigureSize', [32 12], ...
    'FontSize', F.FontSize);
save_fig(fig, fig_dir, ...
    'level11_level03_bilayer_roi_mape_cov_all_freq', ...
    32, 12);

for i = 1:numel(model_order)
    model_name = string(model_order(i));
    fig = adaptive_req.figures.plot_bilayer_frequency_curves( ...
        T_all_roi, model_name, roi_specs, ...
        'Title', readable_model_name(model_name), ...
        'FigureSize', [F.FrequencyWidthCm, F.FrequencyHeightCm], ...
        'FontSize', F.FontSize);
    save_fig(fig, fig_dir, ...
        "level11_level03_frequency_curve_" + sanitize_filename(model_name), ...
        F.FrequencyWidthCm, F.FrequencyHeightCm);
end

end

function name = readable_model_name(name)

name = string(name);
name = replace(name, "TheoryDiffuse3D", "Theory diffuse 3D");
name = replace(name, "GlobalQSingleModel", "Global q model");
name = replace(name, "LocalOnly", "Local only");
name = replace(name, "HybridLocalGlobal", "Hybrid");

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

function save_fig(fig, fig_dir, base_name, width_cm, height_cm)

adaptive_req.templates.paper_export(fig, fig_dir, base_name, ...
    'DPI', 300, ...
    'VectorPDF', true, ...
    'WidthCm', width_cm, ...
    'HeightCm', height_cm);
close(fig);

end
