%% analyze_test_14_level01_resolution_sensitivity.m
% Test 14 Level 01: apply the trained clean model to a dx=dz sweep.
%
% This analysis does not train or retrain any model. It applies the Test 12
% deployment model to the Test 14 resolution-sensitivity dataset and reports
% SWS accuracy/precision as a function of spatial discretization.

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

%% Locate data and output folders

PATHS14 = locate_latest_run(root_dir, 'test_14_dx_dz_resolution_sweep');
PATHS12 = locate_latest_run(root_dir, 'test_12_cs_guess_window_sweep');

fprintf('\nUsing Test 14 run folder:\n%s\n', PATHS14.run_dir);
fprintf('Using Test 12 model source:\n%s\n', PATHS12.run_dir);

analysis_dir = fullfile(PATHS14.run_dir, 'analysis', ...
    'level_01_resolution_sensitivity');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
model_dir = fullfile(analysis_dir, 'models');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(model_dir);

S = load(fullfile(PATHS14.run_dir, 'data', ...
    'test14_dx_dz_resolution_sweep_results.mat'), ...
    'T_all', 'MC', 'CFG');
T_feat = S.T_all;
MC = S.MC;
CFG = S.CFG;

assert(~isempty(T_feat), 'Test 14 feature table is empty.');

%% Load deployment model

deployment_dir = fullfile(PATHS12.run_dir, 'analysis', ...
    'level_01_model_comparison', 'models', 'deployment');

[MODEL_DEPLOY, MODEL_INFO, model_file] = load_resolution_model(deployment_dir);

fprintf('\nLoaded deployment model:\n%s\n', model_file);
disp(struct2table(MODEL_INFO));

%% Predict q and convert to SWS

T_feat = ensure_model_predictors(T_feat, MODEL_DEPLOY);

fprintf('\nApplying deployment model to Test 14 features...\n');
T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_DEPLOY, T_feat, ...
    'ModelType', MODEL_INFO.model_type, ...
    'ModelName', MODEL_INFO.model_name);

T_pred = make_prediction_base(T_feat);
T_pred.model_name = repmat(string(MODEL_INFO.model_name), height(T_feat), 1);
T_pred.feature_set = repmat(string(MODEL_INFO.feature_set), height(T_feat), 1);
T_pred.model_type = repmat(string(MODEL_INFO.model_type), height(T_feat), 1);
T_pred.model_role = repmat("test14_resolution_transfer", height(T_feat), 1);
T_pred.q_pred_raw = T_q.q_pred_raw;
T_pred.q_pred = min(max(T_q.q_pred, 0.001), 0.999);
T_pred.q_true = T_feat.q_theory;
T_pred.q_residual = T_pred.q_pred - T_pred.q_true;
T_pred.q_abs_error = abs(T_pred.q_residual);
T_pred.cs_true = T_feat.cs_true;
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, ...
    T_feat.req_mapping, T_feat.SIM_f0);
T_pred = add_error_metrics(T_pred);

assert(mean(isfinite(T_pred.cs_pred)) > 0.95, ...
    'More than 5%% of cs_pred values are non-finite.');

%% Metrics

T_condition_metrics = summarize_condition_metrics(T_pred);
T_dx_metrics = summarize_group_metrics(T_pred, {'SIM_dx'});
T_dx_M_metrics = summarize_group_metrics(T_pred, {'SIM_dx', 'REQ_M'});
T_ppw_metrics = summarize_group_metrics(T_pred, {'pixels_per_wavelength'});

writetable(remove_cell_columns(T_pred), fullfile(table_dir, ...
    'level14_level01_resolution_predictions.csv'));
writetable(T_condition_metrics, fullfile(table_dir, ...
    'level14_level01_resolution_metrics_by_condition.csv'));
writetable(T_dx_metrics, fullfile(table_dir, ...
    'level14_level01_resolution_metrics_by_dx.csv'));
writetable(T_dx_M_metrics, fullfile(table_dir, ...
    'level14_level01_resolution_metrics_by_M_dx.csv'));
writetable(T_ppw_metrics, fullfile(table_dir, ...
    'level14_level01_resolution_metrics_by_pixels_per_wavelength.csv'));

save(fullfile(model_dir, ...
    'level14_level01_resolution_sensitivity_results.mat'), ...
    'MODEL_DEPLOY', 'MODEL_INFO', 'model_file', 'T_pred', ...
    'T_condition_metrics', 'T_dx_metrics', 'T_dx_M_metrics', ...
    'T_ppw_metrics', 'MC', 'CFG', '-v7.3');

%% Figures

plot_metric_vs_dx(T_condition_metrics, 'MAPE_pct', fig_dir, ...
    'level14_level01_mape_vs_dx.png', 'MAPE vs dx');
plot_metric_vs_dx(T_condition_metrics, 'CoV_pct', fig_dir, ...
    'level14_level01_cov_vs_dx.png', 'CoV vs dx');
plot_mape_vs_pixels_per_wavelength(T_condition_metrics, fig_dir);
plot_mape_heatmap_M_dx(T_condition_metrics, fig_dir);

%% Console summary

fprintf('\nTest 14 Level 01 complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
fprintf('\nMetrics by dx:\n');
disp(T_dx_metrics(:, {'SIM_dx', 'N', 'MAPE_pct', 'bias_pct', ...
    'CoV_pct', 'HighError_gt10_pct', 'HighError_gt20_pct'}));

[~, best_idx] = min(T_dx_metrics.MAPE_pct);
[~, worst_idx] = max(T_dx_metrics.MAPE_pct);
fprintf('\nBest dx by MAPE: %.4g m (MAPE %.2f%%)\n', ...
    T_dx_metrics.SIM_dx(best_idx), T_dx_metrics.MAPE_pct(best_idx));
fprintf('Worst dx by MAPE: %.4g m (MAPE %.2f%%)\n', ...
    T_dx_metrics.SIM_dx(worst_idx), T_dx_metrics.MAPE_pct(worst_idx));

%% Local functions

function PATHS = locate_latest_run(root_dir, test_name)

out_root = fullfile(root_dir, 'outputs', test_name);
runs = dir(fullfile(out_root, test_name + "_*"));
runs = runs([runs.isdir]);
assert(~isempty(runs), 'No run folders found in %s', out_root);

[~, idx] = max([runs.datenum]);
PATHS = struct();
PATHS.run_dir = fullfile(runs(idx).folder, runs(idx).name);

end

function [MODEL_DEPLOY, MODEL_INFO, model_file] = load_resolution_model(deployment_dir)

try
    [MODEL_DEPLOY, MODEL_INFO, model_file] = ...
        adaptive_req.analysis.load_q_model_deployment( ...
        deployment_dir, ...
        'ModelName', 'HybridLocalGlobal', ...
        'FeatureSet', 'WithCsGuess', ...
        'ModelType', 'bagged_trees');
catch
    [MODEL_DEPLOY, MODEL_INFO, model_file] = ...
        adaptive_req.analysis.load_q_model_deployment( ...
        deployment_dir, ...
        'ModelName', 'HybridLocalGlobal', ...
        'FeatureSet', 'CsGuess', ...
        'ModelType', 'bagged_trees');
end

end

function T = ensure_model_predictors(T, MODEL)

required = string({MODEL.encoder.entries.name});
vars = string(T.Properties.VariableNames);
missing = setdiff(required, vars);

if isempty(missing)
    return;
end

can_fill_local_ecum = ismember('req_mapping', vars);

if can_fill_local_ecum
    local_missing = missing(startsWith(missing, ["ecum_", "srad_proxy_"]));
    if ~isempty(local_missing)
        fprintf('Reconstructing %d missing local Ecum/Srad-proxy predictors from req_mapping.\n', ...
            numel(local_missing));
        T = add_local_ecum_features_from_mapping(T, local_missing);
    end
end

vars = string(T.Properties.VariableNames);
missing = setdiff(required, vars);

if ~isempty(missing)
    error('Required deployment predictors are still missing: %s', ...
        strjoin(missing, ', '));
end

end

function T = add_local_ecum_features_from_mapping(T, feature_names)

feature_names = string(feature_names(:));
for fi = 1:numel(feature_names)
    T.(char(feature_names(fi))) = nan(height(T), 1);
end

for i = 1:height(T)
    if isempty(T.req_mapping{i})
        continue;
    end
    feat_i = adaptive_req.quantile.extract_ecum_shape_features(T.req_mapping{i});
    for fi = 1:numel(feature_names)
        name_i = char(feature_names(fi));
        if isfield(feat_i, name_i)
            T.(name_i)(i) = feat_i.(name_i);
        end
    end
end

end

function T = make_prediction_base(T_feat)

vars = [
    "condition_id"
    "condition_position"
    "condition_label"
    "resolution_idx"
    "dx_dz_value"
    "step_idx"
    "aperture_label"
    "realization_idx"
    "patch_idx"
    "SIM_WaveModel"
    "SIM_f0"
    "SIM_cs_bg"
    "SIM_dx"
    "SIM_dz"
    "REQ_M"
    "REQ_cs_guess"
    "REQ_win_size"
    "pixels_per_wavelength"
    "pixels_per_window"
    "k0_over_knyquist"
    "M_eff_guess"
    "M_eff_true_diag"
    "Omega_sr"
    "omega_mean"
    "cs_true"
    ];

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
        mappings{i}, q_pred(i), f0(i));
end

end

function T = add_error_metrics(T)

T.cs_error = T.cs_pred - T.cs_true;
T.cs_abs_error = abs(T.cs_error);
T.cs_error_pct = 100 * T.cs_error ./ T.cs_true;
T.cs_abs_error_pct = abs(T.cs_error_pct);

end

function T_metrics = summarize_condition_metrics(T)

group_vars = {'resolution_idx', 'SIM_dx', 'SIM_dz', 'SIM_f0', ...
    'SIM_cs_bg', 'REQ_M', 'REQ_cs_guess', 'step_idx', ...
    'aperture_label', 'pixels_per_wavelength', 'pixels_per_window', ...
    'k0_over_knyquist', 'M_eff_guess', 'M_eff_true_diag'};
group_vars = group_vars(ismember(group_vars, T.Properties.VariableNames));
T_metrics = summarize_group_metrics(T, group_vars);

end

function T_metrics = summarize_group_metrics(T, group_vars)

if isempty(T)
    T_metrics = table();
    return;
end

[G, T_metrics] = findgroups(T(:, group_vars));
T_metrics.N = splitapply(@numel, T.cs_error_pct, G);
T_metrics.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), ...
    T.cs_error_pct, G);
T_metrics.bias_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.cs_error_pct, G);
T_metrics.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), ...
    T.cs_error_pct, G);
T_metrics.CoV_pct = splitapply(@cov_pct, T.cs_pred, G);
T_metrics.HighError_gt10_pct = splitapply(@(x) 100 * mean(abs(x) > 10, 'omitnan'), ...
    T.cs_error_pct, G);
T_metrics.HighError_gt20_pct = splitapply(@(x) 100 * mean(abs(x) > 20, 'omitnan'), ...
    T.cs_error_pct, G);
T_metrics.P95_abs_error_pct = splitapply(@(x) prctile_finite(abs(x), 95), ...
    T.cs_error_pct, G);

if ismember('pixels_per_wavelength', T.Properties.VariableNames) && ...
        ~ismember('pixels_per_wavelength', string(group_vars))
    T_metrics.pixels_per_wavelength_mean = splitapply(@(x) mean(x, 'omitnan'), ...
        T.pixels_per_wavelength, G);
end

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function plot_metric_vs_dx(T, metric_var, fig_dir, file_name, title_text)

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18 10]);
models_M = unique(T.REQ_M, 'stable');
cs_values = unique(T.SIM_cs_bg, 'stable');
hold on;
for ci = 1:numel(cs_values)
    for mi = 1:numel(models_M)
        idx = T.SIM_cs_bg == cs_values(ci) & T.REQ_M == models_M(mi);
        Ti = sortrows(T(idx, :), 'SIM_dx');
        if isempty(Ti)
            continue;
        end
        plot(Ti.SIM_dx * 1e3, Ti.(metric_var), '-o', ...
            'LineWidth', 1.4, ...
            'DisplayName', sprintf('c_s=%g, M=%g', ...
            cs_values(ci), models_M(mi)));
    end
end
grid on;
xlabel('dx = dz (mm)');
ylabel(strrep(metric_var, '_', '\_'));
title(title_text, 'FontWeight', 'normal');
legend('Location', 'eastoutside');
exportgraphics(fig, fullfile(fig_dir, file_name), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function plot_mape_vs_pixels_per_wavelength(T, fig_dir)

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18 10]);
gscatter(T.pixels_per_wavelength, T.MAPE_pct, ...
    categorical("c_s=" + string(T.SIM_cs_bg) + ", M=" + string(T.REQ_M)));
grid on;
xlabel('pixels per wavelength');
ylabel('MAPE (%)');
title('MAPE vs pixels per wavelength', 'FontWeight', 'normal');
legend('Location', 'eastoutside');
exportgraphics(fig, fullfile(fig_dir, ...
    'level14_level01_mape_vs_pixels_per_wavelength.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function plot_mape_heatmap_M_dx(T, fig_dir)

cs_values = unique(T.SIM_cs_bg, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 9*numel(cs_values) 8]);
tl = tiledlayout(1, numel(cs_values), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for ci = 1:numel(cs_values)
    ax = nexttile(tl);
    Ti = T(T.SIM_cs_bg == cs_values(ci), :);
    M_values = unique(Ti.REQ_M, 'stable');
    dx_values = unique(Ti.SIM_dx, 'stable');
    Z = nan(numel(M_values), numel(dx_values));
    for mi = 1:numel(M_values)
        for di = 1:numel(dx_values)
            idx = Ti.REQ_M == M_values(mi) & Ti.SIM_dx == dx_values(di);
            Z(mi, di) = mean(Ti.MAPE_pct(idx), 'omitnan');
        end
    end
    imagesc(ax, dx_values * 1e3, M_values, Z);
    set(ax, 'YDir', 'normal');
    colorbar(ax);
    xlabel(ax, 'dx = dz (mm)');
    ylabel(ax, 'REQ M');
    title(ax, sprintf('c_s = %g m/s', cs_values(ci)), ...
        'FontWeight', 'normal');
end

title(tl, 'MAPE heatmap by REQ M and dx', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, ...
    'level14_level01_mape_heatmap_by_M_and_dx.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function y = prctile_finite(x, p)

x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = prctile(x, p);
end

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
