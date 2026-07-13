%% analyze_test_07_level13_deep_sws_performance.m
% Level 13: Deep analysis of SWS prediction accuracy, precision, and
% robustness across simulated conditions and Monte Carlo realizations.
%
% Adaptation of benchmark terminology:
%   Intramap: metrics within one condition, aperture step, and realization,
%             using local patches as spatial samples.
%   Interseed: variability of intramap metrics across realizations for one
%              fixed condition and aperture step.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load latest Level 12 lightweight predictions

experiment_name = 'test_08_advanced_angular_features';

[~, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

level12_dir = fullfile(PATHS.analysis_dir, ...
    'level_12_sws_from_trained_q');

prediction_file = fullfile(level12_dir, ...
    'sws_from_trained_q_predictions.csv');

if ~exist(prediction_file, 'file')
    error('Level 12 prediction CSV not found:\n%s', prediction_file);
end

T_all = readtable(prediction_file, 'TextType', 'string');

required_vars = [ ...
    "model_name", "model_type", "condition_id", "step_idx", ...
    "realization_idx", "patch_idx", "SIM_WaveModel", "SIM_f0", ...
    "SIM_cs_bg", "REQ_M", "Omega_sr", "q_true", "q_pred", ...
    "cs_true", "cs_pred_from_q", "cs_error_pct", "cs_abs_error_pct"];

validate_variables(T_all, required_vars);

if ismember("split", string(T_all.Properties.VariableNames))
    T = T_all(T_all.split == "test", :);
else
    T = T_all;
end

T_all.is_diagnostic_model = contains(T_all.model_name, "diagnostic", ...
    'IgnoreCase', true);
T_all.is_operational_model = ~T_all.is_diagnostic_model;

T.is_diagnostic_model = contains(T.model_name, "diagnostic", ...
    'IgnoreCase', true);
T.is_operational_model = ~T.is_diagnostic_model;

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_13_deep_sws_performance');

figure_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');

make_dir(analysis_dir);
make_dir(figure_dir);
make_dir(table_dir);

%% Build metrics

fprintf('\nComputing Level 13 intramap metrics...\n');
T_intramap = build_intramap_metrics(T);

fprintf('Computing Level 13 interseed metrics...\n');
T_interseed = build_interseed_metrics(T_intramap);

fprintf('Computing Level 13 condition metrics...\n');
T_condition = build_condition_metrics(T, T_intramap, T_interseed);

fprintf('Computing Level 13 global model ranking...\n');
T_model = build_model_metrics(T, T_intramap, T_interseed);

T_model_operational = T_model(T_model.is_operational_model, :);
T_model_diagnostic = T_model(T_model.is_diagnostic_model, :);

T_model_operational = sortrows(T_model_operational, ...
    {'MAPE_pct', 'median_intramap_CoV_pct'}, {'ascend', 'ascend'});

T_model_diagnostic = sortrows(T_model_diagnostic, ...
    {'MAPE_pct', 'median_intramap_CoV_pct'}, {'ascend', 'ascend'});

best_model_name = T_model_operational.model_name(1);
best_model_type = T_model_operational.model_type(1);

best_mask = T_condition.model_name == best_model_name & ...
    T_condition.model_type == best_model_type;

T_best_conditions = sortrows(T_condition(best_mask, :), ...
    'MAPE_pct', 'ascend');

T_worst_conditions = sortrows(T_condition(best_mask, :), ...
    'MAPE_pct', 'descend');

T_best_model_predictions = T( ...
    T.model_name == best_model_name & ...
    T.model_type == best_model_type, :);

T_best_model_predictions.q_error = ...
    T_best_model_predictions.q_pred - T_best_model_predictions.q_true;

T_best_model_predictions.q_abs_error = ...
    abs(T_best_model_predictions.q_error);

T_outliers = sortrows(T_best_model_predictions, ...
    'cs_abs_error_pct', 'descend');

outlier_threshold_pct = 20;
T_high_error_outliers = T_outliers( ...
    T_outliers.cs_abs_error_pct >= outlier_threshold_pct, :);

T_aperture_outlier_summary = build_aperture_outlier_summary( ...
    T_best_model_predictions, outlier_threshold_pct);

%% Scenario summaries for every model

factor_vars = [ ...
    "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M", ...
    "step_idx", "Omega_sr"];

T_factor_summary = table();

for i = 1:numel(factor_vars)
    T_i = build_factor_summary(T, factor_vars(i));
    T_factor_summary = [T_factor_summary; T_i]; %#ok<AGROW>
end

%% Descriptive metrics across all 54 conditions

fprintf('Computing descriptive metrics across all conditions...\n');

T_intramap_all = build_intramap_metrics(T_all);
T_interseed_all = build_interseed_metrics(T_intramap_all);
T_condition_all = build_condition_metrics(T_all, T_intramap_all, T_interseed_all);
T_model_all = build_model_metrics(T_all, T_intramap_all, T_interseed_all);

all_condition_dir = fullfile(table_dir, 'all_conditions');
make_dir(all_condition_dir);

writetable(T_intramap_all, fullfile(all_condition_dir, ...
    'intramap_metrics_all_conditions.csv'));
writetable(T_interseed_all, fullfile(all_condition_dir, ...
    'interseed_metrics_all_conditions.csv'));
writetable(T_condition_all, fullfile(all_condition_dir, ...
    'condition_metrics_all_conditions.csv'));
writetable(T_model_all, fullfile(all_condition_dir, ...
    'model_metrics_all_conditions.csv'));

%% Save tables and workspace

writetable(T_intramap, fullfile(table_dir, 'intramap_metrics.csv'));
writetable(T_interseed, fullfile(table_dir, 'interseed_metrics.csv'));
writetable(T_condition, fullfile(table_dir, 'condition_metrics.csv'));
writetable(T_model, fullfile(table_dir, 'model_metrics.csv'));
writetable(T_model_operational, fullfile(table_dir, ...
    'operational_model_ranking.csv'));
writetable(T_model_diagnostic, fullfile(table_dir, ...
    'diagnostic_model_ranking.csv'));
writetable(T_factor_summary, fullfile(table_dir, ...
    'factor_summary.csv'));
writetable(T_best_conditions, fullfile(table_dir, ...
    'best_model_conditions_ranked_best_to_worst.csv'));
writetable(T_worst_conditions, fullfile(table_dir, ...
    'best_model_conditions_ranked_worst_to_best.csv'));
writetable(T_outliers, fullfile(table_dir, ...
    'best_model_predictions_ranked_by_patch_error.csv'));
writetable(T_high_error_outliers, fullfile(table_dir, ...
    'best_model_high_error_outliers.csv'));
writetable(T_aperture_outlier_summary, fullfile(table_dir, ...
    'best_model_outlier_summary_by_aperture.csv'));

save(fullfile(analysis_dir, 'level_13_deep_sws_performance.mat'), ...
    'T_intramap', ...
    'T_interseed', ...
    'T_condition', ...
    'T_model', ...
    'T_model_operational', ...
    'T_model_diagnostic', ...
    'T_factor_summary', ...
    'T_best_conditions', ...
    'T_worst_conditions', ...
    'T_outliers', ...
    'T_high_error_outliers', ...
    'T_aperture_outlier_summary', ...
    'outlier_threshold_pct', ...
    'T_intramap_all', ...
    'T_interseed_all', ...
    'T_condition_all', ...
    'T_model_all', ...
    'best_model_name', ...
    'best_model_type', ...
    'PATHS', ...
    '-v7.3');

%% Figures

plot_model_accuracy_precision(T_model, figure_dir);
plot_best_worst_conditions(T_best_conditions, T_worst_conditions, ...
    best_model_name, best_model_type, figure_dir);
plot_best_model_factors(T, best_model_name, best_model_type, figure_dir);
plot_best_model_condition_mape_boxplots(T_condition, best_model_name, ...
    best_model_type, figure_dir);
plot_best_model_interseed(T_interseed, best_model_name, best_model_type, ...
    figure_dir);
plot_best_model_mape_vs_omega(T_interseed, best_model_name, best_model_type, ...
    figure_dir);
plot_outlier_aperture_diagnostics(T_best_model_predictions, ...
    outlier_threshold_pct, best_model_name, best_model_type, figure_dir);
plot_q_error_sws_error_diagnostics(T_best_model_predictions, ...
    outlier_threshold_pct, best_model_name, best_model_type, figure_dir);
plot_worst_outlier_q_predictions(T_outliers, best_model_name, ...
    best_model_type, figure_dir);

%% Text report

write_summary_report( ...
    fullfile(analysis_dir, 'level_13_summary.txt'), ...
    T_model_operational, ...
    T_model_diagnostic, ...
    T_best_conditions, ...
    T_worst_conditions);

fprintf('\n============================================================\n');
fprintf('LEVEL 13 COMPLETE\n');
fprintf('============================================================\n');
fprintf('Best operational model: %s / %s\n', ...
    best_model_name, best_model_type);
fprintf('Operational MAPE: %.3f %%\n', T_model_operational.MAPE_pct(1));
fprintf('Operational RMSE: %.3f %%\n', T_model_operational.RMSE_pct(1));
fprintf('Median intramap CoV: %.3f %%\n', ...
    T_model_operational.median_intramap_CoV_pct(1));
fprintf('Results saved to:\n%s\n', analysis_dir);

%% Local functions

function T_out = build_intramap_metrics(T)

group_vars = ["model_name", "model_type", "condition_id", ...
    "step_idx", "realization_idx"];

[G, T_group] = findgroups(T(:, cellstr(group_vars)));
ids = unique(G);
rows(numel(ids), 1) = empty_intramap_row();

for i = 1:numel(ids)
    Tg = T(G == ids(i), :);
    cs = Tg.cs_pred_from_q;
    cs_true = first_finite(Tg.cs_true);
    valid = isfinite(cs) & isfinite(cs_true) & cs_true > 0;

    cs_valid = cs(valid);
    err_pct = 100 * (cs_valid - cs_true) / cs_true;
    q_err = Tg.q_pred - Tg.q_true;
    q_err = q_err(isfinite(q_err));

    row = empty_intramap_row();
    row.n_patches = height(Tg);
    row.n_valid = nnz(valid);
    row.valid_fraction = nnz(valid) / max(height(Tg), 1);
    row.cs_true = cs_true;
    row.cs_pred_mean = mean_safe(cs_valid);
    row.cs_pred_median = median_safe(cs_valid);
    row.cs_pred_std = std_safe(cs_valid);
    row.cs_pred_iqr = iqr_safe(cs_valid);
    row.map_APE_mean_pct = abs_pct(row.cs_pred_mean, cs_true);
    row.map_APE_median_pct = abs_pct(row.cs_pred_median, cs_true);
    row.patch_MAPE_pct = mean_safe(abs(err_pct));
    row.patch_RMSE_pct = rmse_safe(err_pct);
    row.bias_pct = mean_safe(err_pct);
    row.CoV_pct = relative_std_pct(cs_valid);
    row.IQR_over_median_pct = relative_iqr_pct(cs_valid);
    row.q_MAE = mean_safe(abs(q_err));
    row.q_bias = mean_safe(q_err);
    row.failed = row.valid_fraction < 1;
    row.SIM_WaveModel = first_string(Tg.SIM_WaveModel);
    row.SIM_f0 = first_finite(Tg.SIM_f0);
    row.SIM_cs_bg = first_finite(Tg.SIM_cs_bg);
    row.REQ_M = first_finite(Tg.REQ_M);
    row.Omega_sr = first_finite(Tg.Omega_sr);
    row.is_diagnostic_model = first_logical(Tg.is_diagnostic_model);
    row.is_operational_model = first_logical(Tg.is_operational_model);
    rows(i) = row;
end

T_metrics = struct2table(rows);
T_out = [T_group, T_metrics];

end

function T_out = build_interseed_metrics(T)

group_vars = ["model_name", "model_type", "condition_id", "step_idx"];
[G, T_group] = findgroups(T(:, cellstr(group_vars)));
ids = unique(G);
rows(numel(ids), 1) = empty_interseed_row();

for i = 1:numel(ids)
    Tg = T(G == ids(i), :);
    ape = Tg.map_APE_mean_pct;
    cov_pct = Tg.CoV_pct;
    map_mean = Tg.cs_pred_mean;

    row = empty_interseed_row();
    row.n_seeds = numel(unique(Tg.realization_idx));
    row.median_APE_pct = median_safe(ape);
    row.IQR_APE_pct = iqr_safe(ape);
    row.std_APE_pct = std_safe(ape);
    row.p90_APE_pct = percentile_safe(ape, 90);
    row.p95_APE_pct = percentile_safe(ape, 95);
    row.worst_APE_pct = max_safe(ape);
    row.robust_relative_APE_variability_pct = robust_variability_pct(ape);
    row.median_CoV_pct = median_safe(cov_pct);
    row.IQR_CoV_pct = iqr_safe(cov_pct);
    row.interseed_cs_CoV_pct = relative_std_pct(map_mean);
    row.failure_rate_pct = 100 * mean(Tg.failed);
    row.SIM_WaveModel = first_string(Tg.SIM_WaveModel);
    row.SIM_f0 = first_finite(Tg.SIM_f0);
    row.SIM_cs_bg = first_finite(Tg.SIM_cs_bg);
    row.REQ_M = first_finite(Tg.REQ_M);
    row.Omega_sr = first_finite(Tg.Omega_sr);
    row.is_diagnostic_model = first_logical(Tg.is_diagnostic_model);
    row.is_operational_model = first_logical(Tg.is_operational_model);
    rows(i) = row;
end

T_metrics = struct2table(rows);
T_out = [T_group, T_metrics];

end

function T_out = build_condition_metrics(T, T_intramap, T_interseed)

group_vars = ["model_name", "model_type", "condition_id"];
[G, T_group] = findgroups(T(:, cellstr(group_vars)));
ids = unique(G);
rows(numel(ids), 1) = empty_condition_row();

for i = 1:numel(ids)
    Tg = T(G == ids(i), :);
    Ti = T_intramap( ...
        T_intramap.model_name == Tg.model_name(1) & ...
        T_intramap.model_type == Tg.model_type(1) & ...
        T_intramap.condition_id == Tg.condition_id(1), :);
    Ts = T_interseed( ...
        T_interseed.model_name == Tg.model_name(1) & ...
        T_interseed.model_type == Tg.model_type(1) & ...
        T_interseed.condition_id == Tg.condition_id(1), :);

    abs_pct = Tg.cs_abs_error_pct;
    err_pct = Tg.cs_error_pct;

    row = empty_condition_row();
    row.n = height(Tg);
    row.n_valid = nnz(isfinite(abs_pct));
    row.valid_fraction = row.n_valid / max(row.n, 1);
    row.MAPE_pct = mean_safe(abs_pct);
    row.MdAPE_pct = median_safe(abs_pct);
    row.RMSE_pct = rmse_safe(err_pct);
    row.bias_pct = mean_safe(err_pct);
    row.p90_abs_error_pct = percentile_safe(abs_pct, 90);
    row.p95_abs_error_pct = percentile_safe(abs_pct, 95);
    row.max_abs_error_pct = max_safe(abs_pct);
    row.median_intramap_APE_pct = median_safe(Ti.map_APE_mean_pct);
    row.median_intramap_CoV_pct = median_safe(Ti.CoV_pct);
    row.median_interseed_IQR_APE_pct = median_safe(Ts.IQR_APE_pct);
    row.median_interseed_cs_CoV_pct = median_safe(Ts.interseed_cs_CoV_pct);
    row.failure_rate_pct = 100 * mean(Ti.failed);
    row.SIM_WaveModel = first_string(Tg.SIM_WaveModel);
    row.SIM_f0 = first_finite(Tg.SIM_f0);
    row.SIM_cs_bg = first_finite(Tg.SIM_cs_bg);
    row.REQ_M = first_finite(Tg.REQ_M);
    row.is_diagnostic_model = first_logical(Tg.is_diagnostic_model);
    row.is_operational_model = first_logical(Tg.is_operational_model);
    rows(i) = row;
end

T_metrics = struct2table(rows);
T_out = [T_group, T_metrics];

end

function T_out = build_model_metrics(T, T_intramap, T_interseed)

group_vars = ["model_name", "model_type"];
[G, T_group] = findgroups(T(:, cellstr(group_vars)));
ids = unique(G);
rows(numel(ids), 1) = empty_model_row();

for i = 1:numel(ids)
    Tg = T(G == ids(i), :);
    Ti = T_intramap( ...
        T_intramap.model_name == Tg.model_name(1) & ...
        T_intramap.model_type == Tg.model_type(1), :);
    Ts = T_interseed( ...
        T_interseed.model_name == Tg.model_name(1) & ...
        T_interseed.model_type == Tg.model_type(1), :);

    row = empty_model_row();
    row.n = height(Tg);
    row.MAPE_pct = mean_safe(Tg.cs_abs_error_pct);
    row.MdAPE_pct = median_safe(Tg.cs_abs_error_pct);
    row.RMSE_pct = rmse_safe(Tg.cs_error_pct);
    row.bias_pct = mean_safe(Tg.cs_error_pct);
    row.p90_abs_error_pct = percentile_safe(Tg.cs_abs_error_pct, 90);
    row.p95_abs_error_pct = percentile_safe(Tg.cs_abs_error_pct, 95);
    row.median_intramap_APE_pct = median_safe(Ti.map_APE_mean_pct);
    row.median_intramap_CoV_pct = median_safe(Ti.CoV_pct);
    row.median_interseed_IQR_APE_pct = median_safe(Ts.IQR_APE_pct);
    row.median_interseed_cs_CoV_pct = median_safe(Ts.interseed_cs_CoV_pct);
    row.failure_rate_pct = 100 * mean(Ti.failed);
    row.is_diagnostic_model = first_logical(Tg.is_diagnostic_model);
    row.is_operational_model = first_logical(Tg.is_operational_model);
    rows(i) = row;
end

T_metrics = struct2table(rows);
T_out = [T_group, T_metrics];
T_out = sortrows(T_out, {'is_diagnostic_model', 'MAPE_pct'}, ...
    {'ascend', 'ascend'});

end

function T_out = build_factor_summary(T, factor_var)

group_vars = ["model_name", "model_type", factor_var];
[G, T_group] = findgroups(T(:, cellstr(group_vars)));
ids = unique(G);

metric_name = repmat(factor_var, numel(ids), 1);
factor_value = strings(numel(ids), 1);
n = zeros(numel(ids), 1);
MAPE_pct = NaN(numel(ids), 1);
MdAPE_pct = NaN(numel(ids), 1);
RMSE_pct = NaN(numel(ids), 1);
bias_pct = NaN(numel(ids), 1);
p90_abs_error_pct = NaN(numel(ids), 1);

for i = 1:numel(ids)
    Tg = T(G == ids(i), :);
    factor_value(i) = value_to_string(Tg.(char(factor_var))(1));
    n(i) = height(Tg);
    MAPE_pct(i) = mean_safe(Tg.cs_abs_error_pct);
    MdAPE_pct(i) = median_safe(Tg.cs_abs_error_pct);
    RMSE_pct(i) = rmse_safe(Tg.cs_error_pct);
    bias_pct(i) = mean_safe(Tg.cs_error_pct);
    p90_abs_error_pct(i) = percentile_safe(Tg.cs_abs_error_pct, 90);
end

T_out = T_group(:, {'model_name', 'model_type'});
T_out.factor_name = metric_name;
T_out.factor_value = factor_value;
T_out.n = n;
T_out.MAPE_pct = MAPE_pct;
T_out.MdAPE_pct = MdAPE_pct;
T_out.RMSE_pct = RMSE_pct;
T_out.bias_pct = bias_pct;
T_out.p90_abs_error_pct = p90_abs_error_pct;

end

function T_out = build_aperture_outlier_summary(T, threshold_pct)

group_vars = ["step_idx", "Omega_sr"];
[G, T_group] = findgroups(T(:, cellstr(group_vars)));

n = splitapply(@numel, T.cs_abs_error_pct, G);
MAPE_pct = splitapply(@mean_safe, T.cs_abs_error_pct, G);
p95_abs_error_pct = splitapply(@(x) percentile_safe(x, 95), ...
    T.cs_abs_error_pct, G);
max_abs_error_pct = splitapply(@max_safe, T.cs_abs_error_pct, G);
q_MAE = splitapply(@mean_safe, T.q_abs_error, G);
q_bias = splitapply(@mean_safe, T.q_error, G);
outlier_count = splitapply(@(x) sum(x >= threshold_pct), ...
    T.cs_abs_error_pct, G);
outlier_rate_pct = 100 * outlier_count ./ n;

T_out = T_group;
T_out.n = n;
T_out.MAPE_pct = MAPE_pct;
T_out.p95_abs_error_pct = p95_abs_error_pct;
T_out.max_abs_error_pct = max_abs_error_pct;
T_out.q_MAE = q_MAE;
T_out.q_bias = q_bias;
T_out.outlier_count = outlier_count;
T_out.outlier_rate_pct = outlier_rate_pct;

end

function plot_model_accuracy_precision(T, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1050 650]);
ax = axes(fig);
hold(ax, 'on');

operational = T.is_operational_model;

scatter(ax, T.MAPE_pct(operational), T.median_intramap_CoV_pct(operational), ...
    85, 'filled', 'DisplayName', 'Operational');
scatter(ax, T.MAPE_pct(~operational), T.median_intramap_CoV_pct(~operational), ...
    85, 'filled', 'Marker', 's', 'DisplayName', 'Diagnostic');

for i = 1:height(T)
    label = short_model_label(T.model_name(i), T.model_type(i));
    text(ax, T.MAPE_pct(i), T.median_intramap_CoV_pct(i), ...
        "  " + label, 'FontSize', 8);
end

xlabel(ax, 'Global SWS MAPE (%)');
ylabel(ax, 'Median intramap CoV (%)');
title(ax, 'Accuracy vs intramap precision');
legend(ax, 'Location', 'best');
grid(ax, 'on');
box(ax, 'on');

save_fig(fig, output_dir, 'model_accuracy_vs_precision');

end

function plot_best_worst_conditions(T_best, T_worst, model_name, model_type, output_dir)

n_show = min(10, height(T_best));
T_plot = [T_best(1:n_show, :); flipud(T_worst(1:n_show, :))];

labels = condition_short_label(T_plot);
colors = [repmat([0.20 0.60 0.35], n_show, 1); ...
    repmat([0.80 0.25 0.25], n_show, 1)];

fig = figure('Color', 'w', 'Position', [100 100 1200 750]);
ax = axes(fig);
b = barh(ax, T_plot.MAPE_pct);
b.FaceColor = 'flat';
b.CData = colors;

yticks(ax, 1:height(T_plot));
yticklabels(ax, labels);
xlabel(ax, 'SWS MAPE (%)');
title(ax, sprintf('Best and worst conditions: %s / %s', ...
    model_name, model_type), 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'on');

save_fig(fig, output_dir, 'best_operational_model_best_worst_conditions');

end

function plot_best_model_factors(T, model_name, model_type, output_dir)

T = T(T.model_name == model_name & T.model_type == model_type, :);

fig = figure('Color', 'w', 'Position', [100 100 1350 800]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

factor_names = ["SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M"];
titles = ["Wave model", "Frequency", "True SWS", "REQ M"];

for i = 1:numel(factor_names)
    ax = nexttile(tl);
    x = categorical(string(T.(char(factor_names(i)))));
    boxchart(ax, x, T.cs_abs_error_pct);
    ylabel(ax, 'Patch-level absolute SWS error (%)');
    title(ax, titles(i));
    grid(ax, 'on');
end

title(tl, sprintf(['Patch-level error distributions with outliers: ', ...
    '%s / %s'], model_name, model_type), ...
    'Interpreter', 'none');

save_fig(fig, output_dir, 'best_operational_model_factor_boxplots');

end

function plot_best_model_condition_mape_boxplots(T, model_name, model_type, output_dir)

T = T(T.model_name == model_name & T.model_type == model_type, :);

fig = figure('Color', 'w', 'Position', [100 100 1350 800]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

factor_names = ["SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M"];
titles = ["Wave model", "Frequency", "True SWS", "REQ M"];

for i = 1:numel(factor_names)
    ax = nexttile(tl);
    x = categorical(string(T.(char(factor_names(i)))));
    boxchart(ax, x, T.MAPE_pct);
    ylabel(ax, 'Condition-aggregated SWS MAPE (%)');
    title(ax, titles(i));
    grid(ax, 'on');
end

title(tl, sprintf('Condition-level MAPE distributions: %s / %s', ...
    model_name, model_type), 'Interpreter', 'none');

save_fig(fig, output_dir, 'best_operational_model_condition_mape_boxplots');

end

function plot_best_model_interseed(T, model_name, model_type, output_dir)

T = T(T.model_name == model_name & T.model_type == model_type, :);

fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
ax = axes(fig);

scatter(ax, T.median_APE_pct, T.IQR_APE_pct, ...
    45, T.Omega_sr, 'filled');

xlabel(ax, 'Median intramap APE across seeds (%)');
ylabel(ax, 'Interseed IQR of intramap APE (%)');
title(ax, 'Accuracy vs interseed robustness');
cb = colorbar(ax);
cb.Label.String = 'Omega (sr)';
grid(ax, 'on');
box(ax, 'on');

save_fig(fig, output_dir, 'best_operational_model_interseed_robustness');

end

function plot_best_model_mape_vs_omega(T, model_name, model_type, output_dir)

T = T(T.model_name == model_name & T.model_type == model_type, :);

fig = figure('Color', 'w', 'Position', [100 100 1000 650]);
ax = axes(fig);
hold(ax, 'on');

M_values = unique(T.REQ_M);
colors = lines(numel(M_values));

for i = 1:numel(M_values)
    Ti = T(T.REQ_M == M_values(i), :);
    [omega, ~, G] = unique(Ti.Omega_sr);
    y = splitapply(@mean_safe, Ti.median_APE_pct, G);
    plot(ax, omega, y, '-o', 'LineWidth', 1.6, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('M = %.3g', M_values(i)));
end

xlabel(ax, 'Omega (sr)');
ylabel(ax, 'Condition SWS MAPE (%)');
title(ax, 'MAPE vs aperture for best operational model');
legend(ax, 'Location', 'best');
grid(ax, 'on');
box(ax, 'on');

save_fig(fig, output_dir, 'best_operational_model_mape_vs_omega');

end

function plot_outlier_aperture_diagnostics(T, threshold_pct, model_name, model_type, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1300 800]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
scatter(ax, T.Omega_sr, T.cs_abs_error_pct, 24, T.REQ_M, 'filled', ...
    'MarkerFaceAlpha', 0.55);
yline(ax, threshold_pct, '--r', 'High-error threshold');
xlabel(ax, 'Omega (sr)');
ylabel(ax, 'Patch-level absolute SWS error (%)');
title(ax, 'Patch error vs aperture');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

ax = nexttile(tl);
[step_values, ~, G] = unique(T.step_idx);
MAPE = splitapply(@mean_safe, T.cs_abs_error_pct, G);
p95 = splitapply(@(x) percentile_safe(x, 95), T.cs_abs_error_pct, G);
max_error = splitapply(@max_safe, T.cs_abs_error_pct, G);
plot(ax, step_values, MAPE, '-o', 'LineWidth', 1.7, ...
    'DisplayName', 'MAPE');
hold(ax, 'on');
plot(ax, step_values, p95, '-s', 'LineWidth', 1.7, ...
    'DisplayName', 'P95');
plot(ax, step_values, max_error, '-^', 'LineWidth', 1.7, ...
    'DisplayName', 'Maximum');
xlabel(ax, 'Aperture step');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Mean, tail, and maximum error');
legend(ax, 'Location', 'best');
grid(ax, 'on');

ax = nexttile(tl);
outlier_count = splitapply(@(x) sum(x >= threshold_pct), ...
    T.cs_abs_error_pct, G);
bar(ax, step_values, outlier_count);
xlabel(ax, 'Aperture step');
ylabel(ax, sprintf('Count with error >= %.0f%%', threshold_pct));
title(ax, 'High-error outlier concentration');
grid(ax, 'on');

ax = nexttile(tl);
q_mae = splitapply(@mean_safe, T.q_abs_error, G);
q_bias = splitapply(@mean_safe, T.q_error, G);
plot(ax, step_values, q_mae, '-o', 'LineWidth', 1.7, ...
    'DisplayName', 'q MAE');
hold(ax, 'on');
plot(ax, step_values, q_bias, '-s', 'LineWidth', 1.7, ...
    'DisplayName', 'q bias');
yline(ax, 0, ':k', 'HandleVisibility', 'off');
xlabel(ax, 'Aperture step');
ylabel(ax, 'Quantile error');
title(ax, 'Quantile prediction error vs aperture');
legend(ax, 'Location', 'best');
grid(ax, 'on');

title(tl, sprintf('Outlier-aperture diagnostics: %s / %s', ...
    model_name, model_type), 'Interpreter', 'none');

save_fig(fig, output_dir, 'best_operational_model_outlier_aperture_diagnostics');

end

function plot_q_error_sws_error_diagnostics(T, threshold_pct, model_name, model_type, output_dir)

is_outlier = T.cs_abs_error_pct >= threshold_pct;

fig = figure('Color', 'w', 'Position', [100 100 1300 600]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
scatter(ax, T.q_error(~is_outlier), T.cs_error_pct(~is_outlier), ...
    20, T.Omega_sr(~is_outlier), 'filled', 'MarkerFaceAlpha', 0.35);
hold(ax, 'on');
scatter(ax, T.q_error(is_outlier), T.cs_error_pct(is_outlier), ...
    55, T.Omega_sr(is_outlier), 'filled', 'MarkerEdgeColor', 'k');
xline(ax, 0, ':k');
yline(ax, 0, ':k');
xlabel(ax, 'q_{pred} - q_{true}');
ylabel(ax, 'Signed SWS error (%)');
title(ax, 'Quantile error propagates into SWS error');
cb = colorbar(ax);
cb.Label.String = 'Omega (sr)';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T.q_true(~is_outlier), T.q_pred(~is_outlier), ...
    20, T.Omega_sr(~is_outlier), 'filled', 'MarkerFaceAlpha', 0.35);
hold(ax, 'on');
scatter(ax, T.q_true(is_outlier), T.q_pred(is_outlier), ...
    55, T.Omega_sr(is_outlier), 'filled', 'MarkerEdgeColor', 'k');
plot(ax, [0 1], [0 1], '--k', 'LineWidth', 1.2);
xlabel(ax, 'q true');
ylabel(ax, 'q predicted');
title(ax, sprintf('High-error patches highlighted (>= %.0f%%)', threshold_pct));
xlim(ax, [0 1]);
ylim(ax, [0 1]);
grid(ax, 'on');

title(tl, sprintf('Quantile-to-SWS failure mechanism: %s / %s', ...
    model_name, model_type), 'Interpreter', 'none');

save_fig(fig, output_dir, 'best_operational_model_q_error_vs_sws_error');

end

function plot_worst_outlier_q_predictions(T, model_name, model_type, output_dir)

n_show = min(20, height(T));
T = T(1:n_show, :);

labels = "C" + string(T.condition_id) + "-S" + string(T.step_idx) + ...
    "-R" + string(T.realization_idx) + "-P" + string(T.patch_idx);

fig = figure('Color', 'w', 'Position', [100 100 1250 750]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
x = 1:n_show;
plot(ax, x, T.q_true, '-o', 'LineWidth', 1.5, 'DisplayName', 'q true');
hold(ax, 'on');
plot(ax, x, T.q_pred, '-s', 'LineWidth', 1.5, 'DisplayName', 'q predicted');
xticks(ax, x);
xticklabels(ax, labels);
xtickangle(ax, 45);
ylabel(ax, 'Quantile');
title(ax, 'Quantiles for worst patch-level SWS outliers');
legend(ax, 'Location', 'best');
grid(ax, 'on');

ax = nexttile(tl);
b = bar(ax, x, T.cs_error_pct);
b.FaceColor = 'flat';
b.CData = turbo(n_show);
xticks(ax, x);
xticklabels(ax, labels);
xtickangle(ax, 45);
ylabel(ax, 'Signed SWS error (%)');
xlabel(ax, 'Condition-step-realization-patch');
title(ax, 'Resulting signed SWS error');
grid(ax, 'on');

title(tl, sprintf('Worst outlier details: %s / %s', ...
    model_name, model_type), 'Interpreter', 'none');

save_fig(fig, output_dir, 'best_operational_model_worst_outlier_q_predictions');

end

function write_summary_report(file_path, T_op, T_diag, T_best, T_worst)

fid = fopen(file_path, 'w');

if fid < 0
    warning('Could not write Level 13 summary report: %s', file_path);
    return;
end

cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'LEVEL 13 DEEP SWS PERFORMANCE SUMMARY\n\n');

fprintf(fid, 'Best operational model\n');
write_model_row(fid, T_op(1, :));

fprintf(fid, '\nBest diagnostic model\n');
write_model_row(fid, T_diag(1, :));

fprintf(fid, '\nFive best conditions for best operational model\n');
write_condition_rows(fid, T_best(1:min(5, height(T_best)), :));

fprintf(fid, '\nFive worst conditions for best operational model\n');
write_condition_rows(fid, T_worst(1:min(5, height(T_worst)), :));

end

function write_model_row(fid, T)

fprintf(fid, '%s / %s\n', T.model_name, T.model_type);
fprintf(fid, 'MAPE = %.4f %%\n', T.MAPE_pct);
fprintf(fid, 'RMSE = %.4f %%\n', T.RMSE_pct);
fprintf(fid, 'Bias = %.4f %%\n', T.bias_pct);
fprintf(fid, 'Median intramap CoV = %.4f %%\n', T.median_intramap_CoV_pct);
fprintf(fid, 'Median interseed IQR APE = %.4f %%\n', ...
    T.median_interseed_IQR_APE_pct);

end

function write_condition_rows(fid, T)

for i = 1:height(T)
    fprintf(fid, ['condition=%d, wave=%s, f0=%.0f, cs=%.3g, M=%.3g, ', ...
        'MAPE=%.4f%%, RMSE=%.4f%%, intramap_CoV=%.4f%%, ', ...
        'interseed_IQR_APE=%.4f%%\n'], ...
        T.condition_id(i), T.SIM_WaveModel(i), T.SIM_f0(i), ...
        T.SIM_cs_bg(i), T.REQ_M(i), T.MAPE_pct(i), T.RMSE_pct(i), ...
        T.median_intramap_CoV_pct(i), T.median_interseed_IQR_APE_pct(i));
end

end

function save_fig(fig, output_dir, base_name)

adaptive_req.analysis.save_analysis_figure( ...
    fig, output_dir, base_name, ...
    'SavePNG', true, ...
    'SavePDF', true, ...
    'CloseAfterSave', true);

end

function row = empty_intramap_row()

row = struct( ...
    'n_patches', 0, 'n_valid', 0, 'valid_fraction', NaN, ...
    'cs_true', NaN, 'cs_pred_mean', NaN, 'cs_pred_median', NaN, ...
    'cs_pred_std', NaN, 'cs_pred_iqr', NaN, ...
    'map_APE_mean_pct', NaN, 'map_APE_median_pct', NaN, ...
    'patch_MAPE_pct', NaN, 'patch_RMSE_pct', NaN, 'bias_pct', NaN, ...
    'CoV_pct', NaN, 'IQR_over_median_pct', NaN, ...
    'q_MAE', NaN, 'q_bias', NaN, 'failed', false, ...
    'SIM_WaveModel', "", 'SIM_f0', NaN, 'SIM_cs_bg', NaN, ...
    'REQ_M', NaN, 'Omega_sr', NaN, ...
    'is_diagnostic_model', false, 'is_operational_model', false);

end

function row = empty_interseed_row()

row = struct( ...
    'n_seeds', 0, 'median_APE_pct', NaN, 'IQR_APE_pct', NaN, ...
    'std_APE_pct', NaN, 'p90_APE_pct', NaN, 'p95_APE_pct', NaN, ...
    'worst_APE_pct', NaN, 'robust_relative_APE_variability_pct', NaN, ...
    'median_CoV_pct', NaN, 'IQR_CoV_pct', NaN, ...
    'interseed_cs_CoV_pct', NaN, 'failure_rate_pct', NaN, ...
    'SIM_WaveModel', "", 'SIM_f0', NaN, 'SIM_cs_bg', NaN, ...
    'REQ_M', NaN, 'Omega_sr', NaN, ...
    'is_diagnostic_model', false, 'is_operational_model', false);

end

function row = empty_condition_row()

row = struct( ...
    'n', 0, 'n_valid', 0, 'valid_fraction', NaN, ...
    'MAPE_pct', NaN, 'MdAPE_pct', NaN, 'RMSE_pct', NaN, ...
    'bias_pct', NaN, 'p90_abs_error_pct', NaN, ...
    'p95_abs_error_pct', NaN, 'max_abs_error_pct', NaN, ...
    'median_intramap_APE_pct', NaN, 'median_intramap_CoV_pct', NaN, ...
    'median_interseed_IQR_APE_pct', NaN, ...
    'median_interseed_cs_CoV_pct', NaN, 'failure_rate_pct', NaN, ...
    'SIM_WaveModel', "", 'SIM_f0', NaN, 'SIM_cs_bg', NaN, ...
    'REQ_M', NaN, ...
    'is_diagnostic_model', false, 'is_operational_model', false);

end

function row = empty_model_row()

row = struct( ...
    'n', 0, 'MAPE_pct', NaN, 'MdAPE_pct', NaN, 'RMSE_pct', NaN, ...
    'bias_pct', NaN, 'p90_abs_error_pct', NaN, ...
    'p95_abs_error_pct', NaN, ...
    'median_intramap_APE_pct', NaN, 'median_intramap_CoV_pct', NaN, ...
    'median_interseed_IQR_APE_pct', NaN, ...
    'median_interseed_cs_CoV_pct', NaN, 'failure_rate_pct', NaN, ...
    'is_diagnostic_model', false, 'is_operational_model', false);

end

function validate_variables(T, required)

missing = setdiff(required, string(T.Properties.VariableNames));

if ~isempty(missing)
    error('Missing required Level 12 variables: %s', strjoin(missing, ', '));
end

end

function make_dir(folder)

if ~exist(folder, 'dir')
    mkdir(folder);
end

end

function y = first_finite(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = x(1); end

end

function y = first_string(x)

x = string(x);
x = x(~ismissing(x));
if isempty(x), y = ""; else, y = x(1); end

end

function y = first_logical(x)

if isempty(x), y = false; else, y = logical(x(1)); end

end

function y = mean_safe(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(x); end

end

function y = median_safe(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = median(x); end

end

function y = std_safe(x)

x = x(isfinite(x));
if numel(x) < 2, y = NaN; else, y = std(x, 0); end

end

function y = iqr_safe(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = iqr(x); end

end

function y = max_safe(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = max(x); end

end

function y = percentile_safe(x, p)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = prctile(x, p); end

end

function y = rmse_safe(x)

x = x(isfinite(x));
if isempty(x), y = NaN; else, y = sqrt(mean(x.^2)); end

end

function y = abs_pct(pred, truth)

if ~isfinite(pred) || ~isfinite(truth) || truth == 0
    y = NaN;
else
    y = 100 * abs(pred - truth) / abs(truth);
end

end

function y = relative_std_pct(x)

mu = mean_safe(x);
sd = std_safe(x);
if ~isfinite(mu) || abs(mu) <= eps, y = NaN; else, y = 100 * sd / abs(mu); end

end

function y = relative_iqr_pct(x)

med = median_safe(x);
iqr_val = iqr_safe(x);
if ~isfinite(med) || abs(med) <= eps, y = NaN; else, y = 100 * iqr_val / abs(med); end

end

function y = robust_variability_pct(x)

med = median_safe(x);
iqr_val = iqr_safe(x);
if ~isfinite(med) || abs(med) <= eps, y = NaN; else, y = 100 * iqr_val / abs(med); end

end

function txt = value_to_string(value)

if isnumeric(value) || islogical(value)
    txt = string(value);
else
    txt = string(value);
end

end

function label = short_model_label(model_name, model_type)

label = replace(model_name, ...
    ["ModelA_local_entropy_only", "ModelB_entropy_known_params", ...
     "ModelC_rich_spectral_known_params", "ModelD_diagnostic_with_Meff"], ...
    ["A", "B", "C", "D"]);
label = label + "/" + replace(model_type, "_trees", "");

end

function labels = condition_short_label(T)

labels = "C" + string(T.condition_id) + " | " + ...
    T.SIM_WaveModel + " | f" + string(T.SIM_f0) + ...
    " | cs" + string(T.SIM_cs_bg) + " | M" + string(T.REQ_M);

end
