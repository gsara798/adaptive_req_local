%% analyze_test_08_level17_model_comparison_plots.m
% Level 17: visual comparison of Level 15 and Level 16 SWS models.
%
% This script does not retrain models. It reads Level 15/16 prediction CSVs
% and creates model-comparison plots with cleaner labels/subscripts.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load outputs

experiment_name = 'test_08_advanced_angular_features';

[~, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

level15_dir = fullfile(PATHS.analysis_dir, ...
    'level_15_ecum_features_and_plane_audit');
level16_dir = fullfile(PATHS.analysis_dir, ...
    'level_16_error_learning_feature_audit');
analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_17_model_comparison_plots');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

level15_file = fullfile(level15_dir, 'level15_sws_predictions.csv');
level16_file = fullfile(level16_dir, 'level16_sws_predictions.csv');

if ~exist(level15_file, 'file')
    error('Level 15 predictions not found:\n%s', level15_file);
end
if ~exist(level16_file, 'file')
    error('Level 16 predictions not found:\n%s', level16_file);
end

T15 = readtable(level15_file, 'TextType', 'string');
T16 = readtable(level16_file, 'TextType', 'string');

T15.analysis_level = repmat("Level 15", height(T15), 1);
T16.analysis_level = repmat("Level 16", height(T16), 1);

T_all = [standardize_sws_table(T15); standardize_sws_table(T16)];
T_all = T_all(T_all.split == "test", :);
T_all.model_display = model_display_name(T_all.model_name);

selected_models = [
    "ModelC_baseline"
    "ModelH_ecum_shape"
    "ModelI_angular_ecum_shape"
    "ModelJ_ecum_by_M"
    "Level16_base_ecum_srad_proxy"
    "Level16_residual_corrected"
];

T = T_all(ismember(T_all.model_name, selected_models) & ...
    T_all.model_type == "bagged_trees", :);

%% Summaries

T_metrics = summarize_sws_metrics(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "model_role"]);
T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');

T_by_M = summarize_sws_metrics(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "REQ_M"]);
T_by_cs = summarize_sws_metrics(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "SIM_cs_bg"]);
T_by_aperture = summarize_sws_metrics(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "step_idx", "Omega_sr"]);
T_by_wave = summarize_sws_metrics(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "SIM_WaveModel"]);

threshold_pct = 20;
T_high_error = summarize_high_error_rates(T, ...
    ["analysis_level", "model_name", "model_display", "model_type"], ...
    threshold_pct);
T_high_error_by_M = summarize_high_error_rates(T, ...
    ["analysis_level", "model_name", "model_display", "model_type", ...
     "REQ_M"], threshold_pct);

T_paired = build_level16_paired_delta(T);
T_outliers = build_worst_outliers(T, 30);

%% Save tables

writetable(T_metrics, fullfile(analysis_dir, ...
    'level17_model_comparison_metrics.csv'));
writetable(T_by_M, fullfile(analysis_dir, ...
    'level17_model_comparison_by_M.csv'));
writetable(T_by_cs, fullfile(analysis_dir, ...
    'level17_model_comparison_by_cs_bg.csv'));
writetable(T_by_aperture, fullfile(analysis_dir, ...
    'level17_model_comparison_by_aperture.csv'));
writetable(T_by_wave, fullfile(analysis_dir, ...
    'level17_model_comparison_by_wave_model.csv'));
writetable(T_high_error, fullfile(analysis_dir, ...
    'level17_high_error_rate_by_model.csv'));
writetable(T_high_error_by_M, fullfile(analysis_dir, ...
    'level17_high_error_rate_by_model_M.csv'));
writetable(T_paired, fullfile(analysis_dir, ...
    'level17_level16_paired_error_delta.csv'));
writetable(T_outliers, fullfile(analysis_dir, ...
    'level17_worst_outliers.csv'));

%% Figures

plot_model_metric_dashboard(T_metrics, T_high_error, threshold_pct, ...
    analysis_dir);
plot_pred_vs_true_grid(T, T_metrics, analysis_dir);
plot_error_boxplots(T, analysis_dir);
plot_error_by_M_and_cs(T_by_M, T_by_cs, analysis_dir);
plot_error_vs_aperture(T_by_aperture, T_high_error_by_M, ...
    threshold_pct, analysis_dir);
plot_q_diagnostics(T, T_metrics, analysis_dir);
plot_level16_paired_improvement(T_paired, analysis_dir);
plot_worst_outliers(T_outliers, analysis_dir);

save(fullfile(analysis_dir, 'level_17_model_comparison_plots.mat'), ...
    'T', 'T_metrics', 'T_by_M', 'T_by_cs', 'T_by_aperture', ...
    'T_by_wave', 'T_high_error', 'T_high_error_by_M', ...
    'T_paired', 'T_outliers', '-v7.3');

fprintf('\n============================================================\n');
fprintf('LEVEL 17 MODEL COMPARISON, TEST SPLIT\n');
fprintf('============================================================\n');
disp(T_metrics(:, {'model_display', 'n', 'MAE_pct', 'RMSE_pct', ...
    'p95_APE_pct', 'max_APE_pct'}));
fprintf('\nLevel 17 complete. Results saved to:\n%s\n', analysis_dir);

%% Local helpers

function T_out = standardize_sws_table(T_in)

n = height(T_in);
T_out = table();

string_vars = {'analysis_level', 'model_name', 'model_type', 'model_role', ...
    'split', 'SIM_WaveModel'};
numeric_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx', ...
    'SIM_f0', 'SIM_cs_bg', 'REQ_M', 'Omega_sr', 'q_true', 'q_pred', ...
    'q_base_pred', 'q_residual_pred', 'q_error', 'q_abs_error', ...
    'cs_true', 'cs_pred_from_q', 'cs_error_pct', 'cs_abs_error_pct'};

for i = 1:numel(string_vars)
    T_out.(string_vars{i}) = get_string_column(T_in, string_vars{i}, n);
end
for i = 1:numel(numeric_vars)
    T_out.(numeric_vars{i}) = get_numeric_column(T_in, numeric_vars{i}, n);
end

end

function x = get_string_column(T, name, n)

if ismember(name, T.Properties.VariableNames)
    x = string(T.(name));
else
    x = strings(n, 1);
end

end

function x = get_numeric_column(T, name, n)

if ismember(name, T.Properties.VariableNames)
    x = double(T.(name));
else
    x = NaN(n, 1);
end

end

function label = model_display_name(model_name)

model_name = string(model_name);
label = model_name;
label(model_name == "ModelC_baseline") = "C baseline";
label(model_name == "ModelH_ecum_shape") = "H Ecum";
label(model_name == "ModelI_angular_ecum_shape") = "I Angular + Ecum";
label(model_name == "ModelJ_ecum_by_M") = "J Ecum by M";
label(model_name == "Level16_base_ecum_srad_proxy") = ...
    "L16 base Ecum + Srad";
label(model_name == "Level16_residual_corrected") = ...
    "L16 residual corrected";

end

function T_metrics = summarize_sws_metrics(T, group_vars)

T = T(isfinite(T.cs_error_pct) & isfinite(T.cs_pred_from_q), :);
[G, T_metrics] = findgroups(T(:, cellstr(group_vars)));
T_metrics.n = splitapply(@numel, T.cs_error_pct, G);
T_metrics.MAE_pct = splitapply(@mean_abs_finite, T.cs_error_pct, G);
T_metrics.MAPE_pct = T_metrics.MAE_pct;
T_metrics.RMSE_pct = splitapply(@rmse_finite, T.cs_error_pct, G);
T_metrics.bias_pct = splitapply(@mean_finite, T.cs_error_pct, G);
T_metrics.median_APE_pct = splitapply(@median_abs_finite, ...
    T.cs_error_pct, G);
T_metrics.p95_APE_pct = splitapply(@p95_abs_finite, ...
    T.cs_error_pct, G);
T_metrics.max_APE_pct = splitapply(@max_abs_finite, ...
    T.cs_error_pct, G);
T_metrics.mean_abs_q_error = splitapply(@mean_finite, ...
    T.q_abs_error, G);

end

function T_rate = summarize_high_error_rates(T, group_vars, threshold_pct)

[G, T_rate] = findgroups(T(:, cellstr(group_vars)));
is_high = T.cs_abs_error_pct > threshold_pct;
T_rate.n = splitapply(@numel, T.cs_abs_error_pct, G);
T_rate.n_high_error = splitapply(@sum, double(is_high), G);
T_rate.high_error_pct = 100 * T_rate.n_high_error ./ T_rate.n;
T_rate.threshold_pct = repmat(threshold_pct, height(T_rate), 1);

end

function T_pair = build_level16_paired_delta(T)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_base = T(T.model_name == "Level16_base_ecum_srad_proxy", :);
T_corr = T(T.model_name == "Level16_residual_corrected", :);

if isempty(T_base) || isempty(T_corr)
    T_pair = table();
    return;
end

T_base = T_base(:, [key_vars, {'REQ_M', 'Omega_sr', 'SIM_cs_bg', ...
    'cs_abs_error_pct', 'q_abs_error'}]);
T_corr = T_corr(:, [key_vars, {'cs_abs_error_pct', 'q_abs_error', ...
    'q_base_pred', 'q_residual_pred'}]);
T_base.Properties.VariableNames{'cs_abs_error_pct'} = 'base_abs_sws_error_pct';
T_base.Properties.VariableNames{'q_abs_error'} = 'base_abs_q_error';
T_corr.Properties.VariableNames{'cs_abs_error_pct'} = 'corrected_abs_sws_error_pct';
T_corr.Properties.VariableNames{'q_abs_error'} = 'corrected_abs_q_error';

T_pair = innerjoin(T_base, T_corr, 'Keys', key_vars);
T_pair.delta_abs_sws_error_pct = ...
    T_pair.corrected_abs_sws_error_pct - T_pair.base_abs_sws_error_pct;
T_pair.delta_abs_q_error = ...
    T_pair.corrected_abs_q_error - T_pair.base_abs_q_error;

end

function T_out = build_worst_outliers(T, n_show)

T_best = T(T.model_name == "Level16_residual_corrected", :);
if isempty(T_best)
    T_best = T(T.model_name == T.model_name(1), :);
end
T_out = sortrows(T_best, 'cs_abs_error_pct', 'descend');
n_show = min(n_show, height(T_out));
T_out = T_out(1:n_show, :);

end

function plot_model_metric_dashboard(T_metrics, T_high_error, ...
    threshold_pct, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'descend');
labels = categorical(T_metrics.model_display);
labels = reordercats(labels, cellstr(T_metrics.model_display));

T_high = innerjoin(T_metrics(:, {'model_name', 'model_display'}), ...
    T_high_error(:, {'model_name', 'high_error_pct'}), 'Keys', 'model_name');

fig = figure('Color', 'w', 'Position', [100 100 1450 720]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
barh(ax, labels, T_metrics.RMSE_pct);
xlabel(ax, 'SWS RMSE (%)');
title(ax, 'Overall accuracy');
grid(ax, 'on');

ax = nexttile(tl);
barh(ax, labels, T_metrics.MAE_pct);
xlabel(ax, 'SWS MAPE (%)');
title(ax, 'Average absolute error');
grid(ax, 'on');

ax = nexttile(tl);
barh(ax, labels, T_metrics.p95_APE_pct);
xlabel(ax, '95th percentile APE (%)');
title(ax, 'Tail error');
grid(ax, 'on');

ax = nexttile(tl);
labels_high = categorical(T_high.model_display);
labels_high = reordercats(labels_high, cellstr(T_high.model_display));
barh(ax, labels_high, T_high.high_error_pct);
xlabel(ax, sprintf('APE > %.0f%% rate', threshold_pct));
title(ax, 'Large-error rate');
grid(ax, 'on');

title(tl, 'Level 17 model comparison dashboard', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'level17_model_metric_dashboard');

end

function plot_pred_vs_true_grid(T, T_metrics, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
n_show = min(6, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1450 850]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_show
    Ti = select_model(T, T_metrics.model_name(i));
    ax = nexttile(tl);
    scatter(ax, Ti.cs_true, Ti.cs_pred_from_q, 12, Ti.REQ_M, ...
        'filled', 'MarkerFaceAlpha', 0.35);
    hold(ax, 'on');
    lims = [min([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan'), ...
        max([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')];
    plot(ax, lims, lims, '--k', 'LineWidth', 1.1);
    xlim(ax, lims);
    ylim(ax, lims);
    axis(ax, 'square');
    xlabel(ax, 'True c_s (m/s)', 'Interpreter', 'tex');
    ylabel(ax, 'Predicted c_s (m/s)', 'Interpreter', 'tex');
    title(ax, sprintf('%s\nRMSE = %.2f%%', ...
        T_metrics.model_display(i), T_metrics.RMSE_pct(i)), ...
        'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = 'REQ M';
    grid(ax, 'on');
end

title(tl, 'Predicted vs true SWS', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'level17_sws_pred_vs_true_grid');

end

function plot_error_boxplots(T, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1450 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
boxchart(ax, categorical(T.model_display), T.cs_abs_error_pct);
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Patch-level error by model');
grid(ax, 'on');
xtickangle(ax, 35);

ax = nexttile(tl);
boxchart(ax, categorical(T.model_display), T.cs_error_pct);
yline(ax, 0, '--k');
ylabel(ax, 'Signed SWS error (%)');
title(ax, 'Bias distribution by model');
grid(ax, 'on');
xtickangle(ax, 35);

T_best = select_model(T, "Level16_residual_corrected");
ax = nexttile(tl);
boxchart(ax, categorical(string(T_best.REQ_M)), T_best.cs_abs_error_pct);
xlabel(ax, 'REQ M');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Best model by M');
grid(ax, 'on');

ax = nexttile(tl);
boxchart(ax, categorical(string(T_best.SIM_cs_bg)), T_best.cs_abs_error_pct);
xlabel(ax, 'True c_s (m/s)', 'Interpreter', 'tex');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Best model by true SWS');
grid(ax, 'on');

title(tl, 'Level 17 error distributions', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'level17_error_boxplots');

end

function plot_error_by_M_and_cs(T_by_M, T_by_cs, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1350 600]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
plot_grouped_metric(ax, T_by_M, "REQ_M", "RMSE_pct");
xlabel(ax, 'REQ M');
ylabel(ax, 'SWS RMSE (%)');
title(ax, 'Model comparison by M');
grid(ax, 'on');

ax = nexttile(tl);
plot_grouped_metric(ax, T_by_cs, "SIM_cs_bg", "RMSE_pct");
xlabel(ax, 'True c_s (m/s)', 'Interpreter', 'tex');
ylabel(ax, 'SWS RMSE (%)');
title(ax, 'Model comparison by true SWS');
grid(ax, 'on');

save_fig(fig, output_dir, 'level17_rmse_by_M_and_cs');

end

function plot_error_vs_aperture(T_by_aperture, T_high_error_by_M, ...
    threshold_pct, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1450 650]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
models = unique(T_by_aperture.model_name, 'stable');
hold(ax, 'on');
for i = 1:numel(models)
    Ti = T_by_aperture(T_by_aperture.model_name == models(i), :);
    Ti = sortrows(Ti, 'Omega_sr');
    plot(ax, Ti.Omega_sr, Ti.RMSE_pct, '-o', 'LineWidth', 1.2, ...
        'DisplayName', char(Ti.model_display(1)));
end
xlabel(ax, '\Omega (sr)', 'Interpreter', 'tex');
ylabel(ax, 'SWS RMSE (%)');
title(ax, 'Error vs aperture');
legend(ax, 'Location', 'best', 'Interpreter', 'none');
grid(ax, 'on');

ax = nexttile(tl);
plot_grouped_metric(ax, T_high_error_by_M, "REQ_M", "high_error_pct");
xlabel(ax, 'REQ M');
ylabel(ax, sprintf('APE > %.0f%% rate', threshold_pct));
title(ax, 'Large-error rate by M');
grid(ax, 'on');

save_fig(fig, output_dir, 'level17_aperture_and_high_error');

end

function plot_q_diagnostics(T, T_metrics, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
n_show = min(4, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1400 850]);
tl = tiledlayout(fig, 2, n_show, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:n_show
    Ti = select_model(T, T_metrics.model_name(i));

    ax = nexttile(tl, i);
    scatter(ax, Ti.q_true, Ti.q_pred, 12, Ti.REQ_M, ...
        'filled', 'MarkerFaceAlpha', 0.35);
    hold(ax, 'on');
    plot(ax, [0 1], [0 1], '--k');
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    xlabel(ax, 'q_{true}', 'Interpreter', 'tex');
    ylabel(ax, 'q_{pred}', 'Interpreter', 'tex');
    title(ax, T_metrics.model_display(i), 'Interpreter', 'none');
    grid(ax, 'on');

    ax = nexttile(tl, n_show + i);
    scatter(ax, Ti.q_error, Ti.cs_error_pct, 12, Ti.Omega_sr, ...
        'filled', 'MarkerFaceAlpha', 0.35);
    xline(ax, 0, ':k');
    yline(ax, 0, ':k');
    xlabel(ax, 'q_{pred} - q_{true}', 'Interpreter', 'tex');
    ylabel(ax, 'Signed SWS error (%)');
    title(ax, 'q error propagation');
    grid(ax, 'on');
end

title(tl, 'Quantile diagnostics with corrected subscript formatting', ...
    'FontWeight', 'bold');
save_fig(fig, output_dir, 'level17_q_diagnostics');

end

function plot_level16_paired_improvement(T_pair, output_dir)

if isempty(T_pair)
    return;
end

fig = figure('Color', 'w', 'Position', [100 100 1300 620]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
scatter(ax, T_pair.base_abs_sws_error_pct, ...
    T_pair.corrected_abs_sws_error_pct, 14, T_pair.REQ_M, ...
    'filled', 'MarkerFaceAlpha', 0.35);
hold(ax, 'on');
lims = [0, max([T_pair.base_abs_sws_error_pct; ...
    T_pair.corrected_abs_sws_error_pct], [], 'omitnan')];
plot(ax, lims, lims, '--k');
xlim(ax, lims);
ylim(ax, lims);
axis(ax, 'square');
xlabel(ax, 'Base absolute SWS error (%)');
ylabel(ax, 'Corrected absolute SWS error (%)');
title(ax, 'Points below line improved');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T_pair.Omega_sr, T_pair.delta_abs_sws_error_pct, 14, ...
    T_pair.REQ_M, 'filled', 'MarkerFaceAlpha', 0.35);
yline(ax, 0, '--k');
xlabel(ax, '\Omega (sr)', 'Interpreter', 'tex');
ylabel(ax, '\Delta absolute SWS error (%)', 'Interpreter', 'tex');
title(ax, 'Residual correction change vs aperture');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

title(tl, 'Level 16 residual correction: paired comparison', ...
    'FontWeight', 'bold');
save_fig(fig, output_dir, 'level17_level16_paired_improvement');

end

function plot_worst_outliers(T_outliers, output_dir)

if isempty(T_outliers)
    return;
end

n_show = min(20, height(T_outliers));
T = T_outliers(1:n_show, :);
labels = "C" + string(T.condition_id) + "-S" + string(T.step_idx) + ...
    "-R" + string(T.realization_idx) + "-P" + string(T.patch_idx);

fig = figure('Color', 'w', 'Position', [100 100 1300 760]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
x = 1:n_show;
plot(ax, x, T.q_true, '-o', 'LineWidth', 1.4, ...
    'DisplayName', 'q_{true}');
hold(ax, 'on');
plot(ax, x, T.q_pred, '-s', 'LineWidth', 1.4, ...
    'DisplayName', 'q_{pred}');
xticks(ax, x);
xticklabels(ax, labels);
xtickangle(ax, 45);
ylabel(ax, 'Quantile');
title(ax, 'Worst corrected-model outliers: quantiles', ...
    'Interpreter', 'tex');
legend(ax, 'Location', 'best', 'Interpreter', 'tex');
grid(ax, 'on');

ax = nexttile(tl);
bar(ax, x, T.cs_error_pct);
yline(ax, 0, '--k');
xticks(ax, x);
xticklabels(ax, labels);
xtickangle(ax, 45);
ylabel(ax, 'Signed SWS error (%)');
xlabel(ax, 'Condition-step-realization-patch');
title(ax, 'Worst corrected-model outliers: SWS error');
grid(ax, 'on');

save_fig(fig, output_dir, 'level17_worst_outliers');

end

function plot_grouped_metric(ax, T, x_var, y_var)

models = unique(T.model_name, 'stable');
x_values = unique(T.(char(x_var)), 'sorted');
Y = NaN(numel(x_values), numel(models));
legend_labels = strings(numel(models), 1);

for i = 1:numel(x_values)
    for j = 1:numel(models)
        mask = T.(char(x_var)) == x_values(i) & T.model_name == models(j);
        if any(mask)
            Y(i, j) = T.(char(y_var))(find(mask, 1));
            legend_labels(j) = T.model_display(find(mask, 1));
        end
    end
end

bar(ax, categorical(string(x_values)), Y);
legend(ax, legend_labels, 'Interpreter', 'none', 'Location', 'best');

end

function T_model = select_model(T, model_name)

T_model = T(T.model_name == string(model_name), :);

end

function y = mean_finite(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(x); end
end

function y = mean_abs_finite(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(abs(x)); end
end

function y = median_abs_finite(x)
x = abs(x(isfinite(x)));
if isempty(x), y = NaN; else, y = median(x); end
end

function y = p95_abs_finite(x)
x = abs(x(isfinite(x)));
if isempty(x), y = NaN; else, y = prctile(x, 95); end
end

function y = max_abs_finite(x)
x = abs(x(isfinite(x)));
if isempty(x), y = NaN; else, y = max(x); end
end

function y = rmse_finite(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = sqrt(mean(x.^2)); end
end

function save_fig(fig, output_dir, base_name)

exportgraphics(fig, fullfile(output_dir, base_name + ".png"), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, base_name + ".pdf"), ...
    'ContentType', 'vector');
close(fig);

end
