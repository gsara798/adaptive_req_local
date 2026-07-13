%% analyze_test_07_level12_sws_from_trained_q.m
% Level 12: Convert trained q predictions into shear-wave-speed predictions.
%
% This script takes the q predictions from Level 11 and converts them into
% predicted shear wave speed using the local REQ cumulative spectrum curve:
%
%   E_cum(k_q) = q_pred
%   cs_pred = 2*pi*f0/k_q

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load Level 11 output

experiment_name = 'test_08_advanced_angular_features';

[~, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

level11_dir = fullfile(PATHS.analysis_dir, ...
    'level_11_compare_predictor_families');

level11_file = fullfile(level11_dir, ...
    'level_11_compare_predictor_families.mat');

if ~exist(level11_file, 'file')
    error('Level 11 file not found:\n%s', level11_file);
end

S = load(level11_file);

T_pred = S.T_all_predictions;
T_ref = S.T_mc_eff;

mapping_var = choose_mapping_var(T_ref);

if mapping_var == ""
    error(['T_ref does not contain req_mapping or req_curve. ', ...
           'Rerun the MC sweep with StoreReqMapping = true.']);
end

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_12_sws_from_trained_q');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Convert q_pred into cs_pred

T_sws = add_sws_from_q_predictions(T_pred, T_ref, mapping_var);

%% Keep test split for main evaluation

T_test = T_sws(T_sws.split == "test", :);

%% Summarize SWS errors by model

T_sws_metrics = summarize_sws_metrics_by_model(T_test);

T_sws_metrics = sortrows(T_sws_metrics, 'RMSE_cs', 'ascend');

T_sws_metrics_all = summarize_sws_metrics_by_model(T_sws);
T_sws_metrics_all = sortrows(T_sws_metrics_all, 'RMSE_cs', 'ascend');

fprintf('\n============================================================\n');
fprintf('SWS METRICS SORTED BY RMSE_cs, TEST SPLIT\n');
fprintf('============================================================\n');

disp(T_sws_metrics(:, { ...
    'model_name', ...
    'model_type', ...
    'n', ...
    'MAE_cs', ...
    'RMSE_cs', ...
    'bias_cs', ...
    'MAE_pct', ...
    'RMSE_pct', ...
    'bias_pct', ...
    'R2_cs', ...
    'spearman_cs'}));

%% Sanity check using q_true

sanity_mae = mean(abs(T_test.cs_from_q_true - T_test.cs_true), 'omitnan');

fprintf('\nSanity check:\n');
fprintf('Mean abs(cs_from_q_true - cs_true) = %.6g m/s\n', sanity_mae);

%% Save tables

save(fullfile(analysis_dir, 'level_12_sws_from_trained_q.mat'), ...
    'T_sws', ...
    'T_test', ...
    'T_sws_metrics', ...
    'T_sws_metrics_all', ...
    '-v7.3');

T_sws_light = T_sws;

if ismember('req_curve', T_sws_light.Properties.VariableNames)
    T_sws_light.req_curve = [];
end

if ismember('req_mapping', T_sws_light.Properties.VariableNames)
    T_sws_light.req_mapping = [];
end

writetable(T_sws_light, ...
    fullfile(analysis_dir, 'sws_from_trained_q_predictions.csv'));

writetable(T_sws_metrics, ...
    fullfile(analysis_dir, 'sws_metrics_by_model.csv'));

writetable(T_sws_metrics_all, ...
    fullfile(analysis_dir, 'sws_metrics_by_model_all_splits.csv'));

save(fullfile(analysis_dir, 'level_12_sws_from_trained_q.mat'), ...
    'T_sws', ...
    'T_test', ...
    'T_sws_metrics', ...
    'T_sws_metrics_all', ...
    '-v7.3');

%% Make SWS plots

plot_sws_rmse_by_model(T_sws_metrics, analysis_dir);

plot_sws_pred_vs_true_best_models( ...
    T_test, ...
    T_sws_metrics, ...
    analysis_dir, ...
    6);

plot_sws_error_vs_true_speed_best_models( ...
    T_test, ...
    T_sws_metrics, ...
    analysis_dir, ...
    6);

fprintf('\nSaved Level 12 SWS analysis to:\n%s\n', analysis_dir);

fprintf('\nLevel 12 SWS analysis completed.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function T_out = add_sws_from_q_predictions(T_pred, T_ref, mapping_var)

key_vars = choose_key_vars(T_pred, T_ref);

T_curve = T_ref(:, cellstr([key_vars, mapping_var]));

[~, ia] = unique(T_curve(:, cellstr(key_vars)), 'rows', 'stable');
T_curve = T_curve(ia, :);

T_out = innerjoin( ...
    T_pred, ...
    T_curve, ...
    'Keys', cellstr(key_vars));

n = height(T_out);

cs_pred = NaN(n, 1);
k_pred = NaN(n, 1);

cs_from_q_true = NaN(n, 1);
k_from_q_true = NaN(n, 1);

for i = 1:n

    curve_i = get_curve(T_out.(char(mapping_var)), i);

    q_pred_i = T_out.q_pred(i);
    q_true_i = T_out.q_true(i);
    f0_i = T_out.SIM_f0(i);

    [cs_pred(i), k_pred(i)] = cs_from_q_curve(q_pred_i, curve_i, f0_i);
    [cs_from_q_true(i), k_from_q_true(i)] = cs_from_q_curve(q_true_i, curve_i, f0_i);

end

T_out.k_pred_from_q = k_pred;
T_out.cs_pred_from_q = cs_pred;

T_out.k_from_q_true = k_from_q_true;
T_out.cs_from_q_true = cs_from_q_true;

T_out.cs_true = T_out.SIM_cs_bg;

T_out.cs_error = T_out.cs_pred_from_q - T_out.cs_true;
T_out.cs_abs_error = abs(T_out.cs_error);

T_out.cs_rel_error = T_out.cs_error ./ T_out.cs_true;
T_out.cs_abs_rel_error = abs(T_out.cs_rel_error);

T_out.cs_error_pct = 100 * T_out.cs_rel_error;
T_out.cs_abs_error_pct = 100 * T_out.cs_abs_rel_error;

end

function mapping_var = choose_mapping_var(T)

vars = string(T.Properties.VariableNames);

if ismember("req_mapping", vars)
    mapping_var = "req_mapping";
elseif ismember("req_curve", vars)
    mapping_var = "req_curve";
else
    mapping_var = "";
end

end

function key_vars = choose_key_vars(T_pred, T_ref)

candidates = {
    ["condition_id", "step_idx", "realization_idx", "patch_idx"]
    ["condition_position", "step_idx", "realization_idx", "patch_idx"]
    ["step_idx", "realization_idx", "patch_idx", "SIM_f0", "SIM_cs_bg", "REQ_M"]
};

pred_vars = string(T_pred.Properties.VariableNames);
ref_vars = string(T_ref.Properties.VariableNames);

for i = 1:numel(candidates)

    key_i = candidates{i};

    if all(ismember(key_i, pred_vars)) && all(ismember(key_i, ref_vars))
        key_vars = key_i;
        return;
    end
end

error('Could not find common key variables between T_pred and T_ref.');

end

function curve = get_curve(curve_col, idx)

if iscell(curve_col)
    curve = curve_col{idx};
else
    curve = curve_col(idx);
end

if iscell(curve)
    curve = curve{1};
end

if isempty(curve) || ~isstruct(curve)
    error('Invalid req_curve at row %d.', idx);
end

end

function [cs, kq] = cs_from_q_curve(q, curve, f0)

cs = NaN;
kq = NaN;

if ~isfinite(q) || ~isfinite(f0) || f0 <= 0
    return;
end

[Ecum, k_cent] = extract_ecum_and_k(curve);

valid = isfinite(Ecum) & isfinite(k_cent) & k_cent > 0;

Ecum = Ecum(valid);
k_cent = k_cent(valid);

if numel(Ecum) < 2
    return;
end

[k_cent, idx_sort] = sort(k_cent, 'ascend');
Ecum = Ecum(idx_sort);

% Enforce monotonic cumulative curve.
Ecum = cummax(Ecum);

% Remove duplicate Ecum values because interp1 needs unique x values.
[Ecum_unique, idx_unique] = unique(Ecum, 'stable');
k_unique = k_cent(idx_unique);

if numel(Ecum_unique) < 2
    return;
end

q_clamped = min(max(q, min(Ecum_unique)), max(Ecum_unique));

kq = interp1( ...
    Ecum_unique, ...
    k_unique, ...
    q_clamped, ...
    'linear', ...
    'extrap');

if ~isfinite(kq) || kq <= 0
    return;
end

cs = 2*pi*f0/kq;

end

function [Ecum, k_cent] = extract_ecum_and_k(curve)

ecum_names = ["Ecum", "E_cum", "ecum", "cum_energy", "Ecum_norm"];
k_names = ["k_cent", "k_center", "k_centers", "k_radial", "k"];

Ecum = [];
k_cent = [];

for i = 1:numel(ecum_names)
    name_i = ecum_names(i);

    if isfield(curve, name_i)
        Ecum = curve.(name_i);
        break;
    end
end

for i = 1:numel(k_names)
    name_i = k_names(i);

    if isfield(curve, name_i)
        k_cent = curve.(name_i);
        break;
    end
end

if isempty(Ecum)
    error('Could not find cumulative energy field in req_curve.');
end

if isempty(k_cent)
    error('Could not find k-center field in req_curve.');
end

Ecum = double(Ecum(:));
k_cent = double(k_cent(:));

end

function T_metrics = summarize_sws_metrics_by_model(T)

[G, T_group] = findgroups(T(:, {'model_name', 'model_type'}));

n = splitapply(@numel, T.cs_error, G);

MAE_cs = splitapply(@mean_abs_finite, T.cs_error, G);
RMSE_cs = splitapply(@rmse_finite, T.cs_error, G);
bias_cs = splitapply(@mean_finite, T.cs_error, G);

MAE_pct = splitapply(@mean_abs_finite, T.cs_error_pct, G);
RMSE_pct = splitapply(@rmse_finite, T.cs_error_pct, G);
bias_pct = splitapply(@mean_finite, T.cs_error_pct, G);

R2_cs = splitapply(@r2_finite, T.cs_true, T.cs_pred_from_q, G);
spearman_cs = splitapply(@spearman_finite, T.cs_true, T.cs_pred_from_q, G);

T_metrics = T_group;

T_metrics.n = n;
T_metrics.MAE_cs = MAE_cs;
T_metrics.RMSE_cs = RMSE_cs;
T_metrics.bias_cs = bias_cs;

T_metrics.MAE_pct = MAE_pct;
T_metrics.RMSE_pct = RMSE_pct;
T_metrics.bias_pct = bias_pct;

T_metrics.R2_cs = R2_cs;
T_metrics.spearman_cs = spearman_cs;

end

function y = mean_finite(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = mean(x);
end

end

function y = mean_abs_finite(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = mean(abs(x));
end

end

function y = rmse_finite(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = sqrt(mean(x.^2));
end

end

function y = r2_finite(y_true, y_pred)

valid = isfinite(y_true) & isfinite(y_pred);

y_true = y_true(valid);
y_pred = y_pred(valid);

if numel(y_true) < 2
    y = NaN;
    return;
end

ss_res = sum((y_pred - y_true).^2);
ss_tot = sum((y_true - mean(y_true)).^2);

if ss_tot <= 0
    y = NaN;
else
    y = 1 - ss_res/ss_tot;
end

end

function rho = spearman_finite(x, y)

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

if numel(x) < 3
    rho = NaN;
else
    rho = corr(x, y, 'Type', 'Spearman');
end

end

function plot_sws_rmse_by_model(T_metrics, output_dir)

T_plot = sortrows(T_metrics, 'RMSE_cs', 'ascend');

labels = strcat(T_plot.model_name, newline, T_plot.model_type);
labels = categorical(labels);
labels = reordercats(labels, cellstr(string(labels)));

fig = figure('Color', 'w', 'Position', [100 100 1100 500]);
ax = axes(fig);

bar(ax, labels, T_plot.RMSE_cs);

ylabel(ax, 'SWS RMSE (m/s)');
xlabel(ax, 'Model');
title(ax, 'SWS error from trained q predictions');

grid(ax, 'on');
ax.FontSize = 13;
ax.LineWidth = 1.2;
xtickangle(ax, 45);

save_figure(fig, output_dir, 'sws_rmse_by_model');

end

function plot_sws_pred_vs_true_best_models(T_test, T_metrics, output_dir, n_best)

T_metrics = sortrows(T_metrics, 'RMSE_cs', 'ascend');
n_best = min(n_best, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1200 750]);
tl = tiledlayout(fig, 2, ceil(n_best/2), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_best

    model_name_i = T_metrics.model_name(i);
    model_type_i = T_metrics.model_type(i);

    idx = T_test.model_name == model_name_i & ...
          T_test.model_type == model_type_i;

    Ti = T_test(idx, :);

    ax = nexttile(tl);

    scatter(ax, Ti.cs_true, Ti.cs_pred_from_q, 12, 'filled', ...
        'MarkerFaceAlpha', 0.35);

    hold(ax, 'on');

    lims = [
        min([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')
        max([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')
    ];

    plot(ax, lims, lims, 'LineWidth', 1.2);

    xlim(ax, lims);
    ylim(ax, lims);

    xlabel(ax, 'True c_s (m/s)');
    ylabel(ax, 'Predicted c_s (m/s)');

    title(ax, sprintf('%s\n%s, RMSE = %.3f m/s', ...
        model_name_i, ...
        model_type_i, ...
        T_metrics.RMSE_cs(i)), ...
        'Interpreter', 'none');

    grid(ax, 'on');
    axis(ax, 'square');
    ax.FontSize = 11;
    ax.LineWidth = 1.1;
end

title(tl, 'Predicted SWS from trained q', 'FontWeight', 'bold');

save_figure(fig, output_dir, 'sws_pred_vs_true_best_models');

end

function plot_sws_error_vs_true_speed_best_models(T_test, T_metrics, output_dir, n_best)

T_metrics = sortrows(T_metrics, 'RMSE_cs', 'ascend');
n_best = min(n_best, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1200 750]);
tl = tiledlayout(fig, 2, ceil(n_best/2), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_best

    model_name_i = T_metrics.model_name(i);
    model_type_i = T_metrics.model_type(i);

    idx = T_test.model_name == model_name_i & ...
          T_test.model_type == model_type_i;

    Ti = T_test(idx, :);

    ax = nexttile(tl);

    scatter(ax, Ti.cs_true, Ti.cs_error, 12, 'filled', ...
        'MarkerFaceAlpha', 0.35);

    hold(ax, 'on');
    yline(ax, 0, 'LineWidth', 1.2);

    xlabel(ax, 'True c_s (m/s)');
    ylabel(ax, 'SWS error (m/s)');

    title(ax, sprintf('%s\n%s, bias = %.3f m/s', ...
        model_name_i, ...
        model_type_i, ...
        T_metrics.bias_cs(i)), ...
        'Interpreter', 'none');

    grid(ax, 'on');
    ax.FontSize = 11;
    ax.LineWidth = 1.1;
end

title(tl, 'SWS error from trained q', 'FontWeight', 'bold');

save_figure(fig, output_dir, 'sws_error_vs_true_speed_best_models');

end

function save_figure(fig, output_dir, base_name)

png_file = fullfile(output_dir, base_name + ".png");
pdf_file = fullfile(output_dir, base_name + ".pdf");

exportgraphics(fig, png_file, 'Resolution', 300);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');

end
