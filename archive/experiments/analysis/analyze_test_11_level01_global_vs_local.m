%% analyze_test_11_level01_global_vs_local.m
% Level 01 for Test 11: compare local, global, and hybrid REQ predictors.
%
% Models:
%   LocalOnly          : local patch spectral/Ecum features
%   GlobalOnly         : one global FFT/REQ descriptor per realization
%   HybridLocalGlobal  : local + global descriptors
%   GlobalQDiagnostic  : hybrid + q_global_theory, diagnostic only

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load Test 11

experiment_name = 'test_11_global_req_features';

[T_feat, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_01_global_vs_local');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
model_dir = fullfile(analysis_dir, 'models');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(model_dir);

required_vars = ["req_mapping", "global_req_mapping", ...
    "q_global_theory", "q_theory"];
assert(all(ismember(required_vars, string(T_feat.Properties.VariableNames))), ...
    'Test 11 data is missing required global/local REQ variables.');

T_feat = add_local_ecum_features_from_mapping(T_feat);
T_feat.row_key = make_row_key(T_feat);

%% Predictor sets

base_local = [
    "radial_entropy"
    "radial_peak_width"
    "radial_k_peak_norm"
    "radial_centroid_norm"
    "radial_std_norm"
    "ang_entropy"
    "ang_resultant_R"
    "ang_resultant_R2"
    "ang_moment_1"
    "ang_moment_2"
    "ang_moment_4"
    "ang_peak_count_rel"
    "ang_top2_to_top1"
    "ang_peak_separation_deg"
    "ecum_width_50_rel"
    "ecum_width_80_rel"
    "ecum_asymmetry_25_75"
    "ecum_width_ratio_80_50"
    "ecum_increment_entropy"
    "ecum_increment_peak_frac"
    "ecum_increment_gini"
    "ecum_slope_max"
    "ecum_slope_peak_to_mean"
    "ecum_slope_iqr_to_median"
    "srad_proxy_centroid_k_norm"
    "srad_proxy_std_k_norm"
    "srad_proxy_skewness"
    "srad_proxy_kurtosis"
    "srad_proxy_peak_k_norm"
    "srad_proxy_peak_to_centroid"
    "srad_proxy_low_side_frac"
    "srad_proxy_high_side_frac"
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"];

local_predictors = existing_numeric_predictors(T_feat, base_local);
global_predictors = existing_numeric_predictors(T_feat, [
    "global_radial_entropy"
    "global_radial_peak_width"
    "global_radial_k_peak_norm"
    "global_radial_centroid_norm"
    "global_radial_std_norm"
    "global_ang_entropy"
    "global_ang_resultant_R"
    "global_ang_resultant_R2"
    "global_ang_moment_1"
    "global_ang_moment_2"
    "global_ang_moment_4"
    "global_ang_peak_count_rel"
    "global_ang_top2_to_top1"
    "global_ang_peak_separation_deg"
    "global_ecum_width_50_rel"
    "global_ecum_width_80_rel"
    "global_ecum_asymmetry_25_75"
    "global_ecum_width_ratio_80_50"
    "global_ecum_increment_entropy"
    "global_ecum_increment_peak_frac"
    "global_ecum_increment_gini"
    "global_ecum_slope_max"
    "global_ecum_slope_peak_to_mean"
    "global_ecum_slope_iqr_to_median"
    "global_srad_proxy_centroid_k_norm"
    "global_srad_proxy_std_k_norm"
    "global_srad_proxy_skewness"
    "global_srad_proxy_kurtosis"
    "global_srad_proxy_peak_k_norm"
    "global_srad_proxy_peak_to_centroid"
    "global_srad_proxy_low_side_frac"
    "global_srad_proxy_high_side_frac"
    "REQ_M"
    "SIM_f0"
    "global_REQ_Nbins_effective"]);

hybrid_predictors = unique([local_predictors; global_predictors], 'stable');
diagnostic_predictors = unique([hybrid_predictors; "q_global_theory"], ...
    'stable');

model_specs = struct([]);
model_specs(1).name = "LocalOnly";
model_specs(1).predictors = local_predictors;
model_specs(1).role = "operational";
model_specs(2).name = "GlobalOnly";
model_specs(2).predictors = global_predictors;
model_specs(2).role = "operational";
model_specs(3).name = "HybridLocalGlobal";
model_specs(3).predictors = hybrid_predictors;
model_specs(3).role = "operational";
model_specs(4).name = "GlobalQDiagnostic";
model_specs(4).predictors = diagnostic_predictors;
model_specs(4).role = "diagnostic_only";

%% Train q models

T_all_pred = table();
T_q_metrics = table();
MODELS = struct([]);

for i = 1:numel(model_specs)

    fprintf('\n=== Training %s ===\n', model_specs(i).name);

    [MODEL_i, T_pred_i, T_metrics_i] = ...
        adaptive_req.analysis.train_q_model_from_predictors( ...
            T_feat, ...
            model_specs(i).predictors, ...
            'ModelName', model_specs(i).name, ...
            'SplitMode', 'condition', ...
            'TrainFraction', 0.70, ...
            'ModelTypes', ["linear", "boosted_trees", "bagged_trees"], ...
            'NumLearningCycles', 200, ...
            'MinLeafSize', 8, ...
            'Verbose', true);

    T_pred_i.model_role = repmat(model_specs(i).role, height(T_pred_i), 1);
    T_pred_i.row_key = make_row_key(T_pred_i);

    MODELS(i).name = model_specs(i).name;
    MODELS(i).role = model_specs(i).role;
    MODELS(i).predictors = model_specs(i).predictors;
    MODELS(i).model = MODEL_i;

    T_all_pred = concat_tables_with_missing(T_all_pred, T_pred_i);
    T_q_metrics = [T_q_metrics; T_metrics_i]; %#ok<AGROW>
end

%% Convert predicted q to local SWS and compare global q diagnostic

T_sws = add_sws_metrics(T_all_pred, T_feat, "req_mapping", "local");
T_global_q = build_global_q_baseline(T_feat);
T_sws = concat_tables_with_missing(T_sws, T_global_q);
T_sws = attach_local_global_q_gap(T_sws, T_feat);

T_sws_test = T_sws(T_sws.split == "test", :);

T_sws_metrics = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role"]);
T_sws_by_M = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"]);
T_sws_by_aperture = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "step_idx", "Omega_sr"]);
T_sws_by_wave = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "SIM_WaveModel"]);
T_high_error = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role"], 20);
T_high_error_by_M = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"], 20);
T_high_error_M_composition = summarize_high_error_composition(T_sws_test, ...
    ["model_name", "model_type", "model_role"], "REQ_M", 20);

T_q_global_summary = summarize_q_global_gap(T_feat);
T_worst = sortrows(T_sws_test, 'abs_sws_error_pct', 'descend');
T_worst = T_worst(1:min(50, height(T_worst)), :);

%% Save outputs

writetable(T_q_metrics, fullfile(table_dir, ...
    'level11_q_model_metrics.csv'));
writetable(remove_cell_columns(T_sws), fullfile(table_dir, ...
    'level11_sws_predictions.csv'));
writetable(T_sws_metrics, fullfile(table_dir, ...
    'level11_sws_metrics.csv'));
writetable(T_sws_by_M, fullfile(table_dir, ...
    'level11_sws_metrics_by_M.csv'));
writetable(T_sws_by_aperture, fullfile(table_dir, ...
    'level11_sws_metrics_by_aperture.csv'));
writetable(T_sws_by_wave, fullfile(table_dir, ...
    'level11_sws_metrics_by_wave_model.csv'));
writetable(T_high_error, fullfile(table_dir, ...
    'level11_high_error_gt20_by_model.csv'));
writetable(T_high_error_by_M, fullfile(table_dir, ...
    'level11_high_error_gt20_by_model_M.csv'));
writetable(T_high_error_M_composition, fullfile(table_dir, ...
    'level11_high_error_gt20_M_composition.csv'));
writetable(T_q_global_summary, fullfile(table_dir, ...
    'level11_q_local_global_gap_summary.csv'));
writetable(remove_cell_columns(T_worst), fullfile(table_dir, ...
    'level11_worst_sws_predictions.csv'));

save(fullfile(model_dir, 'level11_global_vs_local_models.mat'), ...
    'MODELS', 'T_q_metrics', 'T_sws_metrics', 'T_sws_by_M', ...
    'T_sws_by_aperture', 'T_sws_by_wave', 'T_high_error', ...
    'T_high_error_by_M', 'T_high_error_M_composition', ...
    'T_q_global_summary', ...
    'MC', 'PATHS', '-v7.3');

%% Figures

plot_metric_bar(T_sws_metrics, fig_dir);
plot_q_scatter(T_sws_test, fig_dir);
plot_sws_error_box(T_sws_test, fig_dir);
plot_sws_pred_scatter(T_sws_test, fig_dir);
plot_residual_histograms(T_sws_test, fig_dir);
plot_error_by_aperture(T_sws_by_aperture, fig_dir);
plot_error_by_M(T_sws_by_M, fig_dir);
plot_error_by_wave(T_sws_by_wave, fig_dir);
plot_error_vs_q_gap(T_sws_test, fig_dir);
plot_high_error_rates(T_high_error, T_high_error_by_M, ...
    T_high_error_M_composition, fig_dir);
plot_local_vs_global_q(T_feat, fig_dir);
plot_q_gap_by_M_and_aperture(T_feat, fig_dir);
plot_q_gap_heatmap(T_q_global_summary, fig_dir);

fprintf('\nTest 11 Level 01 analysis complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
disp(T_sws_metrics);

%% Local functions

function T = add_local_ecum_features_from_mapping(T)

if ~ismember('req_mapping', T.Properties.VariableNames)
    return;
end

feat0 = adaptive_req.quantile.extract_ecum_shape_features( ...
    T.req_mapping{find_first_mapping(T.req_mapping)});
names = fieldnames(feat0);

for j = 1:numel(names)
    if ~ismember(names{j}, T.Properties.VariableNames)
        T.(names{j}) = nan(height(T), 1);
    end
end

for i = 1:height(T)
    if isempty(T.req_mapping{i})
        continue;
    end
    feat_i = adaptive_req.quantile.extract_ecum_shape_features( ...
        T.req_mapping{i});
    for j = 1:numel(names)
        T.(names{j})(i) = feat_i.(names{j});
    end
end

end

function idx = find_first_mapping(C)

idx = find(~cellfun(@isempty, C), 1, 'first');
if isempty(idx)
    error('No non-empty REQ mapping found.');
end

end

function predictors = existing_numeric_predictors(T, candidates)

vars = string(T.Properties.VariableNames);
candidates = string(candidates(:));
predictors = strings(0, 1);

for i = 1:numel(candidates)
    name_i = candidates(i);
    if ismember(name_i, vars) && isnumeric(T.(char(name_i)))
        predictors(end + 1, 1) = name_i; %#ok<AGROW>
    end
end

end

function key = make_row_key(T)

parts = strings(height(T), 4);
parts(:, 1) = string(T.condition_id);
parts(:, 2) = string(T.step_idx);
parts(:, 3) = string(T.realization_idx);
parts(:, 4) = string(T.patch_idx);
key = join(parts, "|", 2);

end

function T_sws = add_sws_metrics(T_pred, T_ref, mapping_var, mapping_label)

[tf, loc] = ismember(T_pred.row_key, T_ref.row_key);
if ~all(tf)
    error('Could not match all prediction rows back to Test 11 reference rows.');
end

q_pred = T_pred.q_pred;
cs_pred = q_to_cs(q_pred, T_ref.(char(mapping_var))(loc), ...
    T_ref.SIM_f0(loc));

T_sws = T_pred;
T_sws.mapping_source = repmat(string(mapping_label), height(T_sws), 1);
T_sws.cs_true = T_ref.SIM_cs_bg(loc);
T_sws.cs_pred = cs_pred;
T_sws.sws_error = T_sws.cs_pred - T_sws.cs_true;
T_sws.abs_sws_error = abs(T_sws.sws_error);
T_sws.sws_error_pct = 100 * T_sws.sws_error ./ T_sws.cs_true;
T_sws.abs_sws_error_pct = abs(T_sws.sws_error_pct);

end

function T_sws = attach_local_global_q_gap(T_sws, T_ref)

if ~ismember('row_key', T_sws.Properties.VariableNames)
    T_sws.q_local_minus_global = nan(height(T_sws), 1);
    T_sws.abs_q_local_minus_global = nan(height(T_sws), 1);
    return;
end

[tf, loc] = ismember(string(T_sws.row_key), string(T_ref.row_key));
T_sws.q_local_minus_global = nan(height(T_sws), 1);
T_sws.abs_q_local_minus_global = nan(height(T_sws), 1);

if ismember('q_local_minus_global', T_ref.Properties.VariableNames)
    T_sws.q_local_minus_global(tf) = T_ref.q_local_minus_global(loc(tf));
    T_sws.abs_q_local_minus_global(tf) = ...
        abs(T_sws.q_local_minus_global(tf));
end

end

function T_global = build_global_q_baseline(T_ref)

T_global = table();
T_global.condition_id = T_ref.condition_id;
T_global.step_idx = T_ref.step_idx;
T_global.realization_idx = T_ref.realization_idx;
T_global.patch_idx = T_ref.patch_idx;
T_global.Omega_sr = T_ref.Omega_sr;
T_global.SIM_WaveModel = T_ref.SIM_WaveModel;
T_global.SIM_f0 = T_ref.SIM_f0;
T_global.SIM_cs_bg = T_ref.SIM_cs_bg;
T_global.REQ_M = T_ref.REQ_M;
T_global.row_key = T_ref.row_key;
T_global.model_name = repmat("GlobalQDirect", height(T_ref), 1);
T_global.model_type = repmat("global_q_oracle", height(T_ref), 1);
T_global.model_role = repmat("diagnostic_only", height(T_ref), 1);
T_global.split = repmat("test", height(T_ref), 1);
T_global.mapping_source = repmat("global", height(T_ref), 1);
T_global.q_true = T_ref.q_theory;
T_global.q_pred_raw = T_ref.q_global_theory;
T_global.q_pred = T_ref.q_global_theory;
T_global.residual = T_global.q_pred - T_global.q_true;
T_global.abs_error = abs(T_global.residual);
T_global.cs_true = T_ref.SIM_cs_bg;
T_global.cs_pred = q_to_cs(T_global.q_pred, T_ref.global_req_mapping, ...
    T_ref.SIM_f0);
T_global.sws_error = T_global.cs_pred - T_global.cs_true;
T_global.abs_sws_error = abs(T_global.sws_error);
T_global.sws_error_pct = 100 * T_global.sws_error ./ T_global.cs_true;
T_global.abs_sws_error_pct = abs(T_global.sws_error_pct);

end

function cs = q_to_cs(q, mappings, f0)

q = double(q(:));
f0 = double(f0(:));
cs = nan(numel(q), 1);

for i = 1:numel(q)
    mapping_i = mappings{i};
    if isempty(mapping_i) || ~isfinite(q(i))
        continue;
    end

    k = adaptive_req.quantile.quantile_to_k(mapping_i, q(i));
    cs(i) = 2*pi*f0(i) ./ k;
end

end

function T_sum = summarize_sws_metrics(T, group_vars)

group_vars = string(group_vars);
[G, T_keys] = findgroups(T(:, cellstr(group_vars)));
n = splitapply(@numel, T.abs_sws_error_pct, G);
mape = splitapply(@(x) mean(x, 'omitnan'), T.abs_sws_error_pct, G);
rmse = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.sws_error_pct, G);
medae = splitapply(@(x) median(x, 'omitnan'), T.abs_sws_error_pct, G);
p95 = splitapply(@(x) prctile(x, 95), T.abs_sws_error_pct, G);
maxe = splitapply(@(x) max(x, [], 'omitnan'), T.abs_sws_error_pct, G);
high20 = splitapply(@(x) 100 * mean(x > 20, 'omitnan'), ...
    T.abs_sws_error_pct, G);

T_sum = T_keys;
T_sum.N = n;
T_sum.MAPE_pct = mape;
T_sum.RMSE_pct = rmse;
T_sum.MedAE_pct = medae;
T_sum.P95_abs_error_pct = p95;
T_sum.Max_abs_error_pct = maxe;
T_sum.HighError_gt20_pct = high20;
T_sum = sortrows(T_sum, 'MAPE_pct', 'ascend');

end

function T_sum = summarize_high_error_rates(T, group_vars, threshold_pct)

group_vars = string(group_vars);
[G, T_keys] = findgroups(T(:, cellstr(group_vars)));
is_high = T.abs_sws_error_pct > threshold_pct;

T_sum = T_keys;
T_sum.N = splitapply(@numel, is_high, G);
T_sum.N_high_error = splitapply(@sum, is_high, G);
T_sum.HighError_gt20_pct = 100 * T_sum.N_high_error ./ T_sum.N;
T_sum.MAPE_all_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);

T_sum = sortrows(T_sum, 'HighError_gt20_pct', 'descend');

end

function T_comp = summarize_high_error_composition( ...
    T, base_group_vars, composition_var, threshold_pct)

base_group_vars = string(base_group_vars);
composition_var = string(composition_var);

T_high = T(T.abs_sws_error_pct > threshold_pct, :);

if isempty(T_high)
    T_comp = table();
    return;
end

group_vars = [base_group_vars(:); composition_var];
[G, T_keys] = findgroups(T_high(:, cellstr(group_vars)));

n_high = splitapply(@numel, T_high.abs_sws_error_pct, G);
T_comp = T_keys;
T_comp.N_high_error = n_high;

[Gbase, T_base] = findgroups(T_high(:, cellstr(base_group_vars)));
total_high = splitapply(@numel, T_high.abs_sws_error_pct, Gbase);
base_key = make_group_key(T_base, base_group_vars);
row_key = make_group_key(T_comp, base_group_vars);

T_comp.Total_high_error_model = nan(height(T_comp), 1);
for i = 1:height(T_comp)
    idx = find(base_key == row_key(i), 1);
    T_comp.Total_high_error_model(i) = total_high(idx);
end

T_comp.Pct_of_model_high_errors = ...
    100 * T_comp.N_high_error ./ T_comp.Total_high_error_model;
T_comp = sortrows(T_comp, 'Pct_of_model_high_errors', 'descend');

end

function key = make_group_key(T, vars)

vars = string(vars);
key = strings(height(T), 1);

for i = 1:numel(vars)
    key = key + "|" + string(T.(char(vars(i))));
end

end

function T_gap = summarize_q_global_gap(T)

T.q_abs_local_global_gap = abs(T.q_local_minus_global);
T_gap = summarize_generic(T, ...
    ["REQ_M", "SIM_WaveModel", "step_idx"], ...
    "q_abs_local_global_gap");

end

function T_sum = summarize_generic(T, group_vars, value_var)

[G, T_keys] = findgroups(T(:, cellstr(group_vars)));
x = T.(char(value_var));
T_sum = T_keys;
T_sum.N = splitapply(@numel, x, G);
T_sum.Mean = splitapply(@(v) mean(v, 'omitnan'), x, G);
T_sum.Median = splitapply(@(v) median(v, 'omitnan'), x, G);
T_sum.P95 = splitapply(@(v) prctile(v, 95), x, G);
T_sum.Max = splitapply(@(v) max(v, [], 'omitnan'), x, G);

end

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

end

function T = concat_tables_with_missing(A, B)

vars_all = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');

A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);

T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];

end

function T = add_missing_columns(T, vars_all)

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    name_i = char(vars_all(i));
    if ismember(vars_all(i), vars)
        continue;
    end

    string_like = any(endsWith(vars_all(i), ...
        ["name", "type", "role", "source", "label", "split"])) || ...
        startsWith(vars_all(i), "SIM_WaveModel");

    if string_like
        T.(name_i) = strings(height(T), 1);
    else
        T.(name_i) = nan(height(T), 1);
    end
end

end

function plot_metric_bar(T_metrics, fig_dir)

T = T_metrics(T_metrics.model_type == "bagged_trees", :);
figure('Color', 'w');
bar(categorical(T.model_name), T.MAPE_pct);
ylabel('MAPE SWS (%)');
title('Test 11: local vs global vs hybrid (bagged trees)');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_mape_model_bar.png'));
close(gcf);

end

function plot_q_scatter(T, fig_dir)

T = T(T.model_type == "bagged_trees", :);
figure('Color', 'w');
tiledlayout('flow');
models = unique(T.model_name, 'stable');
for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    scatter(Ti.q_true, Ti.q_pred, 12, Ti.REQ_M, 'filled');
    hold on;
    plot([0 1], [0 1], 'k--');
    axis equal; xlim([0 1]); ylim([0 1]);
    xlabel('q true local');
    ylabel('q predicted');
    title(models(i), 'Interpreter', 'none');
    grid on;
end
saveas(gcf, fullfile(fig_dir, 'level11_q_true_vs_pred.png'));
close(gcf);

end

function plot_sws_error_box(T, fig_dir)

T = T(T.model_type == "bagged_trees" | T.model_type == "global_q_oracle", :);
figure('Color', 'w');
boxchart(categorical(T.model_name), T.abs_sws_error_pct);
yline(20, 'r--', '20%');
ylabel('|SWS error| (%)');
title('Test 11: SWS error distribution');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_sws_error_boxplot.png'));
close(gcf);

end

function plot_sws_pred_scatter(T, fig_dir)

T = T(T.model_type == "bagged_trees" | T.model_type == "global_q_oracle", :);
figure('Color', 'w');
tiledlayout('flow');
models = unique(T.model_name, 'stable');

for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    scatter(Ti.cs_true, Ti.cs_pred, 14, Ti.abs_sws_error_pct, 'filled');
    hold on;
    lo = min([Ti.cs_true; Ti.cs_pred], [], 'omitnan');
    hi = max([Ti.cs_true; Ti.cs_pred], [], 'omitnan');
    plot([lo hi], [lo hi], 'k--');
    axis equal;
    xlabel('True SWS (m/s)');
    ylabel('Predicted SWS (m/s)');
    title(models(i), 'Interpreter', 'none');
    grid on;
    cb = colorbar;
    cb.Label.String = '|error| (%)';
end

saveas(gcf, fullfile(fig_dir, 'level11_sws_true_vs_pred.png'));
close(gcf);

end

function plot_residual_histograms(T, fig_dir)

T = T(T.model_type == "bagged_trees" | T.model_type == "global_q_oracle", :);
figure('Color', 'w');
tiledlayout('flow');
models = unique(T.model_name, 'stable');

for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    histogram(Ti.sws_error_pct, 35, 'Normalization', 'probability');
    xline(0, 'k-');
    xline(-20, 'r--');
    xline(20, 'r--');
    xlabel('Signed SWS error (%)');
    ylabel('Fraction');
    title(models(i), 'Interpreter', 'none');
    grid on;
end

saveas(gcf, fullfile(fig_dir, 'level11_sws_error_histograms.png'));
close(gcf);

end

function plot_error_by_aperture(T_by_aperture, fig_dir)

T = T_by_aperture(T_by_aperture.model_type == "bagged_trees", :);
models = unique(T.model_name, 'stable');
figure('Color', 'w');
hold on;
for i = 1:numel(models)
    Ti = T(T.model_name == models(i), :);
    plot(Ti.Omega_sr, Ti.MAPE_pct, 'o-', ...
        'DisplayName', char(models(i)));
end
xlabel('\Omega (sr)');
ylabel('MAPE SWS (%)');
title('Test 11: error by aperture');
legend('Location', 'best', 'Interpreter', 'none');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_mape_by_aperture.png'));
close(gcf);

end

function plot_error_by_M(T_by_M, fig_dir)

T = T_by_M(T_by_M.model_type == "bagged_trees" | ...
    T_by_M.model_type == "global_q_oracle", :);
models = unique(T.model_name, 'stable');

figure('Color', 'w');
tiledlayout('flow');

for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    bar(categorical(Ti.REQ_M), Ti.MAPE_pct);
    xlabel('M');
    ylabel('MAPE SWS (%)');
    title(models(i), 'Interpreter', 'none');
    grid on;
end

saveas(gcf, fullfile(fig_dir, 'level11_mape_by_M.png'));
close(gcf);

end

function plot_error_by_wave(T_by_wave, fig_dir)

T = T_by_wave(T_by_wave.model_type == "bagged_trees" | ...
    T_by_wave.model_type == "global_q_oracle", :);
models = unique(T.model_name, 'stable');

figure('Color', 'w');
tiledlayout('flow');

for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    bar(categorical(Ti.SIM_WaveModel), Ti.MAPE_pct);
    xlabel('Wave model');
    ylabel('MAPE SWS (%)');
    title(models(i), 'Interpreter', 'none');
    grid on;
end

saveas(gcf, fullfile(fig_dir, 'level11_mape_by_wave_model.png'));
close(gcf);

end

function plot_error_vs_q_gap(T, fig_dir)

T = T((T.model_type == "bagged_trees" | ...
    T.model_type == "global_q_oracle") & ...
    isfinite(T.abs_q_local_minus_global), :);

figure('Color', 'w');
tiledlayout('flow');
models = unique(T.model_name, 'stable');

for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    scatter(Ti.abs_q_local_minus_global, Ti.abs_sws_error_pct, ...
        14, Ti.REQ_M, 'filled');
    yline(20, 'r--');
    xlabel('|q local - q global|');
    ylabel('|SWS error| (%)');
    title(models(i), 'Interpreter', 'none');
    grid on;
    cb = colorbar;
    cb.Label.String = 'M';
end

saveas(gcf, fullfile(fig_dir, 'level11_error_vs_q_local_global_gap.png'));
close(gcf);

end

function plot_high_error_rates(T_high, T_high_M, T_comp, fig_dir)

T = T_high(T_high.model_type == "bagged_trees" | ...
    T_high.model_type == "global_q_oracle", :);

figure('Color', 'w');
bar(categorical(T.model_name), T.HighError_gt20_pct);
ylabel('Rows with |SWS error| > 20 (%)');
title('Test 11: high-error rate by model');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_high_error_gt20_by_model.png'));
close(gcf);

T = T_high_M(T_high_M.model_type == "bagged_trees" | ...
    T_high_M.model_type == "global_q_oracle", :);
models = unique(T.model_name, 'stable');

figure('Color', 'w');
tiledlayout('flow');
for i = 1:numel(models)
    nexttile;
    Ti = T(T.model_name == models(i), :);
    bar(categorical(Ti.REQ_M), Ti.HighError_gt20_pct);
    xlabel('M');
    ylabel('Rows > 20% (%)');
    title(models(i), 'Interpreter', 'none');
    grid on;
end
saveas(gcf, fullfile(fig_dir, 'level11_high_error_gt20_by_M.png'));
close(gcf);

if ~isempty(T_comp)
    T = T_comp(T_comp.model_type == "bagged_trees" | ...
        T_comp.model_type == "global_q_oracle", :);
    models = unique(T.model_name, 'stable');

    figure('Color', 'w');
    tiledlayout('flow');
    for i = 1:numel(models)
        nexttile;
        Ti = T(T.model_name == models(i), :);
        bar(categorical(Ti.REQ_M), Ti.Pct_of_model_high_errors);
        xlabel('M');
        ylabel('Share of model outliers (%)');
        title(models(i), 'Interpreter', 'none');
        grid on;
    end
    saveas(gcf, fullfile(fig_dir, ...
        'level11_high_error_gt20_M_composition.png'));
    close(gcf);
end

end

function plot_local_vs_global_q(T, fig_dir)

figure('Color', 'w');
scatter(T.q_theory, T.q_global_theory, 12, T.REQ_M, 'filled');
hold on;
plot([0 1], [0 1], 'k--');
axis equal; xlim([0 1]); ylim([0 1]);
xlabel('q local');
ylabel('q global');
title('Test 11: local q vs global q');
cb = colorbar;
cb.Label.String = 'M';
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_q_local_vs_global.png'));
close(gcf);

end

function plot_q_gap_by_M_and_aperture(T, fig_dir)

figure('Color', 'w');
boxchart(categorical(T.REQ_M), abs(T.q_local_minus_global));
xlabel('M');
ylabel('|q local - q global|');
title('Test 11: q local/global gap by M');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_q_gap_by_M.png'));
close(gcf);

T.step_cat = categorical(T.step_idx);
figure('Color', 'w');
boxchart(T.step_cat, abs(T.q_local_minus_global));
xlabel('Aperture step');
ylabel('|q local - q global|');
title('Test 11: q local/global gap by aperture');
grid on;
saveas(gcf, fullfile(fig_dir, 'level11_q_gap_by_aperture.png'));
close(gcf);

end

function plot_q_gap_heatmap(T_gap, fig_dir)

models = unique(T_gap.SIM_WaveModel, 'stable');

figure('Color', 'w');
tiledlayout('flow');

for i = 1:numel(models)
    nexttile;
    Ti = T_gap(T_gap.SIM_WaveModel == models(i), :);
    Mvals = unique(Ti.REQ_M, 'stable');
    steps = unique(Ti.step_idx, 'stable');
    Z = nan(numel(Mvals), numel(steps));

    for m = 1:numel(Mvals)
        for s = 1:numel(steps)
            idx = Ti.REQ_M == Mvals(m) & Ti.step_idx == steps(s);
            if any(idx)
                Z(m, s) = mean(Ti.Mean(idx), 'omitnan');
            end
        end
    end

    imagesc(steps, Mvals, Z);
    set(gca, 'YDir', 'normal');
    xlabel('Aperture step');
    ylabel('M');
    title(models(i), 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'Mean |q local - q global|';
end

saveas(gcf, fullfile(fig_dir, 'level11_q_gap_heatmap.png'));
close(gcf);

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end
