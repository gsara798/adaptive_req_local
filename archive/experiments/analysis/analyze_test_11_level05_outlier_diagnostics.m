%% analyze_test_11_level05_outlier_diagnostics.m
% Test 11 Level 05: outlier diagnostics for grouped generalization.
%
% This script consumes Level 04 outputs. It does not regenerate Test 11 and
% does not modify the Level 04 analysis.

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

%% Load Test 11 and Level 04 predictions

[T_feat, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    'test_11_global_req_features', ...
    'RootDir', root_dir, ...
    'Verbose', true);

level04_dir = fullfile(PATHS.analysis_dir, ...
    'level_04_grouped_generalization');
level04_pred_file = fullfile(level04_dir, 'tables', ...
    'level11_level04_grouped_predictions.csv');

if ~exist(level04_pred_file, 'file')
    error(['Level 04 predictions not found:\n%s\n', ...
        'Run experiments/analysis/analyze_test_11_level04_grouped_generalization.m first.'], ...
        level04_pred_file);
end

fprintf('\nLoading Level 04 grouped predictions:\n%s\n', level04_pred_file);
T_pred = readtable(level04_pred_file, 'TextType', 'string');

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_05_outlier_diagnostics');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
model_dir = fullfile(analysis_dir, 'models');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(model_dir);

%% Prepare joined diagnostic table

T_feat.row_key = ensure_row_key(T_feat);
T_pred.row_key = ensure_row_key(T_pred);

required_pred = [
    "q_true"
    "q_pred"
    "generalization_test"
    "heldout_var"
    "heldout_value"
    "model_name"
    "model_type"
    "model_role"
    "row_key"];
assert(all(ismember(required_pred, string(T_pred.Properties.VariableNames))), ...
    'Level 04 prediction table is missing required columns.');

operational_models = ["LocalOnly", "GlobalOnly", "HybridLocalGlobal"];
T = T_pred(ismember(T_pred.model_name, operational_models) & ...
    T_pred.model_role == "operational", :);

T = attach_reference_columns(T, T_feat);
T = ensure_sws_columns(T, T_feat);
T = add_effective_window_diagnostics(T);
T = add_outlier_flags(T);
T = add_diagnostic_bins(T);

assert(~isempty(T), 'Joined Level 05 table is empty.');
assert(all(ismember(operational_models, unique(T.model_name))), ...
    'Not all operational models are present in Level 05 table.');

%% Summaries and correlations

T_worst = sortrows(T, 'abs_sws_error_pct', 'descend');
T_worst = T_worst(1:min(200, height(T_worst)), :);

T_out_by_model = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"]);
T_out_by_M = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "REQ_M"]);
T_out_by_Meff_true = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "M_eff_true_diag_bin"]);
T_out_by_Meff_guess = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "M_eff_guess_bin"]);
T_out_by_frequency = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "SIM_f0"]);
T_out_by_aperture = summarize_outlier_rates(T, aperture_group_vars(T));
T_out_by_wave = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "SIM_WaveModel"]);
T_out_by_q = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "model_type"
    "q_true_bin"]);
T_out_by_condition = summarize_outlier_rates(T, [
    "generalization_test"
    "heldout_value"
    "model_name"
    "model_type"
    "SIM_f0"
    "REQ_M"
    "SIM_WaveModel"
    "step_idx"]);

if ismember('abs_q_local_minus_global', string(T.Properties.VariableNames))
    T_out_by_q_gap = summarize_outlier_rates(T, [
        "generalization_test"
        "model_name"
        "model_type"
        "q_gap_bin"]);
else
    T_out_by_q_gap = table();
end

T_corr = compute_error_feature_correlations(T);

T_test = T(T.split == "test", :);
T_support_by_M = summarize_support_occupancy(T_test, "REQ_M");
T_support_by_Meff_guess = summarize_support_occupancy(T_test, ...
    "M_eff_guess_bin");
T_support_by_Meff_true = summarize_support_occupancy(T_test, ...
    "M_eff_true_diag_bin");
T_support_by_frequency = summarize_support_occupancy(T_test, "SIM_f0");
T_support_by_aperture = summarize_support_occupancy(T_test, "step_idx");
T_support_by_q = summarize_support_occupancy(T_test, "q_true_bin");

%% Save tables

writetable(remove_cell_columns(T), fullfile(table_dir, ...
    'level11_level05_joined_predictions_with_diagnostics.csv'));
writetable(select_existing_columns(remove_cell_columns(T_worst), [
    "generalization_test"
    "heldout_var"
    "heldout_value"
    "model_name"
    "model_type"
    "SIM_f0"
    "REQ_M"
    "M_eff_guess"
    "M_eff_true_diag"
    "SIM_cs_bg"
    "SIM_WaveModel"
    "step_idx"
    "Omega_sr"
    "q_true"
    "q_pred"
    "residual"
    "abs_error"
    "cs_true"
    "cs_pred"
    "sws_error_pct"
    "abs_sws_error_pct"
    "q_local_minus_global"
    "abs_q_local_minus_global"]), ...
    fullfile(table_dir, 'level11_level05_worst_200_predictions.csv'));
writetable(T_out_by_model, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_model.csv'));
writetable(T_out_by_M, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_M.csv'));
writetable(T_out_by_Meff_true, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_M_eff_true_diag.csv'));
writetable(T_out_by_Meff_guess, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_M_eff_guess.csv'));
writetable(T_out_by_frequency, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_frequency.csv'));
writetable(T_out_by_aperture, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_aperture.csv'));
writetable(T_out_by_wave, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_wave_model.csv'));
writetable(T_out_by_q, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_q_true_bin.csv'));
writetable(T_out_by_condition, fullfile(table_dir, ...
    'level11_level05_outlier_rate_by_condition_extra.csv'));
if ~isempty(T_out_by_q_gap)
    writetable(T_out_by_q_gap, fullfile(table_dir, ...
        'level11_level05_outlier_rate_by_q_local_global_gap.csv'));
end
writetable(T_corr, fullfile(table_dir, ...
    'level11_level05_error_feature_correlations.csv'));
writetable(T_support_by_M, fullfile(table_dir, ...
    'level11_level05_support_by_M.csv'));
writetable(T_support_by_Meff_guess, fullfile(table_dir, ...
    'level11_level05_support_by_M_eff_guess.csv'));
writetable(T_support_by_Meff_true, fullfile(table_dir, ...
    'level11_level05_support_by_M_eff_true_diag.csv'));
writetable(T_support_by_frequency, fullfile(table_dir, ...
    'level11_level05_support_by_frequency.csv'));
writetable(T_support_by_aperture, fullfile(table_dir, ...
    'level11_level05_support_by_aperture.csv'));
writetable(T_support_by_q, fullfile(table_dir, ...
    'level11_level05_support_by_q_true_bin.csv'));

save(fullfile(model_dir, 'level11_level05_outlier_diagnostics.mat'), ...
    'T_out_by_model', 'T_out_by_M', 'T_out_by_Meff_true', ...
    'T_out_by_Meff_guess', 'T_out_by_frequency', 'T_out_by_aperture', ...
    'T_out_by_wave', 'T_out_by_q', 'T_out_by_q_gap', ...
    'T_out_by_condition', 'T_corr', ...
    'T_support_by_M', 'T_support_by_Meff_guess', ...
    'T_support_by_Meff_true', 'T_support_by_frequency', ...
    'T_support_by_aperture', 'T_support_by_q', ...
    'MC', 'PATHS', '-v7.3');

%% Figures

T_fig = T(T.model_type == "bagged_trees", :);

plot_error_vs_numeric(T_fig, 'M_eff_true_diag', ...
    'M_{eff,true,diag}', fig_dir, ...
    'level11_level05_error_vs_M_eff_true_diag');
plot_error_vs_numeric(T_fig, 'M_eff_guess', ...
    'M_{eff,guess}', fig_dir, ...
    'level11_level05_error_vs_M_eff_guess');
plot_high_error_by_bin(T_out_by_Meff_true, ...
    'M_eff_true_diag_bin', fig_dir, ...
    'level11_level05_high_error_rate_by_M_eff_true_diag_bin');
plot_high_error_by_bin(T_out_by_Meff_guess, ...
    'M_eff_guess_bin', fig_dir, ...
    'level11_level05_high_error_rate_by_M_eff_guess_bin');
plot_heatmap_nominal_M_vs_Meff(T_fig, fig_dir);
plot_error_Meff_guess_vs_true(T_fig, fig_dir);
plot_error_vs_numeric(T_fig, 'q_true', 'q true', fig_dir, ...
    'level11_level05_error_vs_q_true');
plot_error_vs_numeric(T_fig, 'abs_error', '|q pred - q true|', fig_dir, ...
    'level11_level05_sws_error_vs_q_error');
if ismember('abs_q_local_minus_global', string(T_fig.Properties.VariableNames))
    plot_error_vs_numeric(T_fig, 'abs_q_local_minus_global', ...
        '|q local - q global|', fig_dir, ...
        'level11_level05_error_vs_q_local_global_gap');
end
plot_outlier_composition(T_fig, 'REQ_M', fig_dir, ...
    'level11_level05_outlier_composition_by_M');
plot_outlier_composition(T_fig, 'M_eff_true_diag_bin', fig_dir, ...
    'level11_level05_outlier_composition_by_M_eff_true_diag');
plot_outlier_composition(T_fig, 'M_eff_guess_bin', fig_dir, ...
    'level11_level05_outlier_composition_by_M_eff_guess');
plot_top_feature_correlations(T_corr, fig_dir);
plot_high_error_heatmap_frequency_M(T_fig, fig_dir);
plot_high_error_heatmap_Meff_aperture(T_fig, fig_dir);
plot_error_by_wave_and_generalization(T_fig, fig_dir);
plot_q_error_by_Meff_bins(T_fig, fig_dir);
plot_support_vs_outlier_rate(T_test, "REQ_M", fig_dir, ...
    'level11_level05_support_vs_outlier_rate_by_M');
plot_support_vs_outlier_rate(T_test, "M_eff_true_diag_bin", fig_dir, ...
    'level11_level05_support_vs_outlier_rate_by_M_eff_true_diag');
plot_support_vs_outlier_rate(T_test, "M_eff_guess_bin", fig_dir, ...
    'level11_level05_support_vs_outlier_rate_by_M_eff_guess');
plot_outlier_support_summary(T_test, fig_dir);

%% Console report

fprintf('\nTest 11 Level 05 outlier diagnostics complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);

fprintf('\nTop 10 bagged-trees conditions with highest high-error rate.\n');
T_top = T_out_by_condition(T_out_by_condition.model_type == "bagged_trees", :);
T_top = sortrows(T_top, 'HighError_gt20_pct', 'descend');
disp(T_top(1:min(10, height(T_top)), :));

fprintf('\nPer-model outlier summary for bagged trees.\n');
T_model_report = build_model_report(T_fig, ...
    T_out_by_model, T_out_by_M, T_out_by_Meff_guess, ...
    T_out_by_Meff_true, T_out_by_frequency, T_out_by_aperture, T_out_by_q);
disp(T_model_report);

print_auto_conclusions(T_model_report);

fprintf('\nSupport/occupancy summary for bagged-trees test rows.\n');
T_support_report = build_support_console_report(T_test);
disp(T_support_report);
print_support_conclusions(T_support_report);

fprintf('\nM_eff_true_diag is diagnostic only because it uses the true simulated shear wave speed.\n');
fprintf('M_eff_guess is closer to an operational/nominal effective-window variable.\n');

%% Local functions

function key = ensure_row_key(T)

if ismember('row_key', string(T.Properties.VariableNames)) && ...
        all(strlength(string(T.row_key)) > 0)
    key = string(T.row_key);
    return;
end

required = ["condition_id", "step_idx", "realization_idx", "patch_idx"];
assert(all(ismember(required, string(T.Properties.VariableNames))), ...
    'Cannot create row_key because required id columns are missing.');

parts = strings(height(T), 4);
parts(:, 1) = string(T.condition_id);
parts(:, 2) = string(T.step_idx);
parts(:, 3) = string(T.realization_idx);
parts(:, 4) = string(T.patch_idx);
key = join(parts, "|", 2);

end

function T = attach_reference_columns(T, T_ref)

ref_cols = [
    "row_key"
    "REQ_M"
    "M_eff_guess"
    "M_eff_true_diag"
    "REQ_cs_guess"
    "SIM_f0"
    "SIM_cs_bg"
    "SIM_WaveModel"
    "step_idx"
    "Omega_sr"
    "omega_sr"
    "q_theory"
    "q_local_minus_global"
    "req_mapping"];

ref_cols = ref_cols(ismember(ref_cols, string(T_ref.Properties.VariableNames)));
T_ref_small = T_ref(:, cellstr(ref_cols));

[tf, loc] = ismember(string(T.row_key), string(T_ref_small.row_key));
assert(all(tf), 'Some Level 04 prediction rows do not match T_feat row_key.');

for i = 1:numel(ref_cols)
    name = ref_cols(i);
    if name == "row_key"
        continue;
    end

    if ismember(name, string(T.Properties.VariableNames))
        continue;
    end

    T.(char(name)) = T_ref_small.(char(name))(loc);
end

if ~ismember('Omega_sr', string(T.Properties.VariableNames)) && ...
        ismember('omega_sr', string(T.Properties.VariableNames))
    T.Omega_sr = T.omega_sr;
end

if ~ismember('q_theory', string(T.Properties.VariableNames)) && ...
        ismember('q_true', string(T.Properties.VariableNames))
    T.q_theory = T.q_true;
end

end

function T = ensure_sws_columns(T, T_ref)

required_sws = ["cs_true", "cs_pred", "sws_error_pct", ...
    "abs_sws_error_pct"];
if all(ismember(required_sws, string(T.Properties.VariableNames)))
    return;
end

[tf, loc] = ismember(string(T.row_key), string(T_ref.row_key));
assert(all(tf), 'Cannot reconstruct SWS because row_key matching failed.');

if ~ismember('cs_true', string(T.Properties.VariableNames))
    T.cs_true = T_ref.SIM_cs_bg(loc);
end

if ~ismember('cs_pred', string(T.Properties.VariableNames))
    cs_pred = nan(height(T), 1);
    for i = 1:height(T)
        mapping_i = T_ref.req_mapping{loc(i)};
        if isempty(mapping_i) || ~isfinite(T.q_pred(i))
            continue;
        end
        k_i = adaptive_req.quantile.quantile_to_k(mapping_i, T.q_pred(i));
        cs_pred(i) = 2*pi*T_ref.SIM_f0(loc(i)) ./ k_i;
    end
    T.cs_pred = cs_pred;
end

if ~ismember('sws_error', string(T.Properties.VariableNames))
    T.sws_error = T.cs_pred - T.cs_true;
end
if ~ismember('abs_sws_error', string(T.Properties.VariableNames))
    T.abs_sws_error = abs(T.sws_error);
end
if ~ismember('sws_error_pct', string(T.Properties.VariableNames))
    T.sws_error_pct = 100 * T.sws_error ./ T.cs_true;
end
if ~ismember('abs_sws_error_pct', string(T.Properties.VariableNames))
    T.abs_sws_error_pct = abs(T.sws_error_pct);
end

end

function T = add_effective_window_diagnostics(T)

has_M_eff_true_diag = ismember('M_eff_true_diag', T.Properties.VariableNames);
has_M_eff_guess = ismember('M_eff_guess', T.Properties.VariableNames);

if ~has_M_eff_true_diag
    if ismember('REQ_cs_guess', T.Properties.VariableNames)
        cs_guess_i = T.REQ_cs_guess;
    else
        cs_guess_i = 3.0 * ones(height(T), 1);
    end
    T.M_eff_true_diag = T.REQ_M .* cs_guess_i ./ T.SIM_cs_bg;
end

if ~has_M_eff_guess
    T.M_eff_guess = T.REQ_M;
end

if ismember('REQ_cs_guess', T.Properties.VariableNames)
    T.cs_guess_over_true = T.REQ_cs_guess ./ T.SIM_cs_bg;
    T.k_guess_over_true = T.SIM_cs_bg ./ T.REQ_cs_guess;
end

if ismember('q_local_minus_global', T.Properties.VariableNames)
    T.abs_q_local_minus_global = abs(T.q_local_minus_global);
end

end

function T = add_outlier_flags(T)

T.is_high_error_20 = T.abs_sws_error_pct > 20;
T.is_high_error_10 = T.abs_sws_error_pct > 10;
T.is_top5_error = false(height(T), 1);
T.is_top1_error = false(height(T), 1);

[G, ~] = findgroups(T(:, {'model_name', 'model_type', ...
    'generalization_test'}));

for g = 1:max(G)
    idx = find(G == g);
    x = T.abs_sws_error_pct(idx);
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    p95 = prctile(x, 95);
    p99 = prctile(x, 99);
    T.is_top5_error(idx) = T.abs_sws_error_pct(idx) >= p95;
    T.is_top1_error(idx) = T.abs_sws_error_pct(idx) >= p99;
end

end

function T = add_diagnostic_bins(T)

eff_edges = [0 1.5 2 2.5 3 3.5 4 4.5 5 6 Inf];
q_edges = [0 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0];

T.M_eff_true_diag_bin = discretize_to_labels(T.M_eff_true_diag, eff_edges);
T.M_eff_guess_bin = discretize_to_labels(T.M_eff_guess, eff_edges);
T.q_true_bin = discretize_to_labels(T.q_true, q_edges);
T.aperture_bin = categorical(string(T.step_idx));

if ismember('abs_q_local_minus_global', string(T.Properties.VariableNames))
    finite_gap = T.abs_q_local_minus_global(isfinite(T.abs_q_local_minus_global));
    if isempty(finite_gap) || max(finite_gap) <= 0
        gap_edges = [0 Inf];
    else
        gap_edges = unique([0, prctile(finite_gap, [50 75 90 95]), Inf]);
    end
    T.q_gap_bin = discretize_to_labels(T.abs_q_local_minus_global, gap_edges);
end

end

function bins = discretize_to_labels(x, edges)

x = double(x);
labels = strings(numel(edges) - 1, 1);
for i = 1:numel(labels)
    if isinf(edges(i + 1))
        labels(i) = sprintf('[%.2g, Inf)', edges(i));
    else
        labels(i) = sprintf('[%.2g, %.2g)', edges(i), edges(i + 1));
    end
end

idx = discretize(x, edges);
bins = strings(numel(x), 1);
valid = isfinite(idx);
bins(valid) = labels(idx(valid));
bins(~valid) = "missing";
bins = categorical(bins, [labels; "missing"], 'Ordinal', true);

end

function group_vars = aperture_group_vars(T)

group_vars = [
    "generalization_test"
    "model_name"
    "model_type"
    "step_idx"];
if ismember('Omega_sr', string(T.Properties.VariableNames))
    group_vars = [group_vars; "Omega_sr"];
end

end

function T_sum = summarize_outlier_rates(T, group_vars)

group_vars = string(group_vars);
[G, T_keys] = findgroups(T(:, cellstr(group_vars)));

T_sum = T_keys;
T_sum.N = splitapply(@numel, T.abs_sws_error_pct, G);
T_sum.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.Median_abs_error_pct = splitapply(@(x) median(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.P95_abs_error_pct = splitapply(@(x) prctile(x, 95), ...
    T.abs_sws_error_pct, G);
T_sum.HighError_gt10_pct = splitapply(@(x) 100 * mean(x > 10, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.HighError_gt20_pct = splitapply(@(x) 100 * mean(x > 20, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.Top5Error_pct = splitapply(@(x) 100 * mean(x, 'omitnan'), ...
    T.is_top5_error, G);
T_sum.Top1Error_pct = splitapply(@(x) 100 * mean(x, 'omitnan'), ...
    T.is_top1_error, G);

T_sum = sortrows(T_sum, 'HighError_gt20_pct', 'descend');

end

function T_support = summarize_support_occupancy(T, bin_var, varargin)

p = inputParser;
p.FunctionName = 'summarize_support_occupancy';
addRequired(p, 'T', @istable);
addRequired(p, 'bin_var', @(x) ischar(x) || isstring(x));
addParameter(p, 'BaseVars', ["generalization_test", "model_name", ...
    "model_type"], @(x) isstring(x) || iscellstr(x));
parse(p, T, bin_var, varargin{:});

bin_var = string(bin_var);
base_vars = string(p.Results.BaseVars);

if ~ismember(bin_var, string(T.Properties.VariableNames))
    T_support = table();
    return;
end

T = T(T.split == "test", :);
T = T(isfinite(T.abs_sws_error_pct), :);
if isempty(T)
    T_support = table();
    return;
end

[Gbase, T_base] = findgroups(T(:, cellstr(base_vars)));
rows = struct([]);
row_idx = 0;

for g = 1:max(Gbase)
    idx_base = Gbase == g;
    Tb = T(idx_base, :);
    total_rows = height(Tb);
    total_outliers = sum(Tb.is_high_error_20);
    bin_values = unique(Tb.(char(bin_var)), 'stable');

    for b = 1:numel(bin_values)
        value_b = bin_values(b);
        idx_bin = is_same_group_value(Tb.(char(bin_var)), value_b);
        n_bin = sum(idx_bin);
        n_out = sum(Tb.is_high_error_20(idx_bin));

        row_idx = row_idx + 1;
        for k = 1:numel(base_vars)
            rows(row_idx).(char(base_vars(k))) = ...
                T_base.(char(base_vars(k)))(g);
        end
        rows(row_idx).support_var = bin_var;
        rows(row_idx).(char(bin_var)) = value_b;
        rows(row_idx).N_total_bin = n_bin;
        rows(row_idx).N_outliers_gt20 = n_out;
        rows(row_idx).Pct_of_all_rows = 100 * n_bin / total_rows;
        rows(row_idx).HighError_gt20_pct = 100 * n_out / max(n_bin, 1);
        if total_outliers > 0
            rows(row_idx).Pct_of_all_outliers = 100 * n_out / total_outliers;
        else
            rows(row_idx).Pct_of_all_outliers = 0;
        end
        rows(row_idx).Total_rows_model_test = total_rows;
        rows(row_idx).Total_outliers_model_test = total_outliers;
    end
end

if isempty(rows)
    T_support = table();
else
    T_support = struct2table(rows);
    sort_vars = [cellstr(base_vars(:).'), {'HighError_gt20_pct'}];
    sort_dirs = [repmat({'ascend'}, 1, numel(base_vars)), {'descend'}];
    T_support = sortrows(T_support, sort_vars, sort_dirs);
end

end

function tf = is_same_group_value(x, value)

if iscategorical(x) || isstring(x) || iscellstr(x) || ischar(x)
    tf = string(x) == string(value);
else
    tf = x == value;
end
tf = tf(:);

end

function T_corr = compute_error_feature_correlations(T)

numeric_vars = string(T.Properties.VariableNames(varfun(@isnumeric, T, ...
    'OutputFormat', 'uniform')));

exclude = lower([
    "q_pred"
    "q_pred_raw"
    "cs_pred"
    "sws_error"
    "abs_sws_error"
    "sws_error_pct"
    "abs_sws_error_pct"
    "residual"
    "abs_error"
    "is_high_error_10"
    "is_high_error_20"
    "is_top5_error"
    "is_top1_error"
    "N_train"
    "N_test"
    "condition_id"
    "condition_position"
    "realization_idx"
    "patch_idx"]);

numeric_vars = numeric_vars(~ismember(lower(numeric_vars), exclude));

[G, keys] = findgroups(T(:, {'model_name', 'model_type', ...
    'generalization_test'}));
rows = struct([]);
row_idx = 0;

for g = 1:max(G)
    idx = G == g;
    y = T.abs_sws_error_pct(idx);
    for j = 1:numel(numeric_vars)
        var_j = numeric_vars(j);
        x = T.(char(var_j))(idx);
        valid = isfinite(x) & isfinite(y);
        if sum(valid) < 10 || std(x(valid)) <= eps || std(y(valid)) <= eps
            continue;
        end
        row_idx = row_idx + 1;
        rows(row_idx).model_name = keys.model_name(g);
        rows(row_idx).model_type = keys.model_type(g);
        rows(row_idx).generalization_test = keys.generalization_test(g);
        rows(row_idx).feature_name = var_j;
        rows(row_idx).N = sum(valid);
        rows(row_idx).Spearman = corr(x(valid), y(valid), ...
            'Type', 'Spearman', 'Rows', 'complete');
        rows(row_idx).Pearson = corr(x(valid), y(valid), ...
            'Type', 'Pearson', 'Rows', 'complete');
        rows(row_idx).AbsSpearman = abs(rows(row_idx).Spearman);
    end
end

if isempty(rows)
    T_corr = table();
else
    T_corr = struct2table(rows);
    T_corr = sortrows(T_corr, 'AbsSpearman', 'descend');
end

end

function plot_error_vs_numeric(T, var_name, x_label, fig_dir, out_name)

if ~ismember(var_name, string(T.Properties.VariableNames))
    return;
end

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 390]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    scatter(ax, Ti.(var_name), Ti.abs_sws_error_pct, 8, ...
        double(Ti.REQ_M), 'filled', 'MarkerFaceAlpha', 0.25);
    yline(ax, 20, 'r--', '20%');
    xlabel(ax, x_label, 'Interpreter', 'tex');
    ylabel(ax, '|SWS error| (%)');
    title(ax, models(i), 'Interpreter', 'none');
    grid(ax, 'on');
end

export_figure(fig, fig_dir, out_name);

end

function plot_high_error_by_bin(T_sum, bin_var, fig_dir, out_name)

T = T_sum(T_sum.model_type == "bagged_trees", :);
if isempty(T) || ~ismember(bin_var, string(T.Properties.VariableNames))
    return;
end

plot_grouped_metric(T, bin_var, 'HighError_gt20_pct', ...
    '|SWS error| > 20 (%)');
export_figure(gcf, fig_dir, out_name);

end

function plot_grouped_metric(T, x_var, y_var, y_label)

models = unique(T.model_name, 'stable');
x_vals = unique(string(T.(char(x_var))), 'stable');
Y = nan(numel(x_vals), numel(models));

for i = 1:numel(x_vals)
    for j = 1:numel(models)
        idx = string(T.(char(x_var))) == x_vals(i) & ...
            T.model_name == models(j);
        if any(idx)
            Y(i, j) = mean(T.(char(y_var))(idx), 'omitnan');
        end
    end
end

figure('Color', 'w', 'Position', [100 100 1050 460]);
bar(categorical(x_vals, x_vals), Y);
ylabel(y_label);
xlabel(strrep(string(x_var), '_', '\_'));
legend(models, 'Location', 'best', 'Interpreter', 'none');
xtickangle(25);
grid on;

end

function plot_heatmap_nominal_M_vs_Meff(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 390]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    M_vals = unique(Ti.REQ_M, 'stable');
    bins = categories(Ti.M_eff_true_diag_bin);
    Z = nan(numel(bins), numel(M_vals));
    for r = 1:numel(bins)
        for c = 1:numel(M_vals)
            idx = Ti.REQ_M == M_vals(c) & ...
                string(Ti.M_eff_true_diag_bin) == string(bins{r});
            if any(idx)
                Z(r, c) = mean(Ti.abs_sws_error_pct(idx), 'omitnan');
            end
        end
    end
    imagesc(ax, M_vals, 1:numel(bins), Z);
    set(ax, 'YDir', 'normal', 'YTick', 1:numel(bins), ...
        'YTickLabel', bins);
    xlabel(ax, 'REQ M');
    ylabel(ax, 'M_{eff,true,diag} bin', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = 'MAPE (%)';
end

export_figure(fig, fig_dir, ...
    'level11_level05_error_by_nominal_M_and_M_eff_true_diag');

end

function plot_error_Meff_guess_vs_true(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 390]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    scatter(ax, Ti.M_eff_guess, Ti.M_eff_true_diag, 10, ...
        Ti.abs_sws_error_pct, 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel(ax, 'M_{eff,guess}', 'Interpreter', 'tex');
    ylabel(ax, 'M_{eff,true,diag}', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    grid(ax, 'on');
    cb = colorbar(ax);
    cb.Label.String = '|SWS error| (%)';
end

export_figure(fig, fig_dir, ...
    'level11_level05_error_M_eff_guess_vs_true_diag');

end

function plot_outlier_composition(T, group_var, fig_dir, out_name)

if ~ismember(group_var, string(T.Properties.VariableNames))
    return;
end

T_high = T(T.is_high_error_20, :);
if isempty(T_high)
    return;
end

models = unique(T.model_name, 'stable');
x_vals = unique(string(T.(char(group_var))), 'stable');
Y = nan(numel(x_vals), numel(models));

for i = 1:numel(x_vals)
    for j = 1:numel(models)
        denom = sum(T_high.model_name == models(j));
        idx = T_high.model_name == models(j) & ...
            string(T_high.(char(group_var))) == x_vals(i);
        if denom > 0
            Y(i, j) = 100 * sum(idx) / denom;
        end
    end
end

figure('Color', 'w', 'Position', [100 100 1050 460]);
bar(categorical(x_vals, x_vals), Y);
ylabel('Share of >20% outliers (%)');
xlabel(strrep(string(group_var), '_', '\_'));
legend(models, 'Location', 'best', 'Interpreter', 'none');
xtickangle(25);
grid on;

export_figure(gcf, fig_dir, out_name);

end

function plot_top_feature_correlations(T_corr, fig_dir)

T = T_corr(T_corr.model_type == "bagged_trees", :);
if isempty(T)
    return;
end

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1250 430]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    Ti = sortrows(Ti, 'AbsSpearman', 'descend');
    [~, unique_idx] = unique(Ti.feature_name, 'stable');
    Ti = Ti(unique_idx, :);
    Ti = Ti(1:min(15, height(Ti)), :);
    feature_order = flip(Ti.feature_name);
    barh(ax, categorical(feature_order, feature_order), ...
        flip(Ti.Spearman));
    xlabel(ax, 'Spearman rho');
    title(ax, models(i), 'Interpreter', 'none');
    grid(ax, 'on');
end

export_figure(fig, fig_dir, ...
    'level11_level05_top_error_feature_correlations');

end

function plot_high_error_heatmap_frequency_M(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 390]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    freqs = unique(Ti.SIM_f0, 'stable');
    M_vals = unique(Ti.REQ_M, 'stable');
    Z = nan(numel(M_vals), numel(freqs));
    for r = 1:numel(M_vals)
        for c = 1:numel(freqs)
            idx = Ti.REQ_M == M_vals(r) & Ti.SIM_f0 == freqs(c);
            if any(idx)
                Z(r, c) = 100 * mean(Ti.is_high_error_20(idx), 'omitnan');
            end
        end
    end
    imagesc(ax, freqs, M_vals, Z);
    set(ax, 'YDir', 'normal');
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'REQ M');
    title(ax, models(i), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = '>20% error rate (%)';
end

export_figure(fig, fig_dir, ...
    'level11_level05_high_error_heatmap_frequency_M');

end

function plot_high_error_heatmap_Meff_aperture(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 410]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    steps = unique(Ti.step_idx, 'stable');
    bins = categories(Ti.M_eff_true_diag_bin);
    Z = nan(numel(bins), numel(steps));
    for r = 1:numel(bins)
        for c = 1:numel(steps)
            idx = string(Ti.M_eff_true_diag_bin) == string(bins{r}) & ...
                Ti.step_idx == steps(c);
            if any(idx)
                Z(r, c) = 100 * mean(Ti.is_high_error_20(idx), 'omitnan');
            end
        end
    end
    imagesc(ax, steps, 1:numel(bins), Z);
    set(ax, 'YDir', 'normal', 'YTick', 1:numel(bins), ...
        'YTickLabel', bins);
    xlabel(ax, 'step idx');
    ylabel(ax, 'M_{eff,true,diag} bin', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = '>20% error rate (%)';
end

export_figure(fig, fig_dir, ...
    'level11_level05_high_error_heatmap_Meff_aperture');

end

function plot_error_by_wave_and_generalization(T, fig_dir)

if ~ismember('SIM_WaveModel', string(T.Properties.VariableNames))
    return;
end

Tsum = summarize_outlier_rates(T, [
    "generalization_test"
    "model_name"
    "SIM_WaveModel"]);
plot_grouped_metric(Tsum, 'SIM_WaveModel', 'HighError_gt20_pct', ...
    '|SWS error| > 20 (%)');
title('High-error rate by wave model');
export_figure(gcf, fig_dir, ...
    'level11_level05_high_error_by_wave_model_extra');

end

function plot_q_error_by_Meff_bins(T, fig_dir)

models = unique(T.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 390]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');
for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T(T.model_name == models(i), :);
    boxchart(ax, Ti.M_eff_true_diag_bin, Ti.abs_error);
    ylabel(ax, '|q error|');
    xlabel(ax, 'M_{eff,true,diag} bin', 'Interpreter', 'tex');
    title(ax, models(i), 'Interpreter', 'none');
    xtickangle(ax, 25);
    grid(ax, 'on');
end
export_figure(fig, fig_dir, ...
    'level11_level05_q_error_by_M_eff_true_diag_bin_extra');

end

function plot_support_vs_outlier_rate(T, bin_var, fig_dir, out_name)

T_plot = T(T.model_type == "bagged_trees" & T.split == "test", :);
if isempty(T_plot) || ~ismember(bin_var, string(T_plot.Properties.VariableNames))
    return;
end

T_support = summarize_support_occupancy(T_plot, bin_var, ...
    'BaseVars', ["model_name", "model_type"]);
if isempty(T_support)
    return;
end

models = unique(T_support.model_name, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1250 420]);
tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T_support(T_support.model_name == models(i), :);
    x_labels = string(Ti.(char(bin_var)));
    x_cat = categorical(x_labels, x_labels);

    yyaxis(ax, 'left');
    bar(ax, x_cat, Ti.Pct_of_all_rows, 0.72, ...
        'FaceColor', [0.75 0.75 0.75], ...
        'DisplayName', '% all rows');
    ylabel(ax, 'Rows in bin (%)');

    yyaxis(ax, 'right');
    hold(ax, 'on');
    plot(ax, x_cat, Ti.HighError_gt20_pct, 'o-', ...
        'Color', [0.8500 0.3250 0.0980], ...
        'LineWidth', 2.0, ...
        'DisplayName', '>20% error rate');
    plot(ax, x_cat, Ti.Pct_of_all_outliers, 's--', ...
        'Color', [0 0.4470 0.7410], ...
        'LineWidth', 1.6, ...
        'DisplayName', '% all outliers');
    ylabel(ax, 'Outlier metrics (%)');
    title(ax, models(i), 'Interpreter', 'none');
    xlabel(ax, strrep(string(bin_var), '_', '\_'));
    xtickangle(ax, 25);
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end

export_figure(fig, fig_dir, out_name);

end

function plot_outlier_support_summary(T, fig_dir)

T_plot = T(T.model_type == "bagged_trees" & T.split == "test", :);
models = unique(T_plot.model_name, 'stable');
support_vars = ["REQ_M", "M_eff_guess_bin", "M_eff_true_diag_bin"];
titles = ["REQ M", "M_{eff,guess}", "M_{eff,true,diag}"];

fig = figure('Color', 'w', 'Position', [100 100 1250 900]);
tl = tiledlayout(fig, numel(models), numel(support_vars), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(models)
    for j = 1:numel(support_vars)
        ax = nexttile(tl);
        bin_var = support_vars(j);
        Ti0 = T_plot(T_plot.model_name == models(i), :);
        Ti = summarize_support_occupancy(Ti0, bin_var, ...
            'BaseVars', ["model_name", "model_type"]);
        if isempty(Ti)
            continue;
        end
        x_labels = string(Ti.(char(bin_var)));
        x_cat = categorical(x_labels, x_labels);

        yyaxis(ax, 'left');
        bar(ax, x_cat, Ti.Pct_of_all_rows, ...
            'FaceColor', [0.72 0.72 0.72]);
        ylabel(ax, 'Rows (%)');

        yyaxis(ax, 'right');
        plot(ax, x_cat, Ti.HighError_gt20_pct, 'o-', ...
            'LineWidth', 1.8, ...
            'Color', [0.8500 0.3250 0.0980]);
        ylabel(ax, '>20% error (%)');

        if i == 1
            title(ax, titles(j), 'Interpreter', 'tex');
        end
        if j == 1
            text(ax, -0.18, 0.5, models(i), ...
                'Units', 'normalized', ...
                'Rotation', 90, ...
                'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end
        xtickangle(ax, 25);
        grid(ax, 'on');
    end
end

export_figure(fig, fig_dir, ...
    'level11_level05_outlier_support_summary');

end

function T_report = build_model_report(T_fig, ~, T_by_M, ...
    T_by_Mguess, T_by_Mtrue, T_by_freq, T_by_aperture, T_by_q)

models = unique(T_fig.model_name, 'stable');
rows = struct([]);

for i = 1:numel(models)
    model = models(i);
    Ti = T_fig(T_fig.model_name == model, :);
    rows(i).model_name = model;
    rows(i).MAPE_pct = mean(Ti.abs_sws_error_pct, 'omitnan');
    rows(i).HighError_gt20_pct = ...
        100 * mean(Ti.abs_sws_error_pct > 20, 'omitnan');
    rows(i).worst_REQ_M = string(pick_worst(T_by_M, model, "REQ_M").REQ_M);
    rows(i).worst_M_eff_guess_bin = string(pick_worst( ...
        T_by_Mguess, model, "M_eff_guess_bin").M_eff_guess_bin);
    rows(i).worst_M_eff_true_diag_bin = string(pick_worst( ...
        T_by_Mtrue, model, "M_eff_true_diag_bin").M_eff_true_diag_bin);
    rows(i).worst_frequency = string(pick_worst( ...
        T_by_freq, model, "SIM_f0").SIM_f0);
    worst_ap = pick_worst(T_by_aperture, model, "step_idx");
    rows(i).worst_aperture = string(worst_ap.step_idx);
    rows(i).worst_q_true_bin = string(pick_worst( ...
        T_by_q, model, "q_true_bin").q_true_bin);
end

T_report = struct2table(rows);

end

function T_support_report = build_support_console_report(T_test)

T_plot = T_test(T_test.model_type == "bagged_trees" & ...
    T_test.split == "test", :);
models = unique(T_plot.model_name, 'stable');
rows = struct([]);

for i = 1:numel(models)
    model = models(i);
    Ti = T_plot(T_plot.model_name == model, :);
    S = summarize_support_occupancy(Ti, "M_eff_true_diag_bin", ...
        'BaseVars', ["model_name", "model_type"]);

    total_rows = height(Ti);
    total_outliers = sum(Ti.is_high_error_20);
    overall_rate = 100 * total_outliers / max(total_rows, 1);

    largest_rows = pick_support_row(S, 'Pct_of_all_rows');
    largest_rate = pick_support_row(S, 'HighError_gt20_pct');
    largest_outlier_share = pick_support_row(S, 'Pct_of_all_outliers');

    rows(i).model_name = model;
    rows(i).N_total_rows = total_rows;
    rows(i).N_outliers_gt20 = total_outliers;
    rows(i).Overall_HighError_gt20_pct = overall_rate;
    rows(i).bin_largest_support = string(largest_rows.M_eff_true_diag_bin);
    rows(i).largest_support_Pct_of_all_rows = largest_rows.Pct_of_all_rows;
    rows(i).bin_largest_error_rate = string(largest_rate.M_eff_true_diag_bin);
    rows(i).largest_HighError_gt20_pct = largest_rate.HighError_gt20_pct;
    rows(i).bin_largest_outlier_share = ...
        string(largest_outlier_share.M_eff_true_diag_bin);
    rows(i).largest_Pct_of_all_outliers = ...
        largest_outlier_share.Pct_of_all_outliers;
    rows(i).largest_outlier_share_bin_support_pct = ...
        largest_outlier_share.Pct_of_all_rows;
    rows(i).largest_outlier_share_bin_error_rate_pct = ...
        largest_outlier_share.HighError_gt20_pct;
end

T_support_report = struct2table(rows);

end

function row = pick_support_row(T, metric)

if isempty(T)
    row = table();
    row.M_eff_true_diag_bin = missing;
    row.Pct_of_all_rows = NaN;
    row.HighError_gt20_pct = NaN;
    row.Pct_of_all_outliers = NaN;
    return;
end

T = sortrows(T, metric, 'descend');
row = T(1, :);

end

function row = pick_worst(T, model, required_var)

Ti = T(T.model_name == model & T.model_type == "bagged_trees", :);
if isempty(Ti)
    row = table();
    row.(char(required_var)) = missing;
    row.MAPE_pct = NaN;
    row.HighError_gt20_pct = NaN;
    return;
end
Ti = sortrows(Ti, {'HighError_gt20_pct', 'MAPE_pct'}, ...
    {'descend', 'descend'});
row = Ti(1, :);

end

function print_auto_conclusions(T_report)

fprintf('\nAutomatic mini-conclusions.\n');
for i = 1:height(T_report)
    fprintf(['For %s, high errors are most concentrated in REQ_M=%s, ', ...
        'M_eff_guess_bin=%s, M_eff_true_diag_bin=%s, frequency=%s, ', ...
        'aperture step=%s, q_true_bin=%s.\n'], ...
        T_report.model_name(i), ...
        T_report.worst_REQ_M(i), ...
        T_report.worst_M_eff_guess_bin(i), ...
        T_report.worst_M_eff_true_diag_bin(i), ...
        T_report.worst_frequency(i), ...
        T_report.worst_aperture(i), ...
        T_report.worst_q_true_bin(i));
end

end

function print_support_conclusions(T_support_report)

fprintf('\nSupport-aware automatic conclusions using M_eff_true_diag_bin.\n');
for i = 1:height(T_support_report)
    fprintf(['For %s, %.3g%% of all test rows are high-error outliers. ', ...
        'Most outliers are concentrated in bin %s, which contains %.3g%% ', ...
        'of the test rows and has a high-error rate of %.3g%%. '], ...
        T_support_report.model_name(i), ...
        T_support_report.Overall_HighError_gt20_pct(i), ...
        T_support_report.bin_largest_outlier_share(i), ...
        T_support_report.largest_outlier_share_bin_support_pct(i), ...
        T_support_report.largest_outlier_share_bin_error_rate_pct(i));

    support_pct = T_support_report.largest_outlier_share_bin_support_pct(i);
    rate_pct = T_support_report.largest_outlier_share_bin_error_rate_pct(i);
    overall_pct = T_support_report.Overall_HighError_gt20_pct(i);

    if support_pct >= 40 && rate_pct <= 2 * max(overall_pct, eps)
        reason = "mostly because it contains many rows.";
    elseif support_pct < 20 && rate_pct > 2 * max(overall_pct, eps)
        reason = "mostly because it is a small but high-risk bin.";
    else
        reason = "because it combines substantial support with elevated error rate.";
    end
    fprintf('This suggests the bin dominates outliers %s\n', reason);
end

end

function T = select_existing_columns(T, cols)

cols = string(cols);
cols = cols(ismember(cols, string(T.Properties.VariableNames)));
T = T(:, cellstr(cols));

end

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

end

function export_figure(fig, fig_dir, out_name)

png_path = fullfile(fig_dir, string(out_name) + ".png");
pdf_path = fullfile(fig_dir, string(out_name) + ".pdf");
exportgraphics(fig, png_path, 'Resolution', 300, ...
    'BackgroundColor', 'white');
exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end
