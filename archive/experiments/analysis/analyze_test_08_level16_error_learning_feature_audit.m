%% analyze_test_08_level16_error_learning_feature_audit.m
% Level 16: residual learning, Srad-proxy features, and feature audit.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load Test 08

experiment_name = 'test_08_advanced_angular_features';

[T_mc, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

COL = adaptive_req.analysis.detect_mc_columns(T_mc);

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_16_error_learning_feature_audit');
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Build feature table

T_feat = add_ecum_features_from_mapping(T_mc);
ecum_features = string(fieldnames( ...
    adaptive_req.quantile.extract_ecum_shape_features( ...
        T_feat.req_mapping{find_first_mapping(T_feat)})));
srad_proxy_features = ecum_features(startsWith(ecum_features, "srad_proxy_"));

T_feat = adaptive_req.analysis.add_effective_window_metrics( ...
    T_feat, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', 3.0);

base_predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"
    "width_75_25_rel"
    "width_90_50_rel"
    "width_90_10_rel"
    "lowk_frac_rel"
    "midband_frac_rel"
    "highk_frac_rel"
    "circ_var"
    "dom_dir_frac"
    "window_max_frac"
    "window_cf"
];

full_predictors = keep_existing_predictors(T_feat, [base_predictors; ecum_features]);

%% Base model with Srad-proxy Ecum features

[~, T_base_pred, T_base_metrics] = ...
    adaptive_req.analysis.train_q_model_from_predictors( ...
        T_feat, ...
        full_predictors, ...
        'QVar', COL.q, ...
        'ModelName', "Level16_base_ecum_srad_proxy", ...
        'SplitMode', 'condition', ...
        'ConditionVar', COL.condition, ...
        'TrainFraction', 0.70, ...
        'RandomSeed', 1701, ...
        'ModelTypes', "bagged_trees", ...
        'NumLearningCycles', 300, ...
        'MinLeafSize', 8, ...
        'LearnRate', 0.05, ...
        'ClipPredictions', true, ...
        'Verbose', false);

T_base_pred.model_role = repmat("base", height(T_base_pred), 1);
T_split_ref = build_split_reference(T_base_pred);

%% Error-learning model: predict residual q_true - q_base_pred

T_residual_data = add_base_predictions_to_feature_table(T_feat, T_base_pred);
T_residual_data.q_residual_target = ...
    T_residual_data.q_base_true - T_residual_data.q_base_pred;
train_mask = T_residual_data.split == "train";

residual_predictors = keep_existing_predictors(T_residual_data, ...
    [full_predictors; "q_base_pred"]);

[~, T_residual_pred, T_residual_metrics] = ...
    adaptive_req.analysis.train_q_model_from_predictors( ...
        T_residual_data, ...
        residual_predictors, ...
        'QVar', 'q_residual_target', ...
        'ModelName', "Level16_residual_corrector", ...
        'SplitMode', 'condition', ...
        'ConditionVar', COL.condition, ...
        'TrainFraction', 0.70, ...
        'TrainMask', train_mask, ...
        'RandomSeed', 2701, ...
        'ModelTypes', "bagged_trees", ...
        'NumLearningCycles', 220, ...
        'MinLeafSize', 12, ...
        'LearnRate', 0.05, ...
        'ClipPredictions', false, ...
        'Verbose', false);

T_corrected_pred = build_corrected_prediction_table( ...
    T_base_pred, T_residual_pred);

T_all_predictions = [standardize_predictions(T_base_pred); ...
    standardize_predictions(T_corrected_pred)];
T_all_predictions = sortrows(T_all_predictions, ...
    {'model_name', 'model_type', 'condition_id', 'step_idx', ...
     'realization_idx', 'patch_idx'});

%% Convert to SWS and summarize

T_sws = add_sws_from_q_predictions(T_all_predictions, T_feat, "req_mapping");
T_sws_test = T_sws(T_sws.split == "test", :);

T_sws_metrics = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role"]);
T_sws_metrics = sortrows(T_sws_metrics, 'RMSE_pct', 'ascend');

T_sws_by_M = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"]);
T_sws_by_aperture = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "step_idx", "Omega_sr"]);

high_error_threshold_pct = 20;
T_high_error = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role"], high_error_threshold_pct);
T_high_error_by_M = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"], ...
    high_error_threshold_pct);
T_high_error_by_cs = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role", "SIM_cs_bg"], ...
    high_error_threshold_pct);

%% Feature audit

T_feature_assoc = compute_feature_audit_table(T_feat, full_predictors, COL.q);

T_ablation = run_drop_one_ablation( ...
    T_feat, ...
    full_predictors, ...
    T_split_ref, ...
    COL.q, ...
    COL.condition);

T_ablation_sws = add_sws_from_q_predictions( ...
    standardize_predictions(T_ablation.predictions), T_feat, "req_mapping");
T_ablation_sws_test = T_ablation_sws(T_ablation_sws.split == "test", :);
T_ablation_metrics = summarize_sws_metrics(T_ablation_sws_test, ...
    ["model_name", "model_type", "dropped_feature"]);
T_ablation_metrics = sortrows(T_ablation_metrics, 'RMSE_pct', 'ascend');

T_full_ablation_ref = T_ablation_metrics(T_ablation_metrics.dropped_feature == "NONE", :);
if ~isempty(T_full_ablation_ref)
    full_rmse_pct = T_full_ablation_ref.RMSE_pct(1);
    T_ablation_metrics.delta_RMSE_pct_vs_full = ...
        T_ablation_metrics.RMSE_pct - full_rmse_pct;
end

%% Save tables

writetable(T_base_metrics, fullfile(analysis_dir, ...
    'level16_base_q_metrics.csv'));
writetable(T_residual_metrics, fullfile(analysis_dir, ...
    'level16_residual_q_metrics.csv'));
writetable(remove_heavy_curve_columns(T_sws), fullfile(analysis_dir, ...
    'level16_sws_predictions.csv'));
writetable(T_sws_metrics, fullfile(analysis_dir, ...
    'level16_sws_metrics_by_model.csv'));
writetable(T_sws_by_M, fullfile(analysis_dir, ...
    'level16_sws_error_by_M.csv'));
writetable(T_sws_by_aperture, fullfile(analysis_dir, ...
    'level16_sws_error_by_aperture.csv'));
writetable(T_high_error, fullfile(analysis_dir, ...
    'level16_high_error_rate_by_model.csv'));
writetable(T_high_error_by_M, fullfile(analysis_dir, ...
    'level16_high_error_rate_by_model_M.csv'));
writetable(T_high_error_by_cs, fullfile(analysis_dir, ...
    'level16_high_error_rate_by_model_cs_bg.csv'));
writetable(T_feature_assoc, fullfile(analysis_dir, ...
    'level16_feature_q_error_associations.csv'));
writetable(T_ablation_metrics, fullfile(analysis_dir, ...
    'level16_drop_one_ablation_sws_metrics.csv'));

save(fullfile(analysis_dir, 'level_16_error_learning_feature_audit.mat'), ...
    'T_feat', ...
    'ecum_features', ...
    'srad_proxy_features', ...
    'full_predictors', ...
    'T_base_pred', ...
    'T_residual_pred', ...
    'T_corrected_pred', ...
    'T_sws', ...
    'T_sws_metrics', ...
    'T_sws_by_M', ...
    'T_sws_by_aperture', ...
    'T_high_error', ...
    'T_high_error_by_M', ...
    'T_high_error_by_cs', ...
    'T_feature_assoc', ...
    'T_ablation_metrics', ...
    '-v7.3');

%% Plots

plot_level16_model_comparison(T_sws_metrics, analysis_dir);
plot_level16_sws_scatter(T_sws_test, T_sws_metrics, analysis_dir);
plot_level16_error_learning(T_sws_test, analysis_dir);
plot_level16_high_error(T_high_error_by_M, T_high_error_by_cs, ...
    high_error_threshold_pct, analysis_dir);
plot_level16_feature_audit(T_feature_assoc, analysis_dir);
plot_level16_ablation(T_ablation_metrics, analysis_dir);

fprintf('\n============================================================\n');
fprintf('LEVEL 16 SWS METRICS, TEST SPLIT\n');
fprintf('============================================================\n');
disp(T_sws_metrics(:, {'model_name', 'model_type', 'n', ...
    'MAE_pct', 'RMSE_pct', 'bias_pct', 'p95_APE_pct', 'max_APE_pct'}));
fprintf('\nLevel 16 complete. Results saved to:\n%s\n', analysis_dir);

%% Local helpers

function idx = find_first_mapping(T)

idx = find(~cellfun(@isempty, T.req_mapping), 1, 'first');
if isempty(idx)
    error('No non-empty req_mapping entries were found.');
end

end

function T = add_ecum_features_from_mapping(T)

if ~ismember('req_mapping', T.Properties.VariableNames)
    error('Input table must contain req_mapping.');
end

n = height(T);
feat0 = adaptive_req.quantile.extract_ecum_shape_features( ...
    T.req_mapping{find_first_mapping(T)});
names = fieldnames(feat0);

for j = 1:numel(names)
    T.(names{j}) = NaN(n, 1);
end

for i = 1:n
    mapping_i = T.req_mapping{i};
    if isempty(mapping_i)
        continue;
    end

    feat_i = adaptive_req.quantile.extract_ecum_shape_features(mapping_i);
    for j = 1:numel(names)
        T.(names{j})(i) = feat_i.(names{j});
    end
end

end

function predictors_out = keep_existing_predictors(T, predictors_in)

vars = string(T.Properties.VariableNames);
predictors_in = string(predictors_in);
predictors_out = strings(0, 1);

for i = 1:numel(predictors_in)
    p_i = predictors_in(i);
    candidates = [p_i; p_i + "_mean"];
    if any(ismember(vars, candidates))
        predictors_out(end + 1, 1) = p_i; %#ok<AGROW>
    end
end

end

function T_ref = build_split_reference(T_pred)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_ref = T_pred(:, [key_vars, {'split'}]);
[~, ia] = unique(T_ref(:, key_vars), 'rows', 'stable');
T_ref = T_ref(ia, :);

end

function train_mask = build_train_mask_from_reference(T, T_ref)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_keys = T(:, key_vars);
T_keys.row_id_for_split = (1:height(T_keys)).';
T_join = innerjoin(T_keys, T_ref, 'Keys', key_vars);

if height(T_join) ~= height(T)
    error('Split reference does not cover all rows.');
end

train_mask = false(height(T), 1);
train_mask(T_join.row_id_for_split) = T_join.split == "train";

end

function T_out = add_base_predictions_to_feature_table(T_feat, T_base_pred)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_base = T_base_pred(:, [key_vars, {'split', 'q_true', 'q_pred'}]);
T_base.Properties.VariableNames{'q_true'} = 'q_base_true';
T_base.Properties.VariableNames{'q_pred'} = 'q_base_pred';
[~, ia] = unique(T_base(:, key_vars), 'rows', 'stable');
T_base = T_base(ia, :);
T_out = innerjoin(T_feat, T_base, 'Keys', key_vars);

end

function T_corr = build_corrected_prediction_table(T_base, T_resid)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_r = T_resid(:, [key_vars, {'q_pred'}]);
T_r.Properties.VariableNames{'q_pred'} = 'q_residual_pred';
T_corr = innerjoin(T_base, T_r, 'Keys', key_vars);
T_corr.q_base_pred = T_corr.q_pred;
T_corr.q_pred = min(max(T_corr.q_base_pred + T_corr.q_residual_pred, 0), 1);
T_corr.residual = T_corr.q_pred - T_corr.q_true;
T_corr.abs_error = abs(T_corr.residual);
T_corr.model_name = repmat("Level16_residual_corrected", height(T_corr), 1);
T_corr.model_type = repmat("bagged_trees", height(T_corr), 1);
T_corr.model_role = repmat("error_learning", height(T_corr), 1);

end

function T_std = standardize_predictions(T)

n = height(T);
T_std = table();

string_vars = {'model_name', 'model_type', 'model_role', 'split', ...
    'SIM_WaveModel', 'condition_label', 'dropped_feature'};
numeric_vars = {'condition_id', 'step_idx', 'realization_idx', ...
    'patch_idx', 'condition_position', 'SIM_f0', 'SIM_cs_bg', ...
    'REQ_M', 'Omega_sr', 'omega_mean', 'aperture_value', ...
    'q_true', 'q_pred', 'residual', 'abs_error', ...
    'q_base_pred', 'q_residual_pred'};

for i = 1:numel(string_vars)
    T_std.(string_vars{i}) = get_string_column(T, string_vars{i}, n);
end

for i = 1:numel(numeric_vars)
    T_std.(numeric_vars{i}) = get_numeric_column(T, numeric_vars{i}, n);
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

function T_out = add_sws_from_q_predictions(T_pred, T_ref, mapping_var)

key_vars = choose_key_vars(T_pred, T_ref);
T_curve = T_ref(:, cellstr([key_vars, mapping_var]));
[~, ia] = unique(T_curve(:, cellstr(key_vars)), 'rows', 'stable');
T_curve = T_curve(ia, :);

T_out = innerjoin(T_pred, T_curve, 'Keys', cellstr(key_vars));
n = height(T_out);

cs_pred = NaN(n, 1);
k_pred = NaN(n, 1);
cs_from_q_true = NaN(n, 1);
k_from_q_true = NaN(n, 1);

for i = 1:n
    curve_i = get_curve(T_out.(char(mapping_var)), i);
    [cs_pred(i), k_pred(i)] = cs_from_q_curve( ...
        T_out.q_pred(i), curve_i, T_out.SIM_f0(i));
    [cs_from_q_true(i), k_from_q_true(i)] = cs_from_q_curve( ...
        T_out.q_true(i), curve_i, T_out.SIM_f0(i));
end

T_out.k_pred_from_q = k_pred;
T_out.cs_pred_from_q = cs_pred;
T_out.k_from_q_true = k_from_q_true;
T_out.cs_from_q_true = cs_from_q_true;
T_out.cs_true = T_out.SIM_cs_bg;
T_out.q_error = T_out.q_pred - T_out.q_true;
T_out.q_abs_error = abs(T_out.q_error);
T_out.cs_error = T_out.cs_pred_from_q - T_out.cs_true;
T_out.cs_abs_error = abs(T_out.cs_error);
T_out.cs_rel_error = T_out.cs_error ./ T_out.cs_true;
T_out.cs_abs_rel_error = abs(T_out.cs_rel_error);
T_out.cs_error_pct = 100 * T_out.cs_rel_error;
T_out.cs_abs_error_pct = 100 * T_out.cs_abs_rel_error;

end

function key_vars = choose_key_vars(T_pred, T_ref)

candidates = {
    ["condition_id", "step_idx", "realization_idx", "patch_idx"]
    ["condition_position", "step_idx", "realization_idx", "patch_idx"]
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

error('Could not find common key variables.');

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
    error('Invalid REQ mapping at row %d.', idx);
end

end

function [cs, kq] = cs_from_q_curve(q, curve, f0)

cs = NaN;
kq = NaN;

if ~isfinite(q) || ~isfinite(f0) || f0 <= 0
    return;
end

k = double(curve.k_cent(:));
E = double(curve.Ecum(:));
valid = isfinite(E) & isfinite(k) & k > 0;
E = E(valid);
k = k(valid);

if numel(E) < 2
    return;
end

[k, ord] = sort(k, 'ascend');
E = cummax(E(ord));
[Euniq, ia] = unique(E, 'stable');
kuniq = k(ia);

if numel(Euniq) < 2
    return;
end

q = min(max(q, min(Euniq)), max(Euniq));
kq = interp1(Euniq, kuniq, q, 'linear', 'extrap');
if isfinite(kq) && kq > 0
    cs = 2*pi*f0/kq;
end

end

function T_metrics = summarize_sws_metrics(T, group_vars)

T = T(isfinite(T.cs_error_pct) & isfinite(T.cs_pred_from_q), :);
if isempty(T)
    T_metrics = table();
    return;
end

[G, T_metrics] = findgroups(T(:, cellstr(group_vars)));
T_metrics.n = splitapply(@numel, T.cs_error, G);
T_metrics.MAE_cs = splitapply(@mean_abs_finite, T.cs_error, G);
T_metrics.RMSE_cs = splitapply(@rmse_finite, T.cs_error, G);
T_metrics.bias_cs = splitapply(@mean_finite, T.cs_error, G);
T_metrics.MAE_pct = splitapply(@mean_abs_finite, T.cs_error_pct, G);
T_metrics.MAPE_pct = T_metrics.MAE_pct;
T_metrics.RMSE_pct = splitapply(@rmse_finite, T.cs_error_pct, G);
T_metrics.bias_pct = splitapply(@mean_finite, T.cs_error_pct, G);
T_metrics.median_APE_pct = splitapply(@median_abs_finite, ...
    T.cs_error_pct, G);
T_metrics.p95_APE_pct = splitapply(@p95_abs_finite, T.cs_error_pct, G);
T_metrics.max_APE_pct = splitapply(@max_abs_finite, T.cs_error_pct, G);
T_metrics.mean_q_abs_error = splitapply(@mean_finite, T.q_abs_error, G);

end

function T_rate = summarize_high_error_rates(T, group_vars, threshold_pct)

T = T(isfinite(T.cs_abs_error_pct), :);
[G, T_rate] = findgroups(T(:, cellstr(group_vars)));
is_high = T.cs_abs_error_pct > threshold_pct;
T_rate.n = splitapply(@numel, T.cs_abs_error_pct, G);
T_rate.n_high_error = splitapply(@sum, double(is_high), G);
T_rate.high_error_pct = 100 * T_rate.n_high_error ./ T_rate.n;
T_rate.threshold_pct = repmat(threshold_pct, height(T_rate), 1);

end

function T_assoc = compute_feature_audit_table(T, predictors, q_var)

rows(numel(predictors)) = struct( ...
    'feature_name', "", ...
    'spearman_q', NaN, ...
    'mutual_information_q', NaN, ...
    'missing_pct', NaN, ...
    'unique_values', NaN);
y = T.(char(q_var));

for i = 1:numel(predictors)
    name = predictors(i);
    x = double(T.(char(resolve_variable(T, name))));
    rows(i).feature_name = name;
    rows(i).spearman_q = spearman_finite(x, y);
    rows(i).mutual_information_q = discretized_mi(x, y, 8);
    rows(i).missing_pct = 100 * mean(~isfinite(x));
    rows(i).unique_values = numel(unique(x(isfinite(x))));
end

T_assoc = struct2table(rows);
T_assoc.abs_spearman_q = abs(T_assoc.spearman_q);
T_assoc = sortrows(T_assoc, ...
    {'mutual_information_q', 'abs_spearman_q'}, {'descend', 'descend'});

end

function out = run_drop_one_ablation(T, predictors, T_split_ref, q_var, condition_var)

train_mask = build_train_mask_from_reference(T, T_split_ref);
model_types = "bagged_trees";
pred_tables = table();
metric_tables = table();

sets = ["NONE"; predictors(:)];

for i = 1:numel(sets)
    dropped = sets(i);
    if dropped == "NONE"
        predictors_i = predictors;
        model_name = "Level16_ablation_full";
    else
        predictors_i = predictors(predictors ~= dropped);
        model_name = "Level16_drop_" + dropped;
    end

    fprintf('Level 16 ablation %d/%d: drop %s\n', ...
        i, numel(sets), dropped);

    [~, T_pred_i, T_metrics_i] = ...
        adaptive_req.analysis.train_q_model_from_predictors( ...
            T, ...
            predictors_i, ...
            'QVar', q_var, ...
            'ModelName', model_name, ...
            'SplitMode', 'condition', ...
            'ConditionVar', condition_var, ...
            'TrainFraction', 0.70, ...
            'TrainMask', train_mask, ...
            'RandomSeed', 3701, ...
            'ModelTypes', model_types, ...
            'NumLearningCycles', 160, ...
            'MinLeafSize', 10, ...
            'LearnRate', 0.05, ...
            'ClipPredictions', true, ...
            'Verbose', false);

    T_pred_i.model_role = repmat("drop_one_ablation", height(T_pred_i), 1);
    T_pred_i.dropped_feature = repmat(dropped, height(T_pred_i), 1);
    T_metrics_i.dropped_feature = repmat(dropped, height(T_metrics_i), 1);

    pred_tables = [pred_tables; standardize_predictions(T_pred_i)]; %#ok<AGROW>
    metric_tables = [metric_tables; T_metrics_i]; %#ok<AGROW>
end

out = struct();
out.predictions = pred_tables;
out.q_metrics = metric_tables;

end

function x = resolve_variable(T, name)

vars = string(T.Properties.VariableNames);
if ismember(name, vars)
    x = name;
elseif ismember(name + "_mean", vars)
    x = name + "_mean";
else
    error('Variable not found: %s', name);
end

end

function mi = discretized_mi(x, y, num_bins)

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < num_bins * 3 || range(x) <= 0 || range(y) <= 0
    mi = NaN;
    return;
end

xb = discretize_by_quantiles(x, num_bins);
yb = discretize_by_quantiles(y, num_bins);
valid = ~isnan(xb) & ~isnan(yb);
xb = xb(valid);
yb = yb(valid);

joint = accumarray([xb, yb], 1, [num_bins, num_bins], @sum, 0);
Pxy = joint / sum(joint(:));
Px = sum(Pxy, 2);
Py = sum(Pxy, 1);

mi = 0;
for a = 1:num_bins
    for b = 1:num_bins
        if Pxy(a, b) > 0 && Px(a) > 0 && Py(b) > 0
            mi = mi + Pxy(a, b) * log(Pxy(a, b) / (Px(a) * Py(b)));
        end
    end
end

end

function bin = discretize_by_quantiles(x, num_bins)

edges = quantile(x, linspace(0, 1, num_bins + 1));
edges(1) = -Inf;
edges(end) = Inf;
edges = unique(edges, 'stable');

if numel(edges) < 3
    bin = NaN(size(x));
    return;
end

bin = discretize(x, edges);
bin = min(bin, num_bins);

end

function T_light = remove_heavy_curve_columns(T)

T_light = T;
heavy_vars = {'req_curve', 'req_mapping', 'feat', 'feature_struct'};
for i = 1:numel(heavy_vars)
    if ismember(heavy_vars{i}, T_light.Properties.VariableNames)
        T_light.(heavy_vars{i}) = [];
    end
end

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

function plot_level16_model_comparison(T_metrics, output_dir)

T_plot = sortrows(T_metrics, 'RMSE_pct', 'ascend');
labels = categorical(strcat(T_plot.model_name, newline, T_plot.model_type));
labels = reordercats(labels, cellstr(string(labels)));

fig = figure('Color', 'w', 'Position', [100 100 1100 520]);
bar(labels, T_plot.RMSE_pct);
ylabel('SWS RMSE (%)');
title('Level 16: base vs residual-corrected model');
grid on;
xtickangle(35);
save_fig(fig, output_dir, 'level16_model_sws_rmse');

end

function plot_level16_sws_scatter(T, T_metrics, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
n_show = min(2, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1150 520]);
tl = tiledlayout(fig, 1, n_show, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:n_show
    Ti = T(T.model_name == T_metrics.model_name(i) & ...
        T.model_type == T_metrics.model_type(i), :);
    ax = nexttile(tl);
    scatter(ax, Ti.cs_true, Ti.cs_pred_from_q, 14, Ti.REQ_M, ...
        'filled', 'MarkerFaceAlpha', 0.35);
    hold(ax, 'on');
    lims = [min([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan'), ...
        max([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')];
    plot(ax, lims, lims, '--k');
    axis(ax, 'square');
    xlim(ax, lims);
    ylim(ax, lims);
    xlabel(ax, 'True c_s (m/s)');
    ylabel(ax, 'Predicted c_s (m/s)');
    title(ax, sprintf('%s\nRMSE %.2f%%', T_metrics.model_name(i), ...
        T_metrics.RMSE_pct(i)), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = 'REQ M';
    grid(ax, 'on');
end

title(tl, 'Level 16 SWS predictions', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'level16_sws_pred_vs_true');

end

function plot_level16_error_learning(T, output_dir)

T_corr = T(T.model_name == "Level16_residual_corrected", :);
if isempty(T_corr)
    return;
end

fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
scatter(ax, T_corr.q_base_pred, T_corr.q_residual_pred, 14, ...
    T_corr.REQ_M, 'filled', 'MarkerFaceAlpha', 0.35);
xlabel(ax, 'Base q prediction');
ylabel(ax, 'Predicted residual correction');
title(ax, 'What the error learner tries to add');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T_corr.q_error, T_corr.cs_error_pct, 14, T_corr.Omega_sr, ...
    'filled', 'MarkerFaceAlpha', 0.35);
xline(ax, 0, ':k');
yline(ax, 0, ':k');
xlabel(ax, 'Corrected q error');
ylabel(ax, 'Signed SWS error (%)');
title(ax, 'Corrected q error propagation');
cb = colorbar(ax);
cb.Label.String = 'Omega (sr)';
grid(ax, 'on');

title(tl, 'Level 16 residual-learning diagnostics', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'level16_error_learning_diagnostics');

end

function plot_level16_high_error(T_by_M, T_by_cs, threshold_pct, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
plot_rate_bars(ax, T_by_M, "REQ_M");
ylabel(ax, sprintf('APE > %.0f%% rate', threshold_pct));
xlabel(ax, 'REQ M');
title(ax, 'High-error rate by M');
grid(ax, 'on');

ax = nexttile(tl);
plot_rate_bars(ax, T_by_cs, "SIM_cs_bg");
ylabel(ax, sprintf('APE > %.0f%% rate', threshold_pct));
xlabel(ax, 'True c_s');
title(ax, 'High-error rate by true SWS');
grid(ax, 'on');

save_fig(fig, output_dir, 'level16_high_error_rates');

end

function plot_rate_bars(ax, T, x_var)

model_names = unique(T.model_name, 'stable');
x_values = unique(T.(char(x_var)), 'sorted');
Y = NaN(numel(x_values), numel(model_names));

for i = 1:numel(x_values)
    for j = 1:numel(model_names)
        mask = T.(char(x_var)) == x_values(i) & T.model_name == model_names(j);
        if any(mask)
            Y(i, j) = T.high_error_pct(find(mask, 1));
        end
    end
end

bar(ax, categorical(string(x_values)), Y);
legend(ax, model_names, 'Interpreter', 'none', 'Location', 'best');

end

function plot_level16_feature_audit(T_assoc, output_dir)

n_show = min(20, height(T_assoc));
T = T_assoc(1:n_show, :);
labels = categorical(T.feature_name);
labels = reordercats(labels, cellstr(T.feature_name));

fig = figure('Color', 'w', 'Position', [100 100 980 720]);
barh(labels, T.mutual_information_q);
xlabel('Mutual information with q');
ylabel('Feature');
title('Level 16 feature association audit');
grid on;
save_fig(fig, output_dir, 'level16_feature_association');

end

function plot_level16_ablation(T_ablation, output_dir)

T = T_ablation(T_ablation.dropped_feature ~= "NONE", :);
T = sortrows(T, 'delta_RMSE_pct_vs_full', 'ascend');
n_show = min(25, height(T));
T = T(1:n_show, :);
labels = categorical(T.dropped_feature);
labels = reordercats(labels, cellstr(T.dropped_feature));

fig = figure('Color', 'w', 'Position', [100 100 1100 720]);
barh(labels, T.delta_RMSE_pct_vs_full);
xline(0, '--k');
xlabel('Delta SWS RMSE (%) after dropping feature');
ylabel('Dropped feature');
title(['Drop-one ablation: negative means the feature may hurt; ', ...
    'positive means it helps']);
grid on;
save_fig(fig, output_dir, 'level16_drop_one_ablation');

end

function save_fig(fig, output_dir, base_name)

exportgraphics(fig, fullfile(output_dir, base_name + ".png"), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, base_name + ".pdf"), ...
    'ContentType', 'vector');
close(fig);

end
