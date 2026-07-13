%% analyze_test_08_level15_ecum_features_and_plane_audit.m
% Level 15: REQ cumulative-curve features and in-plane-wave audit.

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

[T_mc, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

COL = adaptive_req.analysis.detect_mc_columns(T_mc);

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_15_ecum_features_and_plane_audit');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Add Ecum shape features from lightweight req_mapping

T_ecum = add_ecum_features_from_mapping(T_mc);
ecum_features = string(fieldnames( ...
    adaptive_req.quantile.extract_ecum_shape_features( ...
        T_ecum.req_mapping{find_first_mapping(T_ecum)})));

%% Fast in-plane coverage audit for current Fibonacci cone schedule

schedule = adaptive_req.simulate.build_aperture_schedule( ...
    'cone', MC.CFG.EXP.num_steps);

T_plane_no_force = adaptive_req.simulate.estimate_fibonacci_cone_plane_coverage( ...
    'Nwaves', MC.CFG.SIM.Nwaves, ...
    'ConeAxis', [1 0 0], ...
    'FullApertureDeg', schedule.values, ...
    'Tolerance', 1e-6, ...
    'ForceInPlaneWave', false);

T_plane_force = adaptive_req.simulate.estimate_fibonacci_cone_plane_coverage( ...
    'Nwaves', MC.CFG.SIM.Nwaves, ...
    'ConeAxis', [1 0 0], ...
    'FullApertureDeg', schedule.values, ...
    'Tolerance', 1e-6, ...
    'ForceInPlaneWave', true);

writetable(T_plane_no_force, fullfile(analysis_dir, ...
    'fibonacci_cone_plane_audit_no_force.csv'));
writetable(T_plane_force, fullfile(analysis_dir, ...
    'fibonacci_cone_plane_audit_force_in_plane.csv'));

%% Feature-q dependency analysis

[T_assoc_global, T_binned_global] = ...
    adaptive_req.analysis.compute_groupwise_feature_q_associations( ...
        T_ecum, ...
        ecum_features, ...
        'QVar', COL.q, ...
        'GroupVars', string.empty(1, 0), ...
        'NumBins', 8, ...
        'MinN', 20, ...
        'Verbose', false);

[T_assoc_by_M, T_binned_by_M] = ...
    adaptive_req.analysis.compute_groupwise_feature_q_associations( ...
        T_ecum, ...
        ecum_features, ...
        'QVar', COL.q, ...
        'GroupVars', "REQ_M", ...
        'NumBins', 8, ...
        'MinN', 20, ...
        'Verbose', false);

T_mi = compute_feature_mi_table(T_ecum, ecum_features, COL.q, 8);
T_assoc_global = outerjoin(T_assoc_global, T_mi, ...
    'Keys', 'feature_name', ...
    'MergeKeys', true);
T_assoc_global.abs_spearman_rho = abs(T_assoc_global.spearman_rho);
T_assoc_global = sortrows(T_assoc_global, ...
    {'mutual_information', 'abs_spearman_rho'}, ...
    {'descend', 'descend'});

writetable(T_assoc_global, fullfile(analysis_dir, ...
    'ecum_feature_q_association_global.csv'));
writetable(T_assoc_by_M, fullfile(analysis_dir, ...
    'ecum_feature_q_association_by_M.csv'));
writetable(T_binned_global, fullfile(analysis_dir, ...
    'ecum_feature_q_binned_global.csv'));
writetable(T_binned_by_M, fullfile(analysis_dir, ...
    'ecum_feature_q_binned_by_M.csv'));

%% Model ablation: Model C vs Model C + Ecum shape

T_mc_eff = adaptive_req.analysis.add_effective_window_metrics( ...
    T_ecum, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', 3.0);

if ~ismember("M_eff", string(T_mc_eff.Properties.VariableNames))
    T_mc_eff.M_eff = T_mc_eff.M_eff_true_diag;
end

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

advanced_angular_predictors = [
    "ang_moment_2"
    "ang_moment_4"
    "ang_peak_count_rel"
    "ang_top2_window_frac"
    "ang_top3_window_frac"
    "ang_top2_to_top1"
    "ang_peak_separation_deg"
];

model_specs = struct([]);
model_specs(1).name = "ModelC_baseline";
model_specs(1).role = "operational_baseline";
model_specs(1).predictors = base_predictors;

model_specs(2).name = "ModelH_ecum_shape";
model_specs(2).role = "operational_ecum";
model_specs(2).predictors = [base_predictors; ecum_features];

model_specs(3).name = "ModelI_angular_ecum_shape";
model_specs(3).role = "operational_angular_ecum";
model_specs(3).predictors = [base_predictors; advanced_angular_predictors; ecum_features];

model_types = ["boosted_trees", "bagged_trees"];

T_all_metrics = table();
T_all_predictions = table();

for sidx = 1:numel(model_specs)

    spec = model_specs(sidx);
    predictors_i = keep_existing_predictors(T_mc_eff, spec.predictors);

    fprintf('\nTraining %s\n', spec.name);

    [~, T_pred_i, T_metrics_i] = ...
        adaptive_req.analysis.train_q_model_from_predictors( ...
            T_mc_eff, ...
            predictors_i, ...
            'QVar', COL.q, ...
            'ModelName', spec.name, ...
            'SplitMode', 'condition', ...
            'ConditionVar', COL.condition, ...
            'TrainFraction', 0.70, ...
            'RandomSeed', 1701, ...
            'ModelTypes', model_types, ...
            'NumLearningCycles', 300, ...
            'MinLeafSize', 8, ...
            'LearnRate', 0.05, ...
            'ClipPredictions', true, ...
            'Verbose', false);

    T_metrics_i.model_role = repmat(spec.role, height(T_metrics_i), 1);
    T_metrics_i.trained_M = NaN(height(T_metrics_i), 1);
    T_pred_i.model_role = repmat(spec.role, height(T_pred_i), 1);
    T_pred_i.trained_M = NaN(height(T_pred_i), 1);
    T_pred_i.predictor_set = repmat(strjoin(predictors_i, ", "), ...
        height(T_pred_i), 1);

    T_all_metrics = [T_all_metrics; T_metrics_i]; %#ok<AGROW>
    T_all_predictions = [T_all_predictions; standardize_predictions(T_pred_i)]; %#ok<AGROW>
end

%% Model ablation: train a separate Ecum model for each REQ M

by_M_predictors = keep_existing_predictors(T_mc_eff, ...
    [base_predictors; ecum_features]);
M_values_for_training = unique(T_mc_eff.REQ_M(isfinite(T_mc_eff.REQ_M)), ...
    'sorted');
global_split_reference = make_global_split_reference(T_all_predictions);

for midx = 1:numel(M_values_for_training)

    M_i = M_values_for_training(midx);
    T_m = T_mc_eff(T_mc_eff.REQ_M == M_i, :);
    train_mask_m = build_train_mask_from_reference(T_m, ...
        global_split_reference);

    fprintf('\nTraining ModelJ_ecum_by_M for M = %.0f\n', M_i);

    [~, T_pred_i, T_metrics_i] = ...
        adaptive_req.analysis.train_q_model_from_predictors( ...
            T_m, ...
            by_M_predictors, ...
            'QVar', COL.q, ...
            'ModelName', "ModelJ_ecum_by_M", ...
            'SplitMode', 'condition', ...
            'ConditionVar', COL.condition, ...
            'TrainFraction', 0.70, ...
            'TrainMask', train_mask_m, ...
            'RandomSeed', 1701 + midx, ...
            'ModelTypes', model_types, ...
            'NumLearningCycles', 300, ...
            'MinLeafSize', 8, ...
            'LearnRate', 0.05, ...
            'ClipPredictions', true, ...
            'Verbose', false);

    T_metrics_i.model_role = repmat("operational_ecum_by_M", ...
        height(T_metrics_i), 1);
    T_metrics_i.trained_M = repmat(M_i, height(T_metrics_i), 1);

    T_pred_i.model_role = repmat("operational_ecum_by_M", ...
        height(T_pred_i), 1);
    T_pred_i.trained_M = repmat(M_i, height(T_pred_i), 1);
    T_pred_i.predictor_set = repmat(strjoin(by_M_predictors, ", "), ...
        height(T_pred_i), 1);

    T_all_metrics = [T_all_metrics; T_metrics_i]; %#ok<AGROW>
    T_all_predictions = [T_all_predictions; standardize_predictions(T_pred_i)]; %#ok<AGROW>
end

T_test_metrics = T_all_metrics(T_all_metrics.split == "test", :);
T_test_metrics = sortrows(T_test_metrics, 'RMSE', 'ascend');

T_by_M = summarize_prediction_errors(T_all_predictions, ...
    ["model_name", "model_type", "model_role", "REQ_M"]);

%% Convert q predictions into SWS and diagnose failures

T_sws = add_sws_from_q_predictions(T_all_predictions, T_ecum, "req_mapping");
T_sws_test = T_sws(T_sws.split == "test", :);

T_sws_metrics = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role"]);
T_sws_metrics = sortrows(T_sws_metrics, 'RMSE_pct', 'ascend');

T_sws_by_M = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"]);
T_sws_by_aperture = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "step_idx", "Omega_sr"]);
T_sws_by_condition = summarize_sws_metrics(T_sws_test, ...
    ["model_name", "model_type", "model_role", "condition_id", ...
     "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M"]);

high_error_threshold_pct = 20;
T_high_error_by_model_M = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role", "REQ_M"], ...
    high_error_threshold_pct);
T_high_error_by_model_cs = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role", "SIM_cs_bg"], ...
    high_error_threshold_pct);
T_high_error_by_model = summarize_high_error_rates(T_sws_test, ...
    ["model_name", "model_type", "model_role"], ...
    high_error_threshold_pct);
T_high_error_by_model_M = add_high_error_share_within_model( ...
    T_high_error_by_model_M);
T_high_error_by_model_cs = add_high_error_share_within_model( ...
    T_high_error_by_model_cs);

[T_failure_corr, T_failure_outliers] = build_failure_diagnostics( ...
    T_sws_test, T_sws_metrics);

T_sws_light = remove_heavy_curve_columns(T_sws);

writetable(T_all_metrics, fullfile(analysis_dir, ...
    'level15_all_model_metrics.csv'));
writetable(T_test_metrics, fullfile(analysis_dir, ...
    'level15_test_metrics_sorted_by_RMSE.csv'));
writetable(T_all_predictions, fullfile(analysis_dir, ...
    'level15_all_model_predictions.csv'));
writetable(T_by_M, fullfile(analysis_dir, ...
    'level15_q_error_by_M.csv'));
writetable(T_sws_light, fullfile(analysis_dir, ...
    'level15_sws_predictions.csv'));
writetable(T_sws_metrics, fullfile(analysis_dir, ...
    'level15_sws_metrics_by_model.csv'));
writetable(T_sws_by_M, fullfile(analysis_dir, ...
    'level15_sws_error_by_M.csv'));
writetable(T_sws_by_aperture, fullfile(analysis_dir, ...
    'level15_sws_error_by_aperture.csv'));
writetable(T_sws_by_condition, fullfile(analysis_dir, ...
    'level15_sws_error_by_condition.csv'));
writetable(T_failure_corr, fullfile(analysis_dir, ...
    'level15_failure_correlations.csv'));
writetable(T_failure_outliers, fullfile(analysis_dir, ...
    'level15_failure_outliers.csv'));
writetable(T_high_error_by_model, fullfile(analysis_dir, ...
    'level15_high_error_rate_by_model.csv'));
writetable(T_high_error_by_model_M, fullfile(analysis_dir, ...
    'level15_high_error_rate_by_model_M.csv'));
writetable(T_high_error_by_model_cs, fullfile(analysis_dir, ...
    'level15_high_error_rate_by_model_cs_bg.csv'));

plot_top_ecum_features(T_ecum, T_assoc_global, COL.q, analysis_dir);
plot_level15_rmse_by_M(T_by_M, analysis_dir);
plot_sws_rmse_by_model(T_sws_metrics, analysis_dir);
plot_sws_pred_vs_true_best_models(T_sws_test, T_sws_metrics, analysis_dir, 6);
plot_sws_error_vs_true_speed_best_models(T_sws_test, T_sws_metrics, analysis_dir, 6);
plot_sws_error_boxplots(T_sws_test, T_sws_metrics, analysis_dir);
plot_best_model_failure_diagnostics(T_sws_test, T_sws_metrics, analysis_dir);
plot_failure_correlation_summary(T_failure_corr, analysis_dir);
plot_high_error_rates(T_high_error_by_model_M, T_high_error_by_model_cs, ...
    high_error_threshold_pct, analysis_dir);

save(fullfile(analysis_dir, 'level_15_ecum_features_and_plane_audit.mat'), ...
    'T_ecum', ...
    'T_plane_no_force', ...
    'T_plane_force', ...
    'T_assoc_global', ...
    'T_assoc_by_M', ...
    'T_mi', ...
    'T_all_metrics', ...
    'T_test_metrics', ...
    'T_by_M', ...
    'T_sws', ...
    'T_sws_metrics', ...
    'T_sws_by_M', ...
    'T_sws_by_aperture', ...
    'T_sws_by_condition', ...
    'T_failure_corr', ...
    'T_failure_outliers', ...
    'T_high_error_by_model', ...
    'T_high_error_by_model_M', ...
    'T_high_error_by_model_cs', ...
    'high_error_threshold_pct', ...
    '-v7.3');

fprintf('\n============================================================\n');
fprintf('LEVEL 15 SWS METRICS, TEST SPLIT\n');
fprintf('============================================================\n');
disp(T_sws_metrics(:, { ...
    'model_name', 'model_type', 'n', ...
    'MAE_pct', 'RMSE_pct', 'bias_pct', 'R2_cs', 'spearman_cs'}));

fprintf('\nLevel 15 complete. Results saved to:\n%s\n', analysis_dir);

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

function T_mi = compute_feature_mi_table(T, features, q_var, num_bins)

rows(numel(features)) = struct();
y = T.(char(q_var));

for i = 1:numel(features)
    x = T.(char(features(i)));
    rows(i).feature_name = features(i);
    rows(i).mutual_information = discretized_mi(x, y, num_bins);
end

T_mi = struct2table(rows);

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

function T_std = standardize_predictions(T)

n = height(T);
T_std = table();

string_vars = {'model_name', 'model_type', 'model_role', 'split', ...
    'SIM_WaveModel', 'condition_label', 'predictor_set'};
numeric_vars = {'condition_id', 'step_idx', 'realization_idx', ...
    'patch_idx', 'condition_position', 'SIM_f0', 'SIM_cs_bg', ...
    'REQ_M', 'Omega_sr', 'omega_mean', 'aperture_value', ...
    'trained_M', 'q_true', 'q_pred', 'residual', 'abs_error', ...
    'ecum_at_k0_rel', 'ecum_slope_at_k0_rel', ...
    'ecum_slope_peak_to_mean', 'ecum_width_50_rel', ...
    'ecum_width_80_rel', 'ecum_increment_entropy', ...
    'ecum_increment_peak_frac', 'ecum_low_tail_energy_rel', ...
    'ecum_asymmetry_25_75', 'ecum_width_ratio_80_50', ...
    'ecum_lower_upper_width_ratio', 'ecum_increment_gini', ...
    'ecum_slope_iqr_to_median'};

for i = 1:numel(string_vars)
    name = string_vars{i};
    T_std.(name) = get_string_column(T, name, n);
end

for i = 1:numel(numeric_vars)
    name = numeric_vars{i};
    T_std.(name) = get_numeric_column(T, name, n);
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

function T_ref = make_global_split_reference(T_pred)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
mask = T_pred.model_name == "ModelH_ecum_shape" & ...
    T_pred.model_type == "bagged_trees";

if ~any(mask)
    error('Could not build global split reference from ModelH bagged trees.');
end

T_ref = T_pred(mask, [key_vars, {'split'}]);
[~, ia] = unique(T_ref(:, key_vars), 'rows', 'stable');
T_ref = T_ref(ia, :);

end

function train_mask = build_train_mask_from_reference(T, T_ref)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_keys = T(:, key_vars);
T_keys.row_id_for_split = (1:height(T_keys)).';

T_join = innerjoin(T_keys, T_ref, 'Keys', key_vars);

if height(T_join) ~= height(T)
    error('Global split reference does not cover all rows for by-M training.');
end

train_mask = false(height(T), 1);
train_mask(T_join.row_id_for_split) = T_join.split == "train";

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
T_out.cs_from_true_q_error_pct = ...
    100 * (T_out.cs_from_q_true - T_out.cs_true) ./ T_out.cs_true;

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

error('Could not find common key variables between predictions and reference table.');

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
    error('Invalid REQ curve/mapping at row %d.', idx);
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
Ecum = cummax(Ecum(idx_sort));
[Ecum_unique, idx_unique] = unique(Ecum, 'stable');
k_unique = k_cent(idx_unique);

if numel(Ecum_unique) < 2
    return;
end

q_clamped = min(max(q, min(Ecum_unique)), max(Ecum_unique));
kq = interp1(Ecum_unique, k_unique, q_clamped, 'linear', 'extrap');

if isfinite(kq) && kq > 0
    cs = 2*pi*f0/kq;
end

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
    error('Could not find cumulative energy field in REQ mapping.');
end

if isempty(k_cent)
    error('Could not find k-center field in REQ mapping.');
end

Ecum = double(Ecum(:));
k_cent = double(k_cent(:));

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
T_metrics.R2_cs = splitapply(@r2_finite, T.cs_true, T.cs_pred_from_q, G);
T_metrics.spearman_cs = splitapply(@spearman_finite, ...
    T.cs_true, T.cs_pred_from_q, G);

end

function T_rate = summarize_high_error_rates(T, group_vars, threshold_pct)

T = T(isfinite(T.cs_abs_error_pct), :);

if isempty(T)
    T_rate = table();
    return;
end

[G, T_rate] = findgroups(T(:, cellstr(group_vars)));
is_high = T.cs_abs_error_pct > threshold_pct;

T_rate.n = splitapply(@numel, T.cs_abs_error_pct, G);
T_rate.n_high_error = splitapply(@sum, double(is_high), G);
T_rate.high_error_pct = 100 * T_rate.n_high_error ./ T_rate.n;
T_rate.threshold_pct = repmat(threshold_pct, height(T_rate), 1);
T_rate.median_APE_pct = splitapply(@median_abs_finite, ...
    T.cs_error_pct, G);
T_rate.p95_APE_pct = splitapply(@p95_abs_finite, T.cs_error_pct, G);
T_rate.max_APE_pct = splitapply(@max_abs_finite, T.cs_error_pct, G);

end

function T_rate = add_high_error_share_within_model(T_rate)

model_vars = {'model_name', 'model_type', 'model_role'};

if isempty(T_rate) || ~all(ismember(model_vars, T_rate.Properties.VariableNames))
    T_rate.high_error_share_within_model_pct = NaN(height(T_rate), 1);
    return;
end

[G, ~] = findgroups(T_rate(:, model_vars));
total_high = splitapply(@sum, T_rate.n_high_error, G);
T_rate.high_error_share_within_model_pct = ...
    100 * T_rate.n_high_error ./ total_high(G);
T_rate.high_error_share_within_model_pct(total_high(G) == 0) = 0;

end

function [T_corr, T_outliers] = build_failure_diagnostics(T, T_metrics)

if isempty(T_metrics)
    T_corr = table();
    T_outliers = table();
    return;
end

best_name = T_metrics.model_name(1);
best_type = T_metrics.model_type(1);
T_best = T(T.model_name == best_name & T.model_type == best_type, :);

T_outliers = sortrows(T_best, 'cs_abs_error_pct', 'descend');
keep_vars = ["model_name", "model_type", "condition_id", "step_idx", ...
    "realization_idx", "patch_idx", "SIM_WaveModel", "SIM_f0", ...
    "SIM_cs_bg", "REQ_M", "Omega_sr", "q_true", "q_pred", ...
    "q_error", "q_abs_error", "cs_true", "cs_pred_from_q", ...
    "cs_error_pct", "cs_abs_error_pct", "ecum_at_k0_rel", ...
    "ecum_slope_at_k0_rel", "ecum_slope_peak_to_mean", ...
    "ecum_width_50_rel", "ecum_width_80_rel", ...
    "ecum_increment_entropy", "ecum_increment_peak_frac", ...
    "ecum_low_tail_energy_rel", "ecum_asymmetry_25_75", ...
    "ecum_width_ratio_80_50", "ecum_lower_upper_width_ratio", ...
    "ecum_increment_gini", "ecum_slope_iqr_to_median"];
keep_vars = keep_vars(ismember(keep_vars, string(T_outliers.Properties.VariableNames)));
T_outliers = T_outliers(:, cellstr(keep_vars));

diagnostic_vars = ["REQ_M", "Omega_sr", "SIM_f0", "SIM_cs_bg", ...
    "q_true", "q_pred", "q_error", "q_abs_error", ...
    "ecum_at_k0_rel", "ecum_slope_at_k0_rel", ...
    "ecum_slope_peak_to_mean", "ecum_width_50_rel", ...
    "ecum_width_80_rel", "ecum_increment_entropy", ...
    "ecum_increment_peak_frac", "ecum_low_tail_energy_rel", ...
    "ecum_asymmetry_25_75", "ecum_width_ratio_80_50", ...
    "ecum_lower_upper_width_ratio", "ecum_increment_gini", ...
    "ecum_slope_iqr_to_median"];
diagnostic_vars = diagnostic_vars(ismember(diagnostic_vars, ...
    string(T_best.Properties.VariableNames)));

rows(numel(diagnostic_vars)) = struct( ...
    'model_name', "", ...
    'model_type', "", ...
    'feature_name', "", ...
    'spearman_abs_sws_error', NaN, ...
    'spearman_signed_sws_error', NaN, ...
    'spearman_abs_q_error', NaN, ...
    'n', NaN);

for i = 1:numel(diagnostic_vars)
    x = T_best.(char(diagnostic_vars(i)));
    rows(i).model_name = best_name;
    rows(i).model_type = best_type;
    rows(i).feature_name = diagnostic_vars(i);
    rows(i).spearman_abs_sws_error = spearman_finite( ...
        x, T_best.cs_abs_error_pct);
    rows(i).spearman_signed_sws_error = spearman_finite( ...
        x, T_best.cs_error_pct);
    rows(i).spearman_abs_q_error = spearman_finite(x, T_best.q_abs_error);
    rows(i).n = sum(isfinite(x) & isfinite(T_best.cs_abs_error_pct));
end

T_corr = struct2table(rows);
valid_corr = isfinite(T_corr.spearman_abs_sws_error) & T_corr.n > 0;
T_corr = T_corr(valid_corr, :);
T_corr.abs_spearman_abs_sws_error = abs(T_corr.spearman_abs_sws_error);
T_corr = sortrows(T_corr, 'abs_spearman_abs_sws_error', 'descend');

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

function y = median_abs_finite(x)

x = abs(x(isfinite(x)));
if isempty(x)
    y = NaN;
else
    y = median(x);
end

end

function y = p95_abs_finite(x)

x = abs(x(isfinite(x)));
if isempty(x)
    y = NaN;
else
    y = prctile(x, 95);
end

end

function y = max_abs_finite(x)

x = abs(x(isfinite(x)));
if isempty(x)
    y = NaN;
else
    y = max(x);
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
    y = 1 - ss_res / ss_tot;
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

function T_summary = summarize_prediction_errors(T, group_vars)

T = T(T.split == "test", :);
[G, T_summary] = findgroups(T(:, cellstr(group_vars)));

T_summary.n = splitapply(@numel, T.residual, G);
T_summary.MAE_q = splitapply(@(x) mean(abs(x), 'omitnan'), T.residual, G);
T_summary.RMSE_q = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.residual, G);
T_summary.bias_q = splitapply(@(x) mean(x, 'omitnan'), T.residual, G);

end

function plot_top_ecum_features(T, T_assoc, q_var, output_dir)

n_show = min(4, height(T_assoc));

fig = figure('Color', 'w', 'Position', [100 100 1200 800]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_show
    ax = nexttile(tl);
    feat = char(T_assoc.feature_name(i));
    scatter(ax, T.(feat), T.(char(q_var)), 12, T.REQ_M, 'filled', ...
        'MarkerFaceAlpha', 0.35);
    xlabel(ax, feat, 'Interpreter', 'none');
    ylabel(ax, char(q_var), 'Interpreter', 'none');
    title(ax, sprintf('MI = %.3f, Spearman = %.3f', ...
        T_assoc.mutual_information(i), T_assoc.spearman_rho(i)));
    grid(ax, 'on');
end

title(tl, 'Top Ecum feature associations');
exportgraphics(fig, fullfile(output_dir, ...
    'top_ecum_feature_q_associations.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, ...
    'top_ecum_feature_q_associations.pdf'), 'ContentType', 'vector');
close(fig);

end

function plot_level15_rmse_by_M(T_by_M, output_dir)

T = T_by_M(T_by_M.model_type == "bagged_trees", :);
model_names = unique(T.model_name, 'stable');
M_values = unique(T.REQ_M, 'sorted');
Y = NaN(numel(M_values), numel(model_names));

for i = 1:numel(M_values)
    for j = 1:numel(model_names)
        mask = T.REQ_M == M_values(i) & T.model_name == model_names(j);

        if any(mask)
            Y(i, j) = T.RMSE_q(find(mask, 1));
        end
    end
end

fig = figure('Color', 'w', 'Position', [100 100 950 550]);
bar(M_values, Y);
xlabel('REQ M');
ylabel('Test RMSE_q');
title('Level 15 bagged-tree ablation by M');
legend(model_names, 'Interpreter', 'none', 'Location', 'best');
grid on;

exportgraphics(fig, fullfile(output_dir, ...
    'level15_bagged_tree_RMSE_by_M.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, ...
    'level15_bagged_tree_RMSE_by_M.pdf'), 'ContentType', 'vector');
close(fig);

end

function plot_sws_rmse_by_model(T_metrics, output_dir)

T_plot = sortrows(T_metrics, 'RMSE_pct', 'ascend');
labels = strcat(T_plot.model_name, newline, T_plot.model_type);
labels = categorical(labels);
labels = reordercats(labels, cellstr(string(labels)));

fig = figure('Color', 'w', 'Position', [100 100 1150 520]);
ax = axes(fig);
bar(ax, labels, T_plot.RMSE_pct);
ylabel(ax, 'SWS RMSE (%)');
xlabel(ax, 'Model');
title(ax, 'Level 15 SWS error from Ecum-feature q models');
grid(ax, 'on');
ax.FontSize = 12;
xtickangle(ax, 45);
save_level15_figure(fig, output_dir, 'level15_sws_rmse_by_model');

end

function plot_sws_pred_vs_true_best_models(T_test, T_metrics, output_dir, n_best)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
n_best = min(n_best, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1250 760]);
tl = tiledlayout(fig, 2, ceil(n_best / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_best
    Ti = select_model_rows(T_test, T_metrics.model_name(i), ...
        T_metrics.model_type(i));

    ax = nexttile(tl);
    scatter(ax, Ti.cs_true, Ti.cs_pred_from_q, 13, Ti.REQ_M, 'filled', ...
        'MarkerFaceAlpha', 0.35);
    hold(ax, 'on');
    lims = [
        min([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')
        max([Ti.cs_true; Ti.cs_pred_from_q], [], 'omitnan')
    ];
    plot(ax, lims, lims, '--k', 'LineWidth', 1.1);
    xlim(ax, lims);
    ylim(ax, lims);
    axis(ax, 'square');
    xlabel(ax, 'True c_s (m/s)');
    ylabel(ax, 'Predicted c_s (m/s)');
    title(ax, sprintf('%s\n%s, RMSE = %.2f%%', ...
        T_metrics.model_name(i), T_metrics.model_type(i), ...
        T_metrics.RMSE_pct(i)), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = 'REQ M';
    grid(ax, 'on');
end

title(tl, 'Predicted SWS from Level 15 q models', 'FontWeight', 'bold');
save_level15_figure(fig, output_dir, 'level15_sws_pred_vs_true_best_models');

end

function plot_sws_error_vs_true_speed_best_models(T_test, T_metrics, ...
    output_dir, n_best)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
n_best = min(n_best, height(T_metrics));

fig = figure('Color', 'w', 'Position', [100 100 1250 760]);
tl = tiledlayout(fig, 2, ceil(n_best / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n_best
    Ti = select_model_rows(T_test, T_metrics.model_name(i), ...
        T_metrics.model_type(i));

    ax = nexttile(tl);
    scatter(ax, Ti.cs_true, Ti.cs_error_pct, 13, Ti.Omega_sr, 'filled', ...
        'MarkerFaceAlpha', 0.35);
    hold(ax, 'on');
    yline(ax, 0, '--k', 'LineWidth', 1.1);
    xlabel(ax, 'True c_s (m/s)');
    ylabel(ax, 'Signed SWS error (%)');
    title(ax, sprintf('%s\n%s, bias = %.2f%%', ...
        T_metrics.model_name(i), T_metrics.model_type(i), ...
        T_metrics.bias_pct(i)), 'Interpreter', 'none');
    cb = colorbar(ax);
    cb.Label.String = 'Omega (sr)';
    grid(ax, 'on');
end

title(tl, 'Level 15 signed SWS error', 'FontWeight', 'bold');
save_level15_figure(fig, output_dir, ...
    'level15_sws_error_vs_true_speed_best_models');

end

function plot_sws_error_boxplots(T_test, T_metrics, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
best_name = T_metrics.model_name(1);
best_type = T_metrics.model_type(1);
T_best = select_model_rows(T_test, best_name, best_type);

fig = figure('Color', 'w', 'Position', [100 100 1350 760]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
boxchart(ax, categorical(T_test.model_name), T_test.cs_abs_error_pct);
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Patch errors by model');
grid(ax, 'on');
xtickangle(ax, 35);

ax = nexttile(tl);
boxchart(ax, categorical(string(T_best.REQ_M)), T_best.cs_abs_error_pct);
xlabel(ax, 'REQ M');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, sprintf('Best model by M: %s / %s', best_name, best_type), ...
    'Interpreter', 'none');
grid(ax, 'on');

ax = nexttile(tl);
boxchart(ax, categorical(string(T_best.step_idx)), T_best.cs_abs_error_pct);
xlabel(ax, 'Aperture step');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Best model by aperture step');
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T_best.Omega_sr, T_best.cs_abs_error_pct, 16, ...
    T_best.REQ_M, 'filled', 'MarkerFaceAlpha', 0.35);
xlabel(ax, 'Omega (sr)');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Best model error vs aperture');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

title(tl, 'Level 15 SWS error distributions', 'FontWeight', 'bold');
save_level15_figure(fig, output_dir, 'level15_sws_error_boxplots');

end

function plot_best_model_failure_diagnostics(T_test, T_metrics, output_dir)

T_metrics = sortrows(T_metrics, 'RMSE_pct', 'ascend');
best_name = T_metrics.model_name(1);
best_type = T_metrics.model_type(1);
T = select_model_rows(T_test, best_name, best_type);
is_outlier = T.cs_abs_error_pct >= prctile(T.cs_abs_error_pct, 95);

fig = figure('Color', 'w', 'Position', [100 100 1400 900]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
scatter(ax, T.q_error(~is_outlier), T.cs_error_pct(~is_outlier), ...
    16, T.Omega_sr(~is_outlier), 'filled', 'MarkerFaceAlpha', 0.35);
hold(ax, 'on');
scatter(ax, T.q_error(is_outlier), T.cs_error_pct(is_outlier), ...
    45, T.Omega_sr(is_outlier), 'filled', 'MarkerEdgeColor', 'k');
xline(ax, 0, ':k');
yline(ax, 0, ':k');
xlabel(ax, 'q_{pred} - q_{true}');
ylabel(ax, 'Signed SWS error (%)');
title(ax, 'q error propagation');
cb = colorbar(ax);
cb.Label.String = 'Omega (sr)';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T.q_true(~is_outlier), T.q_pred(~is_outlier), ...
    16, T.REQ_M(~is_outlier), 'filled', 'MarkerFaceAlpha', 0.35);
hold(ax, 'on');
scatter(ax, T.q_true(is_outlier), T.q_pred(is_outlier), ...
    45, T.REQ_M(is_outlier), 'filled', 'MarkerEdgeColor', 'k');
plot(ax, [0 1], [0 1], '--k', 'LineWidth', 1.1);
xlim(ax, [0 1]);
ylim(ax, [0 1]);
xlabel(ax, 'q true');
ylabel(ax, 'q predicted');
title(ax, 'Worst 5% SWS errors highlighted');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T.ecum_width_50_rel, T.cs_abs_error_pct, 16, ...
    T.REQ_M, 'filled', 'MarkerFaceAlpha', 0.35);
xlabel(ax, 'Ecum width 50 rel');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Error vs cumulative-curve width');
cb = colorbar(ax);
cb.Label.String = 'REQ M';
grid(ax, 'on');

ax = nexttile(tl);
scatter(ax, T.ecum_increment_peak_frac, T.cs_abs_error_pct, 16, ...
    T.Omega_sr, 'filled', 'MarkerFaceAlpha', 0.35);
xlabel(ax, 'Ecum increment peak fraction');
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Error vs spectral concentration');
cb = colorbar(ax);
cb.Label.String = 'Omega (sr)';
grid(ax, 'on');

title(tl, sprintf('Level 15 failure diagnostics: %s / %s', ...
    best_name, best_type), 'Interpreter', 'none', 'FontWeight', 'bold');
save_level15_figure(fig, output_dir, ...
    'level15_best_model_failure_diagnostics');

end

function plot_failure_correlation_summary(T_corr, output_dir)

if isempty(T_corr)
    return;
end

n_show = min(12, height(T_corr));
T = T_corr(1:n_show, :);
labels = categorical(T.feature_name);
labels = reordercats(labels, cellstr(T.feature_name));

fig = figure('Color', 'w', 'Position', [100 100 950 560]);
ax = axes(fig);
barh(ax, labels, T.spearman_abs_sws_error);
xlabel(ax, 'Spearman rho with absolute SWS error');
ylabel(ax, 'Feature');
title(ax, 'Best-model failure associations');
grid(ax, 'on');
save_level15_figure(fig, output_dir, ...
    'level15_failure_correlation_summary');

end

function plot_high_error_rates(T_by_M, T_by_cs, threshold_pct, output_dir)

fig = figure('Color', 'w', 'Position', [100 100 1400 650]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
T_plot = T_by_M(T_by_M.model_type == "bagged_trees", :);
plot_grouped_rate_bars(ax, T_plot, "REQ_M");
xlabel(ax, 'REQ M');
ylabel(ax, sprintf('Predictions with APE > %.0f%%', threshold_pct));
title(ax, 'High-error rate by M');
grid(ax, 'on');

ax = nexttile(tl);
T_plot = T_by_cs(T_by_cs.model_type == "bagged_trees", :);
plot_grouped_rate_bars(ax, T_plot, "SIM_cs_bg");
xlabel(ax, 'True c_s (m/s)');
ylabel(ax, sprintf('Predictions with APE > %.0f%%', threshold_pct));
title(ax, 'High-error rate by true SWS');
grid(ax, 'on');

title(tl, 'Level 15 high-error diagnostics', 'FontWeight', 'bold');
save_level15_figure(fig, output_dir, 'level15_high_error_rates');

end

function plot_grouped_rate_bars(ax, T, x_var)

if isempty(T)
    return;
end

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

function T = select_model_rows(T_all, model_name, model_type)

T = T_all(T_all.model_name == model_name & ...
    T_all.model_type == model_type, :);

end

function save_level15_figure(fig, output_dir, base_name)

exportgraphics(fig, fullfile(output_dir, base_name + ".png"), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, base_name + ".pdf"), ...
    'ContentType', 'vector');
close(fig);

end
