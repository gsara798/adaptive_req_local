%% analyze_test_16_theory_residual_q.m
% Test 16: theory-informed residual q model.
%
% Residual formulation:
%   delta_q_target = q_theory - q_prior
%   q_pred = q_prior + delta_q_pred
%
% This script uses the existing Test 12 dataset. It does not regenerate data
% and does not retrain the previously deployed hybrid baseline.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLegendFontSize', 11);

%% Project setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

%% Settings

source_test_name = 'test_12_cs_guess_window_sweep';
test16_name = 'test_16_theory_residual_q';
model_types = "bagged_trees";
train_fraction = 0.70;
% Keep this aligned with Test 15 so direct-vs-residual comparisons use the
% same grouped condition split.
random_seed = 1515;
num_learning_cycles = 160;
min_leaf_size = 8;
use_parallel = true;

forbidden_predictors = lower([
    "q_theory"
    "q_true"
    "q_global_theory"
    "q_local_minus_global"
    "abs_q_local_minus_global"
    "M_eff_true_diag"
    "cs_true"
    "cs_pred"
    "sws_error"
    "abs_sws_error"
    "sws_error_pct"
    "abs_sws_error_pct"
    "residual"
    "abs_error"
    "aperture_weight"
    "solid_angle_weight"
    "true_aperture_weight"]);

%% Load Test 12

[T_feat, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    source_test_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

OUT = make_test16_dirs(root_dir);
fprintf('\nTest 16 output folder:\n%s\n', OUT.analysis_dir);

model_file = fullfile(OUT.analysis_dir, 'level16_theory_residual_q_models.mat');
pred_file = fullfile(OUT.table_dir, 'level16_predictions.csv');
metrics_file = fullfile(OUT.table_dir, 'level16_sws_metrics.csv');

if exist(model_file, 'file') == 2 && exist(pred_file, 'file') == 2 && ...
        exist(metrics_file, 'file') == 2
    fprintf('\nExisting Test 16 outputs found. Resuming registry/figures without retraining.\n');
    S16 = load(model_file, 'MODELS');
    MODELS = S16.MODELS;
    T_sws = readtable(pred_file, 'TextType', 'string');
    T_q_metrics = readtable(fullfile(OUT.table_dir, 'level16_q_metrics.csv'), 'TextType', 'string');
    T_sws_metrics = readtable(metrics_file, 'TextType', 'string');
    T_model_comparison = readtable(fullfile(OUT.table_dir, ...
        'level16_model_comparison.csv'), 'TextType', 'string');
    T_sws_by_user = readtable(fullfile(OUT.table_dir, ...
        'level16_sws_metrics_by_user_guess.csv'), 'TextType', 'string');
    T_sws_by_M_eff = readtable(fullfile(OUT.table_dir, ...
        'level16_sws_metrics_by_M_eff.csv'), 'TextType', 'string');
    T_prior_type = readtable(fullfile(OUT.table_dir, ...
        'level16_sws_metrics_by_prior_type.csv'), 'TextType', 'string');
    register_test16_models(root_dir, MODELS, OUT, test16_name, source_test_name, model_types);
    register_test16_baselines(root_dir, OUT, test16_name, source_test_name);
    make_level16_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, ...
        T_prior_type, OUT);
    print_summary(T_sws_metrics, T_model_comparison, OUT);
    return;
end

required_vars = [
    "q_theory"
    "req_mapping"
    "condition_id"
    "step_idx"
    "realization_idx"
    "patch_idx"
    "SIM_dx"
    "SIM_dz"
    "SIM_f0"
    "SIM_cs_bg"
    "REQ_M"
    "REQ_cs_guess"];
adaptive_req.analysis.Test12Analysis.requireVars(T_feat, required_vars, 'Test 16 input');

T_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_feat);
T_feat = adaptive_req.analysis.Test12Analysis.addBins(T_feat);

%% Theory and user-prior features

[T_base, T_theory_cache] = adaptive_req.analysis.build_theory_q_features( ...
    T_feat, 'Verbose', true);
T_base = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_base);
T_base = adaptive_req.analysis.Test12Analysis.addBins(T_base);
T_base.row_instance_id = (1:height(T_base))';
T_base.q_prior = T_base.q_theory_mean_all;
T_base.prior_type = repmat("mean_all_unknown", height(T_base), 1);
T_base.delta_q_target = T_base.q_theory - T_base.q_prior;

T_user = adaptive_req.analysis.build_user_field_guess_features( ...
    T_base, 'ExpandGuesses', true);
T_user = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_user);
T_user = adaptive_req.analysis.Test12Analysis.addBins(T_user);
T_user.row_instance_id = (1:height(T_user))';
T_user.q_prior = T_user.q_user_guess_prior;
T_user.prior_type = "user_" + string(T_user.user_field_guess);
T_user.delta_q_target = T_user.q_theory - T_user.q_prior;

%% Predictor sets

specs12 = adaptive_req.analysis.Test12Analysis.buildModelSpecs(T_base);
idx_hybrid = find([specs12.model_name] == "HybridLocalGlobal" & ...
    [specs12.feature_set] == "NoCsGuess", 1);
assert(~isempty(idx_hybrid), 'Could not find HybridLocalGlobal | NoCsGuess predictors.');

hybrid_predictors = string(specs12(idx_hybrid).predictors(:));
hybrid_predictors = filter_operational_predictors(hybrid_predictors, forbidden_predictors);

theory_predictors = [
    "q_theory_dir2D"
    "q_theory_diffuse2D"
    "q_theory_projected3D"
    "q_theory_mean_dir2D_projected3D"
    "q_theory_mean_all"
    "q_prior"
    "REQ_cs_guess"
    "M_eff_guess"];

user_predictors = [
    "q_user_guess_prior"
    "user_field_guess"];

residual_predictors = adaptive_req.analysis.Test12Analysis.existingPredictors( ...
    T_base, unique([hybrid_predictors; theory_predictors], 'stable'));
residual_predictors = filter_operational_predictors(residual_predictors, forbidden_predictors);
assert_no_forbidden(residual_predictors, forbidden_predictors, 'ResidualTheoryDirect');

user_residual_predictors = adaptive_req.analysis.Test12Analysis.existingPredictors( ...
    T_user, unique([hybrid_predictors; theory_predictors; user_predictors], 'stable'));
user_residual_predictors = filter_operational_predictors(user_residual_predictors, forbidden_predictors);
assert_no_forbidden(user_residual_predictors, forbidden_predictors, ...
    'ResidualTheoryPlusUserGuess');

fprintf('\nOperational predictors confirmed: target/error/oracle variables were not used.\n');

%% Grouped condition split

[train_mask_base, test_mask_base] = adaptive_req.analysis.Test12Analysis.conditionSplit( ...
    T_base, train_fraction, random_seed);
row_repeats = repmat((1:height(T_base))', 4, 1);
train_mask_user = train_mask_base(row_repeats);
test_mask_user = test_mask_base(row_repeats);

%% Train residual models

MODELS = struct([]);
T_pred_all = table();
T_delta_metrics_all = table();
model_counter = 0;

model_counter = model_counter + 1;
[MODELS(model_counter).MODEL, T_delta_pred_i, T_delta_i] = ...
    adaptive_req.analysis.train_q_model_fixed_split( ...
    T_base, residual_predictors, train_mask_base, test_mask_base, ...
    'QVar', 'delta_q_target', ...
    'ModelName', 'ResidualTheoryDirect', ...
    'ModelRole', 'operational', ...
    'ModelTypes', model_types, ...
    'NumLearningCycles', num_learning_cycles, ...
    'MinLeafSize', min_leaf_size, ...
    'UseParallel', use_parallel, ...
    'ClipRange', [-1 1], ...
    'Verbose', true);
T_delta_pred_i = attach_row_instance_id(T_delta_pred_i, T_base);
MODELS(model_counter).predictors = residual_predictors;
MODELS(model_counter).feature_set = "ResidualTheoryDirect";
MODELS(model_counter).prior_description = "q_prior = q_theory_mean_all";
T_delta_i.feature_set = repmat("ResidualTheoryDirect", height(T_delta_i), 1);
T_delta_metrics_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_delta_metrics_all, T_delta_i);
T_pred_i = residual_prediction_table(T_base, T_delta_pred_i(T_delta_pred_i.split == "test", :), ...
    'ResidualTheoryDirect', 'ResidualTheoryDirect', 'bagged_trees', 'operational');
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables(T_pred_all, T_pred_i);

model_counter = model_counter + 1;
[MODELS(model_counter).MODEL, T_delta_pred_i, T_delta_i] = ...
    adaptive_req.analysis.train_q_model_fixed_split( ...
    T_user, user_residual_predictors, train_mask_user, test_mask_user, ...
    'QVar', 'delta_q_target', ...
    'ModelName', 'ResidualTheoryPlusUserGuess', ...
    'ModelRole', 'operational', ...
    'ModelTypes', model_types, ...
    'NumLearningCycles', num_learning_cycles, ...
    'MinLeafSize', min_leaf_size, ...
    'UseParallel', use_parallel, ...
    'ClipRange', [-1 1], ...
    'Verbose', true);
T_delta_pred_i = attach_row_instance_id(T_delta_pred_i, T_user);
MODELS(model_counter).predictors = user_residual_predictors;
MODELS(model_counter).feature_set = "ResidualTheoryPlusUserGuess";
MODELS(model_counter).prior_description = "q_prior = q_user_guess_prior";
T_delta_i.feature_set = repmat("ResidualTheoryPlusUserGuess", height(T_delta_i), 1);
T_delta_metrics_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_delta_metrics_all, T_delta_i);
T_pred_i = residual_prediction_table(T_user, T_delta_pred_i(T_delta_pred_i.split == "test", :), ...
    'ResidualTheoryPlusUserGuess', 'ResidualTheoryPlusUserGuess', ...
    'bagged_trees', 'operational');
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables(T_pred_all, T_pred_i);

%% Existing baseline and no-ML baselines

T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, predict_current_hybrid_baseline(T_base, test_mask_base, PATHS));
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, make_prior_only_baseline(T_user, test_mask_user));
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, make_best_theory_oracle(T_base, test_mask_base));

%% Convert q to SWS and metrics

T_sws = adaptive_req.analysis.Test12Analysis.addSwsMetrics(T_pred_all, T_base);
T_sws = add_theory_metadata(T_sws, T_base);
T_sws = adaptive_req.analysis.Test12Analysis.addBins(T_sws);

T_q_metrics = summarize_q_metrics(T_pred_all);
T_sws_metrics = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role"]);
T_sws_by_user = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role", "user_field_guess"]);
T_sws_by_M_eff = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role", "M_eff_true_diag_bin"]);
T_prior_type = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role", "prior_type"]);
T_model_comparison = make_model_comparison(T_sws_metrics);
T_test15_compare = compare_to_test15(root_dir, T_sws_metrics);

%% Save outputs

writetable(adaptive_req.analysis.Test12Analysis.removeCellColumns(T_sws), pred_file);
writetable(T_q_metrics, fullfile(OUT.table_dir, 'level16_q_metrics.csv'));
writetable(T_sws_metrics, metrics_file);
writetable(T_model_comparison, fullfile(OUT.table_dir, 'level16_model_comparison.csv'));
writetable(T_sws_by_user, fullfile(OUT.table_dir, 'level16_sws_metrics_by_user_guess.csv'));
writetable(T_sws_by_M_eff, fullfile(OUT.table_dir, 'level16_sws_metrics_by_M_eff.csv'));
writetable(T_prior_type, fullfile(OUT.table_dir, 'level16_sws_metrics_by_prior_type.csv'));
writetable(T_delta_metrics_all, fullfile(OUT.table_dir, 'level16_delta_q_training_metrics.csv'));
writetable(T_theory_cache, fullfile(OUT.table_dir, 'level16_theory_q_cache.csv'));
if ~isempty(T_test15_compare)
    writetable(T_test15_compare, fullfile(OUT.table_dir, 'level16_vs_test15_summary.csv'));
end

save(model_file, 'MODELS', 'T_q_metrics', 'T_sws_metrics', ...
    'T_model_comparison', 'T_sws_by_user', 'T_sws_by_M_eff', ...
    'T_prior_type', 'T_delta_metrics_all', 'T_test15_compare', ...
    'T_theory_cache', 'MC', 'PATHS', '-v7.3');

%% Register models

register_test16_models(root_dir, MODELS, OUT, test16_name, source_test_name, model_types);
register_test16_baselines(root_dir, OUT, test16_name, source_test_name);

%% Figures

make_level16_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, ...
    T_prior_type, OUT);

%% Console summary

print_summary(T_sws_metrics, T_model_comparison, OUT);
fprintf('\nWarning: BestTheoryOracle is diagnostic_only and must not be used as an operational model.\n');
fprintf('Confirmed: q_theory, M_eff_true_diag, aperture_weight and error variables were not used as operational predictors.\n');
fprintf('Model registry manifest:\n%s\n', fullfile(root_dir, 'outputs', 'model_registry', 'model_manifest.csv'));
fprintf('Test 16 complete:\n%s\n', OUT.analysis_dir);

%% Local functions

function OUT = make_test16_dirs(root_dir)
OUT.analysis_dir = fullfile(root_dir, 'outputs', ...
    'test_16_theory_residual_q', 'analysis');
OUT.fig_dir = fullfile(OUT.analysis_dir, 'figures');
OUT.table_dir = fullfile(OUT.analysis_dir, 'tables');
OUT.model_dir = fullfile(OUT.analysis_dir, 'models');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if ~exist(dirs(i), 'dir')
        mkdir(dirs(i));
    end
end
end

function predictors = filter_operational_predictors(predictors, forbidden)
predictors = string(predictors(:));
lower_names = lower(predictors);
bad = ismember(lower_names, forbidden) | contains(lower_names, "error") | ...
    contains(lower_names, "residual") | contains(lower_names, "target");
predictors = predictors(~bad);
end

function assert_no_forbidden(predictors, forbidden, name)
lower_names = lower(string(predictors(:)));
bad = intersect(lower_names, forbidden);
assert(isempty(bad), '%s contains forbidden predictors: %s', name, strjoin(bad, ', '));
assert(~any(contains(lower_names, "aperture_weight")), ...
    '%s contains aperture_weight-like predictors.', name);
end

function T_pred = attach_row_instance_id(T_pred, T_source)
if height(T_pred) ~= height(T_source)
    error(['Cannot attach row_instance_id because prediction/source heights differ. ', ...
        'Prediction rows: %d, source rows: %d.'], height(T_pred), height(T_source));
end
T_pred.row_instance_id = T_source.row_instance_id;
end

function T = normalize_output_types(T)
vars = string(T.Properties.VariableNames);
if ismember('user_field_guess', vars)
    T.user_field_guess = string(T.user_field_guess);
end
if ismember('prior_type', vars)
    T.prior_type = string(T.prior_type);
end
end

function T = residual_prediction_table(T_ref, T_delta_pred, model_name, feature_set, model_type, model_role)
if ismember('row_instance_id', string(T_delta_pred.Properties.VariableNames)) && ...
        ismember('row_instance_id', string(T_ref.Properties.VariableNames))
    loc = T_delta_pred.row_instance_id;
    assert(all(loc >= 1 & loc <= height(T_ref)), 'Invalid row_instance_id values.');
else
    [tf, loc] = ismember(string(T_delta_pred.row_key), string(T_ref.row_key));
    assert(all(tf), 'Could not map residual predictions back to reference rows.');
end
R = T_ref(loc, :);
keep = ["condition_id", "step_idx", "realization_idx", "patch_idx", ...
    "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M", "REQ_cs_guess", ...
    "M_eff_guess", "M_eff_true_diag", "row_key", "user_field_guess", ...
    "q_user_guess_prior", "q_prior", "prior_type"];
keep = keep(ismember(keep, string(R.Properties.VariableNames)));
T = R(:, cellstr(keep));
T.model_name = repmat(string(model_name), height(T), 1);
T.feature_set = repmat(string(feature_set), height(T), 1);
T.model_type = repmat(string(model_type), height(T), 1);
T.model_role = repmat(string(model_role), height(T), 1);
T.split = repmat("test", height(T), 1);
T.q_true = R.q_theory;
T.delta_q_true = R.delta_q_target;
T.delta_q_pred_raw = T_delta_pred.q_pred_raw;
T.delta_q_pred = T_delta_pred.q_pred;
T.q_pred_raw = R.q_prior + T.delta_q_pred_raw;
T.q_pred = min(max(R.q_prior + T.delta_q_pred, 0.001), 0.999);
T.residual = T.q_pred - T.q_true;
T.abs_error = abs(T.residual);
T = normalize_output_types(T);
end

function T = predict_current_hybrid_baseline(T_base, test_mask, PATHS)
deployment_dir = fullfile(PATHS.analysis_dir, ...
    'level_01_model_comparison', 'models', 'deployment');
try
    [MODEL, INFO] = adaptive_req.analysis.load_q_model_deployment( ...
        deployment_dir, ...
        'ModelName', 'HybridLocalGlobal', ...
        'FeatureSet', 'NoCsGuess', ...
        'ModelType', 'bagged_trees');
catch ME
    warning('Could not load CurrentHybridBaseline: %s', getReport(ME, 'basic'));
    T = table();
    return;
end
T_ref = T_base(test_mask, :);
Tq = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL, T_ref, 'ModelType', INFO.model_type, 'ModelName', 'CurrentHybridBaseline');
T = base_prediction_table(T_ref, Tq.q_pred_raw, min(max(Tq.q_pred, 0.001), 0.999), ...
    'CurrentHybridBaseline', 'CurrentHybridBaseline', 'bagged_trees', ...
    'existing_model_reference', nan(height(T_ref), 1), nan(height(T_ref), 1), ...
    "none", "current_hybrid");
end

function T = make_prior_only_baseline(T_user, test_mask)
T_ref = T_user(test_mask, :);
T = base_prediction_table(T_ref, T_ref.q_prior, T_ref.q_prior, ...
    'PriorOnly', 'PriorOnly', 'baseline_no_ml', 'baseline_no_ml', ...
    T_ref.delta_q_target, zeros(height(T_ref), 1), T_ref.user_field_guess, T_ref.prior_type);
end

function T = make_best_theory_oracle(T_base, test_mask)
T_ref = T_base(test_mask, :);
cands = [T_ref.q_theory_dir2D, T_ref.q_theory_diffuse2D, ...
    T_ref.q_theory_projected3D, T_ref.q_theory_mean_dir2D_projected3D, ...
    T_ref.q_theory_mean_all];
[~, idx] = min(abs(cands - T_ref.q_theory), [], 2, 'omitnan');
q_pred = nan(height(T_ref), 1);
for i = 1:height(T_ref)
    q_pred(i) = cands(i, idx(i));
end
T = base_prediction_table(T_ref, q_pred, q_pred, ...
    'BestTheoryOracle', 'BestTheoryOracle', 'oracle_no_ml', ...
    'diagnostic_only', T_ref.q_theory - q_pred, zeros(height(T_ref), 1), ...
    "oracle", "best_theory_candidate");
end

function T = base_prediction_table(T_ref, q_raw, q_pred, model_name, feature_set, ...
    model_type, model_role, delta_true, delta_pred, user_guess, prior_type)
keep = ["condition_id", "step_idx", "realization_idx", "patch_idx", ...
    "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M", "REQ_cs_guess", ...
    "M_eff_guess", "M_eff_true_diag", "row_key", "user_field_guess", ...
    "q_user_guess_prior", "q_prior", "prior_type"];
keep = keep(ismember(keep, string(T_ref.Properties.VariableNames)));
T = T_ref(:, cellstr(keep));
T.model_name = repmat(string(model_name), height(T), 1);
T.feature_set = repmat(string(feature_set), height(T), 1);
T.model_type = repmat(string(model_type), height(T), 1);
T.model_role = repmat(string(model_role), height(T), 1);
T.split = repmat("test", height(T), 1);
T.q_true = T_ref.q_theory;
T.q_pred_raw = q_raw;
T.q_pred = min(max(q_pred, 0.001), 0.999);
T.delta_q_true = delta_true;
T.delta_q_pred_raw = delta_pred;
T.delta_q_pred = delta_pred;
if ~ismember('user_field_guess', string(T.Properties.VariableNames))
    T.user_field_guess = repmat(string(user_guess), height(T), 1);
end
if ~ismember('prior_type', string(T.Properties.VariableNames))
    T.prior_type = repmat(string(prior_type), height(T), 1);
end
T.residual = T.q_pred - T.q_true;
T.abs_error = abs(T.residual);
T = normalize_output_types(T);
end

function T = add_theory_metadata(T, T_base)
[tf, loc] = ismember(string(T.row_key), string(T_base.row_key));
assert(all(tf), 'Could not map theory metadata.');
names = ["q_theory_dir2D", "q_theory_diffuse2D", "q_theory_projected3D", ...
    "q_theory_mean_dir2D_projected3D", "q_theory_mean_all"];
for i = 1:numel(names)
    if ~ismember(names(i), string(T.Properties.VariableNames))
        T.(char(names(i))) = T_base.(char(names(i)))(loc);
    end
end
if ~ismember('q_prior', string(T.Properties.VariableNames))
    T.q_prior = nan(height(T), 1);
end
if ~ismember('prior_type', string(T.Properties.VariableNames))
    T.prior_type = repmat("not_used", height(T), 1);
end
T = normalize_output_types(T);
end

function Tq = summarize_q_metrics(T)
groups = ["model_name", "feature_set", "model_type", "model_role", "split"];
[G, Tq] = findgroups(T(:, cellstr(groups)));
Tq.N = splitapply(@numel, T.q_true, G);
Tq.MAE_q = splitapply(@(a,b) mean(abs(b-a), 'omitnan'), T.q_true, T.q_pred, G);
Tq.RMSE_q = splitapply(@(a,b) sqrt(mean((b-a).^2, 'omitnan')), T.q_true, T.q_pred, G);
Tq.bias_q = splitapply(@(a,b) mean(b-a, 'omitnan'), T.q_true, T.q_pred, G);
Tq.Pearson_q = splitapply(@(a,b) safe_corr(a,b,'Pearson'), T.q_true, T.q_pred, G);
Tq.Spearman_q = splitapply(@(a,b) safe_corr(a,b,'Spearman'), T.q_true, T.q_pred, G);
end

function r = safe_corr(x, y, type)
valid = isfinite(x) & isfinite(y);
if nnz(valid) < 3 || std(x(valid)) <= eps || std(y(valid)) <= eps
    r = NaN;
else
    r = corr(x(valid), y(valid), 'Type', type);
end
end

function T = make_model_comparison(T_sws)
T = sortrows(T_sws, 'MAPE_pct', 'ascend');
if any(T.model_name == "CurrentHybridBaseline")
    base = T(T.model_name == "CurrentHybridBaseline", :);
    for i = 1:height(T)
        T.Delta_MAPE_vs_CurrentHybridBaseline(i, 1) = T.MAPE_pct(i) - base.MAPE_pct(1);
        T.Delta_HighError_gt20_vs_CurrentHybridBaseline(i, 1) = ...
            T.HighError_gt20_pct(i) - base.HighError_gt20_pct(1);
    end
end
end

function T = compare_to_test15(root_dir, T16)
test15_file = fullfile(root_dir, 'outputs', 'test_15_theory_informed_direct_q', ...
    'analysis', 'tables', 'level15_sws_metrics.csv');
if exist(test15_file, 'file') ~= 2
    T = table();
    warning('Test 15 metrics not found. Skipping Test 15 comparison.');
    return;
end
T15 = readtable(test15_file, 'TextType', 'string');
T16s = T16(:, {'model_name','feature_set','model_type','model_role','MAPE_pct', ...
    'HighError_gt20_pct','P95_abs_error_pct'});
T16s.source = repmat("Test16Residual", height(T16s), 1);
T15s = T15(:, {'model_name','feature_set','model_type','model_role','MAPE_pct', ...
    'HighError_gt20_pct','P95_abs_error_pct'});
T15s.source = repmat("Test15Direct", height(T15s), 1);
T = adaptive_req.analysis.Test12Analysis.concatTables(T16s, T15s);
end

function make_level16_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, T_prior_type, OUT)
plot_metric_bar(T_sws_metrics, 'MAPE_pct', ...
    fullfile(OUT.fig_dir, 'level16_mape_by_model.png'), ...
    'Test 16 MAPE by model');
plot_metric_bar(T_sws_metrics, 'HighError_gt20_pct', ...
    fullfile(OUT.fig_dir, 'level16_high_error_by_model.png'), ...
    'Test 16 high-error >20% by model');
plot_delta_vs_baseline(T_sws_metrics, ...
    fullfile(OUT.fig_dir, 'level16_delta_mape_vs_CurrentHybridBaseline.png'));
plot_q_scatter(T_sws, 'q_true', 'q_pred', ...
    fullfile(OUT.fig_dir, 'level16_q_true_vs_pred.png'), ...
    'Test 16 q true vs predicted');
plot_q_scatter(T_sws, 'delta_q_true', 'delta_q_pred', ...
    fullfile(OUT.fig_dir, 'level16_delta_q_true_vs_pred.png'), ...
    'Test 16 delta q true vs predicted');
plot_group_metric(T_sws_by_user, "user_field_guess", "model_name", "MAPE_pct", ...
    'Test 16 MAPE by user field guess', ...
    fullfile(OUT.fig_dir, 'level16_mape_by_user_field_guess.png'));
plot_group_metric(T_sws_by_M_eff, "M_eff_true_diag_bin", "model_name", "MAPE_pct", ...
    'Test 16 MAPE by M_eff_true_diag diagnostic bin', ...
    fullfile(OUT.fig_dir, 'level16_mape_by_M_eff_true_diag.png'));
plot_group_metric(T_prior_type, "prior_type", "model_name", "MAPE_pct", ...
    'Test 16 residual error by prior type', ...
    fullfile(OUT.fig_dir, 'level16_residual_error_by_prior_type.png'));
end

function plot_metric_bar(T, metric_var, file_path, title_text)
T = T(T.model_type == "bagged_trees" | T.model_type == "baseline_no_ml" | ...
    T.model_type == "oracle_no_ml", :);
T = sortrows(T, metric_var, 'ascend');
figure('Color', 'w', 'Position', [100 100 980 560]);
barh(categorical(short_model_labels(T)), T.(metric_var), 0.68);
xlabel(metric_label(metric_var));
title(title_text, 'Interpreter', 'none', 'FontWeight', 'normal');
grid on;
set(gca, 'YDir', 'reverse', 'FontSize', 11);
export_clean_figure(gcf, file_path);
close(gcf);
end

function plot_delta_vs_baseline(T, file_path)
Tc = make_model_comparison(T);
if ~ismember('Delta_MAPE_vs_CurrentHybridBaseline', string(Tc.Properties.VariableNames))
    return;
end
figure('Color', 'w', 'Position', [100 100 1050 520]);
Tc = sortrows(Tc, 'Delta_MAPE_vs_CurrentHybridBaseline', 'ascend');
barh(categorical(short_model_labels(Tc)), Tc.Delta_MAPE_vs_CurrentHybridBaseline);
xline(0, 'k--');
xlabel('\Delta MAPE vs current hybrid baseline (%)');
title('Negative delta means improvement over current hybrid baseline', ...
    'Interpreter', 'none', 'FontWeight', 'normal');
grid on;
set(gca, 'YDir', 'reverse', 'FontSize', 11);
export_clean_figure(gcf, file_path);
close(gcf);
end

function plot_q_scatter(T, x_var, y_var, file_path, title_text)
T = T(T.model_type == "bagged_trees" | T.model_type == "baseline_no_ml", :);
models = unique(T.model_name, 'stable');
figure('Color', 'w', 'Position', [100 100 1250 760]);
tl = tiledlayout(2, ceil(numel(models)/2), 'TileSpacing', 'compact', 'Padding', 'compact');
rng(1616);
for i = 1:numel(models)
    ax = nexttile(tl);
    idx = T.model_name == models(i) & isfinite(T.(char(x_var))) & isfinite(T.(char(y_var)));
    ii = find(idx);
    if isempty(ii)
        title(ax, short_model_name(models(i)), 'Interpreter', 'none', 'FontWeight', 'normal');
        text(ax, 0.5, 0.5, 'No finite values', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        continue;
    end
    if numel(ii) > 25000
        ii = ii(randperm(numel(ii), 25000));
    end
    scatter(ax, T.(char(x_var))(ii), T.(char(y_var))(ii), 5, ...
        'filled', 'MarkerFaceAlpha', 0.12);
    hold(ax, 'on');
    vals = [T.(char(x_var))(ii); T.(char(y_var))(ii)];
    lo = min(vals, [], 'omitnan');
    hi = max(vals, [], 'omitnan');
    if ~isscalar(lo) || ~isscalar(hi) || ~isfinite(lo) || ~isfinite(hi) || lo == hi
        lo = 0; hi = 1;
    end
    plot(ax, [lo hi], [lo hi], 'k--');
    axis(ax, 'equal'); xlim(ax, [lo hi]); ylim(ax, [lo hi]); grid(ax, 'on');
    title(ax, short_model_name(models(i)), 'Interpreter', 'none', 'FontWeight', 'normal');
    xlabel(ax, strrep(char(x_var), '_', '\_'));
    ylabel(ax, strrep(char(y_var), '_', '\_'));
end
title(tl, title_text, 'Interpreter', 'none', 'FontWeight', 'normal');
export_clean_figure(gcf, file_path);
close(gcf);
end

function plot_group_metric(T, x_var, series_var, metric_var, title_text, file_path)
if isempty(T)
    return;
end
if ismember('model_role', string(T.Properties.VariableNames))
    T = T(T.model_role ~= "diagnostic_only", :);
end
x = string(T.(char(x_var)));
series = string(T.(char(series_var)));
x_values = unique(x, 'stable');
series_values = unique(series, 'stable');
Y = nan(numel(x_values), numel(series_values));
for i = 1:numel(x_values)
    for j = 1:numel(series_values)
        idx = x == x_values(i) & series == series_values(j);
        if any(idx)
            Y(i, j) = mean(T.(char(metric_var))(idx), 'omitnan');
        end
    end
end
figure('Color', 'w', 'Position', [100 100 1120 560]);
bar(categorical(x_values), Y);
ylabel(metric_label(metric_var));
xlabel(strrep(char(x_var), '_', '\_'));
title(title_text, 'Interpreter', 'none', 'FontWeight', 'normal');
legend(arrayfun(@short_model_name, series_values), ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'none');
xtickangle(25);
grid on;
set(gca, 'FontSize', 10);
export_clean_figure(gcf, file_path);
close(gcf);
end

function labels = short_model_labels(T)
labels = strings(height(T), 1);
for i = 1:height(T)
    labels(i) = short_model_name(T.model_name(i));
end
end

function out = short_model_name(name)
name = string(name);
switch name
    case "CurrentHybridBaseline"
        out = "Current hybrid";
    case "ResidualTheoryDirect"
        out = "Residual theory";
    case "ResidualTheoryPlusUserGuess"
        out = "Residual + user guess";
    case "PriorOnly"
        out = "Prior only";
    case "BestTheoryOracle"
        out = "Best theory oracle";
    otherwise
        out = name;
end
end

function out = metric_label(metric_var)
switch string(metric_var)
    case "MAPE_pct"
        out = 'MAPE (%)';
    case "HighError_gt20_pct"
        out = 'High-error >20% (%)';
    otherwise
        out = strrep(char(metric_var), '_', '\_');
end
end

function export_clean_figure(fig, file_path)
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end
drawnow;
exportgraphics(fig, file_path, 'Resolution', 300, 'BackgroundColor', 'white');
end

function register_test16_models(root_dir, MODELS, OUT, test_name, source_test_name, model_types)
for i = 1:numel(MODELS)
    model_i = MODELS(i).MODEL;
    adaptive_req.analysis.register_trained_model( ...
        'RootDir', root_dir, ...
        'ModelObject', model_i, ...
        'ModelId', "test16__" + string(model_i.model_name) + "__" + string(model_types), ...
        'RegistrySubdir', 'test16_theory_residual_q', ...
        'TestName', test_name, ...
        'AnalysisLevel', 'analysis', ...
        'ModelName', model_i.model_name, ...
        'FeatureSet', MODELS(i).feature_set, ...
        'ModelType', model_types, ...
        'ModelRole', model_i.model_role, ...
        'TrainingDataset', source_test_name, ...
        'Target', 'delta_q_target', ...
        'PredictorSummary', strjoin(string(MODELS(i).predictors), ', '), ...
        'MetricsFile', fullfile(OUT.table_dir, 'level16_sws_metrics.csv'), ...
        'Notes', 'Theory-informed residual q model: q_pred = q_prior + delta_q_pred.');
end
end

function register_test16_baselines(root_dir, OUT, test_name, source_test_name)
adaptive_req.analysis.register_trained_model( ...
    'RootDir', root_dir, ...
    'ModelId', 'test16__PriorOnly__baseline_no_ml', ...
    'RegistrySubdir', 'test16_theory_residual_q', ...
    'TestName', test_name, ...
    'AnalysisLevel', 'analysis', ...
    'ModelName', 'PriorOnly', ...
    'FeatureSet', 'PriorOnly', ...
    'ModelType', 'baseline_no_ml', ...
    'ModelRole', 'baseline_no_ml', ...
    'TrainingDataset', source_test_name, ...
    'Target', 'q_theory', ...
    'PredictorSummary', 'q_pred = q_prior', ...
    'MetricsFile', fullfile(OUT.table_dir, 'level16_sws_metrics.csv'), ...
    'Notes', 'No-ML theory/user prior baseline for Test 16.');
end

function print_summary(T_sws_metrics, T_model_comparison, OUT)
fprintf('\nTest 16 SWS metrics:\n');
disp(T_sws_metrics(:, {'model_name','model_type','model_role','MAPE_pct', ...
    'MedAE_pct','bias_pct','HighError_gt10_pct','HighError_gt20_pct'}));

base = T_sws_metrics(T_sws_metrics.model_name == "CurrentHybridBaseline", :);
if ~isempty(base)
    fprintf('\nCurrentHybridBaseline MAPE: %.4f%%\n', base.MAPE_pct(1));
end
residuals = T_sws_metrics(ismember(T_sws_metrics.model_name, ...
    ["ResidualTheoryDirect", "ResidualTheoryPlusUserGuess"]), :);
for i = 1:height(residuals)
    delta = NaN;
    if ~isempty(base)
        delta = residuals.MAPE_pct(i) - base.MAPE_pct(1);
    end
    fprintf('%s MAPE: %.4f%% | Delta vs current hybrid: %.4f%% | HighError>20: %.4f%%\n', ...
        residuals.model_name(i), residuals.MAPE_pct(i), delta, ...
        residuals.HighError_gt20_pct(i));
end

test15_cmp = fullfile(OUT.table_dir, 'level16_vs_test15_summary.csv');
if exist(test15_cmp, 'file') == 2
    fprintf('\nTest 16 vs Test 15 summary table:\n%s\n', test15_cmp);
end
fprintf('\nModel comparison table:\n%s\n', ...
    fullfile(OUT.table_dir, 'level16_model_comparison.csv'));
disp(T_model_comparison(:, {'model_name','MAPE_pct','HighError_gt20_pct', ...
    'Delta_MAPE_vs_CurrentHybridBaseline'}));
end
