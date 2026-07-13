%% analyze_test_15_theory_informed_direct_q.m
% Test 15: theory-informed direct q model.
%
% Direct model, not residual:
%   q_pred = F(spectral_features, theory_q_candidates, optional user prior)
%
% This analysis uses Test 12 as the training dataset, does not modify Test 11
% or Test 12, and does not retrain the existing Local/Global/Hybrid models.

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

test_name = 'test_12_cs_guess_window_sweep';
level_name = 'level_15_theory_informed_direct_q';
model_types = "bagged_trees";
train_fraction = 0.70;
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
    test_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

OUT = adaptive_req.analysis.Test12Analysis.makeOutputDirs(PATHS, level_name);
OUT15 = make_standalone_test15_dirs(root_dir);
fprintf('\nTest 15 output folder:\n%s\n', OUT.analysis_dir);
fprintf('Standalone Test 15 mirror:\n%s\n', OUT15.analysis_dir);

model_file = fullfile(OUT.model_dir, 'level15_theory_informed_direct_q_models.mat');
pred_file = fullfile(OUT.table_dir, 'level15_predictions.csv');
metrics_file = fullfile(OUT.table_dir, 'level15_sws_metrics.csv');
if exist(model_file, 'file') == 2 && exist(pred_file, 'file') == 2 && ...
        exist(metrics_file, 'file') == 2
    fprintf('\nExisting Test 15 outputs found. Resuming registry/figures without retraining.\n');
    S15 = load(model_file, 'MODELS');
    MODELS = S15.MODELS;
    T_sws = readtable(pred_file, 'TextType', 'string');
    T_sws_metrics = readtable(metrics_file, 'TextType', 'string');
    T_sws_by_user = readtable(fullfile(OUT.table_dir, ...
        'level15_sws_metrics_by_user_guess.csv'), 'TextType', 'string');
    T_sws_by_M_eff = readtable(fullfile(OUT.table_dir, ...
        'level15_sws_metrics_by_M_eff.csv'), 'TextType', 'string');

    register_existing_models(root_dir, PATHS);
    register_test15_models(root_dir, MODELS, OUT, level_name, test_name, model_types);
    register_test15_baselines(root_dir, OUT, level_name, test_name);
    make_level15_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, OUT);
    sync_test15_outputs(OUT, OUT15);
    print_summary(T_sws_metrics, T_sws_by_user, T_sws_by_M_eff);
    fprintf('\nConfirmed: aperture_weight / solid_angle_weight / true_aperture_weight were not used as operational predictors.\n');
    fprintf('Model registry manifest:\n%s\n', fullfile(root_dir, 'outputs', 'model_registry', 'model_manifest.csv'));
    fprintf('Test 15 complete:\n%s\n', OUT.analysis_dir);
    fprintf('Standalone Test 15 mirror:\n%s\n', OUT15.analysis_dir);
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
adaptive_req.analysis.Test12Analysis.requireVars(T_feat, required_vars, 'Test 15 input');

T_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_feat);
T_feat = adaptive_req.analysis.Test12Analysis.addBins(T_feat);

%% Theory and user-prior features

[T_base, T_theory_cache] = adaptive_req.analysis.build_theory_q_features( ...
    T_feat, 'Verbose', true);
T_base = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_base);
T_base = adaptive_req.analysis.Test12Analysis.addBins(T_base);

T_user = adaptive_req.analysis.build_user_field_guess_features( ...
    T_base, 'ExpandGuesses', true);
T_user = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_user);
T_user = adaptive_req.analysis.Test12Analysis.addBins(T_user);

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
    "q_theory_mean_all"];

user_predictors = [
    "q_user_guess_prior"
    "user_field_guess"];

direct_predictors = adaptive_req.analysis.Test12Analysis.existingPredictors( ...
    T_base, unique([hybrid_predictors; theory_predictors], 'stable'));
direct_predictors = filter_operational_predictors(direct_predictors, forbidden_predictors);
assert_no_forbidden(direct_predictors, forbidden_predictors, 'TheoryCandidatesDirect');

user_direct_predictors = adaptive_req.analysis.Test12Analysis.existingPredictors( ...
    T_user, unique([hybrid_predictors; theory_predictors; user_predictors], 'stable'));
user_direct_predictors = filter_operational_predictors(user_direct_predictors, forbidden_predictors);
assert_no_forbidden(user_direct_predictors, forbidden_predictors, ...
    'TheoryCandidatesPlusUserGuessDirect');

fprintf('\nOperational predictors confirmed: aperture_weight was not used.\n');

%% Grouped condition split

[train_mask_base, test_mask_base] = adaptive_req.analysis.Test12Analysis.conditionSplit( ...
    T_base, train_fraction, random_seed);
train_mask_user = train_mask_base(repmat((1:height(T_base))', 4, 1));
test_mask_user = test_mask_base(repmat((1:height(T_base))', 4, 1));

%% Train new direct models

MODELS = struct([]);
T_pred_all = table();
T_q_metrics_all = table();
model_counter = 0;

model_counter = model_counter + 1;
[MODELS(model_counter).MODEL, T_pred_i, T_q_i] = ...
    adaptive_req.analysis.train_q_model_fixed_split( ...
    T_base, direct_predictors, train_mask_base, test_mask_base, ...
    'ModelName', 'TheoryCandidatesDirect', ...
    'ModelRole', 'operational', ...
    'ModelTypes', model_types, ...
    'NumLearningCycles', num_learning_cycles, ...
    'MinLeafSize', min_leaf_size, ...
    'UseParallel', use_parallel, ...
    'Verbose', true);
MODELS(model_counter).predictors = direct_predictors;
MODELS(model_counter).feature_set = "TheoryCandidatesDirect";
T_pred_i.feature_set = repmat("TheoryCandidatesDirect", height(T_pred_i), 1);
T_q_i.feature_set = repmat("TheoryCandidatesDirect", height(T_q_i), 1);
T_pred_i = normalize_output_types(T_pred_i);
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables(T_pred_all, T_pred_i(T_pred_i.split == "test", :));
T_q_metrics_all = adaptive_req.analysis.Test12Analysis.concatTables(T_q_metrics_all, T_q_i);

model_counter = model_counter + 1;
[MODELS(model_counter).MODEL, T_pred_i, T_q_i] = ...
    adaptive_req.analysis.train_q_model_fixed_split( ...
    T_user, user_direct_predictors, train_mask_user, test_mask_user, ...
    'ModelName', 'TheoryCandidatesPlusUserGuessDirect', ...
    'ModelRole', 'operational', ...
    'ModelTypes', model_types, ...
    'NumLearningCycles', num_learning_cycles, ...
    'MinLeafSize', min_leaf_size, ...
    'UseParallel', use_parallel, ...
    'Verbose', true);
MODELS(model_counter).predictors = user_direct_predictors;
MODELS(model_counter).feature_set = "TheoryCandidatesPlusUserGuessDirect";
T_pred_i.feature_set = repmat("TheoryCandidatesPlusUserGuessDirect", height(T_pred_i), 1);
T_q_i.feature_set = repmat("TheoryCandidatesPlusUserGuessDirect", height(T_q_i), 1);
T_pred_i = normalize_output_types(T_pred_i);
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables(T_pred_all, T_pred_i(T_pred_i.split == "test", :));
T_q_metrics_all = adaptive_req.analysis.Test12Analysis.concatTables(T_q_metrics_all, T_q_i);

%% Existing baseline and no-ML baselines

T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, predict_current_hybrid_baseline(T_base, test_mask_base, PATHS));
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, make_user_guess_prior_baseline(T_user, test_mask_user));
T_pred_all = adaptive_req.analysis.Test12Analysis.concatTables( ...
    T_pred_all, make_theory_best_candidate_oracle(T_user, test_mask_user));

%% Convert q to SWS and metrics

T_sws = adaptive_req.analysis.Test12Analysis.addSwsMetrics(T_pred_all, T_base);
T_sws = add_test15_metadata(T_sws, T_base, T_user);
T_sws = adaptive_req.analysis.Test12Analysis.addBins(T_sws);

T_q_metrics = summarize_q_metrics(T_pred_all);
T_sws_metrics = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role"]);
T_sws_by_user = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role", "user_field_guess"]);
T_sws_by_M_eff = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws, ["model_name", "feature_set", "model_type", "model_role", "M_eff_true_diag_bin"]);
T_model_comparison = make_model_comparison(T_sws_metrics);

%% Save outputs

writetable(adaptive_req.analysis.Test12Analysis.removeCellColumns(T_sws), ...
    fullfile(OUT.table_dir, 'level15_predictions.csv'));
writetable(T_q_metrics, fullfile(OUT.table_dir, 'level15_q_metrics.csv'));
writetable(T_sws_metrics, fullfile(OUT.table_dir, 'level15_sws_metrics.csv'));
writetable(T_sws_by_user, fullfile(OUT.table_dir, 'level15_sws_metrics_by_user_guess.csv'));
writetable(T_sws_by_M_eff, fullfile(OUT.table_dir, 'level15_sws_metrics_by_M_eff.csv'));
writetable(T_model_comparison, fullfile(OUT.table_dir, 'level15_model_comparison.csv'));
writetable(T_theory_cache, fullfile(OUT.table_dir, 'level15_theory_q_cache.csv'));

save(model_file, 'MODELS', 'T_q_metrics', 'T_sws_metrics', ...
    'T_sws_by_user', 'T_sws_by_M_eff', 'T_model_comparison', ...
    'T_theory_cache', 'MC', 'PATHS', '-v7.3');

%% Register models

register_existing_models(root_dir, PATHS);
register_test15_models(root_dir, MODELS, OUT, level_name, test_name, model_types);
register_test15_baselines(root_dir, OUT, level_name, test_name);

%% Figures

make_level15_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, OUT);
sync_test15_outputs(OUT, OUT15);

%% Console summary

print_summary(T_sws_metrics, T_sws_by_user, T_sws_by_M_eff);
fprintf('\nConfirmed: aperture_weight / solid_angle_weight / true_aperture_weight were not used as operational predictors.\n');
fprintf('Model registry manifest:\n%s\n', fullfile(root_dir, 'outputs', 'model_registry', 'model_manifest.csv'));
fprintf('Test 15 complete:\n%s\n', OUT.analysis_dir);
fprintf('Standalone Test 15 mirror:\n%s\n', OUT15.analysis_dir);

%% Local functions

function predictors = filter_operational_predictors(predictors, forbidden)
predictors = string(predictors(:));
lower_names = lower(predictors);
bad = ismember(lower_names, forbidden) | contains(lower_names, "error") | ...
    contains(lower_names, "residual") | contains(lower_names, "target");
predictors = predictors(~bad);
end

function T = normalize_output_types(T)
if ismember('user_field_guess', string(T.Properties.VariableNames))
    T.user_field_guess = string(T.user_field_guess);
end
end

function OUT15 = make_standalone_test15_dirs(root_dir)
OUT15.analysis_dir = fullfile(root_dir, 'outputs', ...
    'test_15_theory_informed_direct_q', 'analysis');
OUT15.fig_dir = fullfile(OUT15.analysis_dir, 'figures');
OUT15.table_dir = fullfile(OUT15.analysis_dir, 'tables');
OUT15.model_dir = fullfile(OUT15.analysis_dir, 'models');
dirs = string(struct2cell(OUT15));
for i = 1:numel(dirs)
    if ~exist(dirs(i), 'dir')
        mkdir(dirs(i));
    end
end
end

function sync_test15_outputs(OUT, OUT15)
copy_dir_contents(OUT.fig_dir, OUT15.fig_dir);
copy_dir_contents(OUT.table_dir, OUT15.table_dir);
copy_dir_contents(OUT.model_dir, OUT15.model_dir);
end

function copy_dir_contents(src_dir, dst_dir)
if ~exist(src_dir, 'dir')
    return;
end
if ~exist(dst_dir, 'dir')
    mkdir(dst_dir);
end
files = dir(src_dir);
for i = 1:numel(files)
    if files(i).isdir || startsWith(files(i).name, ".")
        continue;
    end
    copyfile(fullfile(files(i).folder, files(i).name), ...
        fullfile(dst_dir, files(i).name));
end
end

function assert_no_forbidden(predictors, forbidden, name)
lower_names = lower(string(predictors(:)));
bad = intersect(lower_names, forbidden);
assert(isempty(bad), '%s contains forbidden predictors: %s', name, strjoin(bad, ', '));
assert(~any(contains(lower_names, "aperture_weight")), ...
    '%s contains aperture_weight-like predictors.', name);
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
    'existing_model_reference');
end

function T = make_user_guess_prior_baseline(T_user, test_mask)
T_ref = T_user(test_mask, :);
T = base_prediction_table(T_ref, T_ref.q_user_guess_prior, T_ref.q_user_guess_prior, ...
    'UserGuessPriorOnly', 'UserGuessPriorOnly', 'baseline_no_ml', 'baseline_no_ml');
end

function T = make_theory_best_candidate_oracle(T_user, test_mask)
T_ref = T_user(test_mask, :);
cands = [T_ref.q_theory_dir2D, T_ref.q_theory_diffuse2D, ...
    T_ref.q_theory_projected3D, T_ref.q_theory_mean_dir2D_projected3D, ...
    T_ref.q_theory_mean_all];
[~, idx] = min(abs(cands - T_ref.q_theory), [], 2, 'omitnan');
q_pred = nan(height(T_ref), 1);
for i = 1:height(T_ref)
    q_pred(i) = cands(i, idx(i));
end
T = base_prediction_table(T_ref, q_pred, q_pred, ...
    'TheoryBestCandidateOracle', 'TheoryBestCandidateOracle', ...
    'oracle_no_ml', 'diagnostic_only');
end

function T = base_prediction_table(T_ref, q_raw, q_pred, model_name, feature_set, model_type, model_role)
keep = ["condition_id", "step_idx", "realization_idx", "patch_idx", ...
    "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "REQ_M", "REQ_cs_guess", ...
    "M_eff_guess", "M_eff_true_diag", "row_key", "user_field_guess", ...
    "q_user_guess_prior"];
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
T.residual = T.q_pred - T.q_true;
T.abs_error = abs(T.residual);
T = normalize_output_types(T);
end

function T = add_test15_metadata(T, T_base, T_user)
if ~ismember('user_field_guess', string(T.Properties.VariableNames))
    T.user_field_guess = categorical(repmat("not_used", height(T), 1));
end
if ~ismember('q_user_guess_prior', string(T.Properties.VariableNames))
    T.q_user_guess_prior = nan(height(T), 1);
end
if ~ismember('q_theory_dir2D', string(T.Properties.VariableNames))
    [~, loc] = ismember(string(T.row_key), string(T_base.row_key));
    T.q_theory_dir2D = T_base.q_theory_dir2D(loc);
    T.q_theory_diffuse2D = T_base.q_theory_diffuse2D(loc);
    T.q_theory_projected3D = T_base.q_theory_projected3D(loc);
    T.q_theory_mean_dir2D_projected3D = T_base.q_theory_mean_dir2D_projected3D(loc);
    T.q_theory_mean_all = T_base.q_theory_mean_all(loc);
end
if any(ismissing(T.user_field_guess)) && ismember('user_field_guess', string(T_user.Properties.VariableNames))
    [~, loc] = ismember(string(T.row_key), string(T_user.row_key));
    ok = loc > 0;
    T.user_field_guess(ok) = T_user.user_field_guess(loc(ok));
end
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

function plot_metric_bar(T, metric_var, file_path, title_text)
T = T(T.model_type == "bagged_trees" | T.model_type == "baseline_no_ml" | ...
    T.model_type == "oracle_no_ml", :);
T = sortrows(T, metric_var, 'ascend');
labels = short_model_labels(T);
figure('Color', 'w', 'Position', [100 100 980 560]);
barh(categorical(labels), T.(metric_var), 0.68);
xlabel(metric_label(metric_var));
title(title_text, 'Interpreter', 'none', 'FontWeight', 'normal');
grid on;
set(gca, 'YDir', 'reverse', 'FontSize', 11);
export_clean_figure(gcf, file_path);
close(gcf);
end

function plot_q_scatter(T, file_path)
T = T(T.model_type == "bagged_trees" | T.model_type == "baseline_no_ml", :);
models = unique(T.model_name, 'stable');
figure('Color', 'w', 'Position', [100 100 1250 760]);
tl = tiledlayout(2, ceil(numel(models)/2), 'TileSpacing', 'compact', 'Padding', 'compact');
rng(1515);
for i = 1:numel(models)
    ax = nexttile(tl);
    idx = T.model_name == models(i);
    ii = find(idx);
    if numel(ii) > 25000
        ii = ii(randperm(numel(ii), 25000));
    end
    scatter(ax, T.q_true(ii), T.q_pred(ii), 5, 'filled', 'MarkerFaceAlpha', 0.12);
    hold(ax, 'on'); plot(ax, [0 1], [0 1], 'k--');
    axis(ax, 'equal'); xlim(ax, [0 1]); ylim(ax, [0 1]); grid(ax, 'on');
    title(ax, short_model_name(models(i)), 'Interpreter', 'none', 'FontWeight', 'normal');
    xlabel(ax, 'q true'); ylabel(ax, 'q predicted');
end
title(tl, 'Test 15 q true vs predicted', 'Interpreter', 'none', 'FontWeight', 'normal');
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
barh(categorical(short_model_labels(Tc)), ...
    Tc.Delta_MAPE_vs_CurrentHybridBaseline);
xline(0, 'k--');
xlabel('\Delta MAPE vs current hybrid baseline (%)');
title('Negative delta means improvement over current hybrid baseline', ...
    'Interpreter', 'none', 'FontWeight', 'normal');
grid on;
set(gca, 'YDir', 'reverse', 'FontSize', 11);
export_clean_figure(gcf, file_path);
close(gcf);
end

function plot_user_prior_vs_ml(T, file_path)
keep = ismember(T.model_name, ["TheoryCandidatesPlusUserGuessDirect", "UserGuessPriorOnly"]);
plot_metric_bar(T(keep, :), 'MAPE_pct', file_path, 'User prior vs ML direct model');
end

function labels = short_model_labels(T)
labels = strings(height(T), 1);
for i = 1:height(T)
    labels(i) = short_model_name(T.model_name(i));
    if ismember('feature_set', string(T.Properties.VariableNames)) && ...
            strlength(string(T.feature_set(i))) > 0 && ...
            string(T.feature_set(i)) ~= string(T.model_name(i))
        labels(i) = labels(i) + " | " + short_feature_name(T.feature_set(i));
    end
end
end

function out = short_model_name(name)
name = string(name);
switch name
    case "CurrentHybridBaseline"
        out = "Current hybrid";
    case "TheoryCandidatesDirect"
        out = "Theory direct";
    case "TheoryCandidatesPlusUserGuessDirect"
        out = "Theory + user guess";
    case "UserGuessPriorOnly"
        out = "User prior only";
    case "TheoryBestCandidateOracle"
        out = "Best theory oracle";
    otherwise
        out = name;
end
end

function out = short_feature_name(name)
name = string(name);
switch name
    case "CurrentHybridBaseline"
        out = "baseline";
    case "TheoryCandidatesDirect"
        out = "features + theory";
    case "TheoryCandidatesPlusUserGuessDirect"
        out = "features + theory + guess";
    case "UserGuessPriorOnly"
        out = "prior";
    case "TheoryBestCandidateOracle"
        out = "oracle";
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

function register_existing_models(root_dir, PATHS)
deploy_dir = fullfile(PATHS.analysis_dir, 'level_01_model_comparison', 'models', 'deployment');
registry = [
    "test12_hybrid_baseline", "HybridLocalGlobal", "NoCsGuess", "bagged_trees"
    "test12_hybrid_baseline", "HybridLocalGlobal", "WithCsGuess", "bagged_trees"
    "test12_hybrid_baseline", "LocalOnly", "NoCsGuess", "bagged_trees"
    "test12_hybrid_baseline", "GlobalOnly", "NoCsGuess", "bagged_trees"];
for i = 1:size(registry, 1)
    file_i = fullfile(deploy_dir, sprintf('qmodel_%s__%s__%s.mat', ...
        registry(i,2), registry(i,3), registry(i,4)));
    adaptive_req.analysis.register_trained_model( ...
        'RootDir', root_dir, ...
        'SourceModelFile', file_i, ...
        'ModelId', "existing__" + registry(i,2) + "__" + registry(i,3) + "__" + registry(i,4), ...
        'RegistrySubdir', registry(i,1), ...
        'TestName', 'test12_cs_guess_window_sweep', ...
        'AnalysisLevel', 'level_01_model_comparison', ...
        'ModelName', registry(i,2), ...
        'FeatureSet', registry(i,3), ...
        'ModelType', registry(i,4), ...
        'ModelRole', 'operational', ...
        'TrainingDataset', 'test_12_cs_guess_window_sweep', ...
        'Target', 'q_theory', ...
        'PredictorSummary', 'See deployment model bundle.', ...
        'MetricsFile', fullfile(PATHS.analysis_dir, 'level_01_model_comparison', 'tables', 'level12_level01_sws_metrics.csv'), ...
        'Notes', 'Existing model registered by Test 15.');
end
end

function register_test15_models(root_dir, MODELS, OUT, level_name, test_name, model_types)
for i = 1:numel(MODELS)
    model_i = MODELS(i).MODEL;
    adaptive_req.analysis.register_trained_model( ...
        'RootDir', root_dir, ...
        'ModelObject', model_i, ...
        'ModelId', "test15__" + string(model_i.model_name) + "__" + string(model_types), ...
        'RegistrySubdir', 'test15_theory_informed_direct_q', ...
        'TestName', 'test15_theory_informed_direct_q', ...
        'AnalysisLevel', level_name, ...
        'ModelName', model_i.model_name, ...
        'FeatureSet', MODELS(i).feature_set, ...
        'ModelType', model_types, ...
        'ModelRole', model_i.model_role, ...
        'TrainingDataset', test_name, ...
        'Target', 'q_theory', ...
        'PredictorSummary', strjoin(string(MODELS(i).predictors), ', '), ...
        'MetricsFile', fullfile(OUT.table_dir, 'level15_sws_metrics.csv'), ...
        'Notes', 'Theory-informed direct q model; not residual.');
end
end

function register_test15_baselines(root_dir, OUT, level_name, test_name)
adaptive_req.analysis.register_trained_model( ...
    'RootDir', root_dir, ...
    'ModelId', 'test15__UserGuessPriorOnly__baseline_no_ml', ...
    'RegistrySubdir', 'test15_theory_informed_direct_q', ...
    'TestName', 'test15_theory_informed_direct_q', ...
    'AnalysisLevel', level_name, ...
    'ModelName', 'UserGuessPriorOnly', ...
    'FeatureSet', 'UserGuessPriorOnly', ...
    'ModelType', 'baseline_no_ml', ...
    'ModelRole', 'baseline_no_ml', ...
    'TrainingDataset', test_name, ...
    'Target', 'q_theory', ...
    'PredictorSummary', 'q_user_guess_prior', ...
    'MetricsFile', fullfile(OUT.table_dir, 'level15_sws_metrics.csv'), ...
    'Notes', 'No-ML prior baseline.');
end

function make_level15_figures(T_sws, T_sws_metrics, T_sws_by_user, T_sws_by_M_eff, OUT)
plot_metric_bar(T_sws_metrics, 'MAPE_pct', ...
    fullfile(OUT.fig_dir, 'level15_mape_by_model.png'), ...
    'Test 15 MAPE by model');
plot_metric_bar(T_sws_metrics, 'HighError_gt20_pct', ...
    fullfile(OUT.fig_dir, 'level15_high_error_by_model.png'), ...
    'Test 15 high-error >20% by model');
plot_q_scatter(T_sws, fullfile(OUT.fig_dir, 'level15_q_true_vs_pred.png'));
plot_delta_vs_baseline(T_sws_metrics, ...
    fullfile(OUT.fig_dir, 'level15_delta_mape_vs_CurrentHybridBaseline.png'));
plot_group_metric(T_sws_by_user, "user_field_guess", "model_name", "MAPE_pct", ...
    'Test 15 MAPE by user field guess', ...
    fullfile(OUT.fig_dir, 'level15_mape_by_user_field_guess.png'));
plot_group_metric(T_sws_by_M_eff, "M_eff_true_diag_bin", "model_name", "MAPE_pct", ...
    'Test 15 MAPE by M_eff_true_diag diagnostic bin', ...
    fullfile(OUT.fig_dir, 'level15_mape_by_M_eff_true_diag.png'));
plot_user_prior_vs_ml(T_sws_metrics, ...
    fullfile(OUT.fig_dir, 'level15_user_guess_prior_vs_ml_direct.png'));
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

function print_summary(T_sws, T_user, T_meff)
Top = T_sws(T_sws.model_role ~= "diagnostic_only", :);
[~, ib] = min(Top.MAPE_pct);
[~, ih] = min(Top.HighError_gt20_pct);
fprintf('\nBest operational model by MAPE:\n');
disp(Top(ib, {'model_name','feature_set','model_type','MAPE_pct','HighError_gt20_pct'}));
fprintf('\nBest operational model by HighError_gt20:\n');
disp(Top(ih, {'model_name','feature_set','model_type','MAPE_pct','HighError_gt20_pct'}));

show_delta(T_sws, "TheoryCandidatesDirect", "CurrentHybridBaseline", ...
    'Theory-q candidates vs current hybrid');
show_delta(T_sws, "TheoryCandidatesPlusUserGuessDirect", "TheoryCandidatesDirect", ...
    'Adding user_field_guess vs theory candidates only');

prior = T_sws(T_sws.model_name == "UserGuessPriorOnly", :);
if ~isempty(prior)
    fprintf('\nUserGuessPriorOnly:\n');
    disp(prior(:, {'MAPE_pct','HighError_gt20_pct','P95_abs_error_pct'}));
end

fprintf('\nBest regimes by user_field_guess:\n');
disp(sortrows(T_user(T_user.model_role ~= "diagnostic_only", ...
    {'model_name','user_field_guess','MAPE_pct','HighError_gt20_pct'}), 'MAPE_pct'));
fprintf('\nWorst regimes by M_eff_true_diag_bin:\n');
disp(sortrows(T_meff(T_meff.model_role ~= "diagnostic_only", ...
    {'model_name','M_eff_true_diag_bin','MAPE_pct','HighError_gt20_pct'}), ...
    'MAPE_pct', 'descend'));
end

function show_delta(T, model_a, model_b, label)
Ta = T(T.model_name == model_a, :);
Tb = T(T.model_name == model_b, :);
if isempty(Ta) || isempty(Tb), return; end
fprintf('\n%s:\n', label);
fprintf('Delta MAPE = %.3f%%; Delta HighError>20 = %.3f%%\n', ...
    Ta.MAPE_pct(1) - Tb.MAPE_pct(1), ...
    Ta.HighError_gt20_pct(1) - Tb.HighError_gt20_pct(1));
end
