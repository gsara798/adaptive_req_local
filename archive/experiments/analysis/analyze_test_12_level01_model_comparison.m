%% analyze_test_12_level01_model_comparison.m
% Test 12 Level 01: model and feature-set comparison.
%
% This analysis does not regenerate Test 12. It compares LocalOnly,
% GlobalOnly, and HybridLocalGlobal models while ablating REQ_cs_guess and
% M_eff_guess. M_eff_true_diag is diagnostic-only.

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

%% Load Test 12

[T_feat, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    'test_12_cs_guess_window_sweep', ...
    'RootDir', root_dir, ...
    'Verbose', true);

OUT = adaptive_req.analysis.Test12Analysis.makeOutputDirs( ...
    PATHS, 'level_01_model_comparison');
deploy_model_dir = fullfile(OUT.model_dir, 'deployment');
if ~exist(deploy_model_dir, 'dir')
    mkdir(deploy_model_dir);
end

required_vars = [
    "q_theory"
    "req_mapping"
    "global_req_mapping"
    "SIM_f0"
    "SIM_cs_bg"
    "REQ_M"
    "REQ_cs_guess"
    "condition_id"
    "step_idx"
    "realization_idx"
    "patch_idx"];

adaptive_req.analysis.Test12Analysis.requireVars(T_feat, required_vars, ...
    'Test 12 Level 01');

if ~ismember('M_eff_guess', string(T_feat.Properties.VariableNames))
    warning('M_eff_guess is missing. WithMeffGuess variants will omit it.');
end
if ~ismember('M_eff_true_diag', string(T_feat.Properties.VariableNames))
    warning('M_eff_true_diag is missing. DiagnosticWithMeffTrue will omit it.');
end

T_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_feat);
T_feat = adaptive_req.analysis.Test12Analysis.addBins(T_feat);

%% Model specs and clean condition split

model_specs = adaptive_req.analysis.Test12Analysis.buildModelSpecs(T_feat);
model_types = ["linear", "boosted_trees", "bagged_trees"];

[train_mask, test_mask] = adaptive_req.analysis.Test12Analysis.conditionSplit( ...
    T_feat, 0.70, 12001);

assert(~any(train_mask & test_mask), 'Train/test masks overlap.');
assert(any(train_mask) && any(test_mask), 'Train/test split is empty.');

%% Train models

T_all_pred = table();
T_q_metrics = table();
MODELS = struct([]);

for i = 1:numel(model_specs)
    spec = model_specs(i);
    fprintf('\n=== Level 01 | %s | %s | %s ===\n', ...
        spec.model_name, spec.feature_set, spec.model_role);

    [MODEL_i, T_pred_i, T_metrics_i] = ...
        adaptive_req.analysis.train_q_model_fixed_split( ...
            T_feat, spec.predictors, train_mask, test_mask, ...
            'ModelName', spec.model_name, ...
            'ModelRole', spec.model_role, ...
            'ModelTypes', model_types, ...
            'NumLearningCycles', 200, ...
            'MinLeafSize', 8, ...
            'Verbose', true);

    T_pred_i = adaptive_req.analysis.Test12Analysis.addModelMetadata( ...
        T_pred_i, spec);
    T_metrics_i = adaptive_req.analysis.Test12Analysis.addModelMetadata( ...
        T_metrics_i, spec);

    MODELS(i).model_name = spec.model_name;
    MODELS(i).feature_set = spec.feature_set;
    MODELS(i).model_role = spec.model_role;
    MODELS(i).predictors = spec.predictors;
    MODELS(i).model = MODEL_i;

    deploy_files_i = adaptive_req.analysis.save_q_model_deployment( ...
        MODEL_i, deploy_model_dir, ...
        'ModelName', spec.model_name, ...
        'FeatureSet', spec.feature_set, ...
        'ModelRole', spec.model_role, ...
        'ModelTypes', model_types, ...
        'Overwrite', true);
    MODELS(i).deployment_files = deploy_files_i;

    T_all_pred = adaptive_req.analysis.Test12Analysis.concatTables( ...
        T_all_pred, T_pred_i);
    T_q_metrics = adaptive_req.analysis.Test12Analysis.concatTables( ...
        T_q_metrics, T_metrics_i);
end

assert(~isempty(T_all_pred), 'No predictions were generated.');
assert(~isempty(T_q_metrics), 'No q metrics were generated.');

%% Convert q to local SWS and summarize

T_sws = adaptive_req.analysis.Test12Analysis.addSwsMetrics(T_all_pred, T_feat);
T_sws = adaptive_req.analysis.Test12Analysis.addBins(T_sws);
T_sws_test = T_sws(T_sws.split == "test", :);

assert(~isempty(T_sws_test), 'No test predictions were generated.');
assert(mean(isfinite(T_sws_test.q_pred)) > 0.95, ...
    'Too many non-finite q_pred values.');
assert(mean(isfinite(T_sws_test.cs_pred)) > 0.95, ...
    'Too many non-finite cs_pred values.');

base_groups = [
    "model_name"
    "feature_set"
    "model_type"
    "model_role"];

T_sws_metrics = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, base_groups);
T_by_feature_set = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, base_groups);
T_by_cs_guess = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, [base_groups; "REQ_cs_guess"]);
T_by_M = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, [base_groups; "REQ_M"]);
T_by_frequency = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, [base_groups; "SIM_f0"]);

T_by_Meff_guess = table();
if ismember('M_eff_guess_bin', string(T_sws_test.Properties.VariableNames))
    T_by_Meff_guess = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
        T_sws_test, [base_groups; "M_eff_guess_bin"]);
end

T_by_Meff_true = table();
if ismember('M_eff_true_diag_bin', string(T_sws_test.Properties.VariableNames))
    T_by_Meff_true = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
        T_sws_test, [base_groups; "M_eff_true_diag_bin"]);
end

T_delta = adaptive_req.analysis.Test12Analysis.deltaVsNoCsGuess( ...
    T_sws_metrics, ["model_name"; "model_type"; "model_role"]);

%% Save tables and models

writetable(T_q_metrics, fullfile(OUT.table_dir, ...
    'level12_level01_q_metrics.csv'));
writetable(T_sws_metrics, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics.csv'));
writetable(adaptive_req.analysis.Test12Analysis.removeCellColumns(T_sws), ...
    fullfile(OUT.table_dir, 'level12_level01_predictions.csv'));
writetable(T_by_feature_set, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_feature_set.csv'));
writetable(T_by_Meff_guess, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_M_eff_guess.csv'));
writetable(T_by_Meff_true, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_M_eff_true_diag.csv'));
writetable(T_by_cs_guess, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_cs_guess.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_M.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, ...
    'level12_level01_sws_metrics_by_frequency.csv'));
writetable(T_delta, fullfile(OUT.table_dir, ...
    'level12_level01_delta_mape_vs_NoCsGuess.csv'));

save(fullfile(OUT.model_dir, ...
    'level12_level01_model_comparison_models.mat'), ...
    'MODELS', 'T_q_metrics', 'T_sws_metrics', 'T_delta', ...
    'T_by_feature_set', 'T_by_Meff_guess', 'T_by_Meff_true', ...
    'T_by_cs_guess', 'T_by_M', 'T_by_frequency', 'MC', 'PATHS', '-v7.3');

%% Figures

T_fig = T_sws_metrics(T_sws_metrics.model_type == "bagged_trees" & ...
    T_sws_metrics.model_role == "operational", :);
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T_fig, "feature_set", "model_name", "MAPE_pct", ...
    'Test 12 Level 01 MAPE by model and feature set', ...
    fullfile(OUT.fig_dir, 'level12_level01_mape_by_model_and_feature_set.png'));

adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T_fig, "feature_set", "model_name", "HighError_gt20_pct", ...
    'Test 12 Level 01 high-error rate by model and feature set', ...
    fullfile(OUT.fig_dir, 'level12_level01_high_error_by_model_and_feature_set.png'));

adaptive_req.analysis.Test12Analysis.plotQScatter( ...
    T_sws_test(T_sws_test.model_role == "operational", :), ...
    'Test 12 Level 01 q true vs predicted', ...
    fullfile(OUT.fig_dir, 'level12_level01_q_true_vs_pred.png'));

adaptive_req.analysis.Test12Analysis.plotErrorBox( ...
    T_sws_test(T_sws_test.model_role == "operational", :), ...
    'Test 12 Level 01 SWS error distribution', ...
    fullfile(OUT.fig_dir, 'level12_level01_sws_error_boxplot.png'));

adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T_by_cs_guess(T_by_cs_guess.model_type == "bagged_trees" & ...
    T_by_cs_guess.model_role == "operational", :), ...
    "REQ_cs_guess", "model_name", "MAPE_pct", ...
    'Test 12 Level 01 MAPE by REQ cs guess', ...
    fullfile(OUT.fig_dir, 'level12_level01_mape_by_cs_guess.png'));

adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T_by_Meff_guess(T_by_Meff_guess.model_type == "bagged_trees" & ...
    T_by_Meff_guess.model_role == "operational", :), ...
    "M_eff_guess_bin", "model_name", "MAPE_pct", ...
    'Test 12 Level 01 MAPE by M eff guess', ...
    fullfile(OUT.fig_dir, 'level12_level01_mape_by_M_eff_guess.png'));

adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T_by_Meff_true(T_by_Meff_true.model_type == "bagged_trees" & ...
    T_by_Meff_true.model_role == "operational", :), ...
    "M_eff_true_diag_bin", "model_name", "MAPE_pct", ...
    'Test 12 Level 01 MAPE by M eff true diagnostic', ...
    fullfile(OUT.fig_dir, 'level12_level01_mape_by_M_eff_true_diag.png'));

plot_delta_mape(T_delta, OUT.fig_dir);

%% Console report

fprintf('\nTest 12 Level 01 complete.\n');
fprintf('Analysis folder:\n%s\n', OUT.analysis_dir);

T_console = sortrows(T_fig(:, {'model_name', 'feature_set', ...
    'MAPE_pct', 'HighError_gt20_pct'}), 'MAPE_pct', 'ascend');
fprintf('\nBagged-trees operational summary.\n');
disp(T_console);

fprintf('\nDelta vs NoCsGuess. Negative deltas mean improvement.\n');
disp(T_delta(T_delta.model_type == "bagged_trees" & ...
    T_delta.model_role == "operational", ...
    {'model_name', 'feature_set', 'Delta_MAPE_vs_NoCsGuess', ...
    'Delta_HighError_gt20_vs_NoCsGuess'}));

%% Local functions

function plot_delta_mape(T_delta, fig_dir)

T = T_delta(T_delta.model_type == "bagged_trees" & ...
    T_delta.model_role == "operational", :);
if isempty(T)
    return;
end

models = unique(T.model_name, 'stable');
features = unique(T.feature_set, 'stable');
Y_mape = nan(numel(features), numel(models));
Y_high = nan(numel(features), numel(models));

for i = 1:numel(features)
    for j = 1:numel(models)
        idx = T.feature_set == features(i) & T.model_name == models(j);
        if any(idx)
            Y_mape(i, j) = mean(T.Delta_MAPE_vs_NoCsGuess(idx), 'omitnan');
            Y_high(i, j) = mean(T.Delta_HighError_gt20_vs_NoCsGuess(idx), 'omitnan');
        end
    end
end

figure('Color', 'w', 'Position', [100 100 1150 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(categorical(features), Y_mape);
yline(0, 'k-');
ylabel('\Delta MAPE vs NoCsGuess (%)');
title('MAPE delta. Negative is better.');
legend(models, 'Location', 'best', 'Interpreter', 'none');
grid on;

nexttile;
bar(categorical(features), Y_high);
yline(0, 'k-');
ylabel('\Delta high-error >20% vs NoCsGuess');
title('High-error delta. Negative is better.');
legend(models, 'Location', 'best', 'Interpreter', 'none');
grid on;

exportgraphics(gcf, fullfile(fig_dir, ...
    'level12_level01_delta_mape_vs_NoCsGuess.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end
