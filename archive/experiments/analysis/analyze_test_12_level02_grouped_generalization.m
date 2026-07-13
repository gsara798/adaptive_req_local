%% analyze_test_12_level02_grouped_generalization.m
% Test 12 Level 02: grouped generalization for cs-guess/window sweep.
%
% This analysis does not regenerate Test 12. It evaluates leave-one-group-out
% splits and compares whether REQ_cs_guess and M_eff_guess help.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLegendFontSize', 11);

%% Runtime options

% Profiles:
%   "fast" : operational feature sets, bagged trees only, fewer trees.
%   "full" : operational + diagnostic feature sets, all requested model types.
%
% Start with "fast". It answers the main Test 12 question much faster:
% do REQ_cs_guess and/or M_eff_guess improve grouped generalization?
ANALYSIS.profile = "fast";
ANALYSIS.resume_from_checkpoint = true;
ANALYSIS.save_checkpoint_each_job = true;
ANALYSIS.use_parallel = true;
ANALYSIS.verbose_training = false;

switch ANALYSIS.profile
    case "fast"
        ANALYSIS.model_types = "bagged_trees";
        ANALYSIS.num_learning_cycles = 80;
        ANALYSIS.min_leaf_size = 12;
        ANALYSIS.include_diagnostic = false;
        ANALYSIS.keep_feature_sets = [
            "NoCsGuess"
            "WithCsGuess"
            "WithMeffGuess"
            "WithCsGuessAndMeffGuess"];
    case "full"
        ANALYSIS.model_types = ["linear", "boosted_trees", "bagged_trees"];
        ANALYSIS.num_learning_cycles = 200;
        ANALYSIS.min_leaf_size = 8;
        ANALYSIS.include_diagnostic = true;
        ANALYSIS.keep_feature_sets = [
            "NoCsGuess"
            "WithCsGuess"
            "WithMeffGuess"
            "WithCsGuessAndMeffGuess"
            "DiagnosticWithMeffTrue"];
    otherwise
        error('Unknown ANALYSIS.profile: %s', ANALYSIS.profile);
end

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
    PATHS, 'level_02_grouped_generalization');

required_vars = [
    "q_theory"
    "req_mapping"
    "SIM_f0"
    "SIM_cs_bg"
    "REQ_M"
    "REQ_cs_guess"
    "condition_id"
    "step_idx"
    "realization_idx"
    "patch_idx"];

adaptive_req.analysis.Test12Analysis.requireVars(T_feat, required_vars, ...
    'Test 12 Level 02');

T_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T_feat);
T_feat = adaptive_req.analysis.Test12Analysis.addBins(T_feat);

%% Model and split specs

model_specs = adaptive_req.analysis.Test12Analysis.buildModelSpecs(T_feat);
model_specs = filter_model_specs(model_specs, ANALYSIS);
model_types = ANALYSIS.model_types;

if ANALYSIS.use_parallel
    start_parallel_pool_if_available();
end

split_specs = struct([]);
split_specs(1).name = "leave_one_frequency";
split_specs(1).var = "SIM_f0";
split_specs(2).name = "leave_one_M";
split_specs(2).var = "REQ_M";
split_specs(3).name = "leave_one_cs_guess";
split_specs(3).var = "REQ_cs_guess";
split_specs(4).name = "leave_one_cs_true";
split_specs(4).var = "SIM_cs_bg";
split_specs(5).name = "leave_one_aperture";
split_specs(5).var = "step_idx";
split_specs(6).name = "leave_one_wave_model";
split_specs(6).var = "SIM_WaveModel";

%% Train grouped models

T_all_pred = table();
T_q_metrics = table();
MODELS = struct([]);
model_counter = 0;
completed_jobs = strings(0, 1);

checkpoint_file = fullfile(OUT.model_dir, ...
    sprintf('level12_level02_%s_checkpoint.mat', ANALYSIS.profile));

if ANALYSIS.resume_from_checkpoint && exist(checkpoint_file, 'file')
    fprintf('\nLoading Level 02 checkpoint:\n%s\n', checkpoint_file);
    S = load(checkpoint_file, 'T_all_pred', 'T_q_metrics', 'MODELS', ...
        'model_counter', 'completed_jobs');
    T_all_pred = S.T_all_pred;
    T_q_metrics = S.T_q_metrics;
    MODELS = S.MODELS;
    model_counter = S.model_counter;
    completed_jobs = string(S.completed_jobs);
    fprintf('Completed jobs in checkpoint: %d\n', numel(completed_jobs));
end

for sidx = 1:numel(split_specs)
    split_name = string(split_specs(sidx).name);
    heldout_var = string(split_specs(sidx).var);

    if ~ismember(heldout_var, string(T_feat.Properties.VariableNames))
        fprintf('Skipping %s because %s is missing.\n', split_name, heldout_var);
        continue;
    end

    heldout_values = unique(T_feat.(char(heldout_var)), 'stable');
    if numel(heldout_values) < 2
        warning('Skipping %s because %s has only one group.', ...
            split_name, heldout_var);
        continue;
    end

    for hidx = 1:numel(heldout_values)
        heldout_value = heldout_values(hidx);
        test_mask = adaptive_req.analysis.Test12Analysis.isGroupValue( ...
            T_feat.(char(heldout_var)), heldout_value);
        train_mask = ~test_mask;

        adaptive_req.analysis.Test12Analysis.assertGroupedSplit( ...
            T_feat, train_mask, test_mask, heldout_var, heldout_value);

        for midx = 1:numel(model_specs)
            spec = model_specs(midx);
            job_key = make_job_key(split_name, heldout_var, ...
                adaptive_req.analysis.Test12Analysis.valueToString(heldout_value), ...
                spec);

            if any(completed_jobs == job_key)
                fprintf('Skipping completed job: %s\n', job_key);
                continue;
            end

            fprintf('\n=== Level 02 | %s | %s = %s | %s | %s ===\n', ...
                split_name, heldout_var, ...
                adaptive_req.analysis.Test12Analysis.valueToString(heldout_value), ...
                spec.model_name, spec.feature_set);

            [MODEL_i, T_pred_i, T_metrics_i] = ...
                adaptive_req.analysis.train_q_model_fixed_split( ...
                    T_feat, spec.predictors, train_mask, test_mask, ...
                    'ModelName', spec.model_name, ...
                    'ModelRole', spec.model_role, ...
                    'ModelTypes', model_types, ...
                    'NumLearningCycles', ANALYSIS.num_learning_cycles, ...
                    'MinLeafSize', ANALYSIS.min_leaf_size, ...
                    'UseParallel', ANALYSIS.use_parallel, ...
                    'Verbose', ANALYSIS.verbose_training);

            T_pred_i = adaptive_req.analysis.Test12Analysis.addModelMetadata( ...
                T_pred_i, spec, ...
                'GeneralizationTest', split_name, ...
                'HeldoutVar', heldout_var, ...
                'HeldoutValue', adaptive_req.analysis.Test12Analysis.valueToString(heldout_value), ...
                'NTrain', sum(train_mask), ...
                'NTest', sum(test_mask));

            T_metrics_i = adaptive_req.analysis.Test12Analysis.addModelMetadata( ...
                T_metrics_i, spec, ...
                'GeneralizationTest', split_name, ...
                'HeldoutVar', heldout_var, ...
                'HeldoutValue', adaptive_req.analysis.Test12Analysis.valueToString(heldout_value), ...
                'NTrain', sum(train_mask), ...
                'NTest', sum(test_mask));

            model_counter = model_counter + 1;
            MODELS(model_counter).generalization_test = split_name;
            MODELS(model_counter).heldout_var = heldout_var;
            MODELS(model_counter).heldout_value = adaptive_req.analysis.Test12Analysis.valueToString(heldout_value);
            MODELS(model_counter).model_name = spec.model_name;
            MODELS(model_counter).feature_set = spec.feature_set;
            MODELS(model_counter).model_role = spec.model_role;
            MODELS(model_counter).predictors = spec.predictors;
            MODELS(model_counter).model = MODEL_i;

            T_pred_test_i = T_pred_i(T_pred_i.split == "test", :);
            T_all_pred = adaptive_req.analysis.Test12Analysis.concatTables( ...
                T_all_pred, T_pred_test_i);
            T_q_metrics = adaptive_req.analysis.Test12Analysis.concatTables( ...
                T_q_metrics, T_metrics_i);

            completed_jobs(end + 1, 1) = job_key; %#ok<SAGROW>

            if ANALYSIS.save_checkpoint_each_job
                save(checkpoint_file, 'T_all_pred', 'T_q_metrics', ...
                    'MODELS', 'model_counter', 'completed_jobs', ...
                    'ANALYSIS', '-v7.3');
            end
        end
    end
end

assert(~isempty(T_all_pred), 'No grouped predictions were generated.');
assert(~isempty(T_q_metrics), 'No grouped q metrics were generated.');

%% Convert q to SWS and summarize

T_sws = adaptive_req.analysis.Test12Analysis.addSwsMetrics(T_all_pred, T_feat);
T_sws = adaptive_req.analysis.Test12Analysis.addBins(T_sws);
T_sws_test = T_sws(T_sws.split == "test", :);

assert(~isempty(T_sws_test), 'No grouped test predictions were generated.');
assert(mean(isfinite(T_sws_test.q_pred)) > 0.95, ...
    'Too many non-finite q_pred values.');
assert(mean(isfinite(T_sws_test.cs_pred)) > 0.95, ...
    'Too many non-finite cs_pred values.');

main_groups = [
    "generalization_test"
    "heldout_var"
    "heldout_value"
    "model_name"
    "feature_set"
    "model_type"
    "model_role"];

T_sws_metrics = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, main_groups);
T_by_heldout = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, ["heldout_var"; "heldout_value"; "model_name"; ...
    "feature_set"; "model_type"; "model_role"]);
T_by_feature_set = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, ["generalization_test"; "model_name"; "feature_set"; ...
    "model_type"; "model_role"]);
T_by_model = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
    T_sws_test, ["generalization_test"; "model_name"; "feature_set"; ...
    "model_type"; "model_role"]);

T_high_error = T_sws_metrics(:, {'generalization_test', 'heldout_var', ...
    'heldout_value', 'model_name', 'feature_set', 'model_type', ...
    'model_role', 'N', 'HighError_gt20_pct', 'MAPE_pct'});
T_high_error = sortrows(T_high_error, 'HighError_gt20_pct', 'descend');

T_delta = adaptive_req.analysis.Test12Analysis.deltaVsNoCsGuess( ...
    T_sws_metrics, ["generalization_test"; "heldout_var"; "heldout_value"; ...
    "model_name"; "model_type"; "model_role"]);

%% Save tables and models

writetable(adaptive_req.analysis.Test12Analysis.removeCellColumns(T_sws), ...
    fullfile(OUT.table_dir, 'level12_level02_grouped_predictions.csv'));
writetable(T_q_metrics, fullfile(OUT.table_dir, ...
    'level12_level02_q_metrics.csv'));
writetable(T_sws_metrics, fullfile(OUT.table_dir, ...
    'level12_level02_sws_metrics.csv'));
writetable(T_by_heldout, fullfile(OUT.table_dir, ...
    'level12_level02_sws_metrics_by_heldout.csv'));
writetable(T_by_feature_set, fullfile(OUT.table_dir, ...
    'level12_level02_sws_metrics_by_feature_set.csv'));
writetable(T_by_model, fullfile(OUT.table_dir, ...
    'level12_level02_sws_metrics_by_model.csv'));
writetable(T_high_error, fullfile(OUT.table_dir, ...
    'level12_level02_high_error_gt20.csv'));
writetable(T_delta, fullfile(OUT.table_dir, ...
    'level12_level02_delta_mape_vs_NoCsGuess.csv'));

save(fullfile(OUT.model_dir, ...
    'level12_level02_grouped_generalization_models.mat'), ...
    'MODELS', 'T_q_metrics', 'T_sws_metrics', 'T_by_heldout', ...
    'T_by_feature_set', 'T_by_model', 'T_high_error', 'T_delta', ...
    'MC', 'PATHS', 'ANALYSIS', '-v7.3');

%% Figures

plot_mape_by_heldout(T_sws_metrics, "leave_one_frequency", ...
    OUT.fig_dir, 'level12_level02_mape_leave_one_frequency.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_M", ...
    OUT.fig_dir, 'level12_level02_mape_leave_one_M.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_cs_guess", ...
    OUT.fig_dir, 'level12_level02_mape_leave_one_cs_guess.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_cs_true", ...
    OUT.fig_dir, 'level12_level02_mape_leave_one_cs_true.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_aperture", ...
    OUT.fig_dir, 'level12_level02_mape_leave_one_aperture.png');
plot_summary_by_test(T_by_feature_set, OUT.fig_dir);
plot_delta_by_test(T_delta, "WithCsGuess", OUT.fig_dir, ...
    'level12_level02_delta_mape_with_cs_guess.png');
plot_delta_by_test(T_delta, "WithMeffGuess", OUT.fig_dir, ...
    'level12_level02_delta_mape_with_M_eff_guess.png');
plot_q_by_test(T_sws_test, OUT.fig_dir);
plot_error_box_by_test(T_sws_test, OUT.fig_dir);

%% Console report

fprintf('\nTest 12 Level 02 grouped generalization complete.\n');
fprintf('Analysis folder:\n%s\n', OUT.analysis_dir);

T_op = T_by_feature_set(T_by_feature_set.model_type == "bagged_trees" & ...
    T_by_feature_set.model_role == "operational", :);
T_best = sortrows(T_op, 'MAPE_pct', 'ascend');

fprintf('\nBest operational grouped results, bagged trees.\n');
disp(T_best(1:min(15, height(T_best)), {'generalization_test', ...
    'model_name', 'feature_set', 'MAPE_pct', 'HighError_gt20_pct'}));

fprintf('\nFeature-set deltas vs NoCsGuess. Negative means improvement.\n');
T_delta_console = T_delta(T_delta.model_type == "bagged_trees" & ...
    T_delta.model_role == "operational", :);
disp(sortrows(T_delta_console(:, {'generalization_test', 'heldout_value', ...
    'model_name', 'feature_set', 'Delta_MAPE_vs_NoCsGuess', ...
    'Delta_HighError_gt20_vs_NoCsGuess'}), ...
    'Delta_MAPE_vs_NoCsGuess', 'ascend'));

%% Local functions

function specs = filter_model_specs(specs, ANALYSIS)

keep = true(numel(specs), 1);
for i = 1:numel(specs)
    keep(i) = any(string(specs(i).feature_set) == ANALYSIS.keep_feature_sets);
    if ~ANALYSIS.include_diagnostic && string(specs(i).model_role) ~= "operational"
        keep(i) = false;
    end
end

specs = specs(keep);

assert(~isempty(specs), 'No model specs remain after applying runtime filters.');

fprintf('\nLevel 02 runtime profile: %s\n', ANALYSIS.profile);
fprintf('Model types          : %s\n', strjoin(string(ANALYSIS.model_types), ', '));
fprintf('Learning cycles      : %d\n', ANALYSIS.num_learning_cycles);
fprintf('Min leaf size        : %d\n', ANALYSIS.min_leaf_size);
fprintf('Use parallel         : %d\n', ANALYSIS.use_parallel);
fprintf('Diagnostic included  : %d\n', ANALYSIS.include_diagnostic);
fprintf('Model specs to train : %d\n', numel(specs));

end

function start_parallel_pool_if_available()

if exist('gcp', 'file') ~= 2 || exist('parpool', 'file') ~= 2
    warning('Parallel Computing Toolbox was not found. Continuing serially.');
    return;
end

try
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool('threads');
        catch
            parpool;
        end
    end
catch ME
    warning(ME.identifier, ...
        'Could not start a parallel pool. Continuing serially. Reason: %s', ...
        ME.message);
end

end

function key = make_job_key(split_name, heldout_var, heldout_value, spec)

key = strjoin([
    string(split_name)
    string(heldout_var)
    string(heldout_value)
    string(spec.model_name)
    string(spec.feature_set)
    string(spec.model_role)], "|");

end

function plot_mape_by_heldout(T_metrics, test_name, fig_dir, file_name)

T = T_metrics(T_metrics.generalization_test == test_name & ...
    T_metrics.model_type == "bagged_trees" & ...
    T_metrics.model_role == "operational", :);
if isempty(T)
    return;
end

T.series = T.model_name + " | " + T.feature_set;
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "heldout_value", "series", "MAPE_pct", ...
    sprintf('Test 12 Level 02 MAPE: %s', test_name), ...
    fullfile(fig_dir, file_name));

end

function plot_summary_by_test(T_metrics, fig_dir)

T = T_metrics(T_metrics.model_type == "bagged_trees" & ...
    T_metrics.model_role == "operational", :);
if isempty(T)
    return;
end
T.series = T.model_name + " | " + T.feature_set;
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "generalization_test", "series", "MAPE_pct", ...
    'Test 12 Level 02 MAPE summary by grouped generalization test', ...
    fullfile(fig_dir, 'level12_level02_mape_summary_by_generalization_test.png'));

end

function plot_delta_by_test(T_delta, feature_set, fig_dir, file_name)

T = T_delta(T_delta.feature_set == feature_set & ...
    T_delta.model_type == "bagged_trees" & ...
    T_delta.model_role == "operational", :);
if isempty(T)
    return;
end
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "generalization_test", "model_name", "Delta_MAPE_vs_NoCsGuess", ...
    sprintf('Delta MAPE vs NoCsGuess: %s', feature_set), ...
    fullfile(fig_dir, file_name));

end

function plot_q_by_test(T_sws_test, fig_dir)

T = T_sws_test(T_sws_test.model_type == "bagged_trees" & ...
    T_sws_test.model_role == "operational", :);
if isempty(T)
    return;
end
T.panel = T.generalization_test + " | " + T.model_name;
panels = unique(T.panel, 'stable');

figure('Color', 'w', 'Position', [100 100 1450 900]);
tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(panels)
    ax = nexttile(tl);
    Ti = T(T.panel == panels(i), :);
    scatter(ax, Ti.q_true, Ti.q_pred, 6, 'filled', 'MarkerFaceAlpha', 0.25);
    hold(ax, 'on');
    plot(ax, [0 1], [0 1], 'k--');
    axis(ax, 'equal');
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    grid(ax, 'on');
    title(ax, panels(i), 'Interpreter', 'none');
    xlabel(ax, 'q true');
    ylabel(ax, 'q predicted');
end
title(tl, 'Test 12 Level 02 q true vs predicted');
exportgraphics(gcf, fullfile(fig_dir, ...
    'level12_level02_q_true_vs_pred_by_test.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end

function plot_error_box_by_test(T_sws_test, fig_dir)

T = T_sws_test(T_sws_test.model_type == "bagged_trees" & ...
    T_sws_test.model_role == "operational", :);
if isempty(T)
    return;
end
labels = T.generalization_test + " | " + T.model_name + " | " + T.feature_set;
figure('Color', 'w', 'Position', [100 100 1500 620]);
boxchart(categorical(labels), T.abs_sws_error_pct);
yline(20, 'r--', '20%');
ylabel('|SWS error| (%)');
title('Test 12 Level 02 SWS error by grouped test');
xtickangle(30);
grid on;
exportgraphics(gcf, fullfile(fig_dir, ...
    'level12_level02_sws_error_box_by_test.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end
