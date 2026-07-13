%% analyze_test_08_level14_compare_advanced_features.m
% Level 14: Compare advanced angular feature families for q prediction.
%
% Operational models do not use WaveModel, true SWS, q-derived quantities,
% or aperture. Aperture and true-M_eff variants are kept as diagnostic
% ceilings because they answer different scientific questions.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load MC results

experiment_name = 'test_08_advanced_angular_features';

[T_mc, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

COL = adaptive_req.analysis.detect_mc_columns(T_mc);

%% Add effective-window variables

cs_guess = 3.0;

T_mc_eff = adaptive_req.analysis.add_effective_window_metrics( ...
    T_mc, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', cs_guess);

if ~ismember("M_eff", string(T_mc_eff.Properties.VariableNames))
    T_mc_eff.M_eff = T_mc_eff.M_eff_true_diag;
end

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_14_compare_advanced_features');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Predictor families

entropy_features = [
    "ang_entropy"
    "radial_entropy"
];

known_params = [
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"
];

rich_spectral_features = [
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

angular_shape_features = [
    "ang_moment_2"
    "ang_moment_4"
    "ang_peak_count_rel"
    "ang_top2_window_frac"
    "ang_top3_window_frac"
    "ang_top2_to_top1"
    "ang_peak_separation_deg"
];

model_specs = struct([]);

model_specs(1).name = "ModelC_rich_spectral_known_params";
model_specs(1).role = "operational_baseline";
model_specs(1).predictors = [
    entropy_features
    known_params
    rich_spectral_features
];

model_specs(2).name = "ModelE_advanced_angular_shape";
model_specs(2).role = "operational";
model_specs(2).predictors = [
    entropy_features
    known_params
    rich_spectral_features
    angular_shape_features
];

model_specs(3).name = "ModelF_operational_with_Meff_guess";
model_specs(3).role = "operational";
model_specs(3).predictors = [
    entropy_features
    known_params
    rich_spectral_features
    angular_shape_features
    "M_eff_guess"
];

model_specs(4).name = "ModelG_diagnostic_with_aperture";
model_specs(4).role = "diagnostic_aperture";
model_specs(4).predictors = [
    entropy_features
    known_params
    rich_spectral_features
    angular_shape_features
    "Omega_sr"
];

model_specs(5).name = "ModelD_diagnostic_with_Meff_true";
model_specs(5).role = "diagnostic_true_sws";
model_specs(5).predictors = [
    entropy_features
    known_params
    rich_spectral_features
    angular_shape_features
    "M_eff"
];

%% Model types

model_types = [
    "boosted_trees"
    "bagged_trees"
    "quadratic"
];

%% Train all models

ALL = struct();

T_all_metrics = table();
T_all_predictions = table();

for sidx = 1:numel(model_specs)

    spec = model_specs(sidx);

    fprintf('\n============================================================\n');
    fprintf('Training %s (%s)\n', spec.name, spec.role);
    fprintf('============================================================\n');

    predictors_i = keep_existing_predictors(T_mc_eff, spec.predictors);

    fprintf('\nPredictors used:\n');
    disp(predictors_i(:));

    if numel(predictors_i) < numel(spec.predictors)
        missing_i = setdiff(spec.predictors, predictors_i, 'stable');
        fprintf('\nMissing predictors skipped:\n');
        disp(missing_i(:));
    end

    [MODEL_i, T_pred_i, T_metrics_i] = ...
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
            'Verbose', true);

    T_metrics_i.model_role = repmat(spec.role, height(T_metrics_i), 1);

    T_pred_i.predictor_set = repmat( ...
        strjoin(predictors_i, ", "), ...
        height(T_pred_i), ...
        1);
    T_pred_i.model_role = repmat(spec.role, height(T_pred_i), 1);

    ALL(sidx).name = spec.name;
    ALL(sidx).role = spec.role;
    ALL(sidx).predictors_requested = spec.predictors;
    ALL(sidx).predictors_used = predictors_i;
    ALL(sidx).MODEL = MODEL_i;
    ALL(sidx).T_pred = T_pred_i;
    ALL(sidx).T_metrics = T_metrics_i;

    T_all_metrics = [T_all_metrics; T_metrics_i]; %#ok<AGROW>
    T_pred_i_std = standardize_prediction_table_for_concat(T_pred_i);
    T_all_predictions = [T_all_predictions; T_pred_i_std]; %#ok<AGROW>

    spec_dir = fullfile(analysis_dir, char(spec.name));

    if ~exist(spec_dir, 'dir')
        mkdir(spec_dir);
    end

    adaptive_req.analysis.plot_q_model_performance( ...
        T_pred_i, ...
        T_metrics_i, ...
        'Split', 'test', ...
        'SaveFigure', true, ...
        'OutputDir', spec_dir, ...
        'BaseName', char(spec.name), ...
        'SavePNG', true, ...
        'SavePDF', true, ...
        'SaveFIG', false, ...
        'CloseAfterSave', true);

    writetable(T_metrics_i, fullfile(spec_dir, ...
        sprintf('%s_metrics.csv', spec.name)));

    writetable(remove_heavy_columns(T_pred_i), fullfile(spec_dir, ...
        sprintf('%s_predictions.csv', spec.name)));
end

%% Compare test metrics

T_test = T_all_metrics(T_all_metrics.split == "test", :);
T_test = sortrows(T_test, {'model_role', 'RMSE'}, {'ascend', 'ascend'});

T_operational = T_test(contains(T_test.model_role, "operational"), :);
T_operational = sortrows(T_operational, {'RMSE'}, {'ascend'});

T_diagnostic = T_test(~contains(T_test.model_role, "operational"), :);
T_diagnostic = sortrows(T_diagnostic, {'RMSE'}, {'ascend'});

T_q_by_M = summarize_prediction_errors(T_all_predictions, ...
    ["model_name", "model_type", "model_role", "REQ_M"]);

T_q_by_wave_model = summarize_prediction_errors(T_all_predictions, ...
    ["model_name", "model_type", "model_role", "SIM_WaveModel"]);

T_q_by_M_wave_model = summarize_prediction_errors(T_all_predictions, ...
    ["model_name", "model_type", "model_role", "REQ_M", "SIM_WaveModel"]);

fprintf('\n============================================================\n');
fprintf('OPERATIONAL TEST METRICS SORTED BY RMSE\n');
fprintf('============================================================\n');
disp(T_operational(:, {'model_name', 'model_type', 'model_role', ...
    'MAE', 'RMSE', 'bias', 'R2', 'spearman_rho'}));

fprintf('\n============================================================\n');
fprintf('DIAGNOSTIC TEST METRICS SORTED BY RMSE\n');
fprintf('============================================================\n');
disp(T_diagnostic(:, {'model_name', 'model_type', 'model_role', ...
    'MAE', 'RMSE', 'bias', 'R2', 'spearman_rho'}));

%% Save final outputs

T_all_predictions_light = remove_heavy_columns(T_all_predictions);

save(fullfile(analysis_dir, 'level_14_compare_advanced_features.mat'), ...
    'T_mc', ...
    'T_mc_eff', ...
    'MC', ...
    'PATHS', ...
    'COL', ...
    'model_specs', ...
    'model_types', ...
    'ALL', ...
    'T_all_metrics', ...
    'T_all_predictions_light', ...
    'T_test', ...
    'T_operational', ...
    'T_diagnostic', ...
    'T_q_by_M', ...
    'T_q_by_wave_model', ...
    'T_q_by_M_wave_model', ...
    '-v7.3');

writetable(T_all_metrics, ...
    fullfile(analysis_dir, 'all_model_metrics.csv'));

writetable(T_test, ...
    fullfile(analysis_dir, 'test_metrics_sorted_by_role_and_RMSE.csv'));

writetable(T_operational, ...
    fullfile(analysis_dir, 'operational_test_metrics_sorted_by_RMSE.csv'));

writetable(T_diagnostic, ...
    fullfile(analysis_dir, 'diagnostic_test_metrics_sorted_by_RMSE.csv'));

writetable(T_q_by_M, ...
    fullfile(analysis_dir, 'q_error_by_M.csv'));

writetable(T_q_by_wave_model, ...
    fullfile(analysis_dir, 'q_error_by_wave_model.csv'));

writetable(T_q_by_M_wave_model, ...
    fullfile(analysis_dir, 'q_error_by_M_and_wave_model.csv'));

writetable(T_all_predictions_light, ...
    fullfile(analysis_dir, 'all_model_predictions.csv'));

plot_operational_rmse_by_M(T_q_by_M, analysis_dir);

fprintf('\nSaved Level 14 comparison outputs to:\n%s\n', analysis_dir);
fprintf('\nLevel 14 advanced-feature comparison completed.\n');

%% Local helpers

function predictors_out = keep_existing_predictors(T, predictors_in)

vars = string(T.Properties.VariableNames);
predictors_in = string(predictors_in);

predictors_out = strings(0, 1);

for i = 1:numel(predictors_in)

    p_i = predictors_in(i);

    candidates = [
        p_i
        p_i + "_mean"
    ];

    if any(ismember(vars, candidates))
        predictors_out(end + 1, 1) = p_i; %#ok<AGROW>
    end
end

end

function T = remove_heavy_columns(T)

heavy_vars = {'req_curve', 'req_mapping', 'feat', 'feature_struct'};

for i = 1:numel(heavy_vars)
    if ismember(heavy_vars{i}, T.Properties.VariableNames)
        T.(heavy_vars{i}) = [];
    end
end

end

function T_std = standardize_prediction_table_for_concat(T_in)

n = height(T_in);

T_std = table();

string_vars = { ...
    'model_name', ...
    'model_type', ...
    'model_role', ...
    'predictor_set', ...
    'split', ...
    'condition_label', ...
    'SIM_WaveModel'};

numeric_vars = { ...
    'condition_id', ...
    'condition_position', ...
    'step_idx', ...
    'realization_idx', ...
    'patch_idx', ...
    'SIM_f0', ...
    'SIM_cs_bg', ...
    'REQ_M', ...
    'Omega_sr', ...
    'M_eff_guess', ...
    'M_eff_true_diag', ...
    'q_true', ...
    'q_pred_raw', ...
    'q_pred', ...
    'residual', ...
    'abs_error'};

for i = 1:numel(string_vars)
    name = string_vars{i};

    if ismember(name, T_in.Properties.VariableNames)
        T_std.(name) = string(T_in.(name));
    else
        T_std.(name) = strings(n, 1);
    end
end

for i = 1:numel(numeric_vars)
    name = numeric_vars{i};

    if ismember(name, T_in.Properties.VariableNames)
        T_std.(name) = double(T_in.(name));
    else
        T_std.(name) = NaN(n, 1);
    end
end

end

function T_summary = summarize_prediction_errors(T, group_vars)

T = T(T.split == "test", :);

[G, T_summary] = findgroups(T(:, cellstr(group_vars)));

T_summary.n = splitapply(@numel, T.residual, G);
T_summary.MAE_q = splitapply(@(x) mean(abs(x), 'omitnan'), T.residual, G);
T_summary.RMSE_q = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.residual, G);
T_summary.bias_q = splitapply(@(x) mean(x, 'omitnan'), T.residual, G);
T_summary.p95_abs_q_error = splitapply(@p95_abs, T.residual, G);

T_summary = sortrows(T_summary, cellstr(group_vars));

end

function value = p95_abs(x)

x = sort(abs(x(isfinite(x))));

if isempty(x)
    value = NaN;
    return;
end

idx = max(1, ceil(0.95 * numel(x)));
value = x(idx);

end

function plot_operational_rmse_by_M(T, output_dir)

T = T(T.model_type == "bagged_trees" & ...
    contains(T.model_role, "operational"), :);

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
title('Operational bagged-tree models by local window size');
legend(model_names, 'Interpreter', 'none', 'Location', 'best');
grid on;

exportgraphics(fig, fullfile(output_dir, ...
    'operational_bagged_tree_RMSE_by_M.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, ...
    'operational_bagged_tree_RMSE_by_M.pdf'), 'ContentType', 'vector');
close(fig);

end
