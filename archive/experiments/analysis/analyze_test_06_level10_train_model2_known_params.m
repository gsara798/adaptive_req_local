%% analyze_test_06_level10_train_model2_known_params.m
% Level 10: Train Model 2.
%
% Model 2:
%
%   q_local = F(local features, known experimental or estimator parameters)
%
% Default predictors:
%
%   ang_entropy
%   radial_entropy
%   REQ_M
%   SIM_f0
%
% This script controls which predictors enter the model. The training
% function is generic.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load MC results

experiment_name = 'test_06_feature_q_baseline';

[T_mc, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

%% Detect columns

COL = adaptive_req.analysis.detect_mc_columns(T_mc);

%% Choose Model 2 predictors here

model2_predictors = [ ...
    "ang_entropy", ...
    "radial_entropy", ...
    "REQ_M", ...
    "SIM_f0"];

% Examples for later:
%
model2_predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
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

model2_predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
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

% Diagnostic only, not experimental as final model:
%
diagnostic_predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
    "SIM_cs_bg"
];

fprintf('\nModel 2 predictors:\n');
disp(model2_predictors(:));

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, 'level_10_model2_known_params');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Train Model 2

[MODEL2, T_pred_model2, T_metrics_model2] = ...
    adaptive_req.analysis.train_q_model_from_predictors( ...
        T_mc, ...
        model2_predictors, ...
        'QVar', COL.q, ...
        'ModelName', 'Model2_features_M_f0', ...
        'SplitMode', 'condition', ...
        'ConditionVar', COL.condition, ...
        'TrainFraction', 0.70, ...
        'RandomSeed', 1701, ...
        'ModelTypes', ["linear", "quadratic", "knn"], ...
        'KNNK', 15, ...
        'ClipPredictions', true, ...
        'Verbose', true);

%% Plot performance

FIGS = adaptive_req.analysis.plot_q_model_performance( ...
    T_pred_model2, ...
    T_metrics_model2, ...
    'Split', 'test', ...
    'SaveFigure', true, ...
    'OutputDir', analysis_dir, ...
    'BaseName', 'model2_features_M_f0', ...
    'SavePNG', true, ...
    'SavePDF', true, ...
    'SaveFIG', false, ...
    'CloseAfterSave', false);

%% Test metrics

fprintf('\nModel 2 test metrics.\n');

T_test_metrics_model2 = T_metrics_model2(T_metrics_model2.split == "test", :);

disp(T_test_metrics_model2);

%% Optional comparison with Model 1

level07_dir = fullfile(PATHS.analysis_dir, 'level_07_model1_features_only');
level07_file = fullfile(level07_dir, 'level_07_model1_features_only.mat');

T_compare = table();

if exist(level07_file, 'file')

    S1 = load(level07_file, 'T_metrics_model1');

    T_model1 = S1.T_metrics_model1(S1.T_metrics_model1.split == "test", :);
    T_model2 = T_test_metrics_model2;

    if ~ismember('model_name', T_model1.Properties.VariableNames)
        T_model1.model_name = repmat("Model1_features_only", height(T_model1), 1);
        T_model1 = movevars(T_model1, 'model_name', 'Before', 1);
    end

    T_compare = [T_model1; T_model2];

    fprintf('\nModel 1 versus Model 2 test metrics.\n');
    disp(T_compare(:, {'model_name', 'model_type', 'split', ...
        'MAE', 'RMSE', 'bias', 'R2', 'spearman_rho'}));

    writetable(T_compare, fullfile(analysis_dir, ...
        'model1_vs_model2_test_metrics.csv'));

end

%% Save outputs

save(fullfile(analysis_dir, 'level_10_model2_known_params.mat'), ...
    'T_mc', ...
    'MC', ...
    'PATHS', ...
    'COL', ...
    'model2_predictors', ...
    'MODEL2', ...
    'T_pred_model2', ...
    'T_metrics_model2', ...
    'T_test_metrics_model2', ...
    'T_compare', ...
    'FIGS', ...
    '-v7.3');

writetable(T_metrics_model2, ...
    fullfile(analysis_dir, 'model2_metrics.csv'));

writetable(T_test_metrics_model2, ...
    fullfile(analysis_dir, 'model2_test_metrics.csv'));

writetable(T_pred_model2, ...
    fullfile(analysis_dir, 'model2_predictions.csv'));

fprintf('\nSaved Model 2 outputs to:\n%s\n', analysis_dir);

fprintf('\nLevel 10 Model 2 training completed.\n');
