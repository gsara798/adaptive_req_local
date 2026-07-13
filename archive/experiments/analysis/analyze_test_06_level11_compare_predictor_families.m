%% analyze_test_06_level11_compare_predictor_families.m
% Level 11: Compare predictor families and model types.
%
% This script compares:
%
%   Model A: local entropy features only
%   Model B: local entropy features + known parameters
%   Model C: rich local spectral features + known parameters
%   Model D: diagnostic model with M_eff
%
% Model D is diagnostic because M_eff uses cs_true.

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

%% Add diagnostic effective-window variable

cs_guess = 3.0;

T_mc_eff = adaptive_req.analysis.add_effective_window_metrics( ...
    T_mc, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', cs_guess);

% Short alias for convenience.
if ~ismember("M_eff", string(T_mc_eff.Properties.VariableNames))
    T_mc_eff.M_eff = T_mc_eff.M_eff_true_diag;
end

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_11_compare_predictor_families');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Define predictor families

model_specs = struct([]);

model_specs(1).name = "ModelA_local_entropy_only";
model_specs(1).predictors = [
    "ang_entropy"
    "radial_entropy"
];

model_specs(2).name = "ModelB_entropy_known_params";
model_specs(2).predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"

];

model_specs(3).name = "ModelC_rich_spectral_known_params";
model_specs(3).predictors = [
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
    "ang_moment_1"
    "ang_moment_2"
    "ang_moment_4"
    "ang_peak_count_rel"
    "ang_top1_window_frac"
    "ang_top2_window_frac"
    "ang_top3_window_frac"
    "ang_top2_to_top1"
    "ang_peak_separation_deg"


];

model_specs(4).name = "ModelD_diagnostic_with_Meff";
model_specs(4).predictors = [
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
    "M_eff"
];

%% Model types to compare

model_types = [
    "linear"
    "quadratic"
    "boosted_trees"
    "bagged_trees"
];

%% Train all models

ALL = struct();

T_all_metrics = table();
T_all_predictions = table();

for sidx = 1:numel(model_specs)

    spec = model_specs(sidx);

    fprintf('\n============================================================\n');
    fprintf('Training %s\n', spec.name);
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
            'KNNK', 15, ...
            'NumLearningCycles', 300, ...
            'MinLeafSize', 8, ...
            'LearnRate', 0.05, ...
            'ClipPredictions', true, ...
            'Verbose', true);

    ALL(sidx).name = spec.name;
    ALL(sidx).predictors_requested = spec.predictors;
    ALL(sidx).predictors_used = predictors_i;
    ALL(sidx).MODEL = MODEL_i;
    ALL(sidx).T_pred = T_pred_i;
    ALL(sidx).T_metrics = T_metrics_i;

    T_all_metrics = [T_all_metrics; T_metrics_i]; %#ok<AGROW>
    T_pred_i.predictor_set = repmat( ...
    strjoin(predictors_i, ", "), ...
    height(T_pred_i), ...
    1);

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
        'CloseAfterSave', false);

    writetable(T_metrics_i, fullfile(spec_dir, ...
        sprintf('%s_metrics.csv', spec.name)));

    writetable(T_pred_i, fullfile(spec_dir, ...
        sprintf('%s_predictions.csv', spec.name)));
end

%% Compare test metrics

T_test = T_all_metrics(T_all_metrics.split == "test", :);

T_test = sortrows(T_test, {'RMSE'}, {'ascend'});

fprintf('\n============================================================\n');
fprintf('TEST METRICS SORTED BY RMSE\n');
fprintf('============================================================\n');

disp(T_test(:, {'model_name', 'model_type', 'split', ...
    'MAE', 'RMSE', 'bias', 'R2', 'spearman_rho'}));

%% Save final outputs

save(fullfile(analysis_dir, 'level_11_compare_predictor_families.mat'), ...
    'T_mc', ...
    'T_mc_eff', ...
    'MC', ...
    'PATHS', ...
    'COL', ...
    'model_specs', ...
    'model_types', ...
    'ALL', ...
    'T_all_metrics', ...
    'T_all_predictions', ...
    'T_test', ...
    '-v7.3');

writetable(T_all_metrics, ...
    fullfile(analysis_dir, 'all_model_metrics.csv'));

writetable(T_test, ...
    fullfile(analysis_dir, 'test_metrics_sorted_by_RMSE.csv'));

writetable(T_all_predictions, ...
    fullfile(analysis_dir, 'all_model_predictions.csv'));

fprintf('\nSaved Level 11 comparison outputs to:\n%s\n', analysis_dir);

fprintf('\nLevel 11 predictor-family comparison completed.\n');

%% Local helper

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

function T_std = standardize_prediction_table_for_concat(T_in)
%STANDARDIZE_PREDICTION_TABLE_FOR_CONCAT Create common prediction table.
%
% Different predictor families generate prediction tables with different
% predictor columns. This function keeps only a fixed set of common metadata,
% prediction, and error columns so tables can be vertically concatenated.

n = height(T_in);

T_std = table();

T_std.model_name = get_string_column(T_in, 'model_name', n);
T_std.model_type = get_string_column(T_in, 'model_type', n);
T_std.predictor_set = get_string_column(T_in, 'predictor_set', n);
T_std.split = get_string_column(T_in, 'split', n);

T_std.condition_id = get_numeric_column(T_in, 'condition_id', n);
T_std.condition_position = get_numeric_column(T_in, 'condition_position', n);
T_std.condition_label = get_string_column(T_in, 'condition_label', n);

T_std.step_idx = get_numeric_column(T_in, 'step_idx', n);
T_std.realization_idx = get_numeric_column(T_in, 'realization_idx', n);
T_std.patch_idx = get_numeric_column(T_in, 'patch_idx', n);

T_std.SIM_WaveModel = get_string_column(T_in, 'SIM_WaveModel', n);
T_std.SIM_f0 = get_numeric_column(T_in, 'SIM_f0', n);
T_std.SIM_cs_bg = get_numeric_column(T_in, 'SIM_cs_bg', n);
T_std.REQ_M = get_numeric_column(T_in, 'REQ_M', n);

T_std.Omega_sr = get_numeric_column(T_in, 'Omega_sr', n);
T_std.omega_mean = get_numeric_column(T_in, 'omega_mean', n);
T_std.aperture_value = get_numeric_column(T_in, 'aperture_value', n);

T_std.ang_entropy = get_numeric_column(T_in, 'ang_entropy', n);
T_std.radial_entropy = get_numeric_column(T_in, 'radial_entropy', n);

T_std.M_eff = get_numeric_column(T_in, 'M_eff', n);
T_std.M_eff_true_diag = get_numeric_column(T_in, 'M_eff_true_diag', n);

T_std.q_true = get_numeric_column(T_in, 'q_true', n);
T_std.q_pred_raw = get_numeric_column(T_in, 'q_pred_raw', n);
T_std.q_pred = get_numeric_column(T_in, 'q_pred', n);
T_std.residual = get_numeric_column(T_in, 'residual', n);
T_std.abs_error = get_numeric_column(T_in, 'abs_error', n);

end

function x = get_numeric_column(T, var_name, n)

if ismember(var_name, T.Properties.VariableNames)
    x = T.(var_name);

    if iscell(x)
        x = cellfun(@double, x);
    end

    x = double(x(:));
else
    x = NaN(n, 1);
end

end

function x = get_string_column(T, var_name, n)

if ismember(var_name, T.Properties.VariableNames)
    x = string(T.(var_name));
    x = x(:);
else
    x = strings(n, 1);
end

end