%% analyze_test_06_level08_model1_error_by_group.m
% Level 8: Analyze where Model 1 fails.
%
% This script analyzes Model 1 prediction errors by:
%
%   WaveModel
%   M
%   f0
%   cs_bg
%   condition_id
%
% The goal is to identify which known variables explain the remaining
% residual structure of the feature-only model.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load latest MC result to locate analysis folder

experiment_name = 'test_06_feature_q_baseline';

[~, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

level07_dir = fullfile(PATHS.analysis_dir, 'level_07_model1_features_only');
level07_file = fullfile(level07_dir, 'level_07_model1_features_only.mat');

if ~exist(level07_file, 'file')
    error(['Model 1 file not found:\n%s\n', ...
           'Run analyze_test_06_level07_train_model1_features_only.m first.'], ...
           level07_file);
end

S = load(level07_file);

T_pred_model1 = S.T_pred_model1;
T_metrics_model1 = S.T_metrics_model1;

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, 'level_08_model1_error_by_group');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Group variables to analyze

group_vars = [ ...
    "SIM_WaveModel", ...
    "REQ_M", ...
    "SIM_f0", ...
    "SIM_cs_bg", ...
    "condition_id"];

metrics_to_plot = ["RMSE", "MAE", "bias", "R2"];

T_all_group_errors = table();

%% Analyze one group variable at a time

for gidx = 1:numel(group_vars)

    group_var = group_vars(gidx);

    if ~ismember(group_var, string(T_pred_model1.Properties.VariableNames))
        fprintf('Skipping missing group variable: %s\n', group_var);
        continue;
    end

    fprintf('\nAnalyzing errors by %s\n', group_var);

    T_err_g = adaptive_req.analysis.summarize_q_model_errors_by_group( ...
        T_pred_model1, ...
        'GroupVars', group_var, ...
        'Split', 'test');

    % Save the original table with the real grouping column.
    writetable(T_err_g, fullfile(analysis_dir, ...
        sprintf('model1_errors_by_%s.csv', group_var)));
    
    % Standardize the table so it can be concatenated with other group summaries.
    T_err_g_std = standardize_group_error_table( ...
        T_err_g, ...
        group_var, ...
        group_var);
    
    T_all_group_errors = [T_all_group_errors; T_err_g_std]; %#ok<AGROW>

    for midx = 1:numel(metrics_to_plot)

        metric_i = metrics_to_plot(midx);

        adaptive_req.analysis.plot_q_model_error_by_group( ...
            T_err_g, ...
            group_var, ...
            'Metric', metric_i, ...
            'Title', sprintf('Model 1 %s by %s', metric_i, group_var), ...
            'SaveFigure', true, ...
            'OutputDir', analysis_dir, ...
            'BaseName', sprintf('model1_%s_by_%s', metric_i, group_var), ...
            'SavePNG', true, ...
            'SavePDF', true, ...
            'SaveFIG', false, ...
            'CloseAfterSave', false);

    end
end

%% Analyze combined groups

combined_group_sets = { ...
    ["SIM_WaveModel", "REQ_M"], ...
    ["SIM_WaveModel", "SIM_f0"], ...
    ["SIM_WaveModel", "SIM_cs_bg"], ...
    ["REQ_M", "SIM_f0"], ...
    ["REQ_M", "SIM_cs_bg"]};

T_combined_errors = table();

for sidx = 1:numel(combined_group_sets)

    group_set = combined_group_sets{sidx};

    valid_group_set = group_set(ismember(group_set, ...
        string(T_pred_model1.Properties.VariableNames)));

    if numel(valid_group_set) ~= numel(group_set)
        continue;
    end

    grouping_label = strjoin(group_set, "_");

    fprintf('\nAnalyzing combined errors by %s\n', grouping_label);

    T_err_s = adaptive_req.analysis.summarize_q_model_errors_by_group( ...
        T_pred_model1, ...
        'GroupVars', group_set, ...
        'Split', 'test');

    % Save the original combined-group table.
    writetable(T_err_s, fullfile(analysis_dir, ...
        sprintf('model1_errors_by_%s.csv', grouping_label)));
    
    % Standardize before concatenation.
    T_err_s_std = standardize_group_error_table( ...
        T_err_s, ...
        group_set, ...
        grouping_label);
    
    T_combined_errors = [T_combined_errors; T_err_s_std]; %#ok<AGROW>

end

%% Print compact summary

fprintf('\nOverall Model 1 test metrics.\n');
disp(T_metrics_model1(T_metrics_model1.split == "test", :));

fprintf('\nError by group summary, sorted by RMSE.\n');

T_sorted = sortrows(T_all_group_errors, 'RMSE', 'descend');

display_vars = intersect( ...
    {'grouping_name', 'group_value', 'model_type', 'split', ...
     'n', 'RMSE', 'MAE', 'bias', 'R2', 'spearman_rho'}, ...
    T_sorted.Properties.VariableNames, ...
    'stable');

disp(T_sorted(:, display_vars));

%% Save outputs

save(fullfile(analysis_dir, 'level_08_model1_error_by_group.mat'), ...
    'T_pred_model1', ...
    'T_metrics_model1', ...
    'T_all_group_errors', ...
    'T_combined_errors', ...
    '-v7.3');

writetable(T_all_group_errors, ...
    fullfile(analysis_dir, 'model1_errors_all_single_groups.csv'));

if ~isempty(T_combined_errors)
    writetable(T_combined_errors, ...
        fullfile(analysis_dir, 'model1_errors_combined_groups.csv'));
end

fprintf('\nSaved Level 8 Model 1 error analysis to:\n%s\n', analysis_dir);

fprintf('\nLevel 8 Model 1 error-by-group analysis completed.\n');

%%
function T_std = standardize_group_error_table(T_in, group_vars, grouping_name)
%STANDARDIZE_GROUP_ERROR_TABLE Convert arbitrary group columns to a common format.
%
% This allows tables grouped by different variables to be vertically
% concatenated.

group_vars = string(group_vars);
grouping_name = string(grouping_name);

available_vars = string(T_in.Properties.VariableNames);

group_vars = group_vars(ismember(group_vars, available_vars));

if isempty(group_vars)
    error('No valid group variables found for standardization.');
end

group_value = strings(height(T_in), 1);

for r = 1:height(T_in)

    parts = strings(1, numel(group_vars));

    for g = 1:numel(group_vars)

        var_g = group_vars(g);
        value_g = T_in.(char(var_g))(r);

        parts(g) = var_g + "=" + value_to_string(value_g);

    end

    group_value(r) = strjoin(parts, " | ");

end

T_std = T_in;

% Remove original group-specific columns so all standardized tables share
% the same variable names.
T_std(:, cellstr(group_vars)) = [];

T_std.grouping_name = repmat(grouping_name, height(T_std), 1);
T_std.group_value = group_value;

T_std = movevars(T_std, {'grouping_name', 'group_value'}, 'Before', 1);

end

function txt = value_to_string(value)

if iscell(value)
    value = value{1};
end

if isnumeric(value) || islogical(value)
    txt = string(sprintf('%.6g', value));
elseif isstring(value)
    txt = value;
elseif ischar(value)
    txt = string(value);
elseif iscategorical(value)
    txt = string(value);
else
    txt = string(value);
end

end