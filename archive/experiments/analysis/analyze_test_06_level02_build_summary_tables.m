%% analyze_test_06_level02_build_summary_tables.m
% Build reusable summary tables for test_06_feature_q_baseline.
%
% This script:
%   1. Loads the latest MC result.
%   2. Detects key columns.
%   3. Detects candidate feature columns.
%   4. Builds a condition-step summary table.
%   5. Builds a condition-level summary table.
%   6. Saves all analysis tables.

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

fprintf('\nDetected columns.\n');
disp(COL);

%% Detect feature variables

feature_vars = adaptive_req.analysis.detect_feature_variables(T_mc);

fprintf('\nDetected feature variables.\n');

if isempty(feature_vars)
    fprintf('No feature variables were detected automatically.\n');
else
    disp(feature_vars(:));
end

%% Build summary tables

T_step = adaptive_req.analysis.make_step_summary_table( ...
    T_mc, ...
    COL, ...
    feature_vars);

T_condition = adaptive_req.analysis.make_condition_summary_table(T_step);

fprintf('\nGenerated summary tables.\n');
fprintf('T_step rows      : %d\n', height(T_step));
fprintf('T_condition rows : %d\n', height(T_condition));

fprintf('\nFirst rows of T_step.\n');
disp(head(T_step, min(10, height(T_step))));

fprintf('\nFirst rows of T_condition.\n');
disp(head(T_condition, min(10, height(T_condition))));

%% Save analysis outputs

analysis_dir = fullfile(PATHS.analysis_dir, 'level_02_summary_tables');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

save(fullfile(analysis_dir, 'level_02_summary_tables.mat'), ...
    'T_mc', ...
    'MC', ...
    'PATHS', ...
    'COL', ...
    'feature_vars', ...
    'T_step', ...
    'T_condition', ...
    '-v7.3');

writetable(T_step, fullfile(analysis_dir, 'step_summary_table.csv'));
writetable(T_condition, fullfile(analysis_dir, 'condition_summary_table.csv'));

T_feature_vars = table(feature_vars(:), 'VariableNames', {'feature_variable'});
writetable(T_feature_vars, fullfile(analysis_dir, 'detected_feature_variables.csv'));

fprintf('\nSaved Level 2 summary tables to:\n%s\n', analysis_dir);

fprintf('\nLevel 2 summary table generation completed.\n');
