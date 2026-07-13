%% run_test_08_advanced_angular_features.m
% Generate the Test 08 dataset with advanced angular spectral features.
%
% This experiment preserves the Test 07 physical sweep but uses the current
% feature extractor, which includes angular moments and lobe descriptors.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
experiments_dir = fileparts(this_file);
root_dir = fileparts(experiments_dir);

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Run complete Test 08 Monte Carlo sweep

[T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
    'test_08_advanced_angular_features', ...
    'RootDir', root_dir, ...
    'MaxConditions', Inf, ...
    'SaveResults', true, ...
    'SaveIntermediate', false, ...
    'ContinueOnError', false, ...
    'Verbose', true);

%% Validate and summarize results

assert(all(MC.T_status.success), ...
    'Test 08 dataset generation completed with failed conditions.');

assert(height(T_mc) == MC.expected_rows_selected, ...
    'Generated row count does not match the expected row count.');

assert(ismember('req_mapping', T_mc.Properties.VariableNames), ...
    'The generated dataset does not contain req_mapping.');

assert(~ismember('req_curve', T_mc.Properties.VariableNames), ...
    'The generated dataset unexpectedly contains the heavy req_curve.');

required_features = { ...
    'ang_moment_1', ...
    'ang_moment_2', ...
    'ang_moment_4', ...
    'ang_peak_count_rel', ...
    'ang_top2_to_top1', ...
    'ang_peak_separation_deg'};

missing_features = setdiff(required_features, T_mc.Properties.VariableNames);

assert(isempty(missing_features), ...
    'Missing advanced angular features: %s', strjoin(missing_features, ', '));

fprintf('\nTest 08 dataset completed successfully.\n');
fprintf('Conditions generated: %d\n', MC.n_conditions_selected);
fprintf('Rows generated      : %d\n', height(T_mc));
fprintf('Output folder:\n%s\n', MC.SAVE.output_dir);
