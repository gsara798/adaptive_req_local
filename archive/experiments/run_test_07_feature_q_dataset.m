%% run_test_07_feature_q_dataset.m
% Generate the complete Test 07 dataset for feature-to-q ML experiments.
%
% The resulting table stores the lightweight req_mapping required to
% convert predicted quantiles into predicted wavenumber and SWS.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
experiments_dir = fileparts(this_file);
root_dir = fileparts(experiments_dir);

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Run complete Test 07 Monte Carlo sweep

[T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
    'test_07_feature_q_baseline', ...
    'RootDir', root_dir, ...
    'MaxConditions', Inf, ...
    'SaveResults', true, ...
    'SaveIntermediate', false, ...
    'ContinueOnError', false, ...
    'Verbose', true);

%% Validate and summarize results

assert(all(MC.T_status.success), ...
    'Test 07 dataset generation completed with failed conditions.');

assert(height(T_mc) == MC.expected_rows_selected, ...
    'Generated row count does not match the expected row count.');

assert(ismember('req_mapping', T_mc.Properties.VariableNames), ...
    'The generated dataset does not contain req_mapping.');

assert(~ismember('req_curve', T_mc.Properties.VariableNames), ...
    'The generated dataset unexpectedly contains the heavy req_curve.');

fprintf('\nTest 07 dataset completed successfully.\n');
fprintf('Conditions generated: %d\n', MC.n_conditions_selected);
fprintf('Rows generated      : %d\n', height(T_mc));
fprintf('Output folder:\n%s\n', MC.SAVE.output_dir);
