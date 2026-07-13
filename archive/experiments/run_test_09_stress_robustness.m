%% run_test_09_stress_robustness.m
% Generate Test 09 stress-screening dataset.
%
% This is intentionally separate from Test 08. It tests robustness across
% SNR, number of waves, wave model, M, and ForceInPlaneWave.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
experiments_dir = fileparts(this_file);
root_dir = fileparts(experiments_dir);

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Run stress sweep

[T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
    'test_09_stress_robustness', ...
    'RootDir', root_dir, ...
    'MaxConditions', Inf, ...
    'SaveResults', true, ...
    'SaveIntermediate', false, ...
    'ContinueOnError', false, ...
    'Verbose', true);

%% Validate

assert(all(MC.T_status.success), ...
    'Test 09 dataset generation completed with failed conditions.');
assert(height(T_mc) == MC.expected_rows_selected, ...
    'Generated row count does not match the expected row count.');
assert(ismember('req_mapping', T_mc.Properties.VariableNames), ...
    'The generated dataset does not contain req_mapping.');

fprintf('\nTest 09 stress robustness dataset completed successfully.\n');
fprintf('Conditions generated: %d\n', MC.n_conditions_selected);
fprintf('Rows generated      : %d\n', height(T_mc));
fprintf('Output folder:\n%s\n', MC.SAVE.output_dir);
