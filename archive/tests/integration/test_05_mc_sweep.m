%% test_05_mc_sweep.m
% Run a small MC parameter sweep using test_07_feature_q_baseline.
%
% For testing, only the first two parameter conditions are run.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Run small MC sweep

[T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
    'test_07_feature_q_baseline', ...
    'RootDir', root_dir, ...
    'MaxConditions', 2, ...
    'SaveResults', false, ...
    'SaveIntermediate', false, ...
    'ContinueOnError', false, ...
    'Verbose', true);

%% Display summary

fprintf('\nTest 05 MC sweep summary.\n');
fprintf('Conditions selected = %d\n', MC.n_conditions_selected);
fprintf('Rows generated = %d\n', height(T_mc));

disp(MC.T_status);

fprintf('\nTest 05 MC sweep completed.\n');
