%% export_test_12_level01_deployment_models.m
% Export lightweight deployment q-model files from the existing Test 12
% Level 01 all-model MAT file.
%
% This does not retrain. It only repackages trained models so later map-level
% analyses can load one model directly instead of opening a multi-GB MAT.

clear; clc;
format compact;

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

[~, ~, PATHS12] = adaptive_req.analysis.load_mc_results( ...
    'test_12_cs_guess_window_sweep', ...
    'RootDir', root_dir, ...
    'Verbose', true);

level01_dir = fullfile(PATHS12.analysis_dir, 'level_01_model_comparison');
model_file = fullfile(level01_dir, 'models', ...
    'level12_level01_model_comparison_models.mat');
deploy_dir = fullfile(level01_dir, 'models', 'deployment');

if ~exist(model_file, 'file')
    error('Level 01 all-model MAT file not found:\n%s', model_file);
end
if ~exist(deploy_dir, 'dir')
    mkdir(deploy_dir);
end

fprintf('\nLoading existing Level 01 all-model file. This may take a while.\n%s\n', ...
    model_file);
S = load(model_file, 'MODELS');
MODELS = S.MODELS;

files = strings(0, 1);
for i = 1:numel(MODELS)
    fprintf('\nExporting deployment bundle %d / %d: %s | %s\n', ...
        i, numel(MODELS), string(MODELS(i).model_name), ...
        string(MODELS(i).feature_set));

    files_i = adaptive_req.analysis.save_q_model_deployment( ...
        MODELS(i).model, deploy_dir, ...
        'ModelName', MODELS(i).model_name, ...
        'FeatureSet', MODELS(i).feature_set, ...
        'ModelRole', MODELS(i).model_role, ...
        'ModelTypes', ["linear", "boosted_trees", "bagged_trees"], ...
        'Overwrite', true);

    files = [files; files_i(:)]; %#ok<AGROW>
end

T_files = table(files, 'VariableNames', {'deployment_model_file'});
writetable(T_files, fullfile(deploy_dir, ...
    'level12_level01_deployment_model_files.csv'));

fprintf('\nExport complete.\nDeployment directory:\n%s\n', deploy_dir);
fprintf('Files exported: %d\n', numel(files));
