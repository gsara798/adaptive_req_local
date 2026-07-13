%% analyze_test_06_level05_interpret_associations.m
% Level 5 automatic interpretation of feature-q associations.
%
% This script reads the Level 4 association tables and prints a compact
% answer about whether each feature is associated with q.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load latest MC result to identify analysis folder

experiment_name = 'test_06_feature_q_baseline';

[~, ~, PATHS] = adaptive_req.analysis.load_mc_results( ...
    experiment_name, ...
    'RootDir', root_dir, ...
    'Verbose', true);

level04_dir = fullfile(PATHS.analysis_dir, 'level_04_feature_q_associations');

level04_file = fullfile(level04_dir, 'level_04_feature_q_associations.mat');

if ~exist(level04_file, 'file')
    error(['Level 4 association file not found:\n%s\n', ...
           'Run analyze_test_06_level04_feature_q_associations.m first.'], ...
           level04_file);
end

S = load(level04_file);

required_tables = { ...
    'T_assoc_global', ...
    'T_assoc_by_model_M', ...
    'T_assoc_controlled'};

for i = 1:numel(required_tables)
    if ~isfield(S, required_tables{i})
        error('Missing table in Level 4 file: %s', required_tables{i});
    end
end

%% Output folder

analysis_dir = fullfile(PATHS.analysis_dir, 'level_05_interpret_associations');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

%% Interpret global evidence

[T_feature_global, T_group_global] = ...
    adaptive_req.analysis.summarize_feature_q_evidence( ...
        S.T_assoc_global, ...
        'Label', 'global', ...
        'MinN', 8, ...
        'MinGroups', 1, ...
        'RequireSignificance', false, ...
        'Verbose', true);

%% Interpret evidence grouped by WaveModel and M

[T_feature_by_model_M, T_group_by_model_M] = ...
    adaptive_req.analysis.summarize_feature_q_evidence( ...
        S.T_assoc_by_model_M, ...
        'Label', 'grouped by WaveModel and M', ...
        'MinN', 8, ...
        'MinGroups', 4, ...
        'RequireSignificance', false, ...
        'Verbose', true);

%% Interpret controlled evidence

[T_feature_controlled, T_group_controlled] = ...
    adaptive_req.analysis.summarize_feature_q_evidence( ...
        S.T_assoc_controlled, ...
        'Label', 'controlled by WaveModel, M, f0, and cs_bg', ...
        'MinN', 6, ...
        'MinGroups', 10, ...
        'RequireSignificance', false, ...
        'Verbose', true);

%% Print final compact answer

fprintf('\n============================================================\n');
fprintf('FINAL AUTOMATIC INTERPRETATION\n');
fprintf('============================================================\n');

fprintf('\nMost important table: controlled evidence.\n');
fprintf('This asks whether features predict q after fixing WaveModel, M, f0, and cs_bg.\n\n');

disp(T_feature_controlled(:, { ...
    'feature_name', ...
    'n_valid_groups', ...
    'median_abs_spearman', ...
    'frac_moderate_or_better', ...
    'frac_strong', ...
    'direction_consistency', ...
    'summary_decision', ...
    'recommendation'}));

%% Save outputs

save(fullfile(analysis_dir, 'level_05_interpret_associations.mat'), ...
    'T_feature_global', ...
    'T_group_global', ...
    'T_feature_by_model_M', ...
    'T_group_by_model_M', ...
    'T_feature_controlled', ...
    'T_group_controlled', ...
    '-v7.3');

writetable(T_feature_global, ...
    fullfile(analysis_dir, 'feature_evidence_global.csv'));

writetable(T_group_global, ...
    fullfile(analysis_dir, 'group_evidence_global.csv'));

writetable(T_feature_by_model_M, ...
    fullfile(analysis_dir, 'feature_evidence_by_WaveModel_M.csv'));

writetable(T_group_by_model_M, ...
    fullfile(analysis_dir, 'group_evidence_by_WaveModel_M.csv'));

writetable(T_feature_controlled, ...
    fullfile(analysis_dir, 'feature_evidence_controlled.csv'));

writetable(T_group_controlled, ...
    fullfile(analysis_dir, 'group_evidence_controlled.csv'));

write_text_report( ...
    fullfile(analysis_dir, 'feature_q_evidence_report.txt'), ...
    T_feature_global, ...
    T_feature_by_model_M, ...
    T_feature_controlled);

fprintf('\nSaved Level 5 interpretation outputs to:\n%s\n', analysis_dir);

fprintf('\nLevel 5 automatic interpretation completed.\n');

%% Local helper

function write_text_report(file_path, T_global, T_model_M, T_controlled)

fid = fopen(file_path, 'w');

if fid < 0
    warning('Could not write report file: %s', file_path);
    return;
end

fprintf(fid, 'Feature-q evidence report\n');
fprintf(fid, '=========================\n\n');

fprintf(fid, 'Global evidence\n');
fprintf(fid, '---------------\n');
write_feature_table(fid, T_global);

fprintf(fid, '\nGrouped by WaveModel and M\n');
fprintf(fid, '--------------------------\n');
write_feature_table(fid, T_model_M);

fprintf(fid, '\nControlled by WaveModel, M, f0, and cs_bg\n');
fprintf(fid, '------------------------------------------\n');
write_feature_table(fid, T_controlled);

fprintf(fid, '\nInterpretation\n');
fprintf(fid, '--------------\n');
fprintf(fid, 'The controlled table is the most important for testing whether local features are associated with q after fixing major simulation and estimator parameters.\n');
fprintf(fid, 'Spearman rho is the main metric because it detects monotonic but non-linear associations.\n');
fprintf(fid, 'Linear R2 is included only as a diagnostic.\n');

fclose(fid);

end

function write_feature_table(fid, T)

for i = 1:height(T)

    fprintf(fid, '%s\n', string(T.feature_name(i)));
    fprintf(fid, '  decision: %s\n', string(T.summary_decision(i)));
    fprintf(fid, '  median |Spearman rho|: %.4f\n', T.median_abs_spearman(i));
    fprintf(fid, '  fraction moderate or better: %.4f\n', T.frac_moderate_or_better(i));
    fprintf(fid, '  fraction strong: %.4f\n', T.frac_strong(i));
    fprintf(fid, '  direction consistency: %.4f\n', T.direction_consistency(i));
    fprintf(fid, '  recommendation: %s\n\n', string(T.recommendation(i)));

end

end
