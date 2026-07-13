%% analyze_test_06_level01_sanity_check.m
% Level 1 sanity check for test_06_feature_q_baseline.
%
% This script checks:
%   1. Dataset size.
%   2. Sweep parameter coverage.
%   3. Row counts per condition.
%   4. Row counts per aperture step.
%   5. Missing values.
%   6. Basic q and Omega ranges.
%   7. Basic feature ranges.
%   8. Within-step variability of q.
%
% The goal is not to interpret the physics yet.
% The goal is to verify that the generated dataset is complete and usable.

clear; clc; close all;
format compact;

%% ========================================================================
% 1. Project setup
% ========================================================================

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% ========================================================================
% 2. Locate and load latest MC result
% ========================================================================

experiment_name = 'test_06_feature_q_baseline';

output_root = fullfile(root_dir, 'outputs', experiment_name);

if ~exist(output_root, 'dir')
    error('Output folder not found: %s', output_root);
end

run_dir = find_latest_run_dir(output_root);

fprintf('\nUsing run folder:\n%s\n', run_dir);

data_dir = fullfile(run_dir, 'data');

mat_file = find_mc_result_file(data_dir);

fprintf('\nLoading MC result:\n%s\n', mat_file);

S = load(mat_file);

if ~isfield(S, 'T_mc')
    error('The selected MAT file does not contain T_mc.');
end

T_mc = S.T_mc;

if isfield(S, 'MC')
    MC = S.MC;
else
    MC = struct();
end

fprintf('\nLoaded table.\n');
fprintf('Rows    : %d\n', height(T_mc));
fprintf('Columns : %d\n', width(T_mc));

%% ========================================================================
% 3. Identify important columns
% ========================================================================

vars = string(T_mc.Properties.VariableNames);

q_var = find_first_existing_var(vars, { ...
    'q_theory', ...
    'q', ...
    'q_val', ...
    'q_reference'});

omega_var = find_first_existing_var(vars, { ...
    'Omega_sr', ...
    'Omega', ...
    'omega_sr', ...
    'solid_angle', ...
    'solid_angle_sr'});

condition_var = find_first_existing_var(vars, { ...
    'condition_id'});

step_var = find_first_existing_var(vars, { ...
    'step_idx', ...
    'step_id', ...
    'aperture_step', ...
    'step'});

realization_var = find_first_existing_var(vars, { ...
    'realization_idx', ...
    'realization_id', ...
    'rep_idx', ...
    'rep', ...
    'irep'});

patch_var = find_first_existing_var(vars, { ...
    'patch_idx', ...
    'patch_id', ...
    'patch'});

fprintf('\nDetected key columns.\n');
fprintf('q column          : %s\n', q_var);
fprintf('Omega column      : %s\n', omega_var);
fprintf('condition column  : %s\n', condition_var);
fprintf('step column       : %s\n', step_var);
fprintf('realization column: %s\n', realization_var);
fprintf('patch column      : %s\n', patch_var);

if q_var == ""
    error('Could not detect q column. Check T_mc.Properties.VariableNames.');
end

if omega_var == ""
    warning('Could not detect Omega column. Some checks will be skipped.');
end

%% ========================================================================
% 4. Expected size check
% ========================================================================

expected_conditions = NaN;
expected_rows = NaN;
rows_per_condition = NaN;

if isfield(MC, 'n_conditions_selected')
    expected_conditions = MC.n_conditions_selected;
end

if isfield(MC, 'expected_rows_selected')
    expected_rows = MC.expected_rows_selected;
end

if isfield(MC, 'rows_per_condition')
    rows_per_condition = MC.rows_per_condition;
end

fprintf('\nDataset size check.\n');

if ~isnan(expected_conditions)
    fprintf('Expected selected conditions: %d\n', expected_conditions);
end

if ~isnan(rows_per_condition)
    fprintf('Expected rows per condition : %d\n', rows_per_condition);
end

if ~isnan(expected_rows)
    fprintf('Expected total rows         : %d\n', expected_rows);
    fprintf('Actual total rows           : %d\n', height(T_mc));

    if height(T_mc) == expected_rows
        fprintf('Status: PASS. Total row count matches expected rows.\n');
    else
        fprintf('Status: WARNING. Total row count does not match expected rows.\n');
    end
else
    fprintf('Expected total rows not available in MC metadata.\n');
end

%% ========================================================================
% 5. Sweep coverage check
% ========================================================================

fprintf('\nSweep parameter coverage.\n');

sweep_vars = vars(startsWith(vars, "SIM_") | startsWith(vars, "REQ_"));

for i = 1:numel(sweep_vars)

    v = sweep_vars(i);
    values_i = unique(T_mc.(v));

    fprintf('\n%s:\n', v);
    disp(values_i);

end

if isfield(MC, 'T_conditions')
    fprintf('\nExpected condition matrix from MC.T_conditions:\n');
    disp(MC.T_conditions);
end

%% ========================================================================
% 6. Row count per condition
% ========================================================================

if condition_var ~= ""

    T_condition_counts = count_by_group(T_mc, condition_var);

    fprintf('\nRows per condition.\n');
    disp(T_condition_counts);

    if ~isnan(rows_per_condition)

        bad_condition_rows = T_condition_counts.GroupCount ~= rows_per_condition;

        if any(bad_condition_rows)
            fprintf('\nWARNING. Some conditions do not have the expected number of rows.\n');
            disp(T_condition_counts(bad_condition_rows, :));
        else
            fprintf('\nPASS. All conditions have the expected number of rows.\n');
        end
    end

else

    T_condition_counts = table();
    warning('Skipping condition row count because condition_id was not found.');

end

%% ========================================================================
% 7. Row count per condition and step
% ========================================================================

if condition_var ~= "" && step_var ~= ""

    T_step_counts = count_by_group(T_mc, [condition_var, step_var]);

    fprintf('\nRows per condition and step.\n');
    disp(head(T_step_counts, min(20, height(T_step_counts))));

    if realization_var ~= "" && patch_var ~= ""

        n_realizations = numel(unique(T_mc.(realization_var)));
        n_patches = numel(unique(T_mc.(patch_var)));

        expected_rows_per_step = n_realizations * n_patches;

        fprintf('\nDetected realizations: %d\n', n_realizations);
        fprintf('Detected patches     : %d\n', n_patches);
        fprintf('Expected rows per condition-step: %d\n', expected_rows_per_step);

        bad_step_rows = T_step_counts.GroupCount ~= expected_rows_per_step;

        if any(bad_step_rows)
            fprintf('\nWARNING. Some condition-step groups do not have expected rows.\n');
            disp(T_step_counts(bad_step_rows, :));
        else
            fprintf('\nPASS. All condition-step groups have expected rows.\n');
        end
    end

else

    T_step_counts = table();
    warning('Skipping condition-step row count because condition or step columns were not found.');

end

%% ========================================================================
% 8. Missing values and non-finite values
% ========================================================================

T_missing = summarize_missing_values(T_mc);

fprintf('\nVariables with missing or non-finite values.\n');

bad_missing = T_missing.n_missing > 0 | T_missing.n_nonfinite > 0;

if any(bad_missing)
    disp(T_missing(bad_missing, :));
else
    fprintf('PASS. No missing or non-finite values detected in numeric/string columns.\n');
end

%% ========================================================================
% 9. q and Omega sanity check
% ========================================================================

q = T_mc.(q_var);

fprintf('\nq sanity check.\n');
fprintf('q min    = %.6g\n', min(q, [], 'omitnan'));
fprintf('q max    = %.6g\n', max(q, [], 'omitnan'));
fprintf('q mean   = %.6g\n', mean(q, 'omitnan'));
fprintf('q median = %.6g\n', median(q, 'omitnan'));
fprintf('q std    = %.6g\n', std(q, 'omitnan'));

if any(q < 0 | q > 1)
    warning('Some q values are outside [0, 1].');
else
    fprintf('PASS. All q values are inside [0, 1].\n');
end

if omega_var ~= ""

    omega = T_mc.(omega_var);

    fprintf('\nOmega sanity check.\n');
    fprintf('Omega min    = %.6g\n', min(omega, [], 'omitnan'));
    fprintf('Omega max    = %.6g\n', max(omega, [], 'omitnan'));
    fprintf('Omega mean   = %.6g\n', mean(omega, 'omitnan'));
    fprintf('Omega median = %.6g\n', median(omega, 'omitnan'));
    fprintf('Omega std    = %.6g\n', std(omega, 'omitnan'));

end

%% ========================================================================
% 10. Feature column detection and summary
% ========================================================================

feature_vars = detect_feature_variables(T_mc);

fprintf('\nDetected candidate feature columns.\n');

if isempty(feature_vars)
    fprintf('No feature columns detected automatically.\n');
else
    disp(feature_vars(:));

    T_feature_summary = summarize_numeric_variables(T_mc, feature_vars);

    fprintf('\nFeature summary.\n');
    disp(T_feature_summary);
end

%% ========================================================================
% 11. Duplicate row check
% ========================================================================

key_vars = string.empty(1, 0);

for v = [condition_var, step_var, realization_var, patch_var]
    if v ~= ""
        key_vars(end + 1) = v; %#ok<SAGROW>
    end
end

if numel(key_vars) >= 3

    T_keys = T_mc(:, cellstr(key_vars));
    [~, unique_idx] = unique(T_keys, 'rows', 'stable');

    n_duplicates = height(T_mc) - numel(unique_idx);

    fprintf('\nDuplicate key check.\n');
    fprintf('Key variables: %s\n', strjoin(key_vars, ', '));
    fprintf('Duplicate rows based on key variables: %d\n', n_duplicates);

    if n_duplicates == 0
        fprintf('PASS. No duplicate rows detected using available key variables.\n');
    else
        warning('Duplicate rows detected using available key variables.');
    end

else

    fprintf('\nDuplicate key check skipped. Not enough key variables detected.\n');

end

%% ========================================================================
% 12. Within-step q variability
% ========================================================================

if condition_var ~= "" && step_var ~= ""

    T_q_within_step = summarize_within_group_numeric( ...
        T_mc, ...
        [condition_var, step_var], ...
        q_var);

    fprintf('\nWithin-step q variability.\n');
    disp(head(T_q_within_step, min(20, height(T_q_within_step))));

    fprintf('\nSummary of within-step q std.\n');
    fprintf('min std(q)    = %.6g\n', min(T_q_within_step.std_value, [], 'omitnan'));
    fprintf('median std(q) = %.6g\n', median(T_q_within_step.std_value, 'omitnan'));
    fprintf('max std(q)    = %.6g\n', max(T_q_within_step.std_value, [], 'omitnan'));

else

    T_q_within_step = table();

end

%% ========================================================================
% 13. Save sanity outputs
% ========================================================================

analysis_dir = fullfile(run_dir, 'analysis', 'level_01_sanity_check');

if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

save(fullfile(analysis_dir, 'level_01_sanity_check.mat'), ...
    'T_mc', ...
    'MC', ...
    'T_condition_counts', ...
    'T_step_counts', ...
    'T_missing', ...
    'feature_vars', ...
    'T_q_within_step', ...
    '-v7.3');

writetable(T_missing, fullfile(analysis_dir, 'missing_values.csv'));

if exist('T_condition_counts', 'var') && ~isempty(T_condition_counts)
    writetable(T_condition_counts, fullfile(analysis_dir, 'condition_counts.csv'));
end

if exist('T_step_counts', 'var') && ~isempty(T_step_counts)
    writetable(T_step_counts, fullfile(analysis_dir, 'step_counts.csv'));
end

if exist('T_feature_summary', 'var') && ~isempty(T_feature_summary)
    writetable(T_feature_summary, fullfile(analysis_dir, 'feature_summary.csv'));
end

if exist('T_q_within_step', 'var') && ~isempty(T_q_within_step)
    writetable(T_q_within_step, fullfile(analysis_dir, 'within_step_q_variability.csv'));
end

fprintf('\nSaved Level 1 sanity check outputs to:\n%s\n', analysis_dir);

fprintf('\nLevel 1 sanity check completed.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function run_dir = find_latest_run_dir(output_root)

D = dir(output_root);
D = D([D.isdir]);

names = string({D.name});
mask = names ~= "." & names ~= "..";

D = D(mask);

if isempty(D)
    error('No run folders found in: %s', output_root);
end

[~, idx] = max([D.datenum]);

run_dir = fullfile(output_root, D(idx).name);

end

function mat_file = find_mc_result_file(data_dir)

if ~exist(data_dir, 'dir')
    error('Data folder not found: %s', data_dir);
end

D = dir(fullfile(data_dir, '*mc*sweep*result*.mat'));

if isempty(D)
    D = dir(fullfile(data_dir, '*.mat'));
end

if isempty(D)
    error('No MAT result file found in: %s', data_dir);
end

[~, idx] = max([D.datenum]);

mat_file = fullfile(data_dir, D(idx).name);

end

function var_name = find_first_existing_var(vars, candidates)

vars = string(vars);
candidates = string(candidates);

var_name = "";

for i = 1:numel(candidates)
    idx = find(strcmpi(vars, candidates(i)), 1);

    if ~isempty(idx)
        var_name = vars(idx);
        return;
    end
end

end

function T_count = count_by_group(T, group_vars)

group_vars = cellstr(string(group_vars));

[G, T_count] = findgroups(T(:, group_vars));

T_count.GroupCount = splitapply( ...
    @numel, ...
    ones(height(T), 1), ...
    G);

end

function T_row = varargin_to_table_row(varargin)

n = nargin;

T_row = table();

for i = 1:n
    value_i = varargin{i};

    if iscell(value_i)
        value_i = value_i{1};
    else
        value_i = value_i(1);
    end

    T_row.(sprintf('Var%d', i)) = value_i;
end

end

function T_missing = summarize_missing_values(T)

vars = string(T.Properties.VariableNames);

n_vars = numel(vars);

var_name = strings(n_vars, 1);
class_name = strings(n_vars, 1);
n_missing = zeros(n_vars, 1);
n_nonfinite = zeros(n_vars, 1);

for i = 1:n_vars

    v = vars(i);
    x = T.(v);

    var_name(i) = v;
    class_name(i) = string(class(x));

    try
        n_missing(i) = sum(ismissing(x));
    catch
        n_missing(i) = NaN;
    end

    if isnumeric(x)
        n_nonfinite(i) = sum(~isfinite(x));
    else
        n_nonfinite(i) = 0;
    end

end

T_missing = table(var_name, class_name, n_missing, n_nonfinite);

end

function feature_vars = detect_feature_variables(T)

vars = string(T.Properties.VariableNames);
vars_lower = lower(vars);

is_numeric = false(size(vars));

for i = 1:numel(vars)
    is_numeric(i) = isnumeric(T.(vars(i)));
end

feature_keywords = [ ...
    "entropy", ...
    "anisotropy", ...
    "width", ...
    "spread", ...
    "sharp", ...
    "peak", ...
    "radial", ...
    "angular", ...
    "coherence", ...
    "energy", ...
    "spectrum", ...
    "spectral"];

metadata_keywords = [ ...
    "condition", ...
    "step", ...
    "patch", ...
    "realization", ...
    "rep", ...
    "seed", ...
    "sim_", ...
    "req_", ...
    "omega", ...
    "q_", ...
    "qtheory", ...
    "q"];

is_feature_name = false(size(vars));

for k = 1:numel(feature_keywords)
    is_feature_name = is_feature_name | contains(vars_lower, feature_keywords(k));
end

is_metadata = false(size(vars));

for k = 1:numel(metadata_keywords)
    is_metadata = is_metadata | startsWith(vars_lower, metadata_keywords(k));
end

feature_vars = vars(is_numeric & is_feature_name & ~is_metadata);

end

function T_summary = summarize_numeric_variables(T, selected_vars)

selected_vars = string(selected_vars);
n = numel(selected_vars);

variable = strings(n, 1);
n_missing = zeros(n, 1);
min_value = zeros(n, 1);
max_value = zeros(n, 1);
mean_value = zeros(n, 1);
median_value = zeros(n, 1);
std_value = zeros(n, 1);

for i = 1:n

    v = selected_vars(i);
    x = T.(v);

    variable(i) = v;
    n_missing(i) = sum(~isfinite(x));

    min_value(i) = min(x, [], 'omitnan');
    max_value(i) = max(x, [], 'omitnan');
    mean_value(i) = mean(x, 'omitnan');
    median_value(i) = median(x, 'omitnan');
    std_value(i) = std(x, 'omitnan');

end

T_summary = table( ...
    variable, ...
    n_missing, ...
    min_value, ...
    max_value, ...
    mean_value, ...
    median_value, ...
    std_value);

end

function T_group = summarize_within_group_numeric(T, group_vars, value_var)

group_vars = cellstr(string(group_vars));
value_var = char(value_var);

[G, T_group] = findgroups(T(:, group_vars));

x = T.(value_var);

T_group.n_value = splitapply(@numel, x, G);

T_group.mean_value = splitapply( ...
    @(y) mean(y, 'omitnan'), ...
    x, ...
    G);

T_group.median_value = splitapply( ...
    @(y) median(y, 'omitnan'), ...
    x, ...
    G);

T_group.std_value = splitapply( ...
    @(y) std(y, 'omitnan'), ...
    x, ...
    G);

T_group.min_value = splitapply( ...
    @(y) min(y, [], 'omitnan'), ...
    x, ...
    G);

T_group.max_value = splitapply( ...
    @(y) max(y, [], 'omitnan'), ...
    x, ...
    G);

end

