%% run_test_14_dx_dz_resolution_sweep.m
% Generate the Test 14 clean dx=dz spatial-resolution sensitivity dataset.
%
% This script intentionally does not train any model. It generates local and
% global REQ features/mappings over a short dx=dz sweep. The analysis script
% applies the already trained Test 12 deployment model.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
experiments_dir = fileparts(this_file);
root_dir = fileparts(experiments_dir);

addpath(root_dir);
root_dir = setup_adaptive_req();

CFG = adaptive_req.config.load_profile_config( ...
    'test_14_dx_dz_resolution_sweep', ...
    'RootDir', root_dir);

if ~isfield(CFG.EXP, 'timestamp') || isempty(CFG.EXP.timestamp)
    CFG.EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
end

%% Output folders

SAVE = struct();
SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
SAVE.run_name = sprintf('%s_%s', CFG.EXP.name, CFG.EXP.timestamp);
SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

make_dir_if_needed(SAVE.output_dir);
make_dir_if_needed(SAVE.table_dir);
make_dir_if_needed(SAVE.data_dir);

%% Paired dx=dz loop

dx_values = CFG.RESOLUTION.dx_dz(:).';
T_all = table();
MC_PARTS = struct([]);
T_run_status = table();
internal_root_dir = tempname(tempdir);
make_dir_if_needed(internal_root_dir);
cleanup_internal = onCleanup(@() cleanup_temp_root(internal_root_dir));

fprintf('\nRunning Test 14 dx=dz resolution sweep.\n');
fprintf('dx=dz values: %s m\n', mat2str(dx_values));

for di = 1:numel(dx_values)
    dx_i = dx_values(di);

    CFG_i = CFG;
    CFG_i.SIM.dx = dx_i;
    CFG_i.SIM.dz = dx_i;
    CFG_i.EXP.name = sprintf('%s_dx_%g', CFG.EXP.name, dx_i);
    CFG_i.EXP.profile_name = CFG.EXP.name;
    CFG_i.EXP.timestamp = CFG.EXP.timestamp;
    CFG_i.EXP.seed_base = CFG.EXP.seed_base + 1000000 * (di - 1);

    fprintf('\n=== dx = dz = %.4g m (%d / %d) ===\n', ...
        dx_i, di, numel(dx_values));

    [T_i, MC_i] = adaptive_req.studies.run_mc_sweep( ...
        CFG_i, ...
        'RootDir', internal_root_dir, ...
        'MaxConditions', Inf, ...
        'SaveResults', false, ...
        'SaveIntermediate', false, ...
        'ContinueOnError', false, ...
        'Verbose', true);

    assert(all(MC_i.T_status.success), ...
        'Test 14 failed for dx=dz=%g.', dx_i);

    T_i.resolution_idx = di * ones(height(T_i), 1);
    T_i.dx_dz_value = dx_i * ones(height(T_i), 1);
    T_i.pixels_per_wavelength = (T_i.SIM_cs_bg ./ T_i.SIM_f0) ./ T_i.SIM_dx;
    T_i.pixels_per_window = T_i.REQ_win_size;
    T_i.k0_over_knyquist = T_i.k0_true ./ (pi ./ T_i.SIM_dx);
    T_i.aperture_label = aperture_label_from_step(T_i.step_idx, CFG.EXP.num_steps);

    T_all = concat_tables(T_all, T_i);

    MC_PARTS(di).dx_dz = dx_i;
    MC_PARTS(di).MC = MC_i;

    T_status_i = MC_i.T_status;
    T_status_i.resolution_idx = di * ones(height(T_status_i), 1);
    T_status_i.dx_dz_value = dx_i * ones(height(T_status_i), 1);
    T_run_status = concat_tables(T_run_status, T_status_i);
end

%% Validate and save

assert(~isempty(T_all), 'Test 14 generated an empty table.');
assert(~ismember('req_curve', T_all.Properties.VariableNames), ...
    'The generated dataset unexpectedly contains the heavy req_curve.');

required_vars = { ...
    'req_mapping', ...
    'global_req_mapping', ...
    'q_theory', ...
    'q_global_theory', ...
    'q_local_minus_global', ...
    'REQ_M', ...
    'REQ_cs_guess', ...
    'SIM_f0', ...
    'SIM_cs_bg', ...
    'SIM_dx', ...
    'SIM_dz', ...
    'pixels_per_wavelength', ...
    'pixels_per_window', ...
    'k0_over_knyquist'};

missing_vars = setdiff(required_vars, T_all.Properties.VariableNames);
assert(isempty(missing_vars), ...
    'Missing Test 14 required variables: %s', strjoin(missing_vars, ', '));

T_csv = remove_heavy_table_columns(T_all);
writetable(T_csv, fullfile(SAVE.table_dir, ...
    'test14_dx_dz_resolution_sweep_table.csv'));
writetable(T_run_status, fullfile(SAVE.table_dir, ...
    'test14_dx_dz_resolution_sweep_status.csv'));

MC = struct();
MC.CFG = CFG;
MC.SAVE = SAVE;
MC.MC_PARTS = MC_PARTS;
MC.T_status = T_run_status;
MC.n_rows = height(T_all);
MC.n_dx = numel(dx_values);
MC.dx_dz_values = dx_values;

save(fullfile(SAVE.data_dir, 'test14_dx_dz_resolution_sweep_results.mat'), ...
    'T_all', 'MC', 'CFG', '-v7.3');

fprintf('\nTest 14 dataset completed successfully.\n');
fprintf('Rows generated: %d\n', height(T_all));
fprintf('Unique dx     : %s\n', mat2str(unique(T_all.SIM_dx).'));
fprintf('Unique cs_bg  : %s\n', mat2str(unique(T_all.SIM_cs_bg).'));
fprintf('Unique REQ_M  : %s\n', mat2str(unique(T_all.REQ_M).'));
fprintf('Output folder:\n%s\n', SAVE.output_dir);

%% Local functions

function labels = aperture_label_from_step(step_idx, num_steps)

step_idx = double(step_idx(:));
labels = strings(numel(step_idx), 1);

if num_steps == 1
    labels(:) = "single";
    return;
end

for i = 1:numel(step_idx)
    t = (step_idx(i) - 1) / max(num_steps - 1, 1);
    if t <= 1/3
        labels(i) = "narrow";
    elseif t <= 2/3
        labels(i) = "mid";
    else
        labels(i) = "wide";
    end
end

end

function T = remove_heavy_table_columns(T)

vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars)
    v = char(vars(i));
    drop(i) = iscell(T.(v)) || isstruct(T.(v));
end
T(:, drop) = [];

end

function T = concat_tables(A, B)

if isempty(A)
    T = B;
    return;
end
if isempty(B)
    T = A;
    return;
end

vars_all = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');
A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);
T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];

end

function T = add_missing_columns(T, vars_all)

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    if ismember(vars_all(i), vars)
        continue;
    end
    name_i = char(vars_all(i));
    string_like = any(endsWith(vars_all(i), ...
        ["name", "label", "type", "role", "set", "model"]));
    if string_like
        T.(name_i) = strings(height(T), 1);
    else
        T.(name_i) = nan(height(T), 1);
    end
end

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end

function cleanup_temp_root(path_i)

if exist(path_i, 'dir')
    try
        rmdir(path_i, 's');
    catch ME
        warning('Could not remove temporary Test 14 working directory:\n%s\n%s', ...
            path_i, ME.message);
    end
end

end
