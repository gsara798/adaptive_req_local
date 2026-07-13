function [T_mc, MC] = run_mc_sweep(CFG_or_profile, varargin)
%RUN_MC_SWEEP Run a Monte Carlo parameter sweep for adaptive_req.
%
% This function loops over the parameter conditions defined by CFG.SWEEP.
% For each condition, it resolves a scalar RUN configuration and calls
% adaptive_req.studies.run_aperture_sweep.
%
% Usage:
%
%   [T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
%       'test_04_mc_baseline', ...
%       'RootDir', root_dir);
%
%   [T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
%       CFG, ...
%       'RootDir', root_dir);
%
% For testing:
%
%   [T_mc, MC] = adaptive_req.studies.run_mc_sweep( ...
%       'test_04_mc_baseline', ...
%       'RootDir', root_dir, ...
%       'MaxConditions', 2);

p = inputParser;
p.FunctionName = 'adaptive_req.studies.run_mc_sweep';

addRequired(p, 'CFG_or_profile', @(x) isstruct(x) || ischar(x) || isstring(x));

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'ConfigDir', '', @(x) ischar(x) || isstring(x));

addParameter(p, 'ConditionIndices', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'MaxConditions', Inf, @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'SaveResults', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveIntermediate', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ContinueOnError', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, CFG_or_profile, varargin{:});

root_dir = char(p.Results.RootDir);
config_dir = char(p.Results.ConfigDir);

save_results = logical(p.Results.SaveResults);
save_intermediate = logical(p.Results.SaveIntermediate);
continue_on_error = logical(p.Results.ContinueOnError);
verbose = logical(p.Results.Verbose);

%% ------------------------------------------------------------------------
% Load configuration
% -------------------------------------------------------------------------

if isstruct(CFG_or_profile)

    CFG = CFG_or_profile;

    if ~isfield(CFG.EXP, 'timestamp') || isempty(CFG.EXP.timestamp)
        CFG.EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
    end

    if ~isfield(CFG.EXP, 'profile_name')
        CFG.EXP.profile_name = 'struct_input';
    end

else

    profile_name = char(CFG_or_profile);

    if isempty(config_dir)
        CFG = adaptive_req.config.load_profile_config( ...
            profile_name, ...
            'RootDir', root_dir);
    else
        CFG = adaptive_req.config.load_profile_config( ...
            profile_name, ...
            'RootDir', root_dir, ...
            'ConfigDir', config_dir);
    end
end

OUTPUT = CFG.OUTPUT;

save_condition_table = getfield_with_default(OUTPUT, ...
    'save_condition_table', true);

save_condition_mat = getfield_with_default(OUTPUT, ...
    'save_condition_mat', false);

save_condition_summary_figures = getfield_with_default(OUTPUT, ...
    'save_condition_summary_figures', false);
%% ------------------------------------------------------------------------
% Build sweep matrix
% -------------------------------------------------------------------------

[T_axes, T_conditions, CFG_base] = ...
    adaptive_req.config.build_sweep_preview(CFG);

n_conditions_total = height(T_conditions);

if isempty(p.Results.ConditionIndices)
    condition_indices = 1:n_conditions_total;
else
    condition_indices = p.Results.ConditionIndices(:).';
end

if any(condition_indices < 1) || any(condition_indices > n_conditions_total)
    error('ConditionIndices contains values outside the valid range.');
end

if isfinite(p.Results.MaxConditions)
    n_keep = min(numel(condition_indices), round(p.Results.MaxConditions));
    condition_indices = condition_indices(1:n_keep);
end

n_to_run = numel(condition_indices);

%% ------------------------------------------------------------------------
% Build output folders
% -------------------------------------------------------------------------

SAVE = CFG.SAVE;

SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
SAVE.run_name = sprintf('%s_%s', CFG.EXP.name, CFG.EXP.timestamp);

SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.condition_dir = fullfile(SAVE.output_dir, 'conditions');
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');

if save_results || save_intermediate
    make_dir_if_needed(SAVE.output_dir);
    make_dir_if_needed(SAVE.condition_dir);
    make_dir_if_needed(SAVE.table_dir);
    make_dir_if_needed(SAVE.data_dir);
end

%% ------------------------------------------------------------------------
% Print sweep summary
% -------------------------------------------------------------------------

rows_per_condition = ...
    CFG.EXP.num_steps * ...
    CFG.EXP.num_realizations * ...
    CFG.EXP.num_patches;

expected_rows_selected = n_to_run * rows_per_condition;

if verbose

    fprintf('\nRunning MC parameter sweep.\n');
    fprintf('Profile            : %s\n', CFG.EXP.profile_name);
    fprintf('Experiment name    : %s\n', CFG.EXP.name);
    fprintf('Total conditions   : %d\n', n_conditions_total);
    fprintf('Selected conditions: %d\n', n_to_run);
    fprintf('Rows per condition : %d\n', rows_per_condition);
    fprintf('Expected rows      : %d\n', expected_rows_selected);

    fprintf('\nActive sweep axes.\n');
    disp(T_axes);
end

%% ------------------------------------------------------------------------
% Main MC loop
% -------------------------------------------------------------------------

T_list = cell(n_to_run, 1);
sweep_list = cell(n_to_run, 1);
run_list = cell(n_to_run, 1);

status_rows(n_to_run) = struct();

t_all = tic;

for pos = 1:n_to_run

    condition_id = condition_indices(pos);

    CFG_i = apply_condition_to_cfg(CFG, T_axes, T_conditions, condition_id);

    CFG_i.EXP.condition_id = condition_id;
    CFG_i.EXP.condition_position = pos;
    CFG_i.EXP.n_conditions_total = n_conditions_total;
    CFG_i.EXP.n_conditions_selected = n_to_run;

    CFG_i.EXP.seed_base = CFG.EXP.seed_base + 100000 * (condition_id - 1);

    condition_label = make_condition_label(T_axes, T_conditions, condition_id);

    condition_output_dir = fullfile( ...
        SAVE.condition_dir, ...
        sprintf('condition_%03d_%s', condition_id, condition_label));

    if verbose
        fprintf('\nCondition %d / %d\n', pos, n_to_run);
        fprintf('  condition_id = %d\n', condition_id);
        fprintf('  label        = %s\n', condition_label);
    end

    t_condition = tic;

    try

        RUN_i = adaptive_req.config.resolve_run_config( ...
            CFG_i, ...
            'RootDir', root_dir, ...
            'OutputDir', condition_output_dir, ...
            'MakeDirs', true);

        [T_i, sweep_i] = run_single_condition(RUN_i);

        T_i.condition_id = repmat(condition_id, height(T_i), 1);
        T_i.condition_position = repmat(pos, height(T_i), 1);
        T_i.condition_label = repmat(string(condition_label), height(T_i), 1);

        for ax_idx = 1:height(T_axes)
            var_name = char(T_axes.variable_name(ax_idx));
            T_i.(var_name) = repmat( ...
                T_conditions.(var_name)(condition_id), ...
                height(T_i), 1);
        end

        T_list{pos} = T_i;
        sweep_list{pos} = sweep_i;
        run_list{pos} = RUN_i;

        % -------------------------------------------------------------
        % Save per-condition outputs
        % -------------------------------------------------------------
        
        if save_condition_table
        
            T_i_csv = remove_heavy_table_columns(T_i);
        
            writetable(T_i_csv, fullfile(RUN_i.SAVE.table_dir, ...
                sprintf('condition_%03d_table.csv', condition_id)));
        
        end
        
        if save_condition_mat
        
            save(fullfile(RUN_i.SAVE.data_dir, ...
                sprintf('condition_%03d_result.mat', condition_id)), ...
                'T_i', ...
                'sweep_i', ...
                'RUN_i', ...
                '-v7.3');
        
        end
        
        if save_condition_summary_figures
        
            RUN_i.PLOT.save_summary_figures = true;
        
            adaptive_req.figures.plot_run_summary( ...
                T_i, sweep_i, RUN_i);
        
        end

        elapsed_i = toc(t_condition);

        status_rows(pos).condition_id = condition_id;
        status_rows(pos).condition_position = pos;
        status_rows(pos).condition_label = string(condition_label);
        status_rows(pos).success = true;
        status_rows(pos).n_rows = height(T_i);
        status_rows(pos).elapsed_s = elapsed_i;
        status_rows(pos).error_message = "";

        if verbose
            fprintf('  completed in %.2f s\n', elapsed_i);
            fprintf('  rows = %d\n', height(T_i));
        end

        if save_intermediate
            save(fullfile(RUN_i.SAVE.data_dir, ...
                sprintf('condition_%03d_intermediate.mat', condition_id)), ...
                'T_i', ...
                'sweep_i', ...
                'RUN_i', ...
                '-v7.3');
        end

    catch ME

        elapsed_i = toc(t_condition);

        status_rows(pos).condition_id = condition_id;
        status_rows(pos).condition_position = pos;
        status_rows(pos).condition_label = string(condition_label);
        status_rows(pos).success = false;
        status_rows(pos).n_rows = 0;
        status_rows(pos).elapsed_s = elapsed_i;
        status_rows(pos).error_message = string(ME.message);

        warning('Condition %d failed: %s', condition_id, ME.message);

        if ~continue_on_error
            rethrow(ME);
        end
    end
end

%% ------------------------------------------------------------------------
% Combine outputs
% -------------------------------------------------------------------------

valid_tables = ~cellfun(@isempty, T_list);

if any(valid_tables)
    T_mc = vertcat(T_list{valid_tables});
else
    T_mc = table();
end

T_status = struct2table(status_rows);

MC = struct();

MC.CFG = CFG;
MC.CFG_base = CFG_base;

MC.T_axes = T_axes;
MC.T_conditions = T_conditions;
MC.T_status = T_status;

MC.SAVE = SAVE;

MC.condition_indices = condition_indices;
MC.n_conditions_total = n_conditions_total;
MC.n_conditions_selected = n_to_run;

MC.rows_per_condition = rows_per_condition;
MC.expected_rows_selected = expected_rows_selected;
MC.actual_rows = height(T_mc);

MC.total_elapsed_s = toc(t_all);

MC.sweep_list = sweep_list;
MC.run_list = run_list;

if verbose
    fprintf('\nMC sweep completed.\n');
    fprintf('Successful conditions = %d / %d\n', sum(T_status.success), height(T_status));
    fprintf('Generated rows = %d\n', height(T_mc));
    fprintf('Total elapsed time = %.2f s\n', MC.total_elapsed_s);
end

%% ------------------------------------------------------------------------
% Save outputs
% -------------------------------------------------------------------------

if save_results

    save(fullfile(SAVE.data_dir, 'mc_sweep_results.mat'), ...
        'T_mc', ...
        'MC', ...
        '-v7.3');

    T_csv = remove_heavy_table_columns(T_mc);

    writetable(T_csv, fullfile(SAVE.table_dir, ...
        'mc_sweep_table.csv'));

    writetable(T_axes, fullfile(SAVE.table_dir, ...
        'sweep_axes.csv'));

    writetable(T_conditions, fullfile(SAVE.table_dir, ...
        'sweep_conditions.csv'));

    writetable(T_status, fullfile(SAVE.table_dir, ...
        'sweep_status.csv'));

    if verbose
        fprintf('\nSaved MC MAT output to:\n%s\n', ...
            fullfile(SAVE.data_dir, 'mc_sweep_results.mat'));

        fprintf('\nSaved MC CSV table to:\n%s\n', ...
            fullfile(SAVE.table_dir, 'mc_sweep_table.csv'));
    end
end

end

%% =========================================================================
% Local helper functions
% =========================================================================

function [T_i, sweep_i] = run_single_condition(RUN_i)

[T_i, sweep_i] = adaptive_req.studies.run_aperture_sweep( ...
    RUN_i.cfg, RUN_i.feat_cfg, ...
    'SamplingMode', RUN_i.EXP.sampling_mode, ...
    'NumSteps', RUN_i.EXP.num_steps, ...
    'StepIndices', RUN_i.EXP.step_indices, ...
    'NumRealizations', RUN_i.EXP.num_realizations, ...
    'NumPatches', RUN_i.EXP.num_patches, ...
    'SeedBase', RUN_i.EXP.seed_base, ...
    'ReqOptions', RUN_i.req_options, ...
    'StoreWavefields', RUN_i.PLOT.store_wavefields, ...
    'StoreReqCurve', RUN_i.OUTPUT.store_req_curve, ...
    'StoreReqMapping', RUN_i.OUTPUT.store_req_mapping, ...
    'ComputeGlobalReq', RUN_i.OUTPUT.compute_global_req, ...
    'StoreGlobalReqMapping', RUN_i.OUTPUT.store_global_req_mapping, ...
    'StoreReqMetadata', RUN_i.OUTPUT.store_req_metadata, ...
    'StoreFeatureStruct', RUN_i.OUTPUT.store_feature_struct, ...
    'PlotStepDiagnostics', RUN_i.PLOT.show_step_diagnostics || RUN_i.PLOT.save_step_diagnostics, ...
    'SaveStepDiagnostics', RUN_i.PLOT.save_step_diagnostics, ...
    'StepDiagnosticPatchIndex', RUN_i.EXP.selected_patch, ...
    'StepDiagnosticDir', RUN_i.SAVE.step_diag_dir, ...
    'StepDiagnosticVisible', RUN_i.PLOT.step_diagnostic_visible, ...
    'CloseStepDiagnosticsAfterSave', RUN_i.PLOT.close_step_diagnostics_after_save, ...
    'SaveDiagnosticPNG', RUN_i.SAVE.save_png, ...
    'SaveDiagnosticPDF', RUN_i.SAVE.save_pdf, ...
    'SaveDiagnosticFIG', RUN_i.SAVE.save_fig, ...
    'DiagnosticResolution', RUN_i.SAVE.png_resolution, ...
    'Verbose', RUN_i.OUTPUT.verbose);

end

function CFG_i = apply_condition_to_cfg(CFG, T_axes, T_conditions, condition_id)

CFG_i = CFG;

for ax_idx = 1:height(T_axes)

    path_i = T_axes.path(ax_idx);
    var_name_i = char(T_axes.variable_name(ax_idx));

    value_i = T_conditions.(var_name_i)(condition_id);

    CFG_i = set_by_path(CFG_i, path_i, value_i);

end

end

function S = set_by_path(S, path, value)

parts = cellstr(split(string(path), '.'));
S = set_by_parts(S, parts, value);

end

function S = set_by_parts(S, parts, value)

field_i = parts{1};

if numel(parts) == 1
    S.(field_i) = value;
else
    S.(field_i) = set_by_parts(S.(field_i), parts(2:end), value);
end

end

function label = make_condition_label(T_axes, T_conditions, condition_id)

parts = strings(1, height(T_axes));

for ax_idx = 1:height(T_axes)

    var_name_i = char(T_axes.variable_name(ax_idx));
    value_i = T_conditions.(var_name_i)(condition_id);

    parts(ax_idx) = string(var_name_i) + "_" + value_to_label(value_i);

end

label = strjoin(parts, "__");
label = regexprep(label, '[^A-Za-z0-9_]+', '_');

end

function label = value_to_label(value)

if islogical(value)
    label = string(value);
    return;
end

if isnumeric(value)

    if isinf(value)
        label = "Inf";
    elseif isnan(value)
        label = "NaN";
    else
        label = string(sprintf('%.6g', value));
    end

    label = strrep(label, '-', 'm');
    label = strrep(label, '.', 'p');

else

    label = string(value);

end

end

function T = remove_heavy_table_columns(T)

vars_to_remove = { ...
    'req_curve', ...
    'req_mapping', ...
    'global_req_mapping', ...
    'feat', ...
    'feature_struct'};

vars_to_remove = vars_to_remove(ismember(vars_to_remove, ...
    T.Properties.VariableNames));

if ~isempty(vars_to_remove)
    T(:, vars_to_remove) = [];
end

end

function make_dir_if_needed(folder_path)

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

end

function val = getfield_with_default(S, field_name, default_val)

if isstruct(S) && isfield(S, field_name)
    val = S.(field_name);
else
    val = default_val;
end

end
