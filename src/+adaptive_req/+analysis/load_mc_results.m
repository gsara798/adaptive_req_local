function [T_mc, MC, PATHS] = load_mc_results(experiment_name, varargin)
%LOAD_MC_RESULTS Load Monte Carlo results from the outputs folder.
%
% Usage
%   [T_mc, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
%       'test_06_feature_q_baseline', ...
%       'RootDir', root_dir);
%
% Optional inputs
%   RootDir:
%       Project root directory.
%
%   RunDir:
%       Specific run folder. If empty, the latest run is used.
%
%   MatFile:
%       Specific MAT file to load. If empty, the latest MC result MAT file
%       inside the selected run folder is used.
%
% Outputs
%   T_mc:
%       Main Monte Carlo table.
%
%   MC:
%       Metadata structure saved by run_mc_sweep.
%
%   PATHS:
%       Structure with useful paths.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.load_mc_results';

addRequired(p, 'experiment_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'RunDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'MatFile', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, experiment_name, varargin{:});

experiment_name = char(p.Results.experiment_name);
root_dir = char(p.Results.RootDir);
run_dir = char(p.Results.RunDir);
mat_file = char(p.Results.MatFile);
verbose = logical(p.Results.Verbose);

output_root = fullfile(root_dir, 'outputs', experiment_name);

if isempty(run_dir)
    run_dir = find_latest_run_dir(output_root);
end

if isempty(mat_file)
    data_dir = fullfile(run_dir, 'data');
    mat_file = find_latest_result_file(data_dir);
end

if verbose
    fprintf('\nLoading MC results.\n');
    fprintf('Experiment: %s\n', experiment_name);
    fprintf('Run folder: %s\n', run_dir);
    fprintf('MAT file  : %s\n', mat_file);
end

S = load(mat_file);

if ~isfield(S, 'T_mc')
    error('Selected file does not contain T_mc: %s', mat_file);
end

T_mc = S.T_mc;

if isfield(S, 'MC')
    MC = S.MC;
else
    MC = struct();
end

PATHS = struct();

PATHS.root_dir = root_dir;
PATHS.experiment_name = experiment_name;
PATHS.output_root = output_root;
PATHS.run_dir = run_dir;
PATHS.data_dir = fullfile(run_dir, 'data');
PATHS.table_dir = fullfile(run_dir, 'tables');
PATHS.figure_dir = fullfile(run_dir, 'figures');
PATHS.analysis_dir = fullfile(run_dir, 'analysis');
PATHS.mat_file = mat_file;

if ~exist(PATHS.analysis_dir, 'dir')
    mkdir(PATHS.analysis_dir);
end

if verbose
    fprintf('\nLoaded T_mc.\n');
    fprintf('Rows    : %d\n', height(T_mc));
    fprintf('Columns : %d\n', width(T_mc));
end

end

function run_dir = find_latest_run_dir(output_root)

if ~exist(output_root, 'dir')
    error('Output root not found: %s', output_root);
end

D = dir(output_root);
D = D([D.isdir]);

names = string({D.name});
mask = names ~= "." & names ~= "..";

D = D(mask);

if isempty(D)
    error('No run folders found in: %s', output_root);
end

has_results = false(size(D));

for i = 1:numel(D)
    data_dir_i = fullfile(output_root, D(i).name, 'data');
    if ~exist(data_dir_i, 'dir')
        continue;
    end

    mat_files = dir(fullfile(data_dir_i, '*mc*sweep*result*.mat'));
    if isempty(mat_files)
        mat_files = dir(fullfile(data_dir_i, '*.mat'));
    end

    has_results(i) = ~isempty(mat_files);
end

D_valid = D(has_results);

if isempty(D_valid)
    error('No completed run folders with data MAT files found in: %s', ...
        output_root);
end

[~, idx] = max([D_valid.datenum]);

run_dir = fullfile(output_root, D_valid(idx).name);

end

function mat_file = find_latest_result_file(data_dir)

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
