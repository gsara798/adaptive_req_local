function [T_assoc, T_binned] = compute_groupwise_feature_q_associations(T, feature_vars, varargin)
%COMPUTE_GROUPWISE_FEATURE_Q_ASSOCIATIONS Quantify feature-q associations.
%
% Purpose
%   Quantifies how strongly q is associated with one or more local features,
%   optionally within controlled groups.
%
%   This function does not assume that the relation is linear. It computes:
%
%       Pearson correlation    linear association
%       Spearman correlation   monotonic association
%       Kendall correlation    rank-based ordinal association
%       Linear fit             diagnostic only
%       Binned trend           non-parametric trend
%
% Usage
%   [T_assoc, T_binned] = adaptive_req.analysis.compute_groupwise_feature_q_associations( ...
%       T_step, ...
%       ["ang_entropy", "radial_entropy"], ...
%       'QVar', 'q_mean', ...
%       'GroupVars', ["SIM_WaveModel", "REQ_M", "SIM_f0", "SIM_cs_bg"]);
%
% Inputs
%   T:
%       Table containing q and feature columns.
%
%   feature_vars:
%       Feature names. Each name can be either the raw feature name
%       such as "ang_entropy", or the summary name such as
%       "ang_entropy_mean".
%
% Name-value options
%   QVar:
%       Target q variable. Default: 'q_mean'.
%
%   GroupVars:
%       Variables used to define controlled groups.
%       Use [] for global association.
%
%   NumBins:
%       Number of bins for non-parametric binned trend.
%
%   BinMode:
%       'equal_count' or 'equal_width'.
%
%   MinN:
%       Minimum number of samples required to compute correlations.
%
% Outputs
%   T_assoc:
%       One row per group and feature.
%
%   T_binned:
%       Binned q-feature trend for each group and feature.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.compute_groupwise_feature_q_associations';

addRequired(p, 'T', @istable);

addRequired(p, 'feature_vars', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'QVar', 'q_mean', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'GroupVars', string.empty(1, 0), ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'NumBins', 6, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2);

addParameter(p, 'BinMode', 'equal_count', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'MinN', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2);

addParameter(p, 'Verbose', true, ...
    @(x) islogical(x) || isnumeric(x));

parse(p, T, feature_vars, varargin{:});

feature_vars = string(feature_vars);
q_var = char(p.Results.QVar);
group_vars = string(p.Results.GroupVars);
num_bins = round(p.Results.NumBins);
bin_mode = lower(string(p.Results.BinMode));
min_n = round(p.Results.MinN);
verbose = logical(p.Results.Verbose);

available_vars = string(T.Properties.VariableNames);

if ~ismember(q_var, available_vars)
    error('QVar not found in table: %s', q_var);
end

group_vars = group_vars(group_vars ~= "");

for i = 1:numel(group_vars)
    if ~ismember(group_vars(i), available_vars)
        error('Group variable not found in table: %s', group_vars(i));
    end
end

resolved_features = strings(numel(feature_vars), 1);

for i = 1:numel(feature_vars)
    resolved_features(i) = resolve_feature_variable(T, feature_vars(i));
end

if verbose
    fprintf('\nFeature-q association analysis.\n');
    fprintf('Q variable: %s\n', q_var);

    if isempty(group_vars)
        fprintf('Grouping: global\n');
    else
        fprintf('Grouping: %s\n', strjoin(group_vars, ', '));
    end

    fprintf('\nResolved features:\n');
    disp(table(feature_vars(:), resolved_features(:), ...
        'VariableNames', {'InputFeature', 'ResolvedVariable'}));
end

if isempty(group_vars)

    G = ones(height(T), 1);
    T_groups = table();
    n_groups = 1;

else

    [G, T_groups] = findgroups(T(:, cellstr(group_vars)));
    n_groups = height(T_groups);

end

T_assoc = table();
T_binned = table();

for gidx = 1:n_groups

    group_mask = G == gidx;

    if isempty(group_vars)
        group_row = table();
    else
        group_row = T_groups(gidx, :);
    end

    for fidx = 1:numel(feature_vars)

        feature_input = feature_vars(fidx);
        feature_var = resolved_features(fidx);

        x_all = T.(char(feature_var))(group_mask);
        y_all = T.(q_var)(group_mask);

        valid = isfinite(x_all) & isfinite(y_all);

        x = x_all(valid);
        y = y_all(valid);

        stats = compute_pair_stats(x, y, min_n);

        assoc_row = group_row;

        assoc_row.feature_name = string(feature_input);
        assoc_row.feature_var = string(feature_var);
        assoc_row.q_var = string(q_var);

        assoc_row.n_total = numel(x_all);
        assoc_row.n_valid = numel(x);

        assoc_row.feature_min = stats.feature_min;
        assoc_row.feature_max = stats.feature_max;
        assoc_row.feature_range = stats.feature_range;
        assoc_row.feature_mean = stats.feature_mean;
        assoc_row.feature_std = stats.feature_std;

        assoc_row.q_min = stats.q_min;
        assoc_row.q_max = stats.q_max;
        assoc_row.q_range = stats.q_range;
        assoc_row.q_mean = stats.q_mean;
        assoc_row.q_std = stats.q_std;

        assoc_row.pearson_r = stats.pearson_r;
        assoc_row.pearson_p = stats.pearson_p;

        assoc_row.spearman_rho = stats.spearman_rho;
        assoc_row.spearman_p = stats.spearman_p;

        assoc_row.kendall_tau = stats.kendall_tau;
        assoc_row.kendall_p = stats.kendall_p;

        assoc_row.linear_slope = stats.linear_slope;
        assoc_row.linear_intercept = stats.linear_intercept;
        assoc_row.linear_R2 = stats.linear_R2;
        assoc_row.linear_RMSE = stats.linear_RMSE;

        T_assoc = [T_assoc; assoc_row]; %#ok<AGROW>

        T_bins_i = compute_binned_trend( ...
            x, ...
            y, ...
            num_bins, ...
            bin_mode);

        if ~isempty(T_bins_i)

            n_bin_rows = height(T_bins_i);

            if isempty(group_vars)

                group_bin = table();

            else

                group_bin = group_row(ones(n_bin_rows, 1), :);

            end

            group_bin.feature_name = repmat(string(feature_input), n_bin_rows, 1);
            group_bin.feature_var = repmat(string(feature_var), n_bin_rows, 1);
            group_bin.q_var = repmat(string(q_var), n_bin_rows, 1);

            T_bins_i = [group_bin T_bins_i];

            T_binned = [T_binned; T_bins_i]; %#ok<AGROW>

        end
    end
end

if verbose
    fprintf('\nGenerated association table rows: %d\n', height(T_assoc));
    fprintf('Generated binned trend rows   : %d\n', height(T_binned));
end

end

% =========================================================================
% Local helper functions
% =========================================================================

function feature_var = resolve_feature_variable(T, feature_name)

vars = string(T.Properties.VariableNames);
feature_name = string(feature_name);

candidates = [
    feature_name
    feature_name + "_mean"
];

idx = find(ismember(vars, candidates), 1);

if isempty(idx)
    error('Could not resolve feature variable: %s', feature_name);
end

feature_var = vars(idx);

end

function stats = compute_pair_stats(x, y, min_n)

stats = empty_stats();

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

n = numel(x);

if n == 0
    return;
end

stats.feature_min = min(x);
stats.feature_max = max(x);
stats.feature_range = stats.feature_max - stats.feature_min;
stats.feature_mean = mean(x);
stats.feature_std = std_safe(x);

stats.q_min = min(y);
stats.q_max = max(y);
stats.q_range = stats.q_max - stats.q_min;
stats.q_mean = mean(y);
stats.q_std = std_safe(y);

has_variation = ...
    n >= min_n && ...
    numel(unique(x)) >= 2 && ...
    numel(unique(y)) >= 2;

if ~has_variation
    return;
end

[stats.pearson_r, stats.pearson_p] = safe_corr(x, y, 'Pearson');
[stats.spearman_rho, stats.spearman_p] = safe_corr(x, y, 'Spearman');
[stats.kendall_tau, stats.kendall_p] = safe_corr(x, y, 'Kendall');

[ ...
    stats.linear_slope, ...
    stats.linear_intercept, ...
    stats.linear_R2, ...
    stats.linear_RMSE] = linear_fit_stats(x, y);

end

function stats = empty_stats()

stats = struct();

stats.feature_min = NaN;
stats.feature_max = NaN;
stats.feature_range = NaN;
stats.feature_mean = NaN;
stats.feature_std = NaN;

stats.q_min = NaN;
stats.q_max = NaN;
stats.q_range = NaN;
stats.q_mean = NaN;
stats.q_std = NaN;

stats.pearson_r = NaN;
stats.pearson_p = NaN;

stats.spearman_rho = NaN;
stats.spearman_p = NaN;

stats.kendall_tau = NaN;
stats.kendall_p = NaN;

stats.linear_slope = NaN;
stats.linear_intercept = NaN;
stats.linear_R2 = NaN;
stats.linear_RMSE = NaN;

end

function s = std_safe(x)

if numel(x) < 2
    s = NaN;
else
    s = std(x, 0);
end

end

function [r, p] = safe_corr(x, y, corr_type)

r = NaN;
p = NaN;

try
    [r, p] = corr(x(:), y(:), ...
        'Type', corr_type, ...
        'Rows', 'complete');
catch
    r = NaN;
    p = NaN;
end

end

function [slope, intercept, R2, RMSE] = linear_fit_stats(x, y)

slope = NaN;
intercept = NaN;
R2 = NaN;
RMSE = NaN;

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

if numel(x) < 2 || numel(unique(x)) < 2
    return;
end

p = polyfit(x, y, 1);

slope = p(1);
intercept = p(2);

y_fit = polyval(p, x);

residuals = y - y_fit;

RMSE = sqrt(mean(residuals.^2));

ss_res = sum(residuals.^2);
ss_tot = sum((y - mean(y)).^2);

if ss_tot > eps
    R2 = 1 - ss_res / ss_tot;
end

end

function T_bins = compute_binned_trend(x, y, num_bins, bin_mode)

T_bins = table();

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

n = numel(x);

if n < 2 || numel(unique(x)) < 2
    return;
end

switch bin_mode

    case "equal_count"

        bin_id = equal_count_bins(x, num_bins);

    case "equal_width"

        bin_id = equal_width_bins(x, num_bins);

    otherwise

        error('Unknown BinMode: %s', bin_mode);

end

valid_bins = unique(bin_id(isfinite(bin_id) & bin_id > 0));

if isempty(valid_bins)
    return;
end

n_bins = numel(valid_bins);

bin_number = zeros(n_bins, 1);
n_bin = zeros(n_bins, 1);

feature_mean = zeros(n_bins, 1);
feature_median = zeros(n_bins, 1);
feature_min = zeros(n_bins, 1);
feature_max = zeros(n_bins, 1);

q_mean = zeros(n_bins, 1);
q_median = zeros(n_bins, 1);
q_std = zeros(n_bins, 1);
q_min = zeros(n_bins, 1);
q_max = zeros(n_bins, 1);

for i = 1:n_bins

    b = valid_bins(i);

    mask = bin_id == b;

    xb = x(mask);
    yb = y(mask);

    bin_number(i) = b;
    n_bin(i) = numel(xb);

    feature_mean(i) = mean(xb);
    feature_median(i) = median(xb);
    feature_min(i) = min(xb);
    feature_max(i) = max(xb);

    q_mean(i) = mean(yb);
    q_median(i) = median(yb);
    q_std(i) = std_safe(yb);
    q_min(i) = min(yb);
    q_max(i) = max(yb);

end

T_bins = table( ...
    bin_number, ...
    n_bin, ...
    feature_mean, ...
    feature_median, ...
    feature_min, ...
    feature_max, ...
    q_mean, ...
    q_median, ...
    q_std, ...
    q_min, ...
    q_max);

end

function bin_id = equal_count_bins(x, num_bins)

n = numel(x);

num_bins = min(num_bins, n);
num_bins = max(num_bins, 1);

[~, sort_idx] = sort(x);

bin_id = zeros(n, 1);

for rank_idx = 1:n

    b = ceil(rank_idx * num_bins / n);
    b = min(max(b, 1), num_bins);

    original_idx = sort_idx(rank_idx);

    bin_id(original_idx) = b;

end

end

function bin_id = equal_width_bins(x, num_bins)

xmin = min(x);
xmax = max(x);

if xmax <= xmin
    bin_id = ones(size(x));
    return;
end

edges = linspace(xmin, xmax, num_bins + 1);

[~, ~, bin_id] = histcounts(x, edges);

bin_id(x == xmax) = num_bins;

end
