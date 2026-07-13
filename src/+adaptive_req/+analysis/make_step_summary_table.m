function T_step = make_step_summary_table(T, COL, feature_vars)
%MAKE_STEP_SUMMARY_TABLE Build a condition-step summary table.
%
% Usage
%   T_step = adaptive_req.analysis.make_step_summary_table( ...
%       T_mc, COL, feature_vars);
%
% The output has one row per:
%
%   condition_id × step_idx
%
% and contains q statistics, Omega statistics, parameter values, and feature
% summaries.

if nargin < 2 || isempty(COL)
    COL = adaptive_req.analysis.detect_mc_columns(T);
end

if nargin < 3 || isempty(feature_vars)
    feature_vars = adaptive_req.analysis.detect_feature_variables(T);
end

if COL.condition == ""
    error('Condition column was not detected.');
end

if COL.step == ""
    error('Step column was not detected.');
end

if COL.q == ""
    error('q column was not detected.');
end

group_vars = cellstr([COL.condition, COL.step]);

[G, T_step] = findgroups(T(:, group_vars));

T_step.n_rows = splitapply(@numel, ones(height(T), 1), G);

% Add parameter columns by taking the first value in each group.
param_vars = string(COL.params);

for i = 1:numel(param_vars)

    v = param_vars(i);

    if ismember(v, string(T_step.Properties.VariableNames))
        continue;
    end

    T_step.(char(v)) = splitapply( ...
        @(x) first_value(x), ...
        T.(char(v)), ...
        G);
end

% Add Omega statistics.
if COL.omega ~= ""

    omega = T.(char(COL.omega));

    T_step.omega_mean = splitapply(@mean_omitnan, omega, G);
    T_step.omega_std = splitapply(@std_omitnan, omega, G);
    T_step.omega_min = splitapply(@min_omitnan, omega, G);
    T_step.omega_max = splitapply(@max_omitnan, omega, G);

end

% Add q statistics.
q = T.(char(COL.q));

T_step.q_mean = splitapply(@mean_omitnan, q, G);
T_step.q_median = splitapply(@median_omitnan, q, G);
T_step.q_std = splitapply(@std_omitnan, q, G);
T_step.q_min = splitapply(@min_omitnan, q, G);
T_step.q_max = splitapply(@max_omitnan, q, G);
T_step.q_range = T_step.q_max - T_step.q_min;

% Add feature summaries.
feature_vars = string(feature_vars);

for i = 1:numel(feature_vars)

    v = feature_vars(i);

    if ~ismember(v, string(T.Properties.VariableNames))
        continue;
    end

    x = T.(char(v));

    if ~isnumeric(x)
        continue;
    end

    base_name = char(v);

    T_step.([base_name '_mean']) = splitapply(@mean_omitnan, x, G);
    T_step.([base_name '_median']) = splitapply(@median_omitnan, x, G);
    T_step.([base_name '_std']) = splitapply(@std_omitnan, x, G);
    T_step.([base_name '_min']) = splitapply(@min_omitnan, x, G);
    T_step.([base_name '_max']) = splitapply(@max_omitnan, x, G);

end

end

function y = first_value(x)

if iscell(x)
    y = x{1};
else
    y = x(1);
end

end

function y = mean_omitnan(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = mean(x);
end

end

function y = median_omitnan(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = median(x);
end

end

function y = std_omitnan(x)

x = x(isfinite(x));

if numel(x) < 2
    y = NaN;
else
    y = std(x, 0);
end

end

function y = min_omitnan(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = min(x);
end

end

function y = max_omitnan(x)

x = x(isfinite(x));

if isempty(x)
    y = NaN;
else
    y = max(x);
end

end
