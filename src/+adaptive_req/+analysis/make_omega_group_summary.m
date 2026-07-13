function T_omega = make_omega_group_summary(T_step, varargin)
%MAKE_OMEGA_GROUP_SUMMARY Average T_step over repeated conditions.
%
% Default grouping:
%   SIM_WaveModel, REQ_M, omega_mean
%
% This is useful for plotting average q versus aperture without connecting
% unrelated condition-step points as if they were one continuous curve.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.make_omega_group_summary';

addRequired(p, 'T_step', @istable);

addParameter(p, 'GroupVars', ...
    ["SIM_WaveModel", "REQ_M", "omega_mean"], ...
    @(x) isstring(x) || iscellstr(x) || ischar(x));

addParameter(p, 'QVar', 'q_mean', ...
    @(x) ischar(x) || isstring(x));

parse(p, T_step, varargin{:});

group_vars = string(p.Results.GroupVars);
q_var = char(p.Results.QVar);

available_vars = string(T_step.Properties.VariableNames);

group_vars = group_vars(ismember(group_vars, available_vars));

if isempty(group_vars)
    error('No valid grouping variables were found.');
end

if ~ismember(q_var, available_vars)
    error('q variable not found: %s', q_var);
end

[G, T_omega] = findgroups(T_step(:, cellstr(group_vars)));

q = T_step.(q_var);

T_omega.n_rows = splitapply(@numel, q, G);
T_omega.q_mean = splitapply(@mean_omitnan, q, G);
T_omega.q_std = splitapply(@std_omitnan, q, G);
T_omega.q_median = splitapply(@median_omitnan, q, G);
T_omega.q_min = splitapply(@min_omitnan, q, G);
T_omega.q_max = splitapply(@max_omitnan, q, G);
T_omega.q_range = T_omega.q_max - T_omega.q_min;

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