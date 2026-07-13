function T_condition = make_condition_summary_table(T_step, varargin)
%MAKE_CONDITION_SUMMARY_TABLE Build a condition-level summary table.
%
% Usage
%   T_condition = adaptive_req.analysis.make_condition_summary_table(T_step);
%
% Input
%   T_step:
%       Table produced by make_step_summary_table.
%
% Output
%   T_condition:
%       One row per condition_id.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.make_condition_summary_table';

addRequired(p, 'T_step', @istable);

addParameter(p, 'ConditionVar', 'condition_id', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'OmegaVar', 'omega_mean', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'QVar', 'q_mean', ...
    @(x) ischar(x) || isstring(x));

parse(p, T_step, varargin{:});

condition_var = char(p.Results.ConditionVar);
omega_var = char(p.Results.OmegaVar);
q_var = char(p.Results.QVar);

vars = string(T_step.Properties.VariableNames);

if ~ismember(condition_var, vars)
    error('Condition variable not found in T_step: %s', condition_var);
end

if ~ismember(q_var, vars)
    error('q variable not found in T_step: %s', q_var);
end

has_omega = ismember(omega_var, vars);

[G, T_condition] = findgroups(T_step(:, {condition_var}));

T_condition.n_steps = splitapply(@numel, ones(height(T_step), 1), G);

param_vars = vars(startsWith(vars, "SIM_") | startsWith(vars, "REQ_"));

for i = 1:numel(param_vars)

    v = param_vars(i);

    if ismember(v, string(T_condition.Properties.VariableNames))
        continue;
    end

    T_condition.(char(v)) = splitapply( ...
        @(x) first_value(x), ...
        T_step.(char(v)), ...
        G);
end

q = T_step.(q_var);

T_condition.q_mean_over_steps = splitapply(@mean_omitnan, q, G);
T_condition.q_median_over_steps = splitapply(@median_omitnan, q, G);
T_condition.q_std_over_steps = splitapply(@std_omitnan, q, G);
T_condition.q_min_over_steps = splitapply(@min_omitnan, q, G);
T_condition.q_max_over_steps = splitapply(@max_omitnan, q, G);
T_condition.q_range_over_steps = ...
    T_condition.q_max_over_steps - T_condition.q_min_over_steps;

if has_omega

    omega = T_step.(omega_var);

    T_condition.omega_min = splitapply(@min_omitnan, omega, G);
    T_condition.omega_max = splitapply(@max_omitnan, omega, G);
    T_condition.omega_range = T_condition.omega_max - T_condition.omega_min;

    T_condition.q_vs_omega_slope = splitapply( ...
        @linear_slope, ...
        omega, ...
        q, ...
        G);

    T_condition.q_vs_omega_r2 = splitapply( ...
        @linear_r2, ...
        omega, ...
        q, ...
        G);

else

    T_condition.omega_min = NaN(height(T_condition), 1);
    T_condition.omega_max = NaN(height(T_condition), 1);
    T_condition.omega_range = NaN(height(T_condition), 1);
    T_condition.q_vs_omega_slope = NaN(height(T_condition), 1);
    T_condition.q_vs_omega_r2 = NaN(height(T_condition), 1);

end

% Summarize feature mean columns if present.
feature_mean_vars = vars(endsWith(vars, "_mean"));
feature_mean_vars = feature_mean_vars(~ismember(feature_mean_vars, ...
    ["omega_mean", "q_mean"]));

for i = 1:numel(feature_mean_vars)

    v = feature_mean_vars(i);
    x = T_step.(char(v));

    if ~isnumeric(x)
        continue;
    end

    out_name = char("condition_" + v);

    T_condition.(out_name) = splitapply(@mean_omitnan, x, G);

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

function slope = linear_slope(x, y)

mask = isfinite(x) & isfinite(y);

x = x(mask);
y = y(mask);

if numel(x) < 2 || numel(unique(x)) < 2
    slope = NaN;
    return;
end

p = polyfit(x, y, 1);
slope = p(1);

end

function r2 = linear_r2(x, y)

mask = isfinite(x) & isfinite(y);

x = x(mask);
y = y(mask);

if numel(x) < 2 || numel(unique(x)) < 2
    r2 = NaN;
    return;
end

p = polyfit(x, y, 1);
yfit = polyval(p, x);

ss_res = sum((y - yfit).^2);
ss_tot = sum((y - mean(y)).^2);

if ss_tot <= eps
    r2 = NaN;
else
    r2 = 1 - ss_res / ss_tot;
end

end
