function T_err = summarize_q_model_errors_by_group(T_pred, varargin)
%SUMMARIZE_Q_MODEL_ERRORS_BY_GROUP Summarize q-prediction errors by group.
%
% Purpose
%   Analyze where a q-prediction model fails by computing error metrics
%   grouped by variables such as WaveModel, M, f0, cs_bg, or condition_id.
%
% Usage
%   T_err = adaptive_req.analysis.summarize_q_model_errors_by_group( ...
%       T_pred_model1, ...
%       'GroupVars', "REQ_M", ...
%       'Split', 'test');
%
% Output
%   T_err:
%       One row per model type and group.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.summarize_q_model_errors_by_group';

addRequired(p, 'T_pred', @istable);

addParameter(p, 'GroupVars', "REQ_M", ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'Split', 'test', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'ModelTypes', string.empty(1, 0), ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

parse(p, T_pred, varargin{:});

group_vars = string(p.Results.GroupVars);
split_name = string(p.Results.Split);
model_types = string(p.Results.ModelTypes);

required_vars = ["model_type", "split", "q_true", "q_pred", "residual", "abs_error"];

for i = 1:numel(required_vars)
    if ~ismember(required_vars(i), string(T_pred.Properties.VariableNames))
        error('Required variable not found in T_pred: %s', required_vars(i));
    end
end

group_vars = group_vars(group_vars ~= "");

for i = 1:numel(group_vars)
    if ~ismember(group_vars(i), string(T_pred.Properties.VariableNames))
        error('Group variable not found in T_pred: %s', group_vars(i));
    end
end

T_use = T_pred(T_pred.split == split_name, :);

if ~isempty(model_types)
    T_use = T_use(ismember(T_use.model_type, model_types), :);
end

if isempty(T_use)
    error('No rows found for split = %s.', split_name);
end

all_group_vars = ["model_type", group_vars];

[G, T_err] = findgroups(T_use(:, cellstr(all_group_vars)));

q_true = T_use.q_true;
q_pred = T_use.q_pred;
residual = T_use.residual;
abs_error = T_use.abs_error;

T_err.n = splitapply(@numel, q_true, G);

T_err.MAE = splitapply(@mean_omitnan, abs_error, G);
T_err.RMSE = splitapply(@rmse_from_residual, residual, G);
T_err.bias = splitapply(@mean_omitnan, residual, G);
T_err.residual_std = splitapply(@std_omitnan, residual, G);

T_err.median_abs_error = splitapply(@median_omitnan, abs_error, G);
T_err.p90_abs_error = splitapply(@p90_omitnan, abs_error, G);
T_err.max_abs_error = splitapply(@max_omitnan, abs_error, G);

T_err.R2 = splitapply(@r2_score, q_true, q_pred, G);
T_err.pearson_r = splitapply(@(x, y) safe_corr(x, y, 'Pearson'), q_true, q_pred, G);
T_err.spearman_rho = splitapply(@(x, y) safe_corr(x, y, 'Spearman'), q_true, q_pred, G);

T_err.q_true_min = splitapply(@min_omitnan, q_true, G);
T_err.q_true_max = splitapply(@max_omitnan, q_true, G);
T_err.q_true_mean = splitapply(@mean_omitnan, q_true, G);

T_err.q_pred_min = splitapply(@min_omitnan, q_pred, G);
T_err.q_pred_max = splitapply(@max_omitnan, q_pred, G);
T_err.q_pred_mean = splitapply(@mean_omitnan, q_pred, G);

T_err.frac_overestimate = splitapply(@(x) mean(x > 0, 'omitnan'), residual, G);
T_err.frac_underestimate = splitapply(@(x) mean(x < 0, 'omitnan'), residual, G);

T_err.split = repmat(split_name, height(T_err), 1);
T_err = movevars(T_err, 'split', 'After', 'model_type');

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

function y = p90_omitnan(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = prctile(x, 90);
end
end

function y = rmse_from_residual(residual)
residual = residual(isfinite(residual));
if isempty(residual)
    y = NaN;
else
    y = sqrt(mean(residual.^2));
end
end

function R2 = r2_score(y_true, y_pred)

valid = isfinite(y_true) & isfinite(y_pred);

y_true = y_true(valid);
y_pred = y_pred(valid);

if numel(y_true) < 2
    R2 = NaN;
    return;
end

ss_res = sum((y_true - y_pred).^2);
ss_tot = sum((y_true - mean(y_true)).^2);

if ss_tot <= eps
    R2 = NaN;
else
    R2 = 1 - ss_res / ss_tot;
end

end

function r = safe_corr(x, y, corr_type)

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

if numel(x) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    r = NaN;
    return;
end

try
    r = corr(x(:), y(:), 'Type', corr_type, 'Rows', 'complete');
catch
    r = NaN;
end

end
