function [MODEL, T_pred, T_metrics] = train_q_model_fixed_split( ...
    T, predictors, train_mask, test_mask, varargin)
%TRAIN_Q_MODEL_FIXED_SPLIT Train q models with explicit train/test masks.
%
% This function is intentionally strict: the caller owns the split. It is
% useful for leave-one-group-out tests where random window-level splits would
% leak repeated global features across train/test rows.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.train_q_model_fixed_split';

addRequired(p, 'T', @istable);
addRequired(p, 'predictors', @(x) ischar(x) || isstring(x) || iscellstr(x));
addRequired(p, 'train_mask', @(x) islogical(x) || isnumeric(x));
addRequired(p, 'test_mask', @(x) islogical(x) || isnumeric(x));

addParameter(p, 'QVar', 'q_theory', @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelName', 'q_model', @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelRole', 'operational', @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelTypes', ["linear", "boosted_trees", "bagged_trees"], ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'NumLearningCycles', 200, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'MinLeafSize', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'LearnRate', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'UseParallel', false, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ClipRange', [0.001 0.999], ...
    @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, T, predictors, train_mask, test_mask, varargin{:});

predictors = string(predictors(:));
q_var = string(p.Results.QVar);
model_name = string(p.Results.ModelName);
model_role = string(p.Results.ModelRole);
clip_range = double(p.Results.ClipRange);

train_mask = logical(train_mask(:));
test_mask = logical(test_mask(:));

if numel(train_mask) ~= height(T) || numel(test_mask) ~= height(T)
    error('train_mask and test_mask must have height(T) elements.');
end
if any(train_mask & test_mask)
    error('train_mask and test_mask must be disjoint.');
end
if ~any(train_mask) || ~any(test_mask)
    error('train_mask and test_mask must both be non-empty.');
end

assert_no_leakage_predictors(predictors, model_role);

[MODEL, T_pred, ~] = adaptive_req.analysis.train_q_model_from_predictors( ...
    T, predictors, ...
    'QVar', q_var, ...
    'ModelName', model_name, ...
    'TrainMask', train_mask, ...
    'ModelTypes', p.Results.ModelTypes, ...
    'NumLearningCycles', p.Results.NumLearningCycles, ...
    'MinLeafSize', p.Results.MinLeafSize, ...
    'LearnRate', p.Results.LearnRate, ...
    'UseParallel', p.Results.UseParallel, ...
    'ClipPredictions', false, ...
    'Verbose', p.Results.Verbose);

T_pred.q_pred = min(max(T_pred.q_pred_raw, clip_range(1)), clip_range(2));
T_pred.residual = T_pred.q_pred - T_pred.q_true;
T_pred.abs_error = abs(T_pred.residual);
T_pred.model_role = repmat(model_role, height(T_pred), 1);

if ~ismember('row_key', string(T_pred.Properties.VariableNames))
    T_pred.row_key = make_row_key(T_pred);
end

MODEL.model_role = model_role;
MODEL.fixed_split = true;
MODEL.clip_range = clip_range;

T_metrics = compute_q_metrics(T_pred);

end

function assert_no_leakage_predictors(predictors, model_role)

if string(model_role) ~= "operational"
    return;
end

predictors = lower(string(predictors(:)));

exact_banned = lower([
    "q_theory"
    "q_global_theory"
    "q_local_minus_global"
    "cs_true"
    "cs_pred"
    "sws_error"
    "abs_sws_error"
    "sws_error_pct"
    "abs_sws_error_pct"
    "q_true"
    "q_pred"
    "q_pred_raw"
    "residual"
    "abs_error"
    "target"]);

bad_exact = intersect(predictors, exact_banned);
if ~isempty(bad_exact)
    error('Operational predictors contain leakage variables: %s', ...
        strjoin(bad_exact, ', '));
end

bad_patterns = contains(predictors, "error") | ...
    contains(predictors, "residual") | ...
    contains(predictors, "target");
if any(bad_patterns)
    error('Operational predictors contain target/error-like variables: %s', ...
        strjoin(predictors(bad_patterns), ', '));
end

end

function T_metrics = compute_q_metrics(T_pred)

groups = ["model_name", "model_type", "model_role", "split"];
[G, T_metrics] = findgroups(T_pred(:, cellstr(groups)));

T_metrics.N = splitapply(@numel, T_pred.q_true, G);
T_metrics.MAE_q = splitapply(@(a, b) mean(abs(b - a), 'omitnan'), ...
    T_pred.q_true, T_pred.q_pred, G);
T_metrics.RMSE_q = splitapply(@(a, b) sqrt(mean((b - a).^2, 'omitnan')), ...
    T_pred.q_true, T_pred.q_pred, G);
T_metrics.bias_q = splitapply(@(a, b) mean(b - a, 'omitnan'), ...
    T_pred.q_true, T_pred.q_pred, G);
T_metrics.Pearson = splitapply(@(a, b) safe_corr(a, b, 'Pearson'), ...
    T_pred.q_true, T_pred.q_pred, G);
T_metrics.Spearman = splitapply(@(a, b) safe_corr(a, b, 'Spearman'), ...
    T_pred.q_true, T_pred.q_pred, G);

end

function r = safe_corr(x, y, corr_type)

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 3 || std(x) <= eps || std(y) <= eps
    r = NaN;
    return;
end

r = corr(x(:), y(:), 'Type', corr_type, 'Rows', 'complete');

end

function key = make_row_key(T)

required = ["condition_id", "step_idx", "realization_idx", "patch_idx"];
if ~all(ismember(required, string(T.Properties.VariableNames)))
    key = strings(height(T), 1);
    return;
end

parts = strings(height(T), 4);
parts(:, 1) = string(T.condition_id);
parts(:, 2) = string(T.step_idx);
parts(:, 3) = string(T.realization_idx);
parts(:, 4) = string(T.patch_idx);
key = join(parts, "|", 2);

end
