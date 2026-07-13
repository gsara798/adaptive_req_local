function [MODEL, T_pred, T_metrics] = train_q_model_features_only(T, feature_vars, varargin)
%TRAIN_Q_MODEL_FEATURES_ONLY Train feature-only models for local q prediction.
%
% Purpose
%   Train models of the form:
%
%       q = F(local spectral features)
%
%   No physical or acquisition parameters are used as predictors.
%
% Supported model types
%   linear
%   quadratic
%   knn
%
% Usage
%   [MODEL, T_pred, T_metrics] = adaptive_req.analysis.train_q_model_features_only( ...
%       T_mc, ...
%       ["ang_entropy", "radial_entropy"], ...
%       'QVar', 'q_theory', ...
%       'SplitMode', 'condition');

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.train_q_model_features_only';

addRequired(p, 'T', @istable);

addRequired(p, 'feature_vars', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'QVar', 'q_theory', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'SplitMode', 'condition', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'ConditionVar', 'condition_id', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'TrainFraction', 0.70, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);

addParameter(p, 'RandomSeed', 1701, ...
    @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'ModelTypes', ["linear", "quadratic", "knn"], ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'KNNK', 15, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'ClipPredictions', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'Verbose', true, ...
    @(x) islogical(x) || isnumeric(x));

parse(p, T, feature_vars, varargin{:});

feature_vars = string(feature_vars);
q_var = char(p.Results.QVar);
split_mode = lower(string(p.Results.SplitMode));
condition_var = char(p.Results.ConditionVar);
train_fraction = p.Results.TrainFraction;
random_seed = p.Results.RandomSeed;
model_types = lower(string(p.Results.ModelTypes));
knn_k = round(p.Results.KNNK);
clip_predictions = logical(p.Results.ClipPredictions);
verbose = logical(p.Results.Verbose);

available_vars = string(T.Properties.VariableNames);

if ~ismember(q_var, available_vars)
    error('QVar not found in table: %s', q_var);
end

resolved_features = strings(numel(feature_vars), 1);

for i = 1:numel(feature_vars)
    resolved_features(i) = resolve_feature_variable(T, feature_vars(i));
end

X_raw = zeros(height(T), numel(resolved_features));

for i = 1:numel(resolved_features)
    X_raw(:, i) = T.(char(resolved_features(i)));
end

y_raw = T.(q_var);

valid = isfinite(y_raw);

for i = 1:size(X_raw, 2)
    valid = valid & isfinite(X_raw(:, i));
end

T_model = T(valid, :);
X = X_raw(valid, :);
y = y_raw(valid);

n = numel(y);

if n < 10
    error('Not enough valid rows for training. Valid rows: %d', n);
end

rng(random_seed);

is_train = make_train_split( ...
    T_model, ...
    split_mode, ...
    condition_var, ...
    train_fraction);

is_test = ~is_train;

if ~any(is_train) || ~any(is_test)
    error('Invalid train-test split. Check SplitMode and TrainFraction.');
end

mu = mean(X(is_train, :), 1, 'omitnan');
sigma = std(X(is_train, :), 0, 1, 'omitnan');

sigma(sigma <= eps | ~isfinite(sigma)) = 1;

Xz = (X - mu) ./ sigma;

if verbose
    fprintf('\nTraining feature-only q model.\n');
    fprintf('Q variable     : %s\n', q_var);
    fprintf('Split mode     : %s\n', split_mode);
    fprintf('Train fraction : %.2f\n', train_fraction);
    fprintf('Valid rows     : %d\n', n);
    fprintf('Train rows     : %d\n', sum(is_train));
    fprintf('Test rows      : %d\n', sum(is_test));

    fprintf('\nFeatures used as predictors:\n');
    disp(resolved_features(:));
end

MODEL = struct();

MODEL.model_family = "features_only";
MODEL.q_var = string(q_var);
MODEL.input_feature_vars = feature_vars(:);
MODEL.resolved_feature_vars = resolved_features(:);
MODEL.split_mode = split_mode;
MODEL.condition_var = string(condition_var);
MODEL.train_fraction = train_fraction;
MODEL.random_seed = random_seed;
MODEL.mu = mu;
MODEL.sigma = sigma;
MODEL.clip_predictions = clip_predictions;
MODEL.models = struct([]);

T_pred = table();
T_metrics = table();

for midx = 1:numel(model_types)

    model_type = model_types(midx);

    model_i = train_one_model(model_type, Xz(is_train, :), y(is_train), knn_k);

    y_pred_raw = predict_one_model(model_i, Xz, Xz(is_train, :), y(is_train));

    if clip_predictions
        y_pred = min(max(y_pred_raw, 0), 1);
    else
        y_pred = y_pred_raw;
    end

    MODEL.models(midx).model_type = model_type;
    MODEL.models(midx).model = model_i;

    T_pred_i = build_prediction_table( ...
        T_model, ...
        resolved_features, ...
        y, ...
        y_pred_raw, ...
        y_pred, ...
        is_train, ...
        model_type);

    T_metrics_i = compute_prediction_metrics_table( ...
        y, ...
        y_pred, ...
        is_train, ...
        model_type);

    T_pred = [T_pred; T_pred_i]; %#ok<AGROW>
    T_metrics = [T_metrics; T_metrics_i]; %#ok<AGROW>

end

if verbose
    fprintf('\nPrediction metrics.\n');
    disp(T_metrics);
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

function is_train = make_train_split(T, split_mode, condition_var, train_fraction)

n = height(T);

switch split_mode

    case "random_rows"

        idx = randperm(n);
        n_train = round(train_fraction * n);

        is_train = false(n, 1);
        is_train(idx(1:n_train)) = true;

    case "condition"

        if ~ismember(condition_var, T.Properties.VariableNames)
            error('ConditionVar not found in table: %s', condition_var);
        end

        condition_values = unique(T.(condition_var), 'stable');
        n_conditions = numel(condition_values);

        idx = randperm(n_conditions);
        n_train_conditions = round(train_fraction * n_conditions);

        train_conditions = condition_values(idx(1:n_train_conditions));

        is_train = ismember(T.(condition_var), train_conditions);

    otherwise

        error('Unknown SplitMode: %s', split_mode);

end

end

function model = train_one_model(model_type, Xtrain, ytrain, knn_k)

model = struct();
model.model_type = model_type;

switch model_type

    case "linear"

        Phi = build_design_matrix(Xtrain, "linear");
        beta = Phi \ ytrain;

        model.beta = beta;
        model.design = "linear";

    case "quadratic"

        Phi = build_design_matrix(Xtrain, "quadratic");
        beta = Phi \ ytrain;

        model.beta = beta;
        model.design = "quadratic";

    case "knn"

        model.k = knn_k;
        model.design = "knn";

    otherwise

        error('Unknown model type: %s', model_type);

end

end

function y_pred = predict_one_model(model, X, Xtrain, ytrain)

switch model.model_type

    case {"linear", "quadratic"}

        Phi = build_design_matrix(X, model.design);
        y_pred = Phi * model.beta;

    case "knn"

        y_pred = knn_predict(Xtrain, ytrain, X, model.k);

    otherwise

        error('Unknown model type: %s', model.model_type);

end

end

function Phi = build_design_matrix(X, design)

n = size(X, 1);
p = size(X, 2);

switch design

    case "linear"

        Phi = [ones(n, 1), X];

    case "quadratic"

        Phi = [ones(n, 1), X, X.^2];

        for i = 1:p
            for j = i+1:p
                Phi = [Phi, X(:, i) .* X(:, j)]; %#ok<AGROW>
            end
        end

    otherwise

        error('Unknown design: %s', design);

end

end

function y_pred = knn_predict(Xtrain, ytrain, Xquery, k)

n_query = size(Xquery, 1);
n_train = size(Xtrain, 1);

k = min(k, n_train);

y_pred = NaN(n_query, 1);

chunk_size = 1000;

for start_idx = 1:chunk_size:n_query

    stop_idx = min(start_idx + chunk_size - 1, n_query);

    Xq = Xquery(start_idx:stop_idx, :);

    D2 = zeros(size(Xq, 1), n_train);

    for d = 1:size(Xtrain, 2)
        D2 = D2 + (Xq(:, d) - Xtrain(:, d)').^2;
    end

    [~, idx_sorted] = sort(D2, 2, 'ascend');

    idx_knn = idx_sorted(:, 1:k);

    for i = 1:size(idx_knn, 1)
        y_pred(start_idx + i - 1) = mean(ytrain(idx_knn(i, :)), 'omitnan');
    end
end

end

function T_pred = build_prediction_table( ...
    T_model, ...
    feature_vars, ...
    y_true, ...
    y_pred_raw, ...
    y_pred, ...
    is_train, ...
    model_type)

keep_candidates = string({ ...
    'condition_id', ...
    'condition_position', ...
    'condition_label', ...
    'step_idx', ...
    'realization_idx', ...
    'patch_idx', ...
    'Omega_sr', ...
    'Omega', ...
    'omega_sr', ...
    'aperture_value', ...
    'SIM_WaveModel', ...
    'SIM_f0', ...
    'SIM_cs_bg', ...
    'REQ_M'});

vars = string(T_model.Properties.VariableNames);

keep_vars = keep_candidates(ismember(keep_candidates, vars));
keep_vars = unique([keep_vars(:); feature_vars(:)], 'stable');

T_pred = T_model(:, cellstr(keep_vars));

T_pred.model_type = repmat(string(model_type), height(T_model), 1);
T_pred.split = strings(height(T_model), 1);
T_pred.split(is_train) = "train";
T_pred.split(~is_train) = "test";

T_pred.q_true = y_true;
T_pred.q_pred_raw = y_pred_raw;
T_pred.q_pred = y_pred;
T_pred.residual = y_pred - y_true;
T_pred.abs_error = abs(T_pred.residual);

end

function T_metrics = compute_prediction_metrics_table(y_true, y_pred, is_train, model_type)

splits = ["train", "test", "all"];

T_metrics = table();

for i = 1:numel(splits)

    split_i = splits(i);

    switch split_i
        case "train"
            mask = is_train;
        case "test"
            mask = ~is_train;
        case "all"
            mask = true(size(y_true));
    end

    yt = y_true(mask);
    yp = y_pred(mask);

    valid = isfinite(yt) & isfinite(yp);

    yt = yt(valid);
    yp = yp(valid);

    row = compute_metrics_row(yt, yp);

    row.model_type = string(model_type);
    row.split = split_i;

    T_metrics = [T_metrics; row]; %#ok<AGROW>
end

T_metrics = movevars(T_metrics, {'model_type', 'split'}, 'Before', 1);

end

function row = compute_metrics_row(y_true, y_pred)

n = numel(y_true);

err = y_pred - y_true;

MAE = mean(abs(err), 'omitnan');
RMSE = sqrt(mean(err.^2, 'omitnan'));
bias = mean(err, 'omitnan');

q_std = std_safe(y_true);
err_std = std_safe(err);

R2 = compute_r2(y_true, y_pred);

pearson_r = safe_corr(y_true, y_pred, 'Pearson');
spearman_rho = safe_corr(y_true, y_pred, 'Spearman');

q_true_min = min(y_true, [], 'omitnan');
q_true_max = max(y_true, [], 'omitnan');
q_pred_min = min(y_pred, [], 'omitnan');
q_pred_max = max(y_pred, [], 'omitnan');

row = table( ...
    n, ...
    MAE, ...
    RMSE, ...
    bias, ...
    q_std, ...
    err_std, ...
    R2, ...
    pearson_r, ...
    spearman_rho, ...
    q_true_min, ...
    q_true_max, ...
    q_pred_min, ...
    q_pred_max);

end

function R2 = compute_r2(y_true, y_pred)

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

r = NaN;

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

if numel(x) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    return;
end

try
    r = corr(x(:), y(:), 'Type', corr_type, 'Rows', 'complete');
catch
    r = NaN;
end

end

function s = std_safe(x)

x = x(isfinite(x));

if numel(x) < 2
    s = NaN;
else
    s = std(x, 0);
end

end
