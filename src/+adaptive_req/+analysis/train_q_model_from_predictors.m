function [MODEL, T_pred, T_metrics] = train_q_model_from_predictors(T, predictor_vars, varargin)
%TRAIN_Q_MODEL_FROM_PREDICTORS Train q-prediction models from selected predictors.
%
% Purpose
%   Train models of the form:
%
%       q = F(selected predictors)
%
%   The predictors are fully controlled by the calling script.
%
% Examples
%   Model 1:
%       predictors = ["ang_entropy", "radial_entropy"];
%
%   Model 2:
%       predictors = ["ang_entropy", "radial_entropy", "REQ_M", "SIM_f0"];
%
%   Diagnostic model:
%       predictors = ["ang_entropy", "radial_entropy", "REQ_M", "SIM_f0", "SIM_cs_bg"];
%
% Supported model types
%   linear
%   quadratic
%   knn

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.train_q_model_from_predictors';

addRequired(p, 'T', @istable);

addRequired(p, 'predictor_vars', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

addParameter(p, 'QVar', 'q_theory', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'ModelName', 'q_model', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'SplitMode', 'condition', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'ConditionVar', 'condition_id', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'TrainFraction', 0.70, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);

addParameter(p, 'TrainMask', [], ...
    @(x) isempty(x) || islogical(x) || isnumeric(x));

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

addParameter(p, 'NumLearningCycles', 200, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'MinLeafSize', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'LearnRate', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'UseParallel', false, ...
    @(x) islogical(x) || isnumeric(x));

parse(p, T, predictor_vars, varargin{:});

predictor_vars = string(predictor_vars);
q_var = char(p.Results.QVar);
model_name = string(p.Results.ModelName);
split_mode = lower(string(p.Results.SplitMode));
condition_var = char(p.Results.ConditionVar);
train_fraction = p.Results.TrainFraction;
train_mask = p.Results.TrainMask;
random_seed = p.Results.RandomSeed;
model_types = lower(string(p.Results.ModelTypes));
knn_k = round(p.Results.KNNK);
clip_predictions = logical(p.Results.ClipPredictions);
verbose = logical(p.Results.Verbose);
num_learning_cycles = round(p.Results.NumLearningCycles);
min_leaf_size = round(p.Results.MinLeafSize);
learn_rate = p.Results.LearnRate;
use_parallel = logical(p.Results.UseParallel);

train_options = struct();
train_options.knn_k = knn_k;
train_options.num_learning_cycles = num_learning_cycles;
train_options.min_leaf_size = min_leaf_size;
train_options.learn_rate = learn_rate;
train_options.use_parallel = use_parallel;

available_vars = string(T.Properties.VariableNames);

if ~ismember(q_var, available_vars)
    error('QVar not found in table: %s', q_var);
end

resolved_predictors = strings(numel(predictor_vars), 1);

for i = 1:numel(predictor_vars)
    resolved_predictors(i) = resolve_table_variable(T, predictor_vars(i));
end

valid = isfinite(T.(q_var));

for i = 1:numel(resolved_predictors)
    valid = valid & is_valid_predictor_column(T.(char(resolved_predictors(i))));
end

T_model = T(valid, :);
y = T_model.(q_var);

if height(T_model) < 10
    error('Not enough valid rows for training. Valid rows: %d', height(T_model));
end

[X_raw, encoder] = encode_predictors(T_model, resolved_predictors);

rng(random_seed);

if isempty(train_mask)
    is_train = make_train_split( ...
        T_model, ...
        split_mode, ...
        condition_var, ...
        train_fraction);
else
    is_train = resolve_train_mask(train_mask, valid, height(T_model));
    split_mode = "provided";
end

is_test = ~is_train;

if ~any(is_train) || ~any(is_test)
    error('Invalid train-test split. Check SplitMode and TrainFraction.');
end

mu = mean(X_raw(is_train, :), 1, 'omitnan');
sigma = std(X_raw(is_train, :), 0, 1, 'omitnan');

sigma(sigma <= eps | ~isfinite(sigma)) = 1;

X = (X_raw - mu) ./ sigma;

if verbose
    fprintf('\nTraining q model from selected predictors.\n');
    fprintf('Model name     : %s\n', model_name);
    fprintf('Q variable     : %s\n', q_var);
    fprintf('Split mode     : %s\n', split_mode);
    fprintf('Train fraction : %.2f\n', train_fraction);
    fprintf('Valid rows     : %d\n', height(T_model));
    fprintf('Train rows     : %d\n', sum(is_train));
    fprintf('Test rows      : %d\n', sum(is_test));

    fprintf('\nInput predictors:\n');
    disp(table(predictor_vars(:), resolved_predictors(:), ...
        'VariableNames', {'InputName', 'ResolvedName'}));

    fprintf('\nEncoded predictor columns:\n');
    disp(encoder.encoded_names(:));
end

MODEL = struct();

MODEL.model_name = model_name;
MODEL.q_var = string(q_var);
MODEL.input_predictor_vars = predictor_vars(:);
MODEL.resolved_predictor_vars = resolved_predictors(:);
MODEL.encoded_predictor_names = encoder.encoded_names(:);
MODEL.encoder = encoder;
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

    model_i = train_one_model( ...
    model_type, ...
    X(is_train, :), ...
    y(is_train), ...
    train_options);

    y_pred_raw = predict_one_model(model_i, X, X(is_train, :), y(is_train));

    if clip_predictions
        y_pred = min(max(y_pred_raw, 0), 1);
    else
        y_pred = y_pred_raw;
    end

    MODEL.models(midx).model_type = model_type;
    MODEL.models(midx).model = model_i;

    T_pred_i = build_prediction_table( ...
        T_model, ...
        resolved_predictors, ...
        y, ...
        y_pred_raw, ...
        y_pred, ...
        is_train, ...
        model_name, ...
        model_type);

    T_metrics_i = compute_prediction_metrics_table( ...
        y, ...
        y_pred, ...
        is_train, ...
        model_name, ...
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
% Variable handling
% =========================================================================

function var_name = resolve_table_variable(T, input_name)

vars = string(T.Properties.VariableNames);
input_name = string(input_name);

candidates = [
    input_name
    input_name + "_mean"
];

idx = find(ismember(vars, candidates), 1);

if isempty(idx)
    error('Could not resolve table variable: %s', input_name);
end

var_name = vars(idx);

end

function valid = is_valid_predictor_column(x)

if isnumeric(x) || islogical(x)
    valid = isfinite(double(x));
elseif isstring(x)
    valid = ~ismissing(x) & strlength(x) > 0;
elseif iscategorical(x)
    valid = ~isundefined(x);
elseif iscellstr(x)
    valid = ~cellfun(@isempty, x);
else
    error('Unsupported predictor type: %s', class(x));
end

valid = valid(:);

end

function [X, encoder] = encode_predictors(T, predictor_vars)

X = [];
encoded_names = strings(0, 1);

encoder = struct();
encoder.predictor_vars = predictor_vars(:);
encoder.entries = struct([]);

for i = 1:numel(predictor_vars)

    var_i = char(predictor_vars(i));
    x = T.(var_i);

    entry = struct();
    entry.name = string(var_i);
    entry.class = string(class(x));

    if isnumeric(x) || islogical(x)

        x_num = double(x(:));

        X = [X, x_num]; %#ok<AGROW>
        encoded_names(end + 1, 1) = string(var_i); %#ok<AGROW>

        entry.type = "numeric";
        entry.encoded_names = string(var_i);

    elseif isstring(x) || iscategorical(x) || iscellstr(x)

        x_str = string(x(:));
        cats = unique(x_str, 'stable');

        if numel(cats) <= 1

            entry.type = "categorical_constant";
            entry.categories = cats;
            entry.encoded_names = strings(0, 1);

        else

            % Reference category is cats(1). Create one-hot columns for the
            % remaining categories.
            encoded_i = zeros(numel(x_str), numel(cats) - 1);
            encoded_i_names = strings(numel(cats) - 1, 1);

            for c = 2:numel(cats)
                encoded_i(:, c - 1) = double(x_str == cats(c));
                encoded_i_names(c - 1) = string(var_i) + "_" + sanitize_name(cats(c));
            end

            X = [X, encoded_i]; %#ok<AGROW>
            encoded_names = [encoded_names; encoded_i_names]; %#ok<AGROW>

            entry.type = "categorical";
            entry.categories = cats;
            entry.reference_category = cats(1);
            entry.encoded_names = encoded_i_names;

        end

    else

        error('Unsupported predictor type for %s: %s', var_i, class(x));

    end

    encoder.entries = append_encoder_entry(encoder.entries, entry);

end

encoder.encoded_names = encoded_names;

if isempty(X)
    error('No usable predictor columns were generated.');
end

end

function entries = append_encoder_entry(entries, entry)

if isempty(entries)
    entries = entry;
    return;
end

all_fields = unique([string(fieldnames(entries)); string(fieldnames(entry))], 'stable');
entries = add_missing_entry_fields(entries, all_fields);
entry = add_missing_entry_fields(entry, all_fields);
entries = [entries; entry];

end

function S = add_missing_entry_fields(S, fields)

existing = string(fieldnames(S));
for i = 1:numel(fields)
    if ismember(fields(i), existing)
        continue;
    end
    [S.(char(fields(i)))] = deal([]);
end

end

function name = sanitize_name(x)

name = regexprep(string(x), '[^A-Za-z0-9_]+', '_');
name = matlab.lang.makeValidName(char(name));
name = string(name);

end

% =========================================================================
% Split
% =========================================================================

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

function is_train = resolve_train_mask(train_mask, valid, n_model)

train_mask = logical(train_mask(:));

if numel(train_mask) == numel(valid)
    is_train = train_mask(valid);
elseif numel(train_mask) == n_model
    is_train = train_mask;
else
    error(['TrainMask must have either height(T) elements or the number ', ...
           'of valid model rows.']);
end

if numel(is_train) ~= n_model
    error('Resolved TrainMask has the wrong number of rows.');
end

end

% =========================================================================
% Models
% =========================================================================

function model = train_one_model(model_type, Xtrain, ytrain, opts)

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

        model.k = opts.knn_k;
        model.design = "knn";

    case {"ensemble", "boosted_trees", "lsboost"}

        if exist('fitrensemble', 'file') ~= 2
            error(['fitrensemble was not found. ', ...
                   'This model requires the Statistics and Machine Learning Toolbox.']);
        end

        learner = templateTree( ...
            'MinLeafSize', opts.min_leaf_size);

        model.ensemble = fitrensemble( ...
            Xtrain, ...
            ytrain, ...
            'Method', 'LSBoost', ...
            'Learners', learner, ...
            'NumLearningCycles', opts.num_learning_cycles, ...
            'LearnRate', opts.learn_rate, ...
            'Options', statset('UseParallel', opts.use_parallel));

        model.design = "boosted_trees";
        model.num_learning_cycles = opts.num_learning_cycles;
        model.min_leaf_size = opts.min_leaf_size;
        model.learn_rate = opts.learn_rate;

    case {"bagged_trees", "bag"}

        if exist('fitrensemble', 'file') ~= 2
            error(['fitrensemble was not found. ', ...
                   'This model requires the Statistics and Machine Learning Toolbox.']);
        end

        learner = templateTree( ...
            'MinLeafSize', opts.min_leaf_size);

        model.ensemble = fitrensemble( ...
            Xtrain, ...
            ytrain, ...
            'Method', 'Bag', ...
            'Learners', learner, ...
            'NumLearningCycles', opts.num_learning_cycles, ...
            'Options', statset('UseParallel', opts.use_parallel));

        model.design = "bagged_trees";
        model.num_learning_cycles = opts.num_learning_cycles;
        model.min_leaf_size = opts.min_leaf_size;

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

    case {"ensemble", "boosted_trees", "lsboost", "bagged_trees", "bag"}

        y_pred = predict(model.ensemble, X);

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

% =========================================================================
% Output tables
% =========================================================================

function T_pred = build_prediction_table( ...
    T_model, ...
    predictor_vars, ...
    y_true, ...
    y_pred_raw, ...
    y_pred, ...
    is_train, ...
    model_name, ...
    model_type)

keep_candidates = [
    "condition_id"
    "condition_position"
    "condition_label"
    "step_idx"
    "realization_idx"
    "patch_idx"
    "Omega_sr"
    "Omega"
    "omega_sr"
    "omega_mean"
    "aperture_value"
    "SIM_WaveModel"
    "SIM_f0"
    "SIM_cs_bg"
    "REQ_M"
    "win_size"
    "M_eff_true_diag"];

vars = string(T_model.Properties.VariableNames);

keep_vars = keep_candidates(ismember(keep_candidates, vars));
keep_vars = unique([keep_vars(:); predictor_vars(:)], 'stable');

T_pred = T_model(:, cellstr(keep_vars));

T_pred.model_name = repmat(string(model_name), height(T_model), 1);
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

function T_metrics = compute_prediction_metrics_table( ...
    y_true, ...
    y_pred, ...
    is_train, ...
    model_name, ...
    model_type)

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

    row.model_name = string(model_name);
    row.model_type = string(model_type);
    row.split = split_i;

    T_metrics = [T_metrics; row]; %#ok<AGROW>
end

T_metrics = movevars(T_metrics, {'model_name', 'model_type', 'split'}, 'Before', 1);

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
