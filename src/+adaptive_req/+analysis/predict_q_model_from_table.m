function T_pred = predict_q_model_from_table(MODEL, T, varargin)
%PREDICT_Q_MODEL_FROM_TABLE Apply a trained q model to a new feature table.
%
% MODEL is returned by adaptive_req.analysis.train_q_model_from_predictors.
% This function uses MODEL.encoder, MODEL.mu, and MODEL.sigma so deployment
% uses the exact same feature encoding as training.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.predict_q_model_from_table';

addRequired(p, 'MODEL', @isstruct);
addRequired(p, 'T', @istable);
addParameter(p, 'ModelType', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ClipPredictions', MODEL.clip_predictions, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ModelName', MODEL.model_name, ...
    @(x) ischar(x) || isstring(x));

parse(p, MODEL, T, varargin{:});

model_type = lower(string(p.Results.ModelType));
clip_predictions = logical(p.Results.ClipPredictions);
model_name = string(p.Results.ModelName);

if model_type == ""
    model_idx = 1;
else
    available = lower(string({MODEL.models.model_type}));
    model_idx = find(available == model_type, 1, 'first');

    if isempty(model_idx)
        error('ModelType not found in MODEL: %s', model_type);
    end
end

model_i = MODEL.models(model_idx);
X_raw = encode_table_with_model(T, MODEL.encoder);
X = (X_raw - MODEL.mu) ./ MODEL.sigma;

y_pred_raw = predict_one_model(model_i.model, model_i.model_type, X);

if clip_predictions
    y_pred = min(max(y_pred_raw, 0), 1);
else
    y_pred = y_pred_raw;
end

T_pred = keep_metadata(T);
T_pred.model_name = repmat(model_name, height(T), 1);
T_pred.model_type = repmat(string(model_i.model_type), height(T), 1);
T_pred.q_pred_raw = y_pred_raw;
T_pred.q_pred = y_pred;

if isfield(MODEL, 'q_var') && ismember(MODEL.q_var, string(T.Properties.VariableNames))
    T_pred.q_true = T.(char(MODEL.q_var));
    T_pred.residual = T_pred.q_pred - T_pred.q_true;
    T_pred.abs_error = abs(T_pred.residual);
end

end

function X = encode_table_with_model(T, encoder)

entries = encoder.entries;
X_parts = cell(numel(entries), 1);

for i = 1:numel(entries)
    entry = entries(i);
    name = char(entry.name);

    if ~ismember(name, T.Properties.VariableNames)
        error('Required predictor is missing from table: %s', name);
    end

    x = T.(name);

    switch string(entry.type)
        case "numeric"
            X_parts{i} = double(x(:));

        case "categorical_constant"
            X_parts{i} = zeros(height(T), 0);

        case "categorical"
            x_str = string(x(:));
            cats = string(entry.categories);
            encoded = zeros(height(T), max(0, numel(cats) - 1));

            for c = 2:numel(cats)
                encoded(:, c - 1) = double(x_str == cats(c));
            end

            X_parts{i} = encoded;

        otherwise
            error('Unsupported encoder entry type: %s', entry.type);
    end
end

X = cat(2, X_parts{:});

end

function y_pred = predict_one_model(model, model_type, X)

switch string(model_type)
    case {"linear", "quadratic"}
        Phi = build_design_matrix(X, model.design);
        y_pred = Phi * model.beta;

    case {"ensemble", "boosted_trees", "lsboost", "bagged_trees", "bag"}
        y_pred = predict(model.ensemble, X);

    otherwise
        error('Prediction is not implemented for model type: %s', model_type);
end

end

function Phi = build_design_matrix(X, design)

n = size(X, 1);
p = size(X, 2);

switch string(design)
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

function T_meta = keep_metadata(T)

candidate_vars = [
    "condition_id"
    "condition_position"
    "condition_label"
    "step_idx"
    "realization_idx"
    "patch_idx"
    "patch_label"
    "SIM_WaveModel"
    "SIM_f0"
    "SIM_cs_bg"
    "REQ_M"
    "Omega_sr"
    "omega_mean"
    "aperture_value"
    "cx"
    "cz"
    "x_center_m"
    "z_center_m"
    "cs_true_patch_mean"
    "cs_true_patch_median"
    "roi_name"];

vars = string(T.Properties.VariableNames);
keep_vars = candidate_vars(ismember(candidate_vars, vars));

if isempty(keep_vars)
    T_meta = table();
else
    T_meta = T(:, cellstr(keep_vars));
end

end
