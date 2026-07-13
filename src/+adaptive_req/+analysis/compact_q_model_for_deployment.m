function MODEL_DEPLOY = compact_q_model_for_deployment(MODEL, varargin)
%COMPACT_Q_MODEL_FOR_DEPLOYMENT Keep only what is needed for q prediction.
%
% The full training output is useful for analysis, but deployment scripts
% should not need to load a multi-GB Level 01 model file. This helper keeps
% the encoder/scaling metadata and only the requested trained learners.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.compact_q_model_for_deployment';

addRequired(p, 'MODEL', @isstruct);
addParameter(p, 'ModelTypes', string.empty(0, 1), ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

parse(p, MODEL, varargin{:});

model_types = lower(string(p.Results.ModelTypes(:)));

MODEL_DEPLOY = struct();

copy_fields = [
    "model_name"
    "q_var"
    "input_predictor_vars"
    "resolved_predictor_vars"
    "encoded_predictor_names"
    "encoder"
    "mu"
    "sigma"
    "clip_predictions"
    "model_role"
    "fixed_split"
    "clip_range"];

for i = 1:numel(copy_fields)
    f = char(copy_fields(i));
    if isfield(MODEL, f)
        MODEL_DEPLOY.(f) = MODEL.(f);
    end
end

if ~isfield(MODEL_DEPLOY, 'clip_predictions')
    MODEL_DEPLOY.clip_predictions = true;
end

if ~isfield(MODEL, 'models') || isempty(MODEL.models)
    error('MODEL does not contain trained learners in MODEL.models.');
end

keep = true(numel(MODEL.models), 1);
if ~isempty(model_types)
    available = lower(string({MODEL.models.model_type}));
    keep = ismember(available, model_types);
end

models = MODEL.models(keep);
if isempty(models)
    error('No requested ModelTypes were found in MODEL.models.');
end

for i = 1:numel(models)
    models(i).model = compact_one_learner(models(i).model);
end

MODEL_DEPLOY.models = models;
MODEL_DEPLOY.deployment_compact = true;
MODEL_DEPLOY.created_at = string(datetime('now'));

end

function model = compact_one_learner(model)

if isstruct(model) && isfield(model, 'ensemble')
    try
        model.ensemble = compact(model.ensemble);
    catch
        % Some learner objects may already be compact or not support compact.
        % Prediction still works with the original object.
    end
end

end
