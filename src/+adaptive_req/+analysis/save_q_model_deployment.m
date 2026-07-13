function files = save_q_model_deployment(MODEL, out_dir, varargin)
%SAVE_Q_MODEL_DEPLOYMENT Save one q-model per deployable model type.
%
% Each output MAT contains a single variable MODEL_DEPLOY plus MODEL_INFO.
% This avoids loading the huge all-model Level 01 MAT file during map-level
% analyses.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.save_q_model_deployment';

addRequired(p, 'MODEL', @isstruct);
addRequired(p, 'out_dir', @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelName', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'FeatureSet', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelRole', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelTypes', string.empty(0, 1), ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'Overwrite', true, @(x) islogical(x) || isnumeric(x));

parse(p, MODEL, out_dir, varargin{:});

out_dir = string(out_dir);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

model_name = string(p.Results.ModelName);
if model_name == "" && isfield(MODEL, 'model_name')
    model_name = string(MODEL.model_name);
end

feature_set = string(p.Results.FeatureSet);
model_role = string(p.Results.ModelRole);
model_types = lower(string(p.Results.ModelTypes(:)));

if isempty(model_types)
    model_types = lower(string({MODEL.models.model_type}));
end

files = strings(0, 1);

for i = 1:numel(model_types)
    model_type = model_types(i);
    MODEL_DEPLOY = adaptive_req.analysis.compact_q_model_for_deployment( ...
        MODEL, 'ModelTypes', model_type);

    MODEL_INFO = struct();
    MODEL_INFO.model_name = model_name;
    MODEL_INFO.feature_set = feature_set;
    MODEL_INFO.model_role = model_role;
    MODEL_INFO.model_type = model_type;
    MODEL_INFO.created_at = string(datetime('now'));
    MODEL_INFO.file_format = "adaptive_req_q_model_deployment_v1";

    file_name = sprintf('qmodel_%s__%s__%s.mat', ...
        sanitize_token(model_name), sanitize_token(feature_set), ...
        sanitize_token(model_type));
    file_path = fullfile(out_dir, file_name);

    if exist(file_path, 'file') && ~logical(p.Results.Overwrite)
        files(end + 1, 1) = string(file_path); %#ok<AGROW>
        continue;
    end

    save(file_path, 'MODEL_DEPLOY', 'MODEL_INFO', '-v7.3');
    files(end + 1, 1) = string(file_path); %#ok<AGROW>
end

end

function token = sanitize_token(x)

token = regexprep(string(x), '[^A-Za-z0-9_]+', '_');
token = matlab.lang.makeValidName(char(token));

end
