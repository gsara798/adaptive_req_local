function [MODEL_DEPLOY, MODEL_INFO, file_path] = load_q_model_deployment(model_dir, varargin)
%LOAD_Q_MODEL_DEPLOYMENT Load one deployment q-model bundle.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.load_q_model_deployment';

addRequired(p, 'model_dir', @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelName', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'FeatureSet', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelType', "bagged_trees", @(x) ischar(x) || isstring(x));
addParameter(p, 'FilePath', "", @(x) ischar(x) || isstring(x));

parse(p, model_dir, varargin{:});

if strlength(string(p.Results.FilePath)) > 0
    file_path = string(p.Results.FilePath);
else
    file_name = sprintf('qmodel_%s__%s__%s.mat', ...
        sanitize_token(p.Results.ModelName), ...
        sanitize_token(p.Results.FeatureSet), ...
        sanitize_token(p.Results.ModelType));
    file_path = fullfile(string(model_dir), file_name);
end

if ~exist(file_path, 'file')
    error('Deployment q-model file not found:\n%s', file_path);
end

S = load(file_path, 'MODEL_DEPLOY', 'MODEL_INFO');
MODEL_DEPLOY = S.MODEL_DEPLOY;
MODEL_INFO = S.MODEL_INFO;

end

function token = sanitize_token(x)

token = regexprep(string(x), '[^A-Za-z0-9_]+', '_');
token = matlab.lang.makeValidName(char(token));

end
