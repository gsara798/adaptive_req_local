function entry = register_trained_model(varargin)
%REGISTER_TRAINED_MODEL Copy/register a trained model in outputs/model_registry.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.register_trained_model';

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'SourceModelFile', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelObject', [], @(x) true);
addParameter(p, 'ModelId', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'RegistrySubdir', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'TestName', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'AnalysisLevel', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelName', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'FeatureSet', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelType', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'ModelRole', "operational", @(x) ischar(x) || isstring(x));
addParameter(p, 'TrainingDataset', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Target', "q_theory", @(x) ischar(x) || isstring(x));
addParameter(p, 'PredictorSummary', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'MetricsFile', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'SplitType', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'PerformanceSummary', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Notes', "", @(x) ischar(x) || isstring(x));

parse(p, varargin{:});
R = p.Results;

root_dir = char(R.RootDir);
registry_root = fullfile(root_dir, 'outputs', 'model_registry');
if strlength(string(R.RegistrySubdir)) > 0
    target_dir = fullfile(registry_root, char(R.RegistrySubdir));
else
    target_dir = registry_root;
end
make_dir_if_needed(registry_root);
make_dir_if_needed(target_dir);

model_id = string(R.ModelId);
if strlength(model_id) == 0
    model_id = matlab.lang.makeValidName(char(strjoin([ ...
        string(R.TestName), string(R.ModelName), string(R.FeatureSet), ...
        string(R.ModelType)], "__")));
end

target_model_file = "";
source_file = string(R.SourceModelFile);
if strlength(source_file) > 0 && exist(source_file, 'file') == 2
    [~, name, ext] = fileparts(source_file);
    target_model_file = string(fullfile(target_dir, name + ext));
    if ~strcmp(char(source_file), char(target_model_file))
        copyfile(source_file, target_model_file);
    end
elseif ~isempty(R.ModelObject)
    target_model_file = string(fullfile(target_dir, model_id + ".mat"));
    MODEL = R.ModelObject;
    save(target_model_file, 'MODEL', '-v7.3');
elseif strlength(source_file) > 0
    warning('Model file not found, skipping copy: %s', source_file);
end

entry = table();
entry.model_id = model_id;
entry.test_name = string(R.TestName);
entry.analysis_level = string(R.AnalysisLevel);
entry.model_name = string(R.ModelName);
entry.feature_set = string(R.FeatureSet);
entry.model_type = string(R.ModelType);
entry.model_role = string(R.ModelRole);
entry.training_dataset = string(R.TrainingDataset);
entry.target = string(R.Target);
entry.predictor_summary = string(R.PredictorSummary);
entry.model_file = target_model_file;
entry.metrics_file = string(R.MetricsFile);
entry.split_type = string(R.SplitType);
entry.performance_summary = string(R.PerformanceSummary);
entry.created_datetime = string(datetime('now'));
entry.notes = string(R.Notes);

manifest_file = fullfile(registry_root, 'model_manifest.csv');
if exist(manifest_file, 'file') == 2
    T_manifest = readtable(manifest_file, ...
        'TextType', 'string', ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');
    T_manifest = normalize_manifest_types(T_manifest);
else
    T_manifest = table();
end
entry = normalize_manifest_types(entry);

expected_vars = string(entry.Properties.VariableNames);
if ~isempty(T_manifest)
    manifest_vars = string(T_manifest.Properties.VariableNames);
    for i = 1:numel(expected_vars)
        if ~ismember(expected_vars(i), manifest_vars)
            T_manifest.(char(expected_vars(i))) = strings(height(T_manifest), 1);
        end
    end
    manifest_vars = string(T_manifest.Properties.VariableNames);
    for i = 1:numel(manifest_vars)
        if ~ismember(manifest_vars(i), expected_vars)
            entry.(char(manifest_vars(i))) = strings(height(entry), 1);
        end
    end
end

if ~isempty(T_manifest) && ismember('model_id', string(T_manifest.Properties.VariableNames))
    T_manifest(T_manifest.model_id == model_id, :) = [];
end
T_manifest = adaptive_req.analysis.Test12Analysis.concatTables(T_manifest, entry);
writetable(T_manifest, manifest_file);

end

function T = normalize_manifest_types(T)

for i = 1:width(T)
    name_i = T.Properties.VariableNames{i};
    if isstring(T.(name_i))
        continue;
    elseif iscellstr(T.(name_i)) || ischar(T.(name_i)) || iscategorical(T.(name_i))
        T.(name_i) = string(T.(name_i));
    elseif isdatetime(T.(name_i))
        T.(name_i) = string(T.(name_i));
    else
        T.(name_i) = string(T.(name_i));
    end
end

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end
