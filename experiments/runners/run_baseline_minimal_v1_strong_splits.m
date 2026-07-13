%% run_baseline_minimal_v1_strong_splits.m
% Strong-split validation for baseline_minimal_v1.
%
% This runner does not recompute REQ and does not use any previous frozen
% q/confidence/correction models. It loads the feature dataset produced by
% experiments/runners/run_baseline_minimal_v1.m, builds grouped folds, and
% retrains the two deployable baselines inside every fold:
%
%   1. q_spectrum_only
%   2. q_spectrum_plus_composition
%
% For q_spectrum_plus_composition, the auxiliary composition models
% (predicted_patch_purity, p_mixed, p_strong_mixed) are also trained only on
% the train rows of the current fold. This is the main anti-leakage rule.
%
% Runtime controls:
%   ADAPTIVE_REQ_STRONG_SPLITS_MODE                  = quick | full
%   ADAPTIVE_REQ_STRONG_SPLITS_VALIDATE_ONLY         = true | false
%   ADAPTIVE_REQ_STRONG_SPLITS_CONFIG                = optional JSON path
%   ADAPTIVE_REQ_STRONG_SPLITS_OUTPUT_NAME           = optional output name
%   ADAPTIVE_REQ_STRONG_SPLITS_MAX_TRAIN_ROWS        = optional row cap
%   ADAPTIVE_REQ_STRONG_SPLITS_MAX_EVAL_ROWS         = optional row cap
%   ADAPTIVE_REQ_STRONG_SPLITS_USE_PARALLEL_TRAINING = true | false
%   ADAPTIVE_REQ_STRONG_SPLITS_SAVE_PATCH_PRED       = true | false
%   ADAPTIVE_REQ_STRONG_SPLITS_SAVE_FOLD_MODELS      = true | false

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'baseline_minimal_v1_strong_splits_configuration.json'));
log_path = fullfile(OUT.table_dir, 'runtime_log.txt');
log_fid = fopen(log_path, 'w');
cleanup_obj = onCleanup(@() fclose_if_open(log_fid)); %#ok<NASGU>

logmsg(log_fid, 'Baseline minimal v1 strong splits');
logmsg(log_fid, 'Mode: %s | validate only: %d', CFG.Mode, CFG.ValidateOnly);
logmsg(log_fid, 'Input dataset: %s', CFG.InputDatasetPath);
logmsg(log_fid, 'Output root: %s', OUT.root_dir);

if CFG.UseParallelTraining
    try
        if isempty(gcp('nocreate'))
            parpool('threads');
        end
    catch ME
        warning('StrongSplits:ParallelPool', 'Could not start parallel pool: %s', ME.message);
        CFG.UseParallelTraining = false;
    end
end

%% Load feature dataset

t_load = tic;
S = load(CFG.InputDatasetPath, 'T_train');
assert(isfield(S, 'T_train'), 'Dataset does not contain T_train: %s', CFG.InputDatasetPath);
T_all = S.T_train;
clear S;
logmsg(log_fid, 'Loaded %d patch rows in %.1f s.', height(T_all), toc(t_load));
T_all = normalize_dataset_table(T_all);
BASE_FEATURES = select_base_predictors(T_all);
assert_no_forbidden_predictors(BASE_FEATURES);
logmsg(log_fid, 'Leakage guard passed: %d base predictors.', numel(BASE_FEATURES));

FOLDS = build_folds(T_all, CFG);
assert(~isempty(FOLDS), 'No folds were created. Check configuration and dataset metadata.');
logmsg(log_fid, 'Candidate folds: %d', height(FOLDS));
writetable(FOLDS, fullfile(OUT.table_dir, 'fold_definitions.csv'));

ASSIGN = build_assignment_table(T_all, FOLDS);
writetable(ASSIGN, fullfile(OUT.table_dir, 'train_test_group_assignments.csv'));

if CFG.ValidateOnly
    logmsg(log_fid, 'Validation-only mode: checking the first fold with tiny samples.');
    CFG.MaxRowsForModelTraining = min(CFG.MaxRowsForModelTraining, 5000);
    CFG.MaxRowsForModelEvaluation = min(CFG.MaxRowsForModelEvaluation, 5000);
    FOLDS = FOLDS(1,:);
end

%% Run folds

summary_fold_parts = {};
summary_split_parts = {};
by_frequency_parts = {};
by_geometry_parts = {};
by_geometry_family_parts = {};
by_field_regime_parts = {};
by_field_family_parts = {};
by_M_parts = {};
by_purity_parts = {};
by_purity_group_parts = {};
by_roi_parts = {};
by_roi_frequency_parts = {};
by_geometry_frequency_parts = {};
q_scatter_parts = {};

for fi = 1:height(FOLDS)
    fold = FOLDS(fi,:);
    fold_id = fold.fold_id;
    fold_key = string(fold.fold_key);
    split_name = string(fold.split_name);
    t_fold = tic;
    logmsg(log_fid, '\n[%d/%d] %s | %s', fi, height(FOLDS), split_name, fold_key);

    [train_mask_full, test_mask_full] = masks_for_fold(T_all, fold);
    assert_no_group_leakage(T_all, train_mask_full, test_mask_full, fold);

    train_mask = sample_mask_by_condition(T_all, train_mask_full, CFG.MaxRowsForModelTraining, CFG.RandomSeed + 1000 + fi);
    test_mask = sample_mask_by_condition(T_all, test_mask_full, CFG.MaxRowsForModelEvaluation, CFG.RandomSeed + 2000 + fi);

    logmsg(log_fid, '  rows: train %d/%d | test %d/%d', ...
        sum(train_mask), sum(train_mask_full), sum(test_mask), sum(test_mask_full));

    T_train = T_all(train_mask,:);
    T_eval = T_all(test_mask,:);
    T_train.is_train_row = true(height(T_train),1);
    T_train.is_heldout_row = false(height(T_train),1);
    T_eval.is_train_row = false(height(T_eval),1);
    T_eval.is_heldout_row = true(height(T_eval),1);

    t_train = tic;
    MODELS = train_models_for_fold(T_train, BASE_FEATURES, CFG);
    train_seconds = toc(t_train);

    t_pred = tic;
    T_pred = apply_models_for_fold(T_eval, MODELS, BASE_FEATURES, fold, CFG);
    predict_seconds = toc(t_pred);

    if CFG.SavePatchPredictions
        pred_file = fullfile(OUT.patch_dir, sprintf('patch_predictions__fold_%03d__%s__%s.csv', ...
            fold_id, sanitize(split_name), sanitize(fold_key)));
        writetable(remove_cell_columns(T_pred), pred_file);
    end

    if CFG.SaveFoldModels
        save(fullfile(OUT.model_dir, sprintf('fold_%03d__%s__%s_models.mat', ...
            fold_id, sanitize(split_name), sanitize(fold_key))), ...
            'MODELS', 'BASE_FEATURES', 'fold', 'CFG', '-v7.3');
    end

    T_fold = summarize_predictions(T_pred, ["split_name","fold_key","model_name"]);
    T_fold.train_seconds = repmat(train_seconds, height(T_fold), 1);
    T_fold.predict_seconds = repmat(predict_seconds, height(T_fold), 1);
    T_fold.total_seconds = repmat(toc(t_fold), height(T_fold), 1);
    T_fold.n_train_rows_used = repmat(sum(train_mask), height(T_fold), 1);
    T_fold.n_eval_rows_used = repmat(sum(test_mask), height(T_fold), 1);
    summary_fold_parts{end+1,1} = T_fold; %#ok<AGROW>
    by_frequency_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","f0"]); %#ok<AGROW>
    by_geometry_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","case_id"]); %#ok<AGROW>
    by_geometry_family_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","geometry_family"]); %#ok<AGROW>
    by_field_regime_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","field_regime_ood"]); %#ok<AGROW>
    by_field_family_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","field_family"]); %#ok<AGROW>
    by_M_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","M"]); %#ok<AGROW>
    by_purity_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","purity_bin"]); %#ok<AGROW>
    by_purity_group_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","purity_group"]); %#ok<AGROW>
    by_roi_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","roi_label"]); %#ok<AGROW>
    by_roi_frequency_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","roi_label","f0"]); %#ok<AGROW>
    by_geometry_frequency_parts{end+1,1} = summarize_predictions(T_pred, ["split_name","fold_key","model_name","case_id","f0"]); %#ok<AGROW>
    q_scatter_parts{end+1,1} = sample_q_scatter(T_pred, CFG, fi); %#ok<AGROW>

    logmsg(log_fid, '  fold complete: train %.1f s | predict %.1f s | total %.1f s', ...
        train_seconds, predict_seconds, toc(t_fold));
end

%% Save summaries

T_by_fold = vertcat(summary_fold_parts{:});
T_by_split = summarize_fold_metrics(T_by_fold, ["split_name","model_name"]);
writetable(T_by_split, fullfile(OUT.table_dir, 'summary_metrics_by_split.csv'));
writetable(T_by_fold, fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'));
writetable(vertcat(by_frequency_parts{:}), fullfile(OUT.table_dir, 'metrics_by_frequency.csv'));
writetable(vertcat(by_geometry_parts{:}), fullfile(OUT.table_dir, 'metrics_by_geometry.csv'));
writetable(vertcat(by_geometry_family_parts{:}), fullfile(OUT.table_dir, 'metrics_by_geometry_family.csv'));
writetable(vertcat(by_field_regime_parts{:}), fullfile(OUT.table_dir, 'metrics_by_field_regime.csv'));
writetable(vertcat(by_field_family_parts{:}), fullfile(OUT.table_dir, 'metrics_by_field_family.csv'));
writetable(vertcat(by_M_parts{:}), fullfile(OUT.table_dir, 'metrics_by_M.csv'));
writetable(vertcat(by_purity_parts{:}), fullfile(OUT.table_dir, 'metrics_by_patch_purity_bin.csv'));
writetable(vertcat(by_purity_group_parts{:}), fullfile(OUT.table_dir, 'metrics_by_purity_group.csv'));
writetable(vertcat(by_roi_parts{:}), fullfile(OUT.table_dir, 'metrics_by_roi.csv'));
writetable(vertcat(by_roi_frequency_parts{:}), fullfile(OUT.table_dir, 'metrics_by_roi_frequency.csv'));
writetable(vertcat(by_geometry_frequency_parts{:}), fullfile(OUT.table_dir, 'metrics_by_geometry_frequency.csv'));
writetable(vertcat(q_scatter_parts{:}), fullfile(OUT.table_dir, 'q_scatter_sample.csv'));

RESULTS = struct('CFG', CFG, 'BASE_FEATURES', BASE_FEATURES, 'FOLDS', FOLDS, ...
    'summary_by_split', T_by_split, 'summary_by_fold', T_by_fold); %#ok<NASGU>
save(fullfile(OUT.data_dir, 'baseline_minimal_v1_strong_splits_results.mat'), ...
    'RESULTS', '-v7.3');

write_runner_readme(OUT, CFG, T_by_split, T_by_fold);
logmsg(log_fid, '\nStrong splits complete. Tables: %s', OUT.table_dir);
logmsg(log_fid, 'Run analysis next: run(''experiments/analysis/analyze_baseline_minimal_v1_strong_splits.m'')');

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_STRONG_SPLITS_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_STRONG_SPLITS_MODE must be quick or full.');

cfg_path = string(getenv('ADAPTIVE_REQ_STRONG_SPLITS_CONFIG'));
if cfg_path == ""
    cfg_path = fullfile(root_dir, 'configs', 'final_training', 'baseline_minimal_v1_strong_splits.json');
end
J = jsondecode(fileread(cfg_path));
CFG = struct();
CFG.Mode = mode;
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_STRONG_SPLITS_VALIDATE_ONLY', false);
CFG.ConfigPath = cfg_path;
CFG.InputDatasetPath = fullfile(root_dir, J.InputDataset);
CFG.OutputName = string_field(J, 'OutputName', 'baseline_minimal_v1_strong_splits');
out_override = lower(strtrim(string(getenv('ADAPTIVE_REQ_STRONG_SPLITS_OUTPUT_NAME'))));
if out_override ~= "", CFG.OutputName = out_override; end
CFG.TrainFraction = numeric_field(J, 'TrainFraction', 0.70);
CFG.RandomSeed = numeric_field(J, 'RandomSeed', 53100);
CFG.TreeLearners = numeric_field(J, 'TreeLearners', 180);
CFG.MinLeafSize = numeric_field(J, 'MinLeafSize', 8);
CFG.GroupedSeeds = double(J.GroupedSeeds(:))';
CFG.Splits = J.Splits;
CFG.QuickGeometryFamilies = string(J.QuickGeometryFamilies(:));
CFG.MaxRowsForModelTraining = numeric_field(J, 'MaxRowsForModelTraining', 200000);
CFG.MaxRowsForModelEvaluation = numeric_field(J, 'MaxRowsForModelEvaluation', 200000);
if mode == "quick"
    CFG.GroupedSeeds = double(J.QuickGroupedSeeds(:))';
    CFG.Splits = J.QuickSplits;
    CFG.MaxRowsForModelTraining = numeric_field(J, 'QuickMaxRowsForModelTraining', 50000);
    CFG.MaxRowsForModelEvaluation = numeric_field(J, 'QuickMaxRowsForModelEvaluation', 50000);
    CFG.TreeLearners = numeric_field(J, 'QuickTreeLearners', 80);
end
CFG.MaxRowsForModelTraining = env_number('ADAPTIVE_REQ_STRONG_SPLITS_MAX_TRAIN_ROWS', CFG.MaxRowsForModelTraining);
CFG.MaxRowsForModelEvaluation = env_number('ADAPTIVE_REQ_STRONG_SPLITS_MAX_EVAL_ROWS', CFG.MaxRowsForModelEvaluation);
CFG.UseParallelTraining = env_true('ADAPTIVE_REQ_STRONG_SPLITS_USE_PARALLEL_TRAINING', logical(J.UseParallelTraining));
CFG.SavePatchPredictions = env_true('ADAPTIVE_REQ_STRONG_SPLITS_SAVE_PATCH_PRED', logical(J.SavePatchPredictions));
CFG.SaveFoldModels = env_true('ADAPTIVE_REQ_STRONG_SPLITS_SAVE_FOLD_MODELS', logical(J.SaveFoldModels));
CFG.MaxQScatterRowsPerFoldModel = 5000;
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', char(CFG.OutputName));
if CFG.Mode == "quick"
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.model_dir = fullfile(OUT.root_dir, 'models');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.patch_dir = fullfile(OUT.table_dir, 'patch_predictions_by_fold');
for d = string(struct2cell(OUT))'
    if exist(d, 'dir') ~= 7, mkdir(d); end
end
end

%% Dataset and folds

function T = normalize_dataset_table(T)
vars = string(T.Properties.VariableNames);
if ~ismember('geometry_family', vars)
    if ismember('case_family', vars)
        T.geometry_family = string(T.case_family);
    else
        T.geometry_family = geometry_family_from_case(string(T.case_id));
    end
else
    T.geometry_family = string(T.geometry_family);
end
if ~ismember('field_family', string(T.Properties.VariableNames))
    T.field_family = field_family_from_regime(string(T.field_regime_ood));
end
if ~ismember('purity_bin', string(T.Properties.VariableNames))
    T.purity_bin = purity_bin(T.patch_purity);
else
    T.purity_bin = string(T.purity_bin);
end
if ~ismember('purity_group', string(T.Properties.VariableNames))
    T.purity_group = purity_group(T.patch_purity);
end
if ~ismember('is_mixed', string(T.Properties.VariableNames))
    T.is_mixed = T.patch_purity < 0.95;
end
if ~ismember('is_strong_mixed', string(T.Properties.VariableNames))
    T.is_strong_mixed = T.patch_purity < 0.75;
end
T.condition_key = string(T.condition_key);
T.case_id = string(T.case_id);
T.field_regime_ood = string(T.field_regime_ood);
T.geometry_family = string(T.geometry_family);
T.field_family = string(T.field_family);
end

function FOLDS = build_folds(T, CFG)
rows = {};
fold_id = 0;
if logical(CFG.Splits.GroupedConditionRepeated)
    groups = unique(string(T.condition_key), 'stable');
    for seed = CFG.GroupedSeeds
        rng(seed);
        order = groups(randperm(numel(groups)));
        ntrain = max(1, round(CFG.TrainFraction*numel(order)));
        fold_id = fold_id + 1;
        rows{end+1,1} = fold_row(fold_id, "grouped_condition_repeated", ...
            "seed_" + string(seed), "condition_key", "random_70_30", ...
            strjoin(order(1:ntrain), '|'), strjoin(order(ntrain+1:end), '|')); %#ok<AGROW>
    end
end
if logical(CFG.Splits.LeaveOneFrequencyOut)
    freqs = unique(T.f0(:))';
    for f = freqs
        fold_id = fold_id + 1;
        rows{end+1,1} = fold_row(fold_id, "leave_one_frequency_out", ...
            "f" + string(f), "f0", string(f), "", ""); %#ok<AGROW>
    end
end
if logical(CFG.Splits.LeaveOneGeometryFamilyOut)
    fams = unique(string(T.geometry_family), 'stable')';
    if CFG.Mode == "quick"
        fams = fams(ismember(fams, CFG.QuickGeometryFamilies));
    end
    for fam = fams
        fold_id = fold_id + 1;
        rows{end+1,1} = fold_row(fold_id, "leave_one_geometry_family_out", ...
            fam, "geometry_family", fam, "", ""); %#ok<AGROW>
    end
end
if logical(CFG.Splits.FieldRegimeWithinFamilyOOD)
    regimes = unique(string(T.field_regime_ood), 'stable');
    candidates = ["directional_2d_angle30", "partial_3d_16src", "diffuse_3d_seed2", "diffuse_3d_seed3"];
    for c = candidates
        if any(lower(regimes) == c)
            fold_id = fold_id + 1;
            rows{end+1,1} = fold_row(fold_id, "field_regime_within_family_ood", ...
                c, "field_regime_ood", c, "", ""); %#ok<AGROW>
        end
    end
end
if logical(CFG.Splits.LeaveFieldFamilyOut)
    fams = unique(string(T.field_family), 'stable')';
    for fam = fams
        fold_id = fold_id + 1;
        rows{end+1,1} = fold_row(fold_id, "leave_field_family_out", ...
            fam, "field_family", fam, "", ""); %#ok<AGROW>
    end
end
if isfield(CFG.Splits, 'LeaveOneMOut') && logical(CFG.Splits.LeaveOneMOut)
    mvals = unique(T.M(:))';
    for m = mvals
        fold_id = fold_id + 1;
        rows{end+1,1} = fold_row(fold_id, "leave_one_M_out", ...
            "M" + string(m), "M", string(m), "", ""); %#ok<AGROW>
    end
end
FOLDS = vertcat(rows{:});
end

function r = fold_row(fold_id, split_name, fold_key, group_var, test_value, train_groups, test_groups)
r = table(fold_id, string(split_name), string(fold_key), string(group_var), ...
    string(test_value), string(train_groups), string(test_groups), ...
    'VariableNames', {'fold_id','split_name','fold_key','group_var', ...
    'test_value','train_groups','test_groups'});
end

function [train_mask, test_mask] = masks_for_fold(T, fold)
split = string(fold.split_name);
gv = string(fold.group_var);
tv = string(fold.test_value);
switch split
    case "grouped_condition_repeated"
        train_groups = split_groups(fold.train_groups);
        test_groups = split_groups(fold.test_groups);
        train_mask = ismember(string(T.condition_key), train_groups);
        test_mask = ismember(string(T.condition_key), test_groups);
    case {"leave_one_frequency_out", "leave_one_M_out"}
        val = str2double(tv);
        test_mask = T.(gv) == val;
        train_mask = ~test_mask;
    otherwise
        test_mask = lower(string(T.(gv))) == lower(tv);
        train_mask = ~test_mask;
end
end

function G = split_groups(s)
if strlength(string(s)) == 0
    G = strings(0,1);
else
    G = string(strsplit(char(s), '|'))';
end
end

function ASSIGN = build_assignment_table(T, FOLDS)
conds = unique(string(T.condition_key), 'stable');
case_id = strings(numel(conds),1);
geom = strings(numel(conds),1);
regime = strings(numel(conds),1);
f0 = nan(numel(conds),1);
M = nan(numel(conds),1);
for i = 1:numel(conds)
    idx = find(string(T.condition_key)==conds(i),1);
    case_id(i) = string(T.case_id(idx));
    geom(i) = string(T.geometry_family(idx));
    regime(i) = string(T.field_regime_ood(idx));
    f0(i) = T.f0(idx);
    M(i) = T.M(idx);
end
base = table(conds, case_id, geom, regime, f0, M, ...
    'VariableNames', {'condition_key','case_id','geometry_family','field_regime_ood','f0','M'});
parts = cell(height(FOLDS),1);
for fi = 1:height(FOLDS)
    [tr, te] = masks_for_fold(T, FOLDS(fi,:));
    trc = unique(string(T.condition_key(tr)));
    tec = unique(string(T.condition_key(te)));
    X = base;
    X.fold_id = repmat(FOLDS.fold_id(fi), height(X), 1);
    X.split_name = repmat(string(FOLDS.split_name(fi)), height(X), 1);
    X.fold_key = repmat(string(FOLDS.fold_key(fi)), height(X), 1);
    X.is_train = ismember(X.condition_key, trc);
    X.is_test = ismember(X.condition_key, tec);
    parts{fi} = X;
end
ASSIGN = vertcat(parts{:});
end

function assert_no_group_leakage(T, train_mask, test_mask, fold)
assert(~any(train_mask & test_mask), 'Train/test row overlap in fold %d.', fold.fold_id);
train_conditions = unique(string(T.condition_key(train_mask)));
test_conditions = unique(string(T.condition_key(test_mask)));
assert(isempty(intersect(train_conditions, test_conditions)), ...
    'Condition leakage in fold %d.', fold.fold_id);
if string(fold.split_name) == "leave_one_frequency_out"
    assert(~any(T.f0(train_mask) == str2double(fold.test_value)), 'Frequency leakage.');
end
if string(fold.split_name) == "leave_one_geometry_family_out"
    assert(~any(lower(string(T.geometry_family(train_mask))) == lower(string(fold.test_value))), 'Geometry family leakage.');
end
if string(fold.split_name) == "leave_field_family_out"
    assert(~any(lower(string(T.field_family(train_mask))) == lower(string(fold.test_value))), 'Field family leakage.');
end
assert(any(train_mask) && any(test_mask), 'Empty train or test split in fold %d.', fold.fold_id);
end

%% Training and prediction

function MODELS = train_models_for_fold(T, base_features, CFG)
Fbase = T(:, cellstr(base_features));
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallelTraining);
train_mask = true(height(T),1);
valid_base = all(isfinite(table2array(Fbase)),2);
train_mask = train_mask & valid_base;
MIX = struct();
MIX.base_features = base_features;
MIX.purity = fitrensemble(Fbase(train_mask,:), T.patch_purity(train_mask), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, 'Options', opts);
MIX.mixed = fitcensemble(Fbase(train_mask,:), logical(T.is_mixed(train_mask)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, 'Options', opts);
MIX.strong_mixed = fitcensemble(Fbase(train_mask,:), logical(T.is_strong_mixed(train_mask)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, 'Options', opts);
P = predict_composition(MIX, T, base_features);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;
Q = struct();
Q.spectrum_only = fit_q_regressor(T_aug, train_mask, base_features, T_aug.q_oracle, CFG);
Q.spectrum_plus_composition = fit_q_regressor(T_aug, train_mask, ...
    [base_features; "predicted_patch_purity"; "p_mixed"; "p_strong_mixed"], T_aug.q_oracle, CFG);
MODELS = struct('composition', MIX, 'q', Q, ...
    'model_names', ["q_spectrum_only","q_spectrum_plus_composition"]);
end

function M = fit_q_regressor(T, train_mask, features, y, CFG)
X = T(:, cellstr(features));
valid = train_mask & all(isfinite(table2array(X)),2) & isfinite(y);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallelTraining);
M = struct();
M.features = features(:);
M.model = fitrensemble(X(valid,:), y(valid), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, 'Options', opts);
M.n_train = sum(valid);
end

function T_out = apply_models_for_fold(T, MODELS, base_features, fold, CFG) %#ok<INUSD>
P = predict_composition(MODELS.composition, T, base_features);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;
parts = cell(numel(MODELS.model_names),1);
keep = intersect(string(T_aug.Properties.VariableNames), ...
    ["dataset","condition_key","case_id","case_family","geometry_family", ...
    "field_regime","field_regime_ood","field_family","f0","M","dx","dz", ...
    "map_iz","map_ix","cx","cz","x_center_m","z_center_m", ...
    "true_SWS","k_true","patch_purity","purity_bin","purity_group", ...
    "is_mixed","is_strong_mixed","distance_to_boundary_mm","distance_bin", ...
    "q_oracle","q_theory_prior","sws_theory","is_train_row","is_heldout_row", ...
    "predicted_patch_purity","p_mixed","p_strong_mixed"], 'stable');
for mi = 1:numel(MODELS.model_names)
    name = MODELS.model_names(mi);
    R = T_aug(:, cellstr(keep));
    R.split_name = repmat(string(fold.split_name), height(R), 1);
    R.fold_key = repmat(string(fold.fold_key), height(R), 1);
    R.fold_id = repmat(fold.fold_id, height(R), 1);
    R.model_name = repmat(name, height(R), 1);
    switch name
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, T_aug);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T_aug);
        otherwise
            error('Unknown model: %s', name);
    end
    sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T_aug.q_oracle;
    R.sws_signed_error_pct = 100*(sws_pred - T_aug.true_SWS)./T_aug.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.sws_abs_error_ms = abs(sws_pred - T_aug.true_SWS);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    R.roi_label = assign_roi_labels(R);
    parts{mi} = R;
end
T_out = vertcat(parts{:});
end

function COMP = predict_composition(MIX, F, features)
X = F(:, cellstr(features));
COMP.predicted_patch_purity = min(max(predict(MIX.purity, X), 0), 1);
[~, score] = predict(MIX.mixed, X);
COMP.p_mixed = positive_score(MIX.mixed, score);
[~, score] = predict(MIX.strong_mixed, X);
COMP.p_strong_mixed = positive_score(MIX.strong_mixed, score);
end

function s = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
s = score(:,idx);
end

function y = predict_q_model(M, T)
X = T(:, cellstr(M.features));
y = clamp01(predict(M.model, X));
end

function y = q_to_sws(mappings, q, f0)
y = nan(numel(q),1);
for i = 1:numel(q)
    if isscalar(f0), fi = f0; else, fi = f0(i); end
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), fi);
end
end

%% Metrics

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T)
    S = table();
    return;
end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G == gi,:);
    rows{gi} = table(height(X), ...
        mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), ...
        mean(X.sws_abs_error_ms,'omitnan'), ...
        sqrt(mean((X.sws_pred - X.true_SWS).^2,'omitnan')), ...
        mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.sws_signed_error_pct > 0,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), ...
        mean(abs(X.q_error),'omitnan'), ...
        mean(X.q_error,'omitnan'), ...
        mean(X.true_SWS,'omitnan'), ...
        mean(X.sws_pred,'omitnan'), ...
        mean(X.patch_purity,'omitnan'), ...
        mean(X.predicted_patch_purity,'omitnan'), ...
        mean(X.p_mixed,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'MAE_sws','RMSE_sws','mean_signed_error_pct', ...
        'median_signed_error_pct','underestimate_pct','overestimate_pct', ...
        'high_error10_pct','high_error20_pct','mean_abs_q_error', ...
        'mean_q_error','mean_true_SWS','mean_pred_SWS', ...
        'mean_patch_purity','mean_predicted_patch_purity','mean_p_mixed'});
end
S = [groups vertcat(rows{:})];
end

function S = summarize_fold_metrics(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G == gi,:);
    rows{gi} = table(height(X), ...
        mean(X.MAPE_pct,'omitnan'), std(X.MAPE_pct,'omitnan'), ...
        mean(X.high_error20_pct,'omitnan'), std(X.high_error20_pct,'omitnan'), ...
        mean(X.mean_signed_error_pct,'omitnan'), std(X.mean_signed_error_pct,'omitnan'), ...
        mean(X.mean_abs_q_error,'omitnan'), std(X.mean_abs_q_error,'omitnan'), ...
        sum(X.N), mean(X.total_seconds,'omitnan'), ...
        'VariableNames', {'N_folds','MAPE_pct_mean','MAPE_pct_std', ...
        'high_error20_pct_mean','high_error20_pct_std', ...
        'bias_pct_mean','bias_pct_std','mean_abs_q_error_mean', ...
        'mean_abs_q_error_std','N_eval_rows_total','mean_fold_seconds'});
end
S = [groups vertcat(rows{:})];
end

function Tq = sample_q_scatter(T, CFG, fold_index)
rng(CFG.RandomSeed + 5000 + fold_index);
keep = ["split_name","fold_key","fold_id","model_name","case_id","geometry_family", ...
    "field_regime_ood","field_family","f0","M","q_oracle","q_pred", ...
    "q_error","sws_abs_error_pct","patch_purity","roi_label"];
keep = intersect(keep, string(T.Properties.VariableNames), 'stable');
Tq = T(:, cellstr(keep));
if height(Tq) > CFG.MaxQScatterRowsPerFoldModel
    [G,~] = findgroups(Tq.model_name);
    parts = {};
    for gi = unique(G)'
        idx = find(G == gi);
        n = min(numel(idx), ceil(CFG.MaxQScatterRowsPerFoldModel/numel(unique(G))));
        parts{end+1,1} = Tq(idx(randperm(numel(idx), n)),:); %#ok<AGROW>
    end
    Tq = vertcat(parts{:});
end
end

%% Predictors and labels

function features = select_base_predictors(T)
vars = string(T.Properties.VariableNames);
numeric = false(size(vars));
for i = 1:numel(vars)
    numeric(i) = isnumeric(T.(vars(i))) || islogical(T.(vars(i)));
end
forbidden_patterns = ["true","oracle","purity","mixed","confidence", ...
    "error","pred","sws","cs_","k_true","q_local","q_pred", ...
    "q_theory","req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition","distance"];
allowed = numeric;
for p = forbidden_patterns
    allowed = allowed & ~contains(lower(vars), lower(p));
end
keep_names = ["REQ_M","M","SIM_f0","f0","dx","dz","REQ_Nbins_effective"];
allowed = allowed | ismember(vars, keep_names);
features = vars(allowed);
features = features(:);
end

function assert_no_forbidden_predictors(features)
bad_patterns = ["true","oracle","purity","mixed","confidence","error", ...
    "pred","sws","cs_","k_true","q_local","q_pred","q_theory", ...
    "req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition","distance"];
features = lower(string(features));
for p = bad_patterns
    hit = features(contains(features, p));
    assert(isempty(hit), 'Forbidden predictor detected: %s', strjoin(hit, ', '));
end
end

function roi = assign_roi_labels(T)
roi = repmat("other", height(T), 1);
fam = string(T.geometry_family);
cs = T.true_SWS;
has_dist = ismember("distance_to_boundary_mm", string(T.Properties.VariableNames));
if has_dist, d = T.distance_to_boundary_mm; else, d = nan(height(T),1); end
is_hom = fam == "homogeneous";
if any(is_hom)
    x0 = median(T.x_center_m(is_hom), 'omitnan');
    z0 = median(T.z_center_m(is_hom), 'omitnan');
    center = is_hom & abs(T.x_center_m - x0) <= 4e-3 & abs(T.z_center_m - z0) <= 4e-3;
    roi(center) = "homogeneous_center";
    roi(is_hom & ~center) = "homogeneous_other";
end
is_hetero = ~is_hom;
if any(is_hetero)
    for ck = unique(T.condition_key(is_hetero), 'stable')'
        idx = T.condition_key == ck;
        vals = unique(round(cs(idx), 4));
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        soft = idx & abs(cs - min(vals)) < 1e-4;
        hard = idx & abs(cs - max(vals)) < 1e-4;
        mid = idx & ~soft & ~hard;
        roi(soft & d > 8) = "soft_core";
        roi(hard & d > 4) = "hard_core";
        roi(mid & d > 4) = "intermediate_core";
        roi(idx & d >= 0 & d <= 1) = "interface_0_1mm";
        roi(idx & d > 1 & d <= 2) = "interface_1_2mm";
        roi(idx & d > 2 & d <= 4) = "interface_2_4mm";
    end
end
end

function labels = purity_bin(p)
labels = repmat("[0.99,1.00]", size(p));
labels(p < 0.99) = "[0.95,0.99)";
labels(p < 0.95) = "[0.90,0.95)";
labels(p < 0.90) = "[0.75,0.90)";
labels(p < 0.75) = "[0.50,0.75)";
labels(p < 0.50) = "[0,0.50)";
labels(~isfinite(p)) = "unknown";
end

function labels = purity_group(p)
labels = repmat("pure_ge_0p95", size(p));
labels(p < 0.95) = "mixed_lt_0p95";
labels(p < 0.75) = "strongly_mixed_lt_0p75";
labels(~isfinite(p)) = "unknown";
end

function fam = geometry_family_from_case(case_id)
case_id = lower(string(case_id));
fam = repmat("complex", size(case_id));
fam(contains(case_id, "homogeneous")) = "homogeneous";
fam(contains(case_id, "bilayer") & ~contains(case_id, "smooth")) = "bilayer";
fam(contains(case_id, "inclusion") | contains(case_id, "ellipse")) = "inclusion";
fam(contains(case_id, "thin") | contains(case_id, "three_material") | contains(case_id, "smooth")) = "complex";
end

function fam = field_family_from_regime(regime)
r = lower(string(regime));
fam = repmat("other", size(r));
fam(contains(r, "directional")) = "directional_2D";
fam(contains(r, "diffuse_2d")) = "diffuse_2D";
fam(contains(r, "diffuse_3d")) = "diffuse_3D";
fam(contains(r, "partial_3d")) = "partial_3D";
end

%% Sampling and utilities

function mask = sample_mask_by_condition(T, base_mask, max_rows, seed)
idx = find(base_mask);
if max_rows <= 0 || numel(idx) <= max_rows
    mask = base_mask;
    return;
end
rng(seed);
conds = unique(string(T.condition_key(idx)), 'stable');
per_cond = max(1, floor(max_rows / max(1,numel(conds))));
keep = [];
for c = conds'
    ci = idx(string(T.condition_key(idx)) == c);
    n = min(numel(ci), per_cond);
    keep = [keep; ci(randperm(numel(ci), n))]; %#ok<AGROW>
end
if numel(keep) < max_rows
    rest = setdiff(idx, keep);
    nadd = min(numel(rest), max_rows - numel(keep));
    if nadd > 0, keep = [keep; rest(randperm(numel(rest), nadd))]; end %#ok<AGROW>
end
mask = false(height(T),1);
mask(keep) = true;
end

function q = clamp01(q)
q = min(max(q, 0), 1);
end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
keep = true(size(vars));
for i = 1:numel(vars)
    keep(i) = ~iscell(T.(vars(i)));
end
T = T(:, cellstr(vars(keep)));
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function logmsg(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
if fid > 0, fprintf(fid, '%s\n', msg); end
end

function fclose_if_open(fid)
if ~isempty(fid) && fid > 0, fclose(fid); end
end

function write_config_json(CFG, path)
try
    txt = jsonencode(CFG, 'PrettyPrint', true);
catch
    txt = jsonencode(CFG);
end
fid = fopen(path, 'w'); assert(fid > 0, 'Cannot write config: %s', path);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function write_runner_readme(OUT, CFG, T_by_split, T_by_fold)
path = fullfile(OUT.root_dir, 'README_results.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Baseline minimal v1 strong splits\n\n');
fprintf(fid, 'Generated by `experiments/runners/run_baseline_minimal_v1_strong_splits.m`.\n\n');
fprintf(fid, '## Scope\n\n');
fprintf(fid, '- No REQ recomputation.\n');
fprintf(fid, '- Models: `q_spectrum_only`, `q_spectrum_plus_composition`.\n');
fprintf(fid, '- Composition auxiliary models are retrained inside each fold from train rows only.\n\n');
fprintf(fid, '## Configuration\n\n');
fprintf(fid, '- Mode: `%s`\n', CFG.Mode);
fprintf(fid, '- Max train rows/fold: `%d`\n', CFG.MaxRowsForModelTraining);
fprintf(fid, '- Max eval rows/fold: `%d`\n', CFG.MaxRowsForModelEvaluation);
fprintf(fid, '- Tree learners: `%d`\n\n', CFG.TreeLearners);
fprintf(fid, '## Summary by split\n\n');
write_markdown_table(fid, T_by_split);
fprintf(fid, '\n## Summary by fold\n\n');
write_markdown_table(fid, T_by_fold(:, intersect(string(T_by_fold.Properties.VariableNames), ...
    ["split_name","fold_key","model_name","N","MAPE_pct","mean_signed_error_pct","high_error20_pct","mean_abs_q_error","total_seconds"], 'stable')));
fprintf(fid, '\n## Tables\n\nAll CSV tables are saved under `tables/`. Patch predictions, when enabled, are saved per fold under `tables/patch_predictions_by_fold/`.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T)
    fprintf(fid, '_No rows available._\n');
    return;
end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
max_rows = min(height(T), 80);
for i = 1:max_rows
    vals = strings(1, numel(names));
    for j = 1:numel(names)
        v = T.(names(j))(i);
        if isnumeric(v)
            vals(j) = sprintf('%.4g', v);
        elseif islogical(v)
            vals(j) = string(v);
        else
            vals(j) = string(v);
        end
    end
    fprintf(fid, '| %s |\n', strjoin(vals, ' | '));
end
if height(T) > max_rows
    fprintf(fid, '\n_Only first %d rows shown._\n', max_rows);
end
end

function tf = env_true(name, default_value)
v = lower(strtrim(string(getenv(name))));
if v == "", tf = default_value; return; end
tf = ismember(v, ["1","true","yes","y","on"]);
end

function x = env_number(name, default_value)
v = strtrim(string(getenv(name)));
if v == "", x = default_value; else, x = str2double(v); end
if ~isfinite(x), x = default_value; end
end

function x = numeric_field(S, name, default_value)
if isfield(S, name), x = double(S.(name)); else, x = default_value; end
end

function s = string_field(S, name, default_value)
if isfield(S, name), s = string(S.(name)); else, s = string(default_value); end
end
