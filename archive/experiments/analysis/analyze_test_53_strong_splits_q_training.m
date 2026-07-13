%% analyze_test_53_strong_splits_q_training.m
% Test 53 strong-split validation for the clean q/SWS baseline.
%
% Purpose
%   Re-train the Test 53 clean q models under stricter train/test splits:
%   repeated grouped condition split, leave-one-frequency-out,
%   leave-one-geometry-family-out, field-regime OOD, and optional
%   leave-one-M-out. Composition predictors are re-trained inside every fold
%   using only that fold's training rows.
%
% Runtime controls
%   ADAPTIVE_REQ_TEST53_STRONG_MODE              = quick | full
%   ADAPTIVE_REQ_TEST53_STRONG_ESTIMATE_ONLY     = true | false
%   ADAPTIVE_REQ_TEST53_STRONG_VALIDATE_ONLY     = true | false
%   ADAPTIVE_REQ_TEST53_STRONG_USE_PARALLEL      = true | false
%   ADAPTIVE_REQ_TEST53_STRONG_SAVE_MODELS       = true | false
%   ADAPTIVE_REQ_TEST53_STRONG_SAVE_PREDICTIONS  = true | false
%   ADAPTIVE_REQ_TEST53_STRONG_MAX_TRAIN_ROWS    = integer
%   ADAPTIVE_REQ_TEST53_STRONG_MAX_EVAL_ROWS     = integer
%   ADAPTIVE_REQ_TEST53_STRONG_BENCH_ROWS        = integer
%   ADAPTIVE_REQ_TEST53_STRONG_MAX_MAPS          = integer
%   ADAPTIVE_REQ_TEST53_STRONG_MAP_MIN_COVERAGE  = fraction, diagnostic only
%   ADAPTIVE_REQ_TEST53_STRONG_MAX_FOLDS         = integer, optional debug cap
%
% Notes
%   - true_SWS, q_oracle, true patch_purity, errors, distances, coordinates,
%     and confidence variables are never used as predictors.
%   - true_SWS/q_oracle are used only as labels/evaluation targets.
%   - true patch_purity/distance are used only for stratified diagnostics.
%   - q_spectrum_plus_composition uses composition models trained only on
%     the fold training rows.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
diary_file = fullfile(OUT.root_dir, 'runtime_log.txt');
if exist(diary_file, 'file') == 2
    delete(diary_file);
end
diary(diary_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\nTest 53 strong-split q baseline validation\n');
fprintf('Mode: %s | estimate only: %d | validate only: %d\n', ...
    CFG.Mode, CFG.EstimateOnly, CFG.ValidateOnly);
fprintf('Dataset: %s\n', CFG.DatasetFile);
fprintf('No full-run is started unless ESTIMATE_ONLY=false.\n');

if CFG.ValidateOnly
    validate_script(CFG, OUT);
    fprintf('Validation-only checks passed.\n');
    return;
end

timer_total = tic;

fprintf('\nLoading Test 53 feature dataset. This can take a while for full.\n');
t_load = tic;
S = load(CFG.DatasetFile, 'T_train', 'CFG_saved');
T = S.T_train;
fprintf('Loaded %d patch rows in %.1f s.\n', height(T), toc(t_load));

T = prepare_dataset(T);
BASE_FEATURES = select_base_predictors(T);
assert_no_forbidden_predictors(BASE_FEATURES);
write_config_json(CFG, OUT, BASE_FEATURES);

FOLDS = build_fold_table(T, CFG);
if isfinite(CFG.MaxFolds) && CFG.MaxFolds > 0 && height(FOLDS) > CFG.MaxFolds
    fprintf('Debug fold cap active: retaining first %d/%d folds.\n', CFG.MaxFolds, height(FOLDS));
    FOLDS = FOLDS(1:CFG.MaxFolds,:);
end
fprintf('Candidate folds: %d\n', height(FOLDS));
disp(groupsummary(FOLDS, "split_name"));

write_fold_definitions(FOLDS, OUT);
write_group_assignments(T, FOLDS, OUT);

BENCH = run_mini_benchmark(T, FOLDS, BASE_FEATURES, CFG, OUT);
write_benchmark_readme(BENCH, FOLDS, T, CFG, OUT);

if CFG.EstimateOnly
    fprintf('\nEstimate-only requested. Not running all folds.\n');
    fprintf('Runtime estimate written to:\n%s\n', fullfile(OUT.root_dir, 'README_results.md'));
    return;
end

all_pred_parts = {};
all_summary_parts = {};
all_roi_parts = {};
fold_runtime = table();
maps_saved = 0;

fprintf('\nRunning selected folds...\n');
for fi = 1:height(FOLDS)
    fold = FOLDS(fi,:);
    fprintf('\n[%d/%d] %s | %s\n', fi, height(FOLDS), ...
        fold.split_name, fold.fold_name);
    t_fold = tic;
    [T_pred, T_summary, T_roi, RUNTIME] = run_fold(T, fold, BASE_FEATURES, CFG, OUT);
    fold_runtime = [fold_runtime; RUNTIME]; %#ok<AGROW>
    all_summary_parts{end+1,1} = T_summary; %#ok<AGROW>
    all_roi_parts{end+1,1} = T_roi; %#ok<AGROW>
    if CFG.MaxExampleMaps > maps_saved
        maps_saved = maps_saved + plot_representative_maps( ...
            T_pred, fold, CFG, OUT, CFG.MaxExampleMaps - maps_saved);
    end

    if CFG.KeepPredictionsInMemory
        all_pred_parts{end+1,1} = T_pred; %#ok<AGROW>
    end

    if CFG.SavePredictions
        pred_file = fullfile(OUT.pred_dir, "predictions__" + fold.fold_key + ".mat");
        save(pred_file, 'T_pred', 'fold', '-v7.3');
    end

    fprintf('Fold complete in %.1f s.\n', toc(t_fold));
end

T_summary_all = vertcat(all_summary_parts{:});
T_roi_all = vertcat_nonempty(all_roi_parts);
if CFG.KeepPredictionsInMemory && ~isempty(all_pred_parts)
    T_predictions_all = vertcat(all_pred_parts{:}); %#ok<NASGU>
else
    T_predictions_all = table(); %#ok<NASGU>
end

write_all_outputs(T_summary_all, T_roi_all, fold_runtime, FOLDS, CFG, OUT);
plot_all_outputs(T_summary_all, T_roi_all, CFG, OUT);
write_results_readme(T_summary_all, T_roi_all, fold_runtime, BENCH, FOLDS, CFG, OUT);

save(fullfile(OUT.root_dir, 'test53_strong_splits_results.mat'), ...
    'T_summary_all', 'T_roi_all', 'fold_runtime', 'FOLDS', 'CFG', ...
    'BASE_FEATURES', 'BENCH', '-v7.3');

fprintf('\nTest 53 strong-splits complete in %.1f min.\n', toc(timer_total)/60);
fprintf('Tables: %s\nFigures: %s\nREADME: %s\n', ...
    OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST53_STRONG_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), ...
    'ADAPTIVE_REQ_TEST53_STRONG_MODE must be quick or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickRun = mode == "quick";
CFG.FullRun = mode == "full";
CFG.EstimateOnly = env_true('ADAPTIVE_REQ_TEST53_STRONG_ESTIMATE_ONLY', true);
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST53_STRONG_VALIDATE_ONLY', false);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST53_STRONG_USE_PARALLEL', true);
CFG.SaveModels = env_true('ADAPTIVE_REQ_TEST53_STRONG_SAVE_MODELS', false);
CFG.SavePredictions = env_true('ADAPTIVE_REQ_TEST53_STRONG_SAVE_PREDICTIONS', false);
CFG.KeepPredictionsInMemory = env_true('ADAPTIVE_REQ_TEST53_STRONG_KEEP_PREDICTIONS', false);
CFG.RandomSeed = 53053;
CFG.TrainFraction = 0.70;
CFG.TreeLearners = 180;
CFG.MinLeafSize = 8;
CFG.BenchmarkRows = env_number('ADAPTIVE_REQ_TEST53_STRONG_BENCH_ROWS', 50000);
CFG.MaxFolds = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_FOLDS', Inf);
CFG.OutputName = "test53_strong_splits";
CFG.DatasetFile = fullfile(root_dir, 'outputs', ...
    'test_53_paper_final_clean_q_training', 'data', ...
    'test38_velocity_field_diverse_dataset.mat');

if CFG.QuickRun
    CFG.GroupedSeeds = [53053 53054];
    CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_TRAIN_ROWS', 200000);
    CFG.MaxEvalRows = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_EVAL_ROWS', 200000);
    CFG.GeometryFamilies = ["bilayer","inclusion"];
    CFG.IncludeFieldWithinOOD = true;
    CFG.IncludeFieldFamilyOOD = true;
    CFG.IncludeLeaveMOut = false;
    CFG.MaxExampleMaps = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_MAPS', 8);
else
    CFG.GroupedSeeds = [53053 53054 53055 53056 53057];
    CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_TRAIN_ROWS', 500000);
    CFG.MaxEvalRows = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_EVAL_ROWS', 500000);
    CFG.GeometryFamilies = ["homogeneous","bilayer","inclusion","complex"];
    CFG.IncludeFieldWithinOOD = true;
    CFG.IncludeFieldFamilyOOD = true;
    CFG.IncludeLeaveMOut = env_true('ADAPTIVE_REQ_TEST53_STRONG_INCLUDE_M_OUT', true);
    CFG.MaxExampleMaps = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAX_MAPS', 20);
end

CFG.Frequencies = [200 300 400 500 600];
CFG.ModelNames = ["theory_discrete","fixed_q_train_median", ...
    "q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"];
CFG.PurityEdges = [0 0.50 0.75 0.90 0.95 0.99 1.0000001];
CFG.PurityLabels = ["0-0.50","0.50-0.75","0.75-0.90", ...
    "0.90-0.95","0.95-0.99","0.99-1.00"];
CFG.ROI.SideMm = 6;
CFG.ROI.HomogeneousSideMm = 8;
CFG.ROI.SoftDistanceMm = 8;
CFG.ROI.HardDistanceMm = 4;
CFG.ROI.InterfaceBandMm = 1;
CFG.ROI.MinPatches = 4;
CFG.MapMinCoverage = env_number('ADAPTIVE_REQ_TEST53_STRONG_MAP_MIN_COVERAGE', 0.35);
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', char(CFG.OutputName));
if CFG.QuickRun
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
else
    OUT.root_dir = fullfile(OUT.root_dir, 'full');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.model_dir = fullfile(OUT.root_dir, 'models');
OUT.pred_dir = fullfile(OUT.root_dir, 'predictions');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_representative');
for d = string(struct2cell(OUT))'
    if exist(d, 'dir') ~= 7
        mkdir(d);
    end
end
end

function write_config_json(CFG, OUT, BASE_FEATURES)
cfg = CFG;
cfg.BaseFeatures = BASE_FEATURES;
txt = jsonencode(cfg, 'PrettyPrint', true);
fid = fopen(fullfile(OUT.root_dir, 'test53_strong_splits_config.json'), 'w');
fprintf(fid, '%s\n', txt);
fclose(fid);
end

%% Data preparation and fold definitions

function T = prepare_dataset(T)
T.paper_geometry_family = map_string_values(T.case_id, @paper_geometry_family);
T.paper_field_family = map_string_values(T.field_regime_ood, @paper_field_family);
T.condition_group = T.condition_key;
T.patch_purity_bin_strong = purity_bin_strong(T.patch_purity);
T.purity_group = purity_group(T.patch_purity);
if ~ismember("sws_theory", string(T.Properties.VariableNames))
    T.sws_theory = q_to_sws(T.req_mapping, T.q_theory_prior, T.f0);
end

function out = map_string_values(values, mapper)
values = string(values);
u = unique(values, 'stable');
mapped = strings(size(u));
for i = 1:numel(u)
    mapped(i) = mapper(u(i));
end
out = strings(size(values));
for i = 1:numel(u)
    out(values == u(i)) = mapped(i);
end
end
end

function FOLDS = build_fold_table(T, CFG)
rows = {};
row = 0;
base_groups = unique(T(:, {'condition_group','case_id','paper_geometry_family', ...
    'field_regime_ood','paper_field_family','f0','M'}), 'rows', 'stable');

for si = 1:numel(CFG.GroupedSeeds)
    seed = CFG.GroupedSeeds(si);
    [train_groups, test_groups] = random_group_split(base_groups.condition_group, ...
        CFG.TrainFraction, seed);
    row = row + 1;
    rows{row,1} = make_fold("grouped_condition_repeated", ...
        "seed_" + string(seed), "grouped_seed_" + string(seed), ...
        train_groups, test_groups, seed, "condition_group");
end

for f = CFG.Frequencies
    train_groups = base_groups.condition_group(base_groups.f0 ~= f);
    test_groups = base_groups.condition_group(base_groups.f0 == f);
    row = row + 1;
    rows{row,1} = make_fold("leave_one_frequency_out", ...
        "test_" + string(f) + "Hz", "lofo_" + string(f), ...
        train_groups, test_groups, CFG.RandomSeed, "f0");
end

for fam = CFG.GeometryFamilies
    train_groups = base_groups.condition_group(base_groups.paper_geometry_family ~= fam);
    test_groups = base_groups.condition_group(base_groups.paper_geometry_family == fam);
    row = row + 1;
    rows{row,1} = make_fold("leave_one_geometry_family_out", ...
        "test_" + fam, "logeo_" + fam, ...
        train_groups, test_groups, CFG.RandomSeed, "paper_geometry_family");
end

if CFG.IncludeFieldWithinOOD
    within = ["directional_2D_angle30","diffuse_3D_seed3","partial_3D_16src"];
    for r = within
        train_groups = base_groups.condition_group(base_groups.field_regime_ood ~= r);
        test_groups = base_groups.condition_group(base_groups.field_regime_ood == r);
        row = row + 1;
        rows{row,1} = make_fold("field_regime_within_family_ood", ...
            "test_" + r, "field_within_" + sanitize_token(r), ...
            train_groups, test_groups, CFG.RandomSeed, "field_regime_ood");
    end
end

if CFG.IncludeFieldFamilyOOD
    families = ["directional_2D","diffuse_2D","diffuse_3D","partial_3D"];
    for fam = families
        train_groups = base_groups.condition_group(base_groups.paper_field_family ~= fam);
        test_groups = base_groups.condition_group(base_groups.paper_field_family == fam);
        row = row + 1;
        rows{row,1} = make_fold("leave_field_family_out", ...
            "test_" + fam, "field_family_" + sanitize_token(fam), ...
            train_groups, test_groups, CFG.RandomSeed, "paper_field_family");
    end
end

if CFG.IncludeLeaveMOut
    for m = unique(base_groups.M)'
        train_groups = base_groups.condition_group(base_groups.M ~= m);
        test_groups = base_groups.condition_group(base_groups.M == m);
        row = row + 1;
        rows{row,1} = make_fold("leave_one_M_out", ...
            "test_M" + string(m), "looM_" + string(m), ...
            train_groups, test_groups, CFG.RandomSeed, "M");
    end
end

FOLDS = struct2table(vertcat(rows{:}));
FOLDS.n_train_groups = cellfun(@numel, FOLDS.train_groups);
FOLDS.n_test_groups = cellfun(@numel, FOLDS.test_groups);
end

function F = make_fold(split_name, fold_name, fold_key, train_groups, test_groups, seed, group_variable)
F = struct();
F.split_name = string(split_name);
F.fold_name = string(fold_name);
F.fold_key = string(fold_key);
F.train_groups = {unique(string(train_groups), 'stable')};
F.test_groups = {unique(string(test_groups), 'stable')};
F.seed = seed;
F.group_variable = string(group_variable);
end

function [train_groups, test_groups] = random_group_split(groups, train_fraction, seed)
rng(seed, 'twister');
groups = unique(string(groups), 'stable');
groups = groups(randperm(numel(groups)));
ntrain = round(train_fraction * numel(groups));
ntrain = min(max(ntrain, 1), numel(groups)-1);
train_groups = groups(1:ntrain);
test_groups = groups((ntrain+1):end);
end

function write_fold_definitions(FOLDS, OUT)
F = FOLDS;
F.train_groups = cellfun(@(x) strjoin(string(x), '|'), F.train_groups, 'UniformOutput', false);
F.test_groups = cellfun(@(x) strjoin(string(x), '|'), F.test_groups, 'UniformOutput', false);
writetable(F, fullfile(OUT.table_dir, 'fold_definitions.csv'));
end

function write_group_assignments(T, FOLDS, OUT)
base = unique(T(:, {'condition_group','condition_key','case_id','paper_geometry_family', ...
    'field_regime_ood','paper_field_family','f0','M'}), 'rows', 'stable');
parts = cell(height(FOLDS),1);
for i = 1:height(FOLDS)
    F = base;
    F.split_name = repmat(FOLDS.split_name(i), height(F), 1);
    F.fold_name = repmat(FOLDS.fold_name(i), height(F), 1);
    F.fold_key = repmat(FOLDS.fold_key(i), height(F), 1);
    F.assignment = repmat("unused", height(F), 1);
    F.assignment(ismember(F.condition_group, FOLDS.train_groups{i})) = "train";
    F.assignment(ismember(F.condition_group, FOLDS.test_groups{i})) = "test";
    parts{i} = F;
end
writetable(vertcat(parts{:}), fullfile(OUT.table_dir, 'train_test_group_assignments.csv'));
end

%% Fold execution

function [T_pred, T_summary, T_roi, RUNTIME] = run_fold(T, fold, BASE_FEATURES, CFG, OUT)
t0 = tic;
train_mask = ismember(T.condition_group, fold.train_groups{1});
test_mask = ismember(T.condition_group, fold.test_groups{1});
validate_fold_masks(T, train_mask, test_mask, fold);

train_fit_mask = sample_mask_by_condition(T, train_mask, CFG.MaxTrainRows, fold.seed + 11);
test_eval_mask = sample_mask_by_condition(T, test_mask, CFG.MaxEvalRows, fold.seed + 29);

T_train = T(train_fit_mask,:);
T_eval = T(test_eval_mask,:);
fprintf('Rows: train %d/%d | test %d/%d\n', ...
    height(T_train), sum(train_mask), height(T_eval), sum(test_mask));

t_train = tic;
MODELS = train_fold_models(T_train, BASE_FEATURES, CFG);
train_seconds = toc(t_train);

t_pred = tic;
T_pred = apply_fold_models(T_eval, MODELS, BASE_FEATURES, fold);
predict_seconds = toc(t_pred);

t_metrics = tic;
T_summary = make_all_summaries(T_pred, fold);
T_roi = compute_roi_metrics(T_pred, fold, CFG);
metric_seconds = toc(t_metrics);

if CFG.SaveModels
    model_file = fullfile(OUT.model_dir, "models__" + fold.fold_key + ".mat");
    save(model_file, 'MODELS', 'BASE_FEATURES', 'fold', 'CFG', '-v7.3');
end

RUNTIME = table(fold.split_name, fold.fold_name, fold.fold_key, ...
    height(T_train), height(T_eval), train_seconds, predict_seconds, ...
    metric_seconds, toc(t0), ...
    'VariableNames', {'split_name','fold_name','fold_key', ...
    'n_train_rows','n_eval_rows','train_seconds','predict_seconds', ...
    'metric_seconds','total_seconds'});
end

function MODELS = train_fold_models(T, BASE_FEATURES, CFG)
Fbase = T(:, cellstr(BASE_FEATURES));
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);

comp = struct();
comp.base_features = BASE_FEATURES;
comp.purity = fitrensemble(Fbase, T.patch_purity, ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
comp.mixed = fitcensemble(Fbase, logical(T.patch_purity < 0.95), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
comp.strong_mixed = fitcensemble(Fbase, logical(T.patch_purity < 0.75), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);

P = predict_composition(comp, T, BASE_FEATURES);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;

q = struct();
q.spectrum_only = fit_q_regressor(T_aug, BASE_FEATURES, T_aug.q_oracle, CFG);
q.spectrum_plus_composition = fit_q_regressor(T_aug, ...
    [BASE_FEATURES; "predicted_patch_purity"; "p_mixed"; "p_strong_mixed"], ...
    T_aug.q_oracle, CFG);
q.spectrum_plus_theory_composition = fit_q_regressor(T_aug, ...
    [BASE_FEATURES; "q_theory_prior"; "predicted_patch_purity"; ...
    "p_mixed"; "p_strong_mixed"], T_aug.q_oracle, CFG);
q.fixed_q_train_median = median(T_aug.q_oracle, 'omitnan');

MODELS = struct();
MODELS.composition = comp;
MODELS.q = q;
MODELS.model_names = CFG.ModelNames;
end

function M = fit_q_regressor(T, features, y, CFG)
valid = all(isfinite(table2array(T(:, cellstr(features)))), 2) & isfinite(y);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
M = struct();
M.features = features(:);
M.n_train = sum(valid);
M.model = fitrensemble(T(valid, cellstr(features)), y(valid), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
end

function T_pred = apply_fold_models(T, MODELS, BASE_FEATURES, fold)
P = predict_composition(MODELS.composition, T, BASE_FEATURES);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;

parts = cell(numel(MODELS.model_names), 1);
for mi = 1:numel(MODELS.model_names)
    model_name = MODELS.model_names(mi);
    R = keep_eval_columns(T_aug);
    R.split_name = repmat(fold.split_name, height(R), 1);
    R.fold_name = repmat(fold.fold_name, height(R), 1);
    R.fold_key = repmat(fold.fold_key, height(R), 1);
    R.model_name = repmat(model_name, height(R), 1);

    switch model_name
        case "theory_discrete"
            q_pred = T_aug.q_theory_prior;
            sws_pred = T_aug.sws_theory;
        case "fixed_q_train_median"
            q_pred = repmat(MODELS.q.fixed_q_train_median, height(T_aug), 1);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_only"
            q_pred = clamp01(predict(MODELS.q.spectrum_only.model, ...
                T_aug(:, cellstr(MODELS.q.spectrum_only.features))));
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_composition"
            q_pred = clamp01(predict(MODELS.q.spectrum_plus_composition.model, ...
                T_aug(:, cellstr(MODELS.q.spectrum_plus_composition.features))));
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = clamp01(predict(MODELS.q.spectrum_plus_theory_composition.model, ...
                T_aug(:, cellstr(MODELS.q.spectrum_plus_theory_composition.features))));
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        otherwise
            error('Unknown model name: %s', model_name);
    end

    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T_aug.q_oracle;
    R.sws_error = sws_pred - T_aug.true_SWS;
    R.sws_signed_error_pct = 100 * R.sws_error ./ T_aug.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{mi} = R;
end
T_pred = vertcat(parts{:});
end

function R = keep_eval_columns(T)
vars = ["condition_key","condition_group","case_id","case_family", ...
    "paper_geometry_family","field_regime","field_regime_ood", ...
    "paper_field_family","f0","M","dx","dz","map_iz","map_ix", ...
    "cx","cz","x_center_m","z_center_m","true_SWS","k_true", ...
    "patch_purity","patch_purity_bin_strong","purity_group", ...
    "distance_to_boundary_mm","distance_bin","q_oracle","q_theory_prior", ...
    "sws_theory","is_mixed","is_strong_mixed", ...
    "predicted_patch_purity","p_mixed","p_strong_mixed"];
vars = intersect(vars, string(T.Properties.VariableNames), 'stable');
R = T(:, cellstr(vars));
end

function validate_fold_masks(T, train_mask, test_mask, fold)
assert(any(train_mask), 'Fold %s has empty train set.', fold.fold_key);
assert(any(test_mask), 'Fold %s has empty test set.', fold.fold_key);
train_conditions = unique(T.condition_key(train_mask));
test_conditions = unique(T.condition_key(test_mask));
assert(isempty(intersect(train_conditions, test_conditions)), ...
    'Condition leakage detected in fold %s.', fold.fold_key);

switch string(fold.split_name)
    case "leave_one_frequency_out"
        assert(isempty(intersect(unique(T.f0(train_mask)), unique(T.f0(test_mask)))), ...
            'Frequency leakage detected in %s.', fold.fold_key);
    case "leave_one_geometry_family_out"
        assert(isempty(intersect(unique(T.paper_geometry_family(train_mask)), ...
            unique(T.paper_geometry_family(test_mask)))), ...
            'Geometry-family leakage detected in %s.', fold.fold_key);
    case "leave_field_family_out"
        assert(isempty(intersect(unique(T.paper_field_family(train_mask)), ...
            unique(T.paper_field_family(test_mask)))), ...
            'Field-family leakage detected in %s.', fold.fold_key);
    case "leave_one_M_out"
        assert(isempty(intersect(unique(T.M(train_mask)), unique(T.M(test_mask)))), ...
            'M leakage detected in %s.', fold.fold_key);
end
end

%% Metrics

function T_summary = make_all_summaries(T_pred, fold)
parts = {};
parts{end+1,1} = summarize_predictions(T_pred, "overall", "all", "all");
parts{end+1,1} = summarize_predictions(T_pred, "frequency", "f0", "all");
parts{end+1,1} = summarize_predictions(T_pred, "M", "M", "all");
parts{end+1,1} = summarize_predictions(T_pred, "geometry", "case_id", "all");
parts{end+1,1} = summarize_predictions(T_pred, "geometry_family", "paper_geometry_family", "all");
parts{end+1,1} = summarize_predictions(T_pred, "field_regime", "field_regime_ood", "all");
parts{end+1,1} = summarize_predictions(T_pred, "field_family", "paper_field_family", "all");
parts{end+1,1} = summarize_predictions(T_pred, "patch_purity_bin", "patch_purity_bin_strong", "all");
parts{end+1,1} = summarize_predictions(T_pred, "purity_group", "purity_group", "all");
T_summary = vertcat(parts{:});
T_summary.split_name = repmat(fold.split_name, height(T_summary), 1);
T_summary.fold_name = repmat(fold.fold_name, height(T_summary), 1);
T_summary.fold_key = repmat(fold.fold_key, height(T_summary), 1);
T_summary = movevars(T_summary, ["split_name","fold_name","fold_key"], 'Before', 1);
end

function S = summarize_predictions(T, group_type, group_var, group_value_name)
if group_type == "overall" || group_var == "overall" || group_var == "all"
    G = ones(height(T),1);
    groups = table("overall", "all", 'VariableNames', {'group_type','group_value'});
else
    group_vals = T.(group_var);
    [G, group_raw] = findgroups(group_vals);
    groups = table(repmat(string(group_type), numel(group_raw), 1), ...
        string(group_raw), 'VariableNames', {'group_type','group_value'});
end

[GM, model_raw] = findgroups(T.model_name);
combo = findgroups(G, GM);
rows = cell(max(combo), 1);
for ci = 1:max(combo)
    idx = combo == ci;
    gi = mode(G(idx));
    mi = mode(GM(idx));
    X = T(idx,:);
    rows{ci} = metric_row(X, groups.group_type(gi), groups.group_value(gi), ...
        model_raw(mi));
end
S = vertcat(rows{:});
if group_value_name ~= "all"
    S.group_value_name = repmat(string(group_value_name), height(S), 1);
else
    S.group_value_name = repmat(string(group_var), height(S), 1);
end
S = movevars(S, 'group_value_name', 'After', 'group_value');
end

function R = metric_row(X, group_type, group_value, model_name)
err = X.sws_error;
err_pct = X.sws_signed_error_pct;
abs_pct = X.sws_abs_error_pct;
qerr = X.q_error;
R = table(string(group_type), string(group_value), string(model_name), ...
    height(X), ...
    mean(abs_pct, 'omitnan'), ...
    median(abs_pct, 'omitnan'), ...
    mean(abs(err), 'omitnan'), ...
    sqrt(mean(err.^2, 'omitnan')), ...
    mean(err_pct, 'omitnan'), ...
    median(err_pct, 'omitnan'), ...
    100*mean(err_pct < 0, 'omitnan'), ...
    100*mean(X.high_error10, 'omitnan'), ...
    100*mean(X.high_error20, 'omitnan'), ...
    mean(abs(qerr), 'omitnan'), ...
    mean(qerr, 'omitnan'), ...
    mean(X.predicted_patch_purity, 'omitnan'), ...
    mean(X.p_mixed, 'omitnan'), ...
    'VariableNames', {'group_type','group_value','model_name','N', ...
    'MAPE_pct','median_abs_error_pct','MAE_sws','RMSE_sws', ...
    'mean_bias_pct','median_bias_pct','underestimate_pct', ...
    'high_error10_pct','high_error20_pct','MAE_q','bias_q', ...
    'mean_predicted_patch_purity','mean_p_mixed'});
end

function T_roi = compute_roi_metrics(T_pred, fold, CFG)
keys = unique(T_pred(:, {'condition_key','case_id','paper_geometry_family', ...
    'field_regime_ood','paper_field_family','f0','M'}), 'rows', 'stable');
parts = {};
for ki = 1:height(keys)
    idx_condition = T_pred.condition_key == keys.condition_key(ki);
    Tcond = T_pred(idx_condition,:);
    base = Tcond(Tcond.model_name == Tcond.model_name(1),:);
    ROI = define_rois_for_condition(base, CFG);
    for ri = 1:numel(ROI)
        mask = ROI(ri).mask;
        if sum(mask) < CFG.ROI.MinPatches
            continue;
        end
        for model_name = unique(Tcond.model_name)'
            idx = Tcond.model_name == model_name;
            Xm = Tcond(idx,:);
            X = Xm(mask,:);
            if isempty(X), continue; end
            R = metric_row(X, "roi", ROI(ri).roi_type, model_name);
            R.roi_type = repmat(ROI(ri).roi_type, height(R), 1);
            R.roi_side_mm = repmat(ROI(ri).side_mm, height(R), 1);
            R.roi_center_x_mm = repmat(ROI(ri).center_x_mm, height(R), 1);
            R.roi_center_z_mm = repmat(ROI(ri).center_z_mm, height(R), 1);
            R.true_SWS_mean = mean(X.true_SWS, 'omitnan');
            R.pred_SWS_mean = mean(X.sws_pred, 'omitnan');
            R.pred_SWS_std = std(X.sws_pred, 'omitnan');
            R.mean_patch_purity = mean(X.patch_purity, 'omitnan');
            meta = keys(ki,:);
            R.condition_key = repmat(meta.condition_key, height(R), 1);
            R.case_id = repmat(meta.case_id, height(R), 1);
            R.paper_geometry_family = repmat(meta.paper_geometry_family, height(R), 1);
            R.field_regime_ood = repmat(meta.field_regime_ood, height(R), 1);
            R.paper_field_family = repmat(meta.paper_field_family, height(R), 1);
            R.f0 = repmat(meta.f0, height(R), 1);
            R.M = repmat(meta.M, height(R), 1);
            parts{end+1,1} = R; %#ok<AGROW>
        end
    end
end
if isempty(parts)
    T_roi = table();
else
    T_roi = vertcat(parts{:});
    T_roi.split_name = repmat(fold.split_name, height(T_roi), 1);
    T_roi.fold_name = repmat(fold.fold_name, height(T_roi), 1);
    T_roi.fold_key = repmat(fold.fold_key, height(T_roi), 1);
    T_roi = movevars(T_roi, ["split_name","fold_name","fold_key"], 'Before', 1);
end
end

function ROI = define_rois_for_condition(T, CFG)
ROI = struct('roi_type',{},'mask',{},'side_mm',{},'center_x_mm',{},'center_z_mm',{});
xmm = T.x_center_m * 1e3;
zmm = T.z_center_m * 1e3;
cx0 = median(xmm, 'omitnan');
cz0 = median(zmm, 'omitnan');
family = string(T.paper_geometry_family(1));
true_vals = unique(round(T.true_SWS, 4));

if family == "homogeneous" || numel(true_vals) == 1
    side = CFG.ROI.HomogeneousSideMm;
    mask = square_mask(xmm, zmm, cx0, cz0, side);
    ROI(end+1) = make_roi("homogeneous_center", mask, side, cx0, cz0); %#ok<AGROW>
    return;
end

soft = min(true_vals);
hard = max(true_vals);
dist = T.distance_to_boundary_mm;
side = CFG.ROI.SideMm;

soft_candidates = abs(T.true_SWS - soft) < 1e-3 & dist >= CFG.ROI.SoftDistanceMm;
if any(soft_candidates)
    [cx, cz] = robust_roi_center(xmm, zmm, soft_candidates);
    mask = square_mask(xmm, zmm, cx, cz, side) & soft_candidates;
    ROI(end+1) = make_roi("soft_core", mask, side, cx, cz); %#ok<AGROW>
end

hard_candidates = abs(T.true_SWS - hard) < 1e-3 & dist >= CFG.ROI.HardDistanceMm;
if any(hard_candidates)
    [cx, cz] = robust_roi_center(xmm, zmm, hard_candidates);
    mask = square_mask(xmm, zmm, cx, cz, side) & hard_candidates;
    ROI(end+1) = make_roi("hard_core", mask, side, cx, cz); %#ok<AGROW>
end

if numel(true_vals) > 2
    mid_vals = true_vals(true_vals ~= soft & true_vals ~= hard);
    for mv = reshape(mid_vals, 1, [])
        mid_candidates = abs(T.true_SWS - mv) < 1e-3 & dist >= CFG.ROI.HardDistanceMm;
        if any(mid_candidates)
            [cx, cz] = robust_roi_center(xmm, zmm, mid_candidates);
            mask = square_mask(xmm, zmm, cx, cz, side) & mid_candidates;
            ROI(end+1) = make_roi("mid_core", mask, side, cx, cz); %#ok<AGROW>
        end
    end
end

interface_mask = dist <= CFG.ROI.InterfaceBandMm;
if any(interface_mask)
    ROI(end+1) = make_roi("interface_0_1mm", interface_mask, NaN, NaN, NaN); %#ok<AGROW>
end
end

function R = make_roi(type, mask, side_mm, cx, cz)
R = struct('roi_type', string(type), 'mask', logical(mask), ...
    'side_mm', side_mm, 'center_x_mm', cx, 'center_z_mm', cz);
end

function mask = square_mask(xmm, zmm, cx, cz, side_mm)
mask = abs(xmm - cx) <= side_mm/2 & abs(zmm - cz) <= side_mm/2;
end

function [cx, cz] = robust_roi_center(xmm, zmm, idx)
cx = median(xmm(idx), 'omitnan');
cz = median(zmm(idx), 'omitnan');
end

%% Benchmark and outputs

function BENCH = run_mini_benchmark(T, FOLDS, BASE_FEATURES, CFG, OUT)
fprintf('\nMini-benchmark for runtime estimate.\n');
fold = FOLDS(1,:);
train_mask = ismember(T.condition_group, fold.train_groups{1});
test_mask = ismember(T.condition_group, fold.test_groups{1});
train_fit_mask = sample_mask_by_condition(T, train_mask, CFG.BenchmarkRows, fold.seed + 101);
test_eval_mask = sample_mask_by_condition(T, test_mask, CFG.BenchmarkRows, fold.seed + 202);
T_train = T(train_fit_mask,:);
T_eval = T(test_eval_mask,:);

t_train = tic;
MODELS = train_fold_models(T_train, BASE_FEATURES, CFG);
train_seconds = toc(t_train);

t_pred = tic;
T_pred = apply_fold_models(T_eval, MODELS, BASE_FEATURES, fold);
predict_seconds = toc(t_pred);

t_met = tic;
T_summary = make_all_summaries(T_pred, fold); %#ok<NASGU>
metric_seconds = toc(t_met);

BENCH = struct();
BENCH.n_train_rows = height(T_train);
BENCH.n_eval_rows = height(T_eval);
BENCH.train_seconds = train_seconds;
BENCH.predict_seconds = predict_seconds;
BENCH.metric_seconds = metric_seconds;
BENCH.total_seconds = train_seconds + predict_seconds + metric_seconds;
BENCH.seconds_per_train_row = train_seconds / max(1, height(T_train));
BENCH.seconds_per_eval_row = (predict_seconds + metric_seconds) / max(1, height(T_eval));
BENCH.fold_count = height(FOLDS);
BENCH.estimated_full_seconds = height(FOLDS) * ( ...
    BENCH.seconds_per_train_row * CFG.MaxTrainRows + ...
    BENCH.seconds_per_eval_row * CFG.MaxEvalRows);

writetable(struct2table(BENCH), fullfile(OUT.table_dir, 'mini_benchmark.csv'));
fprintf('Benchmark rows: train %d, eval %d\n', BENCH.n_train_rows, BENCH.n_eval_rows);
fprintf('Benchmark time: train %.1f s, predict %.1f s, metrics %.1f s\n', ...
    train_seconds, predict_seconds, metric_seconds);
fprintf('Estimated selected-run time: %.1f min for %d folds.\n', ...
    BENCH.estimated_full_seconds/60, height(FOLDS));
end

function write_benchmark_readme(BENCH, FOLDS, T, CFG, OUT)
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w');
fprintf(fid, '# Test 53 Strong-Split Validation\n\n');
fprintf(fid, 'This file was generated after the mini-benchmark. ');
fprintf(fid, 'If `EstimateOnly=true`, no full fold analysis was run.\n\n');
fprintf(fid, '## Dataset\n\n');
fprintf(fid, '- Patch rows loaded: %d\n', height(T));
fprintf(fid, '- Conditions: %d\n', numel(unique(T.condition_key)));
fprintf(fid, '- Cases: %d\n', numel(unique(T.case_id)));
fprintf(fid, '- Frequencies: %s Hz\n', mat2str(unique(T.f0)'));
fprintf(fid, '- M values: %s\n', mat2str(unique(T.M)'));
fprintf(fid, '- Folds selected in this mode: %d\n\n', height(FOLDS));
fprintf(fid, '## Mini-Benchmark\n\n');
fprintf(fid, '- Train rows: %d\n', BENCH.n_train_rows);
fprintf(fid, '- Eval rows: %d\n', BENCH.n_eval_rows);
fprintf(fid, '- Train time: %.1f s\n', BENCH.train_seconds);
fprintf(fid, '- Prediction time: %.1f s\n', BENCH.predict_seconds);
fprintf(fid, '- Metric time: %.1f s\n', BENCH.metric_seconds);
fprintf(fid, '- Estimated selected-run time: %.1f min\n\n', BENCH.estimated_full_seconds/60);
fprintf(fid, '## Recommendation\n\n');
if BENCH.estimated_full_seconds/3600 > 8
    fprintf(fid, 'The selected run is expensive. Run Quick first, then launch Full overnight.\n');
else
    fprintf(fid, 'The selected run appears feasible on this machine.\n');
end
fclose(fid);
end

function write_all_outputs(T_summary, T_roi, fold_runtime, FOLDS, CFG, OUT)
writetable(T_summary, fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'));
writetable(fold_runtime, fullfile(OUT.table_dir, 'fold_runtime.csv'));

S_overall = T_summary(T_summary.group_type == "overall", :);
writetable(aggregate_over_folds(S_overall, ["split_name","model_name"]), ...
    fullfile(OUT.table_dir, 'summary_metrics_by_split.csv'));
writetable(select_group(T_summary, "frequency"), fullfile(OUT.table_dir, 'metrics_by_frequency.csv'));
writetable(select_group(T_summary, "geometry"), fullfile(OUT.table_dir, 'metrics_by_geometry.csv'));
writetable(select_group(T_summary, "geometry_family"), fullfile(OUT.table_dir, 'metrics_by_geometry_family.csv'));
writetable(select_group(T_summary, "field_regime"), fullfile(OUT.table_dir, 'metrics_by_field_regime.csv'));
writetable(select_group(T_summary, "field_family"), fullfile(OUT.table_dir, 'metrics_by_field_family.csv'));
writetable(select_group(T_summary, "M"), fullfile(OUT.table_dir, 'metrics_by_M.csv'));
writetable(select_group(T_summary, "patch_purity_bin"), fullfile(OUT.table_dir, 'metrics_by_patch_purity_bin.csv'));
writetable(select_group(T_summary, "purity_group"), fullfile(OUT.table_dir, 'metrics_by_purity_group.csv'));
if ~isempty(T_roi)
    writetable(T_roi, fullfile(OUT.table_dir, 'metrics_by_roi.csv'));
else
    writetable(table(), fullfile(OUT.table_dir, 'metrics_by_roi.csv'));
end

% Compact compatibility aliases requested by the prompt.
copyfile(fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'), ...
    fullfile(OUT.table_dir, 'summary_metrics_by_fold_all_groups.csv'));
end

function G = select_group(T, name)
G = T(T.group_type == string(name), :);
end

function A = aggregate_over_folds(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
metrics = ["MAPE_pct","median_abs_error_pct","MAE_sws","RMSE_sws", ...
    "mean_bias_pct","high_error10_pct","high_error20_pct","MAE_q","bias_q"];
for gi = 1:max(G)
    X = T(G == gi,:);
    R = table(height(X), 'VariableNames', {'N_folds'});
    for m = metrics
        R.(m + "_mean") = mean(X.(m), 'omitnan');
        R.(m + "_std") = std(X.(m), 'omitnan');
    end
    rows{gi} = R;
end
A = [groups vertcat(rows{:})];
end

function plot_all_outputs(T_summary, T_roi, CFG, OUT)
safe_plot(@() plot_split_box(T_summary, OUT, "MAPE_pct", ...
    "MAPE by model and split", "MAPE (%)", 'test53_strong_mape_by_split.png'));
safe_plot(@() plot_split_box(T_summary, OUT, "high_error20_pct", ...
    "High-error >20% by model and split", "Pixels with error >20% (%)", ...
    'test53_strong_high20_by_split.png'));
safe_plot(@() plot_fold_lines(T_summary, OUT, "leave_one_frequency_out", ...
    'test53_strong_leave_frequency_mape.png'));
safe_plot(@() plot_fold_lines(T_summary, OUT, "leave_one_geometry_family_out", ...
    'test53_strong_leave_geometry_family_mape.png'));
safe_plot(@() plot_fold_lines(T_summary, OUT, "leave_field_family_out", ...
    'test53_strong_leave_field_family_mape.png'));
safe_plot(@() plot_purity(T_summary, OUT, "MAPE_pct", ...
    "MAPE versus patch-purity bin", "MAPE (%)", 'test53_strong_mape_by_purity.png'));
safe_plot(@() plot_purity(T_summary, OUT, "high_error20_pct", ...
    "High-error >20% versus patch-purity bin", "Pixels with error >20% (%)", ...
    'test53_strong_high20_by_purity.png'));
if ~isempty(T_roi)
    safe_plot(@() plot_roi_bars(T_roi, OUT), 'ROI plots');
end
end

function n_saved = plot_representative_maps(T_pred, fold, CFG, OUT, max_to_save)
n_saved = 0;
if max_to_save <= 0 || isempty(T_pred)
    return;
end
main_model = "q_spectrum_plus_composition";
Tmain = T_pred(T_pred.model_name == main_model, :);
if isempty(Tmain)
    return;
end

[G, cond] = findgroups(Tmain.condition_key);
mape = splitapply(@(x) mean(x, 'omitnan'), Tmain.sws_abs_error_pct, G);
fam = splitapply(@(x) string(x(1)), Tmain.paper_geometry_family, G);
case_id = splitapply(@(x) string(x(1)), Tmain.case_id, G);
regime = splitapply(@(x) string(x(1)), Tmain.field_regime_ood, G);
freq = splitapply(@(x) x(1), Tmain.f0, G);
mval = splitapply(@(x) x(1), Tmain.M, G);
C = table(cond, fam, case_id, regime, freq, mval, mape, ...
    'VariableNames', {'condition_key','family','case_id','regime','f0','M','MAPE_pct'});

selected = strings(0,1);
for f = unique(C.family, 'stable')'
    X = C(C.family == f,:);
    if isempty(X), continue; end
    [~, iworst] = max(X.MAPE_pct);
    selected(end+1,1) = X.condition_key(iworst); %#ok<AGROW>
end
[~, ibest] = min(C.MAPE_pct);
[~, iworst] = max(C.MAPE_pct);
selected = unique([selected; C.condition_key(ibest); C.condition_key(iworst)], 'stable');
selected = selected(1:min(numel(selected), max_to_save));

for ci = 1:numel(selected)
    key = selected(ci);
    X = T_pred(T_pred.condition_key == key, :);
    Xm = X(X.model_name == main_model, :);
    if isempty(Xm), continue; end
    try
        save_condition_map(Xm, fold, CFG, OUT);
        n_saved = n_saved + 1;
    catch ME
        warning('Representative map failed for %s: %s', key, ME.message);
    end
end
end

function save_condition_map(T, fold, CFG, OUT)
[Ztrue, xg, zg] = table_to_grid(T, 'true_SWS');
[Zpred, ~, ~] = table_to_grid(T, 'sws_pred');
[Zabs, ~, ~] = table_to_grid(T, 'sws_abs_error_pct');
[Zsigned, ~, ~] = table_to_grid(T, 'sws_signed_error_pct');
[Zq, ~, ~] = table_to_grid(T, 'q_pred');
[Zqo, ~, ~] = table_to_grid(T, 'q_oracle');
[Zpur, ~, ~] = table_to_grid(T, 'patch_purity');
[Zppur, ~, ~] = table_to_grid(T, 'predicted_patch_purity');
[Zpmix, ~, ~] = table_to_grid(T, 'p_mixed');

fig = figure('Color','w','Position',[80 80 1450 980]);
tl = tiledlayout(fig, 3, 3, 'Padding','compact', 'TileSpacing','compact');
map_coverage = nnz(isfinite(Zpred)) / max(numel(Zpred), 1);
coverage_note = sprintf('map coverage %.0f%%', 100*map_coverage);
if map_coverage < CFG.MapMinCoverage
    coverage_note = sprintf('%s, sparse diagnostic', coverage_note);
end
title(tl, sprintf('%s | %s | %s | f=%g Hz | M=%g | %s | %s', ...
    char(T.case_id(1)), char(T.field_regime_ood(1)), char(fold.fold_name), ...
    T.f0(1), T.M(1), char(pretty_model(T.model_name(1))), coverage_note), ...
    'Interpreter','none');

sws_lim = finite_prctile([Ztrue(:); Zpred(:)], [1 99]);
err_lim = [0 finite_prctile(Zabs(:), 99)];
signed_lim = finite_symmetric_limit(Zsigned(:), 98);

plot_map_panel(xg, zg, Ztrue, 'True SWS', 'SWS (m/s)', sws_lim);
hold on; overlay_rois(T, CFG); hold off;
plot_map_panel(xg, zg, Zpred, 'Predicted SWS', 'SWS (m/s)', sws_lim);
plot_map_panel(xg, zg, Zabs, 'Absolute error', 'Error (%)', err_lim);
plot_map_panel(xg, zg, Zsigned, 'Signed error', 'Error (%)', signed_lim);
plot_map_panel(xg, zg, Zq, 'Predicted q', 'REQ quantile q', [0 1]);
plot_map_panel(xg, zg, Zqo, 'Oracle q', 'REQ quantile q', [0 1]);
plot_map_panel(xg, zg, Zpur, 'True patch purity', 'Fraction', [0 1]);
plot_map_panel(xg, zg, Zppur, 'Fold-predicted patch purity', 'Fraction', [0 1]);
plot_map_panel(xg, zg, Zpmix, 'Fold-predicted mixedness', 'Probability', [0 1]);

folder = fullfile(OUT.map_dir, char(fold.split_name), char(sanitize_token(fold.fold_key)), ...
    char(sanitize_token(T.paper_geometry_family(1))));
if exist(folder, 'dir') ~= 7, mkdir(folder); end
fname = sprintf('test53_strong_map__%s__%s__f%g__M%g__%s.png', ...
    char(sanitize_token(T.case_id(1))), char(sanitize_token(T.field_regime_ood(1))), ...
    T.f0(1), T.M(1), char(sanitize_token(fold.fold_key)));
exportgraphics(fig, fullfile(folder, fname), 'Resolution', 220);
savefig(fig, fullfile(folder, replace(fname, ".png", ".fig")));
close(fig);
end

function [Z, xg, zg] = table_to_grid(T, varname)
x = T.x_center_m * 1e3;
z = T.z_center_m * 1e3;
xg = unique(x, 'sorted');
zg = unique(z, 'sorted');
Z = nan(numel(zg), numel(xg));
[~, ix] = ismember(x, xg);
[~, iz] = ismember(z, zg);
v = T.(varname);
for i = 1:height(T)
    if ix(i) > 0 && iz(i) > 0
        Z(iz(i), ix(i)) = v(i);
    end
end
end

function plot_map_panel(xg, zg, Z, ttl, cblabel, limits)
nexttile;
imagesc(xg, zg, Z);
axis image tight;
set(gca, 'YDir','normal');
title(ttl);
xlabel('x (mm)');
ylabel('z (mm)');
cb = colorbar;
cb.Label.String = cblabel;
if nargin >= 6 && numel(limits) == 2 && all(isfinite(limits)) && limits(2) > limits(1)
    clim(limits);
end
set(gca, 'Color', [0.94 0.94 0.94]);
end

function lim = finite_prctile(x, p)
x = x(isfinite(x));
if isempty(x)
    lim = [0 1];
    return;
end
lim = prctile(x, p);
if isscalar(lim)
    lim = double(lim);
else
    lim = double(lim(:)');
end
if numel(lim) == 2 && lim(2) <= lim(1)
    pad = max(abs(lim(1))*0.05, 1e-6);
    lim = [lim(1)-pad lim(2)+pad];
end
end

function lim = finite_symmetric_limit(x, p)
x = abs(x(isfinite(x)));
if isempty(x)
    lim = [-1 1];
    return;
end
a = prctile(x, p);
a = max(double(a), 1);
lim = [-a a];
end

function overlay_rois(T, CFG)
ROI = define_rois_for_condition(T, CFG);
for ri = 1:numel(ROI)
    if isfinite(ROI(ri).side_mm)
        x0 = ROI(ri).center_x_mm - ROI(ri).side_mm/2;
        z0 = ROI(ri).center_z_mm - ROI(ri).side_mm/2;
        rectangle('Position', [x0 z0 ROI(ri).side_mm ROI(ri).side_mm], ...
            'EdgeColor', 'w', 'LineWidth', 1.8, 'LineStyle', '--');
        text(ROI(ri).center_x_mm, ROI(ri).center_z_mm, char(ROI(ri).roi_type), ...
            'Color','w', 'FontWeight','bold', 'HorizontalAlignment','center', ...
            'BackgroundColor',[0 0 0 0.35], 'Interpreter','none');
    end
end
end

function plot_split_box(T_summary, OUT, metric, ttl, ylab, fname)
T = T_summary(T_summary.group_type == "overall" & ...
    T_summary.model_name ~= "theory_discrete", :);
fig = figure('Color','w','Position',[80 80 1200 520]);
boxchart(categorical(T.split_name), T.(metric), 'GroupByColor', categorical(T.model_name));
grid on; ylabel(ylab); title(ttl, 'Interpreter','none');
legend('Location','eastoutside', 'Interpreter','none');
saveas(fig, fullfile(OUT.figure_dir, fname));
close(fig);
end

function plot_fold_lines(T_summary, OUT, split_name, fname)
T = T_summary(T_summary.group_type == "overall" & ...
    T_summary.split_name == split_name & T_summary.model_name ~= "theory_discrete", :);
if isempty(T), return; end
fig = figure('Color','w','Position',[80 80 1100 520]);
models = unique(T.model_name, 'stable');
hold on;
for mi = 1:numel(models)
    X = T(T.model_name == models(mi), :);
    plot(categorical(X.fold_name), X.MAPE_pct, '-o', ...
        'DisplayName', pretty_model(models(mi)), 'LineWidth', 1.4);
end
grid on; ylabel('MAPE (%)'); title("MAPE by fold: " + split_name, 'Interpreter','none');
legend('Location','eastoutside');
saveas(fig, fullfile(OUT.figure_dir, fname));
close(fig);
end

function plot_purity(T_summary, OUT, metric, ttl, ylab, fname)
T = T_summary(T_summary.group_type == "patch_purity_bin" & ...
    ismember(T_summary.model_name, ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"]), :);
if isempty(T), return; end
A = aggregate_over_folds(T, ["model_name","group_value"]);
fig = figure('Color','w','Position',[80 80 1000 520]);
models = unique(A.model_name, 'stable');
hold on;
for mi = 1:numel(models)
    X = A(A.model_name == models(mi), :);
    [~, ord] = sort(string(X.group_value));
    X = X(ord,:);
    plot(categorical(X.group_value), X.(metric + "_mean"), '-o', ...
        'DisplayName', pretty_model(models(mi)), 'LineWidth', 1.4);
end
grid on; ylabel(ylab); xlabel('True patch-purity bin'); title(ttl);
legend('Location','eastoutside');
saveas(fig, fullfile(OUT.figure_dir, fname));
close(fig);
end

function plot_roi_bars(T_roi, OUT)
T = T_roi(ismember(T_roi.model_name, ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"]), :);
if isempty(T), return; end
A = aggregate_over_folds(T, ["model_name","roi_type"]);
roi_types = unique(string(A.roi_type), 'stable');
models = unique(string(A.model_name), 'stable');
Y = nan(numel(roi_types), numel(models));
for ri = 1:numel(roi_types)
    for mi = 1:numel(models)
        idx = string(A.roi_type) == roi_types(ri) & string(A.model_name) == models(mi);
        if any(idx)
            Y(ri, mi) = A.MAPE_pct_mean(find(idx, 1));
        end
    end
end
fig = figure('Color','w','Position',[80 80 1000 520]);
bar(categorical(roi_types), Y);
grid on; ylabel('MAPE (%)'); title('ROI MAPE by model');
legend(arrayfun(@pretty_model, models), 'Location','eastoutside', 'Interpreter','none');
saveas(fig, fullfile(OUT.figure_dir, 'test53_strong_roi_mape.png'));
close(fig);
end

function write_results_readme(T_summary, T_roi, fold_runtime, BENCH, FOLDS, CFG, OUT)
S = aggregate_over_folds(T_summary(T_summary.group_type == "overall", :), ...
    ["split_name","model_name"]);
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w');
fprintf(fid, '# Test 53 Strong-Split Results\n\n');
fprintf(fid, 'Mode: `%s`\n\n', CFG.Mode);
fprintf(fid, '## Splits run\n\n');
for s = unique(FOLDS.split_name)'
    fprintf(fid, '- `%s`: %d folds\n', s, sum(FOLDS.split_name == s));
end
fprintf(fid, '\n## Runtime\n\n');
fprintf(fid, '- Mini-benchmark estimate: %.1f min\n', BENCH.estimated_full_seconds/60);
if ~isempty(fold_runtime)
    fprintf(fid, '- Actual fold runtime: %.1f min\n', sum(fold_runtime.total_seconds)/60);
end
fprintf(fid, '\n## Overall summary\n\n');
fprintf(fid, '| Split | Model | MAPE mean | MAPE std | High >20 mean |\n');
fprintf(fid, '|---|---|---:|---:|---:|\n');
for i = 1:height(S)
    fprintf(fid, '| %s | %s | %.3f | %.3f | %.3f |\n', ...
        S.split_name(i), pretty_model(S.model_name(i)), ...
        S.MAPE_pct_mean(i), S.MAPE_pct_std(i), S.high_error20_pct_mean(i));
end
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- Tables: `%s`\n', OUT.table_dir);
fprintf(fid, '- Figures: `%s`\n', OUT.figure_dir);
fprintf(fid, '- Fold runtime log: `%s`\n', fullfile(OUT.table_dir, 'fold_runtime.csv'));
fprintf(fid, '\n## Interpretation placeholder\n\n');
fprintf(fid, 'Review `summary_metrics_by_split.csv`, purity tables, and ROI tables. ');
fprintf(fid, 'The key paper question is whether composition improves robustly under ');
fprintf(fid, 'frequency, geometry, and field-regime OOD splits without increasing ');
fprintf(fid, 'mixed/interface failures.\n');
fclose(fid);
end

%% Validation-only

function validate_script(CFG, OUT)
feature_dir = fullfile(fileparts(CFG.DatasetFile), 'condition_features');
D = dir(fullfile(feature_dir, 'features__*.mat'));
assert(~isempty(D), 'No condition feature files found for validation.');
S = load(fullfile(D(1).folder, D(1).name), 'T_condition');
T = prepare_dataset(S.T_condition);
features = select_base_predictors(T);
assert_no_forbidden_predictors(features);
assert(all(isfinite(T.q_oracle)), 'q_oracle contains non-finite values.');
assert(~any(contains(lower(features), "purity")), 'Leakage: purity in predictors.');
writetable(table(string(D(1).name), height(T), numel(features), ...
    'VariableNames', {'sample_file','N_rows','N_base_features'}), ...
    fullfile(OUT.table_dir, 'validate_only_summary.csv'));
end

%% Utility functions

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
fprintf('Leakage guard passed: %d base predictors.\n', numel(features));
end

function P = predict_composition(MIX, T, base_features)
X = T(:, cellstr(base_features));
P = struct();
P.predicted_patch_purity = clamp01(predict(MIX.purity, X));
P.p_mixed = class_probability(MIX.mixed, X);
P.p_strong_mixed = class_probability(MIX.strong_mixed, X);
end

function p = class_probability(model, X)
[~, score] = predict(model, X);
classes = string(model.ClassNames);
idx = find(classes == "1" | lower(classes) == "true", 1);
if isempty(idx)
    idx = size(score,2);
end
p = clamp01(score(:, idx));
end

function y = q_to_sws(mapping, q, f0)
y = nan(numel(q), 1);
for i = 1:numel(q)
    if isscalar(f0)
        fi = f0;
    else
        fi = f0(i);
    end
    y(i) = adaptive_req.quantile.quantile_to_cs(mapping{i}, q(i), fi);
end
end

function x = clamp01(x)
x = max(0, min(1, x));
end

function keep = sample_mask_by_condition(T, base_mask, max_rows, seed)
idx_all = find(base_mask);
keep = false(height(T),1);
if max_rows <= 0 || numel(idx_all) <= max_rows
    keep(idx_all) = true;
    return;
end
rng(seed, 'twister');
cond = T.condition_group(idx_all);
[G,~] = findgroups(cond);
ncond = max(G);
per_cond = max(1, floor(max_rows / ncond));
selected = false(numel(idx_all),1);
for gi = 1:ncond
    local = find(G == gi);
    n = min(numel(local), per_cond);
    selected(local(randperm(numel(local), n))) = true;
end
remaining = find(~selected);
n_extra = max_rows - sum(selected);
if n_extra > 0 && ~isempty(remaining)
    extra = remaining(randperm(numel(remaining), min(n_extra, numel(remaining))));
    selected(extra) = true;
end
keep(idx_all(selected)) = true;
end

function fam = paper_geometry_family(case_id)
case_id = string(case_id);
hom = ["homogeneous_cs1p5","homogeneous_cs2","homogeneous_cs2p5", ...
    "homogeneous_cs3","homogeneous_cs3p5","homogeneous_cs4p0"];
bil = ["bilayer_1p5_3","bilayer_2_3","bilayer_2_4","bilayer_2p5_3p5", ...
    "oblique_bilayer_2_4","smooth_bilayer_2_4"];
inc = ["inclusion_1p5_3","inclusion_2_3","inclusion_2_4", ...
    "inclusion_2p5_3p5","offcenter_inclusion_2_4", ...
    "ellipse_2_4","two_inclusions_2_4"];
if ismember(case_id, hom)
    fam = "homogeneous";
elseif ismember(case_id, bil)
    fam = "bilayer";
elseif ismember(case_id, inc)
    fam = "inclusion";
else
    fam = "complex";
end
end

function fam = paper_field_family(regime_id)
r = string(regime_id);
if startsWith(r, "directional_2D")
    fam = "directional_2D";
elseif startsWith(r, "diffuse_2D")
    fam = "diffuse_2D";
elseif startsWith(r, "diffuse_3D")
    fam = "diffuse_3D";
elseif startsWith(r, "partial_3D")
    fam = "partial_3D";
else
    fam = "other";
end
end

function b = purity_bin_strong(p)
labels = strings(size(p));
edges = [0 0.50 0.75 0.90 0.95 0.99 1.0000001];
names = ["0-0.50","0.50-0.75","0.75-0.90","0.90-0.95","0.95-0.99","0.99-1.00"];
for i = 1:numel(names)
    labels(p >= edges(i) & p < edges(i+1)) = names(i);
end
labels(labels == "") = "unknown";
b = categorical(labels, names, 'Ordinal', true);
end

function g = purity_group(p)
g = strings(size(p));
g(p >= 0.95) = "pure_ge_0p95";
g(p < 0.95) = "mixed_lt_0p95";
g(p < 0.75) = "strong_mixed_lt_0p75";
g(g == "") = "unknown";
g = categorical(g);
end

function T = vertcat_nonempty(parts)
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    T = table();
else
    T = vertcat(parts{:});
end
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Plot failed (%s): %s', label, ME.message);
end
end

function name = pretty_model(name)
name = string(name);
name = replace(name, "q_spectrum_only", "q spectrum only");
name = replace(name, "q_spectrum_plus_composition", "q spectrum + composition");
name = replace(name, "q_spectrum_plus_theory_composition", "q spectrum + theory + composition");
name = replace(name, "theory_discrete", "Theory discrete");
name = replace(name, "fixed_q_train_median", "Fixed train-median q");
end

function s = sanitize_token(s)
s = regexprep(string(s), '[^A-Za-z0-9_]+', '_');
end

function x = env_number(name, default_value)
raw = getenv(name);
if isempty(raw)
    x = default_value;
else
    x = str2double(raw);
    if ~isfinite(x), x = default_value; end
end
end

function tf = env_true(name, default_value)
raw = lower(strtrim(string(getenv(name))));
if raw == ""
    tf = logical(default_value);
else
    tf = ismember(raw, ["1","true","yes","y","on"]);
end
end
