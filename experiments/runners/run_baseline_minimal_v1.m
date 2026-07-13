%% run_baseline_minimal_v1.m
% Baseline minimal v1: clean spectral q training for paper closure.
%
% This script is intentionally derived from the Test 38/Test 53 training
% engine, but it trains only the two deployable baselines:
%
%   1. q_spectrum_only
%   2. q_spectrum_plus_composition
%
% The composition block is trained inside this experiment and saved as an
% auxiliary artifact. No frozen T18/old/confidence/correction models are
% used as predictors.
%
% Runtime controls:
%   ADAPTIVE_REQ_BASELINE_MODE                  = quick | medium | full
%   ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY         = true | false
%   ADAPTIVE_REQ_BASELINE_SAVE_ALL_MAPS         = true | false
%   ADAPTIVE_REQ_BASELINE_USE_PARFOR            = true | false
%   ADAPTIVE_REQ_BASELINE_USE_PARALLEL_TRAINING = true | false
%   ADAPTIVE_REQ_BASELINE_OUTPUT_NAME           = optional output folder name
%   ADAPTIVE_REQ_BASELINE_FREQUENCIES_HZ        = optional numeric list
%   ADAPTIVE_REQ_BASELINE_M_LIST                = optional numeric list
%   ADAPTIVE_REQ_BASELINE_CS_GUESS              = optional scalar, default 3
%   ADAPTIVE_REQ_BASELINE_TARGET_STEP_M         = optional scalar, default 2e-3

clear; clc; close all;
format compact;

%% Setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.08);

CFG = default_config();
CFG = apply_json_config_if_requested(CFG, root_dir);
CFG = apply_env_overrides_after_json(CFG);
CFG = setup_parallel_if_requested(CFG);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir,'baseline_minimal_v1_configuration.json'));

fprintf('\nBaseline minimal v1: clean spectral q training\n');
fprintf('Mode: %s | validate only: %d | save all maps: %d\n', ...
    CFG.Mode, CFG.ValidateOnly, CFG.SaveAllMaps);
fprintf('Parallel REQ windows: %d | parallel training: %d\n', ...
    CFG.UseParfor, CFG.UseParallelTraining);
fprintf('Training q_spectrum_only and q_spectrum_plus_composition only.\n');
fprintf('Composition auxiliary models are trained from the same train split.\n');
fprintf('No previous q/confidence/correction predictors.\n');
fprintf('At least one source is forced/aligned in every regime.\n');

if CFG.ValidateOnly
    validate_test38(root_dir, CFG);
    fprintf('Baseline minimal validation-only checks passed.\n');
    return;
end

%% Build/load velocity-diverse training dataset

CASE_SPECS = build_case_specs(CFG);
REGIMES = build_regime_specs(CFG);
fprintf('Training cases: %d | regimes: %d | frequencies: %d | M: %d\n', ...
    numel(CASE_SPECS), numel(REGIMES), numel(CFG.Frequencies), numel(CFG.M));
fprintf('Spatial resolutions dx=dz: %s mm | target REQ step: %.3f mm\n', ...
    mat2str(CFG.DxList*1e3), CFG.TargetStepM*1e3);
if CFG.MaxPatchesPerCondition > 0
    fprintf('REQ rows capped after step-grid extraction: %d patches/condition.\n', ...
        CFG.MaxPatchesPerCondition);
else
    fprintf('REQ rows uncapped: retaining every valid StepX/StepZ center.\n');
end

dataset_file = fullfile(OUT.data_dir, 'baseline_minimal_v1_dataset.mat');
if exist(dataset_file,'file') == 2
    S = load(dataset_file, 'T_train', 'CFG_saved');
    if isfield(S,'CFG_saved') && cache_signature(S.CFG_saved) == cache_signature(CFG)
        T_train = S.T_train;
        fprintf('Reused baseline dataset: %d patches.\n', height(T_train));
    else
        fprintf('Ignoring incompatible baseline dataset cache.\n');
        T_train = build_training_dataset(CASE_SPECS, REGIMES, CFG, OUT);
        CFG_saved = CFG; %#ok<NASGU>
        save(dataset_file, 'T_train', 'CFG_saved', '-v7.3');
    end
else
    T_train = build_training_dataset(CASE_SPECS, REGIMES, CFG, OUT);
    CFG_saved = CFG; %#ok<NASGU>
    save(dataset_file, 'T_train', 'CFG_saved', '-v7.3');
end

assert(~isempty(T_train), 'Baseline minimal dataset is empty.');
[train_mask, test_mask] = condition_split(T_train, CFG);
BASE_FEATURES = select_base_predictors(T_train);
assert_no_forbidden_predictors(BASE_FEATURES);
fprintf('Train patches: %d | held-out patches: %d | conditions: %d | predictors: %d\n', ...
    sum(train_mask), sum(test_mask), numel(unique(T_train.condition_key)), numel(BASE_FEATURES));
T_train.is_train_row = train_mask;
T_train.is_heldout_row = test_mask;

model_train_mask = sample_mask_by_condition(T_train, train_mask, ...
    CFG.MaxRowsForModelTraining, CFG.RandomSeed + 101);
model_eval_mask = sample_mask_by_condition(T_train, test_mask, ...
    CFG.MaxRowsForModelEvaluation, CFG.RandomSeed + 202);
fprintf(['Rows used for model fitting: %d/%d train rows | ', ...
    'rows used for held-out prediction tables: %d/%d held-out rows.\n'], ...
    sum(model_train_mask), sum(train_mask), sum(model_eval_mask), sum(test_mask));
T_model_train = T_train(model_train_mask,:);
T_model_eval = T_train(model_eval_mask,:);
clear T_train;

%% Train minimal model family

MODELS = train_models(T_model_train, true(height(T_model_train),1), BASE_FEATURES, CFG);
clear T_model_train;
T_pred = apply_models(T_model_eval, MODELS, BASE_FEATURES, "baseline_minimal_v1");
clear T_model_eval;
T_held = T_pred(T_pred.is_heldout_row,:);

%% Tables

T_comp = composition_metrics(T_held);
T_overall = summarize_predictions(T_held, "model_name");
T_by_case = summarize_predictions(T_held, ["model_name","case_id"]);
T_by_family = summarize_predictions(T_held, ["model_name","case_family"]);
T_by_regime = summarize_predictions(T_held, ["model_name","field_regime_ood"]);
T_by_frequency = summarize_predictions(T_held, ["model_name","f0"]);
T_by_M = summarize_predictions(T_held, ["model_name","M"]);
T_by_purity = summarize_predictions(T_held, ["model_name","purity_bin"]);
T_by_geometry_frequency = summarize_predictions(T_held, ["model_name","case_id","f0"]);
T_by_roi_frequency = summarize_predictions(T_held, ["model_name","roi_label","f0"]);
T_by_geometry_roi_frequency = summarize_predictions(T_held, ...
    ["model_name","case_id","roi_label","f0"]);
T_fail = top_failure_conditions( ...
    summarize_predictions(T_held, ["model_name","condition_key","case_id", ...
    "case_family","field_regime_ood","f0","M"]), 30);

writetable(remove_cell_columns(T_pred), fullfile(OUT.table_dir, ...
    'baseline_minimal_v1_patch_level_predictions.csv'));
writetable(T_comp, fullfile(OUT.table_dir, 'baseline_minimal_v1_composition_metrics.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_case.csv'));
writetable(T_by_family, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_family.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_regime.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_frequency.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_M.csv'));
writetable(T_by_purity, fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_purity_bin.csv'));
writetable(T_by_geometry_frequency, fullfile(OUT.table_dir, ...
    'baseline_minimal_v1_summary_by_geometry_frequency.csv'));
writetable(T_by_roi_frequency, fullfile(OUT.table_dir, ...
    'baseline_minimal_v1_summary_by_roi_frequency.csv'));
writetable(T_by_geometry_roi_frequency, fullfile(OUT.table_dir, ...
    'baseline_minimal_v1_summary_by_geometry_roi_frequency.csv'));
writetable(T_fail, fullfile(OUT.table_dir, 'baseline_minimal_v1_worst_conditions.csv'));

MODEL_BUNDLE = struct('MODELS', MODELS, 'BASE_FEATURES', BASE_FEATURES, ...
    'CFG', CFG, 'policy', "Baseline minimal v1 clean spectral q");
save(fullfile(OUT.model_dir, 'baseline_minimal_v1_q_models.mat'), ...
    'MODEL_BUNDLE', '-v7.3');
COMPOSITION_AUXILIARY = MODELS.composition; %#ok<NASGU>
save(fullfile(OUT.model_dir, 'composition_auxiliary_models.mat'), ...
    'COMPOSITION_AUXILIARY', 'BASE_FEATURES', 'CFG', '-v7.3');
save(fullfile(OUT.data_dir,'baseline_minimal_v1_results.mat'), ...
    'T_pred','T_held','T_comp','T_overall','T_by_case','T_by_family', ...
    'T_by_regime','T_by_frequency','T_by_M','T_by_purity', ...
    'T_by_geometry_frequency','T_by_roi_frequency', ...
    'T_by_geometry_roi_frequency','T_fail','CFG','-v7.3');

%% Figures

safe_plot(@() plot_composition_diagnostics(T_held, OUT), 'composition diagnostics');
safe_plot(@() plot_model_ranking(T_overall, OUT), 'model ranking');
safe_plot(@() plot_summary_bars(T_overall, T_by_family, T_by_purity, OUT), 'summary bars');
safe_plot(@() plot_heatmap_table(T_by_case, "case_id", OUT, 'baseline_minimal_v1_mape_by_case.png'), 'case heatmap');
safe_plot(@() plot_heatmap_table(T_by_regime, "field_regime_ood", OUT, 'baseline_minimal_v1_mape_by_regime.png'), 'regime heatmap');
safe_plot(@() plot_heatmap_table(T_by_frequency, "f0", OUT, 'baseline_minimal_v1_mape_by_frequency.png'), 'frequency heatmap');
safe_plot(@() plot_heatmap_table(T_by_M, "M", OUT, 'baseline_minimal_v1_mape_by_M.png'), 'M heatmap');
safe_plot(@() plot_error_vs_predicted_purity(T_held, OUT), 'error vs predicted purity');
if CFG.SaveAllMaps
    safe_plot(@() plot_training_maps(T_held, OUT, CFG), 'condition maps');
end

print_interpretation(T_overall, T_by_family, T_by_purity, T_fail, OUT);

fprintf('\nTables: %s\nFigures: %s\nModels: %s\nBaseline minimal v1 complete.\n', ...
    OUT.table_dir, OUT.figure_dir, OUT.model_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MODE'))));
if mode == ""
    mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST38_MODE'))));
end
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","medium","full"]), ...
    'ADAPTIVE_REQ_BASELINE_MODE must be quick, medium, or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.MediumMode = mode == "medium";
CFG.FullMode = mode == "full";
CFG.ValidateOnly = env_true_any(["ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY", ...
    "ADAPTIVE_REQ_TEST38_VALIDATE_ONLY"], false);
CFG.SaveAllMaps = env_true_any(["ADAPTIVE_REQ_BASELINE_SAVE_ALL_MAPS", ...
    "ADAPTIVE_REQ_TEST38_SAVE_ALL_MAPS"], true);
CFG.RandomSeed = 53001;
CFG.DxList = env_number_list_any(["ADAPTIVE_REQ_BASELINE_DX_LIST_MM", ...
    "ADAPTIVE_REQ_TEST38_DX_LIST_MM"], 0.2) * 1e-3;
CFG.dx = CFG.DxList(1);
CFG.dz = CFG.dx;
CFG.Lx = 0.05;
CFG.Lz = 0.05;
CFG.TargetStepM = env_number_any(["ADAPTIVE_REQ_BASELINE_TARGET_STEP_M", ...
    "ADAPTIVE_REQ_TEST38_TARGET_STEP_M"], 2.0e-3);
CFG.cs_guess = env_number_any(["ADAPTIVE_REQ_BASELINE_CS_GUESS", ...
    "ADAPTIVE_REQ_TEST38_CS_GUESS"], 3);
CFG.REQ.Nbins = 'auto';
CFG.REQ.Nbins_auto_oversample = 1;
CFG.REQ.Nbins_min = 16;
CFG.REQ.SmoothSigma = 1;
CFG.REQ.Gamma = 1;
CFG.REQ.PadFactor = 1;
CFG.REQ.EdgeMode = 'valid';
CFG.PurityMixedThreshold = 0.95;
CFG.PurityStrongMixedThreshold = 0.75;
CFG.DistanceEdgesMm = [0 0.5 1 2 4 8 Inf];
CFG.PhysicalSwsRange = [0.5 10];
CFG.UseParfor = env_true_any(["ADAPTIVE_REQ_BASELINE_USE_PARFOR", ...
    "ADAPTIVE_REQ_TEST38_USE_PARFOR"], false);
CFG.UseParallelTraining = env_true_any(["ADAPTIVE_REQ_BASELINE_USE_PARALLEL_TRAINING", ...
    "ADAPTIVE_REQ_TEST38_USE_PARALLEL_TRAINING"], false);
CFG.ParallelPoolType = "threads";
CFG.TrainFraction = 0.70;
CFG.MaxPatchesPerCondition = env_number_any(["ADAPTIVE_REQ_BASELINE_MAX_PATCHES_PER_CONDITION", ...
    "ADAPTIVE_REQ_TEST38_MAX_PATCHES_PER_CONDITION"], 300);
CFG.MaxRowsForModelTraining = env_number_any(["ADAPTIVE_REQ_BASELINE_MAX_MODEL_TRAIN_ROWS", ...
    "ADAPTIVE_REQ_TEST38_MAX_MODEL_TRAIN_ROWS"], 250000);
CFG.MaxRowsForModelEvaluation = env_number_any(["ADAPTIVE_REQ_BASELINE_MAX_MODEL_EVAL_ROWS", ...
    "ADAPTIVE_REQ_TEST38_MAX_MODEL_EVAL_ROWS"], 250000);
if CFG.QuickMode
    CFG.MaxPatchesPerCondition = 550;
    CFG.MaxRowsForModelTraining = env_number_any(["ADAPTIVE_REQ_BASELINE_MAX_MODEL_TRAIN_ROWS", ...
        "ADAPTIVE_REQ_TEST38_MAX_MODEL_TRAIN_ROWS"], 0);
    CFG.MaxRowsForModelEvaluation = env_number_any(["ADAPTIVE_REQ_BASELINE_MAX_MODEL_EVAL_ROWS", ...
        "ADAPTIVE_REQ_TEST38_MAX_MODEL_EVAL_ROWS"], 0);
end
CFG.DiagnosticMapInterpScale = env_number_any(["ADAPTIVE_REQ_BASELINE_MAP_INTERP_SCALE", ...
    "ADAPTIVE_REQ_TEST38_MAP_INTERP_SCALE"], 4);
CFG.TreeLearners = 180;
CFG.MinLeafSize = 8;
if CFG.QuickMode
    CFG.Frequencies = 450;
    CFG.M = 2;
    CFG.RegimeSet = "quick";
    CFG.TreeLearners = 80;
    CFG.MinLeafSize = 10;
elseif CFG.MediumMode
    CFG.Frequencies = [300 500];
    CFG.M = 2;
    CFG.RegimeSet = "full_regimes_medium_freq_M2";
    CFG.TreeLearners = 140;
else
    CFG.Frequencies = [200 300 400 500 600];
    CFG.M = [2 3];
    CFG.DxList = env_number_list_any(["ADAPTIVE_REQ_BASELINE_DX_LIST_MM", ...
        "ADAPTIVE_REQ_TEST38_DX_LIST_MM"], 0.2) * 1e-3;
    CFG.dx = CFG.DxList(1);
    CFG.dz = CFG.dx;
    CFG.MaxPatchesPerCondition = env_number_any( ...
        ["ADAPTIVE_REQ_BASELINE_MAX_PATCHES_PER_CONDITION", ...
        "ADAPTIVE_REQ_TEST38_MAX_PATCHES_PER_CONDITION"], 0);
    CFG.MaxRowsForModelTraining = env_number_any( ...
        ["ADAPTIVE_REQ_BASELINE_MAX_MODEL_TRAIN_ROWS", ...
        "ADAPTIVE_REQ_TEST38_MAX_MODEL_TRAIN_ROWS"], 500000);
    CFG.MaxRowsForModelEvaluation = env_number_any( ...
        ["ADAPTIVE_REQ_BASELINE_MAX_MODEL_EVAL_ROWS", ...
        "ADAPTIVE_REQ_TEST38_MAX_MODEL_EVAL_ROWS"], 500000);
    CFG.RegimeSet = "full";
end
CFG.Frequencies = env_number_list_any(["ADAPTIVE_REQ_BASELINE_FREQUENCIES_HZ", ...
    "ADAPTIVE_REQ_TEST38_FREQUENCIES_HZ"], CFG.Frequencies);
CFG.M = env_number_list_any(["ADAPTIVE_REQ_BASELINE_M_LIST", ...
    "ADAPTIVE_REQ_TEST38_M_LIST"], CFG.M);
CFG.OutputName = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_OUTPUT_NAME'))));
if CFG.OutputName == ""
    CFG.OutputName = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST38_OUTPUT_NAME'))));
end
if CFG.OutputName == ""
    CFG.OutputName = "baseline_minimal_v1";
end
CFG.CaseIDs = strings(0,1);
CFG.RegimeIDs = strings(0,1);
CFG.ConfigPath = "";
end

function CFG = apply_json_config_if_requested(CFG, root_dir)
cfg_path = string(getenv('ADAPTIVE_REQ_BASELINE_CONFIG'));
if cfg_path == ""
    cfg_path = fullfile(root_dir, 'configs', 'final_training', ...
        'baseline_minimal_v1.json');
end
if exist(cfg_path, 'file') ~= 2
    return;
end
J = jsondecode(fileread(cfg_path));
CFG.ConfigPath = cfg_path;
CFG = assign_json_scalar(CFG, J, 'Mode');
CFG.Mode = lower(string(CFG.Mode));
CFG.QuickMode = CFG.Mode == "quick";
CFG.MediumMode = CFG.Mode == "medium";
CFG.FullMode = CFG.Mode == "full";
CFG = assign_json_scalar(CFG, J, 'ValidateOnly');
CFG = assign_json_scalar(CFG, J, 'SaveAllMaps');
CFG = assign_json_scalar(CFG, J, 'RandomSeed');
CFG = assign_json_scalar(CFG, J, 'TrainFraction');
CFG = assign_json_scalar(CFG, J, 'TargetStepM');
CFG = assign_json_scalar(CFG, J, 'cs_guess');
CFG = assign_json_scalar(CFG, J, 'Lx');
CFG = assign_json_scalar(CFG, J, 'Lz');
CFG = assign_json_scalar(CFG, J, 'MaxPatchesPerCondition');
CFG = assign_json_scalar(CFG, J, 'MaxRowsForModelTraining');
CFG = assign_json_scalar(CFG, J, 'MaxRowsForModelEvaluation');
CFG = assign_json_scalar(CFG, J, 'UseParfor');
CFG = assign_json_scalar(CFG, J, 'UseParallelTraining');
CFG = assign_json_scalar(CFG, J, 'TreeLearners');
CFG = assign_json_scalar(CFG, J, 'MinLeafSize');
CFG = assign_json_scalar(CFG, J, 'OutputName');
if isfield(J, 'Frequencies'), CFG.Frequencies = rowvec(J.Frequencies); end
if isfield(J, 'M'), CFG.M = rowvec(J.M); end
if isfield(J, 'DxListMm')
    CFG.DxList = rowvec(J.DxListMm) * 1e-3;
    CFG.dx = CFG.DxList(1); CFG.dz = CFG.dx;
end
if isfield(J, 'CaseIDs'), CFG.CaseIDs = string(J.CaseIDs(:)); end
if isfield(J, 'RegimeIDs'), CFG.RegimeIDs = string(J.RegimeIDs(:)); end
if isfield(J, 'REQ')
    R = J.REQ;
    if isfield(R, 'Nbins'), CFG.REQ.Nbins = R.Nbins; end
    if isfield(R, 'Nbins_auto_oversample'), CFG.REQ.Nbins_auto_oversample = R.Nbins_auto_oversample; end
    if isfield(R, 'Nbins_min'), CFG.REQ.Nbins_min = R.Nbins_min; end
    if isfield(R, 'SmoothSigma'), CFG.REQ.SmoothSigma = R.SmoothSigma; end
    if isfield(R, 'Gamma'), CFG.REQ.Gamma = R.Gamma; end
    if isfield(R, 'PadFactor'), CFG.REQ.PadFactor = R.PadFactor; end
    if isfield(R, 'EdgeMode'), CFG.REQ.EdgeMode = char(R.EdgeMode); end
end
end

function CFG = apply_env_overrides_after_json(CFG)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MODE'))));
if mode ~= ""
    assert(ismember(mode, ["quick","medium","full"]), ...
        'ADAPTIVE_REQ_BASELINE_MODE must be quick, medium, or full.');
    CFG.Mode = mode;
    CFG.QuickMode = mode == "quick";
    CFG.MediumMode = mode == "medium";
    CFG.FullMode = mode == "full";
end
if strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY'))) ~= ""
    CFG.ValidateOnly = env_true('ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY', CFG.ValidateOnly);
end
if strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_SAVE_ALL_MAPS'))) ~= ""
    CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_BASELINE_SAVE_ALL_MAPS', CFG.SaveAllMaps);
end
CFG.Frequencies = env_number_list('ADAPTIVE_REQ_BASELINE_FREQUENCIES_HZ', CFG.Frequencies);
CFG.M = env_number_list('ADAPTIVE_REQ_BASELINE_M_LIST', CFG.M);
CFG.cs_guess = env_number('ADAPTIVE_REQ_BASELINE_CS_GUESS', CFG.cs_guess);
CFG.TargetStepM = env_number('ADAPTIVE_REQ_BASELINE_TARGET_STEP_M', CFG.TargetStepM);
CFG.MaxPatchesPerCondition = env_number('ADAPTIVE_REQ_BASELINE_MAX_PATCHES_PER_CONDITION', ...
    CFG.MaxPatchesPerCondition);
CFG.MaxRowsForModelTraining = env_number('ADAPTIVE_REQ_BASELINE_MAX_MODEL_TRAIN_ROWS', ...
    CFG.MaxRowsForModelTraining);
CFG.MaxRowsForModelEvaluation = env_number('ADAPTIVE_REQ_BASELINE_MAX_MODEL_EVAL_ROWS', ...
    CFG.MaxRowsForModelEvaluation);
out_name = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_OUTPUT_NAME'))));
if out_name ~= "", CFG.OutputName = out_name; end

% A full JSON config is the default paper recipe. When quick/medium is
% requested explicitly from the environment, keep that run genuinely small
% unless the caller also supplied explicit frequencies/M.
freq_explicit = strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_FREQUENCIES_HZ'))) ~= "";
m_explicit = strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_M_LIST'))) ~= "";
train_rows_explicit = strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MAX_MODEL_TRAIN_ROWS'))) ~= "";
eval_rows_explicit = strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MAX_MODEL_EVAL_ROWS'))) ~= "";
patch_rows_explicit = strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MAX_PATCHES_PER_CONDITION'))) ~= "";
if mode == "quick"
    if ~freq_explicit, CFG.Frequencies = 450; end
    if ~m_explicit, CFG.M = 2; end
    if ~patch_rows_explicit, CFG.MaxPatchesPerCondition = 550; end
    if ~train_rows_explicit, CFG.MaxRowsForModelTraining = 0; end
    if ~eval_rows_explicit, CFG.MaxRowsForModelEvaluation = 0; end
    CFG.TreeLearners = min(CFG.TreeLearners, 80);
    CFG.MinLeafSize = max(CFG.MinLeafSize, 10);
    CFG.CaseIDs = ["homogeneous_cs1p5"; "homogeneous_cs3"; ...
        "homogeneous_cs4p0"; "bilayer_2_4"; "inclusion_2_4"; ...
        "three_material_2_3_4"];
    CFG.RegimeIDs = ["directional_2D_angle0"; "diffuse_2D_seed1"; ...
        "diffuse_3D_seed1"; "partial_3D_8src"];
elseif mode == "medium"
    if ~freq_explicit, CFG.Frequencies = [300 500]; end
    if ~m_explicit, CFG.M = 2; end
    if ~patch_rows_explicit, CFG.MaxPatchesPerCondition = 300; end
    CFG.TreeLearners = min(CFG.TreeLearners, 140);
end
end

function CFG = assign_json_scalar(CFG, J, name)
if isfield(J, name)
    CFG.(name) = J.(name);
end
end

function v = rowvec(v)
v = double(v);
v = v(:)';
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', char(CFG.OutputName));
if CFG.QuickMode
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
elseif CFG.MediumMode
    OUT.root_dir = fullfile(OUT.root_dir, 'medium');
end
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.field_dir = fullfile(OUT.data_dir, 'field_cache');
OUT.condition_dir = fullfile(OUT.data_dir, 'condition_features');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.model_dir = fullfile(OUT.root_dir, 'models');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.fail_dir = fullfile(OUT.figure_dir, 'worst_conditions');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function assert_no_forbidden_predictors(features)
bad_patterns = ["true","oracle","purity","mixed","confidence","error", ...
    "pred","sws","cs_","k_true","q_local","q_pred","q_theory", ...
    "req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition"];
features = lower(string(features));
for p = bad_patterns
    hit = features(contains(features, p));
    assert(isempty(hit), 'Forbidden baseline predictor detected: %s', strjoin(hit, ', '));
end
fprintf('Base predictors passed leakage guard (%d predictors).\n', numel(features));
end

function validate_test38(root_dir, CFG)
C = build_case_specs(CFG); R = build_regime_specs(CFG);
CFG.dx = CFG.DxList(1);
CFG.dz = CFG.dx;
cfg = build_sim_cfg(CFG, C(1), R(1), CFG.Frequencies(1), 1);
sim = run_ood_simulation(cfg, C(1));
T = extract_training_condition(sim, cfg, C(1), R(1), CFG.M(1), "validate", CFG);
features = select_base_predictors(T);
assert_no_forbidden_predictors(features);
assert(~isempty(T) && all(isfinite(T.q_oracle)));
fprintf('Validation sample: %s, %d patches, %d clean predictors.\n', ...
    C(1).case_id, height(T), numel(features));
end

function T_all = build_training_dataset(CASE_SPECS, REGIMES, CFG, OUT)
parts = {};
condition_id = 0;
total = numel(CASE_SPECS) * numel(REGIMES) * numel(CFG.Frequencies) * ...
    numel(CFG.DxList) * numel(CFG.M);
for ci = 1:numel(CASE_SPECS)
    C = CASE_SPECS(ci);
    for di = 1:numel(CFG.DxList)
        CFG_dx = CFG;
        CFG_dx.dx = CFG.DxList(di);
        CFG_dx.dz = CFG_dx.dx;
        for fi = 1:numel(CFG.Frequencies)
            f0 = CFG.Frequencies(fi);
            for ri = 1:numel(REGIMES)
                R = REGIMES(ri);
                condition_id = condition_id + 1;
                field_key = sprintf('%s__f%g__%s__dx%gum', ...
                    C.case_id, f0, lower(R.regime_id), round(1e6*CFG_dx.dx));
                field_file = fullfile(OUT.field_dir, "field__" + sanitize(field_key) + ".mat");
                if exist(field_file,'file') == 2
                    S = load(field_file, 'sim', 'cfg_sim');
                    sim = S.sim; cfg_sim = S.cfg_sim;
                    fprintf('[field] reused %s\n', field_key);
                else
                    cfg_sim = build_sim_cfg(CFG_dx, C, R, f0, condition_id);
                    sim = run_ood_simulation(cfg_sim, C);
                    case_spec = C; regime_spec = R; %#ok<NASGU>
                    save(field_file, 'sim', 'cfg_sim', 'case_spec', 'regime_spec', '-v7.3');
                    fprintf('[field] built %s | size %dx%d\n', field_key, ...
                        size(sim.Uxz,1), size(sim.Uxz,2));
                end
                for mi = 1:numel(CFG.M)
                    M = CFG.M(mi);
                    key = string(field_key) + "__M" + string(M);
                    feature_file = fullfile(OUT.condition_dir, ...
                        "features__" + sanitize(key) + ".mat");
                    condition_signature = condition_feature_signature(CFG_dx, M);
                    if exist(feature_file,'file') == 2
                        S = load(feature_file, 'T_condition', 'condition_signature');
                        if isfield(S, 'condition_signature') && ...
                                isequal(S.condition_signature, condition_signature)
                            T_condition = S.T_condition;
                            fprintf('[%d/%d] reused %s (%d patches).\n', ...
                                numel(parts)+1, total, key, height(T_condition));
                        else
                            timer = tic;
                            T_condition = extract_training_condition(sim, cfg_sim, C, R, M, key, CFG_dx);
                            save(feature_file, 'T_condition', 'condition_signature', '-v7.3');
                            fprintf(['[%d/%d] regenerated %s (%d patches, dx %.3f mm, ', ...
                                'StepX %d, StepZ %d) in %.1f s.\n'], ...
                                numel(parts)+1, total, key, height(T_condition), ...
                                cfg_sim.dx*1e3, T_condition.REQ_StepX(1), ...
                                T_condition.REQ_StepZ(1), toc(timer));
                        end
                    else
                        timer = tic;
                        T_condition = extract_training_condition(sim, cfg_sim, C, R, M, key, CFG_dx);
                        save(feature_file, 'T_condition', 'condition_signature', '-v7.3');
                        fprintf(['[%d/%d] built %s (%d patches, dx %.3f mm, ', ...
                            'StepX %d, StepZ %d) in %.1f s.\n'], ...
                            numel(parts)+1, total, key, height(T_condition), ...
                            cfg_sim.dx*1e3, T_condition.REQ_StepX(1), ...
                            T_condition.REQ_StepZ(1), toc(timer));
                    end
                    parts{end+1,1} = T_condition; %#ok<AGROW>
                end
            end
        end
    end
end
T_all = vertcat(parts{:});
T_all.dataset = repmat("baseline_minimal_v1", height(T_all), 1);
T_all = movevars(T_all, 'dataset', 'Before', 1);
end

function T = extract_training_condition(sim, cfg, C, R, M, key, CFG)
feat = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', CFG.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
step_x = max(1, round(CFG.TargetStepM / cfg.dx));
step_z = max(1, round(CFG.TargetStepM / cfg.dz));
O = adaptive_req.estimators.req_estimator_map(sim.Uxz, cfg, feat, ...
    'StepX', step_x, 'StepZ', step_z, ...
    'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins',CFG.REQ.Nbins, ...
    'Nbins_auto_oversample',CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min',CFG.REQ.Nbins_min,'smooth_sigma',CFG.REQ.SmoothSigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseParfor, ...
    'StoreReqCurves', false, 'Verbose', false);
T = O.feature_table;
T.REQ_StepX = step_x * ones(height(T),1);
T.REQ_StepZ = step_z * ones(height(T),1);
T.TargetStepM = CFG.TargetStepM * ones(height(T),1);
if CFG.MaxPatchesPerCondition > 0 && height(T) > CFG.MaxPatchesPerCondition
    rng(CFG.RandomSeed + sum(double(char(key))), 'twister');
    T = T(sort(randperm(height(T), CFG.MaxPatchesPerCondition)), :);
end
T = attach_metadata(T, sim, cfg, C, R, M, key, O.win_size, CFG);
T.is_mixed = T.patch_purity < CFG.PurityMixedThreshold;
T.is_strong_mixed = T.patch_purity < CFG.PurityStrongMixedThreshold;
end

function sig = condition_feature_signature(CFG, M)
sig = struct();
sig.Version = 2;
sig.M = M;
sig.dx = CFG.dx;
sig.dz = CFG.dz;
sig.TargetStepM = CFG.TargetStepM;
sig.MaxPatchesPerCondition = CFG.MaxPatchesPerCondition;
sig.cs_guess = CFG.cs_guess;
sig.REQ_Nbins = CFG.REQ.Nbins;
sig.REQ_Nbins_auto_oversample = CFG.REQ.Nbins_auto_oversample;
sig.REQ_Nbins_min = CFG.REQ.Nbins_min;
sig.REQ_SmoothSigma = CFG.REQ.SmoothSigma;
sig.REQ_Gamma = CFG.REQ.Gamma;
sig.REQ_PadFactor = CFG.REQ.PadFactor;
sig.REQ_EdgeMode = CFG.REQ.EdgeMode;
end

function [train_mask, test_mask] = condition_split(T, CFG)
rng(CFG.RandomSeed, 'twister');
cond = unique(T(:, {'condition_key','case_id','field_regime_ood','f0','M'}), 'rows', 'stable');
train_cond = false(height(cond),1);
groups = unique(cond(:, {'case_id','field_regime_ood'}), 'rows', 'stable');
for gi = 1:height(groups)
    idx = find(cond.case_id == groups.case_id(gi) & ...
        cond.field_regime_ood == groups.field_regime_ood(gi));
    idx = idx(randperm(numel(idx)));
    ntrain = max(1, round(CFG.TrainFraction * numel(idx)));
    if ntrain >= numel(idx) && numel(idx) > 1
        ntrain = numel(idx)-1;
    end
    train_cond(idx(1:ntrain)) = true;
end
train_keys = cond.condition_key(train_cond);
train_mask = ismember(T.condition_key, train_keys);
test_mask = ~train_mask;
if ~any(test_mask)
    idx = randperm(height(T));
    train_mask = false(height(T),1);
    train_mask(idx(1:round(CFG.TrainFraction*height(T)))) = true;
    test_mask = ~train_mask;
end
end

function keep = sample_mask_by_condition(T, base_mask, max_rows, seed)
keep = false(height(T),1);
idx_all = find(base_mask);
if max_rows <= 0 || numel(idx_all) <= max_rows
    keep(idx_all) = true;
    return;
end
rng(seed, 'twister');
if ~ismember("condition_key", string(T.Properties.VariableNames))
    idx = idx_all(randperm(numel(idx_all), max_rows));
    keep(idx) = true;
    return;
end
[G,~] = findgroups(T.condition_key(idx_all));
groups = unique(G(:))';
per_group = max(1, floor(max_rows / numel(groups)));
selected = false(numel(idx_all),1);
for g = groups
    local = find(G == g);
    n = min(per_group, numel(local));
    local = local(randperm(numel(local), n));
    selected(local) = true;
end
sel_idx = find(selected);
if numel(sel_idx) < max_rows
    remaining = find(~selected);
    n_extra = min(max_rows - numel(sel_idx), numel(remaining));
    if n_extra > 0
        extra = remaining(randperm(numel(remaining), n_extra));
        selected(extra) = true;
    end
elseif numel(sel_idx) > max_rows
    selected(:) = false;
    selected(sel_idx(randperm(numel(sel_idx), max_rows))) = true;
end
keep(idx_all(selected)) = true;
end

%% OOD design

function cases = build_case_specs(CFG)
base = [
    homogeneous_case("homogeneous_cs1p5", 1.5)
    homogeneous_case("homogeneous_cs2", 2.0)
    homogeneous_case("homogeneous_cs2p5", 2.5)
    homogeneous_case("homogeneous_cs3", 3.0)
    homogeneous_case("homogeneous_cs3p5", 3.5)
    homogeneous_case("homogeneous_cs4p0", 4.0)
    bilayer_case("bilayer_1p5_3", 1.5, 3.0, 0)
    bilayer_case("bilayer_2_3", 2.0, 3.0, 0)
    bilayer_case("bilayer_2_4", 2.0, 4.0, 0)
    bilayer_case("bilayer_2p5_3p5", 2.5, 3.5, 0)
    inclusion_case("inclusion_1p5_3", 1.5, 3.0, [0.025 0.025], 0.008)
    inclusion_case("inclusion_2_3", 2.0, 3.0, [0.025 0.025], 0.008)
    inclusion_case("inclusion_2_4", 2.0, 4.0, [0.025 0.025], 0.008)
    inclusion_case("inclusion_2p5_3p5", 2.5, 3.5, [0.025 0.025], 0.008)
    ellipse_case("ellipse_2_4", 2.0, 4.0, [0.025 0.025], [0.013 0.006], deg2rad(25))
    inclusion_case("offcenter_inclusion_2_4", 2.0, 4.0, [0.032 0.022], 0.007)
    two_inclusion_case("two_inclusions_2_4", 2.0, 4.0)
    bilayer_case("oblique_bilayer_2_4", 2.0, 4.0, deg2rad(25))
    thin_layer_case("thin_layer_2_4", 2.0, 4.0)
    three_material_case("three_material_2_3_4", 2.0, 3.0, 4.0)
    smooth_bilayer_case("smooth_bilayer_2_4", 2.0, 4.0, 1.2e-3)
    ];
if CFG.QuickMode
    keep = ["homogeneous_cs1p5","homogeneous_cs3","homogeneous_cs4p0", ...
        "bilayer_2_4","inclusion_2_4","three_material_2_3_4"];
    cases = base(ismember([base.case_id], keep));
else
    cases = base;
end
if isfield(CFG, 'CaseIDs') && ~isempty(CFG.CaseIDs)
    cases = cases(ismember([cases.case_id], CFG.CaseIDs));
end
assert(~isempty(cases), 'No baseline cases selected. Check CFG.CaseIDs.');
end

function R = build_regime_specs(CFG)
if CFG.QuickMode
    R = [
        regime("directional_2D_angle0", "directional_2D", true, 1, true, ...
        "ranges", "fibonacci", [0 0], [0 2*pi], [0 pi], 11)
        regime("directional_2D_angle45", "directional_2D", true, 1, true, ...
        "ranges", "fibonacci", [pi/4 pi/4], [0 2*pi], [0 pi], 13)
        regime("diffuse_2D_seed1", "diffuse_2D", true, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 17)
        regime("diffuse_3D_seed1", "diffuse_3D", false, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 19)
        regime("partial_3D_8src", "partial_3D", false, 8, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 23)
        ];
else
    R = [
        regime("directional_2D_angle0", "directional_2D", true, 1, true, ...
        "ranges", "fibonacci", [0 0], [0 2*pi], [0 pi], 11)
        regime("directional_2D_angle30", "directional_2D", true, 1, true, ...
        "ranges", "fibonacci", [pi/6 pi/6], [0 2*pi], [0 pi], 13)
        regime("directional_2D_angle60", "directional_2D", true, 1, true, ...
        "ranges", "fibonacci", [pi/3 pi/3], [0 2*pi], [0 pi], 17)
        regime("diffuse_2D_seed1", "diffuse_2D", true, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 19)
        regime("diffuse_3D_seed1", "diffuse_3D", false, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 23)
        regime("diffuse_3D_seed2", "diffuse_3D", false, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 29)
        regime("diffuse_3D_seed3", "diffuse_3D", false, 128, true, ...
        "ranges", "random", [0 2*pi], [0 2*pi], [0 pi], 31)
        regime("partial_3D_8src", "partial_3D", false, 8, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 37)
        regime("partial_3D_16src", "partial_3D", false, 16, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 41)
        regime("partial_3D_32src", "partial_3D", false, 32, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 43)
        ];
end
if isfield(CFG, 'RegimeIDs') && ~isempty(CFG.RegimeIDs)
    R = R(ismember([R.regime_id], CFG.RegimeIDs));
end
assert(~isempty(R), 'No field regimes selected. Check CFG.RegimeIDs.');
end

function C = homogeneous_case(id, cs)
C = empty_case(id, "homogeneous");
C.cs_values = cs;
C.cs_bg = cs;
C.mask_builder = "homogeneous";
end

function C = bilayer_case(id, cs_low, cs_high, angle)
C = empty_case(id, "bilayer");
C.cs_values = [cs_low cs_high];
C.cs_bg = cs_low;
C.cs_inc = cs_high;
C.angle = angle;
C.offset = 0.025;
C.sigma_edge = 1e-6;
C.mask_builder = "bilayer";
end

function C = smooth_bilayer_case(id, cs_low, cs_high, sigma_edge)
C = bilayer_case(id, cs_low, cs_high, 0);
C.case_family = "smooth_transition";
C.sigma_edge = sigma_edge;
end

function C = inclusion_case(id, cs_bg, cs_inc, center, radius)
C = empty_case(id, "inclusion");
C.cs_values = [cs_bg cs_inc];
C.cs_bg = cs_bg;
C.cs_inc = cs_inc;
C.center = center;
C.radius = radius;
C.sigma_edge = 1e-6;
C.mask_builder = "circle";
end

function C = ellipse_case(id, cs_bg, cs_inc, center, axes_m, angle)
C = empty_case(id, "ellipse");
C.cs_values = [cs_bg cs_inc];
C.cs_bg = cs_bg;
C.cs_inc = cs_inc;
C.center = center;
C.axes = axes_m;
C.angle = angle;
C.mask_builder = "ellipse_custom";
end

function C = two_inclusion_case(id, cs_bg, cs_inc)
C = empty_case(id, "two_inclusions");
C.cs_values = [cs_bg cs_inc];
C.cs_bg = cs_bg;
C.cs_inc = cs_inc;
C.centers = [0.018 0.024; 0.034 0.028];
C.radii = [0.006 0.005];
C.mask_builder = "two_circles";
end

function C = thin_layer_case(id, cs_bg, cs_layer)
C = empty_case(id, "thin_layer");
C.cs_values = [cs_bg cs_layer];
C.cs_bg = cs_bg;
C.cs_inc = cs_layer;
C.center_x = 0.025;
C.thickness = 0.004;
C.mask_builder = "thin_layer_custom";
end

function C = three_material_case(id, cs_bg, cs_mid, cs_high)
C = empty_case(id, "three_material");
C.cs_values = [cs_bg cs_mid cs_high];
C.cs_bg = cs_bg;
C.cs_mid = cs_mid;
C.cs_high = cs_high;
C.mask_builder = "three_material_custom";
end

function C = empty_case(id, fam)
C = struct();
C.case_id = string(id);
C.case_family = string(fam);
C.cs_values = [];
C.cs_bg = NaN;
C.cs_inc = NaN;
C.mask_builder = "";
C.cs_mid = NaN;
C.cs_high = NaN;
C.angle = NaN;
C.offset = NaN;
C.sigma_edge = 1e-6;
C.center = [NaN NaN];
C.radius = NaN;
C.axes = [NaN NaN];
C.centers = nan(0,2);
C.radii = nan(0,1);
C.center_x = NaN;
C.thickness = NaN;
end

function R = regime(id, label, is2d, nwaves, force_plane, sampling, method, angle2d, phi, theta, seed_offset)
R = struct();
R.regime_id = string(id);
R.field_regime = string(label);
R.Is2D = is2d;
R.Nwaves = nwaves;
R.ForceInPlaneWave = force_plane;
R.SourceSampling = string(sampling);
R.AngularSamplingMethod = string(method);
R.AngleRange2D = angle2d;
R.PhiRange = phi;
R.ThetaRange = theta;
R.SeedOffset = seed_offset;
end

%% Simulation and feature extraction

function cfg = build_sim_cfg(CFG, C, R, f0, condition_id)
cfg = struct();
cfg.Lx = CFG.Lx; cfg.Lz = CFG.Lz;
cfg.dx = CFG.dx; cfg.dz = CFG.dz;
cfg.f0 = f0;
cfg.cs_bg = C.cs_bg;
cfg.cs_inc = max(C.cs_values);
cfg.Nwaves = R.Nwaves;
cfg.Is2D = R.Is2D;
cfg.WaveModel = 'spherical';
cfg.AngularSamplingMethod = char(R.AngularSamplingMethod);
cfg.ForceInPlaneWave = logical(R.ForceInPlaneWave);
cfg.SNR = Inf;
cfg.AmpJitter = 0.05;
cfg.DecayAlpha = 0;
cfg.UseParfor = false;
cfg.Seed = CFG.RandomSeed + 1000*condition_id + R.SeedOffset;
cfg.PhiRange = R.PhiRange;
cfg.ThetaRange = R.ThetaRange;
cfg.AngleRange2D = R.AngleRange2D;
cfg.SourceSampling = char(R.SourceSampling);
cfg.MaskConfig = build_mask_config(CFG, C);
end

function mask_cfg = build_mask_config(CFG, C)
switch C.mask_builder
    case "homogeneous"
        masks = {};
    case "bilayer"
        masks = {struct('Type','bilayer','cs_inc',C.cs_inc,'Params', ...
            struct('Bi_Angle',C.angle,'Bi_Offset',C.offset, ...
            'SigmaEdge',C.sigma_edge))};
    case "circle"
        masks = {struct('Type','circle','cs_inc',C.cs_inc,'Params', ...
            struct('Center',C.center,'Radius',C.radius,'SigmaEdge',C.sigma_edge))};
    case "ellipse_custom"
        masks = {struct('Type','custom','cs_inc',C.cs_inc,'Params', ...
            struct('CustomMask',custom_mask(CFG,C,"ellipse"),'SigmaEdge',1e-6))};
    case "two_circles"
        masks = cell(1,size(C.centers,1));
        for i = 1:numel(masks)
            masks{i} = struct('Type','circle','cs_inc',C.cs_inc,'Params', ...
                struct('Center',C.centers(i,:),'Radius',C.radii(i),'SigmaEdge',1e-6));
        end
    case "thin_layer_custom"
        masks = {struct('Type','custom','cs_inc',C.cs_inc,'Params', ...
            struct('CustomMask',custom_mask(CFG,C,"thin_layer"),'SigmaEdge',1e-6))};
    case "three_material_custom"
        masks = {
            struct('Type','bilayer','cs_inc',C.cs_mid,'Params', ...
            struct('Bi_Angle',0,'Bi_Offset',0.025,'SigmaEdge',1e-6))
            struct('Type','custom','cs_inc',C.cs_high,'Params', ...
            struct('CustomMask',custom_mask(CFG,C,"three_circle"),'SigmaEdge',1e-6))
            };
    otherwise
        error('Unknown mask_builder: %s', C.mask_builder);
end
mask_cfg = struct('cs_bg', C.cs_bg, 'CombineMode', 'overlay', 'Masks', {masks});
end

function M = custom_mask(CFG, C, kind)
Nx = round(CFG.Lx / CFG.dx) + 1;
Nz = round(CFG.Lz / CFG.dz) + 1;
x = linspace(0, CFG.Lx, Nx);
z = linspace(0, CFG.Lz, Nz);
[X,Z] = ndgrid(x,z);
switch string(kind)
    case "ellipse"
        ca = cos(C.angle); sa = sin(C.angle);
        X0 = X - C.center(1); Z0 = Z - C.center(2);
        xr = ca*X0 + sa*Z0;
        zr = -sa*X0 + ca*Z0;
        M = (xr./C.axes(1)).^2 + (zr./C.axes(2)).^2 <= 1;
    case "thin_layer"
        M = abs(X - C.center_x) <= C.thickness/2;
    case "three_circle"
        M = hypot(X - 0.017, Z - 0.031) <= 0.006;
    otherwise
        error('Unknown custom mask: %s', kind);
end
end

function sim = run_ood_simulation(cfg, C)
sim = adaptive_req.simulate.run_single_simulation(cfg);
if C.case_family == "smooth_transition"
    % MaskConfig smoothing already creates a smooth material transition.
end
end

function T = evaluate_condition(sim, cfg, C, R, M, key, MODEL, CFG)
feat = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', CFG.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
step_x = max(1, round(CFG.TargetStepM / cfg.dx));
step_z = max(1, round(CFG.TargetStepM / cfg.dz));
O = adaptive_req.estimators.req_estimator_map(sim.Uxz, cfg, feat, ...
    'StepX', step_x, 'StepZ', step_z, ...
    'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins',CFG.REQ.Nbins, ...
    'Nbins_auto_oversample',CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min',CFG.REQ.Nbins_min,'smooth_sigma',CFG.REQ.SmoothSigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseParfor, ...
    'StoreReqCurves', false, 'Verbose', false);
F = O.feature_table;
F = attach_metadata(F, sim, cfg, C, R, M, key, O.win_size, CFG);
F = ensure_predictor_columns(F, MODEL.base_features);
q_pred = predict(MODEL.q_model.model, F(:, cellstr(MODEL.q_model.features)));
q_pred = clamp01(q_pred);
sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
COMP = predict_composition(MODEL.composition, F, MODEL.base_features);
T = F(:, {'condition_key','case_id','case_family','field_regime','field_regime_ood', ...
    'f0','M','dx','dz','map_iz','map_ix','cx','cz','x_center_m','z_center_m', ...
    'true_SWS','patch_purity','purity_bin','distance_to_boundary_mm','distance_bin', ...
    'q_oracle','q_theory_prior'});
T.model_name = repmat(MODEL.model_name, height(T), 1);
T.q_pred = q_pred;
T.sws_pred = sws_pred;
T.q_error = q_pred - T.q_oracle;
T.sws_signed_error_pct = 100*(sws_pred - T.true_SWS) ./ T.true_SWS;
T.sws_abs_error_pct = abs(T.sws_signed_error_pct);
T.high_error10 = T.sws_abs_error_pct > 10;
T.high_error20 = T.sws_abs_error_pct > 20;
T.predicted_patch_purity = COMP.predicted_patch_purity;
T.p_mixed = COMP.p_mixed;
T.p_strong_mixed = COMP.p_strong_mixed;
T.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
T.sws_theory_error_pct = 100*(T.sws_theory - T.true_SWS) ./ T.true_SWS;
end

function F = attach_metadata(F, sim, cfg, C, R, M, key, win_size, CFG)
n = height(F);
F.condition_key = repmat(string(key), n, 1);
F.case_id = repmat(C.case_id, n, 1);
F.case_family = repmat(C.case_family, n, 1);
F.field_regime = repmat(R.field_regime, n, 1);
F.field_regime_ood = repmat(R.regime_id, n, 1);
F.f0 = cfg.f0 * ones(n,1);
F.M = M * ones(n,1);
F.REQ_M = M * ones(n,1);
F.dx = cfg.dx * ones(n,1);
F.dz = cfg.dz * ones(n,1);
F.SIM_f0 = cfg.f0 * ones(n,1);
F.true_SWS = nan(n,1);
F.k_true = nan(n,1);
F.patch_purity = nan(n,1);
F.q_oracle = nan(n,1);
half = floor(win_size/2);
for i = 1:n
    cx = F.cx(i); cz = F.cz(i);
    F.true_SWS(i) = sim.cs_map(cz, cx);
    F.k_true(i) = 2*pi*cfg.f0/F.true_SWS(i);
    patch = sim.cs_map((cz-half):(cz+half), (cx-half):(cx+half));
    F.patch_purity(i) = dominant_fraction(patch);
    F.q_oracle(i) = invert_mapping_to_q(F.req_mapping{i}, F.k_true(i));
end
F.purity_bin = purity_bin(F.patch_purity);
D = distance_to_material_boundary(sim.cs_map, cfg.dx);
F.distance_to_boundary_mm = D(sub2ind(size(D), F.cz, F.cx));
F.distance_bin = distance_bin(F.distance_to_boundary_mm, CFG.DistanceEdgesMm);
F.q_theory_prior = repmat(theory_q_for_regime(cfg, M, R.field_regime), n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
end

function F = ensure_predictor_columns(F, features)
for f = string(features(:))'
    if ~ismember(f, string(F.Properties.VariableNames))
        switch f
            case "dx", F.dx = F.dx;
            case "dz", F.dz = F.dz;
            case "f0", F.f0 = F.f0;
            case "M", F.M = F.M;
            otherwise
                error('Required baseline predictor missing in table: %s', f);
        end
    end
end
end

%% Training and model application

function MODELS = train_models(T, train_mask, base_features, CFG)
Fbase = T(:, cellstr(base_features));
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallelTraining);
fprintf('Training composition models with %d primitive predictors.\n', numel(base_features));
MIX = struct();
MIX.base_features = base_features;
MIX.purity = fitrensemble(Fbase(train_mask,:), T.patch_purity(train_mask), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
MIX.mixed = fitcensemble(Fbase(train_mask,:), logical(T.is_mixed(train_mask)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
MIX.strong_mixed = fitcensemble(Fbase(train_mask,:), logical(T.is_strong_mixed(train_mask)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);

P = predict_composition(MIX, T, base_features);
Ttrain = T;
Ttrain.predicted_patch_purity = P.predicted_patch_purity;
Ttrain.p_mixed = P.p_mixed;
Ttrain.p_strong_mixed = P.p_strong_mixed;

Q = struct();
Q.spectrum_only = fit_q_regressor(Ttrain, train_mask, base_features, T.q_oracle, CFG);
Q.spectrum_plus_composition = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "predicted_patch_purity"; "p_mixed"; "p_strong_mixed"], ...
    T.q_oracle, CFG);

MODELS = struct('composition', MIX, 'q', Q, ...
    'model_names', ["q_spectrum_only","q_spectrum_plus_composition"]);
end

function M = fit_q_regressor(T, train_mask, features, y, CFG)
valid = train_mask & all(isfinite(table2array(T(:, cellstr(features)))),2) & isfinite(y);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallelTraining);
M = struct();
M.features = features(:);
M.model = fitrensemble(T(valid, cellstr(features)), y(valid), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template, ...
    'Options', opts);
M.n_train = sum(valid);
end

function T_out = apply_models(T, MODELS, base_features, dataset_name)
P = predict_composition(MODELS.composition, T, base_features);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;
parts = cell(numel(MODELS.model_names),1);
for mi = 1:numel(MODELS.model_names)
    name = MODELS.model_names(mi);
    keep = intersect(string(T_aug.Properties.VariableNames), ...
        ["dataset","condition_key","case_id","case_family","field_regime", ...
        "field_regime_ood","f0","M","dx","dz","map_iz","map_ix","cx","cz", ...
        "x_center_m","z_center_m","true_SWS","k_true","patch_purity", ...
        "purity_bin","is_mixed","is_strong_mixed","distance_to_boundary_mm", ...
        "distance_bin","q_oracle","q_theory_prior","sws_theory", ...
        "is_train_row","is_heldout_row","predicted_patch_purity", ...
        "p_mixed","p_strong_mixed"], 'stable');
    R = T_aug(:, cellstr(keep));
    R.model_name = repmat(name, height(R), 1);
    R.dataset = repmat(dataset_name, height(R), 1);
    switch name
        case "theory_discrete"
            q_pred = T_aug.q_theory_prior;
            sws_pred = T_aug.sws_theory;
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "delta_q_theory_composition"
            delta = predict_q_model(MODELS.q.delta_q_theory_composition, T_aug);
            q_pred = clamp01(T_aug.q_theory_prior + delta);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "delta_logk_theory_composition"
            delta = predict_q_model(MODELS.q.delta_logk_theory_composition, T_aug);
            k_theory = 2*pi*T_aug.f0 ./ T_aug.sws_theory;
            k_pred = k_theory .* exp(delta);
            sws_pred = 2*pi*T_aug.f0 ./ k_pred;
            q_pred = arrayfun(@(i) invert_mapping_to_q(T_aug.req_mapping{i}, k_pred(i)), ...
                (1:height(T_aug))');
        otherwise
            error('Unknown model: %s', name);
    end
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T_aug.q_oracle;
    R.sws_signed_error_pct = 100*(sws_pred - T_aug.true_SWS)./T_aug.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    R.roi_label = assign_roi_labels(R);
    parts{mi} = R;
end
T_out = vertcat(parts{:});
end

function roi = assign_roi_labels(T)
% Evaluation-only ROI labels. These labels use true material and distance
% information only after prediction, never as model predictors.
roi = repmat("other", height(T), 1);
fam = string(T.case_family);
cs = T.true_SWS;
has_dist = ismember("distance_to_boundary_mm", string(T.Properties.VariableNames));
if has_dist
    d = T.distance_to_boundary_mm;
else
    d = nan(height(T),1);
end

is_hom = fam == "homogeneous";
if any(is_hom)
    x0 = median(T.x_center_m(is_hom), 'omitnan');
    z0 = median(T.z_center_m(is_hom), 'omitnan');
    center = is_hom & abs(T.x_center_m - x0) <= 4e-3 & ...
        abs(T.z_center_m - z0) <= 4e-3;
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

function y = predict_q_model(M, T)
X = T(:, cellstr(M.features));
y = clamp01(predict(M.model, X));
end

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

function T = composition_metrics(T)
if ismember("model_name", string(T.Properties.VariableNames))
    T = T(T.model_name == T.model_name(1), :);
end
y = logical(T.is_mixed);
ys = logical(T.is_strong_mixed);
T = table(height(T), rmse(T.patch_purity, T.predicted_patch_purity), ...
    mean(abs(T.patch_purity - T.predicted_patch_purity), 'omitnan'), ...
    auc_binary(y, T.p_mixed, "roc"), auc_binary(y, T.p_mixed, "pr"), ...
    auc_binary(ys, T.p_strong_mixed, "roc"), auc_binary(ys, T.p_strong_mixed, "pr"), ...
    'VariableNames', {'N','purity_RMSE','purity_MAE','mixed_ROC_AUC', ...
    'mixed_PR_AUC','strong_mixed_ROC_AUC','strong_mixed_PR_AUC'});
end

%% Metrics and plots

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G == gi,:);
    rows{gi} = table(height(X), ...
        mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), ...
        mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), ...
        mean(abs(X.q_error),'omitnan'), ...
        mean(X.predicted_patch_purity,'omitnan'), ...
        mean(X.p_mixed,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','median_signed_error_pct', ...
        'underestimate_pct','high_error10_pct','high_error20_pct', ...
        'mean_abs_q_error','mean_predicted_patch_purity','mean_p_mixed'});
end
S = [groups vertcat(rows{:})];
end

function T = top_failure_conditions(T_condition_summary, n)
T = sortrows(T_condition_summary, 'MAPE_pct', 'descend');
T = T(1:min(n,height(T)),:);
end

function plot_condition_maps(T, OUT)
key = T.condition_key(1);
out_dir = fullfile(OUT.map_dir, sanitize(T.case_id(1)), ...
    sanitize(T.field_regime_ood(1)), "f" + string(T.f0(1)), "M" + string(T.M(1)));
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
out_file = fullfile(out_dir, "baseline_minimal_v1_map__" + sanitize(key) + ".png");
[true_map,nz,nx] = rows_to_grid(T, T.true_SWS);
pred_map = rows_to_grid(T, T.sws_pred, nz, nx);
err_map = rows_to_grid(T, T.sws_abs_error_pct, nz, nx);
signed_map = rows_to_grid(T, T.sws_signed_error_pct, nz, nx);
q_map = rows_to_grid(T, T.q_pred, nz, nx);
qerr_map = rows_to_grid(T, T.q_error, nz, nx);
pur_map = rows_to_grid(T, T.patch_purity, nz, nx);
ppur_map = rows_to_grid(T, T.predicted_patch_purity, nz, nx);
pmix_map = rows_to_grid(T, T.p_mixed, nz, nx);
dist_map = rows_to_grid(T, T.distance_to_boundary_mm, nz, nx);
theory_map = rows_to_grid(T, T.sws_theory, nz, nx);
theory_err = abs(rows_to_grid(T, T.sws_theory_error_pct, nz, nx));

fig = figure('Color','w','Units','centimeters','Position',[1 1 32 22]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), true_map, 'True SWS');
plot_map(nexttile(tl), pred_map, 'q spectrum SWS');
plot_map(nexttile(tl), err_map, 'Abs error %');
plot_map(nexttile(tl), signed_map, 'Signed error %');
plot_map(nexttile(tl), theory_map, 'Theory SWS');
plot_map(nexttile(tl), theory_err, 'Theory abs error %');
plot_map(nexttile(tl), q_map, 'q pred');
plot_map(nexttile(tl), qerr_map, 'q pred - q oracle');
plot_map(nexttile(tl), pur_map, 'True patch purity');
plot_map(nexttile(tl), ppur_map, 'Predicted purity');
plot_map(nexttile(tl), pmix_map, 'Predicted p(mixed)');
plot_map(nexttile(tl), dist_map, 'Distance to boundary (mm)');
title(tl, key, 'Interpreter','none');
export_fig(fig, out_file);
end

function tf = should_plot_representative(CFG, C, R, f0, M)
tf = CFG.QuickMode || (M == 2 && ismember(f0, [300 500]) && ...
    ismember(R.regime_id, ["directional_2D_angle0","diffuse_3D_seed1"]) && ...
    ismember(C.case_family, ["homogeneous","bilayer","inclusion","ellipse", ...
    "three_material","smooth_transition"]));
end

function plot_summary_bars(T_overall, T_family, T_purity, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 30 10]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
T_overall = sortrows(T_overall, 'MAPE_pct', 'ascend');
barh(ax, categorical(T_overall.model_name), T_overall.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Held-out model ranking','FontWeight','normal', ...
    'Interpreter','none'); grid(ax,'on');

ax = nexttile(tl);
plot_group_heatmap(ax, T_family, "case_family", "MAPE_pct", 'MAPE by family');

ax = nexttile(tl);
plot_group_heatmap(ax, T_purity, "purity_bin", "MAPE_pct", 'MAPE by purity');
export_fig(fig, fullfile(OUT.figure_dir,'baseline_minimal_v1_summary_bars.png'));
end

function plot_heatmap_table(T, col, OUT, name)
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 9]);
ax = axes(fig);
plot_group_heatmap(ax, T, col, "MAPE_pct", "MAPE by " + string(col));
export_fig(fig, fullfile(OUT.figure_dir,name));
end

function plot_error_vs_predicted_purity(T, OUT)
keep = ismember(T.model_name, ["q_spectrum_only","q_spectrum_plus_composition"]);
T = T(keep,:);
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 9]);
scatter(T.predicted_patch_purity, T.sws_abs_error_pct, 6, T.p_mixed, ...
    'filled', 'MarkerFaceAlpha', 0.15);
xlabel('Predicted patch purity'); ylabel('Abs SWS error (%)');
title('Held-out error vs predicted composition','FontWeight','normal', ...
    'Interpreter','none');
grid on; colorbar;
export_fig(fig, fullfile(OUT.figure_dir,'baseline_minimal_v1_error_vs_predicted_purity.png'));
end

function plot_composition_diagnostics(T, OUT)
models = unique(T.model_name, 'stable');
T0 = T(T.model_name == models(1),:);
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 8]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
scatter(ax, T0.patch_purity, T0.predicted_patch_purity, 5, T0.p_mixed, ...
    'filled', 'MarkerFaceAlpha', 0.18);
hold(ax,'on'); plot(ax,[0 1],[0 1],'k-','LineWidth',0.8);
xlabel(ax,'True patch purity'); ylabel(ax,'Predicted patch purity');
title(ax,'Composition regression','Interpreter','none','FontWeight','normal');
grid(ax,'on'); colorbar(ax);
ax = nexttile(tl);
edges = linspace(0,1,11);
[~,~,bin] = histcounts(T0.predicted_patch_purity, edges);
mape = accumarray(max(bin,1), T0.sws_abs_error_pct, [numel(edges)-1 1], ...
    @(x) mean(x,'omitnan'), NaN);
plot(ax, edges(1:end-1)+diff(edges)/2, mape, '-o');
xlabel(ax,'Predicted purity bin'); ylabel(ax,'Abs SWS error (%)');
title(ax,'Error vs predicted purity','Interpreter','none','FontWeight','normal');
grid(ax,'on');
ax = nexttile(tl);
histogram(ax, T0.predicted_patch_purity, edges);
xlabel(ax,'Predicted patch purity'); ylabel(ax,'N');
title(ax,'Predicted purity distribution','Interpreter','none','FontWeight','normal');
grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'baseline_minimal_v1_composition_diagnostics.png'));
end

function plot_model_ranking(T_overall, OUT)
T_overall = sortrows(T_overall, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 10]);
yyaxis left
bar(categorical(T_overall.model_name), T_overall.MAPE_pct);
ylabel('MAPE (%)');
yyaxis right
plot(categorical(T_overall.model_name), T_overall.high_error20_pct, 'ko-', ...
    'LineWidth',1.0,'MarkerFaceColor','k');
ylabel('High-error >20% (%)');
xtickangle(35); grid on;
title('Baseline minimal held-out model ranking','Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir,'baseline_minimal_v1_model_ranking.png'));
end

function plot_training_maps(T, OUT, CFG)
preferred = ["q_spectrum_only","q_spectrum_plus_composition"];
keys = unique(T.condition_key, 'stable');
if ~CFG.QuickMode
    keep = contains(keys, ["inclusion","bilayer","three_material","homogeneous_cs4p0"]);
    keys = keys(keep);
end
for ki = 1:numel(keys)
    Xall = T(T.condition_key == keys(ki),:);
    if isempty(Xall), continue; end
    X0 = Xall(Xall.model_name == "q_spectrum_only",:);
    if isempty(X0), X0 = Xall(Xall.model_name == Xall.model_name(1),:); end
    [true_map,nz,nx] = rows_to_diagnostic_grid(X0, X0.true_SWS, CFG, "nearest");
    pur_map = rows_to_diagnostic_grid(X0, X0.patch_purity, CFG, "nearest", nz, nx);
    ppur_map = rows_to_diagnostic_grid(X0, X0.predicted_patch_purity, CFG, "natural", nz, nx);
    pmix_map = rows_to_diagnostic_grid(X0, X0.p_mixed, CFG, "natural", nz, nx);
    dist_map = rows_to_diagnostic_grid(X0, X0.distance_to_boundary_mm, CFG, "natural", nz, nx);

    fig = figure('Color','w','Units','centimeters','Position',[1 1 32 20]);
    tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS');
    plot_map(nexttile(tl), pur_map, 'True patch purity');
    plot_map(nexttile(tl), ppur_map, 'Predicted purity');
    plot_map(nexttile(tl), pmix_map, 'Predicted p(mixed)');

    for mi = 1:numel(preferred)
        M = Xall(Xall.model_name == preferred(mi),:);
        if isempty(M)
            plot_map(nexttile(tl), nan(nz,nx), preferred(mi) + ' SWS');
            continue;
        end
        plot_map(nexttile(tl), rows_to_diagnostic_grid(M, M.sws_pred, CFG, "natural", nz, nx), ...
            preferred(mi) + ' SWS');
    end
    for mi = 1:numel(preferred)
        M = Xall(Xall.model_name == preferred(mi),:);
        if isempty(M)
            plot_map(nexttile(tl), nan(nz,nx), preferred(mi) + ' abs err');
            continue;
        end
        plot_map(nexttile(tl), rows_to_diagnostic_grid(M, M.sws_abs_error_pct, CFG, "natural", nz, nx), ...
            preferred(mi) + ' abs err %');
    end
    for mi = 1:numel(preferred)
        M = Xall(Xall.model_name == preferred(mi),:);
        if isempty(M)
            plot_map(nexttile(tl), nan(nz,nx), preferred(mi) + ' q err');
            continue;
        end
        plot_map(nexttile(tl), rows_to_diagnostic_grid(M, M.q_error, CFG, "natural", nz, nx), ...
            preferred(mi) + ' q error');
    end
    plot_map(nexttile(tl), dist_map, 'Distance to boundary (mm)');
    title(tl, keys(ki), 'Interpreter','none');

    out_dir = fullfile(OUT.map_dir, sanitize(X0.case_id(1)), ...
        sanitize(X0.field_regime_ood(1)), "f" + string(X0.f0(1)), ...
        "M" + string(X0.M(1)));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    out_file = fullfile(out_dir, "baseline_minimal_v1_map__" + sanitize(keys(ki)) + ".png");
    export_fig(fig, out_file);
end
end

function plot_group_heatmap(ax, T, col, value_col, ttl)
models = unique(string(T.model_name), 'stable');
groups = unique(string(T.(col)), 'stable');
Z = nan(numel(models), numel(groups));
for i = 1:numel(models)
    for j = 1:numel(groups)
        idx = string(T.model_name) == models(i) & string(T.(col)) == groups(j);
        if any(idx), Z(i,j) = mean(T.(value_col)(idx), 'omitnan'); end
    end
end
imagesc(ax, Z, 'AlphaData', isfinite(Z));
colorbar(ax); grid(ax,'off');
set(ax,'XTick',1:numel(groups),'XTickLabel',groups, ...
    'YTick',1:numel(models),'YTickLabel',models, ...
    'TickLabelInterpreter','none');
xtickangle(ax,35);
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
end

function plot_worst_condition_gallery(T, T_fail, OUT)
n = min(6, height(T_fail));
for i = 1:n
    X = T(T.condition_key == T_fail.condition_key(i), :);
    if isempty(X), continue; end
    plot_condition_maps(X, struct('map_dir',OUT.fail_dir));
end
end

function print_interpretation(T_overall, T_family, T_purity, T_fail, OUT)
fprintf('\n================ Baseline minimal v1 held-out summary ================\n');
disp(T_overall(:, {'model_name','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct','mean_p_mixed'}));
fprintf('\nWorst held-out conditions:\n');
disp(T_fail(:, {'model_name','condition_key','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct','mean_p_mixed'}));
[~, bi] = min(T_overall.MAPE_pct);
best_model = T_overall.model_name(bi);
fprintf('\nBest held-out model: %s, MAPE %.2f%%, high-error20 %.2f%%.\n', ...
    best_model, T_overall.MAPE_pct(bi), T_overall.high_error20_pct(bi));
BF = T_family(T_family.model_name == best_model,:);
if ~isempty(BF)
    [~, wi] = max(BF.MAPE_pct);
    fprintf('Worst family for best model: %s, MAPE %.2f%%.\n', ...
        BF.case_family(wi), BF.MAPE_pct(wi));
end
BP = T_purity(T_purity.model_name == best_model,:);
if ~isempty(BP)
    [~, pi] = max(BP.MAPE_pct);
    fprintf('Worst purity bin for best model: %s, MAPE %.2f%%.\n', ...
        BP.purity_bin(pi), BP.MAPE_pct(pi));
end
fprintf('==========================================================\n');
fprintf('Figures: %s\n', OUT.figure_dir);
end

%% Numerical helpers

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('BaselineMinimal:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
end
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

function y = q_to_sws(mappings, q, f0)
y = nan(numel(q),1);
for i = 1:numel(q)
    if isscalar(f0), fi = f0; else, fi = f0(i); end
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), fi);
end
end

function q = invert_mapping_to_q(mapping, k_target)
k = double(mapping.k_cent(:)); E = double(mapping.Ecum(:));
valid = isfinite(k) & isfinite(E);
k = k(valid); E = E(valid);
if numel(k) < 2 || ~isfinite(k_target), q = NaN; return; end
[ku, ia] = unique(k, 'stable'); Eu = E(ia);
q = interp1(ku, Eu, k_target, 'linear', 'extrap');
q = clamp01(q);
end

function q = theory_q_for_regime(cfg, M, regime)
feat = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', 3, 'gamma_win', 1, 'pad_factor', 1);
switch string(regime)
    case "directional_2D", type = "SingleWave";
    case "diffuse_2D", type = "Diffuse2D";
    case "partial_3D", type = "Partial3D";
    otherwise, type = "Diffuse3D";
end
try
    if type == "Partial3D"
        q = 0.5*(theory_one(cfg,feat,"Diffuse2D") + theory_one(cfg,feat,"Diffuse3D"));
    else
        q = theory_one(cfg,feat,type);
    end
catch
    q = 0.5;
end
q = clamp01(q);
end

function q = theory_one(cfg, feat, field_type)
out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    cfg.dx, cfg.dz, cfg.f0, feat.cs_guess_used, ...
    'M', feat.M, 'Gamma', feat.gamma_win, 'PadFactor', feat.pad_factor, ...
    'Nbins', 'auto', 'SmoothSigma', 1, 'TheoryMode', 'S2D', ...
    'FieldType', field_type, 'Plot', false);
q = out.q_th;
end

function p = dominant_fraction(patch)
v = round(double(patch(:))*1e5)/1e5;
v = v(isfinite(v));
if isempty(v), p = NaN; return; end
[~,~,ic] = unique(v);
counts = accumarray(ic,1);
p = max(counts)/numel(v);
end

function labels = purity_bin(p)
labels = repmat("strongly_mixed", size(p));
labels(p >= 0.75) = "moderately_mixed";
labels(p >= 0.95) = "near_pure";
labels(p >= 0.99) = "pure";
labels(~isfinite(p)) = "unknown";
end

function Dmm = distance_to_material_boundary(cs_map, dx)
Q = round(double(cs_map)*1e5);
B = false(size(Q));
B(:,1:end-1) = B(:,1:end-1) | Q(:,1:end-1) ~= Q(:,2:end);
B(:,2:end) = B(:,2:end) | Q(:,1:end-1) ~= Q(:,2:end);
B(1:end-1,:) = B(1:end-1,:) | Q(1:end-1,:) ~= Q(2:end,:);
B(2:end,:) = B(2:end,:) | Q(1:end-1,:) ~= Q(2:end,:);
if ~any(B(:))
    Dmm = nan(size(Q));
else
    Dmm = bwdist(B) * dx * 1e3;
end
end

function labels = distance_bin(d, edges)
labels = strings(size(d));
for i = 1:numel(edges)-1
    idx = d >= edges(i) & d < edges(i+1);
    labels(idx) = sprintf('%.1f-%.1fmm', edges(i), edges(i+1));
end
labels(d >= edges(end-1)) = sprintf('>%.1fmm', edges(end-1));
labels(~isfinite(d)) = "no_interface";
end

function y = clamp01(x), y = min(max(x,0),1); end

function e = rmse(y, yhat)
e = sqrt(mean((y - yhat).^2, 'omitnan'));
end

function a = auc_binary(y, score, mode)
y = logical(y(:)); score = score(:);
valid = isfinite(score);
y = y(valid); score = score(valid);
if numel(unique(y)) < 2
    a = NaN; return;
end
try
    switch string(mode)
        case "pr"
            [~,~,~,a] = perfcurve(y, score, true, 'XCrit','reca', 'YCrit','prec');
        otherwise
            [~,~,~,a] = perfcurve(y, score, true);
    end
catch
    a = NaN;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function sig = cache_signature(CFG)
C = CFG;
if isfield(C,'ValidateOnly'), C = rmfield(C,'ValidateOnly'); end
sig = string(jsonencode(C));
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3, nz = max(T.map_iz); nx = max(T.map_ix); end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function [Z,nz,nx] = rows_to_diagnostic_grid(T, values, CFG, method, nz, nx)
if nargin < 5
    scale = max(1, round(CFG.DiagnosticMapInterpScale));
    nz = max(T.map_iz) * scale;
    nx = max(T.map_ix) * scale;
end
values = double(values(:));
x = double(T.x_center_m(:));
z = double(T.z_center_m(:));
valid = isfinite(x) & isfinite(z) & isfinite(values);
if nnz(valid) < 3
    Z = nan(nz,nx);
    return;
end
xq = linspace(min(x(valid)), max(x(valid)), nx);
zq = linspace(min(z(valid)), max(z(valid)), nz);
[Xq,Zq] = meshgrid(xq, zq);
try
    F = scatteredInterpolant(x(valid), z(valid), values(valid), char(method), 'none');
catch
    F = scatteredInterpolant(x(valid), z(valid), values(valid), 'nearest', 'none');
end
Z = F(Xq, Zq);
if any(~isfinite(Z(:)))
    Fallback = scatteredInterpolant(x(valid), z(valid), values(valid), 'nearest', 'nearest');
    Z(~isfinite(Z)) = Fallback(Xq(~isfinite(Z)), Zq(~isfinite(Z)));
end
end

function plot_map(ax, Z, ttl)
imagesc(ax, Z, 'AlphaData', isfinite(Z));
axis(ax,'image'); ax.Color = [0.94 0.94 0.94];
colorbar(ax); title(ax, ttl, 'Interpreter','none', 'FontWeight','normal');
set(ax,'XTick',[],'YTick',[]);
end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
end

function write_config_json(CFG, path)
fid = fopen(path,'w'); assert(fid > 0);
fprintf(fid, '%s', jsonencode(CFG, PrettyPrint=true));
fclose(fid);
end

function CFG = setup_parallel_if_requested(CFG)
if ~(CFG.UseParfor || CFG.UseParallelTraining)
    return;
end
if exist('gcp','file') ~= 2 || exist('parpool','file') ~= 2
    warning('BaselineMinimal:NoParallelToolbox', ...
        'Parallel toolbox not available. Continuing without parallel execution.');
    CFG.UseParfor = false;
    CFG.UseParallelTraining = false;
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool(char(CFG.ParallelPoolType));
        catch
            parpool;
        end
    end
catch ME
    warning('BaselineMinimal:ParallelPoolFailed', ...
        'Could not start parallel pool (%s). Continuing serial.', ME.message);
    CFG.UseParfor = false;
    CFG.UseParallelTraining = false;
end
end

function tf = env_true(name, default)
raw = strtrim(lower(string(getenv(name))));
if raw == "", tf = logical(default); else, tf = ismember(raw, ["1","true","yes","on"]); end
end

function tf = env_true_any(names, default)
tf = logical(default);
for name = string(names)
    raw = strtrim(lower(string(getenv(name))));
    if raw ~= ""
        tf = ismember(raw, ["1","true","yes","on"]);
        return;
    end
end
end

function val = env_number(name, default)
raw = strtrim(string(getenv(name)));
if raw == ""
    val = default;
else
    val = str2double(raw);
    if ~isfinite(val), val = default; end
end
end

function val = env_number_any(names, default)
val = default;
for name = string(names)
    raw = strtrim(string(getenv(name)));
    if raw ~= ""
        x = str2double(raw);
        if isfinite(x), val = x; end
        return;
    end
end
end

function vals = env_number_list(name, default)
raw = strtrim(string(getenv(name)));
if raw == ""
    vals = default;
    return;
end
raw = erase(raw, ["[", "]"]);
tokens = regexp(char(raw), '[,;\s]+', 'split');
tokens = tokens(~cellfun(@isempty, tokens));
vals = str2double(tokens);
vals = vals(isfinite(vals));
if isempty(vals)
    vals = default;
end
vals = vals(:).';
end

function vals = env_number_list_any(names, default)
vals = default;
for name = string(names)
    raw = strtrim(string(getenv(name)));
    if raw ~= ""
        vals = env_number_list(char(name), default);
        return;
    end
end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function export_fig(fig, path)
axs = findall(fig, 'Type', 'Axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end
drawnow;
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
close(fig);
end
