%% analyze_test_68_test66_estimability_mask.m
% Test 68: operational estimability mask for Test 66 Eikonal transfer.
%
% This analysis does not retrain q models and does not correct SWS. It uses
% frozen Test66 predictions from q_spectrum_plus_composition and trains a
% separate, operational risk layer that predicts whether a patch is likely
% non-estimable / high-error. Truth, patch purity, distance to interface, and
% ROI labels are used only as labels or diagnostic grouping variables.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST68_MODE = quick | full
%   ADAPTIVE_REQ_TEST68_TEST66_SOURCE = full_a | quick | /path/to/Test66/run
%   ADAPTIVE_REQ_TEST68_MAX_TRAIN_ROWS = integer
%   ADAPTIVE_REQ_TEST68_MAX_EVAL_ROWS = integer
%   ADAPTIVE_REQ_TEST68_MAX_PATCH_ROWS_TO_SAVE = integer
%   ADAPTIVE_REQ_TEST68_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST68_USE_PARALLEL = true | false

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test68_configuration.json'));

fprintf('\nTest 68: operational estimability mask for Test66\n');
fprintf('Mode: %s | Test66 source: %s\n', CFG.Mode, CFG.Test66RunDir);
fprintf('Frozen q model evaluated: %s\n', CFG.BaseModel);
fprintf('No q, composition, correction, or confidence model from earlier tests is modified.\n');

T = load_test66_rows(CFG);
assert(~isempty(T), 'No Test66 rows matched the Test68 filters.');
fprintf('Loaded %d patch rows after Test68 filters.\n', height(T));

T = add_operational_features(T, CFG);
assert_no_forbidden_predictors([CFG.NumericFeatureVars; CFG.CategoricalFeatureVars]);

SPLITS = build_splits(T, CFG);
fprintf('Risk-model folds to run: %d\n', numel(SPLITS));

all_pred = {};
all_metrics = {};
all_thresh = {};
all_cov = {};
all_cal = {};
MODELS = struct();

for si = 1:numel(SPLITS)
    FOLD = SPLITS(si);
    fprintf('\n[%d/%d] %s | %s\n', si, numel(SPLITS), FOLD.split_name, FOLD.fold_name);
    [P, MET, TH, COV, CAL, MDL] = run_one_fold(T, FOLD, CFG);
    all_pred{end+1,1} = P; %#ok<AGROW>
    all_metrics{end+1,1} = MET; %#ok<AGROW>
    all_thresh{end+1,1} = TH; %#ok<AGROW>
    all_cov{end+1,1} = COV; %#ok<AGROW>
    all_cal{end+1,1} = CAL; %#ok<AGROW>
    MODELS(si).split_name = FOLD.split_name; %#ok<SAGROW>
    MODELS(si).fold_name = FOLD.fold_name;
    MODELS(si).detectors = MDL;
end

PRED = vertcat(all_pred{:});
METRICS = vertcat(all_metrics{:});
THRESH = vertcat(all_thresh{:});
COVERAGE = vertcat(all_cov{:});
CALIB = vertcat(all_cal{:});

PRED_TO_SAVE = cap_patch_table_for_saving(PRED, CFG);
PRED_MAIN = select_main_predictions(PRED, CFG);
SUM = make_grouped_summaries(PRED_MAIN, CFG);

writetable(PRED_TO_SAVE, fullfile(OUT.table_dir, 'test68_patch_level_risk_predictions.csv'));
writetable(METRICS, fullfile(OUT.table_dir, 'test68_detector_metrics.csv'));
writetable(THRESH, fullfile(OUT.table_dir, 'test68_threshold_summary.csv'));
writetable(COVERAGE, fullfile(OUT.table_dir, 'test68_accuracy_coverage.csv'));
writetable(CALIB, fullfile(OUT.table_dir, 'test68_reliability_calibration.csv'));
writetable(SUM.geometry, fullfile(OUT.table_dir, 'test68_summary_by_geometry.csv'));
writetable(SUM.frequency, fullfile(OUT.table_dir, 'test68_summary_by_frequency.csv'));
writetable(SUM.realism, fullfile(OUT.table_dir, 'test68_summary_by_realism_level.csv'));
writetable(SUM.field_regime, fullfile(OUT.table_dir, 'test68_summary_by_field_regime.csv'));
writetable(SUM.roi, fullfile(OUT.table_dir, 'test68_summary_by_roi.csv'));
writetable(SUM.purity, fullfile(OUT.table_dir, 'test68_summary_by_true_patch_purity.csv'));
writetable(SUM.distance, fullfile(OUT.table_dir, 'test68_summary_by_distance_over_window.csv'));

fprintf('\nWriting Test68 figures...\n');
safe_plot(@() plot_auc_summary(METRICS, OUT), 'AUC summary');
safe_plot(@() plot_calibration(CALIB, OUT), 'calibration');
safe_plot(@() plot_accuracy_coverage(COVERAGE, OUT), 'accuracy coverage');
safe_plot(@() plot_threshold_tradeoff(THRESH, OUT), 'threshold tradeoff');
safe_plot(@() plot_diagnostic_risk_vs_interface(PRED_MAIN, OUT), 'risk vs interface diagnostics');
safe_plot(@() plot_risk_by_frequency_regime(SUM.frequency, SUM.field_regime, OUT), 'risk by frequency/regime');
safe_plot(@() plot_region_heatmaps(SUM.roi, SUM.geometry, OUT), 'region heatmaps');
if CFG.SaveAllMaps
    safe_plot(@() plot_representative_maps(PRED_MAIN, OUT, CFG), 'representative maps');
end

write_readme(PRED_MAIN, METRICS, THRESH, SUM, OUT, CFG);
save(fullfile(OUT.data_dir, 'test68_estimability_mask_compact.mat'), ...
    'CFG','METRICS','THRESH','COVERAGE','CALIB','SUM','MODELS','-v7.3');

print_console_summary(METRICS, THRESH, OUT);
fprintf('\nTest 68 complete.\nTables: %s\nFigures: %s\nREADME: %s\n', ...
    OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST68_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST68_MODE must be quick or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.BaseModel = "q_spectrum_plus_composition";
CFG.CsGuess = 3.0;
CFG.TrainFraction = 0.70;
CFG.GroupedSeeds = ternary(CFG.QuickMode, [6801 6802], 6801:6805);
CFG.RiskThresholds = [0.2 0.5 0.8];
CFG.CoverageGrid = [0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00];
CFG.CalibrationEdges = 0:0.1:1;
CFG.TreeLearners = env_number('ADAPTIVE_REQ_TEST68_TREE_LEARNERS', ternary(CFG.QuickMode, 80, 160));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST68_MIN_LEAF_SIZE', ternary(CFG.QuickMode, 12, 18));
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST68_USE_PARALLEL', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST68_SAVE_ALL_MAPS', true);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST68_MAX_MAP_CONDITIONS', ternary(CFG.QuickMode, 8, 36));
CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST68_MAX_TRAIN_ROWS', ternary(CFG.QuickMode, 60000, 400000));
CFG.MaxEvalRows = env_number('ADAPTIVE_REQ_TEST68_MAX_EVAL_ROWS', ternary(CFG.QuickMode, 160000, 600000));
CFG.MaxPatchRowsToSave = env_number('ADAPTIVE_REQ_TEST68_MAX_PATCH_ROWS_TO_SAVE', ternary(CFG.QuickMode, 50000, 200000));

if CFG.QuickMode
    CFG.Geometries = ["homogeneous_cs2","homogeneous_cs4","inclusion_2_4","bilayer_inclusion_2_3_4"];
    CFG.Frequencies = [300 500];
    CFG.RealismLevels = ["clean","readout_medium"];
    CFG.MValues = 2;
    CFG.FieldRegimes = ["single_source_lateral","diffuse_like_8src"];
else
    CFG.Geometries = strings(0,1);
    CFG.Frequencies = [];
    CFG.RealismLevels = strings(0,1);
    CFG.MValues = [];
    CFG.FieldRegimes = strings(0,1);
end

src = env_string('ADAPTIVE_REQ_TEST68_TEST66_SOURCE', "full_a");
base = fullfile(root_dir, 'outputs', 'test_66_eikonal_realistic_transfer_validation');
if isfolder(src)
    CFG.Test66RunDir = char(src);
else
    CFG.Test66RunDir = fullfile(base, char(src));
end
if ~isfolder(CFG.Test66RunDir) && strcmpi(char(src), 'full')
    CFG.Test66RunDir = fullfile(base, 'full_a');
end
CFG.PatchCsv = fullfile(CFG.Test66RunDir, 'tables', 'test66_patch_level_results.csv');
assert(exist(CFG.PatchCsv,'file') == 2, 'Missing Test66 patch-level table: %s', CFG.PatchCsv);

CFG.NumericFeatureVars = ["q_pred";"SWS_pred";"k_pred";"log_sws_pred"; ...
    "M_eff_pred";"q_theory_prior";"sws_theory";"abs_q_minus_theory"; ...
    "rel_sws_minus_theory";"predicted_patch_purity";"p_mixed"; ...
    "p_strong_mixed";"f0";"M";"dx";"dz";"REQ_StepX";"REQ_StepZ"; ...
    "TargetStepM";"local_amplitude";"tracking_snr_db"; ...
    "acoustic_readout_amplitude_factor";"source_to_patch_distance_m"; ...
    "depth_m";"grad_q_pred";"grad_sws_pred";"local_std_q_pred"; ...
    "local_std_sws_pred"];
CFG.CategoricalFeatureVars = "field_regime";
end

function OUT = make_output_dirs(root_dir, CFG)
OUT = struct();
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_68_test66_estimability_mask');
if CFG.QuickMode
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if exist(dirs(i), 'dir') ~= 7, mkdir(dirs(i)); end
end
end

%% Data loading and feature construction

function T = load_test66_rows(CFG)
needed = ["map_iz","map_ix","x_center_m","z_center_m","condition_key", ...
    "geometry","geometry_family","case_id","case_family","realism_level", ...
    "field_regime","field_regime_ood","source_seed","noise_seed","f0","M", ...
    "REQ_M","REQ_StepX","REQ_StepZ","TargetStepM","dx","dz", ...
    "true_SWS","k_true","patch_purity","patch_cs_std","patch_cs_range", ...
    "patch_cs_iqr","q_oracle","distance_to_interface_m", ...
    "distance_to_interface_mm","distance_to_interface_over_window_radius", ...
    "depth_m","local_amplitude","tracking_snr_db","snr_proxy_db", ...
    "shear_attenuation_factor","acoustic_readout_amplitude_factor", ...
    "source_to_patch_distance_m","purity_bin","distance_bin", ...
    "distance_over_window_bin","depth_bin","snr_bin","amplitude_bin", ...
    "roi_region","analysis_region","q_theory_prior","sws_theory", ...
    "predicted_patch_purity","p_mixed","p_strong_mixed","model_name", ...
    "q_pred","k_pred","SWS_pred","sws_pred","q_error","k_error", ...
    "sws_signed_error_pct","sws_abs_error_pct","high_error10","high_error20"];

ds = tabularTextDatastore(CFG.PatchCsv, 'Delimiter', ',');
vars = string(ds.VariableNames);
ds.SelectedVariableNames = cellstr(intersect(needed, vars, 'stable'));
ds.ReadSize = 250000;

parts = {};
while hasdata(ds)
    C = read(ds);
    C = normalize_loaded_types(C);
    m = C.model_name == CFG.BaseModel;
    if ~isempty(CFG.Geometries), m = m & ismember(C.geometry, CFG.Geometries); end
    if ~isempty(CFG.Frequencies), m = m & ismember(C.f0, CFG.Frequencies); end
    if ~isempty(CFG.RealismLevels), m = m & ismember(C.realism_level, CFG.RealismLevels); end
    if ~isempty(CFG.MValues), m = m & ismember(round(C.REQ_M), CFG.MValues); end
    if ~isempty(CFG.FieldRegimes), m = m & ismember(C.field_regime, CFG.FieldRegimes); end
    if any(m)
        parts{end+1,1} = C(m,:); %#ok<AGROW>
    end
end

if isempty(parts)
    T = table();
else
    T = vertcat(parts{:});
end
T = normalize_loaded_types(T);
T.M = T.REQ_M;
if ~ismember("analysis_region", string(T.Properties.VariableNames))
    T.analysis_region = T.roi_region;
end
end

function T = normalize_loaded_types(T)
if isempty(T), return; end
names = string(T.Properties.VariableNames);
string_vars = intersect(names, ["condition_key","geometry","geometry_family", ...
    "case_id","case_family","realism_level","field_regime","field_regime_ood", ...
    "purity_bin","distance_bin","distance_over_window_bin","depth_bin", ...
    "snr_bin","amplitude_bin","roi_region","analysis_region","model_name"], 'stable');
for v = string_vars
    T.(v) = string(T.(v));
end
numeric_vars = setdiff(names, string_vars, 'stable');
for v = numeric_vars
    if ~isnumeric(T.(v)) && ~islogical(T.(v))
        T.(v) = str2double(string(T.(v)));
    end
end
if ismember("high_error10", names), T.high_error10 = logical(T.high_error10); end
if ismember("high_error20", names), T.high_error20 = logical(T.high_error20); end
end

function T = add_operational_features(T, CFG)
T.SWS_pred = first_existing_numeric(T, ["SWS_pred","sws_pred"]);
T.q_pred = clamp(T.q_pred, 0, 1);
T.SWS_pred = clamp(T.SWS_pred, 0.2, 20);
T.k_pred = first_existing_numeric(T, ["k_pred"]);
bad_k = ~isfinite(T.k_pred) | T.k_pred <= 0;
T.k_pred(bad_k) = 2*pi*T.f0(bad_k) ./ T.SWS_pred(bad_k);

T.log_sws_pred = log(T.SWS_pred);
T.lambda_pred_m = T.SWS_pred ./ T.f0;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.M_eff_pred = T.window_length_m ./ max(T.lambda_pred_m, eps);
T.q_theory_prior = ensure_numeric_var(T, "q_theory_prior", nan(height(T),1));
T.sws_theory = ensure_numeric_var(T, "sws_theory", nan(height(T),1));
T.abs_q_minus_theory = abs(T.q_pred - T.q_theory_prior);
T.rel_sws_minus_theory = abs(T.SWS_pred - T.sws_theory) ./ max(abs(T.SWS_pred), eps);

T.predicted_patch_purity = clamp(ensure_numeric_var(T, "predicted_patch_purity", ones(height(T),1)), 0, 1);
T.p_mixed = clamp(ensure_numeric_var(T, "p_mixed", zeros(height(T),1)), 0, 1);
T.p_strong_mixed = clamp(ensure_numeric_var(T, "p_strong_mixed", zeros(height(T),1)), 0, 1);
T.tracking_snr_db = ensure_numeric_var(T, "tracking_snr_db", ensure_numeric_var(T, "snr_proxy_db", nan(height(T),1)));
T.acoustic_readout_amplitude_factor = ensure_numeric_var(T, "acoustic_readout_amplitude_factor", nan(height(T),1));
T.local_amplitude = ensure_numeric_var(T, "local_amplitude", nan(height(T),1));
T.source_to_patch_distance_m = ensure_numeric_var(T, "source_to_patch_distance_m", nan(height(T),1));
T.depth_m = ensure_numeric_var(T, "depth_m", T.z_center_m);

T = add_map_neighborhood_features(T);
if ~ismember("high_error10", string(T.Properties.VariableNames))
    T.high_error10 = T.sws_abs_error_pct > 10;
end
if ~ismember("high_error20", string(T.Properties.VariableNames))
    T.high_error20 = T.sws_abs_error_pct > 20;
end

for v = CFG.NumericFeatureVars(:)'
    if ~ismember(v, string(T.Properties.VariableNames))
        T.(v) = nan(height(T),1);
    end
end
end

function T = add_map_neighborhood_features(T)
T.grad_q_pred = nan(height(T),1);
T.grad_sws_pred = nan(height(T),1);
T.local_std_q_pred = nan(height(T),1);
T.local_std_sws_pred = nan(height(T),1);
conds = unique(T.condition_key, 'stable');
for ci = 1:numel(conds)
    idx = find(T.condition_key == conds(ci));
    [Q, rows, cols] = vector_to_map(T.map_iz(idx), T.map_ix(idx), T.q_pred(idx));
    S = vector_to_map(T.map_iz(idx), T.map_ix(idx), T.SWS_pred(idx));
    Gq = gradient_magnitude(Q);
    Gs = gradient_magnitude(S);
    Lq = local_std_nan(Q, 3);
    Ls = local_std_nan(S, 3);
    lin = sub2ind(size(Q), rows, cols);
    T.grad_q_pred(idx) = Gq(lin);
    T.grad_sws_pred(idx) = Gs(lin);
    T.local_std_q_pred(idx) = Lq(lin);
    T.local_std_sws_pred(idx) = Ls(lin);
end
end

function [A, rr, cc] = vector_to_map(r, c, v)
r = double(r); c = double(c);
ur = unique(r, 'stable'); uc = unique(c, 'stable');
[~, rr] = ismember(r, ur);
[~, cc] = ismember(c, uc);
A = nan(numel(ur), numel(uc));
A(sub2ind(size(A), rr, cc)) = v;
end

function G = gradient_magnitude(A)
Af = fillmissing2(A);
[gz, gx] = gradient(Af);
G = sqrt(gx.^2 + gz.^2);
G(~isfinite(A)) = NaN;
end

function S = local_std_nan(A, win)
Af = fillmissing2(A);
k = ones(win, win);
n = conv2(double(isfinite(A)), k, 'same');
s1 = conv2(Af, k, 'same');
s2 = conv2(Af.^2, k, 'same');
mu = s1 ./ max(n, 1);
S = sqrt(max(s2 ./ max(n, 1) - mu.^2, 0));
S(n <= 1) = 0;
S(~isfinite(A)) = NaN;
end

function B = fillmissing2(A)
B = A;
if all(~isfinite(B(:)))
    B(:) = 0;
    return;
end
med = median(B(isfinite(B)), 'omitnan');
B(~isfinite(B)) = med;
end

%% Splits and training

function SPLITS = build_splits(T, CFG)
SPLITS = struct('split_name', {}, 'fold_name', {}, 'train_mask', {}, 'test_mask', {});
keys = unique(T.condition_key, 'stable');
for si = 1:numel(CFG.GroupedSeeds)
    rng(CFG.GroupedSeeds(si));
    perm = keys(randperm(numel(keys)));
    ntrain = max(1, round(CFG.TrainFraction * numel(perm)));
    tr_keys = perm(1:ntrain);
    SPLITS(end+1) = make_fold("grouped_condition_repeated", ... %#ok<AGROW>
        sprintf('seed_%d', CFG.GroupedSeeds(si)), ismember(T.condition_key, tr_keys), ...
        ~ismember(T.condition_key, tr_keys));
end

if any(T.realism_level == "clean") && any(T.realism_level == "readout_medium")
    SPLITS(end+1) = make_fold("realism_transfer", "train_clean_test_readout_medium", ... %#ok<AGROW>
        T.realism_level == "clean", T.realism_level == "readout_medium");
end

if any(T.field_regime == "single_source_lateral") && any(T.field_regime == "diffuse_like_8src")
    SPLITS(end+1) = make_fold("field_transfer", "train_single_test_diffuse8", ... %#ok<AGROW>
        T.field_regime == "single_source_lateral", T.field_regime == "diffuse_like_8src");
end
end

function F = make_fold(split_name, fold_name, train_mask, test_mask)
F = struct();
F.split_name = string(split_name);
F.fold_name = string(fold_name);
F.train_mask = logical(train_mask);
F.test_mask = logical(test_mask);
end

function [P, MET, TH, COV, CAL, MODEL] = run_one_fold(T, FOLD, CFG)
train_idx0 = find(FOLD.train_mask & isfinite(T.sws_abs_error_pct));
test_idx0 = find(FOLD.test_mask & isfinite(T.sws_abs_error_pct));
assert(~isempty(train_idx0) && ~isempty(test_idx0), 'Fold has empty train or test set.');

train_idx = stratified_sample_idx(train_idx0, T.high_error20(train_idx0), CFG.MaxTrainRows, 100 + numel(train_idx0));
test_idx = stratified_sample_idx(test_idx0, T.high_error20(test_idx0), CFG.MaxEvalRows, 200 + numel(test_idx0));

[Xtr, Xte, DESIGN] = design_matrices(T(train_idx,:), T(test_idx,:), CFG);
MODEL = struct();
MODEL.design = DESIGN;
MODEL.features = [CFG.NumericFeatureVars; "cat_field_regime"];
MODEL.detectors = ["logistic_regression","bagged_trees"];

targets = ["high10","high20"];
pred_parts = {};
metric_parts = {};
threshold_parts = {};
coverage_parts = {};
cal_parts = {};
for target = targets
    ytr = target_values(T(train_idx,:), target);
    yte = target_values(T(test_idx,:), target);
    if numel(unique(ytr)) < 2 || numel(unique(yte)) < 2
        warning('Test68:DegenerateFold', 'Skipping %s/%s because target has one class.', FOLD.split_name, target);
        continue;
    end
    weights = class_balance_weights(ytr);
    fprintf('  target %s: train rows %d, eval rows %d, positive train %.2f%%, positive eval %.2f%%\n', ...
        target, numel(ytr), numel(yte), 100*mean(ytr), 100*mean(yte));

    LOG = fitclinear(Xtr, ytr, 'Learner','logistic', 'Regularization','ridge', ...
        'ClassNames',[false true], 'Weights', weights);
    risk_log = predict_positive(LOG, Xte);
    [pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts] = ...
        collect_outputs(pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts, ...
        T(test_idx,:), yte, risk_log, FOLD, target, "logistic_regression", CFG);
    MODEL.(target).logistic_regression = LOG;

    template = templateTree('MinLeafSize', CFG.MinLeafSize);
    opts = statset('UseParallel', CFG.UseParallel);
    BAG = fitcensemble(Xtr, ytr, 'Method','Bag', 'NumLearningCycles', CFG.TreeLearners, ...
        'Learners', template, 'ClassNames', [false true], 'Weights', weights, 'Options', opts);
    risk_bag = predict_positive(BAG, Xte);
    [pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts] = ...
        collect_outputs(pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts, ...
        T(test_idx,:), yte, risk_bag, FOLD, target, "bagged_trees", CFG);
    MODEL.(target).bagged_trees = BAG;
end

P = vertcat(pred_parts{:});
MET = vertcat(metric_parts{:});
TH = vertcat(threshold_parts{:});
COV = vertcat(coverage_parts{:});
CAL = vertcat(cal_parts{:});

if CFG.MaxPatchRowsToSave > 0 && height(P) > CFG.MaxPatchRowsToSave
    P = P(stratified_sample_idx((1:height(P))', P.high_error20, CFG.MaxPatchRowsToSave, 68000), :);
end
end

function [pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts] = collect_outputs( ...
    pred_parts, metric_parts, threshold_parts, coverage_parts, cal_parts, Te, y, risk, FOLD, target, detector, CFG)
P = prediction_table(Te, risk, FOLD, target, detector);
MET = detector_metrics(P, y, risk);
TH = threshold_summary(P, CFG);
COV = accuracy_coverage(P, CFG);
CAL = calibration_summary(P, CFG);
pred_parts{end+1,1} = P;
metric_parts{end+1,1} = MET;
threshold_parts{end+1,1} = TH;
coverage_parts{end+1,1} = COV;
cal_parts{end+1,1} = CAL;
end

function [Xtr, Xte, DESIGN] = design_matrices(Ttr, Tte, CFG)
Xtr_num = zeros(height(Ttr), numel(CFG.NumericFeatureVars));
Xte_num = zeros(height(Tte), numel(CFG.NumericFeatureVars));
med = zeros(1, numel(CFG.NumericFeatureVars));
for i = 1:numel(CFG.NumericFeatureVars)
    v = CFG.NumericFeatureVars(i);
    a = double(Ttr.(v));
    b = double(Tte.(v));
    med(i) = median(a(isfinite(a)), 'omitnan');
    if ~isfinite(med(i)), med(i) = 0; end
    a(~isfinite(a)) = med(i);
    b(~isfinite(b)) = med(i);
    Xtr_num(:,i) = a;
    Xte_num(:,i) = b;
end

mu = mean(Xtr_num, 1, 'omitnan');
sd = std(Xtr_num, 0, 1, 'omitnan');
sd(~isfinite(sd) | sd <= 0) = 1;
Xtr_num = (Xtr_num - mu) ./ sd;
Xte_num = (Xte_num - mu) ./ sd;

cats = unique(string(Ttr.field_regime), 'stable');
Xtr_cat = onehot(string(Ttr.field_regime), cats);
Xte_cat = onehot(string(Tte.field_regime), cats);
Xtr = [Xtr_num Xtr_cat];
Xte = [Xte_num Xte_cat];

DESIGN = struct('numeric_features', CFG.NumericFeatureVars, 'categorical_features', "field_regime", ...
    'field_regime_categories', cats, 'numeric_median', med, 'numeric_mu', mu, 'numeric_sd', sd);
end

function X = onehot(x, cats)
X = zeros(numel(x), numel(cats));
for i = 1:numel(cats)
    X(:,i) = x == cats(i);
end
end

function y = target_values(T, target)
switch string(target)
    case "high10"
        y = logical(T.high_error10);
    otherwise
        y = logical(T.high_error20);
end
end

function w = class_balance_weights(y)
y = logical(y);
p = mean(y);
if p <= 0 || p >= 1
    w = ones(size(y));
else
    w = ones(size(y));
    w(y) = 0.5 / p;
    w(~y) = 0.5 / (1-p);
    w = w / mean(w);
end
end

function risk = predict_positive(model, X)
[~, score] = predict(model, X);
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
risk = score(:,idx);
risk = clamp(risk, 0, 1);
end

function idx = stratified_sample_idx(idx0, y, max_n, seed)
idx0 = idx0(:);
if max_n <= 0 || numel(idx0) <= max_n
    idx = idx0;
    return;
end
rng(seed);
y = logical(y(:));
pos = idx0(y);
neg = idx0(~y);
n_pos = min(numel(pos), max(1, round(0.35*max_n)));
n_neg = max_n - n_pos;
if numel(neg) < n_neg
    n_neg = numel(neg);
    n_pos = min(numel(pos), max_n - n_neg);
end
pos = pos(randperm(numel(pos), n_pos));
neg = neg(randperm(numel(neg), n_neg));
idx = [pos; neg];
idx = idx(randperm(numel(idx)));
end

%% Prediction summaries

function P = prediction_table(Te, risk, FOLD, target, detector)
vars = ["condition_key","geometry","geometry_family","case_id","case_family", ...
    "realism_level","field_regime","source_seed","noise_seed","f0","M","REQ_M", ...
    "map_iz","map_ix","x_center_m","z_center_m","true_SWS","SWS_pred", ...
    "q_pred","k_pred","predicted_patch_purity","p_mixed","p_strong_mixed", ...
    "local_amplitude","tracking_snr_db","acoustic_readout_amplitude_factor", ...
    "depth_m","patch_purity","distance_to_interface_m", ...
    "distance_to_interface_over_window_radius","purity_bin","distance_bin", ...
    "distance_over_window_bin","roi_region","sws_signed_error_pct", ...
    "sws_abs_error_pct","high_error10","high_error20","grad_q_pred","grad_sws_pred"];
vars = intersect(vars, string(Te.Properties.VariableNames), 'stable');
P = Te(:, cellstr(vars));
P.split_name = repmat(FOLD.split_name, height(P), 1);
P.fold_name = repmat(FOLD.fold_name, height(P), 1);
P.target_name = repmat(string(target), height(P), 1);
P.detector_name = repmat(string(detector), height(P), 1);
P.predicted_risk = risk(:);
P.predicted_reliability = 1 - risk(:);
P.nonestimable_risk_gt_0p2 = risk(:) >= 0.2;
P.nonestimable_risk_gt_0p5 = risk(:) >= 0.5;
P.nonestimable_risk_gt_0p8 = risk(:) >= 0.8;
end

function MET = detector_metrics(P, y, risk)
MET = table(P.split_name(1), P.fold_name(1), P.target_name(1), P.detector_name(1), ...
    height(P), 100*mean(y,'omitnan'), auc_binary(y, risk, "roc"), ...
    auc_binary(y, risk, "pr"), mean((risk - double(y)).^2, 'omitnan'), ...
    mean(P.sws_abs_error_pct,'omitnan'), 100*mean(P.high_error20,'omitnan'), ...
    'VariableNames', {'split_name','fold_name','target_name','detector_name', ...
    'N','positive_rate_pct','ROC_AUC','PR_AUC','Brier_score','MAPE_pct','high20_pct'});
end

function TH = threshold_summary(P, CFG)
parts = {};
for th = CFG.RiskThresholds
    reject = P.predicted_risk >= th;
    keep = ~reject;
    parts{end+1,1} = table(P.split_name(1), P.fold_name(1), P.target_name(1), ...
        P.detector_name(1), th, height(P), 100*mean(reject,'omitnan'), ...
        100*mean(keep,'omitnan'), mean(P.sws_abs_error_pct(keep),'omitnan'), ...
        mean(P.sws_abs_error_pct(reject),'omitnan'), ...
        mean(P.sws_signed_error_pct(keep),'omitnan'), ...
        mean(P.sws_signed_error_pct(reject),'omitnan'), ...
        100*mean(P.high_error20(keep),'omitnan'), ...
        100*mean(P.high_error20(reject),'omitnan'), ...
        100*recall_of_rejected(P.high_error20, reject), ...
        100*precision_of_rejected(P.high_error20, reject), ...
        'VariableNames', {'split_name','fold_name','target_name','detector_name', ...
        'risk_threshold','N','rejected_pct','coverage_pct','MAPE_estimable_pct', ...
        'MAPE_nonestimable_pct','bias_estimable_pct','bias_nonestimable_pct', ...
        'high20_estimable_pct','high20_nonestimable_pct','recall_high20_pct', ...
        'precision_high20_pct'}); %#ok<AGROW>
end
TH = vertcat(parts{:});
end

function COV = accuracy_coverage(P, CFG)
[~, ord] = sort(P.predicted_risk, 'ascend');
P = P(ord,:);
parts = {};
for cov = CFG.CoverageGrid
    n = max(1, round(cov * height(P)));
    X = P(1:n,:);
    parts{end+1,1} = table(P.split_name(1), P.fold_name(1), P.target_name(1), ...
        P.detector_name(1), cov, n, mean(X.sws_abs_error_pct,'omitnan'), ...
        mean(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error20,'omitnan'), ...
        mean(X.predicted_risk,'omitnan'), ...
        'VariableNames', {'split_name','fold_name','target_name','detector_name', ...
        'coverage_fraction','N_kept','MAPE_pct','bias_pct','high20_pct','mean_risk'});
end
COV = vertcat(parts{:});
end

function CAL = calibration_summary(P, CFG)
parts = {};
for i = 1:numel(CFG.CalibrationEdges)-1
    lo = CFG.CalibrationEdges(i); hi = CFG.CalibrationEdges(i+1);
    idx = P.predicted_risk >= lo & P.predicted_risk < hi;
    if i == numel(CFG.CalibrationEdges)-1
        idx = P.predicted_risk >= lo & P.predicted_risk <= hi;
    end
    parts{end+1,1} = table(P.split_name(1), P.fold_name(1), P.target_name(1), ...
        P.detector_name(1), lo, hi, sum(idx), mean(P.predicted_risk(idx),'omitnan'), ...
        100*mean(P.high_error20(idx),'omitnan'), mean(P.sws_abs_error_pct(idx),'omitnan'), ...
        'VariableNames', {'split_name','fold_name','target_name','detector_name', ...
        'risk_bin_low','risk_bin_high','N','mean_predicted_risk', ...
        'observed_high20_pct','MAPE_pct'}); %#ok<AGROW>
end
CAL = vertcat(parts{:});
end

function P = select_main_predictions(PRED, CFG)
m = PRED.target_name == "high20" & PRED.detector_name == "bagged_trees";
g = PRED.split_name == "grouped_condition_repeated";
if any(m & g)
    folds = unique(PRED.fold_name(m & g), 'stable');
    P = PRED(m & g & PRED.fold_name == folds(1), :);
else
    P = PRED(m, :);
end
if CFG.MaxPatchRowsToSave > 0 && height(P) > CFG.MaxPatchRowsToSave
    P = P(stratified_sample_idx((1:height(P))', P.high_error20, CFG.MaxPatchRowsToSave, 68123), :);
end
end

function P = cap_patch_table_for_saving(PRED, CFG)
P = PRED;
if CFG.MaxPatchRowsToSave > 0 && height(P) > CFG.MaxPatchRowsToSave
    P = P(stratified_sample_idx((1:height(P))', P.high_error20, CFG.MaxPatchRowsToSave, 68268), :);
end
end

function SUM = make_grouped_summaries(P, CFG)
SUM = struct();
SUM.geometry = summarize_group(P, ["geometry"]);
SUM.frequency = summarize_group(P, ["f0"]);
SUM.realism = summarize_group(P, ["realism_level"]);
SUM.field_regime = summarize_group(P, ["field_regime"]);
SUM.roi = summarize_group(P, ["roi_region"]);
SUM.purity = summarize_group(P, ["purity_bin"]);
SUM.distance = summarize_group(P, ["distance_over_window_bin"]);

% Add a standard threshold summary by region for the main detector.
SUM.roi = add_group_thresholds(P, SUM.roi, ["roi_region"], CFG);
SUM.geometry = add_group_thresholds(P, SUM.geometry, ["geometry"], CFG);
end

function S = summarize_group(P, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(P(:, group_vars));
rows = {};
for gi = 1:max(G)
    X = P(G == gi,:);
    rows{end+1,1} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), 100*mean(X.high_error20,'omitnan'), ...
        mean(X.predicted_risk,'omitnan'), 100*mean(X.predicted_risk >= 0.5,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_APE_pct','bias_pct', ...
        'high10_pct','high20_pct','mean_risk','rejected_at_0p5_pct'}); %#ok<AGROW>
end
S = [groups vertcat(rows{:})];
end

function S = add_group_thresholds(P, S, group_vars, CFG)
if isempty(S), return; end
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(P(:, group_vars));
for th = CFG.RiskThresholds
    mape_keep = nan(height(groups),1);
    high_keep = nan(height(groups),1);
    rej = nan(height(groups),1);
    for gi = 1:max(G)
        X = P(G==gi,:);
        keep = X.predicted_risk < th;
        mape_keep(gi) = mean(X.sws_abs_error_pct(keep),'omitnan');
        high_keep(gi) = 100*mean(X.high_error20(keep),'omitnan');
        rej(gi) = 100*mean(~keep,'omitnan');
    end
    suffix = strrep(sprintf('%.1f', th), '.', 'p');
    S.("MAPE_estimable_t" + suffix) = mape_keep;
    S.("high20_estimable_t" + suffix) = high_keep;
    S.("rejected_t" + suffix + "_pct") = rej;
end
end

%% Figures

function plot_auc_summary(METRICS, OUT)
X = METRICS(METRICS.target_name=="high20",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 14]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
boxchart(ax, categorical(pretty_detector(X.detector_name)), X.PR_AUC);
ylabel(ax,'PR AUC for high-error >20%'); title(ax,'Risk detector PR AUC','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
boxchart(ax, categorical(pretty_detector(X.detector_name)), X.ROC_AUC);
ylabel(ax,'ROC AUC for high-error >20%'); title(ax,'Risk detector ROC AUC','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_detector_auc_summary.png'));
end

function plot_calibration(CALIB, OUT)
X = CALIB(CALIB.target_name=="high20",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 16]);
tl = tiledlayout(fig,1,1,'Padding','compact');
ax = nexttile(tl); hold(ax,'on');
groups = unique(X.detector_name, 'stable');
for g = groups'
    Y = X(X.detector_name==g & X.N>0,:);
    scatter(ax, Y.mean_predicted_risk, Y.observed_high20_pct/100, 35, 'filled', ...
        'DisplayName', pretty_detector(g));
end
plot(ax,[0 1],[0 1],'k--','DisplayName','ideal');
xlabel(ax,'Predicted high-error probability');
ylabel(ax,'Observed high-error fraction');
title(ax,'Risk calibration for high-error >20%','FontWeight','normal');
legend(ax,'Location','eastoutside'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_risk_calibration.png'));
end

function plot_accuracy_coverage(COVERAGE, OUT)
X = COVERAGE(COVERAGE.target_name=="high20",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on');
for d = unique(X.detector_name,'stable')'
    Y = aggregate_xy(X(X.detector_name==d,:), "coverage_fraction", "MAPE_pct");
    plot(ax, 100*Y.x, Y.y, '-o', 'DisplayName', pretty_detector(d));
end
xlabel(ax,'Kept pixels by lowest risk (%)'); ylabel(ax,'MAPE in kept pixels (%)');
title(ax,'Accuracy-coverage curve','FontWeight','normal'); legend(ax,'Location','best'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for d = unique(X.detector_name,'stable')'
    Y = aggregate_xy(X(X.detector_name==d,:), "coverage_fraction", "high20_pct");
    plot(ax, 100*Y.x, Y.y, '-o', 'DisplayName', pretty_detector(d));
end
xlabel(ax,'Kept pixels by lowest risk (%)'); ylabel(ax,'High-error >20% in kept pixels (%)');
title(ax,'High-error rate versus coverage','FontWeight','normal'); legend(ax,'Location','best'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_accuracy_coverage.png'));
end

function plot_threshold_tradeoff(THRESH, OUT)
X = THRESH(THRESH.target_name=="high20",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 16]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
metrics = ["rejected_pct","recall_high20_pct","precision_high20_pct","MAPE_estimable_pct"];
labels = ["Rejected as non-estimable (%)","Recall of high-error patches (%)", ...
    "Precision among rejected patches (%)","MAPE in estimable patches (%)"];
for i = 1:numel(metrics)
    ax = nexttile(tl); hold(ax,'on');
    for d = unique(X.detector_name,'stable')'
        Y = aggregate_xy(X(X.detector_name==d,:), "risk_threshold", metrics(i));
        plot(ax, Y.x, Y.y, '-o', 'DisplayName', pretty_detector(d));
    end
    xlabel(ax,'Risk threshold'); ylabel(ax, labels(i));
    title(ax, labels(i), 'FontWeight','normal'); grid(ax,'on');
end
legend(nexttile(tl,1), 'Location','best');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_threshold_tradeoff.png'));
end

function plot_diagnostic_risk_vs_interface(P, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 16]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
boxchart(ax, categorical(string(P.distance_over_window_bin)), P.predicted_risk);
ylabel(ax,'Predicted risk'); xlabel(ax,'Distance / window radius bin');
title(ax,'Risk versus distance to interface (diagnostic only)','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
boxchart(ax, categorical(string(P.purity_bin)), P.predicted_risk);
ylabel(ax,'Predicted risk'); xlabel(ax,'True patch purity bin');
title(ax,'Risk versus true patch purity (diagnostic only)','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
scatter(ax, P.distance_to_interface_over_window_radius, P.sws_abs_error_pct, 4, P.predicted_risk, 'filled');
xlabel(ax,'Distance / window radius'); ylabel(ax,'Absolute SWS error (%)');
title(ax,'Error and risk near interfaces','FontWeight','normal'); cb=colorbar(ax); ylabel(cb,'Predicted risk');
grid(ax,'on'); ylim(ax,[0 prctile(P.sws_abs_error_pct,99)]);
ax = nexttile(tl);
scatter(ax, P.patch_purity, P.sws_abs_error_pct, 4, P.predicted_risk, 'filled');
xlabel(ax,'True patch purity'); ylabel(ax,'Absolute SWS error (%)');
title(ax,'Error and risk versus purity','FontWeight','normal'); cb=colorbar(ax); ylabel(cb,'Predicted risk');
grid(ax,'on'); ylim(ax,[0 prctile(P.sws_abs_error_pct,99)]);
export_fig(fig, fullfile(OUT.figure_dir, 'test68_risk_vs_distance_purity_diagnostic.png'));
end

function plot_risk_by_frequency_regime(Sfreq, Sregime, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 14]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
plot(ax, Sfreq.f0, Sfreq.MAPE_pct, '-o'); grid(ax,'on');
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)');
title(ax,'Error by frequency','FontWeight','normal');
ax = nexttile(tl);
bar(ax, categorical(pretty_regime(Sregime.field_regime)), Sregime.MAPE_pct);
xtickangle(ax,20); ylabel(ax,'MAPE (%)'); title(ax,'Error by field regime','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_risk_by_frequency_regime.png'));
end

function plot_region_heatmaps(Sroi, Sgeom, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
bar(ax, categorical(pretty_region(Sroi.roi_region)), Sroi.MAPE_pct);
xtickangle(ax,25); ylabel(ax,'MAPE (%)'); title(ax,'MAPE by diagnostic ROI','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
bar(ax, categorical(string(Sgeom.geometry)), Sgeom.MAPE_pct);
xtickangle(ax,25); ylabel(ax,'MAPE (%)'); title(ax,'MAPE by geometry','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test68_region_geometry_summary.png'));
end

function plot_representative_maps(P, OUT, CFG)
P = P(P.target_name=="high20" & P.detector_name=="bagged_trees", :);
if isempty(P), return; end
[G, groups] = findgroups(P(:, {'condition_key','geometry','f0','realism_level','field_regime'}));
case_mape = splitapply(@(x) mean(x,'omitnan'), P.sws_abs_error_pct, G);
case_risk = splitapply(@(x) mean(x,'omitnan'), P.predicted_risk, G);
CASE = [groups table(case_mape, case_risk, 'VariableNames', {'MAPE_pct','mean_risk'})];
CASE = sortrows(CASE, 'MAPE_pct', 'descend');
keep = [1:min(height(CASE), ceil(CFG.MaxMapConditions/2)), ...
    max(1,height(CASE)-floor(CFG.MaxMapConditions/2)+1):height(CASE)];
keep = unique(keep, 'stable');
CASE = CASE(keep,:);
for i = 1:height(CASE)
    key = CASE.condition_key(i);
    X = P(P.condition_key==key,:);
    out_dir = fullfile(OUT.map_dir, sanitize(CASE.geometry(i)), sanitize(CASE.realism_level(i)), sanitize(CASE.field_regime(i)));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    plot_one_map_condition(X, fullfile(out_dir, "test68_estimability__" + sanitize(key) + ".png"));
end
end

function plot_one_map_condition(X, path)
if isempty(X), return; end
panels = ["true_SWS","SWS_pred","sws_abs_error_pct","sws_signed_error_pct", ...
    "predicted_risk","predicted_reliability","nonestimable_risk_gt_0p5", ...
    "predicted_patch_purity","p_mixed","q_pred","grad_q_pred", ...
    "tracking_snr_db","local_amplitude","patch_purity","distance_to_interface_over_window_radius"];
titles = ["True SWS","Predicted SWS","Absolute error","Signed error", ...
    "Predicted risk","Reliability","Non-estimable mask", ...
    "Predicted purity","Predicted mixedness","Predicted q","|grad q|", ...
    "Tracking SNR","Amplitude","True patch purity","Distance / window radius"];
units = ["m/s","m/s","%","%","probability","probability","0/1", ...
    "probability","probability","q","q/pixel","dB","a.u.","fraction","ratio"];
fig = figure('Color','w','Units','centimeters','Position',[1 1 38 24]);
tl = tiledlayout(fig,3,5,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s | f=%g Hz | %s | %s', X.geometry(1), X.f0(1), ...
    pretty_realism(X.realism_level(1)), pretty_regime(X.field_regime(1))), ...
    'Interpreter','none','FontWeight','bold');
for i = 1:numel(panels)
    ax = nexttile(tl);
    A = map_from_table(X, panels(i));
    imagesc(ax, A); axis(ax,'image'); axis(ax,'off');
    title(ax, titles(i), 'FontWeight','normal');
    cb = colorbar(ax); ylabel(cb, units(i));
end
export_fig(fig, path);
end

function A = map_from_table(T, var)
[A, ~, ~] = vector_to_map(T.map_iz, T.map_ix, double(T.(var)));
end

%% README and console

function write_readme(P, METRICS, THRESH, SUM, OUT, CFG)
main_auc = METRICS(METRICS.target_name=="high20" & METRICS.detector_name=="bagged_trees",:);
main_th = THRESH(THRESH.target_name=="high20" & THRESH.detector_name=="bagged_trees" & THRESH.risk_threshold==0.5,:);
lines = strings(0,1);
lines(end+1) = "# Test 68 Results: Test66 Estimability Mask";
lines(end+1) = "";
lines(end+1) = sprintf("- Mode: `%s`", CFG.Mode);
lines(end+1) = sprintf("- Source Test66 run: `%s`", CFG.Test66RunDir);
lines(end+1) = "- No SWS correction and no q-model retraining were performed.";
lines(end+1) = "- Risk features are operational only; true SWS, q_oracle, true patch purity, distance to interface, and errors are excluded from predictors.";
lines(end+1) = "";
lines(end+1) = "## Main High-Error Detector";
lines(end+1) = sprintf("- Bagged-tree PR AUC high-error >20%%: %.3f ± %.3f", ...
    mean(main_auc.PR_AUC,'omitnan'), std(main_auc.PR_AUC,'omitnan'));
lines(end+1) = sprintf("- Bagged-tree ROC AUC high-error >20%%: %.3f ± %.3f", ...
    mean(main_auc.ROC_AUC,'omitnan'), std(main_auc.ROC_AUC,'omitnan'));
if ~isempty(main_th)
    lines(end+1) = sprintf("- At risk threshold 0.5: rejected %.1f%%, high-error recall %.1f%%, MAPE kept %.2f%%, MAPE rejected %.2f%%.", ...
        mean(main_th.rejected_pct,'omitnan'), mean(main_th.recall_high20_pct,'omitnan'), ...
        mean(main_th.MAPE_estimable_pct,'omitnan'), mean(main_th.MAPE_nonestimable_pct,'omitnan'));
end
lines(end+1) = "";
lines(end+1) = "## Diagnostic Interpretation";
lines(end+1) = "- A useful mask should assign high risk to interface/mixed patches and low risk to homogeneous/core patches.";
lines(end+1) = "- If MAPE in estimable pixels is much lower than in non-estimable pixels, this layer is useful even without correction.";
lines(end+1) = "- Test68 should be paired with Test67: Test67 showed that q_oracle removes most interface error, while Test68 tells us where the frozen q model should not be trusted.";
lines(end+1) = "";
lines(end+1) = "## Output Tables";
lines(end+1) = "- `test68_detector_metrics.csv`";
lines(end+1) = "- `test68_threshold_summary.csv`";
lines(end+1) = "- `test68_accuracy_coverage.csv`";
lines(end+1) = "- grouped summaries by geometry, frequency, realism, field regime, ROI, purity, and distance.";
lines(end+1) = "";
lines(end+1) = "## Suggested Next Step";
lines(end+1) = "- Use this mask first as a reporting/estimability layer. Only after choosing a conservative threshold should M-eff correction be gated by low risk/high estimability.";
lines(end+1) = "";
lines(end+1) = "## Compact Region Summaries";
lines(end+1) = "See CSV tables in `tables/`; the most useful first reads are `test68_summary_by_roi.csv` and `test68_summary_by_distance_over_window.csv`.";
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w');
fprintf(fid, '%s\n', lines);
fclose(fid);
end

function print_console_summary(METRICS, THRESH, OUT)
fprintf('\n================ Test 68 summary ================\n');
M = METRICS(METRICS.target_name=="high20",:);
disp(M(:, {'split_name','fold_name','detector_name','N','positive_rate_pct','PR_AUC','ROC_AUC','Brier_score'}));
T = THRESH(THRESH.target_name=="high20" & THRESH.detector_name=="bagged_trees" & THRESH.risk_threshold==0.5,:);
if ~isempty(T)
    fprintf('\nBagged trees at risk threshold 0.5:\n');
    disp(T(:, {'split_name','fold_name','rejected_pct','MAPE_estimable_pct','MAPE_nonestimable_pct','recall_high20_pct','precision_high20_pct'}));
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

%% Utility functions

function assert_no_forbidden_predictors(features)
low = lower(string(features));
bad_exact = ["patch_purity","true_sws","k_true","q_oracle","roi_region", ...
    "analysis_region","distance_to_interface_m","distance_to_interface_mm", ...
    "distance_to_interface_over_window_radius","sws_abs_error_pct", ...
    "sws_signed_error_pct","high_error10","high_error20"];
bad_patterns = ["oracle","true","error","distance_to_interface","roi","purity_bin", ...
    "distance_bin","analysis_region"];
hit = low(ismember(low, bad_exact));
assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
for p = bad_patterns
    hit = low(contains(low, p));
    assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
end
end

function y = first_existing_numeric(T, names)
for n = names
    if ismember(n, string(T.Properties.VariableNames))
        y = double(T.(n));
        return;
    end
end
error('Missing all candidate variables: %s', strjoin(names, ', '));
end

function y = ensure_numeric_var(T, name, default)
if ismember(name, string(T.Properties.VariableNames))
    y = double(T.(name));
else
    y = default;
end
end

function y = auc_binary(ytrue, score, mode)
ytrue = logical(ytrue(:));
score = double(score(:));
ok = isfinite(score);
ytrue = ytrue(ok); score = score(ok);
if numel(unique(ytrue)) < 2
    y = NaN;
    return;
end
[~, ord] = sort(score, 'descend');
ytrue = ytrue(ord);
tp = cumsum(ytrue);
fp = cumsum(~ytrue);
P = sum(ytrue); N = sum(~ytrue);
rec = tp / max(P, eps);
fpr = fp / max(N, eps);
prec = tp ./ max(tp + fp, eps);
switch string(mode)
    case "roc"
        y = trapz([0; fpr; 1], [0; rec; 1]);
    otherwise
        y = trapz([0; rec], [1; prec]);
end
end

function r = recall_of_rejected(high, rejected)
high = logical(high); rejected = logical(rejected);
r = sum(high & rejected) / max(sum(high), eps);
end

function p = precision_of_rejected(high, rejected)
high = logical(high); rejected = logical(rejected);
p = sum(high & rejected) / max(sum(rejected), eps);
end

function S = aggregate_xy(T, xvar, yvar)
x = unique(T.(xvar), 'stable');
y = nan(numel(x),1);
for i = 1:numel(x)
    y(i) = mean(T.(yvar)(T.(xvar)==x(i)), 'omitnan');
end
S = table(x, y);
end

function y = pretty_detector(x)
x = string(x);
y = x;
y(x=="bagged_trees") = "Bagged trees";
y(x=="logistic_regression") = "Logistic regression";
end

function y = pretty_regime(x)
x = string(x);
y = replace(x, "_", " ");
y = replace(y, "single source lateral", "single source");
y = replace(y, "diffuse like 8src", "diffuse-like 8 sources");
end

function y = pretty_realism(x)
x = string(x);
y = replace(x, "_", " ");
end

function y = pretty_region(x)
x = string(x);
y = replace(x, "_", " ");
y(y=="background far") = "background/core";
y(y=="hard inclusion core") = "hard inclusion core";
y(y=="interface 0 1mm") = "interface 0-1 mm";
end

function s = sanitize(x)
s = regexprep(char(string(x)), '[^A-Za-z0-9_=-]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'unnamed'; end
end

function export_fig(fig, path)
folder = fileparts(path);
if exist(folder,'dir') ~= 7, mkdir(folder); end
exportgraphics(fig, path, 'Resolution', 220);
close(fig);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Test68:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
end
end

function write_config_json(CFG, path)
try
    txt = jsonencode(CFG, 'PrettyPrint', true);
catch
    txt = jsonencode(CFG);
end
fid = fopen(path, 'w');
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function y = ternary(tf, a, b)
if tf
    y = a;
else
    y = b;
end
end

function v = env_string(name, default)
raw = string(getenv(name));
if strlength(raw) == 0
    v = string(default);
else
    v = raw;
end
end

function v = env_number(name, default)
raw = string(getenv(name));
if strlength(raw) == 0
    v = default;
else
    v = str2double(raw);
    if ~isfinite(v), v = default; end
end
end

function tf = env_true(name, default)
raw = lower(strtrim(string(getenv(name))));
if strlength(raw) == 0
    tf = default;
else
    tf = any(raw == ["1","true","yes","on"]);
end
end
