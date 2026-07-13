%% analyze_test_45_candidate_aware_no_theory_correction.m
% Test 45: candidate-aware no-theory correction for
% q_spectrum_plus_composition.
%
% This analysis keeps the base Test38 q model frozen. It trains only a
% conservative operational correction layer that tries to decide:
%   1) whether to keep the base estimate;
%   2) whether SWS should increase or decrease;
%   3) how large the correction should be;
%   4) what final high-error risk remains.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST45_MODE          = quick | medium | full
%   ADAPTIVE_REQ_TEST45_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST45_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST45_USE_PARALLEL  = true | false
%   ADAPTIVE_REQ_TEST45_SOURCE        = test39 | kwave | both

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.08);

CFG = default_config();
CFG = setup_parallel_if_requested(CFG);
OUT = make_output_dirs(root_dir, CFG);
SRC = locate_sources(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test45_configuration.json'));

fprintf('\nTest 45: candidate-aware no-theory correction\n');
fprintf('Mode: %s | source: %s | validate only: %d\n', ...
    CFG.Mode, CFG.Source, CFG.ValidateOnly);
fprintf('Frozen base model: %s\n', CFG.BaseModel);
fprintf('Training correction on Test38 rows only. Base q model is not retrained.\n');

assert(exist(SRC.test38_csv,'file') == 2, 'Missing Test38 table: %s', SRC.test38_csv);
if CFG.UseTest39
    assert(exist(SRC.test39_csv,'file') == 2, 'Missing Test39 table: %s', SRC.test39_csv);
end
if CFG.UseKwave
    assert(exist(SRC.test40_csv,'file') == 2, 'Missing Test40/k-Wave table: %s', SRC.test40_csv);
end

T38_all = standardize_synthetic_table(readtable(SRC.test38_csv), "synthetic_train_source");
T38_all = add_variant_reference_features(T38_all);
T38_base = T38_all(T38_all.model_name == CFG.BaseModel,:);
assert(any(T38_base.is_train_row), 'No Test38 training rows for %s.', CFG.BaseModel);
assert(any(T38_base.is_heldout_row), 'No Test38 held-out rows for %s.', CFG.BaseModel);

eval_parts = {};
eval_parts{end+1,1} = T38_base(T38_base.is_heldout_row,:); %#ok<SAGROW>
eval_parts{end}.source_domain = repmat("synthetic_heldout", height(eval_parts{end}), 1);

if CFG.UseTest39
    T39 = standardize_synthetic_table(readtable(SRC.test39_csv), "test39_ood");
    T39 = add_variant_reference_features(T39);
    T39 = T39(T39.model_name == CFG.BaseModel,:);
    T39.source_domain = repmat("test39_ood", height(T39), 1);
    eval_parts{end+1,1} = T39; %#ok<SAGROW>
end
if CFG.UseKwave
    T40 = standardize_kwave_table(readtable(SRC.test40_csv));
    T40 = add_variant_reference_features(T40);
    T40 = T40(T40.model_name == CFG.BaseModel,:);
    if CFG.QuickMode
        T40 = quick_subset_by_condition(T40, CFG);
    end
    eval_parts{end+1,1} = T40; %#ok<SAGROW>
end
T_eval = vertcat_compatible(eval_parts{:});

T_train = T38_base(T38_base.is_train_row,:);
if CFG.QuickMode
    T_train = stratified_row_sample(T_train, CFG.QuickTrainRows, CFG.RandomSeed);
end
fprintf('Train rows: %d | evaluation rows: %d.\n', height(T_train), height(T_eval));

if CFG.ValidateOnly
    T_train = stratified_row_sample(T_train, min(2500,height(T_train)), CFG.RandomSeed);
    T_eval = stratified_row_sample(T_eval, min(2500,height(T_eval)), CFG.RandomSeed+1);
    fprintf('Validation-only subset: train=%d | eval=%d.\n', height(T_train), height(T_eval));
end

%% Train correction stack

rng(CFG.RandomSeed, 'twister');
[MODEL, TRAIN_INFO] = train_direction_stack(T_train, CFG);

if CFG.ValidateOnly
    validate_predictor_policy(MODEL.feature_names);
end

%% Apply strategies

T_results = apply_all_strategies(T_eval, MODEL, CFG);

if CFG.ValidateOnly
    assert(all(isfinite(T_results.sws_pred_final)), ...
        'Non-finite corrected SWS values found.');
end

%% Summaries and outputs

T_overall = summarize_predictions(T_results, ["strategy_name"]);
T_by_source = summarize_predictions(T_results, ["source_domain","strategy_name"]);
T_by_M = summarize_predictions(T_results, ["source_domain","M","strategy_name"]);
T_by_frequency = summarize_predictions(T_results, ["source_domain","f0","strategy_name"]);
T_by_frequency_M = summarize_predictions(T_results, ["source_domain","M","f0","strategy_name"]);
T_by_regime = summarize_predictions(T_results, ["source_domain","field_regime_ood","strategy_name"]);
T_by_region = summarize_predictions(T_results, ["source_domain","evaluation_region","strategy_name"]);
T_by_predmix = summarize_predictions(T_results, ["source_domain","predicted_mixedness_bin","strategy_name"]);
T_by_truepurity = summarize_predictions(T_results, ["source_domain","true_purity_bin","strategy_name"]);
T_by_distance = summarize_predictions(T_results, ["source_domain","distance_bin_pretty","strategy_name"]);
T_selector = candidate_selector_metrics(T_results);
T_gate = gate_metrics(T_results);
T_cal = reliability_calibration(T_results);
T_best = best_strategy_candidates(T_by_source, T_by_region, CFG);
T_reference = reference_comparison_table(root_dir, T_by_source);

writetable(T_results, fullfile(OUT.table_dir, 'test45_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test45_strategy_summary_overall.csv'));
writetable(T_by_source, fullfile(OUT.table_dir, 'test45_strategy_summary_by_source.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test45_strategy_summary_by_M.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'test45_strategy_summary_by_frequency.csv'));
writetable(T_by_frequency_M, fullfile(OUT.table_dir, 'test45_strategy_summary_by_frequency_M.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test45_strategy_summary_by_regime.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test45_strategy_summary_by_region.csv'));
writetable(T_by_predmix, fullfile(OUT.table_dir, 'test45_strategy_summary_by_predicted_mixedness.csv'));
writetable(T_by_truepurity, fullfile(OUT.table_dir, 'test45_strategy_summary_by_true_purity.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test45_strategy_summary_by_distance.csv'));
writetable(T_selector, fullfile(OUT.table_dir, 'test45_candidate_selector_metrics.csv'));
writetable(T_gate, fullfile(OUT.table_dir, 'test45_gate_metrics.csv'));
writetable(T_cal, fullfile(OUT.table_dir, 'test45_reliability_calibration.csv'));
writetable(T_best, fullfile(OUT.table_dir, 'test45_best_strategy_candidates.csv'));
writetable(T_reference, fullfile(OUT.table_dir, 'test45_reference_comparison_test33_test34.csv'));

save(fullfile(OUT.data_dir, 'test45_compact_results.mat'), ...
    'CFG','SRC','MODEL','TRAIN_INFO','T_results','T_overall','T_by_source', ...
    'T_by_M','T_by_frequency','T_by_frequency_M','T_by_regime','T_by_region','T_by_predmix', ...
    'T_by_truepurity','T_by_distance','T_selector','T_gate','T_cal','T_best','T_reference','-v7.3');

%% Figures

safe_plot(@() plot_strategy_ranking(T_by_source, OUT, "MAPE_pct", ...
    'test45_strategy_ranking_mape.png'), 'strategy ranking MAPE');
safe_plot(@() plot_strategy_ranking(T_by_source, OUT, "high_error20_pct", ...
    'test45_strategy_ranking_high20.png'), 'strategy ranking high20');
safe_plot(@() plot_correction_benefit_harm(T_gate, OUT), 'benefit/harm');
safe_plot(@() plot_direction_confusion(T_results, OUT), 'direction confusion');
safe_plot(@() plot_error_before_after(T_results, OUT, true), 'signed error before/after');
safe_plot(@() plot_error_before_after(T_results, OUT, false), 'absolute error before/after');
safe_plot(@() plot_error_by_region(T_by_region, OUT), 'error by region');
safe_plot(@() plot_error_vs_distance(T_by_distance, OUT), 'error vs distance');
safe_plot(@() plot_error_vs_predicted_mixedness(T_by_predmix, OUT), 'error vs predicted mixedness');
safe_plot(@() plot_correction_magnitude_vs_mixedness(T_results, OUT), 'correction magnitude vs mixedness');
safe_plot(@() plot_reliability_calibration(T_cal, OUT), 'reliability calibration');
safe_plot(@() plot_frequency_dependence(T_by_frequency_M, OUT, 2), 'frequency dependence M2');
safe_plot(@() plot_frequency_dependence(T_by_frequency_M, OUT, 3), 'frequency dependence M3');
safe_plot(@() plot_net_gain_heatmap(T_by_region, OUT, CFG), 'net gain heatmap');
safe_plot(@() plot_reference_comparison(T_reference, OUT), 'Test33/Test34 comparison');
if CFG.SaveAllMaps || CFG.ValidateOnly
    safe_plot(@() plot_condition_maps(T_results, CFG, OUT), 'condition maps');
end

print_console_summary(T_overall, T_by_source, T_by_region, T_gate, T_best, CFG);
fprintf('\nTables: %s\nFigures: %s\nData: %s\nTest 45 complete.\n', ...
    OUT.table_dir, OUT.figure_dir, OUT.data_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST45_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","medium","full"]), ...
    'ADAPTIVE_REQ_TEST45_MODE must be quick, medium, or full.');
source = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST45_SOURCE'))));
if source == "", source = "both"; end
assert(ismember(source, ["test39","kwave","both"]), ...
    'ADAPTIVE_REQ_TEST45_SOURCE must be test39, kwave, or both.');

CFG = struct();
CFG.Mode = mode;
CFG.Source = source;
CFG.QuickMode = mode == "quick";
CFG.MediumMode = mode == "medium";
CFG.FullMode = mode == "full";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST45_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST45_SAVE_ALL_MAPS', true);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST45_USE_PARALLEL', false);
CFG.UseTest39 = ismember(source, ["test39","both"]);
CFG.UseKwave = ismember(source, ["kwave","both"]);
CFG.BaseModel = "q_spectrum_plus_composition";
CFG.RandomSeed = 44001;
CFG.DirectionThresholdPct = env_number('ADAPTIVE_REQ_TEST45_DIRECTION_THRESHOLD_PCT', 5.0);
CFG.GainMarginPct = env_number('ADAPTIVE_REQ_TEST45_GAIN_MARGIN_PCT', 1.0);
CFG.HarmMarginPct = env_number('ADAPTIVE_REQ_TEST45_HARM_MARGIN_PCT', 1.0);
CFG.NumTrees = env_number('ADAPTIVE_REQ_TEST45_NUM_TREES', ternary(mode=="quick", 80, 160));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST45_MIN_LEAF_SIZE', 70);
CFG.PhysicalRangeSWS = [0.5 10.0];
CFG.QuickTrainRows = env_number('ADAPTIVE_REQ_TEST45_QUICK_TRAIN_ROWS', 25000);
CFG.QuickKwaveConditions = env_number('ADAPTIVE_REQ_TEST45_QUICK_KWAVE_CONDITIONS', 8);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST45_MAX_MAP_CONDITIONS', ternary(mode=="quick", 8, 32));
CFG.MapStyle = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST45_MAP_STYLE'))));
if CFG.MapStyle == "", CFG.MapStyle = "interp"; end
assert(ismember(CFG.MapStyle, ["patch","interp"]), ...
    'ADAPTIVE_REQ_TEST45_MAP_STYLE must be patch or interp.');
CFG.MapInterpScale = env_number('ADAPTIVE_REQ_TEST45_MAP_INTERP_SCALE', 4);
CFG.MainStrategies = ["base_q_spectrum_plus_composition", ...
    "candidate_q_spectrum_only", "candidate_region_median", ...
    "candidate_neighbor_median", "candidate_blend_region", ...
    "candidate_selector", "conservative_candidate_selector", ...
    "mixedness_candidate_selector", "oracle_best_candidate", ...
    "posterior_reliability_only"];
CFG.PlotStrategies = ["base_q_spectrum_plus_composition", ...
    "candidate_q_spectrum_only", "candidate_region_median", ...
    "candidate_neighbor_median", "candidate_selector", ...
    "conservative_candidate_selector", "oracle_best_candidate"];
end

function CFG = setup_parallel_if_requested(CFG)
if ~CFG.UseParallel
    CFG.BagOptions = statset('UseParallel', false);
    CFG.BoostOptions = statset('UseParallel', false);
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool), parpool('threads'); end
    CFG.BagOptions = statset('UseParallel', true);
    % MATLAB parallel ensemble fitting is reliable for Bag ensembles here.
    % Keep LSBoost serial to avoid the "Method must be Bag" restriction.
    CFG.BoostOptions = statset('UseParallel', false);
catch ME
    warning('Could not start/use parallel pool (%s). Continuing serial.', ME.message);
    CFG.UseParallel = false;
    CFG.BagOptions = statset('UseParallel', false);
    CFG.BoostOptions = statset('UseParallel', false);
end
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_45_candidate_aware_no_theory_correction');
if CFG.QuickMode
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
elseif CFG.MediumMode
    OUT.root_dir = fullfile(OUT.root_dir, 'medium');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function SRC = locate_sources(root_dir, CFG)
SRC = struct();
SRC.test38_csv = table_path(root_dir, 'test_38_velocity_field_diverse_q_training', ...
    CFG.Mode, 'test38_patch_level_predictions.csv');
SRC.test39_csv = table_path(root_dir, 'test_39_frozen_test38_external_validation', ...
    CFG.Mode, 'test39_external_patch_level_predictions.csv');
SRC.test40_csv = fullfile(root_dir, 'outputs', 'test_40_kwave_latest_model_comparison', ...
    'tables', 'test40_kwave_patch_level_model_comparison.csv');
end

function p = table_path(root_dir, output_name, source, file)
base = fullfile(root_dir, 'outputs', output_name);
if source == "quick"
    p = fullfile(base, 'quick', 'tables', file);
elseif source == "medium"
    p = fullfile(base, 'medium', 'tables', file);
else
    p = fullfile(base, 'tables', file);
end
end

%% Data standardization

function T = standardize_synthetic_table(T, source_domain)
T = convert_strings(T);
T = normalize_logical_columns(T);
T.source_domain = repmat(source_domain, height(T), 1);
T.model_name = string(T.model_name);
T.strategy_name_original = T.model_name;
T.sws_base = T.sws_pred;
T.q_base = T.q_pred;
T.case_name = T.case_id;
T.geometry_type = T.case_family;
if ~ismember("material_region", string(T.Properties.VariableNames))
    T.material_region = repmat("homogeneous", height(T), 1);
    T.material_region(T.case_family ~= "homogeneous" & T.true_SWS <= median(T.true_SWS,'omitnan')) = "soft";
    T.material_region(T.case_family ~= "homogeneous" & T.true_SWS > median(T.true_SWS,'omitnan')) = "hard";
end
if ~ismember("region_zone", string(T.Properties.VariableNames))
    T.region_zone = pretty_distance_zone(T.distance_to_boundary_mm);
end
if ~ismember("material_zone", string(T.Properties.VariableNames))
    T.material_zone = T.material_region + "_" + T.region_zone;
end
if ~ismember("p_strong_mixed", string(T.Properties.VariableNames))
    T.p_strong_mixed = zeros(height(T),1);
end
T.distance_to_interface_mm = T.distance_to_boundary_mm;
T.patch_key = T.condition_key + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
T.evaluation_region = evaluation_region(T);
T.true_purity_bin = true_purity_bin(T.patch_purity, T.case_family);
T.predicted_mixedness_bin = probability_bin(T.p_mixed);
T.distance_bin_pretty = distance_bin_pretty(T.distance_to_interface_mm, T.case_family);
end

function T = standardize_kwave_table(T)
T = convert_strings(T);
T = normalize_logical_columns(T);
T.source_domain = repmat("kwave_transfer", height(T), 1);
T.case_id = "kwave_" + sanitize(T.case_name);
T.case_family = repmat("kwave_inclusion", height(T), 1);
T.field_regime_ood = T.field_regime;
T.model_name = erase(T.strategy_name, "T38_");
T.model_name(startsWith(T.strategy_name, "T34_")) = T.strategy_name(startsWith(T.strategy_name, "T34_"));
T.strategy_name_original = T.strategy_name;
T.sws_base = T.sws_pred;
T.q_base = T.q_pred;
T.patch_key = T.condition_key + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
if ismember("predicted_mixed_probability", string(T.Properties.VariableNames))
    T.p_mixed = T.predicted_mixed_probability;
end
if ~ismember("p_strong_mixed", string(T.Properties.VariableNames))
    T.p_strong_mixed = max(0, min(1, 2*T.p_mixed - 1));
end
T.material_region = T.material_side;
T.region_zone = T.analysis_region;
T.material_zone = T.material_region + "_" + T.region_zone;
T.distance_to_boundary_mm = T.distance_to_interface_mm;
T.evaluation_region = evaluation_region(T);
T.true_purity_bin = true_purity_bin(T.patch_purity, T.case_family);
T.predicted_mixedness_bin = probability_bin(T.p_mixed);
T.distance_bin_pretty = distance_bin_pretty(T.distance_to_interface_mm, T.case_family);
end

function T = add_variant_reference_features(T)
T.q_spectrum_only_sws = lookup_model_numeric(T, "q_spectrum_only", "sws_pred");
T.q_spectrum_only_q = lookup_model_numeric(T, "q_spectrum_only", "q_pred");
T.variant_sws_disagreement = abs(T.sws_base - T.q_spectrum_only_sws) ./ max(T.sws_base, 0.25);
T.variant_q_disagreement = T.q_base - T.q_spectrum_only_q;
end

function y = lookup_model_numeric(T, model_name, value_name)
key = T.patch_key;
if ismember("model_name", string(T.Properties.VariableNames))
    R = T(T.model_name == model_name, cellstr(["patch_key", string(value_name)]));
else
    R = table();
end
if isempty(R)
    y = nan(height(T),1);
else
    [ok,loc] = ismember(key, R.patch_key);
    y = nan(height(T),1);
    y(ok) = R.(value_name)(loc(ok));
end
if all(~isfinite(y))
    if value_name == "sws_pred", y = T.sws_base; else, y = T.q_base; end
else
    fallback = ternary(value_name=="sws_pred", T.sws_base, T.q_base);
    y(~isfinite(y)) = fallback(~isfinite(y));
end
end

function y = lookup_strategy_numeric(T, strategy_name, value_name)
R = T(T.strategy_name == strategy_name, cellstr(["patch_key", string(value_name)]));
if isempty(R)
    y = nan(height(T),1);
else
    [ok,loc] = ismember(T.patch_key, R.patch_key);
    y = nan(height(T),1);
    y(ok) = R.(value_name)(loc(ok));
end
fallback = T.sws_pred;
y(~isfinite(y)) = fallback(~isfinite(y));
end

function T = convert_strings(T)
for v = ["dataset","source_domain","condition_key","case_id","case_family", ...
        "case_name","field_regime","field_regime_ood","purity_bin","distance_bin", ...
        "material_region","region_zone","material_zone","model_name","model_family", ...
        "strategy_name","strategy_name_original","geometry","geometry_type", ...
        "material_side","analysis_region","roi_name","region_label","bundle_id"]
    if ismember(v, string(T.Properties.VariableNames))
        T.(v) = string(T.(v));
    end
end
end

function T = normalize_logical_columns(T)
for v = ["is_train_row","is_heldout_row","is_mixed","is_strong_mixed", ...
        "high_error10","high_error20"]
    if ismember(v, string(T.Properties.VariableNames))
        T.(v) = logical(T.(v));
    end
end
end

function T = quick_subset_by_condition(T, CFG)
C = unique(T(:, {'condition_key','M'}), 'rows', 'stable');
C = C(1:min(height(C), CFG.QuickKwaveConditions),:);
[ok,~] = ismember(T(:, {'condition_key','M'}), C);
T = T(ok,:);
end

function T = stratified_row_sample(T, n, seed)
if height(T) <= n, return; end
rng(seed, 'twister');
if ismember("condition_key", string(T.Properties.VariableNames))
    [G,~] = findgroups(T.condition_key);
    keep = false(height(T),1);
    per = max(20, ceil(n / max(G)));
    for g = unique(G)'
        idx = find(G == g);
        idx = idx(randperm(numel(idx), min(per,numel(idx))));
        keep(idx) = true;
    end
    idx = find(keep);
    if numel(idx) > n, idx = idx(randperm(numel(idx), n)); end
else
    idx = randperm(height(T), n);
end
T = T(idx,:);
end

%% Training and application

function [MODEL, info] = train_direction_stack(T, CFG)
T = add_no_theory_candidates(T, CFG);
[X, feature_names] = build_operational_features(T, CFG);
validate_predictor_policy(feature_names);
[C, candidate_names] = candidate_matrix(T, CFG);
valid = all(isfinite(X),2) & all(isfinite(C),2) & isfinite(T.true_SWS) & T.true_SWS > 0;
X = X(valid,:); C = C(valid,:); Tv = T(valid,:);

err = abs(100*(C - Tv.true_SWS)./Tv.true_SWS);
[best_err, best_idx] = min(err, [], 2);
base_err = err(:,1);
apply_label = best_idx ~= 1 & best_err + CFG.GainMarginPct < base_err;
candidate_delta = C(sub2ind(size(C),(1:size(C,1))',best_idx)) - C(:,1);
harm_label = best_idx ~= 1 & best_err > base_err + CFG.HarmMarginPct;
high20_base = base_err > 20;

tree = templateTree('MinLeafSize', CFG.MinLeafSize, 'MaxNumSplits', 384);
selector_model = fitcensemble(X, categorical(candidate_names(best_idx)), ...
    'Method','Bag','Learners',tree,'NumLearningCycles',max(60,round(CFG.NumTrees/2)), ...
    'ClassNames',categorical(candidate_names),'Options',CFG.BagOptions);
[pred_label, pred_score] = predict(selector_model, X);
[pred_idx, pred_conf] = candidate_label_to_index(string(pred_label), pred_score, selector_model, candidate_names);
pred_sws = C(sub2ind(size(C),(1:size(C,1))',pred_idx));
pred_err = abs(100*(pred_sws - Tv.true_SWS)./Tv.true_SWS);
pred_apply_label = pred_idx ~= 1 & pred_err + CFG.GainMarginPct < base_err;
pred_harm_label = pred_idx ~= 1 & pred_err > base_err + CFG.HarmMarginPct;

gateX = gate_feature_matrix(X, pred_conf, pred_sws, Tv);
gate_names = [feature_names "candidate_confidence" "candidate_delta_sws" "candidate_delta_pct" "heuristic_candidate_need"];
apply_model = fitcensemble(gateX, pred_apply_label, ...
    'Method','Bag','Learners',tree,'NumLearningCycles',max(60,round(CFG.NumTrees/2)), ...
    'ClassNames',[false true],'Options',CFG.BagOptions);
harm_model = fitcensemble(gateX, pred_harm_label, ...
    'Method','Bag','Learners',tree,'NumLearningCycles',max(60,round(CFG.NumTrees/2)), ...
    'ClassNames',[false true],'Options',CFG.BagOptions);
reliability_model = fitcensemble(gateX, high20_base, ...
    'Method','Bag','Learners',tree,'NumLearningCycles',max(60,round(CFG.NumTrees/2)), ...
    'ClassNames',[false true],'Options',CFG.BagOptions);

MODEL = struct();
MODEL.feature_names = feature_names;
MODEL.gate_feature_names = gate_names;
MODEL.candidate_names = candidate_names;
MODEL.selector_model = selector_model;
MODEL.apply_model = apply_model;
MODEL.harm_model = harm_model;
MODEL.reliability_model = reliability_model;

info = table(height(Tv), mean(base_err,'omitnan'), mean(best_err,'omitnan'), ...
    100*mean(best_idx==1,'omitnan'), 100*mean(apply_label,'omitnan'), ...
    100*mean(harm_label,'omitnan'), mean(abs(candidate_delta),'omitnan'), ...
    'VariableNames', {'N_train','train_base_MAPE_pct','train_oracle_candidate_MAPE_pct', ...
    'train_best_is_base_pct','train_apply_label_pct','train_harm_label_pct', ...
    'train_mean_oracle_candidate_delta_mps'});
disp(info);
end

function T_all = apply_all_strategies(T, MODEL, CFG)
T = add_no_theory_candidates(T, CFG);
[X, feature_names] = build_operational_features(T, CFG);
assert(isequal(string(feature_names(:)), string(MODEL.feature_names(:))), 'Feature mismatch.');
[C, candidate_names] = candidate_matrix(T, CFG);
assert(isequal(candidate_names(:), MODEL.candidate_names(:)), 'Candidate mismatch.');

[pred_label, pred_score] = predict(MODEL.selector_model, X);
[pred_idx, pred_conf] = candidate_label_to_index(string(pred_label), pred_score, MODEL.selector_model, candidate_names);
pred_sws = C(sub2ind(size(C),(1:size(C,1))',pred_idx));
gateX = gate_feature_matrix(X, pred_conf, pred_sws, T);
[apply_cls, apply_score] = predict(MODEL.apply_model, gateX);
[harm_cls, harm_score] = predict(MODEL.harm_model, gateX);
[~, rel_score] = predict(MODEL.reliability_model, gateX);
p_apply = positive_score(MODEL.apply_model, apply_score);
p_harm = positive_score(MODEL.harm_model, harm_score);
p_high_error = positive_score(MODEL.reliability_model, rel_score);
apply_pred = class_to_logical(apply_cls);
harm_pred = class_to_logical(harm_cls);
need = heuristic_mixedness_need(T);
base_sws = T.sws_base;

candidate_gate = apply_pred & pred_idx ~= 1 & p_apply >= 0.55 & pred_conf >= 0.35 & ...
    p_harm <= 0.45 & ~harm_pred;
conservative_gate = apply_pred & pred_idx ~= 1 & p_apply >= 0.70 & pred_conf >= 0.50 & ...
    p_harm <= 0.25 & need >= 0.10 & ~harm_pred;
mixed_gate = pred_idx ~= 1 & need >= 0.35 & p_harm <= 0.50;

err = abs(100*(C - T.true_SWS)./T.true_SWS);
[best_err, best_idx] = min(err, [], 2);
base_err = err(:,1);
oracle_gate = best_idx ~= 1 & best_err + CFG.GainMarginPct < base_err;
oracle_sws = C(sub2ind(size(C),(1:size(C,1))',best_idx));

parts = {};
parts{end+1,1} = make_strategy_rows(T, "base_q_spectrum_plus_composition", ...
    base_sws, zeros(height(T),1), false(height(T),1), "base", 1, ...
    p_apply, p_harm, p_high_error, "operational"); %#ok<AGROW>
for ci = 2:numel(candidate_names)
    cname = candidate_names(ci);
    parts{end+1,1} = make_strategy_rows(T, "candidate_" + cname, ...
        C(:,ci), C(:,ci)-base_sws, true(height(T),1), cname, ones(height(T),1), ...
        p_apply, p_harm, p_high_error, "candidate_reference"); %#ok<AGROW>
end
parts{end+1,1} = make_strategy_rows(T, "candidate_selector", ...
    choose_sws(base_sws, pred_sws, candidate_gate), choose_delta(pred_sws-base_sws, candidate_gate), ...
    candidate_gate, string(pred_label), pred_conf, p_apply, p_harm, p_high_error, "operational"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "conservative_candidate_selector", ...
    choose_sws(base_sws, pred_sws, conservative_gate), choose_delta(pred_sws-base_sws, conservative_gate), ...
    conservative_gate, string(pred_label), pred_conf, p_apply, p_harm, p_high_error, "operational_conservative"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "mixedness_candidate_selector", ...
    choose_sws(base_sws, pred_sws, mixed_gate), choose_delta(pred_sws-base_sws, mixed_gate), ...
    mixed_gate, string(pred_label), pred_conf, p_apply, p_harm, p_high_error, "operational_ablation"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "oracle_best_candidate", ...
    choose_sws(base_sws, oracle_sws, oracle_gate), choose_delta(oracle_sws-base_sws, oracle_gate), ...
    oracle_gate, candidate_names(best_idx), ones(height(T),1), p_apply, p_harm, p_high_error, "diagnostic_oracle"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "posterior_reliability_only", ...
    base_sws, zeros(height(T),1), false(height(T),1), "base", 1, ...
    p_apply, p_harm, p_high_error, "reliability_only"); %#ok<AGROW>
T_all = vertcat_compatible(parts{:});
tmp_true_candidate = candidate_names(best_idx_for_rows(T_all, T, candidate_names, CFG));
T_all.true_candidate_label = tmp_true_candidate(:);
T_all.is_harmed = T_all.was_corrected & T_all.sws_abs_error_pct > T_all.base_abs_error_pct + CFG.HarmMarginPct;
T_all.is_improved = T_all.was_corrected & T_all.sws_abs_error_pct + CFG.GainMarginPct < T_all.base_abs_error_pct;
T_all.is_unchanged_like = T_all.was_corrected & ~T_all.is_harmed & ~T_all.is_improved;
T_all.net_abs_error_gain_pct_points = T_all.base_abs_error_pct - T_all.sws_abs_error_pct;
end

function R = make_strategy_rows(T, strategy_name, sws, correction_mps, was_corrected, ...
    candidate_label, candidate_conf, p_apply, p_harm, p_high_error, strategy_type)
R = T;
R.strategy_name = repmat(strategy_name, height(T), 1);
R.strategy_type = repmat(strategy_type, height(T), 1);
R.sws_pred_final = sws;
R.sws_pred = sws;
R.correction_mps = correction_mps;
R.correction_pct = 100*correction_mps ./ max(T.sws_base, 0.25);
R.was_corrected = was_corrected;
R.was_kept = ~was_corrected;
candidate_label = string(candidate_label);
if isscalar(candidate_label)
    candidate_label = repmat(candidate_label, height(T), 1);
end
candidate_label = candidate_label(:);
if isscalar(candidate_conf)
    candidate_conf = repmat(candidate_conf, height(T), 1);
end
candidate_conf = candidate_conf(:);
R.candidate_confidence = candidate_conf;
R.predicted_candidate = candidate_label;
R.predicted_direction = candidate_direction_label(sws, T.sws_base, 0.01);
R.direction_confidence = candidate_conf;
R.p_apply_correction = p_apply;
R.p_harm_correction = p_harm;
R.final_high_error_probability = p_high_error;
R.final_reliability = 1 - p_high_error;
R.base_signed_error_pct = 100*(T.sws_base - T.true_SWS)./T.true_SWS;
R.base_abs_error_pct = abs(R.base_signed_error_pct);
R.sws_signed_error_pct = 100*(sws - T.true_SWS)./T.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
end

function [X, names] = build_operational_features(T, CFG)
names = ["M_norm","f0_norm","regime_code","geometry_code", ...
    "predicted_patch_purity","predicted_impurity","p_mixed","p_strong_mixed", ...
    "log_sws_base","q_base","q_spectrum_only_log_sws", ...
    "q_spectrum_only_q","variant_sws_disagreement","variant_q_disagreement", ...
    "region_median_log_sws","neighbor_median_log_sws","blend_region_log_sws", ...
    "base_minus_region_pct","base_minus_neighbor_pct","base_sws_gradient_proxy", ...
    "mixedness_need_proxy"];
base = max(T.sws_base, 0.1);
qonly_sws = fill_default(T.q_spectrum_only_sws, base);
qonly_q = fill_default(T.q_spectrum_only_q, T.q_base);
region_sws = fill_default(T.candidate_region_median_sws, base);
neighbor_sws = fill_default(T.candidate_neighbor_median_sws, base);
blend_sws = fill_default(T.candidate_blend_region_sws, base);
purity = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-purity);
ps = fill_default(T.p_strong_mixed, 0);
grad_proxy = neighborhood_variation_proxy(T);
need = heuristic_mixedness_need(T);
X = [normalize_like(T.M, 3, 1), normalize_like(T.f0, 500, 150), ...
    regime_code(T), geometry_code(T), purity, 1-purity, pm, ps, ...
    log(base), fill_default(T.q_base, 0.5), log(max(qonly_sws,0.1)), qonly_q, ...
    abs(base-qonly_sws)./max(base,0.25), fill_default(T.q_base,0.5)-qonly_q, ...
    log(max(region_sws,0.1)), log(max(neighbor_sws,0.1)), log(max(blend_sws,0.1)), ...
    (base-region_sws)./max(base,0.25), (base-neighbor_sws)./max(base,0.25), ...
    grad_proxy, need];
X = double(X);
X(~isfinite(X)) = 0;
validate_predictor_policy(names);
end

function G = gate_feature_matrix(X, candidate_conf, candidate_sws, T)
delta = candidate_sws - T.sws_base;
G = [X candidate_conf delta 100*delta./max(T.sws_base,0.25) heuristic_mixedness_need(T)];
G(~isfinite(G)) = 0;
end

function [signv, conf] = direction_to_sign(direction, score, model)
classes = string(model.ClassNames);
conf = max(score, [], 2);
signv = zeros(numel(direction),1);
signv(direction == "increase_sws") = 1;
signv(direction == "decrease_sws") = -1;
for i = 1:numel(direction)
    ci = find(classes == direction(i), 1);
    if ~isempty(ci), conf(i) = score(i,ci); end
end
end

function y = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes=="true" | classes=="1", 1);
if isempty(idx), idx = size(score,2); end
y = score(:,idx);
end

function tf = class_to_logical(cls)
s = lower(string(cls));
tf = s == "true" | s == "1";
end

function y = choose_sws(base, candidate, mask)
y = base;
y(mask) = candidate(mask);
end

function y = choose_delta(delta, mask)
y = zeros(size(delta));
y(mask) = delta(mask);
end

function y = clip_sws(y, CFG)
y = min(max(y, CFG.PhysicalRangeSWS(1)), CFG.PhysicalRangeSWS(2));
end

function T = add_no_theory_candidates(T, CFG)
% Operational candidates only. No TheoryQDiscrete, true SWS, true material
% side, true purity, or distance-to-interface is used here.
base = clip_sws(T.sws_base, CFG);
qonly = clip_sws(fill_default(T.q_spectrum_only_sws, base), CFG);
T.candidate_base_sws = base;
T.candidate_q_spectrum_only_sws = qonly;
T.pseudo_region = pseudo_regions_from_base(T, CFG);
T.candidate_region_median_sws = region_median_candidate(T, base);
T.candidate_neighbor_median_sws = neighbor_median_candidate(T, base, CFG);
need = heuristic_mixedness_need(T);
T.candidate_blend_region_sws = clip_sws((1-0.5*need).*base + ...
    (0.5*need).*T.candidate_region_median_sws, CFG);
end

function [C, names] = candidate_matrix(T, CFG)
names = ["base","q_spectrum_only","region_median","neighbor_median","blend_region"];
C = [T.candidate_base_sws, T.candidate_q_spectrum_only_sws, ...
    T.candidate_region_median_sws, T.candidate_neighbor_median_sws, ...
    T.candidate_blend_region_sws];
C = clip_sws(C, CFG);
end

function [idx, conf] = candidate_label_to_index(label, score, model, names)
classes = string(model.ClassNames);
idx = ones(numel(label),1);
conf = max(score, [], 2);
for i = 1:numel(label)
    ii = find(names == label(i), 1);
    if ~isempty(ii), idx(i) = ii; end
    ci = find(classes == label(i), 1);
    if ~isempty(ci), conf(i) = score(i,ci); end
end
end

function label = candidate_direction_label(sws, base, tol)
label = repmat("base", size(sws));
label(sws > base + tol) = "increase_sws";
label(sws < base - tol) = "decrease_sws";
end

function best_idx = best_idx_for_rows(T_all, T_base, candidate_names, CFG)
% Recompute oracle candidate label for the base rows, then repeat by strategy
% block order. This is diagnostic only.
T_base = add_no_theory_candidates(T_base, CFG);
[C,~] = candidate_matrix(T_base, CFG);
err = abs(100*(C - T_base.true_SWS)./T_base.true_SWS);
[~, idx0] = min(err, [], 2);
best_idx = ones(height(T_all),1);
[ok,loc] = ismember(T_all.patch_key, T_base.patch_key);
best_idx(ok) = idx0(loc(ok));
end

function region = pseudo_regions_from_base(T, CFG)
region = ones(height(T),1);
[G,~] = findgroups(T.condition_key);
for gi = unique(G)'
    idx = find(G == gi);
    y = T.sws_base(idx);
    p = fill_default(T.predicted_patch_purity(idx), 0.95);
    pm = fill_default(T.p_mixed(idx), 1-p);
    good = isfinite(y) & p >= 0.85 & pm <= 0.65;
    if nnz(good) < 20 || range(y(good)) < 0.35
        region(idx) = 1;
        continue;
    end
    thr = median(y(good),'omitnan');
    region(idx) = 1 + double(y > thr);
end
end

function out = region_median_candidate(T, base)
out = base;
[G,~] = findgroups(T.condition_key, T.pseudo_region);
p = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-p);
good = isfinite(base) & p >= 0.85 & pm <= 0.65;
for gi = unique(G)'
    idx = find(G == gi);
    src = idx(good(idx));
    if numel(src) < 5, src = idx(isfinite(base(idx))); end
    if isempty(src), continue; end
    out(idx) = median(base(src),'omitnan');
end
out(~isfinite(out)) = base(~isfinite(out));
end

function out = neighbor_median_candidate(T, base, CFG)
out = base;
[G,~] = findgroups(T.condition_key);
p = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-p);
good_all = isfinite(base) & p >= 0.85 & pm <= 0.65;
radius = 3;
for gi = unique(G)'
    idx = find(G == gi);
    S = T(idx,:);
    for ii = 1:numel(idx)
        k = idx(ii);
        same = S.pseudo_region == T.pseudo_region(k);
        close = abs(S.map_ix - T.map_ix(k)) <= radius & abs(S.map_iz - T.map_iz(k)) <= radius;
        src_local = idx(same & close & good_all(idx));
        if numel(src_local) < 3
            src_local = idx(same & good_all(idx));
        end
        if numel(src_local) >= 3
            out(k) = median(base(src_local),'omitnan');
        end
    end
end
out(~isfinite(out)) = base(~isfinite(out));
end

%% Summaries

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
N = splitapply(@numel, T.sws_abs_error_pct, G);
MAPE_pct = splitapply(@(x) mean(x,'omitnan'), T.sws_abs_error_pct, G);
median_abs_error_pct = splitapply(@(x) median(x,'omitnan'), T.sws_abs_error_pct, G);
mean_signed_error_pct = splitapply(@(x) mean(x,'omitnan'), T.sws_signed_error_pct, G);
median_signed_error_pct = splitapply(@(x) median(x,'omitnan'), T.sws_signed_error_pct, G);
high_error10_pct = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error10, G);
high_error20_pct = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
underestimate_pct = splitapply(@(x) 100*mean(x<0,'omitnan'), T.sws_signed_error_pct, G);
overestimate_pct = splitapply(@(x) 100*mean(x>0,'omitnan'), T.sws_signed_error_pct, G);
corrected_pct = splitapply(@(x) 100*mean(x,'omitnan'), T.was_corrected, G);
kept_pct = splitapply(@(x) 100*mean(x,'omitnan'), T.was_kept, G);
mean_correction_mps = splitapply(@(x) mean(abs(x),'omitnan'), T.correction_mps, G);
median_correction_mps = splitapply(@(x) median(abs(x),'omitnan'), T.correction_mps, G);
harm_rate_corrected_pct = splitapply(@rate_among_corrected, T.is_harmed, T.was_corrected, G);
improvement_rate_corrected_pct = splitapply(@rate_among_corrected, T.is_improved, T.was_corrected, G);
net_gain_MAPE_points = splitapply(@(x) mean(x,'omitnan'), T.net_abs_error_gain_pct_points, G);
S = [groups table(N, MAPE_pct, median_abs_error_pct, mean_signed_error_pct, ...
    median_signed_error_pct, high_error10_pct, high_error20_pct, ...
    underestimate_pct, overestimate_pct, corrected_pct, kept_pct, ...
    mean_correction_mps, median_correction_mps, harm_rate_corrected_pct, ...
    improvement_rate_corrected_pct, net_gain_MAPE_points)];
end

function r = rate_among_corrected(flag, corrected)
if ~any(corrected), r = 0; else, r = 100*mean(flag(corrected),'omitnan'); end
end

function T = candidate_selector_metrics(P)
X = P(P.strategy_name == "candidate_selector",:);
if isempty(X), T = table(); return; end
[G, source_domain] = findgroups(X.source_domain);
groups = table(source_domain);
T = [groups table(splitapply(@numel, X.predicted_candidate, G), ...
    splitapply(@(a,b) 100*mean(a==b,'omitnan'), X.predicted_candidate, X.true_candidate_label, G), ...
    splitapply(@(x) mean(x,'omitnan'), X.candidate_confidence, G), ...
    'VariableNames', {'N','candidate_accuracy_pct','mean_candidate_confidence'})];
end

function T = gate_metrics(P)
X = P(P.was_corrected | ismember(P.strategy_name, ["candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector", ...
    "candidate_q_spectrum_only","candidate_region_median", ...
    "candidate_neighbor_median","candidate_blend_region","oracle_best_candidate"]),:);
if isempty(X), T = table(); return; end
[G, groups] = findgroups(X(:, {'source_domain','strategy_name'}));
T = [groups table(splitapply(@numel, X.was_corrected, G), ...
    splitapply(@(x) 100*mean(x,'omitnan'), X.was_corrected, G), ...
    splitapply(@rate_among_corrected, X.is_improved, X.was_corrected, G), ...
    splitapply(@rate_among_corrected, X.is_harmed, X.was_corrected, G), ...
    splitapply(@rate_among_corrected, X.is_unchanged_like, X.was_corrected, G), ...
    splitapply(@(x) mean(x,'omitnan'), X.net_abs_error_gain_pct_points, G), ...
    'VariableNames', {'N','corrected_pct','improved_among_corrected_pct', ...
    'harmed_among_corrected_pct','unchanged_among_corrected_pct','mean_net_gain_points'})];
end

function T = reliability_calibration(P)
X = P(ismember(P.strategy_name, ["posterior_reliability_only", ...
    "candidate_selector","conservative_candidate_selector"]),:);
if isempty(X), T = table(); return; end
edges = 0:0.1:1;
bin = discretize(X.final_high_error_probability, edges, 'IncludedEdge','right');
labels = strings(size(bin));
for i = 1:numel(bin)
    if isfinite(bin(i))
        labels(i) = sprintf('%.1f-%.1f', edges(bin(i)), edges(bin(i)+1));
    else
        labels(i) = "unknown";
    end
end
X.reliability_bin = labels;
[G, groups] = findgroups(X(:, {'source_domain','strategy_name','reliability_bin'}));
T = [groups table(splitapply(@numel, X.high_error20, G), ...
    splitapply(@(x) mean(x,'omitnan'), X.final_high_error_probability, G), ...
    splitapply(@(x) 100*mean(x,'omitnan'), X.high_error20, G), ...
    'VariableNames', {'N','mean_predicted_high_error_probability', ...
    'observed_high_error20_pct'})];
end

function T = best_strategy_candidates(T_by_source, T_by_region, CFG)
T = table();
for src = unique(T_by_source.source_domain, 'stable')'
    X = T_by_source(T_by_source.source_domain == src,:);
    X = X(~ismember(X.strategy_name, ["oracle_best_candidate"]),:);
    if isempty(X), continue; end
    [~,i1] = min(X.MAPE_pct);
    [~,i2] = min(X.high_error20_pct);
    R = table(src, X.strategy_name(i1), X.MAPE_pct(i1), X.strategy_name(i2), X.high_error20_pct(i2), ...
        'VariableNames', {'source_domain','best_MAPE_strategy','best_MAPE_pct', ...
        'best_high20_strategy','best_high20_pct'});
    T = concat_tables(T, R);
end
if isempty(T_by_region), return; end
base = T_by_region(T_by_region.strategy_name == "base_q_spectrum_plus_composition",:);
main = T_by_region(T_by_region.strategy_name == "conservative_candidate_selector",:);
[ok,loc] = ismember(main(:, {'source_domain','evaluation_region'}), base(:, {'source_domain','evaluation_region'}));
if any(ok)
    main.conservative_gain_vs_base_points = base.MAPE_pct(loc(ok)) - main.MAPE_pct(ok);
end
end

function T = reference_comparison_table(root_dir, T45_source)
T = table();
% Current Test45 rows.
keep45 = ["base_q_spectrum_plus_composition","conservative_candidate_selector", ...
    "candidate_selector","mixedness_candidate_selector","oracle_best_candidate"];
for i = 1:height(T45_source)
    if ~ismember(T45_source.strategy_name(i), keep45), continue; end
    R = table("Test45 " + pretty_source_name(T45_source.source_domain(i)), ...
        T45_source.strategy_name(i), T45_source.MAPE_pct(i), T45_source.high_error20_pct(i), ...
        "test45", 'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
% Test33 synthetic reference.
p33 = fullfile(root_dir, 'outputs', 'test_33_mixedness_aware_q_correction', ...
    'tables', 'test33_strategy_summary_overall_test.csv');
if exist(p33,'file') == 2
    A = readtable(p33);
    A.strategy_name = string(A.strategy_name);
    keep = ["local_baseline","sws_nearest_highconf","test30_theory_region_levels", ...
        "mixedness_logk_corrected","mixedness_q_candidate_selector"];
    A = A(ismember(A.strategy_name, keep),:);
    R = table(repmat("Test33 synthetic reference",height(A),1), A.strategy_name, ...
        A.MAPE, A.high_error_20_pct, repmat("test33",height(A),1), ...
        'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
% Test34 k-Wave reference.
p34 = fullfile(root_dir, 'outputs', 'test_34_kwave_mixedness_transfer', ...
    'tables', 'test34_kwave_strategy_summary_overall.csv');
if exist(p34,'file') == 2
    A = readtable(p34);
    A.strategy_name = string(A.strategy_name);
    keep = ["local_baseline","theory_baseline","test30_theory_region_levels", ...
        "mixedness_logk_corrected","mixedness_q_candidate_selector"];
    A = A(ismember(A.strategy_name, keep),:);
    R = table(repmat("Test34 k-Wave reference",height(A),1), A.strategy_name, ...
        A.MAPE, A.high_error_20_pct, repmat("test34",height(A),1), ...
        'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
end

%% Plots

function plot_strategy_ranking(T, OUT, metric, filename)
T = T(ismember(T.strategy_name, ["base_q_spectrum_plus_composition", ...
    "candidate_q_spectrum_only","candidate_region_median", ...
    "candidate_neighbor_median","candidate_blend_region", ...
    "candidate_selector","conservative_candidate_selector", ...
    "mixedness_candidate_selector","oracle_best_candidate", ...
    "posterior_reliability_only"]),:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 16]);
sources = unique(T.source_domain, 'stable');
tl = tiledlayout(fig,1,numel(sources),'TileSpacing','compact','Padding','compact');
for s = sources(:)'
    ax = nexttile(tl);
    X = sortrows(T(T.source_domain == s,:), metric, 'ascend');
    barh(ax, categorical(arrayfun(@pretty_strategy_name, X.strategy_name)), X.(metric));
    xlabel(ax, pretty_metric_name(metric));
    title(ax, pretty_source_name(s), 'FontWeight','normal');
    grid(ax,'on'); set(ax,'TickLabelInterpreter','none');
end
title(tl, ternary(metric=="MAPE_pct", "Strategy ranking by source", ...
    "High-error rate by source"), 'FontWeight','bold');
export_fig(fig, fullfile(OUT.figure_dir, filename));
end

function plot_reference_comparison(T, OUT)
if isempty(T), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 16]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for panel = 1:2
    ax = nexttile(tl);
    if panel == 1
        mask = contains(T.benchmark, "synthetic", 'IgnoreCase', true) | ...
            contains(T.benchmark, "OOD", 'IgnoreCase', true);
        X = T(mask,:);
        ttl = 'Synthetic / OOD comparison';
    else
        mask = contains(T.benchmark, "k-Wave", 'IgnoreCase', true) | ...
            contains(T.benchmark, "kwave", 'IgnoreCase', true);
        X = T(mask,:);
        ttl = 'k-Wave comparison';
    end
    if isempty(X), continue; end
    X.label = X.benchmark + " | " + arrayfun(@pretty_strategy_name, X.strategy_name);
    X = sortrows(X, 'MAPE_pct', 'ascend');
    barh(ax, categorical(X.label), X.MAPE_pct);
    xlabel(ax, 'MAPE (%)'); title(ax, ttl, 'FontWeight','normal');
    grid(ax,'on'); set(ax,'TickLabelInterpreter','none');
end
title(tl, 'Reference comparison with Test33/Test34', 'FontWeight','bold');
export_fig(fig, fullfile(OUT.figure_dir, 'test45_reference_comparison_test33_test34.png'));
end

function plot_correction_benefit_harm(T, OUT)
if isempty(T), return; end
keep = ["candidate_q_spectrum_only","candidate_region_median", ...
    "candidate_neighbor_median","candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector", ...
    "oracle_best_candidate"];
T = T(ismember(T.strategy_name, keep),:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 16]);
sources = unique(T.source_domain, 'stable');
tl = tiledlayout(fig,1,numel(sources),'TileSpacing','compact','Padding','compact');
for s = sources(:)'
    ax = nexttile(tl);
    X = T(T.source_domain == s,:);
    Y = [X.improved_among_corrected_pct X.harmed_among_corrected_pct X.unchanged_among_corrected_pct];
    bar(ax, categorical(arrayfun(@pretty_strategy_name, X.strategy_name)), Y);
    ylabel(ax, 'Corrected pixels (%)');
    title(ax, pretty_source_name(s), 'FontWeight','normal');
    legend(ax, ["Improved","Harmed","Approximately unchanged"], 'Location','bestoutside');
    grid(ax,'on'); xtickangle(ax,25); set(ax,'TickLabelInterpreter','none');
end
title(tl, 'Correction benefit and harm among corrected pixels', 'FontWeight','bold');
export_fig(fig, fullfile(OUT.figure_dir, 'test45_correction_benefit_harm.png'));
end

function plot_direction_confusion(P, OUT)
X = P(P.strategy_name == "candidate_selector",:);
if isempty(X), return; end
for src = unique(X.source_domain, 'stable')'
    S = X(X.source_domain == src,:);
    classes = ["base","q_spectrum_only","region_median","neighbor_median","blend_region"];
    C = zeros(numel(classes));
    for i = 1:numel(classes)
        idx = S.true_candidate_label == classes(i);
        denom = max(1, nnz(idx));
        for j = 1:numel(classes)
            C(i,j) = 100*nnz(idx & S.predicted_candidate == classes(j))/denom;
        end
    end
    fig = figure('Color','w','Units','centimeters','Position',[1 1 20 16]);
    imagesc(C); axis image; colormap(parula);
    cb = colorbar; ylabel(cb, 'Fraction of true class (%)');
    class_labels = strings(size(classes));
    for ci = 1:numel(classes), class_labels(ci) = pretty_candidate_name(classes(ci)); end
    set(gca,'XTick',1:numel(classes),'XTickLabel',cellstr(class_labels), ...
        'YTick',1:numel(classes),'YTickLabel',cellstr(class_labels), ...
        'TickLabelInterpreter','none');
    xtickangle(35);
    xlabel('Predicted candidate'); ylabel('Best candidate by true error');
    title("Candidate selector confusion matrix: " + pretty_source_name(src), 'FontWeight','normal');
    for i=1:numel(classes), for j=1:numel(classes), text(j,i,sprintf('%.1f',C(i,j)), ...
            'HorizontalAlignment','center','Color','w','FontWeight','bold'); end, end
    export_fig(fig, fullfile(OUT.figure_dir, "test45_candidate_confusion_matrix_" + sanitize(src) + ".png"));
end
end

function plot_error_before_after(P, OUT, signed_mode)
strategy = "conservative_candidate_selector";
X = P(P.strategy_name == strategy,:);
if height(X) > 20000
    X = X(randperm(height(X), 20000),:);
end
for src = unique(X.source_domain, 'stable')'
    S = X(X.source_domain == src,:);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 16 14]);
    if signed_mode
        scatter(S.base_signed_error_pct, S.sws_signed_error_pct, 8, S.p_mixed, 'filled', 'MarkerFaceAlpha',0.35);
        xlabel('Baseline signed SWS error (%)'); ylabel('Corrected signed SWS error (%)');
        title("Signed error before and after correction: " + pretty_source_name(src), 'FontWeight','normal');
        xline(0,'k-'); yline(0,'k-');
        file = "test45_signed_error_before_after_" + sanitize(src) + ".png";
    else
        scatter(S.base_abs_error_pct, S.sws_abs_error_pct, 8, S.p_strong_mixed, 'filled', 'MarkerFaceAlpha',0.35);
        xlabel('Baseline absolute SWS error (%)'); ylabel('Corrected absolute SWS error (%)');
        title("Absolute error before and after correction: " + pretty_source_name(src), 'FontWeight','normal');
        file = "test45_absolute_error_before_after_" + sanitize(src) + ".png";
    end
    hold on; lim = axis; lo = min(lim([1 3])); hi = max(lim([2 4])); plot([lo hi],[lo hi],'k--'); axis([lo hi lo hi]);
    cb = colorbar; ylabel(cb, ternary(signed_mode, 'Predicted mixedness probability', 'Predicted strong-mixed probability'));
    grid on;
    export_fig(fig, fullfile(OUT.figure_dir, file));
end
end

function plot_error_by_region(T, OUT)
keep = ["base_q_spectrum_plus_composition","candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector", ...
    "candidate_neighbor_median","oracle_best_candidate"];
T = T(ismember(T.strategy_name, keep),:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
sources = unique(T.source_domain, 'stable');
tl = tiledlayout(fig,1,numel(sources),'TileSpacing','compact','Padding','compact');
for src = sources(:)'
    ax = nexttile(tl);
    X = T(T.source_domain == src,:);
    grouped_bar(ax, X, "evaluation_region", keep, "MAPE_pct");
    ylabel(ax,'MAPE (%)'); xlabel(ax,'Evaluation region');
    title(ax, pretty_source_name(src), 'FontWeight','normal'); grid(ax,'on');
    legend(ax, arrayfun(@pretty_strategy_name, keep), 'Location','bestoutside');
    xtickangle(ax,25);
end
title(tl, 'Regional SWS error by strategy', 'FontWeight','bold');
export_fig(fig, fullfile(OUT.figure_dir, 'test45_error_by_region.png'));
end

function plot_error_vs_distance(T, OUT)
keep = ["base_q_spectrum_plus_composition","candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector"];
plot_grouped_lines(T, "distance_bin_pretty", keep, "MAPE_pct", ...
    'Distance to interface (mm)', 'MAPE (%)', 'Error versus distance to interface', ...
    fullfile(OUT.figure_dir, 'test45_error_vs_distance.png'));
end

function plot_error_vs_predicted_mixedness(T, OUT)
keep = ["base_q_spectrum_plus_composition","candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector"];
plot_grouped_lines(T, "predicted_mixedness_bin", keep, "MAPE_pct", ...
    'Predicted mixedness probability bin', 'MAPE (%)', 'Error versus predicted mixedness', ...
    fullfile(OUT.figure_dir, 'test45_error_vs_predicted_mixedness.png'));
end

function plot_grouped_lines(T, xvar, strategies, metric, xlabel_txt, ylabel_txt, ttl, file)
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 16]);
sources = unique(T.source_domain, 'stable');
tl = tiledlayout(fig,1,numel(sources),'TileSpacing','compact','Padding','compact');
for src = sources(:)'
    ax = nexttile(tl); hold(ax,'on');
    X = T(T.source_domain == src,:);
    cats = unique(string(X.(xvar)), 'stable');
    for s = strategies
        y = nan(numel(cats),1);
        for i = 1:numel(cats)
            idx = X.strategy_name == s & string(X.(xvar)) == cats(i);
            if any(idx), y(i) = mean(X.(metric)(idx), 'omitnan'); end
        end
        plot(ax, 1:numel(cats), y, '-o', 'LineWidth',1.2, 'MarkerSize',4, ...
            'DisplayName', pretty_strategy_name(s));
    end
    set(ax,'XTick',1:numel(cats),'XTickLabel',cats,'TickLabelInterpreter','none');
    xtickangle(ax,25); xlabel(ax,xlabel_txt); ylabel(ax,ylabel_txt);
    title(ax, pretty_source_name(src), 'FontWeight','normal'); grid(ax,'on');
    legend(ax,'Location','bestoutside');
end
title(tl, ttl, 'FontWeight','bold');
export_fig(fig, file);
end

function plot_correction_magnitude_vs_mixedness(P, OUT)
X = P(ismember(P.strategy_name, ["candidate_selector", ...
    "conservative_candidate_selector","mixedness_candidate_selector"]),:);
X = X(X.was_corrected,:);
if isempty(X), return; end
if height(X) > 20000, X = X(randperm(height(X), 20000),:); end
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 16]);
scatter(X.p_mixed, abs(X.correction_pct), 8, X.final_high_error_probability, 'filled', 'MarkerFaceAlpha',0.35);
xlabel('Predicted mixedness probability'); ylabel('Correction magnitude in SWS (%)');
title('Correction magnitude versus predicted mixedness', 'FontWeight','normal');
cb = colorbar; ylabel(cb, 'Final high-error probability');
grid on;
export_fig(fig, fullfile(OUT.figure_dir, 'test45_correction_magnitude_vs_mixedness.png'));
end

function plot_reliability_calibration(T, OUT)
if isempty(T), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 26 16]);
sources = unique(T.source_domain, 'stable');
tl = tiledlayout(fig,1,numel(sources),'TileSpacing','compact','Padding','compact');
for src = sources(:)'
    ax = nexttile(tl); hold(ax,'on');
    X = T(T.source_domain == src,:);
    strategies = unique(X.strategy_name, 'stable');
    for s = strategies(:)'
        S = X(X.strategy_name == s,:);
        plot(ax, S.mean_predicted_high_error_probability, S.observed_high_error20_pct/100, ...
            '-o', 'DisplayName', pretty_strategy_name(s), 'LineWidth',1.1);
    end
    plot(ax,[0 1],[0 1],'k--','DisplayName','Ideal calibration');
    xlabel(ax,'Predicted high-error probability'); ylabel(ax,'Observed high-error rate');
    title(ax, pretty_source_name(src), 'FontWeight','normal'); grid(ax,'on');
    legend(ax,'Location','bestoutside');
end
title(tl, 'Reliability calibration', 'FontWeight','bold');
export_fig(fig, fullfile(OUT.figure_dir, 'test45_reliability_calibration.png'));
end

function plot_frequency_dependence(T, OUT, M)
X = T(T.M == M,:);
if isempty(X), return; end
keep = ["base_q_spectrum_plus_composition","conservative_candidate_selector", ...
    "candidate_selector","mixedness_candidate_selector"];
plot_grouped_lines(X, "f0", keep, "MAPE_pct", 'Frequency (Hz)', 'MAPE (%)', ...
    sprintf('Frequency dependence of SWS error, M=%d', M), ...
    fullfile(OUT.figure_dir, sprintf('test45_frequency_dependence_M%d.png', M)));
end

function plot_net_gain_heatmap(T, OUT, CFG)
base = T(T.strategy_name == "base_q_spectrum_plus_composition",:);
X = T(ismember(T.strategy_name, ["conservative_candidate_selector","mixedness_candidate_selector", ...
    "candidate_selector","candidate_neighbor_median"]),:);
if isempty(X) || isempty(base), return; end
X.gain = nan(height(X),1);
[ok,loc] = ismember(X(:, {'source_domain','evaluation_region'}), base(:, {'source_domain','evaluation_region'}));
X.gain(ok) = base.MAPE_pct(loc(ok)) - X.MAPE_pct(ok);
rows = unique(X.strategy_name,'stable');
cols = unique(X.source_domain + " | " + X.evaluation_region,'stable');
Z = nan(numel(rows), numel(cols));
for i = 1:numel(rows)
    for j = 1:numel(cols)
        parts = split(cols(j), " | ");
        idx = X.strategy_name == rows(i) & X.source_domain == parts(1) & X.evaluation_region == parts(2);
        if any(idx), Z(i,j) = mean(X.gain(idx),'omitnan'); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 18]);
imagesc(Z); colormap(redblue()); cb = colorbar; ylabel(cb, 'MAPE gain versus baseline (percentage points)');
set(gca,'XTick',1:numel(cols),'XTickLabel',arrayfun(@pretty_combo_name, cols), ...
    'YTick',1:numel(rows),'YTickLabel',arrayfun(@pretty_strategy_name, rows), ...
    'TickLabelInterpreter','none');
xtickangle(35);
title('Net correction gain versus baseline: positive means improvement', 'FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test45_net_gain_heatmap.png'));
end

function plot_condition_maps(P, CFG, OUT)
strategy = "conservative_candidate_selector";
conds = unique(P(:, {'source_domain','case_id','field_regime_ood','f0','M','condition_key'}), 'rows', 'stable');
if isfinite(CFG.MaxMapConditions) && height(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions,:);
end
fprintf('Saving %d Test45 map panels under %s.\n', height(conds), OUT.map_dir);
for ci = 1:height(conds)
    idx = P.source_domain == conds.source_domain(ci) & P.condition_key == conds.condition_key(ci) & ...
        P.M == conds.M(ci);
    Xc = P(idx,:);
    B = Xc(Xc.strategy_name == "base_q_spectrum_plus_composition",:);
    C = Xc(Xc.strategy_name == strategy,:);
    if isempty(B) || isempty(C), continue; end
    [true_map,nz,nx] = rows_to_diagnostic_grid(B, B.true_SWS, CFG, "nearest");
    cand_num = candidate_numeric(C.predicted_candidate);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 42 34]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(B, B.sws_base, CFG, "natural", nz, nx), 'Base predicted SWS', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.sws_pred_final, CFG, "natural", nz, nx), 'Corrected SWS', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.candidate_q_spectrum_only_sws, CFG, "natural", nz, nx), 'q-spectrum-only candidate', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.candidate_region_median_sws, CFG, "natural", nz, nx), 'Region-median candidate', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.candidate_neighbor_median_sws, CFG, "natural", nz, nx), 'Neighbor-median candidate', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.pseudo_region, CFG, "nearest", nz, nx), 'Pseudo-region from base SWS', 'Region label');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.predicted_patch_purity, CFG, "natural", nz, nx), 'Predicted patch purity', 'Probability');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.p_mixed, CFG, "natural", nz, nx), 'Predicted mixedness', 'Probability');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(B, B.base_signed_error_pct, CFG, "natural", nz, nx), 'Baseline signed error', 'Signed error (%)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.sws_signed_error_pct, CFG, "natural", nz, nx), 'Corrected signed error', 'Signed error (%)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.net_abs_error_gain_pct_points, CFG, "natural", nz, nx), 'Absolute-error improvement', 'Improvement (%)');
    plot_candidate_map(nexttile(tl), rows_to_diagnostic_grid(C, cand_num, CFG, "nearest", nz, nx));
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, double(C.was_corrected), CFG, "nearest", nz, nx), 'Correction applied mask', 'Applied (0/1)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.correction_mps, CFG, "natural", nz, nx), 'Correction magnitude', 'SWS change (m/s)');
    plot_map(nexttile(tl), rows_to_diagnostic_grid(C, C.final_high_error_probability, CFG, "natural", nz, nx), 'Final high-error probability', 'Probability');
    title(tl, sprintf('Case: %s, f=%d Hz, M=%d, regime=%s, source=%s', ...
        conds.case_id(ci), conds.f0(ci), conds.M(ci), conds.field_regime_ood(ci), ...
        pretty_source_name(conds.source_domain(ci))), 'Interpreter','none');
    outdir = fullfile(OUT.map_dir, sanitize(conds.source_domain(ci)), sanitize(conds.case_id(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    export_fig(fig, fullfile(outdir, sprintf('test45_map_%s__%s__%s__f%d__M%d.png', ...
        sanitize(CFG.MapStyle), sanitize(conds.source_domain(ci)), sanitize(conds.case_id(ci)), conds.f0(ci), conds.M(ci))));
end
end

function grouped_bar(ax, T, group_var, strategies, metric)
groups = unique(string(T.(group_var)), 'stable');
Y = nan(numel(groups), numel(strategies));
for i = 1:numel(groups)
    for j = 1:numel(strategies)
        idx = string(T.(group_var)) == groups(i) & T.strategy_name == strategies(j);
        if any(idx), Y(i,j) = mean(T.(metric)(idx), 'omitnan'); end
    end
end
bar(ax, categorical(arrayfun(@pretty_region_name, groups)), Y);
set(ax,'TickLabelInterpreter','none');
end

function plot_map(ax, Z, ttl, cb_label)
imagesc(ax, Z); axis(ax,'image'); axis(ax,'off');
cb = colorbar(ax); ylabel(cb, cb_label);
title(ax, ttl, 'FontWeight','normal');
end

function plot_direction_map(ax, Z)
imagesc(ax, Z, [-1 1]); axis(ax,'image'); axis(ax,'off');
colormap(ax, [0.20 0.35 0.90; 0.82 0.82 0.82; 0.90 0.25 0.20]);
cb = colorbar(ax, 'Ticks', [-1 0 1], 'TickLabels', {'Decrease SWS','Keep','Increase SWS'});
ylabel(cb, 'Correction direction');
title(ax, 'Predicted correction direction', 'FontWeight','normal');
end

function plot_candidate_map(ax, Z)
imagesc(ax, Z, [1 5]); axis(ax,'image'); axis(ax,'off');
colormap(ax, lines(5));
cb = colorbar(ax, 'Ticks', 1:5, 'TickLabels', {'Base','q-spectrum only','Region median','Neighbor median','Blend'});
ylabel(cb, 'Selected candidate');
title(ax, 'Predicted correction candidate', 'FontWeight','normal');
end

%% Labels and bins

function s = pretty_strategy_name(name)
name = string(name);
switch name
    case "base_q_spectrum_plus_composition", s = "Base: q spectrum + composition";
    case "candidate_q_spectrum_only", s = "Candidate: q spectrum only";
    case "candidate_region_median", s = "Candidate: region median";
    case "candidate_neighbor_median", s = "Candidate: neighbor median";
    case "candidate_blend_region", s = "Candidate: base-region blend";
    case "candidate_selector", s = "Candidate selector";
    case "conservative_candidate_selector", s = "Conservative candidate selector";
    case "mixedness_candidate_selector", s = "Mixedness-gated candidate selector";
    case "oracle_best_candidate", s = "Oracle best candidate";
    case "posterior_reliability_only", s = "Reliability only";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_region_name(name)
name = string(name);
switch name
    case "homogeneous", s = "Homogeneous";
    case "pure", s = "Pure patches";
    case "near_pure", s = "Near-pure patches";
    case "moderately_mixed", s = "Moderately mixed patches";
    case "strongly_mixed", s = "Strongly mixed patches";
    case "soft", s = "Soft/background";
    case "hard", s = "Hard/inclusion";
    case "interface", s = "Interface/mixed";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_source_name(name)
name = string(name);
switch name
    case "synthetic_heldout", s = "Synthetic held-out";
    case "test39_ood", s = "Synthetic OOD";
    case "kwave_transfer", s = "k-Wave transfer";
    case "synthetic_train_source", s = "Synthetic train source";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_regime_name(name)
s = replace(string(name), "_", " ");
s = replace(s, "directional 2D", "Directional 2D");
s = replace(s, "diffuse 2D", "Diffuse 2D");
s = replace(s, "diffuse 3D", "Diffuse 3D");
s = replace(s, "partial 3D", "Partial 3D");
end

function s = pretty_metric_name(name)
name = string(name);
switch name
    case "MAPE_pct", s = "MAPE (%)";
    case "high_error20_pct", s = "Pixels with error >20% (%)";
    case "mean_signed_error_pct", s = "Mean signed SWS error (%)";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_direction_name(name)
name = string(name);
switch name
    case "increase_sws", s = "Increase SWS";
    case "decrease_sws", s = "Decrease SWS";
    case "keep", s = "Keep";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_candidate_name(name)
name = string(name);
switch name
    case "base", s = "Base";
    case "q_spectrum_only", s = "q-spectrum only";
    case "region_median", s = "Region median";
    case "neighbor_median", s = "Neighbor median";
    case "blend_region", s = "Base-region blend";
    otherwise, s = replace(name, "_", " ");
end
end

function s = pretty_combo_name(name)
parts = split(string(name), " | ");
if numel(parts) >= 2
    s = pretty_source_name(parts(1)) + " | " + pretty_region_name(parts(2));
else
    s = replace(string(name), "_", " ");
end
end

function region = evaluation_region(T)
region = strings(height(T),1);
region(:) = "pure";
hom = contains(T.case_family, "homogeneous");
region(hom) = "homogeneous";
if ismember("patch_purity", string(T.Properties.VariableNames))
    p = T.patch_purity;
    region(~hom & p >= 0.99) = "pure";
    region(~hom & p < 0.99 & p >= 0.90) = "near_pure";
    region(~hom & p < 0.90 & p >= 0.70) = "moderately_mixed";
    region(~hom & p < 0.70) = "strongly_mixed";
end
if ismember("region_zone", string(T.Properties.VariableNames))
    rz = string(T.region_zone);
    region(contains(rz, "interface")) = "interface";
end
end

function bin = true_purity_bin(p, family)
bin = strings(size(p));
hom = contains(string(family), "homogeneous");
bin(hom) = "homogeneous";
bin(~hom & p >= 0.99) = "pure";
bin(~hom & p < 0.99 & p >= 0.90) = "near_pure";
bin(~hom & p < 0.90 & p >= 0.70) = "moderately_mixed";
bin(~hom & p < 0.70) = "strongly_mixed";
bin(bin=="") = "unknown";
end

function bin = probability_bin(p)
edges = [0 0.2 0.4 0.6 0.8 1.00001];
bin_id = discretize(p, edges);
labels = ["0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"];
bin = strings(size(p));
for i = 1:numel(p)
    if isfinite(bin_id(i)), bin(i) = labels(bin_id(i)); else, bin(i) = "unknown"; end
end
end

function bin = distance_bin_pretty(d, family)
bin = strings(size(d));
hom = contains(string(family), "homogeneous");
bin(hom) = "homogeneous";
bin(~hom & d < 1) = "0-1 mm";
bin(~hom & d >= 1 & d < 2) = "1-2 mm";
bin(~hom & d >= 2 & d < 4) = "2-4 mm";
bin(~hom & d >= 4 & d < 8) = "4-8 mm";
bin(~hom & d >= 8) = ">8 mm";
bin(bin=="") = "unknown";
end

function z = pretty_distance_zone(d)
z = strings(size(d));
z(d < 2) = "interface_0_2mm";
z(d >= 2 & d < 4) = "near_interface_2_4mm";
z(d >= 4) = "core_gt4mm";
z(~isfinite(d)) = "unknown";
end

function label = direction_label_from_error(signed_error, CFG)
label = repmat("keep", size(signed_error));
label(signed_error < -CFG.DirectionThresholdPct) = "increase_sws";
label(signed_error > CFG.DirectionThresholdPct) = "decrease_sws";
end

function y = direction_numeric(direction)
y = zeros(size(direction));
y(string(direction) == "increase_sws") = 1;
y(string(direction) == "decrease_sws") = -1;
end

function y = candidate_numeric(candidate)
y = ones(size(candidate));
y(string(candidate) == "q_spectrum_only") = 2;
y(string(candidate) == "region_median") = 3;
y(string(candidate) == "neighbor_median") = 4;
y(string(candidate) == "blend_region") = 5;
end

%% Feature helpers

function need = heuristic_mixedness_need(T)
purity = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-purity);
ps = fill_default(T.p_strong_mixed, 0);
var_dis = abs(fill_default(T.variant_sws_disagreement, 0));
grad = neighborhood_variation_proxy(T);
need = max([1-purity, pm, 0.75*ps, min(var_dis,1), min(grad,1)], [], 2);
need = min(max(need,0),1);
end

function g = neighborhood_variation_proxy(T)
g = zeros(height(T),1);
if ~all(ismember(["condition_key","map_iz","map_ix","sws_base"], string(T.Properties.VariableNames)))
    return;
end
[G,~] = findgroups(T.condition_key);
for gi = unique(G)'
    idx = find(G == gi);
    S = T(idx,:);
    nz = max(S.map_iz); nx = max(S.map_ix);
    Z = nan(nz,nx);
    Z(sub2ind([nz,nx], S.map_iz, S.map_ix)) = S.sws_base;
    K = ones(3); K(2,2) = 0;
    cnt = conv2(double(isfinite(Z)), K, 'same');
    sumv = conv2(fillmissing2(Z), K, 'same');
    mu = sumv ./ max(cnt,1);
    gv = abs(Z - mu) ./ max(Z,0.25);
    g(idx) = gv(sub2ind([nz,nx], S.map_iz, S.map_ix));
end
g(~isfinite(g)) = 0;
end

function Z = fillmissing2(Z)
bad = ~isfinite(Z);
if any(~bad,'all')
    Z(bad) = median(Z(~bad),'omitnan');
else
    Z(:) = 0;
end
end

function c = regime_code(T)
r = string(T.field_regime_ood);
c = zeros(numel(r),1);
c(contains(lower(r),"directional")) = 1;
c(contains(lower(r),"diffuse_2d")) = 2;
c(contains(lower(r),"partial")) = 3;
c(contains(lower(r),"diffuse_3d")) = 4;
c = normalize_like(c, 2.5, 1.2);
end

function c = geometry_code(T)
g = lower(string(T.case_family));
c = zeros(numel(g),1);
c(contains(g,"homogeneous")) = 1;
c(contains(g,"bilayer")) = 2;
c(contains(g,"inclusion")) = 3;
c(contains(g,"kwave")) = 4;
c = normalize_like(c, 2.5, 1.2);
end

function validate_predictor_policy(names)
bad = ["true","oracle","material","side","roi","distance","error", ...
    "signed","abs_error","high_error","cs_true","sws_true","k_true", ...
    "patch_purity_true","purity_bin","region_zone"];
names = lower(string(names));
for b = bad
    hit = names(contains(names, b));
    assert(isempty(hit), 'Oracle/evaluation feature leaked into Test45 predictors: %s', strjoin(hit, ', '));
end
end

%% Console interpretation

function print_console_summary(T_overall, T_by_source, T_by_region, T_gate, T_best, CFG)
fprintf('\n================ Test 45 summary ================\n');
X = T_overall(~ismember(T_overall.strategy_name, ["oracle_best_candidate"]),:);
X = sortrows(X, 'MAPE_pct', 'ascend');
fprintf('Best global non-oracle strategy: %s (MAPE %.2f%%, HE20 %.2f%%).\n', ...
    pretty_strategy_name(X.strategy_name(1)), X.MAPE_pct(1), X.high_error20_pct(1));
for src = unique(T_by_source.source_domain,'stable')'
    S = T_by_source(T_by_source.source_domain == src & ...
        ~ismember(T_by_source.strategy_name, ["oracle_best_candidate"]),:);
    S = sortrows(S, 'MAPE_pct', 'ascend');
    fprintf('Best on %s: %s (MAPE %.2f%%, HE20 %.2f%%).\n', ...
        pretty_source_name(src), pretty_strategy_name(S.strategy_name(1)), S.MAPE_pct(1), S.high_error20_pct(1));
end
base = T_by_region(T_by_region.strategy_name == "base_q_spectrum_plus_composition",:);
cons = T_by_region(T_by_region.strategy_name == "conservative_candidate_selector",:);
if ~isempty(base) && ~isempty(cons)
    [ok,loc] = ismember(cons(:, {'source_domain','evaluation_region'}), base(:, {'source_domain','evaluation_region'}));
    gains = base.MAPE_pct(loc(ok)) - cons.MAPE_pct(ok);
    regs = cons.evaluation_region(ok);
    srcs = cons.source_domain(ok);
    for target = ["interface","strongly_mixed","homogeneous","pure","near_pure"]
        idx = regs == target;
        if any(idx)
            fprintf('Conservative correction gain on %s: %.2f MAPE points averaged over sources.\n', ...
                pretty_region_name(target), mean(gains(idx),'omitnan'));
        end
    end
end
if ~isempty(T_gate)
    G = T_gate(T_gate.strategy_name == "conservative_candidate_selector",:);
    if ~isempty(G)
        fprintf('Conservative correction corrected %.2f%% of pixels on average; harm among corrected %.2f%%.\n', ...
            mean(G.corrected_pct,'omitnan'), mean(G.harmed_among_corrected_pct,'omitnan'));
    end
end
fprintf('Recommendation rule: use correction only if it maintains global/pure/soft regions and improves mixed/interface. Otherwise use baseline plus reliability.\n');
fprintf('=================================================\n');
end

%% Utilities

function y = fill_default(x, d)
y = x;
bad = ~isfinite(y);
if isscalar(d)
    y(bad) = d;
else
    y(bad) = d(bad);
end
end

function y = normalize_like(x, mu, sigma)
y = (double(x)-mu)./sigma;
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
valid = isfinite(T.map_iz) & isfinite(T.map_ix) & T.map_iz>=1 & T.map_ix>=1;
Z(sub2ind([nz,nx], T.map_iz(valid), T.map_ix(valid))) = values(valid);
end

function [Z,nz,nx] = rows_to_diagnostic_grid(T, values, CFG, method, nz, nx)
if nargin < 5
    if CFG.MapStyle == "interp"
        nz = max(2, ceil(max(T.map_iz) * CFG.MapInterpScale));
        nx = max(2, ceil(max(T.map_ix) * CFG.MapInterpScale));
    else
        nz = max(T.map_iz);
        nx = max(T.map_ix);
    end
end
if CFG.MapStyle == "patch"
    Z = rows_to_grid(T, values, nz, nx);
    return;
end
x = double(T.map_ix);
z = double(T.map_iz);
v = double(values);
valid = isfinite(x) & isfinite(z) & isfinite(v);
if nnz(valid) < 3
    Z = nan(nz,nx);
    return;
end
xq = linspace(min(x(valid)), max(x(valid)), nx);
zq = linspace(min(z(valid)), max(z(valid)), nz);
[Xq,Zq] = meshgrid(xq,zq);
try
    F = scatteredInterpolant(x(valid), z(valid), v(valid), char(method), 'nearest');
catch
    F = scatteredInterpolant(x(valid), z(valid), v(valid), 'nearest', 'nearest');
end
Z = F(Xq,Zq);
end

function T = vertcat_compatible(varargin)
vars = strings(1,0);
for i = 1:nargin
    vars = [vars string(varargin{i}.Properties.VariableNames)]; %#ok<AGROW>
end
vars = unique(vars,'stable');
string_vars = ["dataset","source_domain","condition_key","case_id","case_family","case_name", ...
    "field_regime","field_regime_ood","purity_bin","distance_bin","material_region", ...
    "region_zone","material_zone","model_name","model_family","strategy_name", ...
    "strategy_name_original","strategy_type","patch_key","geometry","geometry_type", ...
    "material_side","analysis_region","roi_name","region_label","bundle_id", ...
    "evaluation_region","true_purity_bin","predicted_mixedness_bin","distance_bin_pretty", ...
    "predicted_direction","predicted_candidate","true_candidate_label"];
logical_vars = ["is_train_row","is_heldout_row","is_mixed","is_strong_mixed", ...
    "high_error10","high_error20","was_corrected","was_kept","is_harmed", ...
    "is_improved","is_unchanged_like"];
parts = cell(nargin,1);
for i = 1:nargin
    A = varargin{i};
    for v = vars
        if ismember(v,string(A.Properties.VariableNames)), continue; end
        if ismember(v,string_vars)
            A.(v) = repmat("unknown",height(A),1);
        elseif ismember(v,logical_vars)
            A.(v) = false(height(A),1);
        else
            A.(v) = nan(height(A),1);
        end
    end
    parts{i} = A(:,cellstr(vars));
end
T = vertcat(parts{:});
end

function T = concat_tables(T, R)
if isempty(T), T = R; else, T = vertcat_compatible(T, R); end
end

function export_fig(fig, file)
drawnow;
exportgraphics(fig, file, 'Resolution', 220);
close(fig);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Skipping plot %s: %s', label, ME.message);
end
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
end

function tf = env_true(name, default_value)
v = lower(strtrim(string(getenv(name))));
if v == ""
    tf = default_value;
else
    tf = ismember(v, ["1","true","yes","y","on"]);
end
end

function x = env_number(name, default_value)
v = strtrim(string(getenv(name)));
if v == ""
    x = default_value;
else
    x = str2double(v);
    if ~isfinite(x), x = default_value; end
end
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w');
assert(fid>0, 'Could not write config: %s', file);
fprintf(fid,'%s\n',txt);
fclose(fid);
end

function cmap = redblue()
n = 256;
r = [(0:n/2-1)'/(n/2); ones(n/2,1)];
g = [(0:n/2-1)'/(n/2); flipud((0:n/2-1)'/(n/2))];
b = [ones(n/2,1); flipud((0:n/2-1)'/(n/2))];
cmap = [r g b];
end
