%% analyze_test_46_test33_style_logk_correction.m
% Test 46: Test33-style mixedness-gated residual log-k correction.
%
% Goal:
%   Apply the part of Test33 that behaved well: a smooth residual correction
%   in log-k space, blended by predicted mixedness. This test deliberately
%   avoids pseudo-regions, hard regional medians, nearest-neighbor
%   interpolation, TheoryQDiscrete inputs, and candidate maps that can create
%   geometric artifacts.
%
% Frozen base:
%   q_spectrum_plus_composition.
%
% Operational predictors:
%   base q/SWS, q_spectrum_only reference, predicted composition/mixedness,
%   spectral disagreement proxies, M/frequency/regime/geometry metadata, and
%   local base-map variation. True SWS, true purity, true side, and
%   distance-to-interface are labels/evaluation groups only.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST46_MODE          = quick | medium | full
%   ADAPTIVE_REQ_TEST46_SOURCE        = test39 | kwave | both
%   ADAPTIVE_REQ_TEST46_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST46_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST46_USE_PARALLEL  = true | false
%   ADAPTIVE_REQ_TEST46_MAX_MAP_CONDITIONS = numeric

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.12);

CFG = default_config();
CFG = setup_parallel_if_requested(CFG);
OUT = make_output_dirs(root_dir, CFG);
SRC = locate_sources(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test46_configuration.json'));

fprintf('\nTest 46: Test33-style mixedness-gated log-k correction\n');
fprintf('Mode: %s | source: %s | validate only: %d\n', CFG.Mode, CFG.Source, CFG.ValidateOnly);
fprintf('Frozen base model: %s. No pseudo-regions, no Theory inputs, no interpolation correction.\n', CFG.BaseModel);

assert(exist(SRC.test38_csv,'file') == 2, 'Missing Test38 table: %s', SRC.test38_csv);
if CFG.UseTest39
    assert(exist(SRC.test39_csv,'file') == 2, 'Missing Test39 table: %s', SRC.test39_csv);
end
if CFG.UseKwave
    assert(exist(SRC.test40_csv,'file') == 2, 'Missing Test40/k-Wave table: %s', SRC.test40_csv);
end

T38_all = add_variant_reference_features(standardize_synthetic_table(readtable(SRC.test38_csv), "synthetic_train_source"));
T38_base = T38_all(T38_all.model_name == CFG.BaseModel,:);
assert(any(T38_base.is_train_row), 'No Test38 training rows for %s.', CFG.BaseModel);
assert(any(T38_base.is_heldout_row), 'No Test38 held-out rows for %s.', CFG.BaseModel);

eval_parts = {};
eval_parts{end+1,1} = set_source_domain(T38_base(T38_base.is_heldout_row,:), "synthetic_heldout"); %#ok<SAGROW>
if CFG.UseTest39
    T39 = add_variant_reference_features(standardize_synthetic_table(readtable(SRC.test39_csv), "test39_ood"));
    T39 = T39(T39.model_name == CFG.BaseModel,:);
    eval_parts{end+1,1} = set_source_domain(T39, "test39_ood"); %#ok<SAGROW>
end
if CFG.UseKwave
    T40 = add_variant_reference_features(standardize_kwave_table(readtable(SRC.test40_csv)));
    T40 = T40(T40.model_name == CFG.BaseModel,:);
    if CFG.QuickMode, T40 = quick_subset_by_condition(T40, CFG); end
    eval_parts{end+1,1} = T40; %#ok<SAGROW>
end
T_eval = vertcat_compatible(eval_parts{:});

T_train = T38_base(T38_base.is_train_row,:);
if CFG.QuickMode
    T_train = stratified_row_sample(T_train, CFG.QuickTrainRows, CFG.RandomSeed);
end
if height(T_train) > CFG.MaxTrainRows
    T_train = stratified_row_sample(T_train, CFG.MaxTrainRows, CFG.RandomSeed);
end
if CFG.ValidateOnly
    T_train = stratified_row_sample(T_train, min(2500,height(T_train)), CFG.RandomSeed);
    T_eval = stratified_row_sample(T_eval, min(2500,height(T_eval)), CFG.RandomSeed+1);
end
fprintf('Train rows: %d | evaluation rows: %d.\n', height(T_train), height(T_eval));

rng(CFG.RandomSeed, 'twister');
[MODEL, TRAIN_INFO] = train_logk_stack(T_train, CFG);
if CFG.ValidateOnly, validate_predictor_policy(MODEL.feature_names); end

T_results = apply_test46_strategies(T_eval, MODEL, CFG);
assert(all(isfinite(T_results.sws_pred_final)), 'Non-finite SWS predictions found.');

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
T_gate = gate_metrics(T_results);
T_cal = reliability_calibration(T_results);
T_reference = reference_comparison_table(root_dir, T_by_source);

writetable(T_results, fullfile(OUT.table_dir, 'test46_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test46_strategy_summary_overall.csv'));
writetable(T_by_source, fullfile(OUT.table_dir, 'test46_strategy_summary_by_source.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test46_strategy_summary_by_M.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'test46_strategy_summary_by_frequency.csv'));
writetable(T_by_frequency_M, fullfile(OUT.table_dir, 'test46_strategy_summary_by_frequency_M.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test46_strategy_summary_by_regime.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test46_strategy_summary_by_region.csv'));
writetable(T_by_predmix, fullfile(OUT.table_dir, 'test46_strategy_summary_by_predicted_mixedness.csv'));
writetable(T_by_truepurity, fullfile(OUT.table_dir, 'test46_strategy_summary_by_true_purity.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test46_strategy_summary_by_distance.csv'));
writetable(T_gate, fullfile(OUT.table_dir, 'test46_gate_metrics.csv'));
writetable(T_cal, fullfile(OUT.table_dir, 'test46_reliability_calibration.csv'));
writetable(T_reference, fullfile(OUT.table_dir, 'test46_reference_comparison_test33_test34.csv'));

save(fullfile(OUT.data_dir, 'test46_compact_results.mat'), ...
    'CFG','SRC','MODEL','TRAIN_INFO','T_results','T_overall','T_by_source', ...
    'T_by_M','T_by_frequency','T_by_frequency_M','T_by_regime','T_by_region', ...
    'T_by_predmix','T_by_truepurity','T_by_distance','T_gate','T_cal','T_reference','-v7.3');

safe_plot(@() plot_strategy_ranking(T_by_source, OUT, "MAPE_pct", ...
    'test46_strategy_ranking_mape.png'), 'strategy ranking MAPE');
safe_plot(@() plot_strategy_ranking(T_by_source, OUT, "high_error20_pct", ...
    'test46_strategy_ranking_high20.png'), 'strategy ranking high20');
safe_plot(@() plot_reference_comparison(T_reference, OUT), 'Test33/Test34 comparison');
safe_plot(@() plot_error_by_region(T_by_region, OUT), 'error by region');
safe_plot(@() plot_grouped_lines(T_by_distance, "distance_bin_pretty", CFG.PlotStrategies, ...
    "MAPE_pct", 'Distance to interface (mm)', 'MAPE (%)', ...
    'Error versus distance to interface', fullfile(OUT.figure_dir, 'test46_error_vs_distance.png')), 'distance');
safe_plot(@() plot_grouped_lines(T_by_predmix, "predicted_mixedness_bin", CFG.PlotStrategies, ...
    "MAPE_pct", 'Predicted mixedness probability bin', 'MAPE (%)', ...
    'Error versus predicted mixedness', fullfile(OUT.figure_dir, 'test46_error_vs_predicted_mixedness.png')), 'predicted mixedness');
safe_plot(@() plot_reliability_calibration(T_cal, OUT), 'reliability calibration');
safe_plot(@() plot_frequency_dependence(T_by_frequency_M, OUT, 2), 'frequency dependence M2');
safe_plot(@() plot_frequency_dependence(T_by_frequency_M, OUT, 3), 'frequency dependence M3');
if CFG.SaveAllMaps || CFG.ValidateOnly
    safe_plot(@() plot_condition_maps(T_results, CFG, OUT), 'condition maps');
end

print_console_summary(T_by_source, T_by_region, T_gate);
fprintf('\nTables: %s\nFigures: %s\nData: %s\nTest 46 complete.\n', ...
    OUT.table_dir, OUT.figure_dir, OUT.data_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST46_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","medium","full"]), ...
    'ADAPTIVE_REQ_TEST46_MODE must be quick, medium, or full.');
source = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST46_SOURCE'))));
if source == "", source = "both"; end
assert(ismember(source, ["test39","kwave","both"]), ...
    'ADAPTIVE_REQ_TEST46_SOURCE must be test39, kwave, or both.');

CFG = struct();
CFG.Mode = mode;
CFG.Source = source;
CFG.QuickMode = mode == "quick";
CFG.MediumMode = mode == "medium";
CFG.FullMode = mode == "full";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST46_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST46_SAVE_ALL_MAPS', true);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST46_USE_PARALLEL', false);
CFG.UseTest39 = ismember(source, ["test39","both"]);
CFG.UseKwave = ismember(source, ["kwave","both"]);
CFG.BaseModel = "q_spectrum_plus_composition";
CFG.RandomSeed = 46001;
CFG.NumTrees = env_number('ADAPTIVE_REQ_TEST46_NUM_TREES', ternary(mode=="quick", 80, 160));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST46_MIN_LEAF_SIZE', 70);
CFG.Shrinkage = env_number('ADAPTIVE_REQ_TEST46_SHRINKAGE', 0.05);
CFG.QuickTrainRows = env_number('ADAPTIVE_REQ_TEST46_QUICK_TRAIN_ROWS', 25000);
CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST46_MAX_TRAIN_ROWS', 120000);
CFG.QuickKwaveConditions = env_number('ADAPTIVE_REQ_TEST46_QUICK_KWAVE_CONDITIONS', 8);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST46_MAX_MAP_CONDITIONS', ternary(mode=="quick", 8, 48));
CFG.PhysicalRangeSWS = [0.5 10.0];
CFG.GateFloor = env_number('ADAPTIVE_REQ_TEST46_GATE_FLOOR', 0.0);
CFG.ConservativeGateScale = env_number('ADAPTIVE_REQ_TEST46_CONSERVATIVE_GATE_SCALE', 0.65);
CFG.MainStrategies = ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "residual_logk_all","mixedness_logk_corrected", ...
    "conservative_mixedness_logk","posterior_reliability_only", ...
    "oracle_delta_logk"];
CFG.PlotStrategies = ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "mixedness_logk_corrected","conservative_mixedness_logk"];
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
    CFG.BoostOptions = statset('UseParallel', false);
catch ME
    warning('Could not start/use parallel pool (%s). Continuing serial.', ME.message);
    CFG.UseParallel = false;
    CFG.BagOptions = statset('UseParallel', false);
    CFG.BoostOptions = statset('UseParallel', false);
end
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_46_test33_style_logk_correction');
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

%% Data

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
    T.p_strong_mixed = max(0, min(1, 2*T.p_mixed - 1));
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

function T = set_source_domain(T, src)
T.source_domain = repmat(src, height(T), 1);
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

function T = convert_strings(T)
for v = ["dataset","source_domain","condition_key","case_id","case_family", ...
        "case_name","field_regime","field_regime_ood","purity_bin","distance_bin", ...
        "material_region","region_zone","material_zone","model_name","model_family", ...
        "strategy_name","strategy_name_original","geometry","geometry_type", ...
        "material_side","analysis_region","roi_name","region_label","bundle_id"]
    if ismember(v, string(T.Properties.VariableNames)), T.(v) = string(T.(v)); end
end
end

function T = normalize_logical_columns(T)
for v = ["is_train_row","is_heldout_row","is_mixed","is_strong_mixed", ...
        "high_error10","high_error20"]
    if ismember(v, string(T.Properties.VariableNames)), T.(v) = logical(T.(v)); end
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

function [MODEL, info] = train_logk_stack(T, CFG)
[F, feature_names] = build_operational_features(T);
validate_predictor_policy(feature_names);
base_sws = clip_sws(T.sws_base, CFG);
true_k = 2*pi*T.f0 ./ T.true_SWS;
base_k = 2*pi*T.f0 ./ base_sws;
target_delta_logk = log(true_k) - log(base_k);
valid = all(isfinite(F),2) & isfinite(target_delta_logk) & T.true_SWS > 0;
F = F(valid,:); Tv = T(valid,:); target_delta_logk = target_delta_logk(valid);
tree_reg = templateTree('MinLeafSize', CFG.MinLeafSize, 'MaxNumSplits', 512);
tree_cls = templateTree('MinLeafSize', CFG.MinLeafSize, 'MaxNumSplits', 512);

fprintf('Training Test33-style residual log-k model on %d rows...\n', height(Tv));
delta_model = fitrensemble(F, target_delta_logk, ...
    'Method','LSBoost','Learners',tree_reg,'NumLearningCycles',CFG.NumTrees, ...
    'LearnRate',CFG.Shrinkage,'Options',CFG.BoostOptions);
delta_train = predict(delta_model, F);
gate = mixedness_gate(Tv, CFG, false);
base_err = abs(100*(base_sws(valid)-Tv.true_SWS)./Tv.true_SWS);
sws_corr = apply_delta_logk(Tv, base_sws(valid), delta_train, gate, CFG);
corr_err = abs(100*(sws_corr-Tv.true_SWS)./Tv.true_SWS);
high20_corr = corr_err > 20;

Freliab = [F delta_train gate abs(delta_train) abs(sws_corr-base_sws(valid))./max(base_sws(valid),0.25)];
reliability_model = fitcensemble(Freliab, high20_corr, ...
    'Method','Bag','Learners',tree_cls,'NumLearningCycles',max(60,round(CFG.NumTrees/2)), ...
    'ClassNames',[false true],'Options',CFG.BagOptions);

MODEL = struct();
MODEL.feature_names = feature_names;
MODEL.delta_model = delta_model;
MODEL.reliability_model = reliability_model;
MODEL.reliability_feature_names = [feature_names "delta_logk_pred" "mixedness_gate" ...
    "abs_delta_logk_pred" "relative_correction_magnitude"];
info = table(height(Tv), mean(base_err,'omitnan'), mean(corr_err,'omitnan'), ...
    mean(abs(delta_train),'omitnan'), mean(gate,'omitnan'), ...
    'VariableNames', {'N_train','train_base_MAPE_pct','train_mixedness_logk_MAPE_pct', ...
    'train_mean_abs_delta_logk','train_mean_gate'});
disp(info);
end

function T_all = apply_test46_strategies(T, MODEL, CFG)
[F, feature_names] = build_operational_features(T);
assert(isequal(string(feature_names(:)), string(MODEL.feature_names(:))), 'Feature mismatch.');
base_sws = clip_sws(T.sws_base, CFG);
qonly_sws = clip_sws(fill_default(T.q_spectrum_only_sws, base_sws), CFG);
delta_logk = predict(MODEL.delta_model, F);
gate = mixedness_gate(T, CFG, false);
gate_cons = mixedness_gate(T, CFG, true);
sws_all = apply_delta_logk(T, base_sws, delta_logk, ones(height(T),1), CFG);
sws_mixed = apply_delta_logk(T, base_sws, delta_logk, gate, CFG);
sws_cons = apply_delta_logk(T, base_sws, delta_logk, gate_cons, CFG);
true_delta = log(2*pi*T.f0./T.true_SWS) - log(2*pi*T.f0./base_sws);
sws_oracle = apply_delta_logk(T, base_sws, true_delta, gate, CFG);

Freliab = [F delta_logk gate abs(delta_logk) abs(sws_mixed-base_sws)./max(base_sws,0.25)];
[~, rel_score] = predict(MODEL.reliability_model, Freliab);
p_high20 = positive_score(MODEL.reliability_model, rel_score);

parts = {};
parts{end+1,1} = make_strategy_rows(T, "base_q_spectrum_plus_composition", base_sws, ...
    zeros(height(T),1), false(height(T),1), delta_logk, gate, p_high20, "operational_base"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "q_spectrum_only_reference", qonly_sws, ...
    qonly_sws-base_sws, true(height(T),1), delta_logk, gate, p_high20, "reference"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "residual_logk_all", sws_all, ...
    sws_all-base_sws, true(height(T),1), delta_logk, ones(height(T),1), p_high20, "diagnostic_aggressive"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "mixedness_logk_corrected", sws_mixed, ...
    sws_mixed-base_sws, gate>0.02, delta_logk, gate, p_high20, "operational_test33_style"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "conservative_mixedness_logk", sws_cons, ...
    sws_cons-base_sws, gate_cons>0.02, delta_logk, gate_cons, p_high20, "operational_conservative"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "posterior_reliability_only", base_sws, ...
    zeros(height(T),1), false(height(T),1), delta_logk, gate, p_high20, "reliability_only"); %#ok<AGROW>
parts{end+1,1} = make_strategy_rows(T, "oracle_delta_logk", sws_oracle, ...
    sws_oracle-base_sws, gate>0.02, true_delta, gate, p_high20, "diagnostic_oracle"); %#ok<AGROW>
T_all = vertcat_compatible(parts{:});
T_all.is_harmed = T_all.was_corrected & T_all.sws_abs_error_pct > T_all.base_abs_error_pct + 1;
T_all.is_improved = T_all.was_corrected & T_all.sws_abs_error_pct + 1 < T_all.base_abs_error_pct;
T_all.is_unchanged_like = T_all.was_corrected & ~T_all.is_harmed & ~T_all.is_improved;
T_all.net_abs_error_gain_pct_points = T_all.base_abs_error_pct - T_all.sws_abs_error_pct;
end

function R = make_strategy_rows(T, strategy_name, sws, correction_mps, was_corrected, ...
    delta_logk, gate, p_high20, strategy_type)
R = T;
R.strategy_name = repmat(strategy_name, height(T), 1);
R.strategy_type = repmat(strategy_type, height(T), 1);
R.sws_pred_final = sws;
R.sws_pred = sws;
R.correction_mps = correction_mps;
R.correction_pct = 100*correction_mps ./ max(T.sws_base, 0.25);
R.was_corrected = was_corrected;
R.was_kept = ~was_corrected;
R.delta_logk_pred = delta_logk;
R.mixedness_gate = gate;
R.final_high_error_probability = p_high20;
R.final_reliability = 1 - p_high20;
R.base_signed_error_pct = 100*(T.sws_base - T.true_SWS)./T.true_SWS;
R.base_abs_error_pct = abs(R.base_signed_error_pct);
R.sws_signed_error_pct = 100*(sws - T.true_SWS)./T.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
end

function sws = apply_delta_logk(T, base_sws, delta_logk, gate, CFG)
base_k = 2*pi*T.f0 ./ base_sws;
logk = log(base_k) + gate .* delta_logk;
sws = clip_sws(2*pi*T.f0 ./ exp(logk), CFG);
end

function gate = mixedness_gate(T, CFG, conservative)
p = fill_default(T.p_mixed, 0);
ps = fill_default(T.p_strong_mixed, max(0,2*p-1));
gate = max(p, 0.35*ps);
gate = max(gate, CFG.GateFloor);
if conservative
    gate = CFG.ConservativeGateScale * gate;
    gate(gate < 0.20) = 0;
end
if ismember("case_family", string(T.Properties.VariableNames))
    gate(T.case_family == "homogeneous") = 0;
end
gate = max(0, min(1, gate));
end

function [X, names] = build_operational_features(T)
names = ["M_norm","f0_norm","regime_code","geometry_code", ...
    "predicted_patch_purity","predicted_impurity","p_mixed","p_strong_mixed", ...
    "log_sws_base","q_base","q_spectrum_only_log_sws","q_spectrum_only_q", ...
    "variant_sws_disagreement","variant_q_disagreement","base_sws_gradient_proxy", ...
    "mixedness_need_proxy"];
base = max(T.sws_base, 0.1);
qonly_sws = fill_default(T.q_spectrum_only_sws, base);
qonly_q = fill_default(T.q_spectrum_only_q, T.q_base);
purity = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-purity);
ps = fill_default(T.p_strong_mixed, 0);
X = [normalize_like(T.M, 3, 1), normalize_like(T.f0, 500, 150), ...
    regime_code(T), geometry_code(T), purity, 1-purity, pm, ps, ...
    log(base), fill_default(T.q_base, 0.5), log(max(qonly_sws,0.1)), qonly_q, ...
    abs(base-qonly_sws)./max(base,0.25), fill_default(T.q_base,0.5)-qonly_q, ...
    neighborhood_variation_proxy(T), heuristic_mixedness_need(T)];
X = double(X);
X(~isfinite(X)) = 0;
validate_predictor_policy(names);
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

function T = gate_metrics(P)
X = P(P.was_corrected | ismember(P.strategy_name, ["q_spectrum_only_reference", ...
    "residual_logk_all","mixedness_logk_corrected","conservative_mixedness_logk", ...
    "oracle_delta_logk"]),:);
[G, groups] = findgroups(X(:, {'source_domain','strategy_name'}));
T = [groups table(splitapply(@numel, X.was_corrected, G), ...
    splitapply(@(x) 100*mean(x,'omitnan'), X.was_corrected, G), ...
    splitapply(@rate_among_corrected, X.is_improved, X.was_corrected, G), ...
    splitapply(@rate_among_corrected, X.is_harmed, X.was_corrected, G), ...
    splitapply(@(x) mean(x,'omitnan'), X.net_abs_error_gain_pct_points, G), ...
    'VariableNames', {'N','corrected_pct','improved_among_corrected_pct', ...
    'harmed_among_corrected_pct','mean_net_gain_points'})];
end

function T = reliability_calibration(P)
X = P(P.strategy_name == "mixedness_logk_corrected" | P.strategy_name == "posterior_reliability_only",:);
edges = 0:0.1:1;
bin = discretize(X.final_high_error_probability, edges, 'IncludedEdge','right');
labels = strings(size(bin));
for i = 1:numel(bin)
    if isfinite(bin(i)), labels(i) = sprintf('%.1f-%.1f', edges(bin(i)), edges(bin(i)+1));
    else, labels(i) = "unknown"; end
end
X.reliability_bin = labels;
[G, groups] = findgroups(X(:, {'source_domain','strategy_name','reliability_bin'}));
T = [groups table(splitapply(@numel, X.high_error20, G), ...
    splitapply(@(x) mean(x,'omitnan'), X.final_high_error_probability, G), ...
    splitapply(@(x) 100*mean(x,'omitnan'), X.high_error20, G), ...
    'VariableNames', {'N','mean_predicted_high_error_probability', ...
    'observed_high_error20_pct'})];
end

function T = reference_comparison_table(root_dir, T46_source)
T = table();
keep46 = ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "mixedness_logk_corrected","conservative_mixedness_logk","oracle_delta_logk"];
for i = 1:height(T46_source)
    if ~ismember(T46_source.strategy_name(i), keep46), continue; end
    R = table("Test46 " + pretty_source_name(T46_source.source_domain(i)), ...
        T46_source.strategy_name(i), T46_source.MAPE_pct(i), T46_source.high_error20_pct(i), ...
        "test46", 'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
p33 = fullfile(root_dir, 'outputs', 'test_33_mixedness_aware_q_correction', ...
    'tables', 'test33_strategy_summary_overall_test.csv');
if exist(p33,'file') == 2
    A = readtable(p33); A.strategy_name = string(A.strategy_name);
    keep = ["local_baseline","mixedness_logk_corrected","mixedness_q_candidate_selector"];
    A = A(ismember(A.strategy_name, keep),:);
    R = table(repmat("Test33 synthetic reference",height(A),1), A.strategy_name, ...
        A.MAPE, A.high_error_20_pct, repmat("test33",height(A),1), ...
        'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
p34 = fullfile(root_dir, 'outputs', 'test_34_kwave_mixedness_transfer', ...
    'tables', 'test34_kwave_strategy_summary_overall.csv');
if exist(p34,'file') == 2
    A = readtable(p34); A.strategy_name = string(A.strategy_name);
    keep = ["local_baseline","mixedness_logk_corrected","mixedness_q_candidate_selector"];
    A = A(ismember(A.strategy_name, keep),:);
    R = table(repmat("Test34 k-Wave reference",height(A),1), A.strategy_name, ...
        A.MAPE, A.high_error_20_pct, repmat("test34",height(A),1), ...
        'VariableNames', {'benchmark','strategy_name','MAPE_pct','high_error20_pct','source_family'});
    T = concat_tables(T, R);
end
end

%% Plots

function plot_condition_maps(P, CFG, OUT)
conds = unique(P(:, {'source_domain','case_id','field_regime_ood','f0','M','condition_key'}), 'rows', 'stable');
if isfinite(CFG.MaxMapConditions) && height(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions,:);
end
fprintf('Saving %d Test46 Test33-style map panels under %s.\n', height(conds), OUT.map_dir);
for ci = 1:height(conds)
    idx = P.source_domain == conds.source_domain(ci) & P.condition_key == conds.condition_key(ci) & P.M == conds.M(ci);
    Xc = P(idx,:);
    B = Xc(Xc.strategy_name == "base_q_spectrum_plus_composition",:);
    Q = Xc(Xc.strategy_name == "q_spectrum_only_reference",:);
    A = Xc(Xc.strategy_name == "residual_logk_all",:);
    M = Xc(Xc.strategy_name == "mixedness_logk_corrected",:);
    C = Xc(Xc.strategy_name == "conservative_mixedness_logk",:);
    if isempty(B) || isempty(M), continue; end
    maps = {map_from_rows(B,B.true_SWS), map_from_rows(Q,Q.sws_pred_final), ...
        map_from_rows(B,B.sws_base), map_from_rows(A,A.sws_pred_final), ...
        map_from_rows(M,M.sws_pred_final), map_from_rows(C,C.sws_pred_final), ...
        map_from_rows(B,B.base_abs_error_pct), map_from_rows(M,M.sws_abs_error_pct), ...
        map_from_rows(C,C.sws_abs_error_pct), map_from_rows(M,M.p_mixed), ...
        map_from_rows(M,M.predicted_patch_purity), map_from_rows(M,M.mixedness_gate), ...
        map_from_rows(M,M.delta_logk_pred), map_from_rows(M,M.final_reliability)};
    titles = ["True SWS","q-spectrum only","Base q spectrum + composition", ...
        "Residual log-k all","Mixedness log-k","Conservative log-k", ...
        "Base error","Mixedness log-k error","Conservative error", ...
        "P(mixed)","Predicted purity","Mixedness gate","Delta log-k","Posterior reliability"];
    labels = ["SWS (m/s)","SWS (m/s)","SWS (m/s)","SWS (m/s)","SWS (m/s)","SWS (m/s)", ...
        "Abs error (%)","Abs error (%)","Abs error (%)","Probability","Probability","Gate weight", ...
        "Delta log-k","Reliability"];
    fig = figure('Color','w','Units','centimeters','Position',[1 1 36 31]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for k = 1:numel(maps)
        ax = nexttile(tl);
        imagesc(ax, maps{k}); axis(ax,'image'); axis(ax,'off');
        cb = colorbar(ax); ylabel(cb, labels(k));
        title(ax, titles(k), 'FontWeight','bold');
    end
    title(tl, sprintf('%s__f%d__%s__M%d | %s', conds.case_id(ci), conds.f0(ci), ...
        conds.field_regime_ood(ci), conds.M(ci), pretty_source_name(conds.source_domain(ci))), ...
        'Interpreter','none','FontWeight','bold');
    outdir = fullfile(OUT.map_dir, sanitize(conds.source_domain(ci)), sanitize(conds.case_id(ci)), sanitize(conds.field_regime_ood(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    export_fig(fig, fullfile(outdir, sprintf('test46__%s__f%d__%s__M%d.png', ...
        sanitize(conds.case_id(ci)), conds.f0(ci), sanitize(conds.field_regime_ood(ci)), conds.M(ci))));
end
end

function Z = map_from_rows(T, vals)
[izs,~,iz] = unique(T.map_iz);
[ixs,~,ix] = unique(T.map_ix);
Z = nan(numel(izs), numel(ixs));
for i = 1:height(T), Z(iz(i),ix(i)) = vals(i); end
end

function plot_strategy_ranking(T, OUT, metric, filename)
T = T(ismember(T.strategy_name, ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "residual_logk_all","mixedness_logk_corrected","conservative_mixedness_logk","oracle_delta_logk"]),:);
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
        X = T(contains(T.benchmark, "synthetic", 'IgnoreCase', true) | contains(T.benchmark, "OOD", 'IgnoreCase', true),:);
        ttl = 'Synthetic / OOD comparison';
    else
        X = T(contains(T.benchmark, "k-Wave", 'IgnoreCase', true) | contains(T.benchmark, "kwave", 'IgnoreCase', true),:);
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
export_fig(fig, fullfile(OUT.figure_dir, 'test46_reference_comparison_test33_test34.png'));
end

function plot_error_by_region(T, OUT)
keep = ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "mixedness_logk_corrected","conservative_mixedness_logk","oracle_delta_logk"];
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
export_fig(fig, fullfile(OUT.figure_dir, 'test46_error_by_region.png'));
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
        plot(ax, 1:numel(cats), y, '-o', 'DisplayName', pretty_strategy_name(s));
    end
    set(ax,'XTick',1:numel(cats),'XTickLabel',arrayfun(@pretty_region_name, cats), ...
        'TickLabelInterpreter','none');
    xtickangle(ax,30); ylabel(ax,ylabel_txt); xlabel(ax,xlabel_txt);
    title(ax, pretty_source_name(src), 'FontWeight','normal'); grid(ax,'on');
    legend(ax,'Location','bestoutside');
end
title(tl, ttl, 'FontWeight','bold');
export_fig(fig, file);
end

function plot_frequency_dependence(T, OUT, Mval)
X = T(T.M == Mval & ismember(T.strategy_name, ["base_q_spectrum_plus_composition", ...
    "q_spectrum_only_reference","mixedness_logk_corrected","conservative_mixedness_logk"]),:);
if isempty(X), return; end
plot_grouped_lines(X, "f0", ["base_q_spectrum_plus_composition","q_spectrum_only_reference", ...
    "mixedness_logk_corrected","conservative_mixedness_logk"], "MAPE_pct", ...
    'Frequency (Hz)', 'MAPE (%)', sprintf('Frequency dependence of SWS error, M=%d', Mval), ...
    fullfile(OUT.figure_dir, sprintf('test46_frequency_dependence_M%d.png', Mval)));
end

function plot_reliability_calibration(T, OUT)
if isempty(T), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 18 14]);
hold on;
for s = unique(T.source_domain, 'stable')'
    X = T(T.source_domain == s & T.strategy_name == "mixedness_logk_corrected",:);
    if isempty(X), continue; end
    [~,ord] = sort(X.mean_predicted_high_error_probability);
    plot(X.mean_predicted_high_error_probability(ord), X.observed_high_error20_pct(ord), '-o', ...
        'DisplayName', pretty_source_name(s));
end
plot([0 1],[0 100],'k--','DisplayName','Ideal calibration');
xlabel('Predicted high-error probability'); ylabel('Observed high-error rate (%)');
title('Reliability calibration', 'FontWeight','normal');
grid on; legend('Location','bestoutside');
export_fig(fig, fullfile(OUT.figure_dir, 'test46_reliability_calibration.png'));
end

%% Labels and utilities

function s = pretty_strategy_name(name)
name = string(name);
switch name
    case "base_q_spectrum_plus_composition", s = "Base: q spectrum + composition";
    case "q_spectrum_only_reference", s = "q-spectrum only";
    case "residual_logk_all", s = "Residual log-k applied everywhere";
    case "mixedness_logk_corrected", s = "Mixedness-gated log-k";
    case "conservative_mixedness_logk", s = "Conservative mixedness-gated log-k";
    case "posterior_reliability_only", s = "Reliability only";
    case "oracle_delta_logk", s = "Oracle delta log-k";
    case "local_baseline", s = "Test33 local baseline";
    case "mixedness_q_candidate_selector", s = "Test33 mixedness q selector";
    otherwise, s = strrep(name, "_", " ");
end
end

function s = pretty_source_name(name)
name = string(name);
switch name
    case "synthetic_train_source", s = "Synthetic train source";
    case "synthetic_heldout", s = "Synthetic held-out";
    case "test39_ood", s = "Synthetic OOD";
    case "kwave_transfer", s = "k-Wave transfer";
    otherwise, s = strrep(name, "_", " ");
end
end

function s = pretty_region_name(name)
name = string(name);
switch name
    case "homogeneous", s = "Homogeneous";
    case "pure", s = "Pure";
    case "near_pure", s = "Near-pure";
    case "moderately_mixed", s = "Moderately mixed";
    case "strongly_mixed", s = "Strongly mixed";
    case "interface", s = "Interface/mixed";
    case "0-1", s = "0-1 mm";
    case "1-2", s = "1-2 mm";
    case "2-4", s = "2-4 mm";
    case "4-8", s = "4-8 mm";
    case ">8", s = ">8 mm";
    otherwise, s = strrep(name, "_", " ");
end
end

function s = pretty_metric_name(name)
if name == "MAPE_pct", s = "MAPE (%)";
elseif name == "high_error20_pct", s = "Pixels with error >20% (%)";
else, s = strrep(string(name), "_", " ");
end
end

function r = evaluation_region(T)
r = repmat("pure", height(T), 1);
r(T.case_family == "homogeneous") = "homogeneous";
r(T.case_family ~= "homogeneous" & T.patch_purity < 0.70) = "strongly_mixed";
r(T.case_family ~= "homogeneous" & T.patch_purity >= 0.70 & T.patch_purity < 0.90) = "moderately_mixed";
r(T.case_family ~= "homogeneous" & T.patch_purity >= 0.90 & T.patch_purity < 0.99) = "near_pure";
r(T.case_family ~= "homogeneous" & T.patch_purity >= 0.99) = "pure";
r(T.case_family ~= "homogeneous" & T.distance_to_interface_mm <= 2) = "interface";
end

function b = true_purity_bin(p, family)
b = repmat("pure", numel(p), 1);
b(family == "homogeneous") = "homogeneous";
b(family ~= "homogeneous" & p < 0.70) = "strongly_mixed";
b(family ~= "homogeneous" & p >= 0.70 & p < 0.90) = "moderately_mixed";
b(family ~= "homogeneous" & p >= 0.90 & p < 0.99) = "near_pure";
b(family ~= "homogeneous" & p >= 0.99) = "pure";
end

function b = probability_bin(p)
edges = [0 0.2 0.4 0.6 0.8 1.000001];
labels = ["0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"];
idx = discretize(p, edges);
b = repmat("unknown", numel(p), 1);
ok = isfinite(idx);
b(ok) = labels(idx(ok));
end

function b = distance_bin_pretty(d, family)
b = repmat(">8", numel(d), 1);
b(family == "homogeneous") = "homogeneous";
b(family ~= "homogeneous" & d <= 1) = "0-1";
b(family ~= "homogeneous" & d > 1 & d <= 2) = "1-2";
b(family ~= "homogeneous" & d > 2 & d <= 4) = "2-4";
b(family ~= "homogeneous" & d > 4 & d <= 8) = "4-8";
end

function z = pretty_distance_zone(d)
z = repmat("far", numel(d),1);
z(d <= 2) = "interface";
z(d > 2 & d <= 8) = "near";
end

function code = regime_code(T)
[~,~,code] = unique(T.field_regime_ood);
code = normalize_like(code, median(code,'omitnan'), max(1,std(double(code),'omitnan')));
end

function code = geometry_code(T)
[~,~,code] = unique(T.case_family);
code = normalize_like(code, median(code,'omitnan'), max(1,std(double(code),'omitnan')));
end

function y = normalize_like(x, mu, scale)
y = (double(x) - mu) ./ max(scale, eps);
end

function y = fill_default(y, fallback)
if isscalar(fallback), fallback = repmat(fallback, size(y)); end
y(~isfinite(y)) = fallback(~isfinite(y));
end

function y = heuristic_mixedness_need(T)
p = fill_default(T.p_mixed, 0);
purity = fill_default(T.predicted_patch_purity, 0.95);
dis = fill_default(T.variant_sws_disagreement, 0);
y = max(p, max(0, 1-purity));
y = max(y, min(1, 2*dis));
if ismember("case_family", string(T.Properties.VariableNames))
    y(T.case_family == "homogeneous") = 0;
end
y = max(0, min(1, y));
end

function v = neighborhood_variation_proxy(T)
v = zeros(height(T),1);
[G,~] = findgroups(T.condition_key);
for gi = unique(G)'
    idx = find(G == gi);
    S = T(idx,:);
    base = T.sws_base(idx);
    for ii = 1:numel(idx)
        close = abs(S.map_ix - S.map_ix(ii)) <= 1 & abs(S.map_iz - S.map_iz(ii)) <= 1;
        vals = base(close & isfinite(base));
        if numel(vals) >= 3
            v(idx(ii)) = mad(vals,1) / max(abs(median(vals,'omitnan')),0.25);
        end
    end
end
v(~isfinite(v)) = 0;
v = min(v, 2);
end

function y = clip_sws(y, CFG)
y = min(max(y, CFG.PhysicalRangeSWS(1)), CFG.PhysicalRangeSWS(2));
end

function y = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes=="true" | classes=="1", 1);
if isempty(idx), idx = size(score,2); end
y = score(:,idx);
end

function T = concat_tables(T, R)
if isempty(T), T = R; else, T = [T; R]; end %#ok<AGROW>
end

function T = vertcat_compatible(varargin)
T = varargin{1};
for k = 2:nargin
    U = varargin{k};
    allv = union(string(T.Properties.VariableNames), string(U.Properties.VariableNames), 'stable');
    T = add_missing_vars(T, allv);
    U = add_missing_vars(U, allv);
    T = [T(:,cellstr(allv)); U(:,cellstr(allv))]; %#ok<AGROW>
end
end

function T = add_missing_vars(T, vars)
for v = vars(:)'
    if ~ismember(v, string(T.Properties.VariableNames))
        T.(v) = nan(height(T),1);
    end
end
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Plot failed (%s): %s', label, ME.message);
end
end

function write_config_json(CFG, path)
fid = fopen(path,'w');
if fid < 0, return; end
fprintf(fid, '%s\n', jsonencode(CFG, PrettyPrint=true));
fclose(fid);
end

function tf = env_true(name, default)
v = lower(strtrim(string(getenv(name))));
if v == "", tf = default; else, tf = ismember(v, ["1","true","yes","on"]); end
end

function x = env_number(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; return; end
x = str2double(v);
if ~isfinite(x), x = default; end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = sanitize(x)
s = regexprep(string(x), '[^A-Za-z0-9_\\-]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
end

function validate_predictor_policy(names)
names = lower(string(names));
forbidden = ["true","error","abs_error","signed_error","true_patch_purity", ...
    "material_side","distance","interface","roi","oracle","theory"];
hit = names(contains(names, forbidden));
assert(isempty(hit), 'Oracle/evaluation feature leaked into Test46 predictors: %s', strjoin(hit, ', '));
end

function export_fig(fig, path)
if exist(fileparts(path),'dir') ~= 7, mkdir(fileparts(path)); end
exportgraphics(fig, path, 'Resolution', 180);
close(fig);
end

function print_console_summary(T_by_source, T_by_region, T_gate)
fprintf('\n================ Test 46 summary ================\n');
for src = unique(T_by_source.source_domain, 'stable')'
    X = T_by_source(T_by_source.source_domain == src & T_by_source.strategy_name ~= "oracle_delta_logk",:);
    [~,i] = min(X.MAPE_pct);
    fprintf('Best non-oracle on %s: %s (MAPE %.2f%%, HE20 %.2f%%).\n', ...
        pretty_source_name(src), pretty_strategy_name(X.strategy_name(i)), X.MAPE_pct(i), X.high_error20_pct(i));
end
base = T_by_region(T_by_region.strategy_name=="base_q_spectrum_plus_composition",:);
main = T_by_region(T_by_region.strategy_name=="mixedness_logk_corrected",:);
[ok,loc] = ismember(main(:, {'source_domain','evaluation_region'}), base(:, {'source_domain','evaluation_region'}));
for i = find(ok)'
    gain = base.MAPE_pct(loc(i)) - main.MAPE_pct(i);
    if ismember(main.evaluation_region(i), ["interface","strongly_mixed","pure","homogeneous"])
        fprintf('Mixedness log-k gain on %s / %s: %.2f MAPE points.\n', ...
            pretty_source_name(main.source_domain(i)), pretty_region_name(main.evaluation_region(i)), gain);
    end
end
G = T_gate(T_gate.strategy_name=="mixedness_logk_corrected",:);
for i = 1:height(G)
    fprintf('%s mixedness log-k corrected %.1f%% pixels; improved %.1f%%, harmed %.1f%% among corrected.\n', ...
        pretty_source_name(G.source_domain(i)), G.corrected_pct(i), ...
        G.improved_among_corrected_pct(i), G.harmed_among_corrected_pct(i));
end
fprintf('=================================================\n');
end
