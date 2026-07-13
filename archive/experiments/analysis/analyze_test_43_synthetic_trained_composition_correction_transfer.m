%% analyze_test_43_synthetic_trained_composition_correction_transfer.m
% Test 43: synthetic-trained composition-aware correction transfer.
%
% This test does NOT retrain the Test38 composition/mixedness detector or the
% base q models. It trains only residual correction layers on Test38 synthetic
% training rows, using the already predicted composition variables
% (`predicted_patch_purity`, `p_mixed`, `p_strong_mixed`) and operational
% q/theory/model disagreement features. The frozen correction layers are then
% evaluated on:
%   1) Test38 synthetic held-out rows;
%   2) Test39 OOD synthetic external rows;
%   3) Test40 k-Wave transfer rows.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST43_MODE               = quick | medium | full
%   ADAPTIVE_REQ_TEST43_TRAIN_SOURCE       = quick | medium | full
%   ADAPTIVE_REQ_TEST43_TEST39_SOURCE      = quick | medium | full
%   ADAPTIVE_REQ_TEST43_VALIDATE_ONLY      = true | false
%   ADAPTIVE_REQ_TEST43_SAVE_ALL_MAPS      = true | false
%   ADAPTIVE_REQ_TEST43_MAX_MAP_CONDITIONS = integer or Inf
%   ADAPTIVE_REQ_TEST43_USE_PARALLEL       = true | false

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
write_config_json(CFG, fullfile(OUT.root_dir, 'test43_configuration.json'));
SRC = locate_sources(root_dir, CFG);

fprintf('\nTest 43: synthetic-trained composition-aware correction transfer\n');
fprintf('Mode: %s | train source: %s | Test39 source: %s\n', ...
    CFG.Mode, CFG.TrainSource, CFG.Test39Source);
fprintf('Base models: %s\n', strjoin(CFG.BaseModels, ', '));
fprintf('Parallel training requested: %d\n', CFG.UseParallel);

assert(exist(SRC.test38_csv,'file') == 2, 'Missing Test38 table: %s', SRC.test38_csv);
assert(exist(SRC.test39_csv,'file') == 2, 'Missing Test39 table: %s', SRC.test39_csv);
assert(exist(SRC.test40_csv,'file') == 2, 'Missing Test40 k-Wave table: %s', SRC.test40_csv);

T38 = standardize_synthetic_table(readtable(SRC.test38_csv), "test38_synthetic");
T39 = standardize_synthetic_table(readtable(SRC.test39_csv), "test39_external");
T40 = standardize_kwave_table(readtable(SRC.test40_csv));
fprintf('Loaded rows: Test38=%d | Test39=%d | k-Wave/Test40=%d.\n', ...
    height(T38), height(T39), height(T40));

if CFG.ValidateOnly
    validate_inputs(T38, T39, T40, CFG);
    fprintf('Test 43 validation-only checks passed.\n');
    return;
end

rng(CFG.RandomSeed, 'twister');

%% Train correction layers on Test38 synthetic train rows

parts = {};
model_info = table();
for i = 1:numel(CFG.BaseModels)
    base = CFG.BaseModels(i);
    train_rows = T38(T38.model_name == base & T38.is_train_row,:);
    held_rows = T38(T38.model_name == base & T38.is_heldout_row,:);
    ext_rows = T39(T39.model_name == base,:);
    kw_strategy = "T38_" + base;
    kw_rows = T40(T40.model_name == base,:);

    assert(~isempty(train_rows), 'Missing synthetic train rows for %s.', base);
    assert(~isempty(held_rows), 'Missing synthetic held-out rows for %s.', base);
    assert(~isempty(ext_rows), 'Missing Test39 external rows for %s.', base);
    assert(~isempty(kw_rows), 'Missing Test40/k-Wave rows for %s (%s).', base, kw_strategy);

    fprintf('[%d/%d] Training correction + acceptance gate for %s (%d synthetic train rows)...\n', ...
        i, numel(CFG.BaseModels), base, height(train_rows));
    [MODEL, INFO] = train_correction_stack(train_rows, base, CFG);
    model_info = concat_tables(model_info, INFO);

    parts{end+1,1} = baseline_rows(held_rows, "synthetic_heldout", base); %#ok<AGROW>
    parts{end+1,1} = apply_correction_stack(held_rows, MODEL, "synthetic_heldout", base, CFG); %#ok<AGROW>

    parts{end+1,1} = baseline_rows(ext_rows, "test39_external", base); %#ok<AGROW>
    parts{end+1,1} = apply_correction_stack(ext_rows, MODEL, "test39_external", base, CFG); %#ok<AGROW>

    parts{end+1,1} = baseline_rows(kw_rows, "kwave_transfer", base); %#ok<AGROW>
    parts{end+1,1} = apply_correction_stack(kw_rows, MODEL, "kwave_transfer", base, CFG); %#ok<AGROW>
end

% Add diagnostic references.
for dataset_name = ["synthetic_heldout","test39_external"]
    source = ternary(dataset_name=="synthetic_heldout", T38, T39);
    R = source(source.model_name == "theory_discrete" & source.is_heldout_row,:);
    if dataset_name == "test39_external"
        R = source(source.model_name == "theory_discrete",:);
    end
    if ~isempty(R)
        parts{end+1,1} = baseline_rows(R, dataset_name, "theory_discrete"); %#ok<AGROW>
    end
end
for s = ["T34_local_baseline","T34_theory_baseline","T34_mixedness_logk_corrected"]
    R = T40(T40.strategy_name == s,:);
    if ~isempty(R)
        R.model_name = repmat(erase(s, "T34_"), height(R), 1);
        parts{end+1,1} = baseline_rows(R, "kwave_transfer", "T34_" + R.model_name(1)); %#ok<AGROW>
    end
end

T_all = vertcat_compatible(parts{:});

%% Summaries

T_overall = summarize_predictions(T_all, ["dataset_eval","model_family","strategy_name"]);
T_by_dataset = summarize_predictions(T_all, ["strategy_name","dataset_eval"]);
T_by_case = summarize_predictions(T_all, ["strategy_name","dataset_eval","case_id"]);
T_by_family = summarize_predictions(T_all, ["strategy_name","dataset_eval","case_family"]);
T_by_regime = summarize_predictions(T_all, ["strategy_name","dataset_eval","field_regime_ood"]);
T_by_frequency = summarize_predictions(T_all, ["strategy_name","dataset_eval","f0"]);
T_by_M = summarize_predictions(T_all, ["strategy_name","dataset_eval","M"]);
T_by_purity = summarize_predictions(T_all, ["strategy_name","dataset_eval","purity_bin"]);
T_by_material = summarize_predictions(T_all, ["strategy_name","dataset_eval","case_family","material_region"]);
T_by_zone = summarize_predictions(T_all, ["strategy_name","dataset_eval","case_family","region_zone"]);
T_by_material_zone = summarize_predictions(T_all, ["strategy_name","dataset_eval","case_id","material_zone"]);
T_gain = correction_gain_summary(T_all, CFG);
T_fail = top_failure_conditions(summarize_predictions(T_all, ...
    ["strategy_name","dataset_eval","condition_key","case_id","case_family","field_regime_ood","f0","M"]), 50);

writetable(T_all, fullfile(OUT.table_dir, 'test43_patch_level_results.csv'));
writetable(model_info, fullfile(OUT.table_dir, 'test43_correction_training_summary.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test43_summary_overall.csv'));
writetable(T_by_dataset, fullfile(OUT.table_dir, 'test43_summary_by_dataset.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test43_summary_by_case.csv'));
writetable(T_by_family, fullfile(OUT.table_dir, 'test43_summary_by_family.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test43_summary_by_regime.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'test43_summary_by_frequency.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test43_summary_by_M.csv'));
writetable(T_by_purity, fullfile(OUT.table_dir, 'test43_summary_by_purity_bin.csv'));
writetable(T_by_material, fullfile(OUT.table_dir, 'test43_summary_by_material_region.csv'));
writetable(T_by_zone, fullfile(OUT.table_dir, 'test43_summary_by_region_zone.csv'));
writetable(T_by_material_zone, fullfile(OUT.table_dir, 'test43_summary_by_material_zone.csv'));
writetable(T_gain, fullfile(OUT.table_dir, 'test43_correction_gain_summary.csv'));
writetable(T_fail, fullfile(OUT.table_dir, 'test43_worst_conditions.csv'));

save(fullfile(OUT.data_dir, 'test43_synthetic_trained_composition_correction_transfer.mat'), ...
    'T_all','model_info','T_overall','T_by_dataset','T_by_case','T_by_family', ...
    'T_by_regime','T_by_frequency','T_by_M','T_by_purity','T_by_material', ...
    'T_by_zone','T_by_material_zone','T_gain','T_fail','CFG','SRC','-v7.3');

%% Figures

safe_plot(@() plot_dataset_ranking(T_overall, OUT), 'dataset ranking');
safe_plot(@() plot_transfer_summary(T_by_dataset, OUT), 'transfer summary');
safe_plot(@() plot_region_diagnostics(T_by_material, T_by_zone, CFG, OUT), 'region diagnostics');
safe_plot(@() plot_heatmap_table(T_by_case, "case_id", OUT, 'test43_mape_by_case.png'), 'case heatmap');
safe_plot(@() plot_heatmap_table(T_by_regime, "field_regime_ood", OUT, 'test43_mape_by_regime.png'), 'regime heatmap');
safe_plot(@() plot_heatmap_table(T_by_purity, "purity_bin", OUT, 'test43_mape_by_purity_bin.png'), 'purity heatmap');
safe_plot(@() plot_correction_gain(T_gain, OUT), 'correction gain');
if CFG.SaveAllMaps
    safe_plot(@() plot_condition_maps(T_all, CFG, OUT), 'condition maps');
end

print_interpretation(T_overall, T_by_dataset, T_by_material, T_by_zone, T_gain, T_fail);
fprintf('\nTables: %s\nFigures: %s\nTest 43 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST43_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","medium","full"]), ...
    'ADAPTIVE_REQ_TEST43_MODE must be quick, medium, or full.');
train_src = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST43_TRAIN_SOURCE'))));
if train_src == "", train_src = mode; end
test39_src = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST43_TEST39_SOURCE'))));
if test39_src == "", test39_src = mode; end
CFG = struct();
CFG.Mode = mode;
CFG.TrainSource = train_src;
CFG.Test39Source = test39_src;
CFG.QuickMode = mode == "quick";
CFG.MediumMode = mode == "medium";
CFG.FullMode = mode == "full";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST43_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST43_SAVE_ALL_MAPS', true);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST43_MAX_MAP_CONDITIONS', ternary(CFG.QuickMode, 6, 24));
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST43_USE_PARALLEL', false);
CFG.RandomSeed = 43001;
CFG.NumLearningCycles = env_number('ADAPTIVE_REQ_TEST43_NUM_TREES', ternary(CFG.QuickMode, 100, 180));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST43_MIN_LEAF_SIZE', 60);
CFG.AcceptThresholdPoints = env_number('ADAPTIVE_REQ_TEST43_ACCEPT_THRESHOLD_POINTS', 0.25);
CFG.PhysicalRange = [0.5 10];
if CFG.QuickMode
    CFG.BaseModels = ["q_spectrum_plus_composition","q_spectrum_only"];
else
    CFG.BaseModels = ["q_spectrum_plus_composition","q_spectrum_plus_theory_composition", ...
        "q_spectrum_only","delta_q_theory_composition"];
end
CFG.PlotStrategies = ["q_spectrum_plus_composition", ...
    "T43_q_spectrum_plus_composition_logk_all", ...
    "T43_q_spectrum_plus_composition_accept_gate", ...
    "q_spectrum_only", ...
    "T43_q_spectrum_only_logk_all", ...
    "T43_q_spectrum_only_accept_gate", ...
    "theory_discrete", "T34_mixedness_logk_corrected"];
end

function CFG = setup_parallel_if_requested(CFG)
if ~CFG.UseParallel
    CFG.StatOptions = statset('UseParallel', false);
    CFG.LSBoostOptions = statset('UseParallel', false);
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool), parpool('threads'); end
    CFG.StatOptions = statset('UseParallel', true);
    % MATLAB only supports parallel fitting for Bag ensembles. The residual
    % regressor uses LSBoost, so keep that fit serial even when a pool exists.
    CFG.LSBoostOptions = statset('UseParallel', false);
catch ME
    warning('Could not start/use parallel pool (%s). Continuing serial.', ME.message);
    CFG.UseParallel = false;
    CFG.StatOptions = statset('UseParallel', false);
    CFG.LSBoostOptions = statset('UseParallel', false);
end
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_43_synthetic_trained_composition_correction_transfer');
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
    CFG.TrainSource, 'test38_patch_level_predictions.csv');
SRC.test39_csv = table_path(root_dir, 'test_39_frozen_test38_external_validation', ...
    CFG.Test39Source, 'test39_external_patch_level_predictions.csv');
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

function validate_inputs(T38, T39, T40, CFG)
needed = ["model_name","condition_key","M","f0","q_pred","sws_pred","true_SWS", ...
    "q_theory_ref","sws_theory_ref","predicted_patch_purity","p_mixed", ...
    "is_train_row","is_heldout_row"];
assert(all(ismember(needed, string(T38.Properties.VariableNames))), 'Test38 table missing required columns.');
assert(all(ismember(needed(1:end-2), string(T39.Properties.VariableNames))), 'Test39 table missing required columns.');
assert(all(ismember(["model_name","strategy_name","sws_pred","true_SWS","predicted_patch_purity","p_mixed"], ...
    string(T40.Properties.VariableNames))), 'Test40 table missing required columns.');
for m = CFG.BaseModels
    assert(any(T38.model_name == m & T38.is_train_row), 'No train rows for %s.', m);
    assert(any(T38.model_name == m & T38.is_heldout_row), 'No heldout rows for %s.', m);
    assert(any(T39.model_name == m), 'No Test39 rows for %s.', m);
    assert(any(T40.model_name == m), 'No k-Wave rows for %s.', m);
end
assert_no_oracle_features(operational_feature_names());
end

%% Standardization

function T = standardize_synthetic_table(T, dataset_name)
T = convert_strings(T);
T.dataset_eval = repmat(dataset_name, height(T), 1);
T.strategy_name = T.model_name;
T.model_family = repmat("Test38_frozen_q_baseline", height(T), 1);
T.base_model_name = T.model_name;
T.correction_variant = repmat("none", height(T), 1);
T.sws_theory_ref = T.sws_theory;
T.q_theory_ref = T.q_theory_prior;
if ~ismember("p_strong_mixed", string(T.Properties.VariableNames))
    T.p_strong_mixed = zeros(height(T),1);
end
if ~ismember("material_region", string(T.Properties.VariableNames))
    T.material_region = repmat("unknown", height(T), 1);
end
if ~ismember("region_zone", string(T.Properties.VariableNames))
    T.region_zone = T.distance_bin;
end
if ~ismember("material_zone", string(T.Properties.VariableNames))
    T.material_zone = T.material_region + "_" + T.region_zone;
end
T.patch_key = T.condition_key + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
end

function T = standardize_kwave_table(T)
T = convert_strings(T);
T.dataset_eval = repmat("kwave_transfer", height(T), 1);
T.case_id = "kwave_" + sanitize(T.case_name);
T.case_family = repmat("kwave_inclusion", height(T), 1);
T.field_regime_ood = T.field_regime;
T.model_name = erase(T.strategy_name, "T38_");
T.model_name(startsWith(T.strategy_name, "T34_")) = T.strategy_name(startsWith(T.strategy_name, "T34_"));
T.model_family(startsWith(T.strategy_name, "T38_")) = "Test38_frozen_q_baseline";
T.base_model_name = T.model_name;
T.correction_variant = repmat("none", height(T), 1);
T.patch_key = T.condition_key + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
T.sws_theory_ref = lookup_kwave_strategy(T, "T34_theory_baseline");
T.q_theory_ref = 0.5 * ones(height(T),1);
if ~ismember("p_strong_mixed", string(T.Properties.VariableNames))
    T.p_strong_mixed = zeros(height(T),1);
end
T.material_region = T.material_side;
T.region_zone = T.analysis_region;
T.material_zone = T.material_region + "_" + T.region_zone;
T.distance_to_boundary_mm = T.distance_to_interface_mm;
end

function y = lookup_kwave_strategy(T, strategy)
R = T(T.strategy_name == strategy, {'patch_key','sws_pred'});
[ok,loc] = ismember(T.patch_key, R.patch_key);
y = nan(height(T),1);
y(ok) = R.sws_pred(loc(ok));
y(~isfinite(y)) = T.sws_pred(~isfinite(y));
end

function T = convert_strings(T)
for v = ["dataset","dataset_eval","condition_key","case_id","case_family", ...
        "case_name","field_regime","field_regime_ood","purity_bin","distance_bin", ...
        "material_region","region_zone","material_zone","model_name","model_family", ...
        "strategy_name","base_model_name","correction_variant","geometry","geometry_type", ...
        "material_side","analysis_region","roi_name","region_label","bundle_id"]
    if ismember(v, string(T.Properties.VariableNames))
        T.(v) = string(T.(v));
    end
end
end

%% Training/application

function [MODEL, info] = train_correction_stack(T, model_name, CFG)
[X, names] = build_features(T);
assert_no_oracle_features(names);
valid = all(isfinite(X),2) & isfinite(T.sws_pred) & T.sws_pred > 0 & ...
    isfinite(T.true_SWS) & T.true_SWS > 0;
base_k = 2*pi*T.f0 ./ T.sws_pred;
true_k = 2*pi*T.f0 ./ T.true_SWS;
target = log(true_k) - log(base_k);
tree = templateTree('MinLeafSize', CFG.MinLeafSize, 'MaxNumSplits', 384);
delta_model = fitrensemble(X(valid,:), target(valid), ...
    'Method','LSBoost','Learners',tree,'NumLearningCycles',CFG.NumLearningCycles, ...
    'LearnRate',0.05,'Options',CFG.LSBoostOptions);
delta_train = predict(delta_model, X(valid,:));
base_err = T.sws_abs_error_pct(valid);
sws_corr = clip_sws(2*pi*T.f0(valid) ./ exp(log(base_k(valid)) + delta_train), CFG);
corr_err = abs(100*(sws_corr - T.true_SWS(valid))./T.true_SWS(valid));
accept_label = corr_err + CFG.AcceptThresholdPoints < base_err;
gate_features = [X(valid,:) delta_train abs(delta_train) mixedness_gate(T(valid,:))];
gate_names = [names "delta_logk_pred" "abs_delta_logk_pred" "heuristic_mixedness_gate"];
gate_model = fitcensemble(gate_features, accept_label, ...
    'Method','Bag','Learners',tree,'NumLearningCycles',max(50, round(CFG.NumLearningCycles/2)), ...
    'ClassNames',[false true],'Options',CFG.StatOptions);
MODEL = struct('delta_model',delta_model,'gate_model',gate_model, ...
    'feature_names',names,'gate_feature_names',gate_names,'model_name',model_name);
info = table(model_name, nnz(valid), mean(abs(target(valid)),'omitnan'), ...
    mean(base_err,'omitnan'), mean(corr_err,'omitnan'), 100*mean(accept_label,'omitnan'), ...
    'VariableNames', {'model_name','N_train','mean_abs_delta_logk_label', ...
    'train_baseline_MAPE_pct','train_corrected_MAPE_pct','train_accept_label_pct'});
end

function C = apply_correction_stack(T, MODEL, dataset_eval, model_name, CFG)
[X, names] = build_features(T);
assert(isequal(string(names(:)), string(MODEL.feature_names(:))), 'Feature mismatch.');
delta = predict(MODEL.delta_model, X);
base_k = 2*pi*T.f0 ./ T.sws_pred;
sws_all = clip_sws(2*pi*T.f0 ./ exp(log(base_k) + delta), CFG);
heur_gate = mixedness_gate(T);
sws_mixed = clip_sws(2*pi*T.f0 ./ exp(log(base_k) + heur_gate.*delta), CFG);
gateX = [X delta abs(delta) heur_gate];
[accept_cls, score] = predict(MODEL.gate_model, gateX);
p_accept = positive_score(MODEL.gate_model, score);
accept = string(accept_cls) == "true" | string(accept_cls) == "1";
sws_accept = T.sws_pred;
sws_accept(accept) = sws_all(accept);

C1 = corrected_rows(T, dataset_eval, model_name, "logk_all", sws_all, delta, ones(height(T),1), p_accept);
C2 = corrected_rows(T, dataset_eval, model_name, "mixedness_weighted", sws_mixed, delta, heur_gate, p_accept);
C3 = corrected_rows(T, dataset_eval, model_name, "accept_gate", sws_accept, delta, double(accept), p_accept);
C = vertcat_compatible(C1, C2, C3);
end

function B = baseline_rows(T, dataset_eval, model_name)
B = T;
B.dataset_eval = repmat(dataset_eval, height(T), 1);
B.base_model_name = repmat(model_name, height(T), 1);
B.correction_variant = repmat("none", height(T), 1);
if startsWith(model_name, "T34_")
    B.strategy_name = repmat(model_name, height(T), 1);
    B.model_family = repmat("T34_kwave_reference", height(T), 1);
elseif model_name == "theory_discrete"
    B.strategy_name = repmat("theory_discrete", height(T), 1);
    B.model_family = repmat("theory_reference", height(T), 1);
else
    B.strategy_name = repmat(model_name, height(T), 1);
    B.model_family = repmat("Test38_frozen_q_baseline", height(T), 1);
end
B.delta_logk_pred = nan(height(T),1);
B.correction_weight = zeros(height(T),1);
B.p_accept_correction = nan(height(T),1);
end

function C = corrected_rows(T, dataset_eval, model_name, variant, sws, delta, weight, p_accept)
C = baseline_rows(T, dataset_eval, model_name);
C.model_family = repmat("T43_synthetic_trained_composition_correction", height(T), 1);
C.correction_variant = repmat(variant, height(T), 1);
C.strategy_name = repmat("T43_" + model_name + "_" + variant, height(T), 1);
C.sws_pred = sws;
C.sws_signed_error_pct = 100*(sws - C.true_SWS)./C.true_SWS;
C.sws_abs_error_pct = abs(C.sws_signed_error_pct);
C.high_error10 = C.sws_abs_error_pct > 10;
C.high_error20 = C.sws_abs_error_pct > 20;
C.delta_logk_pred = delta;
C.correction_weight = weight;
C.p_accept_correction = p_accept;
end

function [X, names] = build_features(T)
names = operational_feature_names();
q = fill_default(T.q_pred, 0.5);
qth = fill_default(T.q_theory_ref, 0.5);
theory = fill_default(T.sws_theory_ref, T.sws_pred);
purity = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-purity);
ps = fill_default(T.p_strong_mixed, 0);
X = [normalize_like(T.M,3,1), normalize_like(T.f0,500,150), regime_code(T), ...
    purity, 1-purity, pm, ps, log(max(T.sws_pred,0.1)), q, ...
    log(max(theory,0.1)), qth, ...
    log(max(T.sws_pred,0.1))-log(max(theory,0.1)), q-qth, ...
    abs(T.sws_pred-theory)./max(theory,0.25)];
X = double(X); X(~isfinite(X)) = 0;
end

function names = operational_feature_names()
names = ["M_norm","f0_norm","field_regime_code","predicted_patch_purity", ...
    "predicted_impurity","p_mixed","p_strong_mixed","log_sws_base", ...
    "q_pred","log_sws_theory_ref","q_theory_ref","log_base_minus_theory", ...
    "q_minus_theory","rel_base_theory_disagreement"];
end

function gate = mixedness_gate(T)
purity = fill_default(T.predicted_patch_purity, 0.95);
pm = fill_default(T.p_mixed, 1-purity);
ps = fill_default(T.p_strong_mixed, 0);
gate = max([pm, 1-purity, 0.75*ps], [], 2);
gate = min(max(gate,0),1);
end

function assert_no_oracle_features(names)
bad = ["true","oracle","material","side","roi","distance", ...
    "error","signed","abs_error","high_error","cs_true","sws_true","k_true"];
for b = bad
    hit = names(contains(lower(names), b));
    assert(isempty(hit), 'Oracle/evaluation feature leaked into Test43 correction: %s', strjoin(hit, ', '));
end
end

%% Summaries/plots

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
n = splitapply(@numel, T.sws_abs_error_pct, G);
mape = splitapply(@(x) mean(x,'omitnan'), T.sws_abs_error_pct, G);
signed = splitapply(@(x) mean(x,'omitnan'), T.sws_signed_error_pct, G);
under = splitapply(@(x) 100*mean(x<0,'omitnan'), T.sws_signed_error_pct, G);
he20 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
mean_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
mean_acc = splitapply(@(x) mean(x,'omitnan'), T.p_accept_correction, G);
S = [groups table(n,mape,signed,under,he20,mean_sws,mean_acc, ...
    'VariableNames', {'N','MAPE_pct','mean_signed_error_pct','underestimate_pct', ...
    'high_error20_pct','mean_sws_pred','mean_p_accept_correction'})];
end

function T = correction_gain_summary(P, CFG)
T = table();
for dataset_eval = unique(P.dataset_eval, 'stable')'
    for base = CFG.BaseModels(:)'
        B = P(P.dataset_eval == dataset_eval & P.strategy_name == base,:);
        if isempty(B), continue; end
        for variant = ["logk_all","mixedness_weighted","accept_gate"]
            s = "T43_" + base + "_" + variant;
            C = P(P.dataset_eval == dataset_eval & P.strategy_name == s,:);
            if isempty(C), continue; end
            [ok,loc] = ismember(C.patch_key, B.patch_key);
            R = C(ok, {'dataset_eval','strategy_name','case_family','material_region','region_zone','purity_bin'});
            R.mape_gain_points = B.sws_abs_error_pct(loc(ok)) - C.sws_abs_error_pct(ok);
            T = concat_tables(T, R);
        end
    end
end
if isempty(T), return; end
[G,groups] = findgroups(T(:, {'dataset_eval','strategy_name','case_family','material_region','region_zone','purity_bin'}));
T = [groups table(splitapply(@numel,T.mape_gain_points,G), ...
    splitapply(@(x) mean(x,'omitnan'),T.mape_gain_points,G), ...
    splitapply(@(x) 100*mean(x>0,'omitnan'),T.mape_gain_points,G), ...
    'VariableNames', {'N','mean_mape_gain_points','pct_pixels_improved'})];
end

function T = top_failure_conditions(T, n)
T = sortrows(T, {'MAPE_pct','high_error20_pct'}, {'descend','descend'});
T = T(1:min(n,height(T)),:);
end

function plot_dataset_ranking(T, OUT)
datasets = unique(T.dataset_eval, 'stable');
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tl = tiledlayout(fig, numel(datasets), 1, 'TileSpacing','compact','Padding','compact');
for d = datasets(:)'
    ax = nexttile(tl);
    X = sortrows(T(T.dataset_eval == d,:), 'MAPE_pct', 'ascend');
    barh(ax, categorical(X.strategy_name), X.MAPE_pct);
    set(ax,'TickLabelInterpreter','none'); xlabel(ax,'MAPE (%)'); title(ax,d,'Interpreter','none'); grid(ax,'on');
end
export_fig(fig, fullfile(OUT.figure_dir, 'test43_dataset_model_ranking.png'));
end

function plot_transfer_summary(T, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
strategies = unique(T.strategy_name, 'stable');
grouped_bar(axes(fig), T, "dataset_eval", strategies, "MAPE_pct");
ylabel('MAPE (%)'); title('Transfer summary by dataset','FontWeight','normal');
legend('Location','bestoutside','Interpreter','none'); grid on;
export_fig(fig, fullfile(OUT.figure_dir, 'test43_transfer_summary.png'));
end

function plot_region_diagnostics(Tmat, Tzone, CFG, OUT)
keep = select_plot_strategies([Tmat.strategy_name; Tzone.strategy_name], CFG);
fig = figure('Color','w','Units','centimeters','Position',[1 1 35 20]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for ds = ["test39_external","kwave_transfer"]
    ax = nexttile(tl);
    X = Tmat(Tmat.dataset_eval == ds & ismember(Tmat.strategy_name, keep),:);
    grouped_bar(ax, X, "material_region", keep, "MAPE_pct");
    ylabel(ax,'MAPE (%)'); title(ax, ds + " material MAPE", 'Interpreter','none'); grid(ax,'on');
    legend(ax,'Location','bestoutside','Interpreter','none');
    ax = nexttile(tl);
    Z = Tzone(Tzone.dataset_eval == ds & ismember(Tzone.strategy_name, keep),:);
    grouped_bar(ax, Z, "region_zone", keep, "MAPE_pct");
    ylabel(ax,'MAPE (%)'); title(ax, ds + " zone MAPE", 'Interpreter','none'); grid(ax,'on');
    legend(ax,'Location','bestoutside','Interpreter','none');
end
export_fig(fig, fullfile(OUT.figure_dir, 'test43_region_diagnostics.png'));
end

function plot_heatmap_table(T, group_var, OUT, filename)
datasets = unique(T.dataset_eval, 'stable');
for ds = datasets(:)'
    X = T(T.dataset_eval == ds,:);
    strategies = unique(X.strategy_name, 'stable');
    groups = unique(string(X.(group_var)), 'stable');
    Z = nan(numel(strategies), numel(groups));
    for i = 1:numel(strategies)
        for j = 1:numel(groups)
            idx = X.strategy_name == strategies(i) & string(X.(group_var)) == groups(j);
            if any(idx), Z(i,j) = mean(X.MAPE_pct(idx),'omitnan'); end
        end
    end
    fig = figure('Color','w','Units','centimeters','Position',[1 1 30 16]);
    imagesc(Z); colorbar; colormap(parula);
    set(gca,'XTick',1:numel(groups),'XTickLabel',groups, ...
        'YTick',1:numel(strategies),'YTickLabel',strategies,'TickLabelInterpreter','none');
    xtickangle(35); title(ds + " MAPE by " + group_var, 'Interpreter','none','FontWeight','normal');
    export_fig(fig, fullfile(OUT.figure_dir, erase(filename,'.png') + "__" + sanitize(ds) + ".png"));
end
end

function plot_correction_gain(T, OUT)
if isempty(T), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 16]);
X = T(ismember(T.region_zone, ["interface_0_2mm","near_interface_2_4mm","core_gt4mm","hard_core_gt4mm"]),:);
strategies = unique(X.strategy_name, 'stable');
grouped_bar(axes(fig), X, "dataset_eval", strategies, "mean_mape_gain_points");
yline(0,'k-'); ylabel('MAPE gain points'); title('Correction gain by dataset');
legend('Location','bestoutside','Interpreter','none'); grid on;
export_fig(fig, fullfile(OUT.figure_dir, 'test43_correction_gain_by_dataset.png'));
end

function plot_condition_maps(T, CFG, OUT)
keep = select_plot_strategies(T.strategy_name, CFG);
conds = unique(T(:, {'dataset_eval','case_id','field_regime_ood','f0','M'}), 'rows', 'stable');
if isfinite(CFG.MaxMapConditions) && height(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions,:);
end
fprintf('Saving %d Test43 map panels under %s.\n', height(conds), OUT.map_dir);
for ci = 1:height(conds)
    idx = T.dataset_eval == conds.dataset_eval(ci) & T.case_id == conds.case_id(ci) & ...
        T.field_regime_ood == conds.field_regime_ood(ci) & T.f0 == conds.f0(ci) & T.M == conds.M(ci);
    Xc = T(idx,:);
    base = Xc(Xc.strategy_name == Xc.strategy_name(1),:);
    if isempty(base), continue; end
    [true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 36 24]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS');
    plot_map(nexttile(tl), rows_to_grid(base, base.predicted_patch_purity, nz, nx), 'Predicted purity');
    plot_map(nexttile(tl), rows_to_grid(base, base.p_mixed, nz, nx), 'P(mixed)');
    for si = 1:numel(keep)
        s = keep(si); X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), s + ' SWS');
    end
    err_keep = keep(1:min(5,numel(keep)));
    for si = 1:numel(err_keep)
        s = err_keep(si); X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, nz, nx), s + ' err %');
    end
    title(tl, sprintf('%s | %s | %s | f%d | M%d', conds.dataset_eval(ci), ...
        conds.case_id(ci), conds.field_regime_ood(ci), conds.f0(ci), conds.M(ci)), 'Interpreter','none');
    outdir = fullfile(OUT.map_dir, sanitize(conds.dataset_eval(ci)), sanitize(conds.case_id(ci)), ...
        sanitize(conds.field_regime_ood(ci)), "f"+string(conds.f0(ci)), "M"+string(conds.M(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    export_fig(fig, fullfile(outdir, sprintf('test43_map__%s__%s__%s__f%d__M%d.png', ...
        sanitize(conds.dataset_eval(ci)), sanitize(conds.case_id(ci)), ...
        sanitize(conds.field_regime_ood(ci)), conds.f0(ci), conds.M(ci))));
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
bar(ax, categorical(groups), Y); set(ax,'TickLabelInterpreter','none');
end

function keep = select_plot_strategies(strategies, CFG)
available = unique(string(strategies), 'stable');
preferred = string(CFG.PlotStrategies(:));
keep = preferred(ismember(preferred, available));
if isempty(keep), keep = available(1:min(8,numel(available))); end
end

function print_interpretation(T_overall, T_by_dataset, Tmat, Tzone, Tgain, Tfail)
fprintf('\n================ Test 43 transfer summary ================\n');
for ds = unique(T_overall.dataset_eval, 'stable')'
    X = sortrows(T_overall(T_overall.dataset_eval == ds,:), 'MAPE_pct', 'ascend');
    fprintf('%s best: %s (MAPE %.2f%%, signed %.2f%%, HE20 %.2f%%).\n', ...
        ds, X.strategy_name(1), X.MAPE_pct(1), X.mean_signed_error_pct(1), X.high_error20_pct(1));
end
disp(sortrows(T_by_dataset, {'dataset_eval','MAPE_pct'}));
if ~isempty(Tgain)
    G = sortrows(Tgain, 'mean_mape_gain_points', 'descend');
    fprintf('Largest gain: %s / %s / %s / %s = %.2f points.\n', ...
        G.dataset_eval(1), G.strategy_name(1), G.material_region(1), G.region_zone(1), G.mean_mape_gain_points(1));
end
fprintf('Worst conditions saved: %d rows.\n', height(Tfail));
fprintf('==========================================================\n');
end

%% Utilities

function y = fill_default(x, d)
y = x; bad = ~isfinite(y);
if isscalar(d), y(bad) = d; else, y(bad) = d(bad); end
end

function y = clip_sws(y, CFG)
y = min(max(y, CFG.PhysicalRange(1)), CFG.PhysicalRange(2));
end

function y = normalize_like(x, mu, sigma)
y = (double(x)-mu)./sigma;
end

function code = regime_code(T)
if ismember("field_regime_ood", string(T.Properties.VariableNames)), r = T.field_regime_ood; else, r = T.field_regime; end
r = string(r); code = zeros(numel(r),1);
code(contains(r,"directional")) = 1; code(contains(r,"diffuse_2D")) = 2;
code(contains(r,"partial")) = 3; code(contains(r,"diffuse_3D")) = 4;
code = normalize_like(code,2.5,1.2);
end

function p = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes=="true" | classes=="1", 1);
if isempty(idx), idx = size(score,2); end
p = score(:,idx);
end

function T = vertcat_compatible(varargin)
vars = strings(1,0);
for i = 1:nargin, vars = [vars string(varargin{i}.Properties.VariableNames)]; end %#ok<AGROW>
vars = unique(vars,'stable');
string_vars = ["dataset","dataset_eval","condition_key","case_id","case_family","case_name", ...
    "field_regime","field_regime_ood","purity_bin","distance_bin","material_region", ...
    "region_zone","material_zone","model_name","model_family","strategy_name", ...
    "base_model_name","correction_variant","patch_key","geometry","geometry_type", ...
    "material_side","analysis_region","roi_name","region_label","bundle_id"];
logical_vars = ["is_train_row","is_heldout_row","is_mixed","is_strong_mixed","high_error10","high_error20"];
parts = cell(nargin,1);
for i = 1:nargin
    A = varargin{i};
    for v = vars
        if ismember(v,string(A.Properties.VariableNames)), continue; end
        if ismember(v,string_vars), A.(v) = repmat("unknown",height(A),1);
        elseif ismember(v,logical_vars), A.(v) = false(height(A),1);
        else, A.(v) = nan(height(A),1);
        end
    end
    parts{i} = A(:,cellstr(vars));
end
T = vertcat(parts{:});
end

function T = concat_tables(T, R)
if isempty(T), T = R; else, T = vertcat_compatible(T, R); end
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3, nz = max(T.map_iz); nx = max(T.map_ix); end
Z = nan(nz,nx); Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map(ax, Z, ttl)
imagesc(ax,Z); axis(ax,'image'); axis(ax,'off'); colorbar(ax);
title(ax, ttl, 'Interpreter','none','FontWeight','normal','FontSize',7);
end

function export_fig(fig, file)
drawnow; exportgraphics(fig,file,'Resolution',220); close(fig);
end

function safe_plot(fn, label)
try, fn(); catch ME, warning('Skipping plot %s: %s', label, ME.message); end
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
end

function tf = env_true(name, default_value)
v = lower(strtrim(string(getenv(name))));
if v=="", tf = default_value; else, tf = ismember(v,["1","true","yes","y","on"]); end
end

function x = env_number(name, default_value)
v = strtrim(string(getenv(name)));
if v=="", x = default_value; else, x = str2double(v); if ~isfinite(x), x = default_value; end, end
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w'); assert(fid>0, 'Could not write config: %s', file);
fprintf(fid,'%s\n',txt); fclose(fid);
end
