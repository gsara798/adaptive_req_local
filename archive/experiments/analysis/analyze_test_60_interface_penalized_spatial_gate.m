%% analyze_test_60_interface_penalized_spatial_gate.m
% Test 60: interface-penalized spatial gate.
%
% This analysis combines:
%   - Test56-style M-effective residual correction.
%   - Test57-style M=2 versus M=3 multi-window candidates.
%   - Test58-style high-error risk/estimability gate.
%
% No frozen base q/SWS model is retrained. True SWS/k/error labels are used
% only for supervised residual/selector/risk labels and evaluation.

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
write_config_json(CFG, fullfile(OUT.root_dir,'test60_configuration.json'));

fprintf('\nTest 60: interface-penalized M-eff + spatial gate\n');
fprintf('Mode: %s | source: %s\n', CFG.Mode, CFG.SourceCsv);
fprintf('Base model: %s | M2 baseline, M3 hard candidate.\n', CFG.BaseModel);
fprintf('No base q model retraining. Oracle variables are labels/evaluation only.\n');

T0 = readtable(CFG.SourceCsv);
T0 = T0(string(T0.model_name) == CFG.BaseModel, :);
assert(~isempty(T0), 'No rows found for base model %s.', CFG.BaseModel);
T0 = add_operational_features_long(T0, CFG);
assert_no_forbidden_predictors([CFG.ResidualFeatureVars; CFG.RiskFeatureVars; CFG.SelectorFeatureVars]);

if CFG.ValidateOnly
    T0 = sample_by_base_condition(T0, CFG.ValidateConditions, CFG.RandomSeed);
end

Mvals = unique(T0.M)';
fprintf('Available M values in source: %s\n', mat2str(Mvals));
assert(any(abs(Mvals-2)<1e-8), 'Test60 requires M=2 rows.');
assert(any(abs(Mvals-3)<1e-8), 'Test60 requires M=3 rows. Rerun Test55 with ADAPTIVE_REQ_TEST55_M_LIST including 3.');

W = build_m2_m3_wide(T0, CFG);
fprintf('Aligned M2/M3 pixels: %d\n', height(W));
assert(height(W) > 0, 'No aligned M2/M3 rows found.');

[train_mask, test_mask] = condition_split(W, CFG);
MODEL = train_models(W, train_mask, CFG);
T_strategies = apply_strategies(W(test_mask,:), MODEL, CFG);

T_overall = summarize_predictions(T_strategies, "strategy_name");
T_by_case = summarize_predictions(T_strategies, ["strategy_name","case_id"]);
T_by_freq = summarize_predictions(T_strategies, ["strategy_name","f0"]);
T_by_regime = summarize_predictions(T_strategies, ["strategy_name","field_regime_ood"]);
T_by_roi = summarize_predictions(T_strategies(T_strategies.roi_region ~= "other",:), ...
    ["strategy_name","roi_region"]);
T_by_case_roi = summarize_predictions(T_strategies(T_strategies.roi_region ~= "other",:), ...
    ["strategy_name","case_id","roi_region"]);
T_harm = summarize_harm(T_strategies);
T_success = success_table(T_overall, T_by_roi, CFG);

writetable(remove_cell_columns(T_strategies), fullfile(OUT.table_dir,'test60_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test60_strategy_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test60_strategy_summary_by_case.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir,'test60_strategy_summary_by_frequency.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test60_strategy_summary_by_regime.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir,'test60_strategy_summary_by_roi.csv'));
writetable(T_by_case_roi, fullfile(OUT.table_dir,'test60_strategy_summary_by_case_roi.csv'));
writetable(T_harm, fullfile(OUT.table_dir,'test60_correction_harm_summary.csv'));
writetable(T_success, fullfile(OUT.table_dir,'test60_success_criteria_summary.csv'));

safe_plot(@() plot_summary_figures(T_overall, T_by_case_roi, T_by_freq, T_by_roi, T_harm, OUT, CFG), ...
    'summary figures');
if CFG.SaveAllMaps
    safe_plot(@() plot_all_condition_maps(T_strategies, OUT, CFG), 'condition maps');
end

save(fullfile(OUT.data_dir,'test60_interface_penalized_spatial_gate_results.mat'), ...
    'CFG','MODEL','T_overall','T_by_case','T_by_freq','T_by_regime', ...
    'T_by_roi','T_by_case_roi','T_harm','T_success','-v7.3');

print_summary(T_overall, T_by_roi, T_harm, T_success, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 60 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST60_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST60_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST60_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST60_SAVE_ALL_MAPS', true);
CFG.RandomSeed = 59001;
CFG.SourceCsv = env_string('ADAPTIVE_REQ_TEST60_SOURCE_CSV', ...
    fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep', ...
    'tables','test55_patch_level_results.csv'));
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST60_BASE_MODEL', ...
    "q_spectrum_plus_composition");
CFG.TrainFraction = env_number('ADAPTIVE_REQ_TEST60_TRAIN_FRACTION', 0.70);
CFG.CsGuess = env_number('ADAPTIVE_REQ_TEST60_CS_GUESS', 3.0);
CFG.TreeLearners = env_number('ADAPTIVE_REQ_TEST60_TREE_LEARNERS', ternary(mode=="quick",80,180));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST60_MIN_LEAF_SIZE', 10);
CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST60_MAX_TRAIN_ROWS', ternary(mode=="quick",60000,250000));
CFG.ValidateConditions = env_number('ADAPTIVE_REQ_TEST60_VALIDATE_CONDITIONS', 16);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST60_USE_PARALLEL', false);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST60_MAX_MAP_CONDITIONS', 0);

CFG.RiskThresholdHard = env_number('ADAPTIVE_REQ_TEST60_HARD_RISK_THRESHOLD', 0.20);
CFG.RiskThresholdSelector = env_number('ADAPTIVE_REQ_TEST60_SELECTOR_RISK_THRESHOLD', 0.50);
CFG.PurityThreshold = env_number('ADAPTIVE_REQ_TEST60_PURITY_THRESHOLD', 0.90);
CFG.MixedThreshold = env_number('ADAPTIVE_REQ_TEST60_MIXED_THRESHOLD', 0.25);
CFG.HardSwsRatio = env_number('ADAPTIVE_REQ_TEST60_HARD_SWS_RATIO', 1.05);
CFG.MaxM3RelDiff = env_number('ADAPTIVE_REQ_TEST60_MAX_M3_REL_DIFF', 0.20);
CFG.MaxResidualCorrectionPct = env_number('ADAPTIVE_REQ_TEST60_MAX_RESIDUAL_CORRECTION_PCT', 10);
CFG.LabelGainMarginPct = env_number('ADAPTIVE_REQ_TEST60_LABEL_GAIN_MARGIN_PCT', 0.50);
CFG.AlignToleranceM = env_number('ADAPTIVE_REQ_TEST60_ALIGN_TOLERANCE_M', 0.40e-3);
CFG.InterfaceRiskThreshold = env_number('ADAPTIVE_REQ_TEST60_INTERFACE_RISK_THRESHOLD', 0.30);
CFG.InterfaceMixedThreshold = env_number('ADAPTIVE_REQ_TEST60_INTERFACE_MIXED_THRESHOLD', 0.50);
CFG.InterfacePurityThreshold = env_number('ADAPTIVE_REQ_TEST60_INTERFACE_PURITY_THRESHOLD', 0.80);
CFG.InterfaceExtraGainMarginPct = env_number('ADAPTIVE_REQ_TEST60_INTERFACE_EXTRA_GAIN_MARGIN_PCT', 2.0);
CFG.PenalizedSelectorRiskThreshold = env_number('ADAPTIVE_REQ_TEST60_PENALIZED_SELECTOR_RISK_THRESHOLD', 0.35);
CFG.SpatialGateRadiusPx = env_number('ADAPTIVE_REQ_TEST60_SPATIAL_GATE_RADIUS_PX', 1);
CFG.SpatialGateMinSupport = env_number('ADAPTIVE_REQ_TEST60_SPATIAL_GATE_MIN_SUPPORT', 4);

base_predictors = ["q";"sws";"log_sws";"k";"logk";"M_eff";"lambda_m"; ...
    "window_length_m";"predicted_patch_purity";"p_mixed";"p_strong_mixed"; ...
    "M";"f0";"dx";"dz";"field_code";"case_family_code"];
CFG.ResidualFeatureVars = "M2_" + base_predictors;
CFG.RiskFeatureVars = "M2_" + base_predictors([1:6 9:17]);
CFG.SelectorFeatureVars = ["M2_q";"M2_sws";"M2_M_eff";"M2_predicted_patch_purity"; ...
    "M2_p_mixed";"M2_p_strong_mixed";"M2_risk_high20"; ...
    "M3_q";"M3_sws";"M3_M_eff";"M3_predicted_patch_purity"; ...
    "M3_p_mixed";"M3_p_strong_mixed";"M3_risk_high20"; ...
    "sws_M3_minus_M2";"rel_sws_M3_minus_M2"; ...
    "f0";"dx";"dz";"field_code";"case_family_code"];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_60_interface_penalized_spatial_gate','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_60_interface_penalized_spatial_gate');
end
if CFG.ValidateOnly
    OUT.root_dir = fullfile(OUT.root_dir,'validate');
end
OUT.table_dir = fullfile(OUT.root_dir,'tables');
OUT.figure_dir = fullfile(OUT.root_dir,'figures');
OUT.map_dir = fullfile(OUT.figure_dir,'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir,'data');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

%% Data preparation

function T = add_operational_features_long(T, CFG)
T.sws_pred = clamp(T.sws_pred, 0.2, 20);
T.q_pred = clamp(T.q_pred, 0, 1);
T.k_pred_firstpass = 2*pi*T.f0 ./ T.sws_pred;
T.logk_pred_firstpass = log(T.k_pred_firstpass);
T.log_sws_pred = log(T.sws_pred);
T.lambda_pred_m = T.sws_pred ./ T.f0;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.M_eff_pred = T.window_length_m ./ T.lambda_pred_m;
T.field_code = categorical_code(string(T.field_regime_ood));
T.case_family_code = categorical_code(string(T.case_family));
T.base_condition = regexprep(string(T.condition_key), "__M[0-9.]+$", "");
if all(ismember(["x_center_m","z_center_m"], string(T.Properties.VariableNames)))
    T.x_um = round(1e6*T.x_center_m);
    T.z_um = round(1e6*T.z_center_m);
else
    T.x_um = T.map_ix;
    T.z_um = T.map_iz;
end
T.align_key = T.base_condition + "__xum" + string(T.x_um) + "__zum" + string(T.z_um);
if ~ismember("roi_region", string(T.Properties.VariableNames))
    T.roi_region = repmat("unknown", height(T), 1);
end
if ~ismember("material_region", string(T.Properties.VariableNames))
    T.material_region = repmat("unknown", height(T), 1);
end
end

function W = build_m2_m3_wide(T0, CFG)
% Use the M=3 centers as the evaluation grid and attach the nearest M=2
% prediction within a physical tolerance. Exact center matching silently
% drops many frequencies because StepX/StepZ differs across M/f0.
T2all = sortrows(T0(abs(T0.M-2)<1e-8,:), {'base_condition','x_um','z_um'});
T3all = sortrows(T0(abs(T0.M-3)<1e-8,:), {'base_condition','x_um','z_um'});
conds = intersect(unique(string(T2all.base_condition), 'stable'), ...
    unique(string(T3all.base_condition), 'stable'), 'stable');

parts = cell(numel(conds), 1);
tol_um = max(1, round(1e6*CFG.AlignToleranceM));
for ci = 1:numel(conds)
    T2c = T2all(string(T2all.base_condition)==conds(ci), :);
    T3c = T3all(string(T3all.base_condition)==conds(ci), :);
    if isempty(T2c) || isempty(T3c)
        continue;
    end

    x2 = unique(T2c.x_um);
    z2 = unique(T2c.z_um);
    [~, ix] = min(abs(double(T3c.x_um) - double(x2(:))'), [], 2);
    [~, iz] = min(abs(double(T3c.z_um) - double(z2(:))'), [], 2);
    match_x = x2(ix);
    match_z = z2(iz);
    ok = abs(double(T3c.x_um) - double(match_x)) <= tol_um & ...
        abs(double(T3c.z_um) - double(match_z)) <= tol_um;

    key2 = string(T2c.x_um) + "__" + string(T2c.z_um);
    key_match = string(match_x) + "__" + string(match_z);
    [tf, loc] = ismember(key_match, key2);
    ok = ok & tf;
    if any(ok)
        parts{ci} = assemble_m2_m3_wide_rows(T2c(loc(ok), :), T3c(ok, :));
    end
end

parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    W = table();
    warning('No M=2/M=3 rows could be aligned within %.3f mm.', 1e3*CFG.AlignToleranceM);
    return;
end
W = vertcat(parts{:});
fprintf('Aligned %d M=3 evaluation centers with nearest M=2 centers (tolerance %.3f mm).\n', ...
    height(W), 1e3*CFG.AlignToleranceM);
end

function W = assemble_m2_m3_wide_rows(T2, T3)
keep_meta = ["dataset","base_condition","condition_key","align_key","case_id", ...
    "case_family","seen_status","field_regime","field_regime_ood","f0","dx","dz", ...
    "map_ix","map_iz","x_center_m","z_center_m","x_um","z_um", ...
    "true_SWS","k_true","patch_purity","purity_bin","distance_to_boundary_mm", ...
    "distance_bin","material_region","roi_region","field_code","case_family_code"];
keep_meta = intersect(keep_meta, string(T3.Properties.VariableNames), 'stable');
W = T3(:, cellstr(keep_meta));
W.condition_key_M2 = string(T2.condition_key);
W.condition_key_M3 = string(T3.condition_key);
W.align_key_M2 = string(T2.align_key);
W.align_key_M3 = string(T3.align_key);
W.M2_match_dx_m = T2.x_center_m - T3.x_center_m;
W.M2_match_dz_m = T2.z_center_m - T3.z_center_m;
W.M2_match_distance_m = hypot(W.M2_match_dx_m, W.M2_match_dz_m);

W.M2_q = T2.q_pred;
W.M2_sws = T2.sws_pred;
W.M2_log_sws = T2.log_sws_pred;
W.M2_k = T2.k_pred_firstpass;
W.M2_logk = T2.logk_pred_firstpass;
W.M2_M_eff = T2.M_eff_pred;
W.M2_lambda_m = T2.lambda_pred_m;
W.M2_window_length_m = T2.window_length_m;
W.M2_predicted_patch_purity = T2.predicted_patch_purity;
W.M2_p_mixed = T2.p_mixed;
W.M2_p_strong_mixed = T2.p_strong_mixed;
W.M2_M = T2.M;
W.M2_f0 = T2.f0;
W.M2_dx = T2.dx;
W.M2_dz = T2.dz;
W.M2_field_code = T2.field_code;
W.M2_case_family_code = T2.case_family_code;

W.M3_q = T3.q_pred;
W.M3_sws = T3.sws_pred;
W.M3_log_sws = T3.log_sws_pred;
W.M3_k = T3.k_pred_firstpass;
W.M3_logk = T3.logk_pred_firstpass;
W.M3_M_eff = T3.M_eff_pred;
W.M3_lambda_m = T3.lambda_pred_m;
W.M3_window_length_m = T3.window_length_m;
W.M3_predicted_patch_purity = T3.predicted_patch_purity;
W.M3_p_mixed = T3.p_mixed;
W.M3_p_strong_mixed = T3.p_strong_mixed;
W.M3_M = T3.M;
W.M3_f0 = T3.f0;
W.M3_dx = T3.dx;
W.M3_dz = T3.dz;
W.M3_field_code = T3.field_code;
W.M3_case_family_code = T3.case_family_code;

W.sws_M3_minus_M2 = W.M3_sws - W.M2_sws;
W.rel_sws_M3_minus_M2 = W.sws_M3_minus_M2 ./ max(W.M2_sws, eps);
W.high_error10_M2 = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS) > 10;
W.high_error20_M2 = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS) > 20;
end

function [train_mask, test_mask] = condition_split(W, CFG)
keys = unique(string(W.base_condition), 'stable');
rng(CFG.RandomSeed);
keys = keys(randperm(numel(keys)));
ntrain = max(1, round(CFG.TrainFraction*numel(keys)));
train_keys = keys(1:ntrain);
train_mask = ismember(string(W.base_condition), train_keys);
test_mask = ~train_mask;
end

%% Models

function MODEL = train_models(W, train_mask, CFG)
MODEL = struct();
MODEL.cfg = CFG;
MODEL.residual = train_residual_model(W, train_mask, CFG);
W_with_risk = W;
[MODEL.risk, W_with_risk.M2_risk_high20, W_with_risk.M3_risk_high20] = ...
    train_and_apply_risk_model(W, train_mask, CFG);
MODEL.selector = train_selector_model(W_with_risk, train_mask, MODEL.residual, CFG);
MODEL.interface_penalized_selector = train_interface_penalized_selector_model( ...
    W_with_risk, train_mask, MODEL.residual, CFG);
end

function residual = train_residual_model(W, train_mask, CFG)
target = log(W.k_true) - log(W.M2_k);
X = W(:, cellstr(CFG.ResidualFeatureVars));
valid = train_mask & all(isfinite(table2array(X)),2) & isfinite(target);
idx = downsample_idx(find(valid), CFG.MaxTrainRows, CFG.RandomSeed);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
residual = fitrensemble(W(idx, cellstr(CFG.ResidualFeatureVars)), target(idx), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template,'Options',opts);
fprintf('Trained M-eff residual on %d rows.\n', numel(idx));
end

function [risk_model, risk_M2, risk_M3] = train_and_apply_risk_model(W, train_mask, CFG)
X = W(:, cellstr(CFG.RiskFeatureVars));
valid = train_mask & all(isfinite(table2array(X)),2) & isfinite(W.high_error20_M2);
idx = downsample_idx(find(valid), CFG.MaxTrainRows, CFG.RandomSeed + 7);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
risk_model = fitcensemble(W(idx, cellstr(CFG.RiskFeatureVars)), logical(W.high_error20_M2(idx)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template,'Options',opts);
risk_M2 = predict_risk(risk_model, W, "M2", CFG);
risk_M3 = predict_risk(risk_model, W, "M3", CFG);
fprintf('Trained high-error risk model on %d rows.\n', numel(idx));
end

function selector = train_selector_model(W, train_mask, residual_model, CFG)
delta = predict(residual_model, W(:, cellstr(CFG.ResidualFeatureVars)));
delta = clamp(delta, -0.45, 0.45);
sws_resid = apply_delta_logk(W.M2_k, W.f0, delta);

err_base = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
err_m3 = abs(100*(W.M3_sws - W.true_SWS)./W.true_SWS);
err_resid = abs(100*(sws_resid - W.true_SWS)./W.true_SWS);

labels = repmat("keep_M2", height(W), 1);
best_err = err_base;
improve_m3 = err_m3 + CFG.LabelGainMarginPct < best_err;
labels(improve_m3) = "switch_M3";
best_err(improve_m3) = err_m3(improve_m3);
improve_resid = err_resid + CFG.LabelGainMarginPct < best_err;
labels(improve_resid) = "apply_residual";

valid = train_mask & all(isfinite(table2array(W(:, cellstr(CFG.SelectorFeatureVars)))),2);
idx = downsample_idx(find(valid), CFG.MaxTrainRows, CFG.RandomSeed + 11);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
selector = fitcensemble(W(idx, cellstr(CFG.SelectorFeatureVars)), categorical(labels(idx)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template,'Options',opts);
fprintf('Trained candidate selector on %d rows.\n', numel(idx));
end

function selector = train_interface_penalized_selector_model(W, train_mask, residual_model, CFG)
delta = predict(residual_model, W(:, cellstr(CFG.ResidualFeatureVars)));
delta = clamp(delta, -0.45, 0.45);
sws_resid = apply_delta_logk(W.M2_k, W.f0, delta);

err_base = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
err_m3 = abs(100*(W.M3_sws - W.true_SWS)./W.true_SWS);
err_resid = abs(100*(sws_resid - W.true_SWS)./W.true_SWS);

interface_like = operational_interface_like(W, CFG);
hard_like = operational_hard_pure_like(W, CFG);
extra_margin = CFG.InterfaceExtraGainMarginPct * double(interface_like & ~hard_like);
margin = CFG.LabelGainMarginPct + extra_margin;

labels = repmat("keep_M2", height(W), 1);
best_err = err_base;
improve_m3 = err_m3 + margin < best_err;
labels(improve_m3) = "switch_M3";
best_err(improve_m3) = err_m3(improve_m3);
improve_resid = err_resid + margin < best_err;
labels(improve_resid) = "apply_residual";

valid = train_mask & all(isfinite(table2array(W(:, cellstr(CFG.SelectorFeatureVars)))),2);
idx = downsample_idx(find(valid), CFG.MaxTrainRows, CFG.RandomSeed + 23);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
selector = fitcensemble(W(idx, cellstr(CFG.SelectorFeatureVars)), categorical(labels(idx)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template,'Options',opts);
fprintf('Trained interface-penalized candidate selector on %d rows.\n', numel(idx));
end

%% Strategy application

function T = apply_strategies(W, MODEL, CFG)
W.M2_risk_high20 = predict_risk(MODEL.risk, W, "M2", CFG);
W.M3_risk_high20 = predict_risk(MODEL.risk, W, "M3", CFG);

delta = predict(MODEL.residual, W(:, cellstr(CFG.ResidualFeatureVars)));
delta = clamp(delta, -0.45, 0.45);
sws_resid = apply_delta_logk(W.M2_k, W.f0, delta);
resid_correction_pct = 100*(sws_resid - W.M2_sws)./W.M2_sws;

hard_pure_gate = W.M2_risk_high20 < CFG.RiskThresholdHard & ...
    W.M2_predicted_patch_purity >= CFG.PurityThreshold & ...
    W.M2_p_mixed <= CFG.MixedThreshold & ...
    W.M2_sws >= CFG.HardSwsRatio * CFG.CsGuess;
m3_stable = W.M3_sws >= W.M2_sws & abs(W.M3_sws - W.M2_sws)./max(W.M2_sws,eps) <= CFG.MaxM3RelDiff;
resid_stable = sws_resid >= W.M2_sws & abs(resid_correction_pct) <= CFG.MaxResidualCorrectionPct;

rule_mask = hard_pure_gate & m3_stable;
resid_mask = hard_pure_gate & resid_stable;

T = [
    make_strategy(W, "fixed_M2", W.M2_sws, W.M2_q, 2*ones(height(W),1), ...
        zeros(height(W),1), false(height(W),1), "none", W.M2_risk_high20, W.M2_risk_high20, false(height(W),1));
    make_strategy(W, "fixed_M3", W.M3_sws, W.M3_q, 3*ones(height(W),1), ...
        zeros(height(W),1), true(height(W),1), "switch_M3", W.M2_risk_high20, W.M3_risk_high20, false(height(W),1));
    make_strategy(W, "test56_delta_logk_residual_on_M2", sws_resid, W.M2_q, 2*ones(height(W),1), ...
        delta, true(height(W),1), "residual", W.M2_risk_high20, risk_for_candidate(MODEL.risk,W,sws_resid,W.M2_q,2,CFG), false(height(W),1));
    make_strategy(W, "M2_to_M3_hard_pure_switch", choose(W.M2_sws,W.M3_sws,rule_mask), ...
        choose(W.M2_q,W.M3_q,rule_mask), choose(2*ones(height(W),1),3*ones(height(W),1),rule_mask), ...
        zeros(height(W),1), rule_mask, "rule_M3", W.M2_risk_high20, ...
        risk_for_candidate(MODEL.risk,W,choose(W.M2_sws,W.M3_sws,rule_mask),choose(W.M2_q,W.M3_q,rule_mask),choose(2*ones(height(W),1),3*ones(height(W),1),rule_mask),CFG), false(height(W),1));
    make_strategy(W, "risk_gated_Meff_residual", choose(W.M2_sws,sws_resid,resid_mask), ...
        W.M2_q, 2*ones(height(W),1), delta, resid_mask, "risk_gated_residual", W.M2_risk_high20, ...
        risk_for_candidate(MODEL.risk,W,choose(W.M2_sws,sws_resid,resid_mask),W.M2_q,2,CFG), false(height(W),1));
    apply_learned_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG);
    apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, false);
    apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, true);
    make_oracle(W, sws_resid, delta)
    ];
end

function R = apply_learned_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG)
label = string(predict(MODEL.selector, W(:, cellstr(CFG.SelectorFeatureVars))));
use_m3 = label == "switch_M3" & m3_stable & W.M2_risk_high20 < CFG.RiskThresholdSelector;
use_resid = label == "apply_residual" & resid_stable & W.M2_risk_high20 < CFG.RiskThresholdSelector;
use_resid(use_m3) = false;

sws = W.M2_sws;
q = W.M2_q;
Msrc = 2*ones(height(W),1);
mechanism = repmat("none", height(W), 1);
sws(use_m3) = W.M3_sws(use_m3);
q(use_m3) = W.M3_q(use_m3);
Msrc(use_m3) = 3;
mechanism(use_m3) = "learned_M3";
sws(use_resid) = sws_resid(use_resid);
mechanism(use_resid) = "learned_residual";
was = use_m3 | use_resid;
risk_after = risk_for_candidate(MODEL.risk,W,sws,q,Msrc,CFG);
R = make_strategy(W, "learned_candidate_selector", sws, q, Msrc, delta, was, ...
    mechanism, W.M2_risk_high20, risk_after, false(height(W),1));
end

function R = apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, spatial_gate)
label = string(predict(MODEL.interface_penalized_selector, W(:, cellstr(CFG.SelectorFeatureVars))));
hard_like = operational_hard_pure_like(W, CFG);
interface_like = operational_interface_like(W, CFG);
risk_ok = W.M2_risk_high20 < CFG.PenalizedSelectorRiskThreshold | hard_like;
not_extreme_mixed = W.M2_p_mixed < 0.75 | hard_like;

action = zeros(height(W), 1); % 0 keep M2, 1 switch M3, 2 residual
use_m3 = label == "switch_M3" & m3_stable & risk_ok & not_extreme_mixed;
use_resid = label == "apply_residual" & resid_stable & risk_ok & not_extreme_mixed;
% Be extra conservative in operationally interface-like pixels unless the
% model predicts a stable hard/pure correction.
use_m3(interface_like & ~hard_like) = false;
use_resid(interface_like & ~hard_like & W.M2_risk_high20 >= CFG.InterfaceRiskThreshold) = false;
use_resid(use_m3) = false;
action(use_m3) = 1;
action(use_resid) = 2;

if spatial_gate
    action = spatially_regularize_actions(W, action, CFG);
    strategy_name = "spatial_interface_penalized_selector";
else
    strategy_name = "interface_penalized_selector";
end

use_m3 = action == 1;
use_resid = action == 2;
sws = W.M2_sws;
q = W.M2_q;
Msrc = 2*ones(height(W),1);
mechanism = repmat("none", height(W), 1);
sws(use_m3) = W.M3_sws(use_m3);
q(use_m3) = W.M3_q(use_m3);
Msrc(use_m3) = 3;
mechanism(use_m3) = "penalized_M3";
sws(use_resid) = sws_resid(use_resid);
mechanism(use_resid) = "penalized_residual";
if spatial_gate
    mechanism(use_m3) = "spatial_penalized_M3";
    mechanism(use_resid) = "spatial_penalized_residual";
end
was = use_m3 | use_resid;
risk_after = risk_for_candidate(MODEL.risk,W,sws,q,Msrc,CFG);
R = make_strategy(W, strategy_name, sws, q, Msrc, delta, was, ...
    mechanism, W.M2_risk_high20, risk_after, false(height(W),1));
end

function R = make_oracle(W, sws_resid, delta)
err2 = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
err3 = abs(100*(W.M3_sws - W.true_SWS)./W.true_SWS);
errr = abs(100*(sws_resid - W.true_SWS)./W.true_SWS);
[~, choice] = min([err2 err3 errr], [], 2);
sws = W.M2_sws;
q = W.M2_q;
Msrc = 2*ones(height(W),1);
mech = repmat("oracle_M2", height(W), 1);
idx = choice == 2;
sws(idx) = W.M3_sws(idx); q(idx) = W.M3_q(idx); Msrc(idx) = 3; mech(idx) = "oracle_M3";
idx = choice == 3;
sws(idx) = sws_resid(idx); mech(idx) = "oracle_residual";
was = choice ~= 1;
R = make_strategy(W, "oracle_best_of_M2_M3_residual", sws, q, Msrc, delta, was, mech, ...
    W.M2_risk_high20, NaN(height(W),1), true(height(W),1));
end

function R = make_strategy(W, name, sws, q, source_M, delta, was_corrected, mechanism, risk_before, risk_after, diagnostic_only)
keep = ["dataset","base_condition","condition_key","align_key","case_id","case_family", ...
    "seen_status","field_regime","field_regime_ood","f0","dx","dz","map_ix","map_iz", ...
    "x_center_m","z_center_m","x_um","z_um","true_SWS","k_true","patch_purity", ...
    "purity_bin","distance_to_boundary_mm","distance_bin","material_region","roi_region", ...
    "M2_sws","M3_sws","M2_q","M3_q","M2_predicted_patch_purity","M2_p_mixed", ...
    "M2_p_strong_mixed","M2_M_eff","M3_M_eff"];
keep = intersect(keep, string(W.Properties.VariableNames), 'stable');
R = W(:, cellstr(keep));
R.strategy_name = repmat(string(name), height(W), 1);
R.sws_pred_strategy = sws;
R.q_pred_strategy = q;
R.source_M = source_M;
R.delta_logk_pred = delta;
R.was_corrected = was_corrected;
if isscalar(mechanism), mechanism = repmat(string(mechanism), height(W), 1); end
R.correction_mechanism = string(mechanism);
R.risk_before = risk_before;
R.risk_after = risk_after;
R.reliability_after = 1 - risk_after;
R.diagnostic_only = diagnostic_only;
R.correction_pct = 100*(sws - W.M2_sws)./W.M2_sws;
R.sws_signed_error_pct = 100*(sws - W.true_SWS)./W.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.base_abs_error_pct = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
R.improved_vs_M2 = R.sws_abs_error_pct + 0.5 < R.base_abs_error_pct;
R.harmed_vs_M2 = R.sws_abs_error_pct > R.base_abs_error_pct + 0.5;
end

function tf = operational_hard_pure_like(W, CFG)
tf = W.M2_predicted_patch_purity >= CFG.PurityThreshold & ...
    W.M2_p_mixed <= CFG.MixedThreshold & ...
    W.M2_sws >= CFG.HardSwsRatio * CFG.CsGuess;
end

function tf = operational_interface_like(W, CFG)
tf = W.M2_risk_high20 >= CFG.InterfaceRiskThreshold | ...
    W.M2_p_mixed >= CFG.InterfaceMixedThreshold | ...
    W.M2_predicted_patch_purity <= CFG.InterfacePurityThreshold;
end

function action2 = spatially_regularize_actions(W, action, CFG)
action2 = action;
r = max(0, round(CFG.SpatialGateRadiusPx));
if r == 0 || CFG.SpatialGateMinSupport <= 1
    return;
end
kernel = ones(2*r+1);
conds = unique(string(W.base_condition), 'stable');
for ci = 1:numel(conds)
    idx = find(string(W.base_condition) == conds(ci));
    if isempty(idx), continue; end
    xs = unique(W.x_um(idx));
    zs = unique(W.z_um(idx));
    [~, ix] = ismember(W.x_um(idx), xs);
    [~, iz] = ismember(W.z_um(idx), zs);
    A = zeros(numel(zs), numel(xs));
    lin = sub2ind(size(A), iz, ix);
    A(lin) = action(idx);
    for code = 1:2
        support = conv2(double(A == code), kernel, 'same');
        weak = A == code & support < CFG.SpatialGateMinSupport;
        A(weak) = 0;
    end
    action2(idx) = A(lin);
end
end

%% Risk feature utilities

function risk = predict_risk(model, W, prefix, CFG)
if string(prefix) == "M2"
    X = W(:, cellstr(CFG.RiskFeatureVars));
else
    sws = W.M3_sws;
    q = W.M3_q;
    Msrc = 3*ones(height(W),1);
    X = candidate_risk_features(W, sws, q, Msrc, CFG);
end
[~, score] = predict(model, X);
risk = positive_score(model, score);
end

function risk = risk_for_candidate(model, W, sws, q, Msrc, CFG)
X = candidate_risk_features(W, sws, q, Msrc, CFG);
[~, score] = predict(model, X);
risk = positive_score(model, score);
end

function X = candidate_risk_features(W, sws, q, Msrc, CFG)
if isscalar(Msrc), Msrc = repmat(Msrc, height(W), 1); end
sws = clamp(sws, 0.2, 20);
q = clamp(q, 0, 1);
k = 2*pi*W.f0 ./ sws;
lambda_m = sws ./ W.f0;
window_length_m = Msrc .* CFG.CsGuess ./ W.f0;
M_eff = window_length_m ./ lambda_m;
X = table(q, sws, log(sws), k, log(k), M_eff, ...
    W.M2_predicted_patch_purity, W.M2_p_mixed, W.M2_p_strong_mixed, ...
    Msrc, W.f0, W.dx, W.dz, W.field_code, W.case_family_code, ...
    'VariableNames', cellstr(CFG.RiskFeatureVars));
end

function s = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
s = score(:,idx);
end

function sws = apply_delta_logk(k_base, f0, delta)
k_corr = k_base .* exp(delta);
sws = 2*pi*f0 ./ k_corr;
sws = clamp(sws, 0.2, 20);
end

function y = choose(a, b, mask)
y = a;
y(mask) = b(mask);
end

%% Summaries and plots

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T), S = table(); return; end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), 100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.was_corrected,'omitnan'), mean(abs(X.correction_pct),'omitnan'), ...
        100*mean(X.improved_vs_M2 & X.was_corrected,'omitnan')/max(mean(X.was_corrected,'omitnan'),eps), ...
        100*mean(X.harmed_vs_M2 & X.was_corrected,'omitnan')/max(mean(X.was_corrected,'omitnan'),eps), ...
        mean(X.sws_pred_strategy,'omitnan'), mean(X.true_SWS,'omitnan'), ...
        mean(X.risk_before,'omitnan'), mean(X.risk_after,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','median_signed_error_pct','high_error10_pct', ...
        'high_error20_pct','underestimate_pct','corrected_pct', ...
        'mean_abs_correction_pct','improvement_rate_corrected_pct', ...
        'harm_rate_corrected_pct','mean_pred_sws','mean_true_sws', ...
        'mean_risk_before','mean_risk_after'});
end
S = [groups vertcat(rows{:})];
end

function S = summarize_harm(T)
S = summarize_predictions(T, ["strategy_name","correction_mechanism"]);
end

function S = success_table(T_overall, T_by_roi, CFG)
base_overall = T_overall(T_overall.strategy_name=="fixed_M2",:);
base_hard = T_by_roi(T_by_roi.strategy_name=="fixed_M2" & T_by_roi.roi_region=="hard_core_gt4mm",:);
base_soft = T_by_roi(T_by_roi.strategy_name=="fixed_M2" & T_by_roi.roi_region=="soft_core_gt8mm",:);
base_hom = T_by_roi(T_by_roi.strategy_name=="fixed_M2" & T_by_roi.roi_region=="homogeneous_center_8mm",:);
base_i0 = T_by_roi(T_by_roi.strategy_name=="fixed_M2" & T_by_roi.roi_region=="interface_0_0p5mm",:);
parts = {};
for s = unique(T_overall.strategy_name,'stable')'
    if s == "fixed_M2", continue; end
    O = T_overall(T_overall.strategy_name==s,:);
    H = T_by_roi(T_by_roi.strategy_name==s & T_by_roi.roi_region=="hard_core_gt4mm",:);
    So = T_by_roi(T_by_roi.strategy_name==s & T_by_roi.roi_region=="soft_core_gt8mm",:);
    Ho = T_by_roi(T_by_roi.strategy_name==s & T_by_roi.roi_region=="homogeneous_center_8mm",:);
    I0 = T_by_roi(T_by_roi.strategy_name==s & T_by_roi.roi_region=="interface_0_0p5mm",:);
    hard_gain = pct_gain(base_hard.MAPE_pct, H.MAPE_pct);
    soft_delta = scalar_or_nan(So.MAPE_pct) - scalar_or_nan(base_soft.MAPE_pct);
    hom_delta = scalar_or_nan(Ho.MAPE_pct) - scalar_or_nan(base_hom.MAPE_pct);
    interface_h20_delta = scalar_or_nan(I0.high_error20_pct) - scalar_or_nan(base_i0.high_error20_pct);
    global_h20_delta = scalar_or_nan(O.high_error20_pct) - scalar_or_nan(base_overall.high_error20_pct);
    pass = hard_gain >= 20 & soft_delta <= 0.5 & hom_delta <= 0.5 & ...
        interface_h20_delta <= 0.5 & global_h20_delta <= 0;
    parts{end+1,1} = table(s, hard_gain, soft_delta, hom_delta, interface_h20_delta, ...
        global_h20_delta, pass, ...
        'VariableNames', {'strategy_name','hard_core_mape_gain_pct', ...
        'soft_core_mape_delta_points','homogeneous_mape_delta_points', ...
        'interface_0_0p5_high20_delta_points','global_high20_delta_points', ...
        'passes_success_criteria'}); %#ok<AGROW>
end
S = vertcat(parts{:});
end

function plot_summary_figures(T_overall, T_case_roi, T_freq, T_roi, T_harm, OUT, CFG)
fig = figure('Color','w','Units','centimeters','Position',[1 1 36 22]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl);
X = T_overall(T_overall.strategy_name ~= "oracle_best_of_M2_M3_residual",:);
X = sortrows(X, 'MAPE_pct', 'ascend');
barh_labels(ax, pretty_strategy(X.strategy_name), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Strategy ranking','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl);
X = T_roi(T_roi.roi_region=="hard_core_gt4mm",:);
X = sortrows(X,'MAPE_pct','ascend');
barh_labels(ax, pretty_strategy(X.strategy_name), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Hard core error','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
main = ["fixed_M2","fixed_M3","test56_delta_logk_residual_on_M2", ...
    "M2_to_M3_hard_pure_switch","risk_gated_Meff_residual", ...
    "learned_candidate_selector","interface_penalized_selector", ...
    "spatial_interface_penalized_selector"];
for s = main
    F = sortrows(T_freq(T_freq.strategy_name==s,:), 'f0');
    if ~isempty(F), plot(ax,F.f0,F.MAPE_pct,'-o','DisplayName',pretty_strategy(s)); end
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'Frequency dependence','FontWeight','normal');
legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');

ax = nexttile(tl);
X = T_roi(contains(string(T_roi.roi_region),"interface"),:);
X = X(ismember(X.strategy_name, main),:);
X = sortrows(X,'MAPE_pct','descend');
barh_labels(ax, pretty_strategy(X.strategy_name)+" | "+pretty_region(X.roi_region), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Interface bands','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl);
X = T_overall(ismember(T_overall.strategy_name, main) & T_overall.corrected_pct > 0,:);
barh_labels(ax, pretty_strategy(X.strategy_name), X.harm_rate_corrected_pct);
xlabel(ax,'Harmed among corrected (%)'); title(ax,'Correction harm','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
X = T_overall(ismember(T_overall.strategy_name, main),:);
x = 1:height(X);
plot(ax, x, X.mean_risk_before, 'o-', 'DisplayName','Before');
plot(ax, x, X.mean_risk_after, 'o-', 'DisplayName','After');
set(ax, 'XTick', x, 'XTickLabel', cellstr(pretty_strategy(X.strategy_name)));
ylabel(ax,'High-error risk score'); title(ax,'Reliability before/after','FontWeight','normal');
legend(ax,'Location','best'); grid(ax,'on');
xtickangle(ax,35);

export_fig(fig, fullfile(OUT.figure_dir,'test60_hard_aware_summary.png'));

fig2 = figure('Color','w','Units','centimeters','Position',[1 1 34 16]);
tl2 = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl2);
X = T_case_roi(T_case_roi.roi_region=="hard_core_gt4mm" & ismember(T_case_roi.strategy_name, main),:);
plot_grouped_case_metric(ax, X, 'MAPE_pct', 'Hard core MAPE (%)');
ax = nexttile(tl2);
plot_grouped_case_metric(ax, X, 'mean_signed_error_pct', 'Hard core signed error (%)');
export_fig(fig2, fullfile(OUT.figure_dir,'test60_hard_core_by_case.png'));
end

function plot_all_condition_maps(T, OUT, CFG)
conds = unique(string(T.base_condition), 'stable');
if CFG.MaxMapConditions > 0 && numel(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions);
end
fprintf('Saving %d Test60 condition map panels under %s.\n', numel(conds), OUT.map_dir);
for i = 1:numel(conds)
    X = T(string(T.base_condition)==conds(i), :);
    plot_condition_map_panel(X, conds(i), OUT);
    if mod(i,25)==0 || i==numel(conds)
        fprintf('  saved %d/%d condition maps.\n', i, numel(conds));
    end
end
end

function plot_condition_map_panel(Tc, cond, OUT)
strategies = ["fixed_M2","fixed_M3","learned_candidate_selector", ...
    "interface_penalized_selector","spatial_interface_penalized_selector"];
strategy_vec = string(Tc.strategy_name);
base = Tc(strategy_vec=="fixed_M2",:);
if isempty(base), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 28]);
tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), base, base.true_SWS, 'True SWS', 'SWS (m/s)');
for s = strategies
    X = Tc(strategy_vec==s,:);
    plot_map(nexttile(tl), X, X.sws_pred_strategy, pretty_strategy(s), 'SWS (m/s)');
end
for s = ["fixed_M2","learned_candidate_selector","interface_penalized_selector", ...
        "spatial_interface_penalized_selector"]
    X = Tc(strategy_vec==s,:);
    plot_map(nexttile(tl), X, X.sws_signed_error_pct, pretty_strategy(s)+" signed error", 'Error (%)');
end
X = Tc(strategy_vec=="spatial_interface_penalized_selector",:);
plot_map(nexttile(tl), X, X.risk_before, 'Risk before', 'Risk score');
plot_map(nexttile(tl), X, X.risk_after, 'Risk after', 'Risk score');
plot_map(nexttile(tl), X, double(X.was_corrected), 'Correction mask', '0/1');
plot_map(nexttile(tl), X, X.source_M, 'Selected M', 'M');
title(tl, "Test60: " + cond, 'Interpreter','none');
out = fullfile(OUT.map_dir, sanitize_filename(string(X.case_id(1))), ...
    "test60_map__" + sanitize_filename(cond) + ".png");
export_fig(fig, out);
end

function plot_map(ax, T, values, ttl, cb_label)
[Z, xs, zs] = grid_from_points(T.x_um, T.z_um, values);
imagesc(ax, xs/1000, zs/1000, Z);
axis(ax,'image'); set(ax,'YDir','normal');
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
cb = colorbar(ax); ylabel(cb, cb_label);
end

function [Z, xs, zs] = grid_from_points(x_um, z_um, val)
xs = unique(x_um); zs = unique(z_um);
[~, ix] = ismember(x_um, xs);
[~, iz] = ismember(z_um, zs);
Z = NaN(numel(zs), numel(xs));
idx = sub2ind(size(Z), iz, ix);
Z(idx) = val;
end

%% Helper summaries/plots

function S = aggregate_for_plot(T, group_vars)
S = summarize_predictions(T, group_vars);
end

function plot_grouped_case_metric(ax, T, metric, ylab)
strategy_vec = string(T.strategy_name);
case_vec = string(T.case_id);
strategies = unique(strategy_vec,'stable');
cases = unique(case_vec,'stable');
Y = NaN(numel(cases), numel(strategies));
for i = 1:numel(cases)
    for j = 1:numel(strategies)
        idx = case_vec==cases(i) & strategy_vec==strategies(j);
        if any(idx), Y(i,j) = T.(metric)(find(idx,1)); end
    end
end
bar(ax, 1:numel(cases), Y);
set(ax, 'XTick', 1:numel(cases), 'XTickLabel', cellstr(strrep(string(cases),"_"," ")));
ylabel(ax, ylab); title(ax, ylab, 'FontWeight','normal'); grid(ax,'on');
legend(ax, pretty_strategy(strategies), 'Location','bestoutside','Interpreter','none');
xtickangle(ax,25);
end

function barh_labels(ax, labels, values)
y = 1:numel(values);
barh(ax, y, values);
set(ax, 'YTick', y, 'YTickLabel', cellstr(string(labels)));
end

function print_summary(T_overall, T_by_roi, T_harm, T_success, OUT)
fprintf('\n================ Test 60 summary ================\n');
disp(T_overall(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct','corrected_pct','harm_rate_corrected_pct'}));
fprintf('\nHard core:\n');
H = T_by_roi(T_by_roi.roi_region=="hard_core_gt4mm",:);
disp(H(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
fprintf('\nSuccess criteria:\n');
disp(T_success);
fprintf('Outputs: %s\n', OUT.root_dir);
end

function s = pretty_strategy(x)
s = string(x);
s(s=="fixed_M2") = "Fixed M=2";
s(s=="fixed_M3") = "Fixed M=3";
s(s=="test56_delta_logk_residual_on_M2") = "M-eff residual on M=2";
s(s=="M2_to_M3_hard_pure_switch") = "Hard/pure switch M2->M3";
s(s=="risk_gated_Meff_residual") = "Risk-gated M-eff residual";
s(s=="learned_candidate_selector") = "Learned candidate selector";
s(s=="interface_penalized_selector") = "Interface-penalized selector";
s(s=="spatial_interface_penalized_selector") = "Spatial interface-penalized selector";
s(s=="oracle_best_of_M2_M3_residual") = "Oracle best candidate";
s = strrep(s,"_"," ");
end

function s = pretty_region(x)
s = string(x);
s = strrep(s,"_"," ");
s = strrep(s,"0p5","0.5");
s = strrep(s,"gt",">");
end

function tf = diagnostic_only_or_false(T)
if ismember("diagnostic_only", string(T.Properties.VariableNames))
    tf = logical(T.diagnostic_only);
else
    tf = false(height(T),1);
end
end

%% Generic helpers

function idx = downsample_idx(idx, max_rows, seed)
if max_rows > 0 && numel(idx) > max_rows
    rng(seed);
    idx = idx(randperm(numel(idx), max_rows));
end
end

function T = sample_by_base_condition(T, ncond, seed)
keys = unique(string(T.base_condition), 'stable');
rng(seed);
keys = keys(randperm(numel(keys)));
keys = keys(1:min(ncond,numel(keys)));
T = T(ismember(string(T.base_condition), keys),:);
end

function y = pct_gain(base, val)
base = scalar_or_nan(base); val = scalar_or_nan(val);
y = 100*(base - val)/max(base, eps);
end

function y = scalar_or_nan(x)
if isempty(x), y = NaN; else, y = x(1); end
end

function code = categorical_code(x)
[~,~,ic] = unique(string(x), 'stable');
code = double(ic);
end

function assert_no_forbidden_predictors(features)
low = lower(string(features));
bad_exact = ["patch_purity","true_sws","k_true","q_oracle","material_region", ...
    "roi_region","distance_to_boundary_mm","distance_bin"];
bad_patterns = ["oracle","true","error","high_error"];
hit = low(ismember(low, bad_exact));
assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
for p = bad_patterns
    hit = low(contains(low,p));
    assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
end
end

function y = clamp(x, lo, hi), y = min(max(x,lo),hi); end
function y = ternary(tf, a, b), if tf, y = a; else, y = b; end, end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
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
    loc = "";
    if ~isempty(ME.stack)
        loc = sprintf(' at %s:%d', ME.stack(1).name, ME.stack(1).line);
    end
    warning('Test60:PlotFailed', 'Plot failed (%s)%s: %s', label, loc, ME.message);
end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w'); fprintf(fid,'%s\n',txt); fclose(fid);
end

function tf = env_true(name, default)
v = lower(strtrim(string(getenv(name))));
if v == "", tf = default; else, tf = ismember(v, ["1","true","yes","on"]); end
end

function x = env_number(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x), x = default; end
end

function s = env_string(name, default)
s = string(getenv(name));
if strlength(strtrim(s)) == 0, s = string(default); end
end

function name = sanitize_filename(s)
name = regexprep(char(s), '[^A-Za-z0-9_.-]+', '_');
end
