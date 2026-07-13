%% analyze_test_40_kwave_latest_model_comparison.m
% Test 40: compare Test 34 k-Wave transfer strategies with latest clean q models.
%
% This script does not train anything. It reuses:
%   - Test 34 k-Wave strategy predictions (old correction stack)
%   - Test 35 k-Wave feature table
%   - registered Test 35 and Test 38 frozen model bundles
%
% The goal is to compare old Local/Hybrid/correction strategies against the
% newer clean spectral q models on the same k-Wave inclusion cases.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.08);

%% Configuration

CFG = struct();
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST40_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST40_SAVE_ALL_MAPS', true);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST40_MAX_MAP_CONDITIONS', Inf);
CFG.Test38Source = env_string('ADAPTIVE_REQ_TEST40_TEST38_SOURCE', "full");
CFG.IncludeTest35 = env_true('ADAPTIVE_REQ_TEST40_INCLUDE_TEST35', true);
CFG.Test35ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST40_TEST35_MODELS', ...
    ["q_spectrum_only","delta_logk_theory_composition"]);
CFG.Test38ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST40_TEST38_MODELS', ...
    ["q_spectrum_only","q_spectrum_plus_composition","q_spectrum_plus_theory_composition"]);
CFG.OutputRoot = env_string('ADAPTIVE_REQ_TEST40_OUTPUT_ROOT', ...
    fullfile(root_dir, 'outputs', 'test_40_kwave_latest_model_comparison', ...
    "test38_" + sanitize(CFG.Test38Source)));
CFG.KwaveFeatureSource = "Test35 cached k-Wave REQ feature table";
CFG.OldStrategies = ["local_baseline","hybrid_baseline","theory_baseline", ...
    "test30_theory_region_levels","hybrid_lowconf_region_else_sws_nearest", ...
    "mixedness_logk_corrected","mixedness_q_candidate_selector"];
CFG.NewBundleIds = ["test35_full","test38_" + sanitize(CFG.Test38Source)];
CFG.PlotStrategies = ["T34_local_baseline","T34_hybrid_baseline", ...
    "T34_mixedness_logk_corrected","T34_mixedness_q_candidate_selector", ...
    "T35_q_spectrum_only","T35_delta_logk_theory_composition", ...
    "T38_q_spectrum_only","T38_q_spectrum_plus_composition", ...
    "T38_q_spectrum_plus_theory_composition"];
CFG.CoreDistanceMm = 4;
CFG.InterfaceDistanceMm = 2;

OUT = make_output_dirs(CFG.OutputRoot);
write_config_json(CFG, fullfile(OUT.root_dir, 'test40_configuration.json'));

fprintf('\nTest 40: k-Wave comparison of old Test34 stack vs latest clean q models\n');
fprintf('No training. Output: %s\n', OUT.root_dir);
fprintf('Test38 bundle source: %s\n', CFG.Test38Source);
fprintf('k-Wave features: %s; REQ is not recalculated in this script.\n', CFG.KwaveFeatureSource);

SRC = locate_sources(root_dir, CFG);
assert(exist(SRC.test34_csv,'file') == 2, 'Missing Test34 patch CSV: %s', SRC.test34_csv);
assert(exist(SRC.test35_kwave_features,'file') == 2, 'Missing Test35 k-Wave feature dataset: %s', SRC.test35_kwave_features);
if CFG.IncludeTest35
    assert(exist(SRC.test35_bundle,'file') == 2, 'Missing registered Test35 bundle: %s', SRC.test35_bundle);
end
assert(exist(SRC.test38_bundle,'file') == 2, 'Missing registered Test38 bundle: %s', SRC.test38_bundle);

if CFG.ValidateOnly
    validate_inputs(SRC);
    fprintf('Test 40 validation-only checks passed.\n');
    return;
end

%% Load and evaluate

T_old_wide = readtable(SRC.test34_csv);
T_old = old_test34_to_long(T_old_wide, CFG);

S = load(SRC.test35_kwave_features, 'T_kw');
T_kw = standardize_kwave_feature_table(S.T_kw);
clear S;

fprintf('Loaded Test34 old rows: %d wide patches -> %d strategy rows.\n', ...
    height(T_old_wide), height(T_old));
fprintf('Loaded Test35 k-Wave feature table: %d patches.\n', height(T_kw));

new_parts = {};
if CFG.IncludeTest35
    fprintf('Applying Test35 bundle: %s\n', SRC.test35_bundle);
    T_new35 = apply_registered_bundle(T_kw, SRC.test35_bundle, "T35", "test35_full", ...
        CFG.Test35ModelsToEvaluate);
    new_parts{end+1} = T_new35; %#ok<SAGROW>
    fprintf('Finished Test35 bundle: %d prediction rows.\n', height(T_new35));
else
    fprintf('Skipping Test35 bundle by ADAPTIVE_REQ_TEST40_INCLUDE_TEST35=false.\n');
end

fprintf('Applying Test38 bundle: %s\n', SRC.test38_bundle);
T_new38 = apply_registered_bundle(T_kw, SRC.test38_bundle, "T38", ...
    "test38_" + sanitize(CFG.Test38Source), CFG.Test38ModelsToEvaluate);
new_parts{end+1} = T_new38;
fprintf('Finished Test38 bundle: %d prediction rows.\n', height(T_new38));

T_all = vertcat_compatible(T_old, new_parts{:});
T_all = add_analysis_regions(T_all, CFG);

%% Summaries

T_overall = summarize_long(T_all, ["model_family","strategy_name"]);
T_by_family = summarize_long(T_all, ["model_family","strategy_name"]);
T_by_case = summarize_long(T_all, ["strategy_name","case_name","field_regime"]);
T_by_M = summarize_long(T_all, ["strategy_name","M"]);
T_by_side = summarize_long(T_all, ["strategy_name","material_side"]);
T_by_roi = summarize_long(T_all, ["strategy_name","roi_name"]);
T_by_distance = summarize_long(T_all, ["strategy_name","distance_bin"]);
T_by_region = summarize_long(T_all, ["strategy_name","analysis_region"]);
T_by_case_M = summarize_long(T_all, ["strategy_name","case_name","field_regime","M"]);

writetable(T_all, fullfile(OUT.table_dir, 'test40_kwave_patch_level_model_comparison.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test40_summary_overall.csv'));
writetable(T_by_family, fullfile(OUT.table_dir, 'test40_summary_by_model_family.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test40_summary_by_case.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test40_summary_by_M.csv'));
writetable(T_by_side, fullfile(OUT.table_dir, 'test40_summary_by_soft_hard.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir, 'test40_summary_by_roi.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test40_summary_by_distance.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test40_summary_by_analysis_region.csv'));
writetable(T_by_case_M, fullfile(OUT.table_dir, 'test40_summary_by_case_M.csv'));
save(fullfile(OUT.data_dir, 'test40_kwave_latest_model_comparison.mat'), ...
    'T_all','T_overall','T_by_family','T_by_case','T_by_M','T_by_side', ...
    'T_by_roi','T_by_distance','T_by_region','T_by_case_M','CFG','SRC','-v7.3');

%% Figures

plot_overall_ranking(T_overall, OUT);
plot_group_heatmap(T_by_case, "case_name", OUT, 'test40_mape_by_case.png');
plot_group_heatmap(T_by_M, "M", OUT, 'test40_mape_by_M.png');
plot_soft_hard_roi(T_by_side, T_by_roi, OUT);
plot_distance_curves(T_by_distance, OUT);
plot_family_comparison(T_by_family, OUT);
if CFG.SaveAllMaps
    plot_condition_maps(T_all, CFG, OUT);
end

print_interpretation(T_overall, T_by_case, T_by_side, T_by_roi, T_by_distance);

fprintf('\nTables: %s\nFigures: %s\nTest 40 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Loading and standardization

function OUT = make_output_dirs(root_dir)
OUT.root_dir = root_dir;
OUT.table_dir = fullfile(root_dir, 'tables');
OUT.figure_dir = fullfile(root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(root_dir, 'data');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function SRC = locate_sources(root_dir, CFG)
SRC = struct();
SRC.test34_csv = fullfile(root_dir, 'outputs', 'test_34_kwave_mixedness_transfer', ...
    'tables', 'test34_kwave_patch_level_predictions.csv');
SRC.test35_kwave_features = fullfile(root_dir, 'outputs', ...
    'test_35_spectral_composition_to_q_model', 'data', 'test35_kwave_feature_dataset.mat');
SRC.test35_bundle = fullfile(root_dir, 'outputs', 'model_registry', ...
    'test35_spectral_composition_to_q_model', 'test35__spectral_composition_q__full_bundle.mat');
SRC.test38_bundle = resolve_test38_bundle(root_dir, CFG.Test38Source);
end

function path = resolve_test38_bundle(root_dir, source)
source = string(source);
if exist(source, 'file') == 2
    path = char(source);
    return;
end
switch lower(source)
    case "full"
        path = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
            'models', 'test38_velocity_field_diverse_q_models.mat');
    case "medium"
        path = fullfile(root_dir, 'outputs', 'model_registry', ...
            'test38_velocity_field_diverse_q_training', ...
            'test38__velocity_field_diverse_q__medium_bundle.mat');
    case "quick"
        path = fullfile(root_dir, 'outputs', 'model_registry', ...
            'test38_velocity_field_diverse_q_training', ...
            'test38__velocity_field_diverse_q__quick_bundle.mat');
    otherwise
        error('Unknown ADAPTIVE_REQ_TEST40_TEST38_SOURCE: %s. Use full, medium, quick, or a .mat path.', source);
end
end

function validate_inputs(SRC)
T = readtable(SRC.test34_csv);
assert(ismember("sws_local_baseline", string(T.Properties.VariableNames)));
assert(ismember("abs_error_mixedness_logk_corrected", string(T.Properties.VariableNames)));
S = load(SRC.test35_kwave_features, 'T_kw');
assert(istable(S.T_kw) && height(S.T_kw) > 0);
needed = ["req_mapping","q_oracle","q_theory_prior","true_SWS","case_name","M"];
assert(all(ismember(needed, string(S.T_kw.Properties.VariableNames))), ...
    'T_kw missing one or more required variables.');
B35 = load(SRC.test35_bundle, 'MODEL_BUNDLE');
B38 = load(SRC.test38_bundle, 'MODEL_BUNDLE');
assert(isfield(B35.MODEL_BUNDLE,'MODELS') && isfield(B38.MODEL_BUNDLE,'MODELS'));
end

function T = old_test34_to_long(W, CFG)
parts = cell(numel(CFG.OldStrategies),1);
base_vars = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "M","f0","dx","dz","map_iz","map_ix","x","z","true_SWS","cs_true", ...
    "roi_name","material_side","region_label","distance_to_interface_mm", ...
    "distance_bin","purity_bin","confidence","predicted_mixed_probability", ...
    "predicted_patch_purity"];
base_vars = base_vars(ismember(base_vars, string(W.Properties.VariableNames)));
for i = 1:numel(CFG.OldStrategies)
    s = CFG.OldStrategies(i);
    R = W(:, cellstr(base_vars));
    R.model_family = repmat("T34_old_stack", height(W), 1);
    R.bundle_id = repmat("test34", height(W), 1);
    R.strategy_name = repmat("T34_" + s, height(W), 1);
    R.sws_pred = W.("sws_" + s);
    R.sws_signed_error_pct = W.("signed_error_" + s);
    R.sws_abs_error_pct = W.("abs_error_" + s);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    R.q_pred = nan(height(W),1);
    R.q_error = nan(height(W),1);
    if ~ismember("patch_purity", string(R.Properties.VariableNames))
        R.patch_purity = nan(height(W),1);
    end
    if ~ismember("p_mixed", string(R.Properties.VariableNames))
        if ismember("predicted_mixed_probability", string(R.Properties.VariableNames))
            R.p_mixed = R.predicted_mixed_probability;
        else
            R.p_mixed = nan(height(W),1);
        end
    end
    parts{i} = R;
end
T = vertcat(parts{:});
end

function T = standardize_kwave_feature_table(T)
if ~ismember("true_SWS", string(T.Properties.VariableNames))
    T.true_SWS = T.cs_true;
end
if ~ismember("dx", string(T.Properties.VariableNames)), T.dx = 0.5e-3*ones(height(T),1); end
if ~ismember("dz", string(T.Properties.VariableNames)), T.dz = T.dx; end
if ~ismember("f0", string(T.Properties.VariableNames)), T.f0 = T.frequency_hz; end
if ~ismember("M", string(T.Properties.VariableNames)), T.M = T.REQ_M; end
if ~ismember("material_side", string(T.Properties.VariableNames))
    T.material_side = repmat("soft", height(T), 1);
    T.material_side(T.true_SWS > 2.5) = "hard";
end
if ~ismember("field_regime", string(T.Properties.VariableNames))
    T.field_regime = arrayfun(@canonical_kwave_regime, string(T.case_name));
end
if ~ismember("geometry", string(T.Properties.VariableNames))
    T.geometry = repmat("kwave_circular_inclusion_2_3", height(T), 1);
end
if ~ismember("geometry_type", string(T.Properties.VariableNames))
    T.geometry_type = repmat("inclusion", height(T), 1);
end
if ~ismember("distance_to_interface_mm", string(T.Properties.VariableNames))
    T.distance_to_interface_mm = abs(hypot(T.x_center_m - 0.025, T.z_center_m - 0.025) - 0.010) * 1e3;
end
if ~ismember("distance_bin", string(T.Properties.VariableNames))
    T.distance_bin = distance_bin(T.distance_to_interface_mm);
end
if ~ismember("roi_name", string(T.Properties.VariableNames))
    T.roi_name = repmat("unknown", height(T), 1);
end
if ~ismember("condition_key", string(T.Properties.VariableNames))
    T.condition_key = "kwave__" + sanitize(T.case_name) + "__M" + string(T.M);
end
end

%% Frozen model application

function T_out = apply_registered_bundle(T_kw, bundle_file, prefix, bundle_id, models_to_evaluate)
if nargin < 5 || isempty(models_to_evaluate)
    models_to_evaluate = strings(0,1);
else
    models_to_evaluate = string(models_to_evaluate(:));
end
fprintf('  Loading bundle %s...\n', bundle_file);
B = load(bundle_file, 'MODEL_BUNDLE');
MODEL_BUNDLE = B.MODEL_BUNDLE;
assert_no_forbidden_predictors(MODEL_BUNDLE.BASE_FEATURES);
F = ensure_bundle_predictor_columns(T_kw, MODEL_BUNDLE.BASE_FEATURES);
fprintf('  Predicting composition features...\n');
P = predict_composition(MODEL_BUNDLE.MODELS.composition, F, MODEL_BUNDLE.BASE_FEATURES);
F.predicted_patch_purity = P.predicted_patch_purity;
F.p_mixed = P.p_mixed;
F.p_strong_mixed = P.p_strong_mixed;

parts = {};
model_names = string(MODEL_BUNDLE.MODELS.model_names(:));
if ~isempty(models_to_evaluate)
    model_names = model_names(ismember(model_names, models_to_evaluate));
end
for name = model_names(:)'
    if name == "theory_discrete"
        % Keep Test34 theory baseline as the shared old reference.
        continue;
    end
    fprintf('  Predicting %s_%s...\n', prefix, name);
    R = base_new_prediction_table(F, prefix, bundle_id, name);
    switch name
        case "q_spectrum_only"
            q_pred = predict_q_model(MODEL_BUNDLE.MODELS.q.spectrum_only, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_theory"
            q_pred = predict_q_model(MODEL_BUNDLE.MODELS.q.spectrum_plus_theory, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODEL_BUNDLE.MODELS.q.spectrum_plus_composition, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODEL_BUNDLE.MODELS.q.spectrum_plus_theory_composition, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "delta_q_theory_composition"
            delta = predict_q_model(MODEL_BUNDLE.MODELS.q.delta_q_theory_composition, F, false);
            q_pred = clamp01(F.q_theory_prior + delta);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "delta_logk_theory_composition"
            delta = predict_q_model(MODEL_BUNDLE.MODELS.q.delta_logk_theory_composition, F, false);
            k_theory = 2*pi*F.f0 ./ F.sws_theory;
            k_pred = k_theory .* exp(delta);
            sws_pred = 2*pi*F.f0 ./ k_pred;
            q_pred = arrayfun(@(i) invert_mapping_to_q(F.req_mapping{i}, k_pred(i)), (1:height(F))');
        otherwise
            warning('Skipping unknown model %s in %s.', name, bundle_id);
            continue;
    end
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - F.q_oracle;
    R.sws_signed_error_pct = 100*(sws_pred - F.true_SWS) ./ F.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{end+1,1} = R; %#ok<AGROW>
end
T_out = vertcat(parts{:});
end

function R = base_new_prediction_table(F, prefix, bundle_id, model_name)
keep = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "M","f0","dx","dz","map_iz","map_ix","true_SWS","cs_true","roi_name", ...
    "material_side","distance_to_interface_mm","distance_bin","purity_bin", ...
    "patch_purity","predicted_patch_purity","p_mixed"];
keep = keep(ismember(keep, string(F.Properties.VariableNames)));
R = F(:, cellstr(keep));
if ismember("x_center_m", string(F.Properties.VariableNames))
    R.x = F.x_center_m; R.z = F.z_center_m;
elseif ismember("x", string(F.Properties.VariableNames))
    R.x = F.x; R.z = F.z;
else
    R.x = nan(height(F),1); R.z = nan(height(F),1);
end
if ~ismember("cs_true", string(R.Properties.VariableNames)), R.cs_true = R.true_SWS; end
if ~ismember("purity_bin", string(R.Properties.VariableNames)), R.purity_bin = repmat("unknown",height(F),1); end
if ~ismember("predicted_patch_purity", string(R.Properties.VariableNames)), R.predicted_patch_purity = nan(height(F),1); end
if ~ismember("confidence", string(R.Properties.VariableNames)), R.confidence = nan(height(F),1); end
if ~ismember("predicted_mixed_probability", string(R.Properties.VariableNames))
    if ismember("p_mixed", string(R.Properties.VariableNames))
        R.predicted_mixed_probability = R.p_mixed;
    else
        R.predicted_mixed_probability = nan(height(F),1);
    end
end
R.model_family = repmat(prefix + "_latest_clean_q", height(F), 1);
R.bundle_id = repmat(bundle_id, height(F), 1);
R.strategy_name = repmat(prefix + "_" + model_name, height(F), 1);
end

function F = ensure_bundle_predictor_columns(F, features)
for f = string(features(:))'
    if ~ismember(f, string(F.Properties.VariableNames))
        switch f
            case "dx", F.dx = F.dx;
            case "dz", F.dz = F.dz;
            case "f0", F.f0 = F.f0;
            case "M", F.M = F.M;
            case "REQ_M", F.REQ_M = F.M;
            case "SIM_f0", F.SIM_f0 = F.f0;
            case "REQ_StepX"
                F.REQ_StepX = infer_req_step(F, "x");
            case "REQ_StepZ"
                F.REQ_StepZ = infer_req_step(F, "z");
            case "TargetStepM"
                F.TargetStepM = infer_target_step_m(F);
            otherwise
                error('Required predictor missing for registered bundle: %s', f);
        end
    end
end

function step = infer_req_step(F, axis_name)
n = height(F);
if axis_name == "x"
    coord_name = "x_center_m"; pixel_name = "map_ix"; spacing = F.dx;
else
    coord_name = "z_center_m"; pixel_name = "map_iz"; spacing = F.dz;
end
step_value = NaN;
if ismember(coord_name, string(F.Properties.VariableNames))
    u = unique(F.(coord_name));
    u = sort(u(isfinite(u)));
    if numel(u) > 1
        d = median(diff(u), 'omitnan');
        if isfinite(d) && d > 0
            step_value = max(1, round(d / median(spacing, 'omitnan')));
        end
    end
elseif ismember(pixel_name, string(F.Properties.VariableNames))
    u = unique(F.(pixel_name));
    u = sort(u(isfinite(u)));
    if numel(u) > 1
        step_value = max(1, round(median(diff(u), 'omitnan')));
    end
end
if ~isfinite(step_value), step_value = 1; end
step = step_value * ones(n,1);
end

function target_step_m = infer_target_step_m(F)
if ~ismember("REQ_StepX", string(F.Properties.VariableNames))
    F.REQ_StepX = infer_req_step(F, "x");
end
if ~ismember("REQ_StepZ", string(F.Properties.VariableNames))
    F.REQ_StepZ = infer_req_step(F, "z");
end
target_step_m = max(F.REQ_StepX .* F.dx, F.REQ_StepZ .* F.dz);
bad = ~isfinite(target_step_m) | target_step_m <= 0;
target_step_m(bad) = median(F.dx, 'omitnan');
end
if ~ismember("sws_theory", string(F.Properties.VariableNames))
    F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
end
end

function P = predict_composition(MIX, T, base_features)
X = T(:, cellstr(base_features));
P = struct();
P.predicted_patch_purity = min(max(predict(MIX.purity, X), 0), 1);
[~,score] = predict(MIX.mixed, X);
P.p_mixed = positive_score(MIX.mixed, score);
[~,score] = predict(MIX.strong_mixed, X);
P.p_strong_mixed = positive_score(MIX.strong_mixed, score);
end

function y = predict_q_model(M, T, do_clamp)
if nargin < 3, do_clamp = true; end
X = T(:, cellstr(M.features));
y = predict(M.model, X);
if do_clamp, y = clamp01(y); end
end

%% Summaries

function T = add_analysis_regions(T, CFG)
T.analysis_region = repmat("other", height(T), 1);
T.analysis_region(T.distance_to_interface_mm <= CFG.InterfaceDistanceMm) = "interface_0_2mm";
T.analysis_region(T.distance_to_interface_mm > CFG.InterfaceDistanceMm & ...
    T.distance_to_interface_mm <= CFG.CoreDistanceMm) = "near_interface_2_4mm";
T.analysis_region(T.distance_to_interface_mm > CFG.CoreDistanceMm & T.material_side == "soft") = "soft_core_gt4mm";
T.analysis_region(T.distance_to_interface_mm > CFG.CoreDistanceMm & T.material_side == "hard") = "hard_core_gt4mm";
end

function S = summarize_long(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
n = splitapply(@numel, T.sws_abs_error_pct, G);
mape = splitapply(@(x) mean(x,'omitnan'), T.sws_abs_error_pct, G);
medae = splitapply(@(x) median(x,'omitnan'), T.sws_abs_error_pct, G);
signed = splitapply(@(x) mean(x,'omitnan'), T.sws_signed_error_pct, G);
medsigned = splitapply(@(x) median(x,'omitnan'), T.sws_signed_error_pct, G);
under = splitapply(@(x) 100*mean(x < 0,'omitnan'), T.sws_signed_error_pct, G);
he10 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error10, G);
he20 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
mean_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
std_sws = splitapply(@std_omitnan, T.sws_pred, G);
S = [groups table(n,mape,medae,signed,medsigned,under,he10,he20,mean_sws,std_sws, ...
    'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
    'median_signed_error_pct','underestimate_pct','high_error10_pct','high_error20_pct', ...
    'mean_sws_pred','std_sws_pred'})];
end

%% Plots

function plot_overall_ranking(T, OUT)
T = sortrows(T, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 13]);
bar(categorical(T.strategy_name), T.MAPE_pct);
ylabel('MAPE (%)'); grid on; xtickangle(35);
title('Test 40 k-Wave overall MAPE','Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test40_overall_model_ranking.png'));
end

function plot_group_heatmap(T, group_var, OUT, filename)
strategies = unique(T.strategy_name, 'stable');
groups = unique(string(T.(group_var)), 'stable');
Z = nan(numel(strategies), numel(groups));
for i = 1:numel(strategies)
    for j = 1:numel(groups)
        idx = T.strategy_name == strategies(i) & string(T.(group_var)) == groups(j);
        if any(idx), Z(i,j) = mean(T.MAPE_pct(idx), 'omitnan'); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 16]);
imagesc(Z); colorbar; colormap(parula);
set(gca,'XTick',1:numel(groups),'XTickLabel',groups,'YTick',1:numel(strategies), ...
    'YTickLabel',strategies,'TickLabelInterpreter','none');
xtickangle(35); title("MAPE by " + group_var, 'Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, filename));
end

function plot_soft_hard_roi(T_side, T_roi, OUT)
keep = select_plot_strategies(T_side.strategy_name);
S = T_side(ismember(T_side.strategy_name, keep),:);
R = T_roi(ismember(T_roi.strategy_name, keep),:);
side_groups = unique(string(S.material_side), 'stable');
roi_groups = unique(string(R.roi_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 13]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar(categorical(side_groups), grouped_values(S, "material_side", keep));
legend(keep,'Interpreter','none','Location','bestoutside'); ylabel('MAPE (%)'); grid on;
title('Soft/hard MAPE','FontWeight','normal');
nexttile;
bar(categorical(roi_groups), grouped_values(R, "roi_name", keep));
legend(keep,'Interpreter','none','Location','bestoutside'); ylabel('MAPE (%)'); grid on;
title('ROI MAPE','FontWeight','normal'); xtickangle(35);
export_fig(fig, fullfile(OUT.figure_dir, 'test40_soft_hard_roi_mape.png'));
end

function plot_distance_curves(T_dist, OUT)
keep = select_plot_strategies(T_dist.strategy_name);
order = ["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
fig = figure('Color','w','Units','centimeters','Position',[1 1 22 13]);
hold on;
for si = 1:numel(keep)
    s = keep(si);
    y = nan(size(order));
    for i = 1:numel(order)
        idx = T_dist.strategy_name == s & string(T_dist.distance_bin) == order(i);
        if any(idx), y(i) = mean(T_dist.MAPE_pct(idx), 'omitnan'); end
    end
    plot(1:numel(order), y, '-o', 'LineWidth', 1.0, 'DisplayName', s);
end
set(gca,'XTick',1:numel(order),'XTickLabel',order); xtickangle(25);
ylabel('MAPE (%)'); xlabel('Distance to interface'); grid on;
legend('Interpreter','none','Location','bestoutside');
title('Error vs distance to interface','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test40_mape_vs_distance.png'));
end

function plot_family_comparison(T_family, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 22 11]);
families = unique(T_family.model_family, 'stable');
y = nan(numel(families),1);
for i = 1:numel(families)
    y(i) = min(T_family.MAPE_pct(T_family.model_family == families(i)));
end
bar(categorical(families), y); ylabel('Best MAPE in family (%)'); grid on;
title('Best strategy/model per family','Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test40_best_by_model_family.png'));
end

function plot_condition_maps(T, CFG, OUT)
keep = string(CFG.PlotStrategies(:));
conds = unique(T(:, {'case_name','field_regime','M'}), 'rows', 'stable');
if isfinite(CFG.MaxMapConditions) && height(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions,:);
end
fprintf('Saving %d Test 40 map panels under %s.\n', height(conds), OUT.map_dir);
for ci = 1:height(conds)
    idxc = T.case_name == conds.case_name(ci) & T.field_regime == conds.field_regime(ci) & T.M == conds.M(ci);
    Xc = T(idxc,:);
    if isempty(Xc), continue; end
    base = Xc(Xc.strategy_name == Xc.strategy_name(1),:);
    [true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 34 24]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS');
    plot_map(nexttile(tl), rows_to_grid(base, base.distance_to_interface_mm, nz, nx), 'Distance mm');
    for si = 1:numel(keep)
        s = keep(si);
        X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), s + ' SWS');
    end
    err_keep = keep(1:min(5,numel(keep)));
    for si = 1:numel(err_keep)
        s = err_keep(si);
        X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, nz, nx), s + ' abs err %');
    end
    title(tl, sprintf('%s | %s | M%d', conds.case_name(ci), conds.field_regime(ci), conds.M(ci)), ...
        'Interpreter','none');
    outdir = fullfile(OUT.map_dir, sanitize(conds.case_name(ci)), sanitize(conds.field_regime(ci)), ...
        "M" + string(conds.M(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    export_fig(fig, fullfile(outdir, sprintf('test40_map__%s__%s__M%d.png', ...
        sanitize(conds.case_name(ci)), sanitize(conds.field_regime(ci)), conds.M(ci))));
end
end

function keep = select_plot_strategies(strategies)
preferred = ["T34_local_baseline","T34_hybrid_baseline","T34_mixedness_logk_corrected", ...
    "T34_mixedness_q_candidate_selector","T35_q_spectrum_only", ...
    "T38_q_spectrum_only","T38_q_spectrum_plus_composition", ...
    "T38_q_spectrum_plus_theory_composition"];
keep = preferred(ismember(preferred, unique(strategies)));
end

function Y = grouped_values(T, group_var, strategies)
groups = unique(string(T.(group_var)), 'stable');
Y = nan(numel(groups), numel(strategies));
for i = 1:numel(groups)
    for j = 1:numel(strategies)
        idx = string(T.(group_var)) == groups(i) & T.strategy_name == strategies(j);
        if any(idx), Y(i,j) = mean(T.MAPE_pct(idx), 'omitnan'); end
    end
end
end

function print_interpretation(T_overall, T_by_case, T_by_side, T_by_roi, T_by_distance)
T = sortrows(T_overall, 'MAPE_pct', 'ascend');
fprintf('\n================ Test 40 k-Wave comparison summary ================\n');
disp(T(:, {'strategy_name','model_family','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
fprintf('Best overall: %s (MAPE %.2f%%, HE20 %.2f%%).\n', ...
    T.strategy_name(1), T.MAPE_pct(1), T.high_error20_pct(1));
old = T(startsWith(T.strategy_name,"T34_"),:);
new = T(~startsWith(T.strategy_name,"T34_"),:);
if ~isempty(old) && ~isempty(new)
    old = sortrows(old,'MAPE_pct','ascend'); new = sortrows(new,'MAPE_pct','ascend');
    fprintf('Best old Test34 strategy: %s, MAPE %.2f%%.\n', old.strategy_name(1), old.MAPE_pct(1));
    fprintf('Best latest clean q model: %s, MAPE %.2f%%.\n', new.strategy_name(1), new.MAPE_pct(1));
end
fprintf('Tables by case/M/side/ROI/distance are saved for detailed inspection.\n');
fprintf('====================================================================\n');
end

%% Utilities

function q = clamp01(q)
q = min(max(q, 0.001), 0.999);
end

function T = vertcat_compatible(varargin)
vars = strings(1,0);
for i = 1:nargin
    vars = [vars string(varargin{i}.Properties.VariableNames)]; %#ok<AGROW>
end
vars = unique(vars, 'stable');
string_vars = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "roi_name","material_side","region_label","distance_bin","purity_bin", ...
    "model_family","bundle_id","strategy_name","analysis_region"];
logical_vars = ["high_error10","high_error20"];
parts = cell(nargin,1);
for i = 1:nargin
    A = varargin{i};
    for v = vars
        if ismember(v, string(A.Properties.VariableNames)), continue; end
        if ismember(v, string_vars)
            A.(v) = repmat("unknown", height(A), 1);
        elseif ismember(v, logical_vars)
            A.(v) = false(height(A), 1);
        else
            A.(v) = nan(height(A), 1);
        end
    end
    parts{i} = A(:, cellstr(vars));
end
T = vertcat(parts{:});
end

function y = q_to_sws(mappings, q, f0)
n = numel(q); y = nan(n,1);
for i = 1:n
    if isscalar(f0), f = f0; else, f = f0(i); end
    if isfinite(q(i))
        y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), f);
    end
end
end

function q = invert_mapping_to_q(mapping, k_target)
k = double(mapping.k_cent(:));
E = double(mapping.Ecum(:));
valid = isfinite(k) & isfinite(E);
k = k(valid); E = E(valid);
if numel(k) < 2 || ~isfinite(k_target), q = NaN; return; end
[kuniq, ia] = unique(k, 'stable');
Euniq = E(ia);
if numel(kuniq) < 2
    q = NaN; return;
elseif k_target <= min(kuniq)
    q = Euniq(1);
elseif k_target >= max(kuniq)
    q = Euniq(end);
else
    q = interp1(kuniq, Euniq, k_target, 'linear');
end
q = clamp01(q);
end

function p = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
p = score(:,idx);
end

function assert_no_forbidden_predictors(features)
bad = ["true","oracle","purity","mixed","confidence","error","pred", ...
    "sws","cs_","k_true","q_local","q_pred","req_mapping","patch_idx", ...
    "map_ix","map_iz","cx","cz","x_center","z_center","condition"];
features = lower(string(features));
for b = bad
    hit = features(contains(features,b));
    assert(isempty(hit), 'Forbidden predictor in registered model bundle: %s', strjoin(hit, ', '));
end
end

function bin = distance_bin(d)
bin = repmat("gt_8mm", numel(d), 1);
bin(d <= 8) = "4_8mm";
bin(d <= 4) = "2_4mm";
bin(d <= 2) = "1_2mm";
bin(d <= 1) = "0_1mm";
end

function label = canonical_kwave_regime(case_name)
switch string(case_name)
    case "Directional 2D", label = "directional_2D";
    case "Diffuse 2D", label = "diffuse_2D";
    case {"Diffuse 3D","Projected diffuse 3D"}, label = "diffuse_3D";
    otherwise, label = "partial_3D";
end
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
end

function s = std_omitnan(x)
s = std(x, 'omitnan');
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map(ax, Z, ttl)
imagesc(ax, Z); axis(ax,'image'); axis(ax,'off'); colorbar(ax);
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
end

function export_fig(fig, path)
drawnow;
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
close(fig);
end

function write_config_json(CFG, path)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(path,'w'); fwrite(fid, txt); fclose(fid);
end

function tf = env_true(name, default)
v = string(getenv(name));
if v == "", tf = default; return; end
tf = ismember(lower(v), ["1","true","yes","y","on"]);
end

function x = env_number(name, default)
v = string(getenv(name));
if v == "", x = default; return; end
x = str2double(v);
if isnan(x), x = default; end
end

function x = env_string(name, default)
v = string(getenv(name));
if v == ""
    x = string(default);
else
    x = v;
end
end

function xs = env_string_list(name, default)
v = string(getenv(name));
if v == ""
    xs = string(default);
    return;
end
parts = split(v, [",",";"]);
parts = strip(parts);
xs = parts(parts ~= "");
end
