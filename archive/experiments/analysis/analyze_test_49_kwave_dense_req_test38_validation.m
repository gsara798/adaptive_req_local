%% analyze_test_49_kwave_dense_req_test38_validation.m
% Test 49: dense k-Wave REQ validation for frozen Test38 q models.
%
% This script re-extracts REQ/features from the raw k-Wave harmonic fields
% using the Test38 training REQ profile:
%   Nbins='auto', Nbins_auto_oversample=1, Nbins_min=16, smooth_sigma=1,
%   gamma_win=1, pad_factor=1, EdgeMode='valid'.
%
% No q model is trained. Frozen Test38 models are loaded and applied to the
% newly extracted k-Wave REQ mappings/features. k-Wave truth is used only for
% q_oracle/error/ROI diagnostics, never as a predictor.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

%% Configuration

CFG = struct();
CFG.Mode = lower(env_string('ADAPTIVE_REQ_TEST49_MODE', "quick"));
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST49_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST49_SAVE_ALL_MAPS', true);
CFG.ModelSource = env_string('ADAPTIVE_REQ_TEST49_MODEL_SOURCE', "full");
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST49_TARGET_STEP_M', 0.5e-3);
CFG.UseWindowParfor = env_true('ADAPTIVE_REQ_TEST49_USE_PARFOR', true);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST49_MAX_MAP_CONDITIONS', Inf);
CFG.ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST49_MODELS', ...
    ["q_spectrum_only","q_spectrum_plus_composition","q_spectrum_plus_theory_composition"]);
CFG.ModelLabelPrefix = "T38full";
CFG.cs_guess = 3.0;
CFG.KW = default_kwave_config();
CFG.REQ = struct('Gamma',1,'PadFactor',1,'Nbins',"auto", ...
    'Nbins_auto_oversample',1,'Nbins_min',16,'SmoothSigma',1,'EdgeMode',"valid");

if CFG.Mode == "quick"
    CFG.KW.cases = CFG.KW.cases([1 3]); % Directional 2D and projected diffuse 3D.
    CFG.KW.REQ.M_list = [2];
elseif CFG.Mode == "medium"
    CFG.KW.cases = CFG.KW.cases([1 2 3]);
    CFG.KW.REQ.M_list = [2 3];
elseif CFG.Mode == "full"
    CFG.KW.REQ.M_list = [2 3 4];
else
    error('Unknown ADAPTIVE_REQ_TEST49_MODE=%s. Use quick, medium, or full.', CFG.Mode);
end

OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test49_configuration.json'));

fprintf('\nTest 49: dense k-Wave REQ validation for frozen Test38 models\n');
fprintf('Mode: %s | model source: %s | target step %.3f mm | parfor windows: %d\n', ...
    CFG.Mode, CFG.ModelSource, 1e3*CFG.TargetStepM, CFG.UseWindowParfor);
fprintf('REQ profile: Nbins=%s, oversample=%g, min=%g, smooth_sigma=%g, gamma=%g, pad=%g.\n', ...
    CFG.REQ.Nbins, CFG.REQ.Nbins_auto_oversample, CFG.REQ.Nbins_min, ...
    CFG.REQ.SmoothSigma, CFG.REQ.Gamma, CFG.REQ.PadFactor);
fprintf('No training. k-Wave truth is evaluation-only.\n');

BUNDLE = load_test38_bundle(root_dir, CFG);
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);
fprintf('Loaded frozen Test38 bundle with %d primitive predictors.\n', numel(BASE_FEATURES));

if CFG.ValidateOnly
    validate_sources(root_dir, CFG, BUNDLE, BASE_FEATURES);
    fprintf('Test 49 validation-only checks passed.\n');
    return;
end

%% Main dense k-Wave loop

T_all = table();
condition_count = numel(CFG.KW.cases) * numel(CFG.KW.REQ.M_list);
counter = 0;
for ci = 1:numel(CFG.KW.cases)
    case_i = CFG.KW.cases(ci);
    [field, dinf, cfg_sim] = load_kwave_field(root_dir, CFG, case_i);
    fprintf('\nLoaded k-Wave case %s | field %s | dx %.3f mm | dz %.3f mm\n', ...
        case_i.case_name, mat2str(size(field)), 1e3*cfg_sim.dx, 1e3*cfg_sim.dz);
    for mi = 1:numel(CFG.KW.REQ.M_list)
        M = CFG.KW.REQ.M_list(mi);
        counter = counter + 1;
        key = sprintf('kwave_dense__%s__M%g__step%gum', ...
            sanitize(case_i.case_name), M, round(1e6*CFG.TargetStepM));
        fprintf('[%d/%d] Extract/evaluate %s\n', counter, condition_count, key);
        t0 = tic;
        [T_pred, F] = evaluate_kwave_condition(field, dinf, cfg_sim, case_i, M, ...
            key, BUNDLE, BASE_FEATURES, CFG, OUT);
        fprintf('  completed %d windows x %d models in %.1f s.\n', ...
            height(F), numel(unique(T_pred.strategy_name)), toc(t0));
        T_all = vertcat_compatible(T_all, T_pred);
        if CFG.SaveAllMaps
            plot_condition_maps(T_pred, F, CFG, OUT);
        end
    end
end

T_all = add_analysis_regions(T_all, CFG);

%% Summaries and figures

T_overall = summarize_long(T_all, "strategy_name");
T_by_case = summarize_long(T_all, ["strategy_name","case_name","field_regime"]);
T_by_M = summarize_long(T_all, ["strategy_name","M"]);
T_by_side = summarize_long(T_all, ["strategy_name","material_side"]);
T_by_roi = summarize_long(T_all, ["strategy_name","roi_name"]);
T_by_distance = summarize_long(T_all, ["strategy_name","distance_bin"]);
T_by_region = summarize_long(T_all, ["strategy_name","analysis_region"]);
T_by_case_M = summarize_long(T_all, ["strategy_name","case_name","field_regime","M"]);

writetable(T_all, fullfile(OUT.table_dir, 'test49_kwave_dense_patch_level_predictions.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test49_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test49_summary_by_case.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test49_summary_by_M.csv'));
writetable(T_by_side, fullfile(OUT.table_dir, 'test49_summary_by_soft_hard.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir, 'test49_summary_by_roi.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test49_summary_by_distance.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test49_summary_by_analysis_region.csv'));
writetable(T_by_case_M, fullfile(OUT.table_dir, 'test49_summary_by_case_M.csv'));
save(fullfile(OUT.data_dir, 'test49_kwave_dense_req_test38_validation.mat'), ...
    'T_all','T_overall','T_by_case','T_by_M','T_by_side','T_by_roi', ...
    'T_by_distance','T_by_region','T_by_case_M','CFG','-v7.3');

plot_overall_ranking(T_overall, OUT);
plot_group_heatmap(T_by_case, "case_name", OUT, 'test49_mape_by_case.png');
plot_group_heatmap(T_by_M, "M", OUT, 'test49_mape_by_M.png');
plot_soft_hard_roi(T_by_side, T_by_roi, OUT);
plot_distance_curves(T_by_distance, OUT);
plot_summary_by_case_M(T_by_case_M, OUT);

print_summary(T_overall, T_by_side, T_by_roi, T_by_distance);
fprintf('\nTables: %s\nFigures: %s\nTest 49 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration helpers

function KW = default_kwave_config()
KW = struct();
KW.f0 = 500;
KW.cs_bg = 2.0;
KW.cs_inc = 3.0;
KW.inclusion_center_m = [0.025, 0.025];
KW.inclusion_radius_m = 0.010;
KW.REQ.M_list = [2 3 4];
KW.ROI(1).name = "inclusion_roi";
KW.ROI(1).true_cs = KW.cs_inc;
KW.ROI(1).xlim_m = [0.021, 0.029];
KW.ROI(1).zlim_m = [0.021, 0.029];
KW.ROI(2).name = "background_roi";
KW.ROI(2).true_cs = KW.cs_bg;
KW.ROI(2).xlim_m = [0.006, 0.014];
KW.ROI(2).zlim_m = [0.021, 0.029];
KW.cases = [
    struct('case_name',"Directional 2D", 'folder',"2D-SS", 'wave_model',"SingleWave")
    struct('case_name',"Diffuse 2D", 'folder',"2D-diffuse", 'wave_model',"Diffuse2D")
    struct('case_name',"Projected diffuse 3D", 'folder',"3D-diffuse", 'wave_model',"Diffuse3D")
    struct('case_name',"3D-rev", 'folder',"3D-rev", 'wave_model',"Rev3D")
    ];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.Mode == "full"
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_49_kwave_dense_req_test38_validation');
else
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_49_kwave_dense_req_test38_validation', CFG.Mode);
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.cache_dir = fullfile(OUT.data_dir, 'condition_cache');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function B = load_test38_bundle(root_dir, CFG)
src = string(CFG.ModelSource);
if exist(src, 'file') == 2
    file = char(src);
else
    switch lower(src)
        case "full"
            file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
                'models', 'test38_velocity_field_diverse_q_models.mat');
        case "medium"
            file = fullfile(root_dir, 'outputs', 'model_registry', ...
                'test38_velocity_field_diverse_q_training', ...
                'test38__velocity_field_diverse_q__medium_bundle.mat');
        case "quick"
            file = fullfile(root_dir, 'outputs', 'model_registry', ...
                'test38_velocity_field_diverse_q_training', ...
                'test38__velocity_field_diverse_q__quick_bundle.mat');
        otherwise
            error('Unknown ADAPTIVE_REQ_TEST49_MODEL_SOURCE=%s.', src);
    end
end
assert(exist(file,'file') == 2, 'Missing Test38 bundle: %s', file);
fprintf('Loading frozen bundle: %s\n', file);
S = load(file, 'MODEL_BUNDLE');
B = S.MODEL_BUNDLE;
B.bundle_file = file;
end

function validate_sources(root_dir, CFG, BUNDLE, BASE_FEATURES)
assert(isfield(BUNDLE,'MODELS') && isfield(BUNDLE.MODELS,'composition'));
for ci = 1:min(1,numel(CFG.KW.cases))
    case_i = CFG.KW.cases(ci);
    data_file = fullfile(root_dir, 'data', 'k-wave', char(case_i.folder), ...
        sprintf('data_%dHz.mat', CFG.KW.f0));
    assert(exist(data_file,'file') == 2, 'Missing k-Wave field: %s', data_file);
end
assert_no_forbidden_predictors(BASE_FEATURES);
end

%% Dense REQ extraction and model application

function [field, dinf, cfg] = load_kwave_field(root_dir, CFG, case_i)
data_file = fullfile(root_dir, 'data', 'k-wave', char(case_i.folder), ...
    sprintf('data_%dHz.mat', CFG.KW.f0));
assert(exist(data_file,'file') == 2, 'Missing k-Wave data file: %s', data_file);
S = load(data_file, 'Vz_mg_ph', 'dinf');
field = S.Vz_mg_ph;
dinf = S.dinf;
cfg = struct();
cfg.f0 = CFG.KW.f0;
cfg.cs_bg = CFG.KW.cs_bg;
cfg.dx = dinf.dx;
cfg.dz = dinf.dy;
cfg.UseParfor = false;
end

function [T_pred, F] = evaluate_kwave_condition(field, dinf, cfg, case_i, M, key, ...
    BUNDLE, BASE_FEATURES, CFG, OUT)
cache_file = fullfile(OUT.cache_dir, sprintf('%s__req_%s_os%g_min%g_smooth%g.mat', ...
    key, CFG.REQ.Nbins, CFG.REQ.Nbins_auto_oversample, CFG.REQ.Nbins_min, CFG.REQ.SmoothSigma));
if exist(cache_file,'file') == 2
    S = load(cache_file, 'F', 'cache_sig');
    if isfield(S,'cache_sig') && isequaln(S.cache_sig, cache_signature(cfg, case_i, M, CFG))
        F = S.F;
        fprintf('  reused dense REQ cache (%d windows).\n', height(F));
    else
        F = extract_dense_features(field, dinf, cfg, case_i, M, key, CFG);
        cache_sig = cache_signature(cfg, case_i, M, CFG); %#ok<NASGU>
        save(cache_file, 'F', 'cache_sig', '-v7.3');
    end
else
    F = extract_dense_features(field, dinf, cfg, case_i, M, key, CFG);
    cache_sig = cache_signature(cfg, case_i, M, CFG); %#ok<NASGU>
    save(cache_file, 'F', 'cache_sig', '-v7.3');
end
F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
[T_pred, F] = apply_test38_models(F, BUNDLE.MODELS, BASE_FEATURES, CFG);
end

function F = extract_dense_features(field, dinf, cfg, case_i, M, key, CFG)
feat = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', CFG.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
step_x = max(1, round(CFG.TargetStepM / cfg.dx));
step_z = max(1, round(CFG.TargetStepM / cfg.dz));
req_options = {'Nbins', CFG.REQ.Nbins, ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, 'smooth_sigma', CFG.REQ.SmoothSigma};
fprintf('  REQ extraction: M=%g, StepX=%d, StepZ=%d, field=%s...\n', ...
    M, step_x, step_z, mat2str(size(field)));
t_req = tic;
O = adaptive_req.estimators.req_estimator_map(field, cfg, feat, ...
    'StepX', step_x, 'StepZ', step_z, ...
    'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseWindowParfor, ...
    'StoreReqCurves', false, 'Verbose', false);
F = O.feature_table;
fprintf('  REQ/features completed in %.1f s for %d windows.\n', toc(t_req), height(F));

F.condition_key = repmat(string(key), height(F), 1);
F.case_name = repmat(string(case_i.case_name), height(F), 1);
F.geometry = repmat("kwave_circular_inclusion_2_3", height(F), 1);
F.geometry_type = repmat("inclusion", height(F), 1);
F.field_regime = repmat(map_kwave_regime(case_i.wave_model), height(F), 1);
F.wave_model = repmat(string(case_i.wave_model), height(F), 1);
F.f0 = cfg.f0 * ones(height(F),1);
F.SIM_f0 = F.f0;
F.dx = cfg.dx * ones(height(F),1);
F.dz = cfg.dz * ones(height(F),1);
F.M = M * ones(height(F),1);
F.REQ_M = F.M;
F.REQ_cs_guess = feature_cs_guess(feat, CFG.cs_guess) * ones(height(F),1);
F.REQ_StepX = step_x * ones(height(F),1);
F.REQ_StepZ = step_z * ones(height(F),1);
F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
F.true_SWS = true_cs_at_points(F.x_center_m, F.z_center_m, CFG.KW);
F.cs_true = F.true_SWS;
F.k_true = 2*pi*cfg.f0 ./ F.true_SWS;
F.q_theory_prior = repmat(theory_q_kwave(F, case_i.wave_model, cfg.f0), height(F), 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
F.q_oracle = nan(height(F),1);
for i = 1:height(F)
    F.q_oracle(i) = invert_mapping_to_q(F.req_mapping{i}, F.k_true(i));
end
F.sws_oracle = q_to_sws(F.req_mapping, F.q_oracle, F.f0);
F.distance_to_interface_mm = abs(hypot(F.x_center_m - CFG.KW.inclusion_center_m(1), ...
    F.z_center_m - CFG.KW.inclusion_center_m(2)) - CFG.KW.inclusion_radius_m) * 1e3;
F.distance_bin = distance_bin(F.distance_to_interface_mm);
F.material_side = repmat("soft", height(F), 1);
F.material_side(F.true_SWS > (CFG.KW.cs_bg + CFG.KW.cs_inc)/2) = "hard";
F.roi_name = classify_kwave_roi(F.x_center_m, F.z_center_m, CFG.KW.ROI);
F.patch_purity = nan(height(F),1);
F.purity_bin = repmat("unknown_kWave", height(F), 1);
F.is_mixed = false(height(F),1);
F.is_strong_mixed = false(height(F),1);
F.REQ_Nbins_effective = ensure_nbins_effective(F);
end

function sig = cache_signature(cfg, case_i, M, CFG)
sig = struct();
sig.case_name = string(case_i.case_name);
sig.wave_model = string(case_i.wave_model);
sig.f0 = cfg.f0;
sig.dx = cfg.dx;
sig.dz = cfg.dz;
sig.M = M;
sig.TargetStepM = CFG.TargetStepM;
sig.REQ = CFG.REQ;
sig.cs_guess = CFG.cs_guess;
end

function F = ensure_predictor_columns(F, features, CFG)
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
                F.REQ_StepX = max(1, round(CFG.TargetStepM ./ F.dx));
            case "REQ_StepZ"
                F.REQ_StepZ = max(1, round(CFG.TargetStepM ./ F.dz));
            case "TargetStepM"
                F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
            otherwise
                error('Required Test38 predictor missing in dense k-Wave table: %s', f);
        end
    end
end
end

function [T_out, F] = apply_test38_models(F, MODELS, BASE_FEATURES, CFG)
P = predict_composition(MODELS.composition, F, BASE_FEATURES);
F.predicted_patch_purity = P.predicted_patch_purity;
F.p_mixed = P.p_mixed;
F.p_strong_mixed = P.p_strong_mixed;

parts = {};
for name = string(CFG.ModelsToEvaluate(:))'
    if ~ismember(name, string(MODELS.model_names(:)))
        warning('Skipping unavailable Test38 model: %s', name);
        continue;
    end
    R = base_prediction_table(F, CFG.ModelLabelPrefix + "_" + name);
    switch name
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_theory"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory_composition, F);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "delta_q_theory_composition"
            delta = predict_q_model(MODELS.q.delta_q_theory_composition, F, false);
            q_pred = clamp01(F.q_theory_prior + delta);
            sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
        case "delta_logk_theory_composition"
            delta = predict_q_model(MODELS.q.delta_logk_theory_composition, F, false);
            k_theory = 2*pi*F.f0 ./ F.sws_theory;
            k_pred = k_theory .* exp(delta);
            sws_pred = 2*pi*F.f0 ./ k_pred;
            q_pred = arrayfun(@(i) invert_mapping_to_q(F.req_mapping{i}, k_pred(i)), (1:height(F))');
        otherwise
            warning('Skipping unsupported Test38 model in Test49: %s', name);
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

function R = base_prediction_table(F, strategy_name)
keep = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "wave_model","M","REQ_M","f0","dx","dz","REQ_StepX","REQ_StepZ", ...
    "TargetStepM","map_iz","map_ix","x_center_m","z_center_m","true_SWS", ...
    "cs_true","roi_name","material_side","distance_to_interface_mm", ...
    "distance_bin","purity_bin","patch_purity","predicted_patch_purity", ...
    "p_mixed","p_strong_mixed","q_theory_prior","sws_theory"];
keep = keep(ismember(keep, string(F.Properties.VariableNames)));
R = F(:, cellstr(keep));
R.x = F.x_center_m;
R.z = F.z_center_m;
R.model_family = repmat("Test38_frozen_dense_kWave", height(F), 1);
R.bundle_id = repmat("test38_" + strategy_name, height(F), 1);
R.strategy_name = repmat(string(strategy_name), height(F), 1);
end

%% Predictors

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

%% Summaries and plots

function T = add_analysis_regions(T, CFG)
T.analysis_region = repmat("other", height(T), 1);
T.analysis_region(T.distance_to_interface_mm <= 1) = "interface_0_1mm";
T.analysis_region(T.distance_to_interface_mm > 1 & T.distance_to_interface_mm <= 2) = "interface_1_2mm";
T.analysis_region(T.distance_to_interface_mm > 2 & T.distance_to_interface_mm <= 4) = "near_interface_2_4mm";
T.analysis_region(T.distance_to_interface_mm > 4 & T.material_side == "soft") = "soft_core_gt4mm";
T.analysis_region(T.distance_to_interface_mm > 4 & T.material_side == "hard") = "hard_core_gt4mm";
T.analysis_region(T.roi_name == "background_roi") = "background_roi";
T.analysis_region(T.roi_name == "inclusion_roi") = "inclusion_roi";
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
over = splitapply(@(x) 100*mean(x > 0,'omitnan'), T.sws_signed_error_pct, G);
he10 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error10, G);
he20 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
mean_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
std_sws = splitapply(@std_omitnan, T.sws_pred, G);
S = [groups table(n,mape,medae,signed,medsigned,under,over,he10,he20,mean_sws,std_sws, ...
    'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
    'median_signed_error_pct','underestimate_pct','overestimate_pct','high_error10_pct', ...
    'high_error20_pct','mean_sws_pred','std_sws_pred'})];
end

function plot_overall_ranking(T, OUT)
T = sortrows(T, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 12]);
bar(categorical(T.strategy_name), T.MAPE_pct);
ylabel('MAPE (%)'); grid on; xtickangle(35);
title('Test 49 dense k-Wave REQ: overall MAPE','Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test49_overall_model_ranking.png'));
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
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13]);
imagesc(Z); cb = colorbar; ylabel(cb,'MAPE (%)'); colormap(parula);
set(gca,'XTick',1:numel(groups),'XTickLabel',groups,'YTick',1:numel(strategies), ...
    'YTickLabel',strategies,'TickLabelInterpreter','none');
xtickangle(30); title("Dense k-Wave MAPE by " + group_var, 'Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, filename));
end

function plot_soft_hard_roi(T_side, T_roi, OUT)
strategies = unique(T_side.strategy_name, 'stable');
side_groups = unique(string(T_side.material_side), 'stable');
roi_groups = unique(string(T_roi.roi_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 13]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
nexttile;
bar(categorical(side_groups), grouped_values(T_side, "material_side", strategies));
legend(strategies,'Interpreter','none','Location','bestoutside'); ylabel('MAPE (%)'); grid on;
title('Soft/hard MAPE','FontWeight','normal');
nexttile;
bar(categorical(roi_groups), grouped_values(T_roi, "roi_name", strategies));
legend(strategies,'Interpreter','none','Location','bestoutside'); ylabel('MAPE (%)'); grid on;
title('ROI MAPE','FontWeight','normal'); xtickangle(30);
export_fig(fig, fullfile(OUT.figure_dir, 'test49_soft_hard_roi_mape.png'));
end

function plot_distance_curves(T_dist, OUT)
strategies = unique(T_dist.strategy_name, 'stable');
order = ["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
fig = figure('Color','w','Units','centimeters','Position',[1 1 22 12]);
hold on;
for si = 1:numel(strategies)
    y = nan(size(order));
    for i = 1:numel(order)
        idx = T_dist.strategy_name == strategies(si) & string(T_dist.distance_bin) == order(i);
        if any(idx), y(i) = mean(T_dist.MAPE_pct(idx), 'omitnan'); end
    end
    plot(1:numel(order), y, '-o', 'LineWidth', 1.1, 'DisplayName', strategies(si));
end
set(gca,'XTick',1:numel(order),'XTickLabel',order); xtickangle(25);
ylabel('MAPE (%)'); xlabel('Distance to interface'); grid on;
legend('Interpreter','none','Location','bestoutside');
title('Dense k-Wave error versus distance to interface','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test49_mape_vs_distance.png'));
end

function plot_summary_by_case_M(T_by_case_M, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 13]);
strategies = unique(T_by_case_M.strategy_name, 'stable');
cases = unique(T_by_case_M.case_name, 'stable');
tiledlayout(fig,1,numel(cases),'TileSpacing','compact','Padding','compact');
for ci = 1:numel(cases)
    ax = nexttile;
    hold(ax,'on');
    for si = 1:numel(strategies)
        X = T_by_case_M(T_by_case_M.case_name == cases(ci) & ...
            T_by_case_M.strategy_name == strategies(si),:);
        if isempty(X), continue; end
        X = sortrows(X,'M');
        plot(ax, X.M, X.MAPE_pct, '-o', 'DisplayName', strategies(si), 'LineWidth', 1.0);
    end
    title(ax, cases(ci), 'Interpreter','none','FontWeight','normal');
    xlabel(ax,'REQ M'); ylabel(ax,'MAPE (%)'); grid(ax,'on');
end
legend('Interpreter','none','Location','bestoutside');
export_fig(fig, fullfile(OUT.figure_dir, 'test49_mape_vs_M_by_case.png'));
end

function plot_condition_maps(T_pred, F, CFG, OUT)
conds = unique(T_pred(:, {'case_name','field_regime','M'}), 'rows', 'stable');
persistent n_maps_saved
if isempty(n_maps_saved), n_maps_saved = 0; end
for ci = 1:height(conds)
    if n_maps_saved >= CFG.MaxMapConditions, return; end
    idxc = T_pred.case_name == conds.case_name(ci) & ...
        T_pred.field_regime == conds.field_regime(ci) & T_pred.M == conds.M(ci);
    Xc = T_pred(idxc,:);
    Fc = F(F.case_name == conds.case_name(ci) & F.field_regime == conds.field_regime(ci) & F.M == conds.M(ci),:);
    if isempty(Xc) || isempty(Fc), continue; end
    [true_map,nz,nx] = rows_to_grid(Fc, Fc.true_SWS);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 34 26]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_grid(Fc, Fc.sws_theory, nz, nx), 'Theory SWS', 'SWS (m/s)');
    plot_map(nexttile(tl), rows_to_grid(Fc, Fc.predicted_patch_purity, nz, nx), 'Predicted patch purity', 'Probability');
    plot_map(nexttile(tl), rows_to_grid(Fc, Fc.p_mixed, nz, nx), 'Predicted mixedness', 'Probability');
    strategies = unique(Xc.strategy_name, 'stable');
    for si = 1:numel(strategies)
        X = Xc(Xc.strategy_name == strategies(si),:);
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), strategies(si) + ' SWS', 'SWS (m/s)');
    end
    for si = 1:numel(strategies)
        X = Xc(Xc.strategy_name == strategies(si),:);
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, nz, nx), ...
            strategies(si) + ' abs error', 'Abs error (%)');
    end
    plot_map(nexttile(tl), rows_to_grid(Fc, Fc.distance_to_interface_mm, nz, nx), ...
        'Distance to interface', 'Distance (mm)');
    plot_map(nexttile(tl), rows_to_grid(Fc, double(Fc.material_side == "hard"), nz, nx), ...
        'Hard-region mask', 'Mask');
    title(tl, sprintf('%s | %s | M=%g | dense REQ step %.3f mm', ...
        conds.case_name(ci), conds.field_regime(ci), conds.M(ci), 1e3*CFG.TargetStepM), ...
        'Interpreter','none');
    outdir = fullfile(OUT.map_dir, sanitize(conds.case_name(ci)), ...
        sanitize(conds.field_regime(ci)), "M" + string(conds.M(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    out_file = fullfile(outdir, sprintf('test49_dense_kwave__%s__%s__M%g.png', ...
        sanitize(conds.case_name(ci)), sanitize(conds.field_regime(ci)), conds.M(ci)));
    export_fig(fig, out_file);
    n_maps_saved = n_maps_saved + 1;
end
end

function print_summary(T_overall, T_by_side, T_by_roi, T_by_distance)
T = sortrows(T_overall, 'MAPE_pct', 'ascend');
fprintf('\n================ Test 49 dense k-Wave summary ================\n');
disp(T(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
fprintf('Best dense k-Wave model: %s (MAPE %.2f%%, signed %.2f%%, HE20 %.2f%%).\n', ...
    T.strategy_name(1), T.MAPE_pct(1), T.mean_signed_error_pct(1), T.high_error20_pct(1));
hard = T_by_side(T_by_side.material_side == "hard",:);
soft = T_by_side(T_by_side.material_side == "soft",:);
if ~isempty(hard)
    hard = sortrows(hard,'MAPE_pct','ascend');
    fprintf('Best hard-side MAPE: %s %.2f%% (signed %.2f%%).\n', ...
        hard.strategy_name(1), hard.MAPE_pct(1), hard.mean_signed_error_pct(1));
end
if ~isempty(soft)
    soft = sortrows(soft,'MAPE_pct','ascend');
    fprintf('Best soft-side MAPE: %s %.2f%% (signed %.2f%%).\n', ...
        soft.strategy_name(1), soft.MAPE_pct(1), soft.mean_signed_error_pct(1));
end
fprintf('ROI/distance summaries saved for core/interface inspection.\n');
fprintf('==============================================================\n');
end

%% k-Wave diagnostics

function cs = true_cs_at_points(x, z, KW)
r = hypot(x - KW.inclusion_center_m(1), z - KW.inclusion_center_m(2));
cs = KW.cs_bg * ones(numel(x), 1);
cs(r <= KW.inclusion_radius_m) = KW.cs_inc;
end

function roi = classify_kwave_roi(x, z, rois)
roi = repmat("outside_roi", numel(x), 1);
for ri = 1:numel(rois)
    in_roi = x >= rois(ri).xlim_m(1) & x <= rois(ri).xlim_m(2) & ...
        z >= rois(ri).zlim_m(1) & z <= rois(ri).zlim_m(2);
    roi(in_roi) = string(rois(ri).name);
end
end

function r = map_kwave_regime(w)
switch string(w)
    case "SingleWave", r = "directional_2D";
    case "Diffuse2D", r = "diffuse_2D";
    case "Diffuse3D", r = "diffuse_3D";
    otherwise, r = "partial_3D";
end
end

function q = theory_q_kwave(T, wave_model, f0)
switch string(wave_model)
    case "SingleWave", type = "SingleWave";
    case "Diffuse2D", type = "Diffuse2D";
    case "Diffuse3D", type = "Diffuse3D";
    otherwise, type = "Partial3D";
end
if type == "Partial3D"
    q = 0.5*(theory_one(T,"Diffuse2D",f0) + theory_one(T,"Diffuse3D",f0));
else
    q = theory_one(T,type,f0);
end
q = clamp01(q);
end

function q = theory_one(T, type, f0)
o = adaptive_req.theory.q_theory_REQ_discrete_shearUZ(T.dx(1), T.dz(1), ...
    f0, T.REQ_cs_guess(1), 'M', T.REQ_M(1), 'Gamma', 1, ...
    'PadFactor', 1, 'Nbins', 'auto', 'SmoothSigma', 1, ...
    'TheoryMode', 'S2D', 'FieldType', type, 'Plot', false);
q = o.q_th;
end

function c = feature_cs_guess(feat, fallback)
if isfield(feat, 'cs_guess_used')
    c = feat.cs_guess_used;
elseif isfield(feat, 'cs_guess')
    c = feat.cs_guess;
else
    c = fallback;
end
end

%% Numeric utilities

function q = clamp01(q)
q = min(max(q, 0.001), 0.999);
end

function p = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
p = score(:,idx);
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

function n = ensure_nbins_effective(F)
n = nan(height(F),1);
if ismember("REQ_Nbins_effective", string(F.Properties.VariableNames))
    n = F.REQ_Nbins_effective;
    return;
end
for i = 1:height(F)
    try
        n(i) = F.req_mapping{i}.Nbins_effective;
    catch
        try
            n(i) = numel(F.req_mapping{i}.k_cent);
        catch
            n(i) = NaN;
        end
    end
end
end

function bin = distance_bin(d)
bin = repmat("gt_8mm", numel(d), 1);
bin(d <= 8) = "4_8mm";
bin(d <= 4) = "2_4mm";
bin(d <= 2) = "1_2mm";
bin(d <= 1) = "0_1mm";
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map(ax, Z, ttl, cb_label)
imagesc(ax, Z); axis(ax,'image'); axis(ax,'off');
cb = colorbar(ax); ylabel(cb, cb_label);
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
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

function T = vertcat_compatible(varargin)
nonempty = {};
for i = 1:nargin
    if ~isempty(varargin{i}), nonempty{end+1} = varargin{i}; end %#ok<AGROW>
end
if isempty(nonempty), T = table(); return; end
vars = strings(1,0);
for i = 1:numel(nonempty)
    vars = [vars string(nonempty{i}.Properties.VariableNames)]; %#ok<AGROW>
end
vars = unique(vars, 'stable');
string_vars = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "wave_model","roi_name","material_side","distance_bin","purity_bin", ...
    "model_family","bundle_id","strategy_name","analysis_region"];
logical_vars = ["high_error10","high_error20"];
parts = cell(numel(nonempty),1);
for i = 1:numel(nonempty)
    A = nonempty{i};
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

function assert_no_forbidden_predictors(features)
bad = ["true","oracle","purity","mixed","confidence","error","pred", ...
    "sws","cs_","k_true","q_local","q_pred","req_mapping","patch_idx", ...
    "map_ix","map_iz","cx","cz","x_center","z_center","condition", ...
    "distance"];
features = lower(string(features));
for b = bad
    hit = features(contains(features,b));
    assert(isempty(hit), 'Forbidden predictor in Test38 bundle: %s', strjoin(hit, ', '));
end
end

function s = std_omitnan(x)
s = std(x, 'omitnan');
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
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
if v == "", x = string(default); else, x = v; end
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
