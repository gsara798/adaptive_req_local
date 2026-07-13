%% analyze_test_62_dense_test61_benefit_harm_validation.m
% Test 62: dense REQ validation of the frozen Test 61 benefit/harm gate.
%
% This script recalculates dense REQ maps on representative conditions,
% applies the frozen Test 38 q_spectrum_plus_composition estimator at M=2
% and M=3, then applies the already-trained Test 61 correction/reliability
% gate. No base q model and no correction gate is retrained here.

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
write_config_json(CFG, fullfile(OUT.root_dir,'test62_configuration.json'));

fprintf('\nTest 62: dense validation of frozen Test 61 benefit/harm gate\n');
fprintf('Mode: %s | target REQ step %.3f mm | save all maps: %d\n', ...
    CFG.Mode, 1e3*CFG.TargetStepM, CFG.SaveAllMaps);
fprintf('No training. Frozen model source: %s\n', CFG.ModelBundleFile);
fprintf('Frozen Test 61 gate source: %s\n', CFG.Test61ResultFile);
fprintf('REQ profile: %s | Nbins=%s | oversample=%g | Nbins min=%g\n', ...
    CFG.REQ.Profile, string(CFG.REQ.Nbins), CFG.REQ.Nbins_auto_oversample, CFG.REQ.Nbins_min);
fprintf('At least one source is forced/aligned in every regime.\n');

BUNDLE = load_test38_bundle(CFG);
MODELS = BUNDLE.MODELS;
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);
GATE = load_test61_gate(CFG);

if CFG.ValidateOnly
    CONDITIONS = build_condition_specs(CFG);
    CONDITIONS = CONDITIONS(1:min(1,numel(CONDITIONS)));
else
    CONDITIONS = build_condition_specs(CFG);
end
fprintf('Dense validation conditions: %d | M values: %s\n', ...
    numel(CONDITIONS), mat2str(CFG.M));

strategy_parts = {};
base_parts = {};
for ci = 1:numel(CONDITIONS)
    CND = CONDITIONS(ci);
    condition_timer = tic;
    long_parts = {};
    for mi = 1:numel(CFG.M)
        M = CFG.M(mi);
        key = sprintf('%s__f%g__%s__dx%gum__M%g', CND.case.case_id, ...
            CND.f0, CND.regime.regime_id, round(1e6*CFG.dx), M);
        [T_pred, T_base] = evaluate_dense_condition(CND, M, key, MODELS, ...
            BASE_FEATURES, CFG, OUT);
        T_pred = add_roi_labels(T_pred, CFG);
        T_pred = add_operational_features_long(T_pred, GATE.CFG);
        long_parts{end+1,1} = T_pred; %#ok<AGROW>
        base_parts{end+1,1} = remove_cell_columns(T_pred); %#ok<AGROW>
        fprintf('    M=%g: %d dense patches.\n', M, height(T_base));
    end
    T_long = vertcat(long_parts{:});
    W = build_m2_m3_wide(T_long, GATE.CFG);
    assert(~isempty(W), 'Could not align M=2/M=3 dense rows for %s.', CND.case.case_id);
    T_strat = apply_strategies(W, GATE.MODEL, GATE.CFG);
    strategy_parts{end+1,1} = remove_cell_columns(T_strat); %#ok<AGROW>

    fprintf('[%d/%d] %s f=%g %s: %d aligned dense pixels in %.1f s.\n', ...
        ci, numel(CONDITIONS), CND.case.case_id, CND.f0, CND.regime.regime_id, ...
        height(W), toc(condition_timer));
    if CFG.SaveAllMaps && should_save_map(ci, CFG)
        plot_test62_maps(T_strat, T_long, CND, OUT, CFG);
        plot_roi_overlay_pair(T_long(abs(T_long.M-2)<1e-8,:), string(W.base_condition(1)), OUT, CFG);
    end
end

T_all = vertcat(strategy_parts{:});
T_base_long = vertcat(base_parts{:});
T_overall = summarize_strategy_predictions(T_all, "strategy_name");
T_by_case = summarize_strategy_predictions(T_all, ["strategy_name","case_id","seen_status"]);
T_by_regime = summarize_strategy_predictions(T_all, ["strategy_name","field_regime_ood"]);
T_by_frequency = summarize_strategy_predictions(T_all, ["strategy_name","f0"]);
T_by_roi = summarize_strategy_predictions(T_all(T_all.roi_region ~= "other",:), ...
    ["strategy_name","case_id","roi_region"]);
T_harm = summarize_harm(T_all);

writetable(T_all, fullfile(OUT.table_dir,'test62_dense_strategy_patch_level_results.csv'));
writetable(T_base_long, fullfile(OUT.table_dir,'test62_dense_base_M2_M3_long_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test62_strategy_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test62_strategy_summary_by_case.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test62_strategy_summary_by_regime.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir,'test62_strategy_summary_by_frequency.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir,'test62_strategy_summary_by_roi.csv'));
writetable(T_harm, fullfile(OUT.table_dir,'test62_correction_harm_summary.csv'));

safe_plot(@() plot_test62_summary_figures(T_overall, T_by_case, T_by_frequency, T_by_roi, T_harm, OUT), ...
    'summary figures');

save(fullfile(OUT.data_dir,'test62_compact_results.mat'), ...
    'CFG','T_overall','T_by_case','T_by_regime','T_by_frequency', ...
    'T_by_roi','T_harm','-v7.3');

print_test62_summary(T_overall, T_by_roi, T_harm, OUT);

fprintf('\nTables: %s\nFigures: %s\nTest 62 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST62_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST62_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST62_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST62_SAVE_ALL_MAPS', true);
CFG.UseParfor = env_true('ADAPTIVE_REQ_TEST62_USE_PARFOR', false);
CFG.TargetStepM = env_num('ADAPTIVE_REQ_TEST62_TARGET_STEP_M', 0.5e-3);
CFG.MaxMapConditions = env_num('ADAPTIVE_REQ_TEST62_MAX_MAP_CONDITIONS', 0);
CFG.ModelSource = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST62_MODEL_SOURCE'))));
if CFG.ModelSource == ""
    if CFG.QuickMode || CFG.ValidateOnly
        CFG.ModelSource = "medium";
    else
        CFG.ModelSource = "full";
    end
end
CFG.ModelBundleFile = resolve_model_bundle(root_dir, CFG.ModelSource);
CFG.Test61Source = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST62_TEST61_SOURCE'))));
if CFG.Test61Source == ""
    CFG.Test61Source = "full";
end
CFG.Test61ResultFile = resolve_test61_results(root_dir, CFG.Test61Source);
CFG.dx = 0.2e-3;
CFG.dz = 0.2e-3;
CFG.Lx = 0.05;
CFG.Lz = 0.05;
CFG.cs_guess = 3;
CFG.M = [2 3];
CFG.RandomSeed = 48001;
CFG.DistanceEdgesMm = [0 0.5 1 2 4 8 Inf];
CFG.CenterRoiSizeMm = 8;
CFG.SoftCoreMinDistanceMm = 8;
CFG.HardCoreMinDistanceMm = 4;
CFG.ModelsToPlot = "q_spectrum_plus_composition";
CFG.ModelsToEvaluate = "q_spectrum_plus_composition";
CFG.StrategiesToPlot = ["fixed_M2","fixed_M3","interface_penalized_selector", ...
    "learned_benefit_harm_gate","conservative_benefit_harm_gate"];
CFG.REQ.Profile = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST62_REQ_PROFILE'))));
if CFG.REQ.Profile == "", CFG.REQ.Profile = "test38_training"; end
CFG.REQ.Gamma = 1;
CFG.REQ.PadFactor = 1;
CFG.REQ.EdgeMode = 'valid';
CFG.REQ.Nbins = 'auto';
switch CFG.REQ.Profile
    case {"test38_training","training","matched"}
        CFG.REQ.Nbins_auto_oversample = 1;
        CFG.REQ.Nbins_min = 16;
    case {"dense_default","old_test48"}
        CFG.REQ.Nbins_auto_oversample = 2;
        CFG.REQ.Nbins_min = 128;
    otherwise
        error('Unknown ADAPTIVE_REQ_TEST62_REQ_PROFILE: %s', CFG.REQ.Profile);
end
CFG.REQ.Nbins_auto_oversample = env_num('ADAPTIVE_REQ_TEST62_NBINS_OVERSAMPLE', CFG.REQ.Nbins_auto_oversample);
CFG.REQ.Nbins_min = env_num('ADAPTIVE_REQ_TEST62_NBINS_MIN', CFG.REQ.Nbins_min);
CFG.REQ.SmoothSigma = 1;
if CFG.ValidateOnly
    CFG.M = [2 3];
    CFG.MaxMapConditions = min_positive_or_default(CFG.MaxMapConditions, 1);
end
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_62_dense_test61_benefit_harm_validation','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_62_dense_test61_benefit_harm_validation');
end
if CFG.ValidateOnly
    OUT.root_dir = fullfile(OUT.root_dir, 'validate');
end
OUT.table_dir = fullfile(OUT.root_dir,'tables');
OUT.figure_dir = fullfile(OUT.root_dir,'figures');
OUT.map_dir = fullfile(OUT.figure_dir,'maps_by_condition');
OUT.roi_dir = fullfile(OUT.figure_dir,'roi_overlays');
OUT.data_dir = fullfile(OUT.root_dir,'data');
OUT.field_dir = fullfile(OUT.data_dir,'field_cache');
OUT.req_dir = fullfile(OUT.data_dir,'dense_req_features');
dirs = string(struct2cell(OUT));
for d = dirs(:)'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function file = resolve_test61_results(root_dir, source)
source = string(source);
if isfile(source)
    file = char(source);
    return;
end
switch lower(source)
    case "quick"
        file = fullfile(root_dir, 'outputs', ...
            'test_61_learned_benefit_harm_gate', 'quick', 'data', ...
            'test61_learned_benefit_harm_gate_results.mat');
    otherwise
        file = fullfile(root_dir, 'outputs', ...
            'test_61_learned_benefit_harm_gate', 'data', ...
            'test61_learned_benefit_harm_gate_results.mat');
end
end

function file = resolve_model_bundle(root_dir, source)
source = string(source);
if isfile(source)
    file = char(source);
    return;
end
switch lower(source)
    case "quick"
        file = fullfile(root_dir, 'outputs', ...
            'test_38_velocity_field_diverse_q_training', 'quick', 'models', ...
            'test38_velocity_field_diverse_q_models.mat');
    case "medium"
        file = fullfile(root_dir, 'outputs', ...
            'test_38_velocity_field_diverse_q_training', 'medium', 'models', ...
            'test38_velocity_field_diverse_q_models.mat');
    otherwise
        file = fullfile(root_dir, 'outputs', ...
            'test_38_velocity_field_diverse_q_training', 'models', ...
            'test38_velocity_field_diverse_q_models.mat');
end
end

function B = load_test38_bundle(CFG)
assert(exist(CFG.ModelBundleFile,'file') == 2, ...
    'Missing Test38 model bundle: %s', CFG.ModelBundleFile);
fprintf('Loading frozen Test 38 bundle (%s): %s\n', CFG.ModelSource, CFG.ModelBundleFile);
S = load(CFG.ModelBundleFile, 'MODEL_BUNDLE');
B = S.MODEL_BUNDLE;
fprintf('Loaded bundle with %d model entries.\n', numel(B.MODELS.model_names));
end

function G = load_test61_gate(CFG)
assert(exist(CFG.Test61ResultFile,'file') == 2, ...
    'Missing Test61 trained gate results: %s', CFG.Test61ResultFile);
fprintf('Loading frozen Test61 correction gate: %s\n', CFG.Test61ResultFile);
S = load(CFG.Test61ResultFile, 'CFG', 'MODEL');
G = struct();
G.CFG = S.CFG;
G.MODEL = S.MODEL;
required = ["residual","risk","interface_penalized_selector","benefit_harm"];
for r = required
    assert(isfield(G.MODEL, r), 'Test61 MODEL is missing field: %s', r);
end
fprintf('Loaded Test61 gate. Residual/risk/benefit-harm models are frozen.\n');
end

function assert_no_forbidden_predictors(features)
bad_patterns = ["true","oracle","purity","mixed","confidence","error", ...
    "pred","sws","cs_","k_true","q_local","q_pred","q_theory", ...
    "req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition"];
features = lower(string(features));
for p = bad_patterns
    hit = features(contains(features, p));
    assert(isempty(hit), 'Forbidden Test38 predictor detected: %s', strjoin(hit, ', '));
end
fprintf('Base predictors passed leakage guard (%d predictors).\n', numel(features));
end

%% Conditions

function COND = build_condition_specs(CFG)
cases_seen = [
    homogeneous_case("homogeneous_cs2", 2, "seen")
    homogeneous_case("homogeneous_cs3", 3, "seen")
    bilayer_case("bilayer_2_3", 2, 3, 0, "seen")
    inclusion_case("circular_inclusion_2_3", 2, 3, [0.025 0.025], 0.008, "seen")
    thin_layer_case("thin_layer_2_4", 2, 4, "seen")
    two_inclusion_case("two_inclusions_2_4", 2, 4, "seen")
    three_material_case("three_material_2_3_4", 2, 3, 4, "seen")
    ];
cases_unseen = [
    bilayer_case("ood_bilayer_2p25_3p75", 2.25, 3.75, 0, "unseen_velocity")
    inclusion_case("ood_inclusion_2p25_3p75", 2.25, 3.75, [0.025 0.025], 0.008, "unseen_velocity")
    ellipse_case("ood_ellipse_offcenter_2_4", 2, 4, [0.030 0.021], [0.009 0.005], pi/7, "unseen_geometry")
    ];
regimes = [
    regime("directional_2D_angle15","directional_2D",true,1,true,"ranges","fibonacci",[deg2rad(15) deg2rad(15)],[0 2*pi],[0 pi],15)
    regime("diffuse_2D_seed17","diffuse_2D",true,128,true,"ranges","random",[0 2*pi],[0 2*pi],[0 pi],17)
    regime("partial_3D_12src","partial_3D",false,12,true,"ranges","fibonacci",[0 2*pi],[0 2*pi],[0 pi],19)
    regime("diffuse_3D_seed23","diffuse_3D",false,128,true,"ranges","random",[0 2*pi],[0 2*pi],[0 pi],23)
    ];
if CFG.QuickMode
    cases = [cases_seen([1 2 3 4 5 7]); cases_unseen([1 2])];
    freqs = [375 475 575 675 475 575 475 675];
    reg_idx = [1 4 1 4 2 3 1 4];
    COND = repmat(empty_condition(), 0, 1);
    for i = 1:numel(cases)
        COND(end+1,1) = condition(cases(i), regimes(reg_idx(i)), freqs(i)); %#ok<AGROW>
    end
else
    cases = [cases_seen; cases_unseen];
    freqs = [375 475 575 675];
    COND = repmat(empty_condition(), 0, 1);
    for ci = 1:numel(cases)
        for ri = 1:numel(regimes)
            f0 = freqs(1 + mod(ci+ri-2, numel(freqs)));
            COND(end+1,1) = condition(cases(ci), regimes(ri), f0); %#ok<AGROW>
        end
    end
end
end

function X = empty_condition()
X = struct('case', empty_case("", ""), 'regime', regime("", "", true, 1, true, ...
    "ranges", "fibonacci", 0, [0 2*pi], [0 pi], 0), 'f0', NaN);
end

function X = condition(C, R, f0)
X = struct('case', C, 'regime', R, 'f0', f0);
end

function C = homogeneous_case(id, cs, seen)
C = empty_case(id, "homogeneous");
C.cs_values = cs; C.cs_bg = cs; C.seen_status = string(seen); C.mask_builder = "homogeneous";
end

function C = bilayer_case(id, cs_low, cs_high, angle, seen)
C = empty_case(id, "bilayer");
C.cs_values = [cs_low cs_high]; C.cs_bg = cs_low; C.cs_inc = cs_high;
C.angle = angle; C.offset = 0.025; C.sigma_edge = 1e-6;
C.seen_status = string(seen); C.mask_builder = "bilayer";
end

function C = inclusion_case(id, cs_bg, cs_inc, center, radius, seen)
C = empty_case(id, "inclusion");
C.cs_values = [cs_bg cs_inc]; C.cs_bg = cs_bg; C.cs_inc = cs_inc;
C.center = center; C.radius = radius; C.sigma_edge = 1e-6;
C.seen_status = string(seen); C.mask_builder = "circle";
end

function C = ellipse_case(id, cs_bg, cs_inc, center, axes_m, angle, seen)
C = empty_case(id, "ellipse");
C.cs_values = [cs_bg cs_inc]; C.cs_bg = cs_bg; C.cs_inc = cs_inc;
C.center = center; C.axes = axes_m; C.angle = angle;
C.seen_status = string(seen); C.mask_builder = "ellipse_custom";
end

function C = two_inclusion_case(id, cs_bg, cs_inc, seen)
C = empty_case(id, "two_inclusions");
C.cs_values = [cs_bg cs_inc]; C.cs_bg = cs_bg; C.cs_inc = cs_inc;
C.centers = [0.018 0.024; 0.034 0.028]; C.radii = [0.006 0.005];
C.seen_status = string(seen); C.mask_builder = "two_circles";
end

function C = thin_layer_case(id, cs_bg, cs_layer, seen)
C = empty_case(id, "thin_layer");
C.cs_values = [cs_bg cs_layer]; C.cs_bg = cs_bg; C.cs_inc = cs_layer;
C.center_x = 0.025; C.thickness = 0.004;
C.seen_status = string(seen); C.mask_builder = "thin_layer_custom";
end

function C = three_material_case(id, cs_bg, cs_mid, cs_high, seen)
C = empty_case(id, "three_material");
C.cs_values = [cs_bg cs_mid cs_high]; C.cs_bg = cs_bg; C.cs_mid = cs_mid; C.cs_high = cs_high;
C.seen_status = string(seen); C.mask_builder = "three_material_custom";
end

function C = empty_case(id, fam)
C = struct('case_id',string(id),'case_family',string(fam),'seen_status',"seen", ...
    'cs_values',[],'cs_bg',NaN,'cs_inc',NaN,'cs_mid',NaN,'cs_high',NaN, ...
    'mask_builder',"",'angle',NaN,'offset',NaN,'sigma_edge',1e-6, ...
    'center',[NaN NaN],'radius',NaN,'axes',[NaN NaN], ...
    'centers',nan(0,2),'radii',nan(0,1),'center_x',NaN,'thickness',NaN);
end

function R = regime(id, label, is2d, nwaves, force_plane, sampling, method, angle2d, phi, theta, seed_offset)
R = struct('regime_id',string(id),'field_regime',string(label),'Is2D',is2d, ...
    'Nwaves',nwaves,'ForceInPlaneWave',force_plane,'SourceSampling',string(sampling), ...
    'AngularSamplingMethod',string(method),'AngleRange2D',angle2d, ...
    'PhiRange',phi,'ThetaRange',theta,'SeedOffset',seed_offset);
end

%% Dense evaluation

function [T_pred, F] = evaluate_dense_condition(CND, M, key, MODELS, BASE_FEATURES, CFG, OUT)
C = CND.case; R = CND.regime;
field_key = sprintf('%s__f%g__%s__dx%gum', C.case_id, CND.f0, R.regime_id, round(1e6*CFG.dx));
field_file = fullfile(OUT.field_dir, "field__" + sanitize(field_key) + ".mat");
if exist(field_file,'file') == 2
    S = load(field_file, 'sim', 'cfg_sim');
    sim = S.sim; cfg_sim = S.cfg_sim;
else
    cfg_sim = build_sim_cfg(CFG, C, R, CND.f0, 1000 + sum(double(char(field_key))));
    sim = adaptive_req.simulate.run_single_simulation(cfg_sim);
    case_spec = C; regime_spec = R; %#ok<NASGU>
    save(field_file, 'sim', 'cfg_sim', 'case_spec', 'regime_spec', '-v7.3');
end

req_file = fullfile(OUT.req_dir, "dense_req__" + sanitize(key) + ...
    sprintf('__step%gum__req_%s_os%g_min%g.mat', round(1e6*CFG.TargetStepM), ...
    sanitize(CFG.REQ.Profile), CFG.REQ.Nbins_auto_oversample, CFG.REQ.Nbins_min));
if exist(req_file,'file') == 2
    S = load(req_file, 'F');
    F = S.F;
else
    feat = adaptive_req.config.default_feature_config('M', M, ...
        'cs_guess', CFG.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
        'pad_factor', CFG.REQ.PadFactor);
    step_x = max(1, round(CFG.TargetStepM / cfg_sim.dx));
    step_z = max(1, round(CFG.TargetStepM / cfg_sim.dz));
    O = adaptive_req.estimators.req_estimator_map(sim.Uxz, cfg_sim, feat, ...
        'StepX', step_x, 'StepZ', step_z, ...
        'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
        'ReqOptions', {'Nbins',CFG.REQ.Nbins, ...
        'Nbins_auto_oversample',CFG.REQ.Nbins_auto_oversample, ...
        'Nbins_min',CFG.REQ.Nbins_min,'smooth_sigma',CFG.REQ.SmoothSigma}, ...
        'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
        'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseParfor, ...
        'StoreReqCurves', false, 'Verbose', false);
    F = O.feature_table;
    F = attach_metadata(F, sim, cfg_sim, C, R, M, key, O.win_size, CFG, step_x, step_z);
    save(req_file, 'F', '-v7.3');
end
F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
T_pred = apply_models(F, MODELS, BASE_FEATURES, "test62_dense_test61_validation", CFG.ModelsToEvaluate);
end

function cfg = build_sim_cfg(CFG, C, R, f0, condition_id)
cfg = struct();
cfg.Lx = CFG.Lx; cfg.Lz = CFG.Lz; cfg.dx = CFG.dx; cfg.dz = CFG.dz;
cfg.f0 = f0; cfg.cs_bg = C.cs_bg; cfg.cs_inc = max(C.cs_values);
cfg.Nwaves = R.Nwaves; cfg.Is2D = R.Is2D; cfg.WaveModel = 'spherical';
cfg.AngularSamplingMethod = char(R.AngularSamplingMethod);
cfg.ForceInPlaneWave = logical(R.ForceInPlaneWave);
cfg.SNR = Inf; cfg.AmpJitter = 0.05; cfg.DecayAlpha = 0; cfg.UseParfor = false;
cfg.Seed = CFG.RandomSeed + condition_id + R.SeedOffset;
cfg.PhiRange = R.PhiRange; cfg.ThetaRange = R.ThetaRange;
cfg.AngleRange2D = R.AngleRange2D; cfg.SourceSampling = char(R.SourceSampling);
cfg.MaskConfig = build_mask_config(CFG, C);
end

function mask_cfg = build_mask_config(CFG, C)
switch C.mask_builder
    case "homogeneous"
        masks = {};
    case "bilayer"
        masks = {struct('Type','bilayer','cs_inc',C.cs_inc,'Params', ...
            struct('Bi_Angle',C.angle,'Bi_Offset',C.offset,'SigmaEdge',C.sigma_edge))};
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
        xr = ca*X0 + sa*Z0; zr = -sa*X0 + ca*Z0;
        M = (xr./C.axes(1)).^2 + (zr./C.axes(2)).^2 <= 1;
    case "thin_layer"
        M = abs(X - C.center_x) <= C.thickness/2;
    case "three_circle"
        M = hypot(X - 0.017, Z - 0.031) <= 0.006;
    otherwise
        error('Unknown custom mask: %s', kind);
end
end

function F = attach_metadata(F, sim, cfg, C, R, M, key, win_size, CFG, step_x, step_z)
n = height(F);
F.dataset = repmat("test62_dense_test61_validation", n, 1);
F.condition_key = repmat(string(key), n, 1);
F.case_id = repmat(C.case_id, n, 1);
F.case_family = repmat(C.case_family, n, 1);
F.seen_status = repmat(C.seen_status, n, 1);
F.field_regime = repmat(R.field_regime, n, 1);
F.field_regime_ood = repmat(R.regime_id, n, 1);
F.f0 = cfg.f0 * ones(n,1);
F.M = M * ones(n,1);
F.REQ_M = M * ones(n,1);
F.REQ_StepX = step_x * ones(n,1);
F.REQ_StepZ = step_z * ones(n,1);
F.TargetStepM = CFG.TargetStepM * ones(n,1);
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
F.is_mixed = F.patch_purity < 0.95;
F.is_strong_mixed = F.patch_purity < 0.75;
F.purity_bin = purity_bin(F.patch_purity);
D = distance_to_material_boundary(sim.cs_map, cfg.dx);
F.distance_to_boundary_mm = D(sub2ind(size(D), F.cz, F.cx));
F.distance_bin = distance_bin(F.distance_to_boundary_mm, CFG.DistanceEdgesMm);
F.q_theory_prior = repmat(theory_q_for_regime(cfg, M, R.field_regime), n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
end

function F = ensure_predictor_columns(F, features, CFG)
for f = string(features(:))'
    if ~ismember(f, string(F.Properties.VariableNames))
        switch f
            case "dx", F.dx = F.dx;
            case "dz", F.dz = F.dz;
            case "f0", F.f0 = F.f0;
            case "M", F.M = F.M;
            case "REQ_StepX"
                F.REQ_StepX = max(1, round(CFG.TargetStepM ./ F.dx));
            case "REQ_StepZ"
                F.REQ_StepZ = max(1, round(CFG.TargetStepM ./ F.dz));
            case "TargetStepM"
                F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
            otherwise
                error('Required Test38 predictor missing in dense table: %s', f);
        end
    end
end
end

function T_out = apply_models(T, MODELS, base_features, dataset_name, keep_models)
P = predict_composition(MODELS.composition, T, base_features);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;
model_names = MODELS.model_names(ismember(MODELS.model_names, keep_models));
parts = cell(numel(model_names),1);
for mi = 1:numel(model_names)
    name = model_names(mi);
    keep = intersect(string(T_aug.Properties.VariableNames), ...
        ["dataset","condition_key","case_id","case_family","seen_status", ...
        "field_regime","field_regime_ood","f0","M","dx","dz","map_iz","map_ix", ...
        "cx","cz","x_center_m","z_center_m","true_SWS","k_true","patch_purity", ...
        "purity_bin","is_mixed","is_strong_mixed","distance_to_boundary_mm", ...
        "distance_bin","q_oracle","q_theory_prior","sws_theory", ...
        "predicted_patch_purity","p_mixed","p_strong_mixed"], 'stable');
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
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        otherwise
            continue;
    end
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T_aug.q_oracle;
    R.sws_signed_error_pct = 100*(sws_pred - T_aug.true_SWS)./T_aug.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{mi} = R;
end
parts = parts(~cellfun(@isempty, parts));
T_out = vertcat(parts{:});
end

function y = predict_q_model(M, T)
X = T(:, cellstr(M.features));
y = clamp01(predict(M.model, X));
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

%% Labels, summaries, plots

function T = add_roi_labels(T, CFG)
T.material_region = repmat("unknown", height(T), 1);
T.roi_region = repmat("other", height(T), 1);
[G, keys] = findgroups(T.condition_key);
for gi = 1:numel(keys)
    idx = G == gi;
    X = T(idx,:);
    fam = string(X.case_family(1));
    cs = X.true_SWS;
    cmin = min(cs, [], 'omitnan'); cmax = max(cs, [], 'omitnan'); crange = cmax - cmin;
    mat = repmat("unknown", height(X), 1);
    if fam == "homogeneous" || crange < 0.05
        mat(:) = "homogeneous";
    elseif contains(fam, ["inclusion","ellipse","two_inclusions"])
        mat(:) = "soft_background"; mat(cs >= cmin + 0.5*crange) = "hard_inclusion";
    elseif fam == "three_material"
        mat(:) = "soft";
        mat(cs >= cmin + crange/3 & cs < cmin + 2*crange/3) = "mid";
        mat(cs >= cmin + 2*crange/3) = "hard";
    else
        mat(:) = "soft"; mat(cs >= cmin + 0.5*crange) = "hard";
    end
    d = X.distance_to_boundary_mm;
    roi = repmat("other", height(X), 1);
    x0 = median(X.x_center_m,'omitnan'); z0 = median(X.z_center_m,'omitnan');
    half = 0.5*CFG.CenterRoiSizeMm/1e3;
    center_roi = abs(X.x_center_m-x0) <= half & abs(X.z_center_m-z0) <= half;
    if fam == "homogeneous"
        roi(center_roi) = "homogeneous_center_8mm";
    else
        roi(d <= 0.5) = "interface_0_0p5mm";
        roi(d > 0.5 & d <= 1) = "interface_0p5_1mm";
        roi(d > 1 & d <= 2) = "interface_1_2mm";
        roi(d > 2 & d <= 4) = "interface_2_4mm";
        roi(ismember(mat,["soft","soft_background"]) & d > CFG.SoftCoreMinDistanceMm) = "soft_core_gt8mm";
        roi(ismember(mat,["hard","hard_inclusion"]) & d > CFG.HardCoreMinDistanceMm) = "hard_core_gt4mm";
        roi(mat == "mid" & d > CFG.HardCoreMinDistanceMm) = "mid_core_gt4mm";
    end
    T.material_region(idx) = mat;
    T.roi_region(idx) = roi;
end
end

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
T.x_um = round(1e6*T.x_center_m);
T.z_um = round(1e6*T.z_center_m);
T.align_key = T.base_condition + "__xum" + string(T.x_um) + "__zum" + string(T.z_um);
end

function W = build_m2_m3_wide(T0, CFG)
T2all = sortrows(T0(abs(T0.M-2)<1e-8,:), {'base_condition','x_um','z_um'});
T3all = sortrows(T0(abs(T0.M-3)<1e-8,:), {'base_condition','x_um','z_um'});
conds = intersect(unique(string(T2all.base_condition), 'stable'), ...
    unique(string(T3all.base_condition), 'stable'), 'stable');

parts = cell(numel(conds), 1);
tol_um = max(1, round(1e6*CFG.AlignToleranceM));
for ci = 1:numel(conds)
    T2c = T2all(string(T2all.base_condition)==conds(ci), :);
    T3c = T3all(string(T3all.base_condition)==conds(ci), :);
    if isempty(T2c) || isempty(T3c), continue; end

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
        mean(X.sws_pred,'omitnan'), std_omitnan(X.sws_pred), ...
        mean(X.true_SWS,'omitnan'), mean(abs(X.q_error),'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
        'median_signed_error_pct','high_error10_pct','high_error20_pct', ...
        'underestimate_pct','mean_pred_sws','std_pred_sws','mean_true_sws', ...
        'mean_abs_q_error'});
end
S = [groups vertcat(rows{:})];
end

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
    apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, false);
    apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, true);
    apply_benefit_harm_gate(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, false);
    apply_benefit_harm_gate(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, true);
    make_oracle(W, sws_resid, delta)
    ];
end

function R = apply_interface_penalized_selector(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, spatial_gate)
label = string(predict(MODEL.interface_penalized_selector, W(:, cellstr(CFG.SelectorFeatureVars))));
hard_like = operational_hard_pure_like(W, CFG);
interface_like = operational_interface_like(W, CFG);
risk_ok = W.M2_risk_high20 < CFG.PenalizedSelectorRiskThreshold | hard_like;
not_extreme_mixed = W.M2_p_mixed < 0.75 | hard_like;

action = zeros(height(W), 1);
use_m3 = label == "switch_M3" & m3_stable & risk_ok & not_extreme_mixed;
use_resid = label == "apply_residual" & resid_stable & risk_ok & not_extreme_mixed;
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

function R = apply_benefit_harm_gate(W, MODEL, sws_resid, delta, m3_stable, resid_stable, CFG, conservative)
[C, M] = build_candidate_feature_table(W, MODEL.risk, sws_resid, delta, CFG);
X = C(:, cellstr(CFG.CandidateFeatureVars));
expected_gain = predict(MODEL.benefit_harm.gain_model, X);
[~, harm_score] = predict(MODEL.benefit_harm.harm_model, X);
harm_prob = positive_score(MODEL.benefit_harm.harm_model, harm_score);

n = height(W);
stable = false(2*n,1);
stable(M.candidate_code == 1) = m3_stable;
stable(M.candidate_code == 2) = resid_stable;

if conservative
    min_gain = CFG.ConservativeMinGainPct;
    max_harm = CFG.ConservativeMaxHarmProb;
    harm_penalty = CFG.ConservativeHarmPenaltyPct;
    max_candidate_risk = CFG.ConservativeMaxCandidateRisk;
    strategy_name = "conservative_benefit_harm_gate";
else
    min_gain = CFG.BenefitMinGainPct;
    max_harm = CFG.BenefitMaxHarmProb;
    harm_penalty = CFG.BenefitHarmPenaltyPct;
    max_candidate_risk = CFG.BenefitMaxCandidateRisk;
    strategy_name = "learned_benefit_harm_gate";
end

score = expected_gain - harm_penalty * harm_prob;
accept = stable & expected_gain >= min_gain & harm_prob <= max_harm & ...
    M.candidate_risk_high20 <= max_candidate_risk & score > 0;

interface_like = operational_interface_like(W, CFG);
hard_like = operational_hard_pure_like(W, CFG);
row_interface = interface_like(M.row_idx) & ~hard_like(M.row_idx);
block_interface_m3 = row_interface & M.candidate_code == 1;
block_high_risk_interface_residual = row_interface & M.candidate_code == 2 & ...
    W.M2_risk_high20(M.row_idx) >= CFG.InterfaceRiskThreshold;
accept(block_interface_m3 | block_high_risk_interface_residual) = false;
accept(row_interface) = accept(row_interface) & ...
    expected_gain(row_interface) >= (min_gain + CFG.InterfaceExtraGainMarginPct);

score(~accept) = -Inf;
score_mat = reshape(score, n, 2);
[best_score, best_choice] = max(score_mat, [], 2);
best_choice(~isfinite(best_score)) = 0;

use_m3 = best_choice == 1;
use_resid = best_choice == 2;
sws = W.M2_sws;
q = W.M2_q;
Msrc = 2*ones(n,1);
mechanism = repmat("none", n, 1);
sws(use_m3) = W.M3_sws(use_m3);
q(use_m3) = W.M3_q(use_m3);
Msrc(use_m3) = 3;
mechanism(use_m3) = "benefit_harm_M3";
sws(use_resid) = sws_resid(use_resid);
mechanism(use_resid) = "benefit_harm_residual";
if conservative
    mechanism(use_m3) = "conservative_benefit_harm_M3";
    mechanism(use_resid) = "conservative_benefit_harm_residual";
end
was = use_m3 | use_resid;
risk_after = risk_for_candidate(MODEL.risk,W,sws,q,Msrc,CFG);
R = make_strategy(W, strategy_name, sws, q, Msrc, delta, was, mechanism, ...
    W.M2_risk_high20, risk_after, false(n,1));

gain_mat = reshape(expected_gain, n, 2);
harm_mat = reshape(harm_prob, n, 2);
score_all_mat = reshape(expected_gain - harm_penalty * harm_prob, n, 2);
chosen_idx = sub2ind([n 2], (1:n)', max(best_choice,1));
R.expected_gain_pct = NaN(n,1);
R.predicted_harm_prob = NaN(n,1);
R.accept_score = NaN(n,1);
valid_choice = best_choice > 0;
R.expected_gain_pct(valid_choice) = gain_mat(chosen_idx(valid_choice));
R.predicted_harm_prob(valid_choice) = harm_mat(chosen_idx(valid_choice));
R.accept_score(valid_choice) = score_all_mat(chosen_idx(valid_choice));
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
R.expected_gain_pct = NaN(height(W),1);
R.predicted_harm_prob = NaN(height(W),1);
R.accept_score = NaN(height(W),1);
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
if r == 0 || CFG.SpatialGateMinSupport <= 1, return; end
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

function [C, M] = build_candidate_feature_table(W, risk_model, sws_resid, delta, CFG)
n = height(W);
row_idx = [(1:n)'; (1:n)'];
candidate_code = [ones(n,1); 2*ones(n,1)];
candidate_source_M = [3*ones(n,1); 2*ones(n,1)];
candidate_sws = [W.M3_sws; sws_resid];
candidate_q = [W.M3_q; W.M2_q];
candidate_risk_high20 = [W.M3_risk_high20; ...
    risk_for_candidate(risk_model, W, sws_resid, W.M2_q, 2, CFG)];

candidate_sws = clamp(candidate_sws, 0.2, 20);
candidate_q = clamp(candidate_q, 0, 1);
candidate_k = 2*pi*W.f0(row_idx) ./ candidate_sws;
candidate_log_sws = log(candidate_sws);
candidate_logk = log(candidate_k);
candidate_M_eff = (candidate_source_M .* CFG.CsGuess ./ W.f0(row_idx)) ./ ...
    (candidate_sws ./ W.f0(row_idx));
candidate_delta_sws = candidate_sws - W.M2_sws(row_idx);
candidate_rel_delta_sws = candidate_delta_sws ./ max(W.M2_sws(row_idx), eps);
candidate_abs_rel_delta_sws = abs(candidate_rel_delta_sws);

C = table(candidate_code, candidate_q, candidate_sws, candidate_log_sws, ...
    candidate_k, candidate_logk, candidate_M_eff, candidate_source_M, ...
    candidate_risk_high20, candidate_delta_sws, candidate_rel_delta_sws, ...
    candidate_abs_rel_delta_sws, W.M2_q(row_idx), W.M2_sws(row_idx), ...
    W.M2_M_eff(row_idx), W.M2_risk_high20(row_idx), ...
    W.M2_predicted_patch_purity(row_idx), W.M2_p_mixed(row_idx), ...
    W.M2_p_strong_mixed(row_idx), W.M3_sws(row_idx), W.M3_M_eff(row_idx), ...
    W.sws_M3_minus_M2(row_idx), W.rel_sws_M3_minus_M2(row_idx), ...
    W.f0(row_idx), W.dx(row_idx), W.dz(row_idx), W.field_code(row_idx), ...
    W.case_family_code(row_idx), ...
    'VariableNames', cellstr(CFG.CandidateFeatureVars));

M = table(row_idx, candidate_code, candidate_source_M, candidate_sws, ...
    candidate_q, candidate_risk_high20, delta(row_idx), ...
    abs(100*(candidate_sws - W.true_SWS(row_idx))./W.true_SWS(row_idx)), ...
    'VariableNames', {'row_idx','candidate_code','candidate_source_M', ...
    'candidate_sws','candidate_q','candidate_risk_high20','delta_logk', ...
    'candidate_abs_error_pct'});
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

function S = summarize_strategy_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T), S = table(); return; end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    corrected_rate = mean(X.was_corrected,'omitnan');
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), 100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*corrected_rate, mean(abs(X.correction_pct),'omitnan'), ...
        100*mean(X.improved_vs_M2 & X.was_corrected,'omitnan')/max(corrected_rate,eps), ...
        100*mean(X.harmed_vs_M2 & X.was_corrected,'omitnan')/max(corrected_rate,eps), ...
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
S = summarize_strategy_predictions(T, ["strategy_name","correction_mechanism"]);
end

%% Test 33 reference comparison

function T33 = load_and_summarize_test33_reference(root_dir, CFG, OUT)
T33 = struct();
source = resolve_test33_table(root_dir, CFG.Test33Source);
if exist(source,'file') ~= 2
    warning('Test48:MissingTest33', 'Missing Test33 reference table: %s', source);
    return;
end
fprintf('Loading Test 33 reference predictions (%s): %s\n', CFG.Test33Source, source);
P = readtable(source);
if ismember("split_role", string(P.Properties.VariableNames))
    P = P(string(P.split_role) == "test", :);
end
Tlong = test33_to_long(P, CFG.Test33Strategies);
if isempty(Tlong)
    warning('Test48:NoTest33Strategies', 'No requested Test33 strategies were found.');
    return;
end
T33.overall = summarize_predictions(Tlong, "model_name");
T33.by_geometry = summarize_predictions(Tlong, ["model_name","case_id"]);
T33.by_frequency = summarize_predictions(Tlong, ["model_name","case_id","f0"]);
T33.by_M = summarize_predictions(Tlong, ["model_name","M"]);
T33.by_roi = summarize_predictions(Tlong(Tlong.roi_region ~= "other",:), ...
    ["model_name","case_id","roi_region"]);
T33.by_distance = summarize_predictions(Tlong, ["model_name","case_id","distance_bin"]);
T33.by_side = summarize_predictions(Tlong(ismember(Tlong.material_region,["soft","hard"]),:), ...
    ["model_name","case_id","material_region"]);
T33.notes = table("test33_reference", string(source), height(Tlong), ...
    "These summaries use cached Test 33 test predictions. They are not regenerated on the dense Test 48 fields because Test 33 depends on Test 31/32 intermediate operational maps.", ...
    'VariableNames', {'source_name','source_file','N_long_rows','note'});
writetable(T33.overall, fullfile(OUT.table_dir,'test48_test33_reference_summary_overall.csv'));
writetable(T33.by_geometry, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_geometry.csv'));
writetable(T33.by_frequency, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_frequency.csv'));
writetable(T33.by_M, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_M.csv'));
writetable(T33.by_roi, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_roi.csv'));
writetable(T33.by_distance, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_distance.csv'));
writetable(T33.by_side, fullfile(OUT.table_dir,'test48_test33_reference_summary_by_side.csv'));
writetable(T33.notes, fullfile(OUT.table_dir,'test48_test33_reference_notes.csv'));
fprintf('Test33 reference strategies summarized: %s\n', strjoin(unique(Tlong.model_name,'stable'), ', '));
end

function file = resolve_test33_table(root_dir, source)
base = fullfile(root_dir, 'outputs', 'test_33_mixedness_aware_q_correction');
switch string(source)
    case "quick"
        file = fullfile(base, 'quick', 'tables', 'test33_patch_level_predictions.csv');
    otherwise
        file = fullfile(base, 'tables', 'test33_patch_level_predictions.csv');
end
end

function Tlong = test33_to_long(P, strategies)
parts = {};
for s = string(strategies(:))'
    sws_var = "sws_" + s;
    if s == "mixedness_logk_corrected"
        sws_var = "sws_mixedness_logk_corrected";
    elseif s == "mixedness_q_candidate_selector"
        sws_var = "sws_mixedness_q_candidate_selector";
    end
    if ~ismember(sws_var, string(P.Properties.VariableNames)), continue; end
    R = table();
    R.dataset = repmat("test33_cached_reference", height(P), 1);
    R.condition_key = string(P.condition_key);
    R.case_id = string(P.geometry);
    if ismember("geometry_type", string(P.Properties.VariableNames))
        R.case_family = string(P.geometry_type);
    else
        R.case_family = string(P.geometry);
    end
    R.seen_status = repmat("test33_original_domain", height(P), 1);
    R.field_regime = string(P.field_regime);
    R.field_regime_ood = string(P.field_regime);
    R.f0 = P.f0;
    R.M = P.M;
    R.dx = P.dx;
    R.dz = P.dz;
    R.map_iz = P.map_iz;
    R.map_ix = P.map_ix;
    R.x_center_m = P.x;
    R.z_center_m = P.z;
    R.true_SWS = P.true_SWS;
    R.patch_purity = P.patch_purity;
    R.purity_bin = string(P.purity_bin);
    R.distance_to_boundary_mm = P.distance_to_interface_mm;
    R.distance_bin = string(P.distance_bin);
    R.material_region = string(P.material_side);
    R.roi_region = roi_from_test33(P);
    R.model_name = repmat("T33_" + s, height(P), 1);
    R.q_pred = nan(height(P),1);
    R.sws_pred = P.(sws_var);
    R.q_error = nan(height(P),1);
    R.sws_signed_error_pct = 100*(R.sws_pred - R.true_SWS)./R.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{end+1,1} = R; %#ok<AGROW>
end
if isempty(parts), Tlong = table(); else, Tlong = vertcat(parts{:}); end
end

function roi = roi_from_test33(P)
roi = repmat("other", height(P), 1);
side = string(P.material_side);
d = P.distance_to_interface_mm;
geo = string(P.geometry_type);
hom = geo == "homogeneous" | contains(string(P.geometry), "homogeneous");
roi(hom) = "homogeneous_center_8mm";
roi(~hom & d <= 0.5) = "interface_0_0p5mm";
roi(~hom & d > 0.5 & d <= 1) = "interface_0p5_1mm";
roi(~hom & d > 1 & d <= 2) = "interface_1_2mm";
roi(~hom & d > 2 & d <= 4) = "interface_2_4mm";
roi(~hom & side == "soft" & d > 8) = "soft_core_gt8mm";
roi(~hom & side == "hard" & d > 4) = "hard_core_gt4mm";
end

function plot_test33_reference_comparison(T48_overall, T33_overall, T33_roi, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 15]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
A = T48_overall(T48_overall.model_name ~= "theory_discrete",:);
A.source_group = repmat("Test48 dense, matched REQ", height(A), 1);
B = T33_overall;
B.source_group = repmat("Test33 cached test domain", height(B), 1);
names = [pretty_model(A.model_name); pretty_model(B.model_name)];
vals = [A.MAPE_pct; B.MAPE_pct];
grp = [A.source_group; B.source_group];
barh(ax, categorical(names), vals);
xlabel(ax,'MAPE (%)'); title(ax,'Overall comparison, not same dataset','FontWeight','normal');
grid(ax,'on');
for i = 1:numel(vals)
    text(ax, vals(i), i, "  " + grp(i), 'Interpreter','none', 'FontSize',7);
end
ax = nexttile(tl);
R = T33_roi(ismember(T33_roi.roi_region, ...
    ["soft_core_gt8mm","hard_core_gt4mm","interface_0_0p5mm","interface_0p5_1mm","interface_1_2mm"]),:);
if ~isempty(R)
    R = aggregate_for_plot(R, ["model_name","roi_region"]);
    labels = pretty_model(R.model_name) + " | " + pretty_region(R.roi_region);
    [~,ord] = sort(R.MAPE_pct,'descend');
    ord = ord(1:min(numel(ord),20));
    barh(ax, categorical(labels(ord)), R.MAPE_pct(ord));
    xlabel(ax,'MAPE (%)'); title(ax,'Test33 reference ROI errors','FontWeight','normal'); grid(ax,'on');
else
    text(ax,0.1,0.5,'No Test33 ROI rows available');
end
export_fig(fig, fullfile(OUT.figure_dir,'test48_test33_reference_comparison.png'));
end

function R = summarize_roi_table(T)
R = summarize_predictions(T(T.roi_region ~= "other",:), ...
    ["model_name","condition_key","case_id","case_family","seen_status", ...
    "field_regime_ood","f0","M","roi_region","material_region"]);
R.roi_size_mm = repmat("", height(R), 1);
R.roi_size_mm(R.roi_region=="homogeneous_center_8mm") = "8 x 8 mm";
R.roi_size_mm(R.roi_region=="soft_core_gt8mm") = "material side, distance >8 mm";
R.roi_size_mm(R.roi_region=="hard_core_gt4mm" | R.roi_region=="mid_core_gt4mm") = ...
    "material side, distance >4 mm";
R.roi_size_mm(contains(R.roi_region,"interface")) = "distance band";
end

function plot_dense_maps(T, key, OUT, CFG)
base = T(T.model_name == CFG.ModelsToPlot(1),:);
if isempty(base), base = T(T.model_name == T.model_name(1),:); end
[true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
pur_map = rows_to_grid(base, base.patch_purity, nz, nx);
ppur_map = rows_to_grid(base, base.predicted_patch_purity, nz, nx);
pmix_map = rows_to_grid(base, base.p_mixed, nz, nx);
dist_map = rows_to_grid(base, base.distance_to_boundary_mm, nz, nx);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 25]);
tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), true_map, 'True SWS (m/s)');
plot_map(nexttile(tl), pur_map, 'True patch purity');
plot_map(nexttile(tl), ppur_map, 'Predicted patch purity');
plot_map(nexttile(tl), pmix_map, 'Predicted mixedness probability');
for m = CFG.ModelsToPlot
    X = T(T.model_name == m,:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), pretty_model(m) + " SWS");
end
plot_map(nexttile(tl), rows_to_grid(T(T.model_name=="theory_discrete",:), ...
    T.sws_pred(T.model_name=="theory_discrete"), nz, nx), 'Theory discrete SWS');
for m = CFG.ModelsToPlot
    X = T(T.model_name == m,:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, nz, nx), pretty_model(m) + " abs error (%)");
end
for m = CFG.ModelsToPlot
    X = T(T.model_name == m,:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_signed_error_pct, nz, nx), pretty_model(m) + " signed error (%)");
end
plot_map(nexttile(tl), rows_to_grid(base, base.q_pred, nz, nx), 'Primary q prediction');
plot_map(nexttile(tl), dist_map, 'Distance to interface (mm)');
title(tl, key, 'Interpreter','none');
out_dir = fullfile(OUT.map_dir, sanitize(base.case_id(1)), sanitize(base.field_regime_ood(1)));
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
export_fig(fig, fullfile(out_dir, "test48_map__" + sanitize(key) + ".png"));
end

function plot_roi_overlay_pair(T, key, OUT, CFG)
base = T(T.model_name == "q_spectrum_plus_composition",:);
if isempty(base), base = T(T.model_name == T.model_name(1),:); end
[true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
pred_map = rows_to_grid(base, base.sws_pred, nz, nx);
roi_map = rows_to_grid(base, double(categorical(base.roi_region)), nz, nx);
dist_map = rows_to_grid(base, base.distance_to_boundary_mm, nz, nx);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 12]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); imagesc(ax,true_map); axis(ax,'image'); colorbar(ax);
title(ax,'True SWS + ROI overlay','FontWeight','normal'); hold(ax,'on'); draw_roi_overlay(ax, base, dist_map, CFG);
ax = nexttile(tl); imagesc(ax,pred_map); axis(ax,'image'); colorbar(ax);
title(ax,'Primary predicted SWS + ROI overlay','FontWeight','normal'); hold(ax,'on'); draw_roi_overlay(ax, base, dist_map, CFG);
ax = nexttile(tl); imagesc(ax,roi_map); axis(ax,'image'); colorbar(ax);
title(ax,'ROI label map','FontWeight','normal');
title(tl, key, 'Interpreter','none');
out_dir = fullfile(OUT.roi_dir, sanitize(base.case_id(1)), sanitize(base.field_regime_ood(1)));
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
export_fig(fig, fullfile(out_dir, "test62_roi_overlay__" + sanitize(key) + ".png"));
end

function draw_roi_overlay(ax, X, dist_map, CFG)
if string(X.case_family(1)) == "homogeneous"
    ix0 = median(X.map_ix,'omitnan'); iz0 = median(X.map_iz,'omitnan');
    dx = median(diff(unique(sort(X.x_center_m))),'omitnan');
    dz = median(diff(unique(sort(X.z_center_m))),'omitnan');
    if ~isfinite(dx) || dx <= 0, dx = median(X.dx,'omitnan'); end
    if ~isfinite(dz) || dz <= 0, dz = median(X.dz,'omitnan'); end
    w = CFG.CenterRoiSizeMm/1e3/dx; h = CFG.CenterRoiSizeMm/1e3/dz;
    rectangle(ax,'Position',[ix0-w/2 iz0-h/2 w h],'EdgeColor','w', ...
        'LineWidth',1.8,'LineStyle','--');
    text(ax, ix0, iz0, '8 x 8 mm', 'Color','w','FontWeight','bold', ...
        'HorizontalAlignment','center','BackgroundColor',[0 0 0 0.35]);
else
    if any((dist_map(:) <= 0.5) & isfinite(dist_map(:)))
        contour(ax, dist_map <= 0.5, [1 1], 'w', 'LineWidth',1.4);
    end
    if any((dist_map(:) > 4) & isfinite(dist_map(:)))
        contour(ax, dist_map > 4, [1 1], 'k', 'LineWidth',1.0);
    end
    text(ax, 3, 5, 'white: 0-0.5 mm interface | black: core >4 mm', ...
        'Color','w','FontWeight','bold','BackgroundColor',[0 0 0 0.35]);
end
end

function plot_summary_figures(T_overall, T_case, T_freq, T_M, T_roi, OUT)
main = T_overall(T_overall.model_name ~= "theory_discrete",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 16]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
main = sortrows(main,'MAPE_pct','ascend');
barh(ax, categorical(pretty_model(main.model_name)), main.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Dense representative model ranking','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for m = unique(T_freq.model_name,'stable')'
    if m == "theory_discrete", continue; end
    X = sortrows(T_freq(T_freq.model_name==m,:), 'f0');
    plot(ax, X.f0, X.MAPE_pct, '-o', 'DisplayName', pretty_model(m));
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'MAPE vs frequency','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','best','Interpreter','none');
ax = nexttile(tl); hold(ax,'on');
for m = unique(T_M.model_name,'stable')'
    if m == "theory_discrete", continue; end
    X = sortrows(T_M(T_M.model_name==m,:), 'M');
    plot(ax, X.M, X.MAPE_pct, '-o', 'DisplayName', pretty_model(m));
end
xlabel(ax,'REQ M'); ylabel(ax,'MAPE (%)'); title(ax,'MAPE vs M','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','best','Interpreter','none');
ax = nexttile(tl);
P = T_roi(T_roi.model_name=="q_spectrum_plus_composition",:);
P = aggregate_for_plot(P, "roi_region");
P = sortrows(P,'MAPE_pct','descend');
barh(ax, categorical(pretty_region(P.roi_region)), P.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Primary model error by ROI','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test48_summary_diagnostics.png'));
end

function plot_test62_summary_figures(T_overall, T_case, T_freq, T_roi, T_harm, OUT)
deployable = T_overall(~contains(string(T_overall.strategy_name),"oracle"),:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl);
A = sortrows(deployable,'MAPE_pct','ascend');
barh(ax, categorical(pretty_strategy(A.strategy_name)), A.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Dense strategy ranking','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl);
A = sortrows(deployable,'high_error20_pct','ascend');
barh(ax, categorical(pretty_strategy(A.strategy_name)), A.high_error20_pct);
xlabel(ax,'Pixels with error >20% (%)'); title(ax,'High-error ranking','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl); hold(ax,'on');
main_strategies = ["fixed_M2","fixed_M3","interface_penalized_selector", ...
    "learned_benefit_harm_gate","conservative_benefit_harm_gate"];
for s = main_strategies
    X = sortrows(T_freq(T_freq.strategy_name==s,:), 'f0');
    if isempty(X), continue; end
    plot(ax, X.f0, X.MAPE_pct, '-o', 'DisplayName', pretty_strategy(s));
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'MAPE vs frequency','FontWeight','normal');
grid(ax,'on'); legend(ax,'Location','bestoutside','Interpreter','none');

ax = nexttile(tl);
R = T_roi(ismember(T_roi.roi_region, ["homogeneous_center_8mm","soft_core_gt8mm", ...
    "hard_core_gt4mm","interface_0_0p5mm","interface_0p5_1mm"]),:);
R = R(ismember(R.strategy_name, main_strategies),:);
if ~isempty(R)
    R = aggregate_strategy_for_plot(R, ["strategy_name","roi_region"]);
    labels = pretty_strategy(R.strategy_name) + " | " + pretty_region(R.roi_region);
    [~,ord] = sort(R.MAPE_pct,'descend');
    ord = ord(1:min(24,numel(ord)));
    barh(ax, categorical(labels(ord)), R.MAPE_pct(ord));
    xlabel(ax,'MAPE (%)'); title(ax,'ROI/core/interface errors','FontWeight','normal'); grid(ax,'on');
else
    text(ax,0.1,0.5,'No ROI rows available');
end

ax = nexttile(tl);
H = T_harm(~contains(string(T_harm.strategy_name),"oracle"),:);
H = H(H.correction_mechanism ~= "none",:);
if ~isempty(H)
    H = sortrows(H,'harm_rate_corrected_pct','ascend');
    labels = pretty_strategy(H.strategy_name) + " | " + string(H.correction_mechanism);
    barh(ax, categorical(labels), H.harm_rate_corrected_pct);
    xlabel(ax,'Harmed corrected pixels (%)'); title(ax,'Correction harm among corrected pixels','FontWeight','normal'); grid(ax,'on');
else
    text(ax,0.1,0.5,'No corrections applied');
end

ax = nexttile(tl);
C = T_case(ismember(T_case.strategy_name, main_strategies),:);
if ~isempty(C)
    [~,ord] = sort(C.MAPE_pct,'descend');
    ord = ord(1:min(20,numel(ord)));
    labels = pretty_strategy(C.strategy_name(ord)) + " | " + string(C.case_id(ord));
    barh(ax, categorical(labels), C.MAPE_pct(ord));
    xlabel(ax,'MAPE (%)'); title(ax,'Worst dense cases','FontWeight','normal'); grid(ax,'on');
else
    text(ax,0.1,0.5,'No case summaries available');
end
export_fig(fig, fullfile(OUT.figure_dir,'test62_summary_diagnostics.png'));
end

function plot_test62_maps(T, Tlong, CND, OUT, CFG)
base = T(T.strategy_name=="fixed_M2",:);
[true_map,nz,nx] = rows_to_grid_strategy(base, base.true_SWS);
ppur_map = rows_to_grid_strategy(base, base.M2_predicted_patch_purity, nz, nx);
pmix_map = rows_to_grid_strategy(base, base.M2_p_mixed, nz, nx);
risk_before = rows_to_grid_strategy(base, base.risk_before, nz, nx);

plot_strategies = CFG.StrategiesToPlot(ismember(CFG.StrategiesToPlot, unique(T.strategy_name)));
fig = figure('Color','w','Units','centimeters','Position',[1 1 36 28]);
tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
plot_map_with_label(nexttile(tl), true_map, 'True SWS', 'SWS (m/s)');
plot_map_with_label(nexttile(tl), ppur_map, 'Predicted patch purity', 'Probability');
plot_map_with_label(nexttile(tl), pmix_map, 'Predicted mixedness', 'Probability');
plot_map_with_label(nexttile(tl), risk_before, 'Risk before correction', 'High-error probability');

for s = plot_strategies
    X = T(T.strategy_name==s,:);
    plot_map_with_label(nexttile(tl), rows_to_grid_strategy(X, X.sws_pred_strategy, nz, nx), ...
        pretty_strategy(s) + " SWS", 'SWS (m/s)');
end
for s = plot_strategies(1:min(3,numel(plot_strategies)))
    X = T(T.strategy_name==s,:);
    plot_map_with_label(nexttile(tl), rows_to_grid_strategy(X, X.sws_signed_error_pct, nz, nx), ...
        pretty_strategy(s) + " signed error", 'Signed error (%)');
end

sel = first_existing_strategy(T, ["learned_benefit_harm_gate","conservative_benefit_harm_gate", ...
    "interface_penalized_selector"]);
X = T(T.strategy_name==sel,:);
plot_map_with_label(nexttile(tl), rows_to_grid_strategy(X, X.was_corrected, nz, nx), ...
    pretty_strategy(sel) + " correction mask", 'Corrected (0/1)');
plot_map_with_label(nexttile(tl), rows_to_grid_strategy(X, X.correction_pct, nz, nx), ...
    pretty_strategy(sel) + " correction", 'SWS change (%)');
plot_map_with_label(nexttile(tl), rows_to_grid_strategy(X, X.risk_after, nz, nx), ...
    pretty_strategy(sel) + " risk after", 'High-error probability');

T2 = Tlong(abs(Tlong.M-2)<1e-8,:);
if ~isempty(T2)
    plot_map_with_label(nexttile(tl), rows_to_grid(T2, T2.q_pred), ...
        'M=2 q prediction', 'REQ quantile q');
else
    axis(nexttile(tl),'off');
end

title(tl, sprintf('Case: %s, f=%g Hz, regime=%s', ...
    CND.case.case_id, CND.f0, CND.regime.regime_id), 'Interpreter','none');
out_dir = fullfile(OUT.map_dir, sanitize(CND.case.case_id), sanitize(CND.regime.regime_id));
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
key = sprintf('%s__f%g__%s__dx%gum', CND.case.case_id, CND.f0, CND.regime.regime_id, round(1e6*median(T.dx)));
export_fig(fig, fullfile(out_dir, "test62_map__" + sanitize(key) + ".png"));
end

function s = first_existing_strategy(T, candidates)
s = candidates(1);
for c = candidates
    if any(T.strategy_name==c), s = c; return; end
end
s = T.strategy_name(1);
end

function [Z,nz,nx] = rows_to_grid_strategy(T, values, nz, nx)
if nargin < 3
    nx = max(T.map_ix); nz = max(T.map_iz);
end
Z = nan(nz,nx);
if isempty(T), return; end
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
Z(idx) = values;
end

function S = aggregate_strategy_for_plot(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
S = groups;
S.N = splitapply(@(x) sum(x,'omitnan'), T.N, G);
S.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), T.MAPE_pct, T.N, G);
if ismember("high_error20_pct", string(T.Properties.VariableNames))
    S.high_error20_pct = splitapply(@(x,w) weighted_mean(x,w), T.high_error20_pct, T.N, G);
end
if ismember("harm_rate_corrected_pct", string(T.Properties.VariableNames))
    S.harm_rate_corrected_pct = splitapply(@(x,w) weighted_mean(x,w), T.harm_rate_corrected_pct, T.N, G);
end
end

%% Numerical and plotting helpers

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
if ~any(B(:)), Dmm = nan(size(Q)); else, Dmm = bwdist(B) * dx * 1e3; end
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

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
if isempty(T), return; end
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map(ax, Z, ttl)
imagesc(ax, Z, 'AlphaData', isfinite(Z));
axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
ax.Color = [0.94 0.94 0.94];
colorbar(ax); title(ax, ttl, 'Interpreter','none','FontWeight','normal');
end

function plot_map_with_label(ax, Z, ttl, label)
plot_map(ax, Z, ttl);
cb = colorbar(ax); ylabel(cb, label);
end

function S = aggregate_for_plot(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
S = groups;
S.N = splitapply(@(x) sum(x,'omitnan'), T.N, G);
S.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), T.MAPE_pct, T.N, G);
end

function y = weighted_mean(x,w)
ok = isfinite(x) & isfinite(w) & w > 0;
if ~any(ok), y = NaN; else, y = sum(x(ok).*w(ok))/sum(w(ok)); end
end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
end

function y = std_omitnan(x)
x = x(isfinite(x));
if numel(x) <= 1, y = NaN; else, y = std(x); end
end

function y = clamp01(x), y = min(max(x,0),1); end

function s = pretty_model(x)
x = string(x);
s = x;
s(x=="q_spectrum_plus_composition") = "q spectrum + composition";
s(x=="q_spectrum_only") = "q spectrum only";
s(x=="q_spectrum_plus_theory_composition") = "q spectrum + theory + composition";
s(x=="theory_discrete") = "Theory discrete";
s(x=="T33_local_baseline") = "T33 LocalOnly";
s(x=="T33_mixedness_logk_corrected") = "T33 mixedness log-k";
s(x=="T33_mixedness_q_candidate_selector") = "T33 mixedness q-selector";
s = strrep(s, "_", " ");
end

function s = pretty_region(x)
s = string(x);
s = strrep(s, "_", " ");
s = strrep(s, "0p5", "0.5");
s = strrep(s, "gt", ">");
end

function s = sanitize(x)
s = regexprep(string(x), '[^A-Za-z0-9_=-]+', '_');
s = char(s);
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
    warning('Test48:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
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

function x = env_num(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x), x = default; end
end

function tf = should_save_map(i, CFG)
tf = CFG.SaveAllMaps && (CFG.MaxMapConditions <= 0 || i <= CFG.MaxMapConditions);
end

function x = min_positive_or_default(x, default)
if ~isfinite(x) || x <= 0, x = default; end
end

function code = categorical_code(x)
[~,~,code] = unique(string(x), 'stable');
code = double(code);
end

function y = clamp(x, lo, hi)
y = min(max(x,lo),hi);
end

function s = pretty_strategy(x)
x = string(x);
s = x;
s(x=="fixed_M2") = "Fixed M=2";
s(x=="fixed_M3") = "Fixed M=3";
s(x=="test56_delta_logk_residual_on_M2") = "M-eff residual";
s(x=="M2_to_M3_hard_pure_switch") = "Hard/pure M2→M3 rule";
s(x=="risk_gated_Meff_residual") = "Risk-gated M-eff residual";
s(x=="interface_penalized_selector") = "Interface-penalized selector";
s(x=="spatial_interface_penalized_selector") = "Spatial interface-penalized selector";
s(x=="learned_benefit_harm_gate") = "Learned benefit/harm gate";
s(x=="conservative_benefit_harm_gate") = "Conservative benefit/harm gate";
s(x=="oracle_best_of_M2_M3_residual") = "Oracle best candidate";
s = strrep(s, "_", " ");
end

function print_test62_summary(T_overall, T_roi, T_harm, OUT)
fprintf('\n================ Test 62 dense Test61 summary ================\n');
cols = {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct', ...
    'corrected_pct','harm_rate_corrected_pct'};
disp(T_overall(:, cols));
deploy = T_overall(~contains(string(T_overall.strategy_name),"oracle"),:);
if ~isempty(deploy)
    deploy = sortrows(deploy,'MAPE_pct','ascend');
    fprintf('Best deployable dense strategy: %s, MAPE %.2f%%, high>20 %.2f%%.\n', ...
        pretty_strategy(deploy.strategy_name(1)), deploy.MAPE_pct(1), deploy.high_error20_pct(1));
end
hard = T_roi(T_roi.roi_region=="hard_core_gt4mm" & ...
    ~contains(string(T_roi.strategy_name),"oracle"),:);
if ~isempty(hard)
    hard = sortrows(hard,'MAPE_pct','ascend');
    fprintf('Best hard-core dense strategy: %s on %s, MAPE %.2f%%.\n', ...
        pretty_strategy(hard.strategy_name(1)), string(hard.case_id(1)), hard.MAPE_pct(1));
end
if ~isempty(T_harm)
    H = T_harm(T_harm.correction_mechanism ~= "none" & ...
        ~contains(string(T_harm.strategy_name),"oracle"),:);
    if ~isempty(H)
        H = sortrows(H,'harm_rate_corrected_pct','ascend');
        fprintf('Lowest correction harm among active mechanisms: %s / %s, harm %.2f%%.\n', ...
            pretty_strategy(H.strategy_name(1)), string(H.correction_mechanism(1)), ...
            H.harm_rate_corrected_pct(1));
    end
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

function print_summary(T_overall, T_case, T_roi, OUT)
fprintf('\n================ Test 48 dense summary ================\n');
disp(T_overall(:, {'model_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
P = T_overall(T_overall.model_name=="q_spectrum_plus_composition",:);
if ~isempty(P)
    fprintf('Primary model dense representative MAPE %.2f%%, high>20 %.2f%%.\n', ...
        P.MAPE_pct, P.high_error20_pct);
end
X = T_case(T_case.model_name=="q_spectrum_plus_composition",:);
if ~isempty(X)
    X = sortrows(X,'MAPE_pct','descend');
    fprintf('Worst dense primary case: %s (%s), MAPE %.2f%%.\n', ...
        X.case_id(1), X.seen_status(1), X.MAPE_pct(1));
end
R = T_roi(T_roi.model_name=="q_spectrum_plus_composition",:);
if ~isempty(R)
    R = aggregate_for_plot(R, "roi_region");
    R = sortrows(R,'MAPE_pct','descend');
    fprintf('Worst ROI/zone: %s, MAPE %.2f%%.\n', pretty_region(R.roi_region(1)), R.MAPE_pct(1));
end
fprintf('Outputs: %s\n', OUT.root_dir);
end
