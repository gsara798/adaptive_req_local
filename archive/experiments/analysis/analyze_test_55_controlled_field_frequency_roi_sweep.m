%% analyze_test_55_controlled_field_frequency_roi_sweep.m
% Test 55: controlled field-frequency ROI sweep for frozen Test 38/53 models.
%
% This script recalculates REQ on a small balanced design: selected
% homogeneous, bilayer, and inclusion geometries; all requested frequencies;
% all four field regimes. It never trains or modifies the saved model bundle.

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
write_config_json(CFG, fullfile(OUT.root_dir,'test55_configuration.json'));

fprintf('\nTest 55: controlled field-frequency ROI sweep\n');
fprintf('Mode: %s | target REQ step %.3f mm | save all maps: %d\n', ...
    CFG.Mode, 1e3*CFG.TargetStepM, CFG.SaveAllMaps);
fprintf('No training. Frozen model source: %s\n', CFG.ModelBundleFile);
fprintf('REQ profile: %s | Nbins=%s | oversample=%g | Nbins min=%g\n', ...
    CFG.REQ.Profile, string(CFG.REQ.Nbins), CFG.REQ.Nbins_auto_oversample, CFG.REQ.Nbins_min);
fprintf('At least one source is forced/aligned in every regime.\n');

BUNDLE = load_test38_bundle(CFG);
MODELS = BUNDLE.MODELS;
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);

if CFG.ValidateOnly
    CONDITIONS = build_condition_specs(CFG);
    CONDITIONS = CONDITIONS(1:min(1,numel(CONDITIONS)));
else
    CONDITIONS = build_condition_specs(CFG);
end
fprintf('Balanced validation conditions: %d | M values: %s\n', ...
    numel(CONDITIONS), mat2str(CFG.M));

all_parts = {};
roi_parts = {};
for ci = 1:numel(CONDITIONS)
    CND = CONDITIONS(ci);
    for mi = 1:numel(CFG.M)
        M = CFG.M(mi);
        key = sprintf('%s__f%g__%s__dx%gum__M%g', CND.case.case_id, ...
            CND.f0, CND.regime.regime_id, round(1e6*CFG.dx), M);
        timer = tic;
        [T_pred, T_base] = evaluate_dense_condition(CND, M, key, MODELS, ...
            BASE_FEATURES, CFG, OUT);
        T_pred = add_roi_labels(T_pred, CFG);
        R = summarize_roi_table(T_pred);
        all_parts{end+1,1} = remove_cell_columns(T_pred); %#ok<AGROW>
        roi_parts{end+1,1} = R; %#ok<AGROW>
        fprintf('[%d/%d] %s: %d dense patches in %.1f s.\n', ...
            (ci-1)*numel(CFG.M)+mi, numel(CONDITIONS)*numel(CFG.M), ...
            key, height(T_base), toc(timer));
        if CFG.SaveAllMaps
            plot_dense_maps(T_pred, key, OUT, CFG);
            plot_roi_overlay_pair(T_pred, key, OUT, CFG);
        end
    end
end

T_all = vertcat(all_parts{:});
T_roi = vertcat(roi_parts{:});
T_overall = summarize_predictions(T_all, "model_name");
T_by_case = summarize_predictions(T_all, ["model_name","case_id","seen_status"]);
T_by_regime = summarize_predictions(T_all, ["model_name","field_regime_ood"]);
T_by_frequency = summarize_predictions(T_all, ["model_name","f0"]);
T_by_M = summarize_predictions(T_all, ["model_name","M"]);
T_by_roi = summarize_predictions(T_all(T_all.roi_region ~= "other",:), ...
    ["model_name","case_id","roi_region"]);

T_roi_balanced = summarize_balanced_roi(T_all, CFG);

writetable(T_all, fullfile(OUT.table_dir,'test55_patch_level_results.csv'));
writetable(T_roi, fullfile(OUT.table_dir,'test55_roi_summary.csv'));
writetable(T_roi_balanced, fullfile(OUT.table_dir,'test55_roi_frequency_field_summary.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test55_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test55_summary_by_case.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test55_summary_by_regime.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir,'test55_summary_by_frequency.csv'));
writetable(T_by_M, fullfile(OUT.table_dir,'test55_summary_by_M.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir,'test55_summary_by_roi.csv'));

safe_plot(@() plot_summary_figures(T_overall, T_by_case, T_by_frequency, T_by_M, T_by_roi, OUT), ...
    'summary figures');
safe_plot(@() plot_balanced_roi_frequency_field(T_roi_balanced, OUT, CFG), ...
    'balanced ROI frequency-field figures');

save(fullfile(OUT.data_dir,'test55_compact_results.mat'), ...
    'CFG','T_overall','T_by_case','T_by_regime','T_by_frequency','T_by_M', ...
    'T_by_roi','T_roi','T_roi_balanced','-v7.3');

print_summary(T_overall, T_by_case, T_by_roi, OUT);

fprintf('\nTables: %s\nFigures: %s\nTest 55 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST55_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST55_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST55_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST55_SAVE_ALL_MAPS', false);
CFG.UseParfor = env_true('ADAPTIVE_REQ_TEST55_USE_PARFOR', false);
CFG.TargetStepM = env_num('ADAPTIVE_REQ_TEST55_TARGET_STEP_M', 0.75e-3);
CFG.ModelSource = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST55_MODEL_SOURCE'))));
if CFG.ModelSource == ""
    CFG.ModelSource = "full";
end
CFG.ModelBundleFile = resolve_model_bundle(root_dir, CFG.ModelSource);
CFG.dx = 0.2e-3;
CFG.dz = 0.2e-3;
CFG.Lx = 0.05;
CFG.Lz = 0.05;
CFG.cs_guess = 3;
CFG.M = env_num_list('ADAPTIVE_REQ_TEST55_M_LIST', [2]);
CFG.RandomSeed = 48001;
CFG.DistanceEdgesMm = [0 0.5 1 2 4 8 Inf];
CFG.CenterRoiSizeMm = 8;
CFG.SoftCoreMinDistanceMm = 8;
CFG.HardCoreMinDistanceMm = 4;
CFG.ModelsToPlot = ["q_spectrum_plus_composition","q_spectrum_only"];
CFG.ModelsToEvaluate = ["q_spectrum_only","q_spectrum_plus_composition"];
CFG.PrimaryModel = "q_spectrum_plus_composition";
CFG.REQ.Profile = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST55_REQ_PROFILE'))));
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
        error('Unknown ADAPTIVE_REQ_TEST55_REQ_PROFILE: %s', CFG.REQ.Profile);
end
CFG.REQ.Nbins_auto_oversample = env_num('ADAPTIVE_REQ_TEST55_NBINS_OVERSAMPLE', CFG.REQ.Nbins_auto_oversample);
CFG.REQ.Nbins_min = env_num('ADAPTIVE_REQ_TEST55_NBINS_MIN', CFG.REQ.Nbins_min);
CFG.REQ.SmoothSigma = 1;
if CFG.ValidateOnly
    CFG.M = 2;
end
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep');
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

function file = resolve_model_bundle(root_dir, source)
source = string(source);
if isfile(source)
    file = char(source);
    return;
end
switch lower(source)
    case {"test53","paper","paper_final","full"}
        f53 = fullfile(root_dir, 'outputs', ...
            'test_53_paper_final_clean_q_training', 'models', ...
            'test38_velocity_field_diverse_q_models.mat');
        if exist(f53,'file') == 2
            file = f53;
        else
            file = fullfile(root_dir, 'outputs', ...
                'test_38_velocity_field_diverse_q_training', 'models', ...
                'test38_velocity_field_diverse_q_models.mat');
        end
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
cases = [
    homogeneous_case("homogeneous_cs2", 2, "selected")
    homogeneous_case("homogeneous_cs3", 3, "selected")
    homogeneous_case("homogeneous_cs4", 4, "selected")
    bilayer_case("bilayer_2_3", 2, 3, 0, "selected")
    bilayer_case("bilayer_2_4", 2, 4, 0, "selected")
    inclusion_case("inclusion_2_3", 2, 3, [0.025 0.025], 0.008, "selected")
    inclusion_case("inclusion_2_4", 2, 4, [0.025 0.025], 0.008, "selected")
    ];
regimes = [
    regime("directional_2D_angle15","directional_2D",true,1,true,"ranges","fibonacci",[deg2rad(15) deg2rad(15)],[0 2*pi],[0 pi],15)
    regime("diffuse_2D_seed17","diffuse_2D",true,128,true,"ranges","random",[0 2*pi],[0 2*pi],[0 pi],17)
    regime("partial_3D_12src","partial_3D",false,12,true,"ranges","fibonacci",[0 2*pi],[0 2*pi],[0 pi],19)
    regime("diffuse_3D_seed23","diffuse_3D",false,128,true,"ranges","random",[0 2*pi],[0 2*pi],[0 pi],23)
    ];
freqs = env_num_list('ADAPTIVE_REQ_TEST55_FREQUENCIES_HZ', [200 300 400 500 600]);
if CFG.QuickMode
    cases = cases([1 3 4 6]);
    freqs = [200 500];
else
    % full keeps all selected cases and all frequencies.
end
COND = repmat(empty_condition(), 0, 1);
for ci = 1:numel(cases)
    for fi = 1:numel(freqs)
        for ri = 1:numel(regimes)
            COND(end+1,1) = condition(cases(ci), regimes(ri), freqs(fi)); %#ok<AGROW>
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
T_pred = apply_models(F, MODELS, BASE_FEATURES, "test55_controlled_field_frequency_roi_sweep", CFG.ModelsToEvaluate);
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
F.dataset = repmat("test55_controlled_field_frequency_roi_sweep", n, 1);
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

function S = summarize_balanced_roi(T, CFG)
X = T(T.model_name == CFG.PrimaryModel & T.roi_region ~= "other", :);
keep = ismember(X.roi_region, ["homogeneous_center_8mm", ...
    "soft_core_gt8mm","hard_core_gt4mm","mid_core_gt4mm", ...
    "interface_0_0p5mm","interface_0p5_1mm","interface_1_2mm"]);
X = X(keep,:);
if isempty(X)
    S = table();
    return;
end
field_type = strings(height(X),1);
for ii = 1:height(X)
    field_type(ii) = field_type_from_regime(X.field_regime_ood(ii));
end
X.field_type = field_type;
[G, groups] = findgroups(X(:, {'model_name','case_id','case_family', ...
    'field_type','field_regime_ood','f0','M','roi_region','material_region'}));
S = groups;
S.N = splitapply(@numel, X.sws_pred, G);
S.true_sws_mean = splitapply(@(x) mean(x,'omitnan'), X.true_SWS, G);
S.pred_sws_mean = splitapply(@(x) mean(x,'omitnan'), X.sws_pred, G);
S.pred_sws_std = splitapply(@std_omitnan, X.sws_pred, G);
S.MAPE_pct = splitapply(@(x) mean(x,'omitnan'), X.sws_abs_error_pct, G);
S.mean_signed_error_pct = splitapply(@(x) mean(x,'omitnan'), X.sws_signed_error_pct, G);
S.high_error20_pct = 100*splitapply(@(x) mean(x,'omitnan'), double(X.high_error20), G);
S.mean_patch_purity = splitapply(@(x) mean(x,'omitnan'), X.patch_purity, G);
end

function plot_balanced_roi_frequency_field(S, OUT, CFG)
if isempty(S), return; end
out_dir = fullfile(OUT.figure_dir, 'roi_frequency_by_field');
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
cases = unique(S.case_id, 'stable');
for ci = 1:numel(cases)
    X = S(S.case_id == cases(ci), :);
    rois = unique(X.roi_region, 'stable');
    % Keep the most interpretable regions first.
    order = ["homogeneous_center_8mm","soft_core_gt8mm","hard_core_gt4mm", ...
        "mid_core_gt4mm","interface_0_0p5mm","interface_0p5_1mm","interface_1_2mm"];
    rois = order(ismember(order, rois));
    if isempty(rois), continue; end
    fields = ["directional 2D","diffuse 2D","partial 3D","diffuse 3D"];
    fields = fields(ismember(fields, unique(X.field_type)));
    colors = lines(numel(fields));
    fig = figure('Color','w','Units','centimeters', ...
        'Position',[1 1 max(22,9*numel(rois)) 14]);
    tl = tiledlayout(fig, 2, numel(rois), 'TileSpacing','compact','Padding','compact');
    for ri = 1:numel(rois)
        Rroi = X(X.roi_region == rois(ri),:);
        ax = nexttile(tl, ri); hold(ax,'on');
        truth = groupsummary(Rroi, 'f0', 'mean', 'true_sws_mean');
        plot(ax, truth.f0, truth.mean_true_sws_mean, 'k--', 'LineWidth',1.1, ...
            'DisplayName','true SWS');
        for fi = 1:numel(fields)
            R = sortrows(Rroi(Rroi.field_type == fields(fi),:), 'f0');
            if isempty(R), continue; end
            errorbar(ax, R.f0, R.pred_sws_mean, R.pred_sws_std, '-o', ...
                'Color',colors(fi,:), 'MarkerFaceColor',colors(fi,:), ...
                'LineWidth',1.1, 'DisplayName', fields(fi));
        end
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'ROI SWS (m/s)');
        title(ax, pretty_region(rois(ri)), 'FontWeight','normal');
        grid(ax,'on');
        if ri == numel(rois)
            legend(ax,'Location','bestoutside','Interpreter','none');
        end

        ax = nexttile(tl, numel(rois)+ri); hold(ax,'on');
        for fi = 1:numel(fields)
            R = sortrows(Rroi(Rroi.field_type == fields(fi),:), 'f0');
            if isempty(R), continue; end
            plot(ax, R.f0, R.MAPE_pct, '-o', 'Color',colors(fi,:), ...
                'MarkerFaceColor',colors(fi,:), 'LineWidth',1.1, ...
                'DisplayName', fields(fi));
        end
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'ROI MAPE (%)');
        title(ax,'Error by field type','FontWeight','normal');
        grid(ax,'on');
    end
    title(tl, sprintf('%s: balanced frequency x field ROI sweep, M=%g', ...
        pretty_region(cases(ci)), CFG.M(1)), 'FontWeight','normal');
    export_fig(fig, fullfile(out_dir, "test55_roi_frequency_field__" + sanitize(cases(ci)) + ".png"));
end
end

function plot_dense_maps(T, key, OUT, CFG)
base = T(T.model_name == CFG.ModelsToPlot(1),:);
if isempty(base), base = T(T.model_name == T.model_name(1),:); end
[true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
pur_map = rows_to_grid(base, base.patch_purity, nz, nx);
ppur_map = rows_to_grid(base, base.predicted_patch_purity, nz, nx);
pmix_map = rows_to_grid(base, base.p_mixed, nz, nx);
dist_map = rows_to_grid(base, base.distance_to_boundary_mm, nz, nx);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 22]);
tl = tiledlayout(fig,4,3,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), true_map, 'True SWS (m/s)');
plot_map(nexttile(tl), pur_map, 'True patch purity');
plot_map(nexttile(tl), ppur_map, 'Predicted patch purity');
plot_map(nexttile(tl), pmix_map, 'Predicted mixedness probability');
for m = CFG.ModelsToPlot
    X = T(T.model_name == m,:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), pretty_model(m) + " SWS");
end
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
export_fig(fig, fullfile(out_dir, "test55_map__" + sanitize(key) + ".png"));
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
export_fig(fig, fullfile(out_dir, "test55_roi_overlay__" + sanitize(key) + ".png"));
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
xlabel(ax,'MAPE (%)'); title(ax,'Controlled sweep model ranking','FontWeight','normal'); grid(ax,'on');
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
export_fig(fig, fullfile(OUT.figure_dir,'test55_summary_diagnostics.png'));
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
    warning('Test55:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
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

function vals = env_num_list(name, default)
v = strtrim(string(getenv(name)));
if v == ""
    vals = default;
    return;
end
txt = regexprep(char(v), '[\[\]\(\)]', ' ');
parts = regexp(txt, '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
vals = str2double(parts);
vals = vals(isfinite(vals));
if isempty(vals), vals = default; end
end

function t = field_type_from_regime(regime)
r = lower(string(regime));
if startsWith(r, "directional")
    t = "directional 2D";
elseif startsWith(r, "diffuse_2d")
    t = "diffuse 2D";
elseif startsWith(r, "partial_3d")
    t = "partial 3D";
elseif startsWith(r, "diffuse_3d")
    t = "diffuse 3D";
else
    t = "unknown";
end
end

function print_summary(T_overall, T_case, T_roi, OUT)
fprintf('\n================ Test 55 controlled sweep summary ================\n');
disp(T_overall(:, {'model_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
P = T_overall(T_overall.model_name=="q_spectrum_plus_composition",:);
if ~isempty(P)
    fprintf('Primary model controlled-sweep MAPE %.2f%%, high>20 %.2f%%.\n', ...
        P.MAPE_pct, P.high_error20_pct);
end
X = T_case(T_case.model_name=="q_spectrum_plus_composition",:);
if ~isempty(X)
    X = sortrows(X,'MAPE_pct','descend');
    fprintf('Worst controlled-sweep primary case: %s (%s), MAPE %.2f%%.\n', ...
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
