%% analyze_test_37_ood_q_spectrum_validation.m
% Test 37: OOD validation of the frozen Test 35 q_spectrum_only model.
%
% Goal:
%   Stress-test the simplest clean Test 35 model on simulations outside the
%   original 2/3 m/s training distribution. This script does not retrain.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST37_MODE          = quick | full
%   ADAPTIVE_REQ_TEST37_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST37_SAVE_ALL_MAPS = true | false
%
% Quick mode intentionally uses a small but diagnostic subset:
%   - unseen homogeneous speeds,
%   - unseen bilayer/inclusion contrasts,
%   - new ellipse/off-center/two-inclusion/oblique/thin/three-material cases,
%   - one unseen frequency, M=2, and two new field regimes.
%
% Full mode expands frequencies, regimes, cases, and M values at dx=0.2 mm.

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
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir,'test37_configuration.json'));

fprintf('\nTest 37: OOD validation of frozen Test 35 q_spectrum_only\n');
fprintf('Mode: %s | validate only: %d | save all maps: %d\n', ...
    CFG.Mode, CFG.ValidateOnly, CFG.SaveAllMaps);
fprintf('No training. At least one source is forced/aligned in every regime.\n');

MODEL = load_test35_model(root_dir);
assert_no_forbidden_predictors(MODEL.q_model.features);

if CFG.ValidateOnly
    validate_test37(root_dir, CFG, MODEL);
    fprintf('Test 37 validation-only checks passed.\n');
    return;
end

%% Run OOD design

CASE_SPECS = build_case_specs(CFG);
REGIMES = build_regime_specs(CFG);
fprintf('OOD cases: %d | regimes: %d | frequencies: %d | M: %d\n', ...
    numel(CASE_SPECS), numel(REGIMES), numel(CFG.Frequencies), numel(CFG.M));

parts = {};
condition_rows = {};
condition_id = 0;
for ci = 1:numel(CASE_SPECS)
    C = CASE_SPECS(ci);
    for fi = 1:numel(CFG.Frequencies)
        f0 = CFG.Frequencies(fi);
        for ri = 1:numel(REGIMES)
            R = REGIMES(ri);
            condition_id = condition_id + 1;
            field_key = sprintf('%s__f%g__%s__dx%gum', ...
                C.case_id, f0, lower(R.regime_id), round(1e6*CFG.dx));
            field_file = fullfile(OUT.field_dir, "field__" + sanitize(field_key) + ".mat");
            if exist(field_file,'file') == 2
                S = load(field_file, 'sim', 'cfg_sim', 'case_spec', 'regime_spec');
                sim = S.sim; cfg_sim = S.cfg_sim;
                fprintf('[field] reused %s\n', field_key);
            else
                cfg_sim = build_sim_cfg(CFG, C, R, f0, condition_id);
                sim = run_ood_simulation(cfg_sim, C);
                case_spec = C; regime_spec = R; %#ok<NASGU>
                save(field_file, 'sim', 'cfg_sim', 'case_spec', 'regime_spec', '-v7.3');
                fprintf('[field] built %s | size %dx%d\n', ...
                    field_key, size(sim.Uxz,1), size(sim.Uxz,2));
            end
            for mi = 1:numel(CFG.M)
                M = CFG.M(mi);
                key = string(field_key) + "__M" + string(M);
                result_file = fullfile(OUT.condition_dir, ...
                    "result__" + sanitize(key) + ".mat");
                if exist(result_file,'file') == 2
                    S = load(result_file, 'T_condition', 'T_condition_summary');
                    T_condition = S.T_condition;
                    T_condition_summary = S.T_condition_summary;
                    fprintf('[%s] reused %d patches.\n', key, height(T_condition));
                else
                    timer = tic;
                    T_condition = evaluate_condition(sim, cfg_sim, C, R, M, key, MODEL, CFG);
                    T_condition_summary = summarize_predictions(T_condition, ...
                        ["condition_key","case_id","case_family","field_regime","f0","M"]);
                    save(result_file, 'T_condition', 'T_condition_summary', '-v7.3');
                    fprintf('[%s] evaluated %d patches in %.1f s.\n', ...
                        key, height(T_condition), toc(timer));
                end
                parts{end+1,1} = T_condition; %#ok<AGROW>
                condition_rows{end+1,1} = T_condition_summary; %#ok<AGROW>
                if CFG.SaveAllMaps || should_plot_representative(CFG, C, R, f0, M)
                    plot_condition_maps(T_condition, OUT);
                end
            end
        end
    end
end

T_all = vertcat(parts{:});
T_condition_summary = vertcat(condition_rows{:});

%% Tables

T_overall = summarize_predictions(T_all, "model_name");
T_by_case = summarize_predictions(T_all, ["model_name","case_id"]);
T_by_family = summarize_predictions(T_all, ["model_name","case_family"]);
T_by_regime = summarize_predictions(T_all, ["model_name","field_regime"]);
T_by_frequency = summarize_predictions(T_all, ["model_name","f0"]);
T_by_M = summarize_predictions(T_all, ["model_name","M"]);
T_by_purity = summarize_predictions(T_all, ["model_name","purity_bin"]);
T_by_distance = summarize_predictions(T_all, ["model_name","distance_bin"]);
T_fail = top_failure_conditions(T_condition_summary, 30);

writetable(remove_cell_columns(T_all), fullfile(OUT.table_dir, ...
    'test37_ood_patch_level_predictions.csv'));
writetable(T_condition_summary, fullfile(OUT.table_dir, ...
    'test37_ood_condition_summary.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test37_ood_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test37_ood_summary_by_case.csv'));
writetable(T_by_family, fullfile(OUT.table_dir, 'test37_ood_summary_by_family.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test37_ood_summary_by_regime.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'test37_ood_summary_by_frequency.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test37_ood_summary_by_M.csv'));
writetable(T_by_purity, fullfile(OUT.table_dir, 'test37_ood_summary_by_purity_bin.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test37_ood_summary_by_distance_bin.csv'));
writetable(T_fail, fullfile(OUT.table_dir, 'test37_ood_worst_conditions.csv'));
save(fullfile(OUT.data_dir,'test37_ood_results.mat'), ...
    'T_all','T_condition_summary','T_overall','T_by_case','T_by_family', ...
    'T_by_regime','T_by_frequency','T_by_M','T_by_purity','T_by_distance', ...
    'T_fail','CFG','-v7.3');

%% Figures

plot_summary_bars(T_overall, T_by_family, T_by_purity, OUT);
plot_heatmap_table(T_by_case, "case_id", OUT, 'test37_mape_by_case.png');
plot_heatmap_table(T_by_regime, "field_regime", OUT, 'test37_mape_by_regime.png');
plot_heatmap_table(T_by_frequency, "f0", OUT, 'test37_mape_by_frequency.png');
plot_heatmap_table(T_by_M, "M", OUT, 'test37_mape_by_M.png');
plot_error_vs_predicted_purity(T_all, OUT);
plot_worst_condition_gallery(T_all, T_fail, OUT);

print_interpretation(T_overall, T_by_family, T_by_purity, T_fail, OUT);

fprintf('\nTables: %s\nFigures: %s\nTest 37 complete.\n', ...
    OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST37_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), ...
    'ADAPTIVE_REQ_TEST37_MODE must be quick or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST37_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST37_SAVE_ALL_MAPS', true);
CFG.RandomSeed = 37001;
CFG.dx = 0.2e-3;
CFG.dz = 0.2e-3;
CFG.Lx = 0.05;
CFG.Lz = 0.05;
CFG.TargetStepM = 1.0e-3;
CFG.cs_guess = 3;
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
CFG.UseParfor = false;
if CFG.QuickMode
    CFG.Frequencies = 450;
    CFG.M = 2;
    CFG.RegimeSet = "quick";
else
    CFG.Frequencies = [350 450 550 650];
    CFG.M = [2 3 4];
    CFG.RegimeSet = "full";
end
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_37_ood_q_spectrum_validation');
if CFG.QuickMode
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
end
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.field_dir = fullfile(OUT.data_dir, 'field_cache');
OUT.condition_dir = fullfile(OUT.data_dir, 'condition_results');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.fail_dir = fullfile(OUT.figure_dir, 'worst_conditions');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function MODEL = load_test35_model(root_dir)
file = fullfile(root_dir, 'outputs', 'test_35_spectral_composition_to_q_model', ...
    'models', 'test35_spectral_composition_q_models.mat');
assert(exist(file,'file') == 2, 'Missing Test 35 full model bundle: %s', file);
S = load(file, 'MODEL_BUNDLE');
B = S.MODEL_BUNDLE;
MODEL = struct();
MODEL.bundle_file = file;
MODEL.base_features = B.BASE_FEATURES;
MODEL.q_model = B.MODELS.q.spectrum_only;
MODEL.composition = B.MODELS.composition;
MODEL.model_name = "q_spectrum_only_Test35_frozen";
end

function assert_no_forbidden_predictors(features)
bad_patterns = ["true","oracle","purity","mixed","confidence","error", ...
    "pred","sws","cs_","k_true","q_local","q_pred","q_theory", ...
    "req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition"];
features = lower(string(features));
for p = bad_patterns
    hit = features(contains(features, p));
    assert(isempty(hit), 'Forbidden Test37 predictor detected: %s', strjoin(hit, ', '));
end
fprintf('q_spectrum_only predictors passed leakage guard (%d predictors).\n', numel(features));
end

function validate_test37(root_dir, CFG, MODEL)
assert(exist(MODEL.bundle_file,'file') == 2);
C = build_case_specs(CFG); R = build_regime_specs(CFG);
cfg = build_sim_cfg(CFG, C(1), R(1), CFG.Frequencies(1), 1);
sim = run_ood_simulation(cfg, C(1));
T = evaluate_condition(sim, cfg, C(1), R(1), CFG.M(1), "validate", MODEL, CFG);
assert(~isempty(T) && all(isfinite(T.sws_pred)));
assert(any(strcmp('predicted_patch_purity', T.Properties.VariableNames)));
fprintf('Validation sample: %s, %d patches, MAPE %.2f%%.\n', ...
    C(1).case_id, height(T), mean(T.sws_abs_error_pct,'omitnan'));
end

%% OOD design

function cases = build_case_specs(CFG)
base = [
    homogeneous_case("homogeneous_cs1p5", 1.5)
    homogeneous_case("homogeneous_cs4p0", 4.0)
    bilayer_case("bilayer_1p5_3", 1.5, 3.0, 0)
    bilayer_case("bilayer_2_4", 2.0, 4.0, 0)
    inclusion_case("inclusion_1p5_3", 1.5, 3.0, [0.025 0.025], 0.008)
    inclusion_case("inclusion_2_4", 2.0, 4.0, [0.025 0.025], 0.008)
    ellipse_case("ellipse_2_4", 2.0, 4.0, [0.025 0.025], [0.013 0.006], deg2rad(25))
    inclusion_case("offcenter_inclusion_2_4", 2.0, 4.0, [0.032 0.022], 0.007)
    two_inclusion_case("two_inclusions_2_4", 2.0, 4.0)
    bilayer_case("oblique_bilayer_2_4", 2.0, 4.0, deg2rad(25))
    thin_layer_case("thin_layer_2_4", 2.0, 4.0)
    three_material_case("three_material_2_3_4", 2.0, 3.0, 4.0)
    smooth_bilayer_case("smooth_bilayer_2_4", 2.0, 4.0, 1.2e-3)
    ];
if CFG.QuickMode
    keep = ["homogeneous_cs1p5","homogeneous_cs4p0","bilayer_2_4", ...
        "inclusion_2_4","ellipse_2_4","offcenter_inclusion_2_4", ...
        "oblique_bilayer_2_4","three_material_2_3_4"];
    cases = base(ismember([base.case_id], keep));
else
    more = [
        homogeneous_case("homogeneous_cs2p5", 2.5)
        homogeneous_case("homogeneous_cs3p5", 3.5)
        bilayer_case("bilayer_2p5_3p5", 2.5, 3.5, 0)
        inclusion_case("inclusion_2p5_3p5", 2.5, 3.5, [0.025 0.025], 0.008)
        ];
    cases = [base; more];
end
end

function R = build_regime_specs(CFG)
R = [
    regime("directional_2D_new_angle", "directional_2D", true, 1, true, ...
    "ranges", "fibonacci", [pi/6 pi/6], [0 2*pi], [0 pi], 11)
    regime("diffuse_3D_new_seed", "diffuse_3D", false, 128, true, ...
    "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 29)
    ];
if ~CFG.QuickMode
    R = [R
        regime("diffuse_2D_new_seed", "diffuse_2D", true, 128, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 31)
        regime("partial_3D_16src", "partial_3D", false, 16, true, ...
        "ranges", "fibonacci", [0 2*pi], [0 2*pi], [0 pi], 37)];
end
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
F.patch_purity = nan(n,1);
F.q_oracle = nan(n,1);
half = floor(win_size/2);
for i = 1:n
    cx = F.cx(i); cz = F.cz(i);
    F.true_SWS(i) = sim.cs_map(cz, cx);
    patch = sim.cs_map((cz-half):(cz+half), (cx-half):(cx+half));
    F.patch_purity(i) = dominant_fraction(patch);
    F.q_oracle(i) = invert_mapping_to_q(F.req_mapping{i}, 2*pi*cfg.f0/F.true_SWS(i));
end
F.purity_bin = purity_bin(F.patch_purity);
D = distance_to_material_boundary(sim.cs_map, cfg.dx);
F.distance_to_boundary_mm = D(sub2ind(size(D), F.cz, F.cx));
F.distance_bin = distance_bin(F.distance_to_boundary_mm, CFG.DistanceEdgesMm);
F.q_theory_prior = repmat(theory_q_for_regime(cfg, M, R.field_regime), n, 1);
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
                error('Required Test35 predictor missing in Test37 table: %s', f);
        end
    end
end
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
out_file = fullfile(out_dir, "test37_map__" + sanitize(key) + ".png");
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
tf = CFG.QuickMode || (M == 2 && ismember(f0, [350 650]) && ...
    ismember(R.regime_id, ["directional_2D_new_angle","diffuse_3D_new_seed"]) && ...
    ismember(C.case_family, ["homogeneous","bilayer","inclusion","ellipse", ...
    "three_material","smooth_transition"]));
end

function plot_summary_bars(T_overall, T_family, T_purity, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 8]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact');
ax = nexttile(tl); bar(ax, categorical(T_overall.model_name), T_overall.MAPE_pct);
ylabel(ax,'MAPE (%)'); title(ax,'Overall','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl); bar(ax, categorical(T_family.case_family), T_family.MAPE_pct);
xtickangle(ax,30); title(ax,'By family','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl); bar(ax, categorical(T_purity.purity_bin), T_purity.MAPE_pct);
xtickangle(ax,30); title(ax,'By true patch purity','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test37_summary_bars.png'));
end

function plot_heatmap_table(T, col, OUT, name)
groups = unique(string(T.(col)), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 22 5]);
bar(categorical(groups), T.MAPE_pct);
xtickangle(30); ylabel('MAPE (%)'); title("MAPE by " + string(col), 'FontWeight','normal');
grid on; export_fig(fig, fullfile(OUT.figure_dir,name));
end

function plot_error_vs_predicted_purity(T, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 9]);
scatter(T.predicted_patch_purity, T.sws_abs_error_pct, 6, T.p_mixed, ...
    'filled', 'MarkerFaceAlpha', 0.15);
xlabel('Predicted patch purity'); ylabel('Abs SWS error (%)');
title('OOD error vs predicted composition','FontWeight','normal');
grid on; colorbar;
export_fig(fig, fullfile(OUT.figure_dir,'test37_error_vs_predicted_purity.png'));
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
fprintf('\n================ Test 37 OOD summary ================\n');
disp(T_overall(:, {'model_name','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct','mean_p_mixed'}));
fprintf('\nWorst OOD conditions:\n');
disp(T_fail(:, {'condition_key','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct','mean_p_mixed'}));
[~, wi] = max(T_family.MAPE_pct);
fprintf('\nWorst family: %s, MAPE %.2f%%.\n', ...
    T_family.case_family(wi), T_family.MAPE_pct(wi));
[~, pi] = max(T_purity.MAPE_pct);
fprintf('Worst purity bin: %s, MAPE %.2f%%.\n', ...
    T_purity.purity_bin(pi), T_purity.MAPE_pct(pi));
fprintf('=====================================================\n');
fprintf('Worst-condition maps: %s\n', OUT.fail_dir);
end

%% Numerical helpers

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

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3, nz = max(T.map_iz); nx = max(T.map_ix); end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
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

function tf = env_true(name, default)
raw = strtrim(lower(string(getenv(name))));
if raw == "", tf = logical(default); else, tf = ismember(raw, ["1","true","yes","on"]); end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
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
