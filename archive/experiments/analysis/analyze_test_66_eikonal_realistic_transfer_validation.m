%% analyze_test_66_eikonal_realistic_transfer_validation.m
% Test 66: zero-shot Eikonal realistic transfer validation.
%
% This test evaluates frozen Test53/Test38 q models on realistic 2.5D
% Eikonal/readout fields from wave_sim_project_clean. Nothing is retrained:
% q_spectrum_only, q_spectrum_plus_composition, and the frozen composition
% auxiliaries are loaded from the saved model bundle and applied as-is.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST66_MODE          = validate | quick | full_a | full_b
%   ADAPTIVE_REQ_TEST66_WAVE_SIM_PATH = /Users/sara/Documents/wave_sim_project_clean
%   ADAPTIVE_REQ_TEST66_MODEL_SOURCE  = test53 | full | /path/to/bundle.mat
%   ADAPTIVE_REQ_TEST66_M_VALUES      = comma-separated M list, e.g. 2 or 2,3
%   ADAPTIVE_REQ_TEST66_TARGET_STEP_M = 0.0005
%   ADAPTIVE_REQ_TEST66_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST66_USE_PARFOR    = true | false
%   ADAPTIVE_REQ_TEST66_PLOT_ONLY     = true | false
%   ADAPTIVE_REQ_TEST66_EXPERIMENT    = default | large_inclusion_readout
%   ADAPTIVE_REQ_TEST66_INCLUSION_RADIUS_M = override simple inclusion radius
%   ADAPTIVE_REQ_TEST66_RUN_LABEL     = optional output suffix
%   ADAPTIVE_REQ_TEST66_GEOMETRIES    = optional comma-separated filter
%   ADAPTIVE_REQ_TEST66_FREQUENCIES   = optional comma-separated filter
%   ADAPTIVE_REQ_TEST66_REALISM_LEVELS = optional comma-separated filter
%   ADAPTIVE_REQ_TEST66_FIELD_REGIMES = optional comma-separated filter
%   ADAPTIVE_REQ_TEST66_SOURCE_SEEDS  = optional comma-separated filter
%   ADAPTIVE_REQ_TEST66_NOISE_SEEDS   = optional comma-separated filter
%
% Important: truth maps, patch purity, distance to interface, ROI labels, and
% error variables are diagnostic only. They are never used as model inputs.

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
write_config_json(CFG, fullfile(OUT.root_dir, 'test66_configuration.json'));

fprintf('\nTest 66: zero-shot Eikonal realistic transfer validation\n');
fprintf('Mode: %s | wave_sim_project_clean: %s\n', CFG.Mode, CFG.WaveSimPath);
fprintf('Model source: %s\n', CFG.ModelBundleFile);
fprintf('REQ profile: M=%s, cs_guess=%.2f, target step %.3f mm, Nbins=%s, min=%g, smooth sigma=%g\n', ...
    mat2str(CFG.M), CFG.cs_guess, 1e3*CFG.TargetStepM, string(CFG.REQ.Nbins), ...
    CFG.REQ.Nbins_min, CFG.REQ.SmoothSigma);
fprintf('No q, composition, risk, or correction model is trained in this script.\n');

if CFG.PlotOnly
    fprintf('Plot-only mode: reading existing CSV tables and regenerating figures.\n');
    [T_existing, SUM_existing] = read_existing_results(OUT);
    plot_summary_figures(T_existing, SUM_existing, OUT, CFG);
    fprintf('Plot-only complete. Figures: %s\n', OUT.figure_dir);
    return;
end

assert(exist(fullfile(CFG.WaveSimPath, 'src'), 'dir') == 7, ...
    'Missing wave_sim_project_clean src folder: %s', CFG.WaveSimPath);
addpath(fullfile(CFG.WaveSimPath, 'src'));

BUNDLE = load_model_bundle(CFG);
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);

if CFG.ValidateOnly
    CONDITIONS = build_conditions(CFG);
    CONDITIONS = CONDITIONS(1:min(numel(CONDITIONS),2));
    fprintf('Validation-only conditions: %d\n', numel(CONDITIONS));
    T_validate = run_conditions(CONDITIONS, BUNDLE, BASE_FEATURES, CFG, OUT);
    assert(~isempty(T_validate) && all(isfinite(T_validate.sws_pred)), ...
        'Validation-only produced empty or non-finite predictions.');
    writetable(remove_cell_columns(T_validate), fullfile(OUT.table_dir, 'test66_validation_patch_results.csv'));
    fprintf('Validation-only complete. Rows: %d\n', height(T_validate));
    return;
end

CONDITIONS = build_conditions(CFG);
fprintf('Eikonal transfer conditions: %d | M values: %s\n', numel(CONDITIONS), mat2str(CFG.M));
assert(~isempty(CONDITIONS), ['No Test66 conditions survived the selected filters. ' ...
    'For clean/shear-only realism levels, noise_seed is 0; readout levels use the requested noise seeds.']);
RTE = estimate_runtime(CONDITIONS, BUNDLE, BASE_FEATURES, CFG, OUT);
writetable(RTE, fullfile(OUT.table_dir, 'test66_runtime_estimate.csv'));

T_all = run_conditions(CONDITIONS, BUNDLE, BASE_FEATURES, CFG, OUT);
T_csv = remove_cell_columns(T_all);
writetable(T_csv, fullfile(OUT.table_dir, 'test66_patch_level_results.csv'));

SUM = write_summaries(T_all, OUT);
plot_summary_figures(T_all, SUM, OUT, CFG);
save(fullfile(OUT.data_dir, 'test66_compact_results.mat'), ...
    'CFG','SUM','RTE','-v7.3');

print_summary(SUM.overall, SUM.by_realism_level, SUM.by_roi, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 66 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST66_MODE'))));
if mode == "", mode = "quick"; end
aliases = containers.Map(["validation","validate_only","full"], ["validate","validate","full_a"]);
if isKey(aliases, mode), mode = aliases(mode); end
assert(ismember(mode, ["validate","quick","full_a","full_b"]), ...
    'ADAPTIVE_REQ_TEST66_MODE must be validate, quick, full_a, or full_b.');

CFG = struct();
CFG.Mode = mode;
CFG.ValidateOnly = mode == "validate" || env_true('ADAPTIVE_REQ_TEST66_VALIDATE_ONLY', false);
CFG.QuickMode = mode == "quick";
CFG.FullA = mode == "full_a";
CFG.FullB = mode == "full_b";
CFG.WaveSimPath = char(env_string('ADAPTIVE_REQ_TEST66_WAVE_SIM_PATH', ...
    "/Users/sara/Documents/wave_sim_project_clean"));
CFG.ModelSource = env_string('ADAPTIVE_REQ_TEST66_MODEL_SOURCE', "test53");
CFG.ModelBundleFile = resolve_model_bundle(root_dir, CFG.ModelSource);
CFG.Experiment = lower(env_string('ADAPTIVE_REQ_TEST66_EXPERIMENT', "default"));
CFG.RunLabel = env_string('ADAPTIVE_REQ_TEST66_RUN_LABEL', "");
CFG.FallbackReqCacheDir = env_string('ADAPTIVE_REQ_TEST66_FALLBACK_REQ_CACHE_DIR', "");
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST66_SAVE_ALL_MAPS', true);
CFG.UseParfor = env_true('ADAPTIVE_REQ_TEST66_USE_PARFOR', false);
CFG.PlotOnly = env_true('ADAPTIVE_REQ_TEST66_PLOT_ONLY', false);
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST66_TARGET_STEP_M', 0.5e-3);
CFG.dx = env_number('ADAPTIVE_REQ_TEST66_DX', 0.2e-3);
CFG.dz = env_number('ADAPTIVE_REQ_TEST66_DZ', 0.2e-3);
CFG.Lx = 0.050;
CFG.Lz = 0.050;
CFG.cs_guess = 3.0;
CFG.RandomSeed = 66001;
CFG.MinSourceRoiDistanceM = 4e-3;
CFG.InclusionRadiusM = env_number('ADAPTIVE_REQ_TEST66_INCLUSION_RADIUS_M', 0.007);
CFG.HardInclusionRadiusM = env_number('ADAPTIVE_REQ_TEST66_HARD_INCLUSION_RADIUS_M', 0.0065);
CFG.GeometryTag = "";
if CFG.Experiment == "large_inclusion_readout"
    CFG.InclusionRadiusM = env_number('ADAPTIVE_REQ_TEST66_INCLUSION_RADIUS_M', 0.012);
    CFG.HardInclusionRadiusM = env_number('ADAPTIVE_REQ_TEST66_HARD_INCLUSION_RADIUS_M', 0.012);
    CFG.GeometryTag = sprintf("__r%gmm", round(1e3*CFG.InclusionRadiusM));
end
CFG.ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST66_MODELS', ...
    ["q_spectrum_only","q_spectrum_plus_composition"]);
CFG.FilterGeometries = env_string_list('ADAPTIVE_REQ_TEST66_GEOMETRIES', strings(0,1));
CFG.FilterFrequencies = env_number_list('ADAPTIVE_REQ_TEST66_FREQUENCIES', []);
CFG.FilterRealismLevels = env_string_list('ADAPTIVE_REQ_TEST66_REALISM_LEVELS', strings(0,1));
CFG.FilterFieldRegimes = env_string_list('ADAPTIVE_REQ_TEST66_FIELD_REGIMES', strings(0,1));
CFG.FilterSourceSeeds = env_number_list('ADAPTIVE_REQ_TEST66_SOURCE_SEEDS', []);
CFG.FilterNoiseSeeds = env_number_list('ADAPTIVE_REQ_TEST66_NOISE_SEEDS', []);
CFG.REQ = struct('Gamma',1,'PadFactor',1,'EdgeMode',"valid", ...
    'Nbins',"auto",'Nbins_auto_oversample',1,'Nbins_min',16,'SmoothSigma',1);
CFG.DistanceEdgesMm = [0 1 2 4 8 Inf];
CFG.DistanceOverWindowEdges = [0 0.25 0.5 1 2 Inf];
CFG.DepthEdgesMm = [0 10 20 30 40 Inf];
CFG.SnrEdgesDb = [-Inf 5 10 15 20 30 40 Inf];
CFG.AmplitudeQuantileBins = 5;

if CFG.ValidateOnly
    CFG.M = 2;
elseif CFG.QuickMode
    CFG.M = 2;
else
    CFG.M = [2 3];
end
CFG.M = env_number_list('ADAPTIVE_REQ_TEST66_M_VALUES', CFG.M);
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.FullA || CFG.FullB
    base_root_dir = fullfile(root_dir, 'outputs', 'test_66_eikonal_realistic_transfer_validation', char(CFG.Mode));
else
    base_root_dir = fullfile(root_dir, 'outputs', 'test_66_eikonal_realistic_transfer_validation', char(CFG.Mode));
end
OUT.root_dir = base_root_dir;
if CFG.Experiment ~= "default"
    OUT.root_dir = fullfile(OUT.root_dir, char(CFG.Experiment));
    base_root_dir = OUT.root_dir;
end
if strlength(string(CFG.RunLabel)) > 0
    OUT.root_dir = fullfile(OUT.root_dir, char(sanitize("run_" + string(CFG.RunLabel))));
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.field_dir = fullfile(OUT.data_dir, 'field_cache');
OUT.req_dir = fullfile(OUT.data_dir, 'req_cache');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end

% A labeled run is often a subset/re-analysis of an already completed Test66
% run. Reusing the parent REQ cache prevents paying the expensive
% req_estimator_map cost again for conditions that already exist.
OUT.req_fallback_dirs = strings(0,1);
if strlength(string(CFG.FallbackReqCacheDir)) > 0
    OUT.req_fallback_dirs = string(CFG.FallbackReqCacheDir);
elseif strlength(string(CFG.RunLabel)) > 0
    parent_req_dir = string(fullfile(base_root_dir, 'data', 'req_cache'));
    if parent_req_dir ~= string(OUT.req_dir)
        OUT.req_fallback_dirs = parent_req_dir;
    end
end
end

function file = resolve_model_bundle(root_dir, source)
src = string(source);
if exist(src, 'file') == 2
    file = char(src); return;
end
switch lower(src)
    case {"test53","paper","final"}
        file = fullfile(root_dir, 'outputs', 'test_53_paper_final_clean_q_training', ...
            'models', 'test38_velocity_field_diverse_q_models.mat');
    case "test53_quick"
        file = fullfile(root_dir, 'outputs', 'test_53_paper_final_clean_q_training', ...
            'quick', 'models', 'test38_velocity_field_diverse_q_models.mat');
    case "medium"
        file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
            'medium', 'models', 'test38_velocity_field_diverse_q_models.mat');
    case "quick"
        file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
            'quick', 'models', 'test38_velocity_field_diverse_q_models.mat');
    otherwise
        file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
            'models', 'test38_velocity_field_diverse_q_models.mat');
end
end

function B = load_model_bundle(CFG)
assert(exist(CFG.ModelBundleFile,'file') == 2, ...
    'Missing model bundle: %s', CFG.ModelBundleFile);
S = load(CFG.ModelBundleFile, 'MODEL_BUNDLE');
B = S.MODEL_BUNDLE;
B.bundle_file = CFG.ModelBundleFile;
assert(isfield(B.MODELS,'composition'), 'Model bundle lacks frozen composition models.');
assert(isfield(B.MODELS,'q') && isfield(B.MODELS.q,'spectrum_only') && ...
    isfield(B.MODELS.q,'spectrum_plus_composition'), ...
    'Model bundle lacks q_spectrum_only or q_spectrum_plus_composition.');
fprintf('Loaded frozen bundle with %d model entries.\n', numel(B.MODELS.model_names));
end

%% Conditions

function C = build_conditions(CFG)
if CFG.Experiment == "large_inclusion_readout"
    geos = ["inclusion_2_3","inclusion_2_4"];
    freqs = env_number_list('ADAPTIVE_REQ_TEST66_EXPERIMENT_FREQS', [200 300 500]);
    levels = "readout_medium";
    regimes = "diffuse_like_8src";
    source_seed_candidates = env_number_list('ADAPTIVE_REQ_TEST66_EXPERIMENT_SOURCE_SEEDS', 1);
    noise_seeds = env_number_list('ADAPTIVE_REQ_TEST66_EXPERIMENT_NOISE_SEEDS', 1);
elseif CFG.ValidateOnly
    geos = ["homogeneous_cs2","bilayer_inclusion_2_3_4"];
    freqs = 300;
    levels = "clean";
    regimes = "single_source_lateral";
    source_seed_candidates = 1;
    noise_seeds = 1;
elseif CFG.QuickMode
    geos = ["homogeneous_cs2","homogeneous_cs4","inclusion_2_4","bilayer_inclusion_2_3_4"];
    freqs = [300 500];
    levels = ["clean","readout_medium"];
    regimes = ["single_source_lateral","diffuse_like_8src"];
    source_seed_candidates = 1:2;
    noise_seeds = 1:2;
else
    geos = ["homogeneous_cs2","homogeneous_cs3","homogeneous_cs4", ...
        "bilayer_2_3","bilayer_2_4","inclusion_2_3","inclusion_2_4","bilayer_inclusion_2_3_4"];
    freqs = [200 300 400 500 600];
    regimes = ["single_source_lateral","single_source_diagonal","diffuse_like_8src"];
    source_seed_candidates = 1:2;
    noise_seeds = 1:3;
    if CFG.FullB
        levels = ["clean","shear_attenuation_only","readout_soft","readout_medium","readout_hard"];
        regimes = [regimes "diffuse_like_16src"];
    else
        levels = ["clean","shear_attenuation_only","readout_medium"];
    end
end

C = repmat(empty_condition(), 0, 1);
for gi = 1:numel(geos)
    G = geometry_spec(geos(gi));
    for fi = 1:numel(freqs)
        for li = 1:numel(levels)
            for ri = 1:numel(regimes)
                source_seeds = source_seed_list(regimes(ri), source_seed_candidates);
                for si = 1:numel(source_seeds)
                    ns = noise_seed_list(levels(li), noise_seeds);
                    for ni = 1:numel(ns)
                        X = empty_condition();
                        X.geometry = G.geometry;
                        X.geometry_family = G.geometry_family;
                        X.cs_bg = G.cs_bg;
                        X.cs_inc = G.cs_inc;
                        X.cs_mid = G.cs_mid;
                        X.cs_hard = G.cs_hard;
                        X.f0 = freqs(fi);
                        X.realism_level = levels(li);
                        X.field_regime = regimes(ri);
                        X.source_seed = source_seeds(si);
                        X.noise_seed = ns(ni);
                        source_layout_tag = "";
                        if startsWith(string(X.field_regime), "diffuse_like")
                            source_layout_tag = "__layoutv2";
                        end
                        X.condition_key = sprintf('%s%s__f%g__%s__%s__src%d__noise%d%s', ...
                            X.geometry, CFG.GeometryTag, X.f0, X.realism_level, X.field_regime, ...
                            X.source_seed, X.noise_seed, source_layout_tag);
                        C(end+1,1) = X; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
C = filter_conditions(C, CFG);
end

function C = filter_conditions(C, CFG)
if isempty(C), return; end
if ~isempty(CFG.FilterGeometries)
    C = C(ismember([C.geometry], CFG.FilterGeometries));
end
if ~isempty(CFG.FilterFrequencies)
    C = C(ismember([C.f0], CFG.FilterFrequencies));
end
if ~isempty(CFG.FilterRealismLevels)
    C = C(ismember([C.realism_level], CFG.FilterRealismLevels));
end
if ~isempty(CFG.FilterFieldRegimes)
    C = C(ismember([C.field_regime], CFG.FilterFieldRegimes));
end
if ~isempty(CFG.FilterSourceSeeds)
    C = C(ismember([C.source_seed], CFG.FilterSourceSeeds));
end
if ~isempty(CFG.FilterNoiseSeeds)
    C = C(ismember([C.noise_seed], CFG.FilterNoiseSeeds));
end
end

function ss = source_seed_list(regime, candidates)
% Deterministic single-source regimes do not use source_seed, so running
% seed 1:2 only creates duplicate clean/attenuation fields. Multi-source
% diffuse-like regimes use source_seed to sample independent source layouts.
if startsWith(string(regime), "diffuse_like")
    ss = candidates;
else
    ss = candidates(1);
end
end

function ns = noise_seed_list(level, candidates)
if startsWith(string(level), "readout")
    ns = candidates;
else
    ns = 0;
end
end

function G = geometry_spec(name)
G = struct('geometry', string(name), 'geometry_family', "", ...
    'cs_bg', NaN, 'cs_inc', NaN, 'cs_mid', NaN, 'cs_hard', NaN);
switch string(name)
    case "homogeneous_cs2"
        G.geometry_family = "homogeneous"; G.cs_bg = 2; G.cs_inc = 2;
    case "homogeneous_cs3"
        G.geometry_family = "homogeneous"; G.cs_bg = 3; G.cs_inc = 3;
    case "homogeneous_cs4"
        G.geometry_family = "homogeneous"; G.cs_bg = 4; G.cs_inc = 4;
    case "bilayer_2_3"
        G.geometry_family = "bilayer"; G.cs_bg = 2; G.cs_inc = 3;
    case "bilayer_2_4"
        G.geometry_family = "bilayer"; G.cs_bg = 2; G.cs_inc = 4;
    case "inclusion_2_3"
        G.geometry_family = "inclusion"; G.cs_bg = 2; G.cs_inc = 3;
    case "inclusion_2_4"
        G.geometry_family = "inclusion"; G.cs_bg = 2; G.cs_inc = 4;
    case "bilayer_inclusion_2_3_4"
        G.geometry_family = "bilayer_inclusion_three_material";
        G.cs_bg = 2; G.cs_mid = 3; G.cs_hard = 4; G.cs_inc = 4;
    otherwise
        error('Unknown Test66 geometry: %s', name);
end
end

function X = empty_condition()
X = struct('geometry',"",'geometry_family',"",'cs_bg',NaN,'cs_inc',NaN, ...
    'cs_mid',NaN,'cs_hard',NaN,'f0',NaN,'realism_level',"", ...
    'field_regime',"",'source_seed',NaN,'noise_seed',NaN,'condition_key',"");
end

%% Main loop

function T_all = run_conditions(CONDITIONS, BUNDLE, BASE_FEATURES, CFG, OUT)
parts = {};
for ci = 1:numel(CONDITIONS)
    C = CONDITIONS(ci);
    fprintf('\n[%d/%d] %s\n', ci, numel(CONDITIONS), C.condition_key);
    t0 = tic;
    [SIM, sim_cfg] = get_or_build_field(C, CFG, OUT);
    fprintf('  field size %s | min/max SNR %.1f/%.1f dB | %.1f s\n', ...
        mat2str(size(SIM.U)), min(SIM.tracking_snr_db(:),[],'omitnan'), ...
        max(SIM.tracking_snr_db(:),[],'omitnan'), toc(t0));
    for mi = 1:numel(CFG.M)
        M = CFG.M(mi);
        key = sprintf('%s__M%g__step%gum', C.condition_key, M, round(1e6*CFG.TargetStepM));
        t1 = tic;
        F = get_or_extract_req(SIM, sim_cfg, C, M, key, CFG, OUT);
        F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
        validate_features_for_prediction(F, BASE_FEATURES);
        T = apply_models(F, BUNDLE.MODELS, BASE_FEATURES, CFG);
        parts{end+1,1} = T; %#ok<AGROW>
        fprintf('  M=%g: %d windows x %d models in %.1f s.\n', ...
            M, height(F), numel(unique(T.model_name)), toc(t1));
        if CFG.SaveAllMaps
            plot_condition_maps(T, F, SIM, key, OUT, CFG);
        end
    end
end
T_all = vertcat(parts{:});
end

function [SIM, cfg] = get_or_build_field(C, CFG, OUT)
file = fullfile(OUT.field_dir, "field__" + sanitize(C.condition_key) + ".mat");
if exist(file,'file') == 2
    S = load(file, 'SIM', 'cfg');
    SIM = S.SIM; cfg = S.cfg; return;
end
cfg = build_simcore_cfg(C, CFG);
if startsWith(C.field_regime, "diffuse_like")
    out = simulate_diffuse_like(cfg, C, CFG);
else
    out = simcore.pipelines.runREQReadySimulation(cfg);
    out.source_positions_xyz = cfg.Source.PositionXYZ;
    out.source_amplitudes = cfg.Source.Amplitude;
end
SIM = compact_sim_output(out, C, CFG);
save(file, 'SIM', 'cfg', 'C', '-v7.3');
end

function cfg = build_simcore_cfg(C, CFG)
cfg = simcore.config.defaultConfig();
cfg.Name = C.condition_key;
cfg.f0 = C.f0;
cfg.Lx = CFG.Lx; cfg.Lz = CFG.Lz;
cfg.dx = CFG.dx; cfg.dz = CFG.dz;
cfg.fUS_MHz = 5.0;
cfg.Measurement.Axis = [0 0 1];
cfg.Source.Type = "finite_lateral";
cfg.Source.Axis = [0 0 1];
cfg.Source.MotionAxis = [0 0 1];
cfg.Source.ApertureMM = 4;
cfg.Source.NumPoints = 11;
cfg.Source.Amplitude = 1;
cfg.Shear.PhaseModel = "eikonal";
cfg.Shear.GeometricDecayPower = 0.5;
cfg.Shear.Attenuation.PathModel = "same_as_phase";
cfg.Shear.Attenuation.IntegrationMethod = "trapezoid";
cfg.Shear.InterfaceRT.Enabled = false;
cfg.Shear.InterfaceShadow.Enabled = false;
cfg.Ultrasound.PathModel = "compound_steering";
cfg.Ultrasound.SteeringAnglesDeg = [-10 0 10];
cfg.Ultrasound.NoiseCouplingMode = "relative_snr_from_echo";
cfg.Ultrasound.ReferenceSNRdB = 45;
cfg.Noise.AmplitudeNormalization = "fixed_reference";
cfg.Noise.ReferenceAmplitude = 1;
cfg.Noise.Seed = 100000 + 100*C.source_seed + C.noise_seed + round(C.f0);
cfg.Shear.Noise.Seed = 200000 + 100*C.source_seed + C.noise_seed + round(C.f0);
cfg = apply_geometry(cfg, C, CFG);
cfg = apply_source_regime(cfg, C, CFG);
cfg = apply_realism(cfg, C.realism_level);
end

function cfg = apply_geometry(cfg, C, CFG)
soft = material_powerlaw("soft", C.cs_bg, 4, 0.50, 0);
inc = material_powerlaw("inclusion", C.cs_inc, 18, 1.30, 5);
switch C.geometry_family
    case "homogeneous"
        cfg.Geometry.Type = "homogeneous";
        cfg.Materials.Background = soft;
        cfg.Materials.Inclusion = soft;
    case "bilayer"
        cfg.Geometry.Type = "bilayer";
        cfg.Geometry.Params.InterfaceX = 0.025;
        cfg.Materials.Background = soft;
        cfg.Materials.Inclusion = material_powerlaw("hard_layer", C.cs_inc, 16, 1.20, 4);
    case "inclusion"
        cfg.Geometry.Type = "inclusion";
        cfg.Geometry.Params.Center = [0.030 0.025];
        cfg.Geometry.Params.Radius = CFG.InclusionRadiusM;
        cfg.Materials.Background = soft;
        cfg.Materials.Inclusion = inc;
    case "bilayer_inclusion_three_material"
        cfg.Geometry.Type = "bilayer_inclusion_three_material";
        cfg.Geometry.Params.InterfaceX = 0.025;
        cfg.Geometry.Params.InclusionCenter = [0.036 0.025];
        cfg.Geometry.Params.InclusionRadius = CFG.HardInclusionRadiusM;
        cfg.Materials.Background = material_powerlaw("soft_background", C.cs_bg, 4, 0.50, 0);
        cfg.Materials.Inclusion = material_powerlaw("intermediate_layer", C.cs_mid, 12, 0.90, 2);
        cfg.Materials.HardInclusion = material_powerlaw("hard_inclusion", C.cs_hard, 24, 1.50, 6);
    otherwise
        error('Unsupported geometry family: %s', C.geometry_family);
end
end

function mat = material_powerlaw(name, cs, alpha0, usAlpha0, backscatterDB)
mat = simcore.materials.makeMaterial(name, 'Model','empirical_power_law', ...
    'c', cs, 'cRef', cs, 'fRef', 400, 'beta', 0.0, ...
    'alpha0', alpha0, 'alphaFRef', 400, 'alphaExponent', 1.1, ...
    'USAlpha0', usAlpha0, 'USAlphaPower', 1.0, ...
    'BackscatterDB', backscatterDB);
end

function cfg = apply_source_regime(cfg, C, CFG)
rng(C.source_seed + 66100);
switch C.field_regime
    case "single_source_lateral"
        cfg.Source.PositionXYZ = [-0.010 0 0.025];
    case "single_source_diagonal"
        cfg.Source.PositionXYZ = [-0.010 0 -0.010];
    otherwise
        cfg.Source.PositionXYZ = [-0.010 0 0.025];
end
cfg.Source.Amplitude = 1;
cfg.Source.NumPoints = 11;
cfg.Source.ApertureMM = 4;
cfg.Source.MotionAxis = [0 0 1];
cfg.Measurement.Axis = [0 0 1];
cfg.Shear.Noise.ReferenceDistanceM = hypot(CFG.Lx, CFG.Lz);
end

function cfg = apply_realism(cfg, level)
cfg.Shear.Attenuation.Enabled = false;
cfg.Shear.Noise.Enabled = false;
cfg.Shear.Noise.PropagationNoiseEnabled = false;
cfg.Noise.Enabled = false;
cfg.Noise.IncludeShearAmplitudeSNR = false;
switch string(level)
    case "clean"
        return;
    case "shear_attenuation_only"
        cfg.Shear.Attenuation.Enabled = true;
    case {"readout_soft","readout_medium","readout_hard"}
        preset = extractAfter(string(level), "readout_");
        cfg.Shear.Attenuation.Enabled = true;
        cfg = simcore.readout.noisePreset(preset, cfg);
        cfg.Noise.Enabled = true;
        cfg.Noise.IncludeShearAmplitudeSNR = true;
        cfg.Noise.ReferenceAmplitude = 1;
    otherwise
        error('Unsupported realism level: %s', level);
end
end

function out = simulate_diffuse_like(cfg, C, CFG)
N = sscanf(C.field_regime, 'diffuse_like_%dsrc');
if isempty(N), N = 8; end
S = source_set(N, CFG.Lx, CFG.Lz, C.source_seed);
components = cell(N,1);
for k = 1:N
    ck = cfg;
    ck.Source.PositionXYZ = S.positions(k,:);
    ck.Source.Amplitude = abs(S.amps(k));
    ck.Noise.Enabled = false;
    ck.Shear.Noise.Enabled = false;
    components{k} = simcore.shear.simulate2p5D(ck);
end
shearOut = components{1};
U = zeros(size(shearOut.U_clean));
ampSum = zeros(size(U));
attInt = zeros(size(U));
attFac = zeros(size(U));
measW = zeros(size(U));
wSum = 0;
for k = 1:N
    Ck = components{k};
    ph = S.amps(k) ./ max(abs(S.amps(k)), eps);
    w = abs(S.amps(k));
    U = U + ph .* Ck.U_clean;
    ampSum = ampSum + abs(ph).*Ck.component_amplitude_sum;
    attInt = attInt + w.*Ck.shear_attenuation_integral;
    attFac = attFac + w.*Ck.shear_attenuation_factor;
    measW = measW + w.*Ck.measured_component_weight;
    wSum = wSum + w;
end
wSum = max(wSum, eps);
shearOut.U_clean = U;
shearOut.U = U;
shearOut.amplitude = abs(U);
shearOut.phase = angle(U);
shearOut.component_amplitude_sum = ampSum;
shearOut.shear_attenuation_integral = attInt ./ wSum;
shearOut.shear_attenuation_factor = attFac ./ wSum;
shearOut.measured_component_weight = measW ./ wSum;
shearOut.diag.MultiSourceName = char(C.field_regime);
shearOut.diag.NumIndependentSources = N;
shearOut.multi_source_positions_xyz = S.positions;
shearOut.multi_source_amplitudes = S.amps;
out = simcore.readout.applyReadout(shearOut, cfg);
out.cs_map = out.medium.cs_map;
out.true_SWS = out.medium.cs_map;
out.alpha_shear_map = out.medium.alpha_shear_map;
out.alpha_us_map = out.medium.alpha_us_map;
out.label_map = out.medium.label_map;
out.material_fraction = out.medium.material_fraction;
out.signed_distance_m = out.medium.signed_distance_m;
out.abs_distance_m = out.medium.abs_distance_m;
out.x = out.medium.x;
out.z = out.medium.z;
out.dx = median(diff(out.x));
out.dz = median(diff(out.z));
out.f0 = cfg.f0;
out.fUS_MHz = cfg.fUS_MHz;
out.amplitude = abs(out.U);
out.phase = angle(out.U);
out.U_measured = out.U;
out.U_shear_clean = out.U_clean;
out.U_shear_physical = out.U_shear_physical;
out.amplitude_clean = abs(out.U_clean);
out.phase_clean = angle(out.U_clean);
out.tracking_snr_db = out.echo.tracking_snr_db;
out.readout_noise_weight = out.noise_weight;
out.roi_masks = simcore.geometry.roiMasks(out.medium);
out.source_positions_xyz = S.positions;
out.source_amplitudes = S.amps;
out.config = cfg;
out.adapter_name = 'simcore_clean_test66_diffuse_like_v1';
end

function S = source_set(N, Lx, Lz, seed)
rng(7000 + seed + N);
% Keep two sources exactly in the measurement plane so every diffuse-like
% condition has visible in-plane contributions. The remaining sources are
% sampled around the computational box with nonzero y offsets, which gives
% a more genuinely 3D source distribution than the older fixed template
% that placed four sources in-plane for N=8.
P = zeros(N, 3);
n_in_plane = min(N, 2);
if n_in_plane >= 1
    P(1,:) = [-0.010, 0, Lz/2];
end
if n_in_plane >= 2
    P(2,:) = [Lx + 0.010, 0, Lz/2];
end
y_min = 8e-3;
y_max = 20e-3;
outside = 10e-3;
for k = (n_in_plane+1):N
    side = randi(4);
    y = (y_min + (y_max-y_min)*rand) * sign(rand-0.5);
    switch side
        case 1 % left
            P(k,:) = [-outside, y, Lz*rand];
        case 2 % right
            P(k,:) = [Lx+outside, y, Lz*rand];
        case 3 % bottom
            P(k,:) = [Lx*rand, y, -outside];
        case 4 % top
            P(k,:) = [Lx*rand, y, Lz+outside];
    end
end
phi = 2*pi*rand(N,1);
amp = (0.8 + 0.4*rand(N,1)) ./ sqrt(N);
amp(1:min(2,N)) = 1/sqrt(N);
S = struct('positions', P, 'amps', amp .* exp(1i*phi));
end

function SIM = compact_sim_output(out, C, CFG)
SIM = struct();
SIM.U = out.U;
SIM.Uxz = out.U;
SIM.U_clean = out.U_shear_clean;
SIM.U_shear_physical = out.U_shear_physical;
SIM.x = out.x; SIM.z = out.z;
SIM.dx = out.dx; SIM.dz = out.dz;
SIM.f0 = out.f0;
SIM.cs_map = out.true_SWS;
SIM.label_map = double(out.label_map);
if isfield(out,'medium') && isfield(out.medium,'material_label_map')
    SIM.label_map = double(out.medium.material_label_map);
end
SIM.signed_distance_m = out.signed_distance_m;
SIM.abs_distance_m = out.abs_distance_m;
SIM.alpha_shear_map = out.alpha_shear_map;
SIM.alpha_us_map = out.alpha_us_map;
SIM.amp_map = abs(out.U);
SIM.phase_map = angle(out.U);
SIM.amplitude_clean = out.amplitude_clean;
SIM.shear_attenuation_factor = getfield_or(out, 'shear_attenuation_factor', ones(size(out.U)));
SIM.acoustic_readout_amplitude_factor = getfield_or(out.echo, 'echo_amplitude', ones(size(out.U)));
SIM.tracking_snr_db = out.tracking_snr_db;
SIM.readout_noise_weight = out.readout_noise_weight;
SIM.roi_masks = out.roi_masks;
SIM.source_positions_xyz = getfield_or(out, 'source_positions_xyz', out.config.Source.PositionXYZ);
SIM.source_amplitudes = getfield_or(out, 'source_amplitudes', out.config.Source.Amplitude);
SIM.source_motion_axis = getfield_or(out.config.Source, 'MotionAxis', [0 0 1]);
SIM.source_axis = getfield_or(out.config.Source, 'Axis', [0 0 1]);
SIM.measurement_axis = getfield_or(out.config.Measurement, 'Axis', [0 0 1]);
SIM.source_to_roi_min_m = source_roi_distance(SIM, CFG);
SIM.condition_key = C.condition_key;
SIM.geometry = C.geometry;
SIM.realism_level = C.realism_level;
end

function dmin = source_roi_distance(SIM, CFG)
roi = false(size(SIM.cs_map));
names = fieldnames(SIM.roi_masks);
for i = 1:numel(names)
    if contains(names{i}, 'core') || contains(names{i}, 'center')
        roi = roi | SIM.roi_masks.(names{i});
    end
end
if ~any(roi(:)), roi = true(size(SIM.cs_map)); end
[rr,cc] = find(roi);
xr = SIM.x(cc); zr = SIM.z(rr);
P = SIM.source_positions_xyz;
dmin = inf;
for k = 1:size(P,1)
    d = hypot(xr - P(k,1), zr - P(k,3));
    dmin = min(dmin, min(d,[],'omitnan'));
end
if dmin < CFG.MinSourceRoiDistanceM
    warning('Test66:sourceNearROI', 'A source is %.2f mm from a central ROI.', 1e3*dmin);
end
end

%% REQ, metadata, and model application

function F = get_or_extract_req(SIM, cfg_sim, C, M, key, CFG, OUT)
file = fullfile(OUT.req_dir, "req__" + sanitize(key) + ".mat");
if exist(file,'file') == 2
    S = load(file, 'F'); F = S.F; return;
end
if isfield(OUT, 'req_fallback_dirs') && ~isempty(OUT.req_fallback_dirs)
    for fd = string(OUT.req_fallback_dirs(:))'
        fallback_file = fullfile(fd, "req__" + sanitize(key) + ".mat");
        if exist(fallback_file, 'file') == 2
            S = load(fallback_file, 'F'); F = S.F; return;
        end
    end
end
cfg_req = struct('dx',cfg_sim.dx,'dz',cfg_sim.dz,'f0',cfg_sim.f0, ...
    'cs_bg',CFG.cs_guess,'WaveModel','simcore_clean_eikonal');
feat = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', CFG.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
step_x = max(1, round(CFG.TargetStepM / cfg_sim.dx));
step_z = max(1, round(CFG.TargetStepM / cfg_sim.dz));
O = adaptive_req.estimators.req_estimator_map(SIM.U, cfg_req, feat, ...
    'StepX', step_x, 'StepZ', step_z, ...
    'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins',CFG.REQ.Nbins, ...
    'Nbins_auto_oversample',CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min',CFG.REQ.Nbins_min,'smooth_sigma',CFG.REQ.SmoothSigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseParfor, ...
    'StoreReqCurves', false, 'Verbose', false);
F = O.feature_table;
F = attach_metadata(F, SIM, cfg_sim, C, M, key, O.win_size, CFG, step_x, step_z);
save(file, 'F', '-v7.3');
end

function F = attach_metadata(F, SIM, cfg, C, M, key, win_size, CFG, step_x, step_z)
n = height(F);
F.dataset = repmat("test66_eikonal_realistic_transfer", n, 1);
F.condition_key = repmat(string(key), n, 1);
F.geometry = repmat(C.geometry, n, 1);
F.geometry_family = repmat(C.geometry_family, n, 1);
F.case_id = F.geometry;
F.case_family = F.geometry_family;
F.realism_level = repmat(C.realism_level, n, 1);
F.field_regime = repmat(C.field_regime, n, 1);
F.field_regime_ood = F.field_regime;
F.source_seed = C.source_seed * ones(n,1);
F.noise_seed = C.noise_seed * ones(n,1);
F.f0 = cfg.f0 * ones(n,1);
F.SIM_f0 = F.f0;
F.M = M * ones(n,1);
F.REQ_M = F.M;
F.REQ_StepX = step_x * ones(n,1);
F.REQ_StepZ = step_z * ones(n,1);
F.TargetStepM = CFG.TargetStepM * ones(n,1);
F.dx = cfg.dx * ones(n,1);
F.dz = cfg.dz * ones(n,1);
F.true_SWS = nan(n,1);
F.k_true = nan(n,1);
F.patch_purity = nan(n,1);
F.patch_cs_std = nan(n,1);
F.patch_cs_range = nan(n,1);
F.patch_cs_iqr = nan(n,1);
F.q_oracle = nan(n,1);
F.distance_to_interface_m = nan(n,1);
F.distance_to_interface_mm = nan(n,1);
F.distance_to_boundary_mm = nan(n,1);
F.depth_m = nan(n,1);
F.local_amplitude = nan(n,1);
F.local_phase_rad = nan(n,1);
F.tracking_snr_db = nan(n,1);
F.snr_proxy_db = nan(n,1);
F.shear_attenuation_factor = nan(n,1);
F.acoustic_readout_amplitude_factor = nan(n,1);
F.source_to_roi_min_m = SIM.source_to_roi_min_m * ones(n,1);
F.source_to_patch_distance_m = nan(n,1);
half = floor(win_size/2);
window_radius_m = max(half,1) * mean([cfg.dx cfg.dz]);
F.window_radius_m = window_radius_m * ones(n,1);
for i = 1:n
    cx = F.cx(i); cz = F.cz(i);
    rr = max(1,cz-half):min(size(SIM.cs_map,1),cz+half);
    cc = max(1,cx-half):min(size(SIM.cs_map,2),cx+half);
    patch_cs = SIM.cs_map(rr,cc);
    patch_lab = SIM.label_map(rr,cc);
    F.true_SWS(i) = SIM.cs_map(cz, cx);
    F.k_true(i) = 2*pi*cfg.f0/F.true_SWS(i);
    F.patch_purity(i) = dominant_fraction(patch_lab);
    F.patch_cs_std(i) = std(patch_cs(:), 'omitnan');
    F.patch_cs_range(i) = max(patch_cs(:),[],'omitnan') - min(patch_cs(:),[],'omitnan');
    F.patch_cs_iqr(i) = iqr_omitnan(patch_cs(:));
    F.q_oracle(i) = invert_mapping_to_q(F.req_mapping{i}, F.k_true(i));
    F.distance_to_interface_m(i) = SIM.abs_distance_m(cz,cx);
    F.distance_to_interface_mm(i) = 1e3 * F.distance_to_interface_m(i);
    F.distance_to_boundary_mm(i) = F.distance_to_interface_mm(i);
    F.depth_m(i) = SIM.z(cz);
    F.local_amplitude(i) = SIM.amp_map(cz,cx);
    F.local_phase_rad(i) = SIM.phase_map(cz,cx);
    F.tracking_snr_db(i) = SIM.tracking_snr_db(cz,cx);
    F.snr_proxy_db(i) = F.tracking_snr_db(i);
    F.shear_attenuation_factor(i) = SIM.shear_attenuation_factor(cz,cx);
    F.acoustic_readout_amplitude_factor(i) = SIM.acoustic_readout_amplitude_factor(cz,cx);
    F.source_to_patch_distance_m(i) = min(hypot(SIM.source_positions_xyz(:,1)-SIM.x(cx), ...
        SIM.source_positions_xyz(:,3)-SIM.z(cz)), [], 'omitnan');
end
F.distance_to_interface_over_window_radius = F.distance_to_interface_m ./ F.window_radius_m;
F.is_mixed = F.patch_purity < 0.95;
F.is_strong_mixed = F.patch_purity < 0.75;
F.purity_bin = purity_bin(F.patch_purity);
F.distance_bin = distance_bin(F.distance_to_interface_mm, CFG.DistanceEdgesMm);
F.distance_over_window_bin = ratio_bin(F.distance_to_interface_over_window_radius, CFG.DistanceOverWindowEdges);
F.depth_bin = distance_bin(1e3*F.depth_m, CFG.DepthEdgesMm);
F.snr_bin = snr_bin(F.tracking_snr_db, CFG.SnrEdgesDb);
F.amplitude_bin = amplitude_bins(F.local_amplitude, CFG.AmplitudeQuantileBins);
F.roi_region = assign_roi_region(F, SIM);
F.analysis_region = F.roi_region;
F.q_theory_prior = repmat(theory_q_for_regime(cfg, M, C.field_regime), n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
end

function T_out = apply_models(T, MODELS, base_features, CFG)
P = predict_composition(MODELS.composition, T, base_features);
T.predicted_patch_purity = P.predicted_patch_purity;
T.p_mixed = P.p_mixed;
T.p_strong_mixed = P.p_strong_mixed;
model_names = CFG.ModelsToEvaluate(ismember(CFG.ModelsToEvaluate, MODELS.model_names));
parts = cell(numel(model_names),1);
for mi = 1:numel(model_names)
    name = model_names(mi);
    keep = intersect(string(T.Properties.VariableNames), ...
        ["dataset","condition_key","geometry","geometry_family","case_id","case_family", ...
        "realism_level","field_regime","field_regime_ood","source_seed","noise_seed", ...
        "f0","SIM_f0","M","REQ_M","dx","dz","REQ_StepX","REQ_StepZ","TargetStepM", ...
        "map_iz","map_ix","cx","cz","x_center_m","z_center_m","true_SWS","k_true", ...
        "patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr","purity_bin", ...
        "is_mixed","is_strong_mixed","distance_to_interface_m","distance_to_interface_mm", ...
        "distance_to_interface_over_window_radius","distance_bin","distance_over_window_bin", ...
        "depth_m","depth_bin","local_amplitude","amplitude_bin","local_phase_rad", ...
        "tracking_snr_db","snr_proxy_db","snr_bin","shear_attenuation_factor", ...
        "acoustic_readout_amplitude_factor","source_to_roi_min_m","source_to_patch_distance_m", ...
        "roi_region","analysis_region","q_oracle","q_theory_prior","sws_theory", ...
        "predicted_patch_purity","p_mixed","p_strong_mixed"], 'stable');
    R = T(:, cellstr(keep));
    R.model_name = repmat(name, height(R), 1);
    switch name
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, T);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T);
        otherwise
            continue;
    end
    sws_pred = q_to_sws(T.req_mapping, q_pred, T.f0);
    R.q_pred = q_pred;
    R.k_pred = 2*pi*T.f0 ./ sws_pred;
    R.SWS_pred = sws_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T.q_oracle;
    R.k_error = R.k_pred - T.k_true;
    R.sws_signed_error_pct = 100*(sws_pred - T.true_SWS)./T.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{mi} = R;
end
T_out = vertcat(parts{~cellfun(@isempty,parts)});
end

function F = ensure_predictor_columns(F, features, CFG)
for f = string(features(:))'
    if ~ismember(f, string(F.Properties.VariableNames))
        switch f
            case "REQ_StepX", F.REQ_StepX = max(1, round(CFG.TargetStepM ./ F.dx));
            case "REQ_StepZ", F.REQ_StepZ = max(1, round(CFG.TargetStepM ./ F.dz));
            case "TargetStepM", F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
            otherwise
                error('Required Test38 predictor missing in Test66 table: %s', f);
        end
    end
end
end

function validate_features_for_prediction(F, features)
X = F(:, cellstr(features));
for i = 1:numel(features)
    v = X.(features(i));
    if isnumeric(v)
        assert(all(isfinite(v) | isnan(v)), 'Predictor %s has Inf values.', features(i));
    end
end
end

%% Summaries and plots

function SUM = write_summaries(T, OUT)
SUM = struct();
SUM.overall = summarize_predictions(T, "model_name");
SUM.condition_averaged_overall = condition_average(T, "model_name");
SUM.by_frequency = summarize_predictions(T, ["model_name","f0"]);
SUM.by_geometry = summarize_predictions(T, ["model_name","geometry"]);
SUM.by_geometry_frequency = summarize_predictions(T, ["model_name","geometry","f0"]);
SUM.by_realism_level = summarize_predictions(T, ["model_name","realism_level"]);
SUM.by_frequency_realism = summarize_predictions(T, ["model_name","f0","realism_level"]);
SUM.by_field_regime = summarize_predictions(T, ["model_name","field_regime"]);
SUM.by_snr_bin = summarize_predictions(T, ["model_name","snr_bin"]);
SUM.by_depth_bin = summarize_predictions(T, ["model_name","depth_bin"]);
SUM.by_amplitude_bin = summarize_predictions(T, ["model_name","amplitude_bin"]);
SUM.by_patch_purity_bin = summarize_predictions(T, ["model_name","purity_bin"]);
SUM.by_distance_to_interface_bin = summarize_predictions(T, ["model_name","distance_bin"]);
SUM.by_distance_over_window_radius_bin = summarize_predictions(T, ["model_name","distance_over_window_bin"]);
SUM.by_roi = summarize_predictions(T, ["model_name","roi_region"]);
SUM.by_roi_frequency = summarize_predictions(T, ["model_name","roi_region","f0"]);
SUM.condition_summary = summarize_predictions(T, ["model_name","condition_key","geometry", ...
    "realism_level","field_regime","f0","M","source_seed","noise_seed"]);

names = fieldnames(SUM);
files = ["test66_summary_overall.csv","test66_condition_averaged_overall.csv", ...
    "test66_summary_by_frequency.csv","test66_summary_by_geometry.csv", ...
    "test66_summary_by_geometry_frequency.csv","test66_summary_by_realism_level.csv", ...
    "test66_summary_by_frequency_realism.csv", ...
    "test66_summary_by_field_regime.csv","test66_summary_by_SNR_bin.csv", ...
    "test66_summary_by_depth_bin.csv","test66_summary_by_amplitude_bin.csv", ...
    "test66_summary_by_patch_purity_bin.csv","test66_summary_by_distance_to_interface_bin.csv", ...
    "test66_summary_by_distance_over_window_radius_bin.csv","test66_summary_by_roi.csv", ...
    "test66_summary_by_roi_frequency.csv","test66_condition_summary.csv"];
for i = 1:numel(names)
    writetable(SUM.(names{i}), fullfile(OUT.table_dir, files(i)));
end
end

function [T, SUM] = read_existing_results(OUT)
patch_file = fullfile(OUT.table_dir, 'test66_patch_level_results.csv');
assert(exist(patch_file,'file') == 2, ...
    'Plot-only mode requires existing patch-level results: %s', patch_file);
T = readtable(patch_file, 'TextType','string');
SUM = struct();
SUM.overall = read_summary_table(OUT, 'test66_summary_overall.csv');
SUM.condition_averaged_overall = read_summary_table(OUT, 'test66_condition_averaged_overall.csv');
SUM.by_frequency = read_summary_table(OUT, 'test66_summary_by_frequency.csv');
SUM.by_geometry = read_summary_table(OUT, 'test66_summary_by_geometry.csv');
SUM.by_geometry_frequency = read_summary_table(OUT, 'test66_summary_by_geometry_frequency.csv');
SUM.by_realism_level = read_summary_table(OUT, 'test66_summary_by_realism_level.csv');
SUM.by_frequency_realism = read_summary_table(OUT, 'test66_summary_by_frequency_realism.csv');
SUM.by_field_regime = read_summary_table(OUT, 'test66_summary_by_field_regime.csv');
SUM.by_snr_bin = read_summary_table(OUT, 'test66_summary_by_SNR_bin.csv');
SUM.by_depth_bin = read_summary_table(OUT, 'test66_summary_by_depth_bin.csv');
SUM.by_amplitude_bin = read_summary_table(OUT, 'test66_summary_by_amplitude_bin.csv');
SUM.by_patch_purity_bin = read_summary_table(OUT, 'test66_summary_by_patch_purity_bin.csv');
SUM.by_distance_to_interface_bin = read_summary_table(OUT, 'test66_summary_by_distance_to_interface_bin.csv');
SUM.by_distance_over_window_radius_bin = read_summary_table(OUT, 'test66_summary_by_distance_over_window_radius_bin.csv');
SUM.by_roi = read_summary_table(OUT, 'test66_summary_by_roi.csv');
SUM.by_roi_frequency = read_summary_table(OUT, 'test66_summary_by_roi_frequency.csv');
SUM.condition_summary = read_summary_table(OUT, 'test66_condition_summary.csv');
end

function T = read_summary_table(OUT, name)
file = fullfile(OUT.table_dir, name);
if exist(file,'file') == 2
    T = readtable(file, 'TextType','string');
else
    T = table();
end
end

function S = summarize_predictions(T, group_vars)
if isempty(T), S = table(); return; end
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    err = X.sws_pred - X.true_SWS;
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(abs(err),'omitnan'), ...
        sqrt(mean(err.^2,'omitnan')), mean(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), 100*mean(X.high_error20,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        mean(abs(X.q_error),'omitnan'), mean(X.q_error,'omitnan'), ...
        mean(X.tracking_snr_db,'omitnan'), mean(X.local_amplitude,'omitnan'), ...
        mean(X.patch_purity,'omitnan'), mean(X.predicted_patch_purity,'omitnan'), ...
        mean(X.p_mixed,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_APE_pct','MAE_SWS', ...
        'RMSE_SWS','bias_pct','high_error10_pct','high_error20_pct', ...
        'underestimate_pct','mean_abs_q_error','mean_q_bias', ...
        'mean_tracking_snr_db','mean_amplitude','mean_patch_purity', ...
        'mean_predicted_patch_purity','mean_p_mixed'});
end
S = [groups vertcat(rows{:})];
end

function S = condition_average(T, group_vars)
C = summarize_predictions(T, [string(group_vars) "condition_key"]);
if isempty(C), S = table(); return; end
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(C(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = C(G==gi,:);
    rows{gi} = table(height(X), mean(X.MAPE_pct,'omitnan'), std(X.MAPE_pct,'omitnan'), ...
        mean(X.high_error20_pct,'omitnan'), std(X.high_error20_pct,'omitnan'), ...
        'VariableNames', {'N_conditions','condition_mean_MAPE_pct', ...
        'condition_std_MAPE_pct','condition_mean_high20_pct','condition_std_high20_pct'});
end
S = [groups vertcat(rows{:})];
end

function plot_summary_figures(T, SUM, OUT, CFG)
safe_plot(@() plot_basic_rankings(SUM, OUT), 'basic rankings');
safe_plot(@() plot_curves(SUM, OUT), 'summary curves');
safe_plot(@() plot_heatmaps(SUM, OUT), 'heatmaps');
safe_plot(@() plot_roi_figures(T, SUM, OUT), 'ROI figures');
safe_plot(@() plot_transfer_diagnostics(SUM, OUT), 'transfer diagnostics');
safe_plot(@() plot_condition_diagnostics(SUM, OUT), 'condition diagnostics');
end

function plot_basic_rankings(SUM, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
nexttile;
S = sortrows(SUM.overall, 'MAPE_pct');
barh(categorical(pretty_model(S.model_name)), S.MAPE_pct); grid on;
xlabel('MAPE (%)'); title('Global MAPE');
nexttile;
barh(categorical(pretty_model(S.model_name)), S.high_error20_pct); grid on;
xlabel('Pixels with error >20% (%)'); title('Global high-error fraction');
export_fig(fig, fullfile(OUT.figure_dir,'test66_model_ranking.png'));
end

function plot_curves(SUM, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 24]);
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
plot_group_curve(nexttile, SUM.by_frequency, 'f0', 'Frequency (Hz)', 'MAPE vs frequency');
plot_group_bar(nexttile, SUM.by_realism_level, 'realism_level', 'Realism level', 'MAPE vs realism');
plot_group_bar(nexttile, SUM.by_field_regime, 'field_regime', 'Field regime', 'MAPE vs field regime');
plot_group_bar(nexttile, SUM.by_patch_purity_bin, 'purity_bin', 'Patch purity bin', 'MAPE vs patch purity');
plot_group_bar(nexttile, SUM.by_distance_over_window_radius_bin, 'distance_over_window_bin', ...
    'Distance/window radius bin', 'MAPE vs normalized interface distance');
plot_group_bar(nexttile, SUM.by_snr_bin, 'snr_bin', 'SNR bin', 'MAPE vs SNR');
export_fig(fig, fullfile(OUT.figure_dir,'test66_summary_curves.png'));
end

function plot_heatmaps(SUM, OUT)
S = SUM.by_geometry_frequency(SUM.by_geometry_frequency.model_name=="q_spectrum_plus_composition",:);
S = aggregate_for_heatmap(S, ["geometry","f0"]);
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 12]);
draw_matrix_heatmap(gca, string(S.f0), pretty_geometry(S.geometry), S.MAPE_pct, ...
    'Frequency (Hz)', 'Geometry', ...
    'q spectrum + composition: MAPE by geometry and frequency');
export_fig(fig, fullfile(OUT.figure_dir,'test66_geometry_frequency_heatmap.png'));

if isfield(SUM, 'by_frequency_realism')
    S = SUM.by_frequency_realism(SUM.by_frequency_realism.model_name=="q_spectrum_plus_composition",:);
    if ~isempty(S)
        S = aggregate_for_heatmap(S, ["realism_level","f0"]);
        fig = figure('Color','w','Units','centimeters','Position',[1 1 22 12]);
        draw_matrix_heatmap(gca, string(S.f0), pretty_realism(S.realism_level), S.MAPE_pct, ...
            'Frequency (Hz)', 'Realism level', ...
            'q spectrum + composition: MAPE by frequency and realism');
        export_fig(fig, fullfile(OUT.figure_dir,'test66_frequency_realism_heatmap.png'));
    end
end
end

function draw_matrix_heatmap(ax, x_labels, y_labels, values, xlab, ylab, ttl)
xu = unique(string(x_labels), 'stable');
yu = unique(string(y_labels), 'stable');
M = nan(numel(yu), numel(xu));
for i = 1:numel(yu)
    for j = 1:numel(xu)
        idx = string(y_labels)==yu(i) & string(x_labels)==xu(j);
        if any(idx), M(i,j) = mean(values(idx), 'omitnan'); end
    end
end
imagesc(ax, M);
colormap(ax, parula);
set(ax, 'XTick', 1:numel(xu), 'XTickLabel', xu, ...
    'YTick', 1:numel(yu), 'YTickLabel', yu, ...
    'TickLabelInterpreter','none');
cb = colorbar(ax); ylabel(cb, 'MAPE (%)');
xlabel(ax, xlab); ylabel(ax, ylab);
title(ax, ttl, 'FontWeight','normal');
axis(ax, 'tight');
for i = 1:numel(yu)
    for j = 1:numel(xu)
        if isfinite(M(i,j))
            text(ax, j, i, sprintf('%.1f', M(i,j)), ...
                'HorizontalAlignment','center', 'FontWeight','bold', ...
                'FontSize',7, 'Color','w');
        end
    end
end
end

function S = aggregate_for_heatmap(S, group_vars)
if isempty(S), return; end
[G, groups] = findgroups(S(:, cellstr(group_vars)));
vals = splitapply(@(x) mean(x,'omitnan'), S.MAPE_pct, G);
S = [groups table(vals, 'VariableNames', {'MAPE_pct'})];
end

function plot_roi_figures(T, SUM, OUT)
S = SUM.by_roi_frequency(SUM.by_roi_frequency.model_name=="q_spectrum_plus_composition",:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on');
rois = unique(S.roi_region,'stable')';
for r = rois
    X = sortrows(S(S.roi_region==r,:), 'f0');
    if height(X) < 1, continue; end
    plot(ax, X.f0, X.MAPE_pct, '-o', 'DisplayName', pretty_region(r));
end
grid(ax,'on'); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)');
title(ax,'ROI x frequency, q spectrum + composition');
legend(ax,'Location','eastoutside');
ax = nexttile;
P = T(T.model_name=="q_spectrum_plus_composition",:);
R = summarize_predictions(P, ["roi_region","condition_key","true_SWS"]);
scatter(ax, R.true_SWS, R.MAPE_pct, 24, R.bias_pct, 'filled');
grid(ax,'on'); xlabel(ax,'True ROI/patch SWS (m/s)'); ylabel(ax,'MAPE (%)');
title(ax,'Bias-colored error by true SWS level'); colorbar(ax);
export_fig(fig, fullfile(OUT.figure_dir,'test66_roi_frequency_and_bias.png'));
end

function plot_group_curve(ax, S, xvar, xlab, ttl)
hold(ax,'on');
for m = unique(S.model_name,'stable')'
    X = sortrows(S(S.model_name==m,:), xvar);
    plot(ax, X.(xvar), X.MAPE_pct, '-o', 'DisplayName', pretty_model(m));
end
grid(ax,'on'); xlabel(ax,xlab); ylabel(ax,'MAPE (%)'); title(ax,ttl);
legend(ax,'Location','best');
end

function plot_group_bar(ax, S, xvar, xlab, ttl)
models = unique(S.model_name,'stable');
cats = unique(string(S.(xvar)),'stable');
Y = nan(numel(cats), numel(models));
for i = 1:numel(cats)
    for j = 1:numel(models)
        idx = string(S.(xvar))==cats(i) & S.model_name==models(j);
        if any(idx), Y(i,j) = S.MAPE_pct(find(idx,1)); end
    end
end
bar(ax, categorical(pretty_category(cats, xvar)), Y); grid(ax,'on');
set(ax, 'TickLabelInterpreter','none');
xlabel(ax,xlab); ylabel(ax,'MAPE (%)'); title(ax,ttl);
legend(ax, pretty_model(models), 'Location','best');
xtickangle(ax,30);
end

function plot_transfer_diagnostics(SUM, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 36 24]);
tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
plot_group_bar(nexttile, SUM.by_geometry, 'geometry', 'Geometry', 'MAPE by geometry');
plot_group_bar(nexttile, SUM.by_frequency_realism, 'realism_level', 'Realism level', 'MAPE by realism');
plot_group_bar(nexttile, SUM.by_field_regime, 'field_regime', 'Field regime', 'MAPE by field regime');
plot_group_bar(nexttile, SUM.by_depth_bin, 'depth_bin', 'Depth bin', 'MAPE by depth');
plot_group_bar(nexttile, SUM.by_amplitude_bin, 'amplitude_bin', 'Amplitude bin', 'MAPE by amplitude');
plot_group_bar(nexttile, SUM.by_snr_bin, 'snr_bin', 'Tracking SNR bin', 'MAPE by SNR');
export_fig(fig, fullfile(OUT.figure_dir,'test66_transfer_diagnostics.png'));

fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_group_bar(nexttile, SUM.by_patch_purity_bin, 'purity_bin', 'Patch purity bin', 'MAPE by patch purity');
plot_group_bar(nexttile, SUM.by_distance_over_window_radius_bin, ...
    'distance_over_window_bin', 'Distance / window radius', 'MAPE by normalized distance to interface');
export_fig(fig, fullfile(OUT.figure_dir,'test66_interface_and_purity_diagnostics.png'));
end

function plot_condition_diagnostics(SUM, OUT)
if ~isfield(SUM, 'condition_summary') || isempty(SUM.condition_summary), return; end
S = SUM.condition_summary(SUM.condition_summary.model_name=="q_spectrum_plus_composition",:);
if isempty(S), return; end
S = sortrows(S, 'MAPE_pct', 'descend');
N = min(20, height(S));
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 18]);
barh(categorical(pretty_condition(S.condition_key(1:N))), S.MAPE_pct(1:N)); grid on;
xlabel('MAPE (%)'); ylabel('Condition');
title('Worst conditions, q spectrum + composition');
export_fig(fig, fullfile(OUT.figure_dir,'test66_worst_conditions.png'));

fig = figure('Color','w','Units','centimeters','Position',[1 1 28 16]);
scatter(S.mean_tracking_snr_db, S.MAPE_pct, 36, S.bias_pct, 'filled');
grid on; xlabel('Mean tracking SNR (dB)'); ylabel('Condition MAPE (%)');
title('Condition error vs tracking SNR, q spectrum + composition');
cb = colorbar; ylabel(cb,'Bias (%)');
export_fig(fig, fullfile(OUT.figure_dir,'test66_condition_error_vs_snr.png'));
end

function plot_condition_maps(T, F, SIM, key, OUT, CFG)
primary = "q_spectrum_plus_composition";
if ~any(T.model_name == primary), primary = T.model_name(1); end
X = T(T.model_name == primary, :);
[pred_map,nz,nx] = rows_to_grid(X, X.sws_pred);
true_req_map = rows_to_grid(X, X.true_SWS, nz, nx);
pur_map = rows_to_grid(X, X.patch_purity, nz, nx);
snr_map = rows_to_grid(X, X.tracking_snr_db, nz, nx);
att_map = rows_to_grid(X, X.shear_attenuation_factor, nz, nx);
signed_map = rows_to_grid(X, X.sws_signed_error_pct, nz, nx);
abs_map = rows_to_grid(X, X.sws_abs_error_pct, nz, nx);
fig = figure('Color','w','Units','centimeters','Position',[1 1 36 22]);
tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
plot_img(SIM.cs_map, 'True SWS', 'm/s');
plot_img(SIM.label_map, 'Material label map', 'label');
plot_img(SIM.amp_map, 'Measured amplitude', 'a.u.');
plot_img(SIM.tracking_snr_db, 'Tracking SNR', 'dB');
plot_img(SIM.shear_attenuation_factor, 'Shear attenuation factor', 'relative');
plot_img(SIM.acoustic_readout_amplitude_factor, 'Acoustic/readout amplitude factor', 'relative');
plot_map_with_label(nexttile, true_req_map, 'True SWS at REQ centers', 'm/s');
plot_map_with_label(nexttile, pur_map, 'Patch purity', 'fraction');
plot_map_with_label(nexttile, pred_map, 'Predicted SWS', 'm/s');
plot_map_with_label(nexttile, signed_map, 'Signed error', '%');
plot_map_with_label(nexttile, abs_map, 'Absolute error', '%');
plot_map_with_label(nexttile, snr_map, 'REQ-center SNR', 'dB');
sgtitle(strrep(string(key), '_', '\_'));
out_dir = fullfile(OUT.map_dir, sanitize(X.geometry(1)), sanitize(X.realism_level(1)), sanitize(X.field_regime(1)));
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
export_fig(fig, fullfile(out_dir, "test66_map__" + sanitize(key) + ".png"));
field_key = field_key_from_map_key(key);
plot_field_phase_amplitude_diagnostic(SIM, field_key, out_dir);
plot_source_geometry_diagnostic(SIM, field_key, out_dir);
end

function field_key = field_key_from_map_key(key)
% Field/source diagnostics are independent of REQ M, so use a key without
% the map-specific M/step suffix. This avoids saving duplicate diagnostics
% for M=2 and M=3.
field_key = regexprep(string(key), '__M[^_]+__step[0-9]+um$', '');
end

function plot_img(A, ttl, cb)
nexttile; imagesc(A); axis image off; title(ttl, 'FontWeight','normal');
c = colorbar; ylabel(c, cb);
end

function plot_field_phase_amplitude_diagnostic(SIM, key, out_dir)
diag_dir = fullfile(out_dir, 'field_diagnostics');
if exist(diag_dir,'dir') ~= 7, mkdir(diag_dir); end
xmm = 1e3 * SIM.x;
zmm = 1e3 * SIM.z;
real_field = real(SIM.U);
phase_sign = sign(cos(angle(SIM.U)));
amp = abs(SIM.U);
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 10]);
tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile;
disable_axes_toolbar(ax);
imagesc(ax, xmm, zmm, real_field); axis(ax,'image'); set(ax,'YDir','normal');
title(ax,'Real axial harmonic field','FontWeight','normal');
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
lim = max(abs(real_field(:)), [], 'omitnan');
if isfinite(lim) && lim > 0, clim(ax, [-lim lim]); end
colormap(ax, redblue_colormap(256));
cb = colorbar(ax); ylabel(cb,'a.u.');
ax = nexttile;
disable_axes_toolbar(ax);
imagesc(ax, xmm, zmm, phase_sign); axis(ax,'image'); set(ax,'YDir','normal');
title(ax,'Phase sign: sign(cos phase)','FontWeight','normal');
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
clim(ax,[-1 1]); colormap(ax, [0.05 0.15 0.85; 0.9 0.05 0.05]);
cb = colorbar(ax); ylabel(cb,'blue=-1, red=+1');
ax = nexttile;
disable_axes_toolbar(ax);
imagesc(ax, xmm, zmm, amp); axis(ax,'image'); set(ax,'YDir','normal');
title(ax,'Harmonic amplitude |U|','FontWeight','normal');
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
colormap(ax, parula); cb = colorbar(ax); ylabel(cb,'a.u.');
sgtitle(strrep(string(key), '_', '\_'));
png_file = fullfile(diag_dir, "test66_field_phase_amplitude__" + sanitize(key) + ".png");
fig_file = fullfile(diag_dir, "test66_field_phase_amplitude__" + sanitize(key) + ".fig");
exportgraphics(fig, png_file, 'Resolution', 220);
savefig(fig, fig_file);
close(fig);
end

function plot_source_geometry_diagnostic(SIM, key, out_dir)
src_dir = fullfile(out_dir, 'source_geometry');
if exist(src_dir,'dir') ~= 7, mkdir(src_dir); end
P = double(SIM.source_positions_xyz);
if isvector(P), P = reshape(P,1,[]); end
if isempty(P) || size(P,2) < 3, return; end
amp = SIM.source_amplitudes;
if isempty(amp), amp = ones(size(P,1),1); end
amp = abs(amp(:));
if numel(amp) == 1 && size(P,1) > 1, amp = repmat(amp, size(P,1), 1); end
amp = amp ./ max(max(amp), eps);
motion = double(getfield_or(SIM, 'source_motion_axis', [0 0 1]));
if numel(motion) ~= 3 || norm(motion) <= eps, motion = [0 0 1]; end
motion = motion(:)' ./ norm(motion);
xmm = 1e3 * SIM.x;
zmm = 1e3 * SIM.z;
Pmm = 1e3 * P;
xlim_mm = [min(xmm) max(xmm)];
zlim_mm = [min(zmm) max(zmm)];
ylim_abs = max(abs([Pmm(:,2); 10]));
ylim_mm = [-ylim_abs ylim_abs];
margin = 6;
fig = figure('Color','w','Units','centimeters','Position',[1 1 42 16]);
tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile; hold(ax,'on');
disable_axes_toolbar(ax);
draw_domain_box(ax, xlim_mm, ylim_mm, zlim_mm);
scatter3(ax, Pmm(:,1), Pmm(:,2), Pmm(:,3), 70 + 90*amp, ...
    [0.95 0.15 0.05], 'filled', 'MarkerEdgeColor','k', ...
    'DisplayName','Sources');
arrow_len_mm = 6;
quiver3(ax, Pmm(:,1), Pmm(:,2), Pmm(:,3), ...
    arrow_len_mm*motion(1)*ones(size(P,1),1), ...
    arrow_len_mm*motion(2)*ones(size(P,1),1), ...
    arrow_len_mm*motion(3)*ones(size(P,1),1), ...
    0, 'Color',[0.05 0.05 0.05], 'LineWidth',1.8, ...
    'MaxHeadSize',0.9, 'DisplayName','Motion');
for i = 1:size(Pmm,1)
    text(ax, Pmm(i,1), Pmm(i,2), Pmm(i,3), sprintf('  S%d',i), ...
        'FontSize',8, 'Color','k', 'FontWeight','bold');
end
grid(ax,'on'); axis(ax,'equal');
xlim(ax, [min([xlim_mm Pmm(:,1)'])-margin, max([xlim_mm Pmm(:,1)'])+margin]);
ylim(ax, [ylim_mm(1)-margin, ylim_mm(2)+margin]);
zlim(ax, [min([zlim_mm Pmm(:,3)'])-margin, max([zlim_mm Pmm(:,3)'])+margin]);
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); zlabel(ax,'z (mm)');
title(ax,'3D source geometry','FontWeight','normal');
view(ax, 38, 20);
legend(ax,'Location','northeastoutside');

ax = nexttile; hold(ax,'on');
disable_axes_toolbar(ax);
plot(ax, [xlim_mm(1) xlim_mm(2)], [0 0], '-', ...
    'Color',[0.1 0.1 0.8], 'LineWidth',2.0, ...
    'DisplayName','Measurement plane y=0');
plot(ax, [xlim_mm(1) xlim_mm(2) xlim_mm(2) xlim_mm(1) xlim_mm(1)], ...
    [ylim_mm(1) ylim_mm(1) ylim_mm(2) ylim_mm(2) ylim_mm(1)], ...
    'k-', 'LineWidth',1.0, 'HandleVisibility','off');
in_plane = abs(Pmm(:,2)) < 1e-6;
if any(~in_plane)
    scatter(ax, Pmm(~in_plane,1), Pmm(~in_plane,2), 60 + 70*amp(~in_plane), ...
        'w', 'filled', 'MarkerEdgeColor','k', 'DisplayName','Off-plane sources');
end
if any(in_plane)
    scatter(ax, Pmm(in_plane,1), Pmm(in_plane,2), 80 + 80*amp(in_plane), ...
        [0.95 0.15 0.05], 'filled', 'MarkerEdgeColor','k', 'DisplayName','In-plane sources');
end
quiver(ax, Pmm(:,1), Pmm(:,2), ...
    arrow_len_mm*motion(1)*ones(size(P,1),1), ...
    arrow_len_mm*motion(2)*ones(size(P,1),1), ...
    0, 'Color','k', 'LineWidth',1.3, 'MaxHeadSize',0.8, ...
    'DisplayName','Motion projection');
for i = 1:size(Pmm,1)
    text(ax, Pmm(i,1), Pmm(i,2), sprintf(' S%d',i), ...
        'FontSize',8, 'Color','k', 'FontWeight','bold');
end
axis(ax,'equal'); grid(ax,'on');
xlim(ax, [min([xlim_mm Pmm(:,1)'])-margin, max([xlim_mm Pmm(:,1)'])+margin]);
ylim(ax, [ylim_mm(1)-margin, ylim_mm(2)+margin]);
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)');
title(ax,'x-y source view: plane is y=0','FontWeight','normal');
legend(ax,'Location','southoutside','Orientation','horizontal');

ax = nexttile; hold(ax,'on');
disable_axes_toolbar(ax);
imagesc(ax, xmm, zmm, SIM.cs_map); axis(ax,'image'); set(ax,'YDir','normal');
colormap(ax, parula); cb = colorbar(ax); ylabel(cb,'True SWS (m/s)');
plot(ax, [xlim_mm(1) xlim_mm(2) xlim_mm(2) xlim_mm(1) xlim_mm(1)], ...
    [zlim_mm(1) zlim_mm(1) zlim_mm(2) zlim_mm(2) zlim_mm(1)], ...
    'k-', 'LineWidth',1.2, 'HandleVisibility','off');
in_plane = abs(Pmm(:,2)) < 1e-6;
if any(~in_plane)
    scatter(ax, Pmm(~in_plane,1), Pmm(~in_plane,3), 60 + 70*amp(~in_plane), ...
        'w', 'filled', 'MarkerEdgeColor','k', 'DisplayName','Off-plane sources');
end
if any(in_plane)
    scatter(ax, Pmm(in_plane,1), Pmm(in_plane,3), 80 + 80*amp(in_plane), ...
        [0.95 0.15 0.05], 'filled', 'MarkerEdgeColor','k', 'DisplayName','In-plane sources');
end
quiver(ax, Pmm(:,1), Pmm(:,3), ...
    arrow_len_mm*motion(1)*ones(size(P,1),1), ...
    arrow_len_mm*motion(3)*ones(size(P,1),1), ...
    0, 'Color','k', 'LineWidth',1.3, 'MaxHeadSize',0.8, ...
    'DisplayName','Motion projection');
for i = 1:size(Pmm,1)
    text(ax, Pmm(i,1), Pmm(i,3), sprintf(' S%d',i), ...
        'FontSize',8, 'Color','k', 'FontWeight','bold');
end
xlim(ax, [min([xlim_mm Pmm(:,1)'])-margin, max([xlim_mm Pmm(:,1)'])+margin]);
ylim(ax, [min([zlim_mm Pmm(:,3)'])-margin, max([zlim_mm Pmm(:,3)'])+margin]);
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
title(ax,'x-z projection on true SWS','FontWeight','normal');
legend(ax,'Location','southoutside','Orientation','horizontal');
sgtitle(strrep(string(key), '_', '\_'));
png_file = fullfile(src_dir, "test66_sources_3d__" + sanitize(key) + ".png");
fig_file = fullfile(src_dir, "test66_sources_3d__" + sanitize(key) + ".fig");
exportgraphics(fig, png_file, 'Resolution', 220);
savefig(fig, fig_file);
close(fig);
end

function draw_domain_box(ax, xlim_mm, ylim_mm, zlim_mm)
x = xlim_mm; y = ylim_mm; z = zlim_mm;
corners = [
    x(1) y(1) z(1); x(2) y(1) z(1); x(2) y(2) z(1); x(1) y(2) z(1);
    x(1) y(1) z(2); x(2) y(1) z(2); x(2) y(2) z(2); x(1) y(2) z(2)];
edges = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
for e = 1:size(edges,1)
    p = corners(edges(e,:),:);
    plot3(ax, p(:,1), p(:,2), p(:,3), 'Color',[0.25 0.25 0.25], ...
        'LineWidth',1.0, 'HandleVisibility','off');
end
patch(ax, [x(1) x(2) x(2) x(1)], [0 0 0 0], [z(1) z(1) z(2) z(2)], ...
    [0.85 0.9 1.0], 'FaceAlpha',0.22, 'EdgeColor',[0.1 0.1 0.7], ...
    'LineWidth',1.0, 'DisplayName','Measurement plane y=0');
plot3(ax, [x(1) x(2)], [0 0], [z(1) z(1)], '-', ...
    'Color',[0.05 0.05 0.8], 'LineWidth',2.0, 'HandleVisibility','off');
plot3(ax, [x(1) x(2)], [0 0], [z(2) z(2)], '-', ...
    'Color',[0.05 0.05 0.8], 'LineWidth',2.0, 'HandleVisibility','off');
text(ax, mean(x), 0, z(2), '  y=0', ...
    'Color',[0.05 0.05 0.8], 'FontWeight','bold', 'FontSize',9, ...
    'HandleVisibility','off');
end

function print_summary(T_overall, T_by_level, T_by_roi, OUT)
fprintf('\nInterpretive summary:\n');
if isempty(T_overall), return; end
T = sortrows(T_overall, 'MAPE_pct');
fprintf('  Best global model: %s (MAPE %.2f%%, high>20 %.2f%%).\n', ...
    T.model_name(1), T.MAPE_pct(1), T.high_error20_pct(1));
P = T_by_level(T_by_level.model_name=="q_spectrum_plus_composition",:);
if ~isempty(P)
    [~,ii] = max(P.MAPE_pct);
    fprintf('  Worst realism level for primary model: %s (MAPE %.2f%%).\n', ...
        P.realism_level(ii), P.MAPE_pct(ii));
end
R = T_by_roi(T_by_roi.model_name=="q_spectrum_plus_composition",:);
if ~isempty(R)
    [~,ii] = max(R.MAPE_pct);
    fprintf('  Worst ROI for primary model: %s (MAPE %.2f%%).\n', ...
        R.roi_region(ii), R.MAPE_pct(ii));
end
fprintf('  Outputs: %s\n', OUT.root_dir);
end

%% Runtime estimate

function RTE = estimate_runtime(CONDITIONS, BUNDLE, BASE_FEATURES, CFG, OUT)
benchN = min(numel(CONDITIONS), 1);
t0 = tic;
for i = 1:benchN
    [SIM, sim_cfg] = get_or_build_field(CONDITIONS(i), CFG, OUT);
    key = sprintf('%s__benchmark_M%g', CONDITIONS(i).condition_key, CFG.M(1));
    F = get_or_extract_req(SIM, sim_cfg, CONDITIONS(i), CFG.M(1), key, CFG, OUT);
    F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
    apply_models(F, BUNDLE.MODELS, BASE_FEATURES, CFG);
end
sec_per_M = toc(t0) / max(benchN,1);
nM = numel(CFG.M);
est_quick = sec_per_M * numel(CONDITIONS) * nM / 60;
RTE = table(numel(CONDITIONS), nM, sec_per_M, est_quick, ...
    'VariableNames', {'N_conditions','N_M','benchmark_sec_per_condition_M', ...
    'estimated_minutes_for_selected_run'});
fprintf('Runtime estimate: %.1f sec per condition/M, %.1f min for selected run.\n', ...
    sec_per_M, est_quick);
end

%% Shared helpers

function assert_no_forbidden_predictors(features)
bad_patterns = ["true","oracle","purity","mixed","confidence","error", ...
    "pred","sws","cs_","k_true","q_local","q_pred","q_theory", ...
    "req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition"];
features = lower(string(features));
for p = bad_patterns
    hit = features(contains(features, p));
    assert(isempty(hit), 'Forbidden predictor detected: %s', strjoin(hit, ', '));
end
fprintf('Base predictors passed leakage guard (%d predictors).\n', numel(features));
end

function y = predict_q_model(M, T)
y = clamp01(predict(M.model, T(:, cellstr(M.features))));
end

function COMP = predict_composition(MIX, F, features)
X = F(:, cellstr(features));
COMP.predicted_patch_purity = min(max(predict(MIX.purity, X), 0), 1);
[~, score] = predict(MIX.mixed, X); COMP.p_mixed = positive_score(MIX.mixed, score);
[~, score] = predict(MIX.strong_mixed, X); COMP.p_strong_mixed = positive_score(MIX.strong_mixed, score);
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
if contains(string(regime), "single_source")
    type = "SingleWave";
elseif contains(string(regime), "diffuse_like")
    type = "Diffuse3D";
else
    type = "SingleWave";
end
try
    out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, feat.cs_guess_used, ...
        'M', feat.M, 'Gamma', feat.gamma_win, 'PadFactor', feat.pad_factor, ...
        'Nbins', 'auto', 'SmoothSigma', 1, 'TheoryMode', 'S2D', ...
        'FieldType', type, 'Plot', false);
    q = out.q_th;
catch
    q = 0.5;
end
q = clamp01(q);
end

function p = dominant_fraction(patch)
v = double(patch(:));
v = v(isfinite(v));
if isempty(v), p = NaN; return; end
[~,~,ic] = unique(v);
p = max(accumarray(ic,1)) / numel(v);
end

function y = iqr_omitnan(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = prctile(x,75) - prctile(x,25); end
end

function labels = purity_bin(p)
labels = repmat("strongly_mixed", size(p));
labels(p >= 0.50) = "mixed_0p50_0p75";
labels(p >= 0.75) = "mixed_0p75_0p90";
labels(p >= 0.90) = "near_pure_0p90_0p95";
labels(p >= 0.95) = "near_pure_0p95_0p99";
labels(p >= 0.99) = "pure_0p99_1p00";
labels(~isfinite(p)) = "unknown";
end

function labels = distance_bin(d, edges)
labels = strings(size(d));
for i = 1:numel(edges)-1
    idx = d >= edges(i) & d < edges(i+1);
    if isinf(edges(i+1))
        lab = sprintf('>%gmm', edges(i));
    else
        lab = sprintf('%g_%gmm', edges(i), edges(i+1));
    end
    labels(idx) = string(lab);
end
labels(~isfinite(d)) = "no_interface";
end

function labels = ratio_bin(r, edges)
labels = strings(size(r));
for i = 1:numel(edges)-1
    idx = r >= edges(i) & r < edges(i+1);
    if isinf(edges(i+1))
        lab = sprintf('>%g_window_radius', edges(i));
    else
        lab = sprintf('%g_%g_window_radius', edges(i), edges(i+1));
    end
    labels(idx) = string(strrep(lab,'.','p'));
end
labels(~isfinite(r)) = "no_interface";
end

function labels = snr_bin(s, edges)
labels = strings(size(s));
for i = 1:numel(edges)-1
    idx = s >= edges(i) & s < edges(i+1);
    labels(idx) = sprintf('%g_%g_dB', edges(i), edges(i+1));
end
labels(~isfinite(s)) = "unknown";
end

function labels = amplitude_bins(a, nb)
labels = repmat("unknown", size(a));
ok = isfinite(a);
if sum(ok) < nb, return; end
edges = quantile(a(ok), linspace(0,1,nb+1));
edges(1) = -Inf; edges(end) = Inf;
idx = discretize(a, edges);
for i = 1:nb
    labels(idx==i) = sprintf('amplitude_Q%d', i);
end
end

function roi = assign_roi_region(F, SIM)
roi = repmat("other", height(F), 1);
names = fieldnames(SIM.roi_masks);
for j = 1:numel(names)
    mask = SIM.roi_masks.(names{j});
    inside = mask(sub2ind(size(mask), F.cz, F.cx));
    roi(inside) = string(names{j});
end
interface = isfinite(F.distance_to_interface_mm);
roi(interface & F.distance_to_interface_mm >= 0 & F.distance_to_interface_mm < 1) = "interface_0_1mm";
roi(interface & F.distance_to_interface_mm >= 1 & F.distance_to_interface_mm < 2) = "interface_1_2mm";
roi(interface & F.distance_to_interface_mm >= 2 & F.distance_to_interface_mm < 4 & roi=="other") = "interface_2_4mm";
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
if isempty(T), return; end
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map_with_label(ax, Z, ttl, label)
imagesc(ax, Z, 'AlphaData', isfinite(Z));
axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
ax.Color = [0.94 0.94 0.94];
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
cb = colorbar(ax); ylabel(cb, label);
end

function disable_axes_toolbar(ax)
try
    axtoolbar(ax, {});
catch
end
try
    ax.Toolbar.Visible = 'off';
catch
end
try
    disableDefaultInteractivity(ax);
catch
end
end

function cmap = redblue_colormap(n)
if nargin < 1, n = 256; end
n1 = floor(n/2);
n2 = n - n1;
blue = [linspace(0.05,1,n1)' linspace(0.15,1,n1)' linspace(0.85,1,n1)'];
red = [linspace(1,0.9,n2)' linspace(1,0.05,n2)' linspace(1,0.05,n2)'];
cmap = [blue; red];
end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
end

function y = clamp01(y)
y = min(max(y, 0), 1);
end

function v = getfield_or(S, name, default_value)
if isstruct(S) && isfield(S, name), v = S.(name); else, v = default_value; end
end

function s = pretty_model(x)
x = string(x);
s = x;
s(x=="q_spectrum_plus_composition") = "q spectrum + composition";
s(x=="q_spectrum_only") = "q spectrum only";
s = strrep(s, "_", " ");
end

function s = pretty_region(x)
s = string(x);
s = replace_many(s, [
    "homogeneous_center", "homogeneous center";
    "background_far", "background far";
    "soft_core", "soft core";
    "hard_core", "hard core";
    "hard_inclusion_core", "hard inclusion core";
    "intermediate_layer_core", "intermediate layer core";
    "interface_0_1mm", "interface 0-1 mm";
    "interface_1_2mm", "interface 1-2 mm";
    "interface_2_4mm", "interface 2-4 mm";
    "other", "other"]);
end

function s = pretty_geometry(x)
s = string(x);
s = replace_many(s, [
    "homogeneous_cs2", "homogeneous cs=2 m/s";
    "homogeneous_cs3", "homogeneous cs=3 m/s";
    "homogeneous_cs4", "homogeneous cs=4 m/s";
    "bilayer_2_3", "bilayer 2/3 m/s";
    "inclusion_2_3", "inclusion 2/3 m/s";
    "inclusion_2_4", "inclusion 2/4 m/s";
    "bilayer_inclusion_2_3_4", "bilayer + inclusion 2/3/4 m/s"]);
end

function s = pretty_realism(x)
s = string(x);
s = replace_many(s, [
    "clean", "clean";
    "shear_attenuation_only", "shear attenuation only";
    "readout_soft", "readout soft";
    "readout_medium", "readout medium";
    "readout_hard", "readout hard"]);
end

function s = pretty_regime(x)
s = string(x);
s = replace_many(s, [
    "single_source_lateral", "single lateral source";
    "single_source_diagonal", "single diagonal source";
    "diffuse_like_8src", "diffuse-like, 8 sources";
    "diffuse_like_16src", "diffuse-like, 16 sources"]);
end

function s = pretty_bin(x)
s = string(x);
s = replace_many(s, [
    "pure_0p99_1p00", "pure 0.99-1.00";
    "near_pure_0p95_0p99", "near-pure 0.95-0.99";
    "near_pure_0p90_0p95", "near-pure 0.90-0.95";
    "mixed_0p75_0p90", "mixed 0.75-0.90";
    "mixed_0p50_0p75", "mixed 0.50-0.75";
    "strongly_mixed", "strongly mixed";
    "0_0p25_window_radius", "0-0.25 window radii";
    "0p25_0p5_window_radius", "0.25-0.5 window radii";
    "0p5_1_window_radius", "0.5-1 window radii";
    "1_2_window_radius", "1-2 window radii";
    ">2_window_radius", ">2 window radii";
    "no_interface", "no interface";
    "0_10mm", "0-10 mm";
    "10_20mm", "10-20 mm";
    "20_30mm", "20-30 mm";
    "30_40mm", "30-40 mm";
    ">40mm", ">40 mm";
    "amplitude_Q1", "amplitude Q1";
    "amplitude_Q2", "amplitude Q2";
    "amplitude_Q3", "amplitude Q3";
    "amplitude_Q4", "amplitude Q4";
    "amplitude_Q5", "amplitude Q5";
    "-Inf_5_dB", "<5 dB";
    "5_10_dB", "5-10 dB";
    "10_15_dB", "10-15 dB";
    "15_20_dB", "15-20 dB";
    "20_30_dB", "20-30 dB";
    "30_40_dB", "30-40 dB";
    "40_Inf_dB", ">40 dB"]);
s = strrep(s, "_window_radius", " window radii");
s = strrep(s, "window_radius", "window radii");
s = strrep(s, "_dB", " dB");
s = strrep(s, "mm", " mm");
s = strrep(s, "_", "-");
s = strrep(s, "-Inf", "-Inf");
s = strrep(s, "Inf", "Inf");
end

function s = pretty_category(x, xvar)
x = string(x);
switch string(xvar)
    case "geometry"
        s = pretty_geometry(x);
    case "realism_level"
        s = pretty_realism(x);
    case "field_regime"
        s = pretty_regime(x);
    case "roi_region"
        s = pretty_region(x);
    otherwise
        s = pretty_bin(x);
end
end

function s = pretty_condition(x)
s = string(x);
s = strrep(s, "__", " | ");
s = strrep(s, "_", " ");
for old = ["homogeneous cs2","homogeneous cs4","inclusion 2 4", ...
        "bilayer inclusion 2 3 4", "single source lateral", ...
        "diffuse like 8src", "readout medium"]
    new = old;
    new = strrep(new, "cs", "c_s=");
    s = strrep(s, old, new);
end
end

function s = replace_many(s, pairs)
for i = 1:size(pairs,1)
    s(s == pairs(i,1)) = pairs(i,2);
end
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

function write_config_json(CFG, path)
folder = fileparts(path);
if exist(folder,'dir') ~= 7, mkdir(folder); end
try
    txt = jsonencode(CFG, 'PrettyPrint', true);
catch
    txt = jsonencode(CFG);
end
fid = fopen(path, 'w');
assert(fid > 0, 'Could not open config JSON for writing: %s', path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', txt);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Test66:plotFailed', 'Could not create %s: %s', label, ME.message);
end
end

function s = env_string(name, default)
s = string(getenv(name));
if s == "", s = string(default); end
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

function xs = env_number_list(name, default)
v = strtrim(string(getenv(name)));
if v == ""
    xs = default;
    return;
end
tokens = strtrim(split(v, {',',';',' '})).';
tokens(tokens=="") = [];
xs = str2double(tokens);
xs = xs(isfinite(xs));
if isempty(xs)
    xs = default;
end
end

function xs = env_string_list(name, default)
v = strtrim(string(getenv(name)));
if v == "", xs = string(default); return; end
xs = strtrim(split(v, {',',';'})).';
xs(xs=="") = [];
end
