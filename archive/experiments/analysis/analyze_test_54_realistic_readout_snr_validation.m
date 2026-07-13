%% analyze_test_54_realistic_readout_snr_validation.m
% Test 54: realistic ultrasound/readout/SNR stress-test for frozen q models.
%
% The clean paper model is not retrained here. This script generates
% ultrasound-like 2D/2.5D fields with /Users/sara/Documents/wave_sim_project,
% re-extracts REQ using the Test38 profile, applies a frozen Test38/Test53
% q-model bundle, and reports degradation from clean to moderate/hard readout.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST54_MODE          = validate | quick | full
%   ADAPTIVE_REQ_TEST54_WAVE_SIM_PATH = /Users/sara/Documents/wave_sim_project
%   ADAPTIVE_REQ_TEST54_MODEL_SOURCE  = full | test53 | /path/to/bundle.mat
%   ADAPTIVE_REQ_TEST54_TARGET_STEP_M = 0.0005
%   ADAPTIVE_REQ_TEST54_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST54_USE_PARFOR    = true | false

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
write_config_json(CFG, fullfile(OUT.root_dir, 'test54_configuration.json'));

fprintf('\nTest 54: realistic readout/SNR validation for frozen q models\n');
fprintf('Mode: %s | wave_sim_project: %s\n', CFG.Mode, CFG.WaveSimPath);
fprintf('Model source: %s\n', CFG.ModelBundleFile);
fprintf('REQ profile: Nbins=%s, oversample=%g, min=%g, smooth_sigma=%g\n', ...
    string(CFG.REQ.Nbins), CFG.REQ.Nbins_auto_oversample, ...
    CFG.REQ.Nbins_min, CFG.REQ.SmoothSigma);
fprintf('No q model is trained. Truth/noise diagnostics are evaluation-only.\n');

assert(exist(fullfile(CFG.WaveSimPath, 'src'), 'dir') == 7, ...
    'Missing wave_sim_project src folder: %s', CFG.WaveSimPath);
addpath(fullfile(CFG.WaveSimPath, 'src'));

if CFG.ValidateOnly
    C = build_conditions(CFG);
    assert(~isempty(C), 'No Test54 conditions were built.');
    sim_cfg = build_wave_sim_cfg(C(1), CFG);
    out = simcore.simulateTimeReadout2D(sim_cfg);
    assert(isfield(out,'U') && isfield(out,'cs_map'));
    fprintf('Validation-only generated one field: %s, size %s.\n', ...
        C(1).condition_key, mat2str(size(out.U)));
    return;
end

BUNDLE = load_model_bundle(CFG);
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);

CONDITIONS = build_conditions(CFG);
fprintf('Realistic readout conditions: %d | M values: %s\n', ...
    numel(CONDITIONS), mat2str(CFG.M));

parts = {};
for ci = 1:numel(CONDITIONS)
    C = CONDITIONS(ci);
    fprintf('\n[%d/%d] Generating %s\n', ci, numel(CONDITIONS), C.condition_key);
    t0 = tic;
    [SIM, sim_cfg] = get_or_build_field(C, CFG, OUT);
    fprintf('  field size %s | effective frequency SNR %.2f dB | %.1f s\n', ...
        mat2str(size(SIM.U)), get_snr_db(SIM), toc(t0));
    for mi = 1:numel(CFG.M)
        M = CFG.M(mi);
        key = sprintf('%s__M%g__step%gum', C.condition_key, M, round(1e6*CFG.TargetStepM));
        t1 = tic;
        F = get_or_extract_req(SIM, sim_cfg, C, M, key, CFG, OUT);
        F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
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
T_overall = summarize_predictions(T_all, "model_name");
T_by_level = summarize_predictions(T_all, ["model_name","realism_level"]);
T_by_case = summarize_predictions(T_all, ["model_name","case_id","realism_level"]);
T_by_freq = summarize_predictions(T_all, ["model_name","f0","realism_level"]);
T_by_M = summarize_predictions(T_all, ["model_name","M","realism_level"]);
T_by_region = summarize_predictions(T_all, ["model_name","analysis_region","realism_level"]);
T_degrade = degradation_table(T_by_level);

writetable(T_all, fullfile(OUT.table_dir, 'test54_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test54_summary_overall.csv'));
writetable(T_by_level, fullfile(OUT.table_dir, 'test54_summary_by_realism_level.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test54_summary_by_case_level.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir, 'test54_summary_by_frequency_level.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test54_summary_by_M_level.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test54_summary_by_region_level.csv'));
writetable(T_degrade, fullfile(OUT.table_dir, 'test54_clean_to_realistic_degradation.csv'));

plot_summary(T_by_level, T_degrade, T_by_region, OUT);
save(fullfile(OUT.data_dir, 'test54_compact_results.mat'), ...
    'CFG','T_overall','T_by_level','T_by_case','T_by_freq','T_by_M', ...
    'T_by_region','T_degrade','-v7.3');

print_summary(T_overall, T_degrade, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 54 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST54_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["validate","quick","full"]), ...
    'ADAPTIVE_REQ_TEST54_MODE must be validate, quick, or full.');
CFG = struct();
CFG.Mode = mode;
CFG.ValidateOnly = mode == "validate" || env_true('ADAPTIVE_REQ_TEST54_VALIDATE_ONLY', false);
CFG.QuickMode = mode == "quick";
CFG.FullMode = mode == "full";
CFG.WaveSimPath = char(env_string('ADAPTIVE_REQ_TEST54_WAVE_SIM_PATH', ...
    "/Users/sara/Documents/wave_sim_project"));
CFG.ModelSource = env_string('ADAPTIVE_REQ_TEST54_MODEL_SOURCE', "full");
CFG.ModelBundleFile = resolve_model_bundle(root_dir, CFG.ModelSource);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST54_SAVE_ALL_MAPS', true);
CFG.UseParfor = env_true('ADAPTIVE_REQ_TEST54_USE_PARFOR', false);
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST54_TARGET_STEP_M', 0.5e-3);
CFG.dx = 0.2e-3; CFG.dz = 0.2e-3;
CFG.Lx = 0.05; CFG.Lz = 0.05;
CFG.cs_guess = 3.0;
CFG.RandomSeed = 54001;
CFG.M = [2 3];
if CFG.QuickMode || CFG.ValidateOnly
    CFG.M = 2;
end
if CFG.ValidateOnly
    CFG.dx = 0.5e-3; CFG.dz = 0.5e-3;
    CFG.Lx = 0.02; CFG.Lz = 0.02;
end
CFG.ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST54_MODELS', ...
    ["q_spectrum_only","q_spectrum_plus_composition","q_spectrum_plus_theory_composition"]);
CFG.REQ = struct('Gamma',1,'PadFactor',1,'EdgeMode',"valid", ...
    'Nbins',"auto",'Nbins_auto_oversample',1,'Nbins_min',16,'SmoothSigma',1);
CFG.DistanceEdgesMm = [0 1 2 4 8 Inf];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.FullMode
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_54_realistic_readout_snr_validation');
else
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_54_realistic_readout_snr_validation', char(CFG.Mode));
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
end

function file = resolve_model_bundle(root_dir, source)
src = string(source);
if exist(src, 'file') == 2
    file = char(src); return;
end
switch lower(src)
    case "test53"
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
fprintf('Loaded frozen bundle with %d model entries.\n', numel(B.MODELS.model_names));
end

%% Conditions and wave_sim_project configs

function C = build_conditions(CFG)
cases = [
    case_spec("homogeneous_cs2", "homogeneous", 2, 2)
    case_spec("inclusion_2_3", "inclusion", 2, 3)
    case_spec("bilayer_2_3", "bilayer", 2, 3)
    ];
if CFG.QuickMode || CFG.ValidateOnly
    cases = cases(1:2);
    freqs = [300 500];
else
    freqs = [200 300 400 500 600];
end
levels = ["clean","moderate_realistic","hard_realistic"];
C = repmat(empty_condition(), 0, 1);
for ci = 1:numel(cases)
    for fi = 1:numel(freqs)
        for li = 1:numel(levels)
            X = empty_condition();
            X.case_id = cases(ci).case_id;
            X.case_family = cases(ci).case_family;
            X.cs_bg = cases(ci).cs_bg;
            X.cs_inc = cases(ci).cs_inc;
            X.f0 = freqs(fi);
            X.realism_level = levels(li);
            X.condition_key = sprintf('%s__f%g__%s', X.case_id, X.f0, X.realism_level);
            C(end+1,1) = X; %#ok<AGROW>
        end
    end
end
end

function C = case_spec(id, family, cs_bg, cs_inc)
C = struct('case_id', string(id), 'case_family', string(family), ...
    'cs_bg', cs_bg, 'cs_inc', cs_inc);
end

function X = empty_condition()
X = struct('case_id',"",'case_family',"",'cs_bg',NaN,'cs_inc',NaN, ...
    'f0',NaN,'realism_level',"",'condition_key',"");
end

function cfg = build_wave_sim_cfg(C, CFG)
cfg = struct();
cfg.Lx = CFG.Lx; cfg.Lz = CFG.Lz; cfg.dx = CFG.dx; cfg.dz = CFG.dz;
cfg.f0 = C.f0; cfg.cs_bg = C.cs_bg;
cfg.WaveModel = 'spherical';
cfg.Nwaves = 8;
if CFG.ValidateOnly, cfg.Nwaves = 1; end
cfg.Is2D = true;
cfg.SourceSampling = 'cone';
cfg.AngularSamplingMethod = 'fibonacci';
cfg.ConeAxis = [-1 0 0];
cfg.ConeHalfAngleDeg = 30;
cfg.ForceInPlaneWave = true;
cfg.AngleRange2D = [0 0];
cfg.DecayAlpha = 0.5;
cfg.SNR = Inf;
cfg.Seed = CFG.RandomSeed + round(C.f0) + sum(double(char(C.condition_key)));
cfg.Measurement.Mode = 'projected';
cfg.Measurement.Axis = [0 0 1];
cfg.Measurement.Polarization = 'inplane_sv';
cfg.Acquisition.Nt = 30;
cfg.Acquisition.PRF = max(4000, 8*C.f0);
cfg.Acquisition.StoreTimeSeries = false;
cfg.Medium = build_medium_config(C);
cfg = apply_realism_level(cfg, C.realism_level);
end

function M = build_medium_config(C)
M = struct();
M.cs_bg = C.cs_bg;
M.Material = material(C.cs_bg, 3, 0.5, 0);
switch C.case_family
    case "homogeneous"
        M.Masks = {};
    case "inclusion"
        M.Masks = {struct('Type','circle','Params',struct( ...
            'Center',[0.025 0.025],'Radius',0.008,'SigmaEdge',0), ...
            'Material',material(C.cs_inc, 30, 1.5, 6))};
    case "bilayer"
        M.Masks = {struct('Type','bilayer','Params',struct( ...
            'Bi_Angle',0,'Bi_Offset',0.025,'SigmaEdge',0), ...
            'Material',material(C.cs_inc, 30, 1.5, 6))};
    otherwise
        M.Masks = {};
end
end

function S = material(c, alpha_np_m, us_alpha, backscatter_db)
S = struct('Type','isotropic','c',c, ...
    'Attenuation',struct('Model','power_law','alpha0',alpha_np_m, ...
    'Units','Np/m','fRef',400,'exponent',1.0), ...
    'USAttenuation',struct('alpha0',us_alpha,'Units','dB/cm/MHz'), ...
    'Backscatter',struct('GainDB',backscatter_db));
end

function cfg = apply_realism_level(cfg, level)
cfg.Shear.Attenuation.Enabled = false;
cfg.Shear.Attenuation.Mode = 'ray_integral';
cfg.Shear.Attenuation.Model = 'none';
cfg.Shear.Attenuation.RaySamples = 48;
cfg.Ultrasound.FrequencyMHz = 7;
cfg.Ultrasound.AcousticAttenuation.Enabled = false;
cfg.Ultrasound.ReadoutNoise.Enabled = false;
switch string(level)
    case "clean"
        return;
    case "moderate_realistic"
        cfg.Shear.Attenuation.Enabled = true;
        cfg.Ultrasound.ReadoutNoise.Enabled = true;
        cfg.Ultrasound.ReadoutNoise.Model = 'boukraa_time_domain';
        cfg.Ultrasound.ReadoutNoise.Level = 'medium';
        cfg.Ultrasound.ReadoutNoise.DepthProfile = 'exponential_0_to_1';
        cfg.Ultrasound.ReadoutNoise.DepthExponent = 3;
        cfg.Ultrasound.ReadoutNoise.CompensateTemporalExtraction = true;
        cfg.Ultrasound.ReadoutNoise.FrequencyNoiseGain = 1.0;
        cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.Enabled = true;
        cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaX_mm = 0.45;
        cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaZ_mm = 0.45;
    case "hard_realistic"
        cfg = apply_realism_level(cfg, "moderate_realistic");
        cfg.Ultrasound.ReadoutNoise.Level = 'high';
        cfg.Ultrasound.ReadoutNoise.FrequencyNoiseGain = 1.5;
        cfg.Ultrasound.AcousticAttenuation.Enabled = true;
        cfg.Ultrasound.AcousticAttenuation.Mode = 'beam_integral';
        cfg.Ultrasound.AcousticAttenuation.FrequencyMHz = 7;
        cfg.Ultrasound.AcousticAttenuation.TwoWay = true;
        cfg.Ultrasound.AcousticAttenuation.BeamDirection = [0 0 1];
        cfg.Ultrasound.AcousticAttenuation.CouplingToReadoutNoise = true;
        cfg.Ultrasound.AcousticAttenuation.IncludeBackscatter = true;
        cfg.Ultrasound.AcousticAttenuation.MaxNoiseGain = 8;
        cfg.Ultrasound.AcousticAttenuation.MinNoiseGain = 0.25;
    otherwise
        error('Unknown realism level: %s', level);
end
end

function [SIM, cfg] = get_or_build_field(C, CFG, OUT)
file = fullfile(OUT.field_dir, "field__" + sanitize(C.condition_key) + ".mat");
if exist(file,'file') == 2
    S = load(file, 'SIM', 'cfg');
    SIM = S.SIM; cfg = S.cfg; return;
end
cfg = build_wave_sim_cfg(C, CFG);
out = simcore.simulateTimeReadout2D(cfg);
SIM = struct();
SIM.U = out.U;
SIM.Uxz = out.U;
SIM.U_clean = getfield_or(out, 'U_clean', out.U); %#ok<GFLD>
SIM.x = out.x; SIM.z = out.z;
SIM.cs_map = out.cs_map;
SIM.k_map = out.k_map;
SIM.alpha_map = getfield_or(out, 'alpha_map', zeros(size(out.cs_map)));
SIM.amp_map = abs(out.U);
SIM.phase_map = angle(out.U);
SIM.diag = out.diag;
save(file, 'SIM', 'cfg', 'C', '-v7.3');
end

%% REQ, metadata, and model application

function F = get_or_extract_req(SIM, cfg_sim, C, M, key, CFG, OUT)
file = fullfile(OUT.req_dir, "req__" + sanitize(key) + ".mat");
if exist(file,'file') == 2
    S = load(file, 'F'); F = S.F; return;
end
cfg_req = struct('dx',cfg_sim.dx,'dz',cfg_sim.dz,'f0',cfg_sim.f0, ...
    'cs_bg',C.cs_bg,'WaveModel','wave_sim_project');
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
F.dataset = repmat("test54_realistic_readout", n, 1);
F.condition_key = repmat(string(key), n, 1);
F.case_id = repmat(C.case_id, n, 1);
F.case_family = repmat(C.case_family, n, 1);
F.realism_level = repmat(C.realism_level, n, 1);
F.field_regime = repmat("directional_2D", n, 1);
F.field_regime_ood = repmat("wave_sim_project_directional_side", n, 1);
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
    F.true_SWS(i) = SIM.cs_map(cz, cx);
    F.k_true(i) = 2*pi*cfg.f0/F.true_SWS(i);
    patch = SIM.cs_map((cz-half):(cz+half), (cx-half):(cx+half));
    F.patch_purity(i) = dominant_fraction(patch);
    F.q_oracle(i) = invert_mapping_to_q(F.req_mapping{i}, F.k_true(i));
end
F.is_mixed = F.patch_purity < 0.95;
F.is_strong_mixed = F.patch_purity < 0.75;
F.purity_bin = purity_bin(F.patch_purity);
D = distance_to_boundary(SIM.cs_map, cfg.dx);
F.distance_to_boundary_mm = D(sub2ind(size(D), F.cz, F.cx));
F.distance_bin = distance_bin(F.distance_to_boundary_mm, CFG.DistanceEdgesMm);
F.q_theory_prior = repmat(theory_q_for_regime(cfg, M, "directional_2D"), n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
F.local_amplitude = SIM.amp_map(sub2ind(size(SIM.amp_map), F.cz, F.cx));
F.local_phase_rad = SIM.phase_map(sub2ind(size(SIM.phase_map), F.cz, F.cx));
F.frequency_snr_db = get_snr_db(SIM) * ones(n,1);
F.analysis_region = analysis_region(F);
end

function T_out = apply_models(T, MODELS, base_features, CFG)
P = predict_composition(MODELS.composition, T, base_features);
T.predicted_patch_purity = P.predicted_patch_purity;
T.p_mixed = P.p_mixed;
T.p_strong_mixed = P.p_strong_mixed;
model_names = MODELS.model_names(ismember(MODELS.model_names, CFG.ModelsToEvaluate));
parts = cell(numel(model_names),1);
for mi = 1:numel(model_names)
    name = model_names(mi);
    keep = intersect(string(T.Properties.VariableNames), ...
        ["dataset","condition_key","case_id","case_family","realism_level", ...
        "field_regime","field_regime_ood","f0","M","dx","dz","map_iz","map_ix", ...
        "cx","cz","x_center_m","z_center_m","true_SWS","k_true","patch_purity", ...
        "purity_bin","is_mixed","is_strong_mixed","distance_to_boundary_mm", ...
        "distance_bin","analysis_region","q_oracle","q_theory_prior","sws_theory", ...
        "local_amplitude","local_phase_rad","frequency_snr_db", ...
        "predicted_patch_purity","p_mixed","p_strong_mixed"], 'stable');
    R = T(:, cellstr(keep));
    R.model_name = repmat(name, height(R), 1);
    switch name
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, T);
            sws_pred = q_to_sws(T.req_mapping, q_pred, T.f0);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T);
            sws_pred = q_to_sws(T.req_mapping, q_pred, T.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory_composition, T);
            sws_pred = q_to_sws(T.req_mapping, q_pred, T.f0);
        otherwise
            continue;
    end
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.q_error = q_pred - T.q_oracle;
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
                error('Required predictor missing in Test54 table: %s', f);
        end
    end
end
end

%% Summaries and plots

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), 100*mean(X.high_error20,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), mean(abs(X.q_error),'omitnan'), ...
        mean(X.frequency_snr_db,'omitnan'), mean(X.predicted_patch_purity,'omitnan'), ...
        mean(X.p_mixed,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','underestimate_pct','mean_abs_q_error', ...
        'mean_frequency_snr_db','mean_predicted_patch_purity','mean_p_mixed'});
end
S = [groups vertcat(rows{:})];
end

function D = degradation_table(S)
if isempty(S), D = table(); return; end
models = unique(S.model_name, 'stable');
rows = {};
for i = 1:numel(models)
    m = models(i);
    clean = S(S.model_name == m & S.realism_level == "clean", :);
    if isempty(clean), continue; end
    base = clean.MAPE_pct(1);
    levels = S(S.model_name == m, :);
    for j = 1:height(levels)
        rows{end+1,1} = table(m, levels.realism_level(j), levels.MAPE_pct(j), ...
            levels.MAPE_pct(j) - base, levels.high_error20_pct(j) - clean.high_error20_pct(1), ...
            'VariableNames', {'model_name','realism_level','MAPE_pct', ...
            'delta_MAPE_vs_clean_pctpt','delta_high20_vs_clean_pctpt'}); %#ok<AGROW>
    end
end
D = vertcat(rows{:});
end

function plot_summary(T_by_level, T_degrade, T_by_region, OUT)
fig = figure('Color','w','Position',[100 100 1200 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;
bar_labeled(compose("%s | %s", T_by_level.model_name, T_by_level.realism_level), ...
    T_by_level.MAPE_pct); grid on; ylabel('MAPE (%)');
title('Error by realism level'); xtickangle(35);
nexttile;
bar_labeled(compose("%s | %s", T_degrade.model_name, T_degrade.realism_level), ...
    T_degrade.delta_MAPE_vs_clean_pctpt); grid on; ylabel('\Delta MAPE vs clean (points)');
title('Clean-to-realistic degradation'); xtickangle(35);
nexttile;
keep = T_by_region.analysis_region ~= "other";
bar_labeled(compose("%s | %s", T_by_region.model_name(keep), T_by_region.analysis_region(keep)), ...
    T_by_region.MAPE_pct(keep)); grid on; ylabel('MAPE (%)');
title('Regional error'); xtickangle(35);
saveas(fig, fullfile(OUT.figure_dir, 'test54_realistic_readout_summary.png'));
close(fig);
end

function bar_labeled(labels, values)
labels = string(labels);
bar(1:numel(values), values);
set(gca, 'XTick', 1:numel(values), 'XTickLabel', cellstr(labels));
end

function plot_condition_maps(T, F, SIM, key, OUT, CFG)
primary = "q_spectrum_plus_composition";
if ~any(T.model_name == primary), primary = T.model_name(1); end
X = T(T.model_name == primary, :);
fig = figure('Color','w','Position',[80 80 1500 780]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
plot_img(SIM.cs_map, 'True SWS', 'SWS (m/s)');
plot_img(abs(SIM.U), 'Amplitude', '|velocity|');
plot_img(angle(SIM.U), 'Phase', 'rad');
plot_scatter_map(F, F.local_amplitude, 'Window amplitude', 'amplitude');
plot_scatter_map(X, X.sws_pred, 'Predicted SWS', 'SWS (m/s)');
plot_scatter_map(X, X.sws_abs_error_pct, 'Absolute error', 'error (%)');
plot_scatter_map(X, X.sws_signed_error_pct, 'Signed error', 'error (%)');
plot_scatter_map(X, X.q_pred, 'Predicted q', 'REQ quantile q');
plot_scatter_map(X, X.predicted_patch_purity, 'Predicted patch purity', 'probability');
plot_scatter_map(X, X.p_mixed, 'Predicted mixedness', 'probability');
plot_scatter_map(X, X.distance_to_boundary_mm, 'Distance to boundary', 'mm');
plot_scatter_map(X, X.frequency_snr_db, 'Frequency-domain SNR', 'dB');
sgtitle(strrep(string(key), '_', '\_'));
file = fullfile(OUT.map_dir, "test54_map__" + sanitize(key) + ".png");
saveas(fig, file);
close(fig);
end

function plot_img(A, ttl, cb)
nexttile; imagesc(A); axis image off; title(ttl); colorbar; ylabel(colorbar, cb);
end

function plot_scatter_map(T, v, ttl, cb)
nexttile;
scatter(T.x_center_m*1e3, T.z_center_m*1e3, 16, v, 'filled');
axis image ij; xlabel('x (mm)'); ylabel('z (mm)'); title(ttl);
c = colorbar; ylabel(c, cb);
end

function print_summary(T_overall, T_degrade, OUT)
fprintf('\nInterpretive summary:\n');
if isempty(T_overall), return; end
T = sortrows(T_overall, 'MAPE_pct');
fprintf('  Best global model: %s (MAPE %.2f%%, high>20 %.2f%%).\n', ...
    T.model_name(1), T.MAPE_pct(1), T.high_error20_pct(1));
if ~isempty(T_degrade)
    D = sortrows(T_degrade(T_degrade.realism_level ~= "clean",:), ...
        'delta_MAPE_vs_clean_pctpt', 'descend');
    if ~isempty(D)
        fprintf('  Largest realistic-readout degradation: %s / %s %+0.2f MAPE points.\n', ...
            D.model_name(1), D.realism_level(1), D.delta_MAPE_vs_clean_pctpt(1));
    end
end
fprintf('  Use this as robustness evidence, not as retraining data yet.\n');
fprintf('  Outputs: %s\n', OUT.root_dir);
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
switch string(regime)
    case "directional_2D", type = "SingleWave";
    case "diffuse_2D", type = "Diffuse2D";
    otherwise, type = "Diffuse3D";
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
v = round(double(patch(:))*1e5)/1e5;
v = v(isfinite(v));
if isempty(v), p = NaN; return; end
[~,~,ic] = unique(v);
p = max(accumarray(ic,1)) / numel(v);
end

function D = distance_to_boundary(cs_map, dx)
if max(cs_map(:)) - min(cs_map(:)) < 1e-8
    D = inf(size(cs_map)); return;
end
BW = cs_map > (min(cs_map(:)) + max(cs_map(:))) / 2;
B = bwperim(BW);
D = bwdist(B) * dx * 1e3;
end

function bin = purity_bin(p)
edges = [0 0.7 0.9 0.99 1.00001];
labels = ["strongly_mixed","moderately_mixed","near_pure","pure"];
idx = discretize(p, edges);
bin = labels(max(1,min(numel(labels),idx))).';
end

function bin = distance_bin(d, edges)
labels = strings(1,numel(edges)-1);
for i = 1:numel(labels)
    if isinf(edges(i+1)), labels(i) = sprintf('>%gmm',edges(i));
    else, labels(i) = sprintf('%g_%gmm',edges(i),edges(i+1));
    end
end
idx = discretize(d, edges);
idx(isnan(idx)) = numel(labels);
bin = labels(max(1,min(numel(labels),idx))).';
end

function r = analysis_region(F)
r = repmat("other", height(F), 1);
r(F.case_family == "homogeneous") = "homogeneous";
r(F.distance_to_boundary_mm <= 2) = "interface_0_2mm";
r(F.patch_purity >= 0.99 & F.distance_to_boundary_mm > 4) = "pure_core";
r(F.patch_purity < 0.9) = "mixed_patch";
end

function v = getfield_or(S, name, default_value)
if isfield(S, name), v = S.(name); else, v = default_value; end
end

function snr = get_snr_db(SIM)
snr = Inf;
try
    if isfield(SIM.diag, 'frequencySNRdB')
        snr = SIM.diag.frequencySNRdB;
    elseif isfield(SIM.diag, 'timeReadoutNoise') && isfield(SIM.diag.timeReadoutNoise, 'SNRdB_time')
        snr = SIM.diag.timeReadoutNoise.SNRdB_time;
    end
catch
    snr = Inf;
end
end

function y = clamp01(y)
y = min(max(y, 0), 1);
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^\w\d-]+', '_');
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
end

function vals = env_string_list(name, default)
v = strtrim(string(getenv(name)));
if v == "", vals = string(default); return; end
vals = string(strtrim(split(v, {',',';',' '})));
vals = vals(vals ~= "");
end

function write_config_json(CFG, file)
fid = fopen(file, 'w');
assert(fid > 0, 'Could not write config JSON: %s', file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(CFG, 'PrettyPrint', true));
end
