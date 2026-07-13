%% Test 71: Eikonal cs_guess sweep diagnostic
% Diagnose whether hard-material underestimation improves when REQ is
% re-extracted with different operational cs_guess values.
%
% This script does not retrain q, composition, confidence, or correction
% models. It only loads frozen Test38/Test53 models and applies them to
% cached Eikonal fields from Test66/Test70.

clear; clc; close all;
format compact;

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(ROOT);
if exist('setup_adaptive_req','file') == 2
    setup_adaptive_req;
    if exist('PROJECT_ROOT','var')
        ROOT = PROJECT_ROOT;
    end
end

CFG = local_config(ROOT);
OUT = make_outputs(ROOT, CFG);

fprintf('\nTest 71: Eikonal cs_guess sweep diagnostic\n');
fprintf('Mode: %s | target step %.3f mm | cs_guess=%s | M=%s\n', ...
    CFG.Mode, 1e3*CFG.TargetStepM, mat2str(CFG.CsGuesses), mat2str(CFG.MValues));
fprintf('Source field cache: %s\n', CFG.FieldCache);
fprintf('No training. Frozen q/composition models only.\n\n');

MODELS = load_models(ROOT, CFG);
BASE_FEATURES = string(MODELS.base_features(:));
assert_no_leakage_predictors(BASE_FEATURES);
fprintf('Loaded frozen model bundle: %s\n', MODELS.bundle_path);
fprintf('Frozen predictors required: %d\n', numel(BASE_FEATURES));

FIELD_FILES = select_field_files(CFG);
fprintf('Selected cached fields: %d\n', numel(FIELD_FILES));
if isempty(FIELD_FILES)
    error('No cached Eikonal fields matched the Test71 filters.');
end

t_all = tic;
ALL = table();
for iField = 1:numel(FIELD_FILES)
    S = load(FIELD_FILES(iField), 'SIM', 'C', 'cfg');
    SIM = S.SIM;
    C = S.C;
    cond_key = condition_key(C, SIM, FIELD_FILES(iField));
    fprintf('\n[%d/%d] %s\n', iField, numel(FIELD_FILES), cond_key);

    for M = CFG.MValues
        for cs_guess = CFG.CsGuesses
            t_cond = tic;
            try
                F = extract_or_load_features(SIM, C, cond_key, M, cs_guess, CFG, OUT);
                F = ensure_predictors(F, BASE_FEATURES, CFG);
                R = apply_models(F, MODELS, CFG);
                R.condition_key = repmat(string(cond_key), height(R), 1);
                R.elapsed_condition_s = repmat(toc(t_cond), height(R), 1);
                ALL = [ALL; R]; %#ok<AGROW>
                fprintf('  M=%.1f cs_guess=%.2f: %d patches x %d models in %.1f s\n', ...
                    M, cs_guess, height(F), numel(unique(R.model_name)), toc(t_cond));
            catch ME
                warning('Test71:ConditionFailed', ...
                    'Failed %s | M=%.1f | cs_guess=%.2f: %s', ...
                    cond_key, M, cs_guess, ME.message);
            end
        end
    end
end

if isempty(ALL)
    error('Test71 produced no rows.');
end

fprintf('\nWriting Test71 outputs...\n');
writetable(ALL, fullfile(OUT.tables, 'test71_patch_level_results.csv'));

SUM = struct();
SUM.overall = summarize_metrics(ALL, ["model_name","cs_guess","REQ_M"]);
SUM.geometry = summarize_metrics(ALL, ["model_name","cs_guess","geometry","REQ_M"]);
SUM.frequency = summarize_metrics(ALL, ["model_name","cs_guess","f0","REQ_M"]);
SUM.geometry_frequency = summarize_metrics(ALL, ["model_name","cs_guess","geometry","f0","REQ_M"]);
SUM.field_regime = summarize_metrics(ALL, ["model_name","cs_guess","field_regime","REQ_M"]);
SUM.realism = summarize_metrics(ALL, ["model_name","cs_guess","realism_level","REQ_M"]);
SUM.roi = summarize_metrics(ALL, ["model_name","cs_guess","roi_region","REQ_M"]);
SUM.roi_frequency = summarize_metrics(ALL, ["model_name","cs_guess","roi_region","f0","REQ_M"]);
SUM.purity = summarize_metrics(ALL, ["model_name","cs_guess","purity_bin","REQ_M"]);
SUM.distance = summarize_metrics(ALL, ["model_name","cs_guess","distance_over_window_bin","REQ_M"]);

writetable(SUM.overall, fullfile(OUT.tables, 'test71_summary_overall.csv'));
writetable(SUM.geometry, fullfile(OUT.tables, 'test71_summary_by_geometry.csv'));
writetable(SUM.frequency, fullfile(OUT.tables, 'test71_summary_by_frequency.csv'));
writetable(SUM.geometry_frequency, fullfile(OUT.tables, 'test71_summary_by_geometry_frequency.csv'));
writetable(SUM.field_regime, fullfile(OUT.tables, 'test71_summary_by_field_regime.csv'));
writetable(SUM.realism, fullfile(OUT.tables, 'test71_summary_by_realism.csv'));
writetable(SUM.roi, fullfile(OUT.tables, 'test71_summary_by_roi.csv'));
writetable(SUM.roi_frequency, fullfile(OUT.tables, 'test71_summary_by_roi_frequency.csv'));
writetable(SUM.purity, fullfile(OUT.tables, 'test71_summary_by_patch_purity.csv'));
writetable(SUM.distance, fullfile(OUT.tables, 'test71_summary_by_distance_over_window.csv'));

save(fullfile(OUT.data, 'test71_compact_results.mat'), 'CFG', 'SUM', '-v7.3');

plot_summary_figures(SUM, OUT);
if CFG.SaveMaps
    plot_representative_maps(ALL, OUT, CFG);
end
write_readme(ALL, SUM, OUT, CFG, toc(t_all));

fprintf('\nTest 71 complete in %.1f min.\n', toc(t_all)/60);
fprintf('Tables : %s\n', OUT.tables);
fprintf('Figures: %s\n', OUT.figures);

%% Configuration and I/O

function CFG = local_config(ROOT)
CFG = struct();
CFG.Mode = lower(string(getenv_default('ADAPTIVE_REQ_TEST71_MODE', 'validate')));
CFG.ModelSource = string(getenv_default('ADAPTIVE_REQ_TEST71_MODEL_SOURCE', 'auto'));
CFG.SourceRun = string(getenv_default('ADAPTIVE_REQ_TEST71_SOURCE_RUN', ...
    fullfile(ROOT, 'outputs', 'test_66_eikonal_realistic_transfer_validation', ...
    'full_a', 'run_test70_m2m3_limited')));
CFG.FieldCache = fullfile(CFG.SourceRun, 'data', 'field_cache');
CFG.CsGuesses = env_num_list('ADAPTIVE_REQ_TEST71_CS_GUESSES', [2 3 4]);
CFG.MValues = env_num_list('ADAPTIVE_REQ_TEST71_M_VALUES', 2);
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST71_TARGET_STEP_M', 0.5e-3);
CFG.SaveMaps = env_true('ADAPTIVE_REQ_TEST71_SAVE_ALL_MAPS', true);
CFG.UseParfor = env_true('ADAPTIVE_REQ_TEST71_USE_PARFOR', false);
CFG.IncludeDiagnostics = env_true('ADAPTIVE_REQ_TEST71_INCLUDE_ORACLE_THEORY', true);
CFG.MaxExampleMaps = env_number('ADAPTIVE_REQ_TEST71_MAX_EXAMPLE_MAPS', 24);

CFG.REQ = struct();
CFG.REQ.Nbins = 'auto';
CFG.REQ.NbinsMin = 16;
CFG.REQ.NbinsAutoOversample = 1;
CFG.REQ.SmoothSigma = 1;
CFG.REQ.Gamma = 1;
CFG.REQ.PadFactor = 1;
CFG.REQ.EdgeMode = 'valid';

switch CFG.Mode
    case "validate"
        CFG.Geometries = ["homogeneous_cs4","inclusion_2_4"];
        CFG.Frequencies = 300;
        CFG.RealismLevels = "clean";
        CFG.FieldRegimes = "single_source_lateral";
        CFG.SourceSeeds = 1;
        CFG.NoiseSeeds = 0;
        CFG.CsGuesses = CFG.CsGuesses(1:min(2,numel(CFG.CsGuesses)));
        CFG.MaxFields = 2;
    case "quick"
        CFG.Geometries = ["homogeneous_cs4","bilayer_2_4","inclusion_2_4","bilayer_inclusion_2_3_4"];
        CFG.Frequencies = [300 500];
        CFG.RealismLevels = ["clean","readout_medium"];
        CFG.FieldRegimes = ["single_source_lateral","diffuse_like_8src"];
        CFG.SourceSeeds = 1;
        CFG.NoiseSeeds = [0 1];
        CFG.MaxFields = inf;
    otherwise
        CFG.Geometries = string.empty;
        CFG.Frequencies = [];
        CFG.RealismLevels = string.empty;
        CFG.FieldRegimes = string.empty;
        CFG.SourceSeeds = [];
        CFG.NoiseSeeds = [];
        CFG.MaxFields = inf;
end
end

function OUT = make_outputs(ROOT, CFG)
OUT.root = fullfile(ROOT, 'outputs', 'test_71_eikonal_cs_guess_sweep_diagnostic', char(CFG.Mode));
OUT.tables = fullfile(OUT.root, 'tables');
OUT.figures = fullfile(OUT.root, 'figures');
OUT.data = fullfile(OUT.root, 'data');
OUT.req_cache = fullfile(OUT.data, 'req_cache');
mkdir(OUT.tables); mkdir(OUT.figures); mkdir(OUT.data); mkdir(OUT.req_cache);
end

function MODELS = load_models(ROOT, CFG)
paths = strings(0,1);
if CFG.ModelSource ~= "auto"
    if isfile(CFG.ModelSource)
        paths(end+1) = CFG.ModelSource;
    else
        paths(end+1) = model_source_path(ROOT, CFG.ModelSource);
    end
end
paths = [paths; ...
    fullfile(ROOT, 'outputs', 'test_53_paper_final_clean_q_training', 'models', 'test38_velocity_field_diverse_q_models.mat'); ...
    fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'models', 'test38_velocity_field_diverse_q_models.mat'); ...
    fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'medium', 'models', 'test38_velocity_field_diverse_q_models.mat'); ...
    fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'quick', 'models', 'test38_velocity_field_diverse_q_models.mat')];
paths = unique(paths, 'stable');

last_error = "";
for p = paths(:)'
    if ~isfile(p), continue; end
    try
        S = load(p);
        MODELS = normalize_bundle(S);
        MODELS.bundle_path = string(p);
        return;
    catch ME
        last_error = ME.message;
        warning('Test71:BundleLoadFailed', 'Could not load %s: %s', p, ME.message);
    end
end
error('No usable frozen model bundle found. Last error: %s', last_error);
end

function p = model_source_path(ROOT, name)
switch lower(string(name))
    case "test53"
        p = fullfile(ROOT, 'outputs', 'test_53_paper_final_clean_q_training', 'models', 'test38_velocity_field_diverse_q_models.mat');
    case "full"
        p = fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'models', 'test38_velocity_field_diverse_q_models.mat');
    case "medium"
        p = fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'medium', 'models', 'test38_velocity_field_diverse_q_models.mat');
    otherwise
        p = fullfile(ROOT, 'outputs', 'test_38_velocity_field_diverse_q_training', 'quick', 'models', 'test38_velocity_field_diverse_q_models.mat');
end
end

function MODELS = normalize_bundle(S)
if isfield(S, 'MODEL_BUNDLE')
    B = S.MODEL_BUNDLE;
elseif isfield(S, 'MODELS')
    B = S.MODELS;
else
    B = S;
end

MODELS = struct();
MODELS.q_spectrum_only = find_q_model(B, ["q_spectrum_only","spectrum_only"]);
MODELS.q_spectrum_plus_composition = find_q_model(B, ["q_spectrum_plus_composition","spectrum_plus_composition"]);
MODELS.composition = find_composition(B);

features = string(MODELS.q_spectrum_only.features(:));
features = union(features, string(MODELS.q_spectrum_plus_composition.features(:)), 'stable');
if ~isempty(MODELS.composition.features)
    features = union(features, string(MODELS.composition.features(:)), 'stable');
end
MODELS.base_features = features;
end

function M = find_q_model(B, names)
for name = names
    candidates = {B};
    if isfield(B, 'q_models'), candidates{end+1} = B.q_models; end
    if isfield(B, 'models'), candidates{end+1} = B.models; end
    for ii = 1:numel(candidates)
        C = candidates{ii};
        if isstruct(C) && isfield(C, name)
            M = C.(name);
            return;
        end
    end
end
error('Missing frozen q model: %s', strjoin(names, ', '));
end

function C = find_composition(B)
C = struct('purity', [], 'mixed', [], 'strong_mixed', [], 'features', string.empty);
roots = {B};
if isfield(B, 'composition'), roots{end+1} = B.composition; end
if isfield(B, 'composition_models'), roots{end+1} = B.composition_models; end
if isfield(B, 'models') && isfield(B.models, 'composition'), roots{end+1} = B.models.composition; end
for ii = 1:numel(roots)
    R = roots{ii};
    C.purity = first_field(R, ["purity","patch_purity","predicted_patch_purity","purity_model"]);
    C.mixed = first_field(R, ["mixed","p_mixed","is_mixed","mixed_model"]);
    C.strong_mixed = first_field(R, ["strong_mixed","p_strong_mixed","is_strong_mixed","strong_mixed_model"]);
    if ~isempty(C.purity) && ~isempty(C.mixed) && ~isempty(C.strong_mixed)
        C.features = get_model_features(C.purity);
        return;
    end
end
error('Missing frozen composition/purity models.');
end

function v = first_field(S, names)
v = [];
if ~isstruct(S), return; end
for name = names
    if isfield(S, name)
        v = S.(name);
        return;
    end
end
end

function f = get_model_features(M)
if isstruct(M) && isfield(M, 'features')
    f = string(M.features(:));
elseif isprop(M, 'PredictorNames')
    f = string(M.PredictorNames(:));
else
    f = string.empty;
end
end

%% Field selection and REQ extraction

function files = select_field_files(CFG)
raw = dir(fullfile(CFG.FieldCache, 'field__*.mat'));
meta = table();
for ii = 1:numel(raw)
    p = fullfile(raw(ii).folder, raw(ii).name);
    try
        S = load(p, 'C', 'SIM');
        r = table();
        r.file = string(p);
        r.geometry = string(get_meta(S.C, S.SIM, 'geometry', 'unknown'));
        r.f0 = double(get_meta(S.C, S.SIM, 'f0', NaN));
        r.realism_level = string(get_meta(S.C, S.SIM, 'realism_level', 'unknown'));
        r.field_regime = string(get_meta(S.C, S.SIM, 'field_regime', 'unknown'));
        r.source_seed = double(get_meta(S.C, S.SIM, 'source_seed', NaN));
        r.noise_seed = double(get_meta(S.C, S.SIM, 'noise_seed', 0));
        meta = [meta; r]; %#ok<AGROW>
    catch
    end
end
if isempty(meta)
    files = strings(0,1);
    return;
end
keep = true(height(meta),1);
if ~isempty(CFG.Geometries), keep = keep & ismember(meta.geometry, CFG.Geometries); end
if ~isempty(CFG.Frequencies), keep = keep & ismember(meta.f0, CFG.Frequencies); end
if ~isempty(CFG.RealismLevels), keep = keep & ismember(meta.realism_level, CFG.RealismLevels); end
if ~isempty(CFG.FieldRegimes), keep = keep & ismember(meta.field_regime, CFG.FieldRegimes); end
if ~isempty(CFG.SourceSeeds), keep = keep & ismember(meta.source_seed, CFG.SourceSeeds); end
if ~isempty(CFG.NoiseSeeds), keep = keep & ismember(meta.noise_seed, CFG.NoiseSeeds); end
meta = sortrows(meta(keep,:), {'geometry','f0','realism_level','field_regime','source_seed','noise_seed'});
if isfinite(CFG.MaxFields)
    meta = meta(1:min(height(meta), CFG.MaxFields), :);
end
files = meta.file;
end

function F = extract_or_load_features(SIM, C, cond_key, M, cs_guess, CFG, OUT)
step_x = max(1, round(CFG.TargetStepM / SIM.dx));
step_z = max(1, round(CFG.TargetStepM / SIM.dz));
cache_file = fullfile(OUT.req_cache, sprintf('req__%s__M%s__cs%s__step%dpx.mat', ...
    sanitize(cond_key), num2str(M), num2str(cs_guess), step_x));
if isfile(cache_file)
    S = load(cache_file, 'F');
    F = S.F;
    return;
end

feat = adaptive_req.config.default_feature_config('M', M, 'cs_guess', cs_guess, ...
    'gamma_win', CFG.REQ.Gamma, 'pad_factor', CFG.REQ.PadFactor);
cfg_req = struct('dx', SIM.dx, 'dz', SIM.dz, 'f0', SIM.f0, ...
    'cs_bg', cs_guess, 'WaveModel', 'simcore_eikonal_test71');

[~, F] = adaptive_req.estimators.req_estimator_map(SIM.U, cfg_req, feat, ...
    'StepX', step_x, 'StepZ', step_z, 'EdgeMode', CFG.REQ.EdgeMode, ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins', CFG.REQ.Nbins, ...
        'Nbins_auto_oversample', CFG.REQ.NbinsAutoOversample, ...
        'Nbins_min', CFG.REQ.NbinsMin, ...
        'smooth_sigma', CFG.REQ.SmoothSigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseParfor, ...
    'StoreReqCurves', true, 'Verbose', false);

F = attach_eval_diagnostics(F, SIM, C, M, cs_guess, step_x, step_z, CFG);
save(cache_file, 'F', '-v7.3');
end

function F = attach_eval_diagnostics(F, SIM, C, M, cs_guess, step_x, step_z, CFG)
n = height(F);
ix = max(1, min(size(SIM.cs_map,2), F.map_ix));
iz = max(1, min(size(SIM.cs_map,1), F.map_iz));
idx = sub2ind(size(SIM.cs_map), iz, ix);

F.geometry = repmat(string(SIM.geometry), n, 1);
F.case_id = F.geometry;
F.geometry_family = geometry_family(F.geometry);
F.case_family = F.geometry_family;
F.realism_level = repmat(string(SIM.realism_level), n, 1);
F.field_regime = repmat(string(get_meta(C, SIM, 'field_regime', 'unknown')), n, 1);
F.field_regime_ood = F.field_regime;
F.source_seed = repmat(double(get_meta(C, SIM, 'source_seed', NaN)), n, 1);
F.noise_seed = repmat(double(get_meta(C, SIM, 'noise_seed', 0)), n, 1);
F.f0 = repmat(double(SIM.f0), n, 1);
F.SIM_f0 = F.f0;
F.M = repmat(M, n, 1);
F.REQ_M = F.M;
F.cs_guess = repmat(cs_guess, n, 1);
F.cs_guess_used = F.cs_guess;
F.dx = repmat(double(SIM.dx), n, 1);
F.dz = repmat(double(SIM.dz), n, 1);
F.TargetStepM = repmat(CFG.TargetStepM, n, 1);
F.REQ_StepX = repmat(step_x, n, 1);
F.REQ_StepZ = repmat(step_z, n, 1);
F.REQ_Nbins_effective = repmat(NaN, n, 1);

F.x_m = SIM.x(ix(:));
F.z_m = SIM.z(iz(:));
F.depth_m = F.z_m;
F.true_SWS = SIM.cs_map(idx);
F.SWS_true = F.true_SWS;
F.k_true = 2*pi*F.f0 ./ F.true_SWS;

F.patch_purity = patch_purity(SIM.label_map, iz, ix, step_x, step_z);
F.patch_cs_std = patch_cs_stat(SIM.cs_map, iz, ix, step_x, step_z, "std");
F.patch_cs_range = patch_cs_stat(SIM.cs_map, iz, ix, step_x, step_z, "range");
F.patch_cs_iqr = patch_cs_stat(SIM.cs_map, iz, ix, step_x, step_z, "iqr");
F.is_mixed = F.patch_purity < 0.95;
F.is_strong_mixed = F.patch_purity < 0.75;
F.purity_bin = purity_bin(F.patch_purity);

if isfield(SIM, 'abs_distance_m') && ~isempty(SIM.abs_distance_m)
    F.distance_to_interface_m = SIM.abs_distance_m(idx);
else
    F.distance_to_interface_m = nan(n,1);
end
F.distance_to_interface_mm = 1e3 * F.distance_to_interface_m;
F.window_radius_m = repmat(0.5*max(step_x*SIM.dx, step_z*SIM.dz), n, 1);
F.distance_to_interface_over_window_radius = F.distance_to_interface_m ./ max(F.window_radius_m, eps);
F.distance_over_window_bin = distance_window_bin(F.distance_to_interface_over_window_radius);

F.local_amplitude = sample_map(SIM, 'amp_map', idx, abs(SIM.U(idx)));
F.local_phase_rad = sample_map(SIM, 'phase_map', idx, angle(SIM.U(idx)));
F.tracking_snr_db = sample_map(SIM, 'tracking_snr_db', idx, nan(n,1));
F.snr_proxy_db = F.tracking_snr_db;
F.shear_attenuation_factor = sample_map(SIM, 'shear_attenuation_factor', idx, nan(n,1));
F.acoustic_readout_amplitude_factor = sample_map(SIM, 'acoustic_readout_amplitude_factor', idx, nan(n,1));
F.roi_region = roi_labels(SIM, idx, n);

F.q_oracle = nan(n,1);
if ismember('req_mapping', string(F.Properties.VariableNames))
    for ii = 1:n
        F.q_oracle(ii) = invert_mapping_to_q(F.req_mapping{ii}, F.k_true(ii));
    end
end
qt = theory_q(SIM, M, cs_guess, F.field_regime(1), CFG);
F.q_theory_prior = repmat(qt, n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
end

%% Prediction and metrics

function F = ensure_predictors(F, features, CFG)
for f = string(features(:))'
    if ismember(f, string(F.Properties.VariableNames)), continue; end
    switch f
        case "REQ_StepX"
            F.REQ_StepX = max(1, round(CFG.TargetStepM ./ F.dx));
        case "REQ_StepZ"
            F.REQ_StepZ = max(1, round(CFG.TargetStepM ./ F.dz));
        case "TargetStepM"
            F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
        case "REQ_Nbins_effective"
            F.REQ_Nbins_effective = nan(height(F),1);
        otherwise
            error('Required frozen-model predictor missing: %s', f);
    end
end
end

function R = apply_models(F, MODELS, CFG)
COMP = predict_composition(MODELS.composition, F);
F.predicted_patch_purity = COMP.predicted_patch_purity;
F.p_mixed = COMP.p_mixed;
F.p_strong_mixed = COMP.p_strong_mixed;

R = table();
R = [R; rows_for_q(F, "q_spectrum_only", predict_q(MODELS.q_spectrum_only, F))]; %#ok<AGROW>
R = [R; rows_for_q(F, "q_spectrum_plus_composition", predict_q(MODELS.q_spectrum_plus_composition, F))]; %#ok<AGROW>
if CFG.IncludeDiagnostics
    R = [R; rows_for_q(F, "oracle_q", F.q_oracle)]; %#ok<AGROW>
    R = [R; rows_for_q(F, "theory_q_discrete", F.q_theory_prior)]; %#ok<AGROW>
end
end

function q = predict_q(M, F)
X = F(:, cellstr(string(M.features(:))));
q = clamp01(predict(M.model, X));
end

function COMP = predict_composition(C, F)
X = F(:, cellstr(string(C.features(:))));
COMP.predicted_patch_purity = clamp01(predict(C.purity, X));
[~, score] = predict(C.mixed, X);
COMP.p_mixed = positive_class_score(C.mixed, score);
[~, score] = predict(C.strong_mixed, X);
COMP.p_strong_mixed = positive_class_score(C.strong_mixed, score);
end

function T = rows_for_q(F, model_name, q)
T = output_base(F);
sws = q_to_sws(F.req_mapping, q, F.f0);
T.model_name = repmat(string(model_name), height(F), 1);
T.q_pred = q;
T.SWS_pred = sws;
T.k_pred = 2*pi*T.f0 ./ T.SWS_pred;
T.sws_signed_error_pct = 100*(T.SWS_pred - T.true_SWS) ./ T.true_SWS;
T.sws_abs_error_pct = abs(T.sws_signed_error_pct);
T.high_error_gt10 = T.sws_abs_error_pct > 10;
T.high_error_gt20 = T.sws_abs_error_pct > 20;
T.q_error = T.q_pred - T.q_oracle;
T.q_abs_error = abs(T.q_error);
end

function T = output_base(F)
vars = ["geometry","geometry_family","case_id","case_family","realism_level", ...
    "field_regime","field_regime_ood","source_seed","noise_seed","f0","SIM_f0", ...
    "M","REQ_M","cs_guess","cs_guess_used","dx","dz","TargetStepM", ...
    "REQ_StepX","REQ_StepZ","map_ix","map_iz","x_m","z_m","true_SWS", ...
    "SWS_true","k_true","q_oracle","q_theory_prior","sws_theory", ...
    "patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr", ...
    "is_mixed","is_strong_mixed","purity_bin","distance_to_interface_m", ...
    "distance_to_interface_mm","window_radius_m", ...
    "distance_to_interface_over_window_radius","distance_over_window_bin", ...
    "depth_m","local_amplitude","local_phase_rad","tracking_snr_db", ...
    "snr_proxy_db","shear_attenuation_factor","acoustic_readout_amplitude_factor", ...
    "roi_region","predicted_patch_purity","p_mixed","p_strong_mixed"];
vars = vars(ismember(vars, string(F.Properties.VariableNames)));
T = F(:, cellstr(vars));
end

function S = summarize_metrics(T, groups)
valid = isfinite(T.SWS_pred) & isfinite(T.true_SWS) & T.true_SWS > 0;
T = T(valid,:);
if isempty(T), S = table(); return; end
[G, keys] = findgroups(T(:, cellstr(groups)));
S = keys;
S.N = splitapply(@numel, T.sws_abs_error_pct, G);
S.MAPE_pct = splitapply(@omitmean, T.sws_abs_error_pct, G);
S.MedianAPE_pct = splitapply(@omitmedian, T.sws_abs_error_pct, G);
S.Bias_pct = splitapply(@omitmean, T.sws_signed_error_pct, G);
S.MedianSignedError_pct = splitapply(@omitmedian, T.sws_signed_error_pct, G);
S.MAE_sws = splitapply(@(p,t) omitmean(abs(p-t)), T.SWS_pred, T.true_SWS, G);
S.RMSE_sws = splitapply(@(p,t) sqrt(omitmean((p-t).^2)), T.SWS_pred, T.true_SWS, G);
S.High10_pct = 100*splitapply(@omitmean, double(T.high_error_gt10), G);
S.High20_pct = 100*splitapply(@omitmean, double(T.high_error_gt20), G);
S.Under_pct = 100*splitapply(@omitmean, double(T.sws_signed_error_pct < 0), G);
S.MeanQ = splitapply(@omitmean, T.q_pred, G);
S.MeanQOracle = splitapply(@omitmean, T.q_oracle, G);
S.MeanQError = splitapply(@omitmean, T.q_error, G);
S.MeanPurity = splitapply(@omitmean, T.patch_purity, G);
S.MeanSNR_dB = splitapply(@omitmean, T.snr_proxy_db, G);
end

%% Plotting and report

function plot_summary_figures(SUM, OUT)
try
    plot_global_bar(SUM.overall, OUT);
    plot_frequency(SUM.frequency, OUT);
    plot_roi_bias(SUM.roi, OUT);
    plot_purity(SUM.purity, OUT);
catch ME
    warning('Test71:PlotFailed', 'Summary plotting failed: %s', ME.message);
end
end

function plot_global_bar(S, OUT)
if isempty(S), return; end
S = S(ismember(S.model_name, ["q_spectrum_only","q_spectrum_plus_composition","oracle_q"]), :);
figure('Color','w','Position',[100 100 1050 420]);
bar(categorical(strcat(string(S.model_name), " | cs=", string(S.cs_guess))), S.MAPE_pct);
grid on; ylabel('MAPE (%)'); title('Test 71 global MAPE by cs guess'); xtickangle(35);
saveas(gcf, fullfile(OUT.figures, 'test71_global_mape_by_cs_guess.png')); close(gcf);
end

function plot_frequency(S, OUT)
if isempty(S), return; end
figure('Color','w','Position',[100 100 1050 520]); hold on;
for m = ["q_spectrum_only","q_spectrum_plus_composition","oracle_q"]
    for cs = unique(S.cs_guess)'
        rows = sortrows(S(S.model_name == m & S.cs_guess == cs, :), 'f0');
        if isempty(rows), continue; end
        plot(rows.f0, rows.MAPE_pct, '-o', 'DisplayName', sprintf('%s cs=%.1f', m, cs));
    end
end
grid on; xlabel('Frequency (Hz)'); ylabel('MAPE (%)');
title('MAPE versus frequency by cs guess'); legend('Location','eastoutside', 'Interpreter','none');
saveas(gcf, fullfile(OUT.figures, 'test71_mape_vs_frequency_by_cs_guess.png')); close(gcf);
end

function plot_roi_bias(S, OUT)
if isempty(S), return; end
S = S(S.model_name == "q_spectrum_plus_composition", :);
figure('Color','w','Position',[100 100 1150 450]);
bar(categorical(strcat(string(S.roi_region), " | cs=", string(S.cs_guess))), S.Bias_pct);
grid on; ylabel('Mean signed SWS error (%)'); title('Composition model bias by ROI and cs guess'); xtickangle(35);
saveas(gcf, fullfile(OUT.figures, 'test71_roi_bias_composition_by_cs_guess.png')); close(gcf);
end

function plot_purity(S, OUT)
if isempty(S), return; end
S = S(S.model_name == "q_spectrum_plus_composition", :);
figure('Color','w','Position',[100 100 1050 520]); hold on;
bins = unique(string(S.purity_bin), 'stable');
for cs = unique(S.cs_guess)'
    rows = S(S.cs_guess == cs, :);
    [~, x] = ismember(string(rows.purity_bin), bins);
    plot(x, rows.MAPE_pct, '-o', 'DisplayName', sprintf('cs=%.1f', cs));
end
grid on; xticks(1:numel(bins)); xticklabels(bins); xtickangle(30);
xlabel('Patch purity bin'); ylabel('MAPE (%)'); title('Composition error versus patch purity');
legend('Location','best');
saveas(gcf, fullfile(OUT.figures, 'test71_mape_vs_patch_purity_by_cs_guess.png')); close(gcf);
end

function plot_representative_maps(T, OUT, CFG)
map_dir = fullfile(OUT.figures, 'maps_by_condition');
mkdir(map_dir);
T = T(T.model_name == "q_spectrum_plus_composition", :);
keys = unique(T.condition_key, 'stable');
keys = keys(1:min(numel(keys), CFG.MaxExampleMaps));
for k = 1:numel(keys)
    R = T(T.condition_key == keys(k), :);
    cs_values = unique(R.cs_guess)';
    ncols = 2 + numel(cs_values)*3;
    figure('Color','w','Position',[50 50 290*ncols 780]);
    tiledlayout(3, ncols, 'TileSpacing','compact', 'Padding','compact');
    nexttile; map_panel(R, R.true_SWS, 'True SWS', 'm/s');
    nexttile; map_panel(R, R.patch_purity, 'Patch purity', '');
    for cs = cs_values
        r = R(R.cs_guess == cs, :);
        nexttile; map_panel(r, r.SWS_pred, sprintf('Pred SWS cs=%.1f', cs), 'm/s');
    end
    nexttile; map_panel(R, R.distance_to_interface_over_window_radius, 'Distance/window radius', '');
    nexttile; map_panel(R, R.snr_proxy_db, 'SNR proxy', 'dB');
    for cs = cs_values
        r = R(R.cs_guess == cs, :);
        nexttile; map_panel(r, r.sws_signed_error_pct, sprintf('Signed error cs=%.1f', cs), '%');
    end
    nexttile; map_panel(R, R.q_oracle, 'q oracle', 'q');
    nexttile; map_panel(R, R.local_amplitude, 'Amplitude', 'a.u.');
    for cs = cs_values
        r = R(R.cs_guess == cs, :);
        nexttile; map_panel(r, r.q_pred, sprintf('q pred cs=%.1f', cs), 'q');
    end
    sgtitle(strrep(keys(k), '_', '\_'));
    saveas(gcf, fullfile(map_dir, sprintf('test71_maps__%s.png', sanitize(keys(k)))));
    close(gcf);
end
end

function map_panel(T, values, ttl, units)
nz = max(T.map_iz); nx = max(T.map_ix);
G = nan(nz, nx);
G(sub2ind(size(G), T.map_iz, T.map_ix)) = values;
imagesc(G); axis image off; title(ttl, 'Interpreter','none');
cb = colorbar; if strlength(string(units)) > 0, ylabel(cb, units); end
end

function write_readme(T, SUM, OUT, CFG, runtime_s)
fid = fopen(fullfile(OUT.root, 'README_results.md'), 'w');
fprintf(fid, '# Test 71: Eikonal cs_guess sweep diagnostic\n\n');
fprintf(fid, '- Mode: `%s`\n', CFG.Mode);
fprintf(fid, '- Runtime: %.1f min\n', runtime_s/60);
fprintf(fid, '- Rows: %d\n', height(T));
fprintf(fid, '- cs_guess values: `%s`\n', mat2str(CFG.CsGuesses));
fprintf(fid, '- M values: `%s`\n', mat2str(CFG.MValues));
fprintf(fid, '- Field cache: `%s`\n\n', CFG.FieldCache);
fprintf(fid, '## Overall metrics\n\n');
S = SUM.overall;
if ~isempty(S)
    S = sortrows(S, {'model_name','cs_guess'});
    fprintf(fid, '| Model | cs_guess | M | N | MAPE %% | Bias %% | High20 %% |\n');
    fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
    for ii = 1:height(S)
        fprintf(fid, '| %s | %.2f | %.1f | %d | %.2f | %.2f | %.2f |\n', ...
            S.model_name(ii), S.cs_guess(ii), S.REQ_M(ii), S.N(ii), ...
            S.MAPE_pct(ii), S.Bias_pct(ii), S.High20_pct(ii));
    end
end
fprintf(fid, '\n## How to read this test\n\n');
fprintf(fid, '- If `oracle_q` improves with a different `cs_guess`, the REQ mapping/window choice matters.\n');
fprintf(fid, '- If `oracle_q` improves but the ML q models do not, the frozen model does not transfer to that altered mapping.\n');
fprintf(fid, '- If neither improves, the hard/interfacial error is likely field physics or REQ spectral ambiguity rather than only `cs_guess`.\n');
fclose(fid);
end

%% Small utilities

function assert_no_leakage_predictors(features)
forbidden = ["true_SWS","SWS_true","patch_purity","q_oracle","k_true", ...
    "sws_abs_error_pct","sws_signed_error_pct","distance_to_interface_m", ...
    "distance_to_interface_mm","confidence","risk","high_error_gt20"];
bad = intersect(string(features(:)), forbidden);
if ~isempty(bad)
    error('Leakage guard failed. Forbidden predictors: %s', strjoin(bad, ', '));
end
end

function key = condition_key(C, SIM, file)
if isfield(C, 'condition_key')
    key = string(C.condition_key);
elseif isfield(SIM, 'condition_key')
    key = string(SIM.condition_key);
else
    [~, name] = fileparts(file);
    key = erase(string(name), "field__");
end
end

function v = get_meta(C, SIM, name, fallback)
if isstruct(C) && isfield(C, name)
    v = C.(name);
elseif isstruct(SIM) && isfield(SIM, name)
    v = SIM.(name);
else
    v = fallback;
end
end

function family = geometry_family(g)
family = repmat("other", size(g));
family(contains(g, "homogeneous")) = "homogeneous";
family(contains(g, "bilayer")) = "bilayer";
family(contains(g, "inclusion")) = "inclusion";
family(contains(g, "three")) = "complex";
end

function p = patch_purity(label_map, iz, ix, step_x, step_z)
n = numel(iz); p = nan(n,1);
rx = max(1, round(step_x/2)); rz = max(1, round(step_z/2));
for ii = 1:n
    rr = max(1,iz(ii)-rz):min(size(label_map,1),iz(ii)+rz);
    cc = max(1,ix(ii)-rx):min(size(label_map,2),ix(ii)+rx);
    v = label_map(rr,cc);
    p(ii) = mean(v(:) == label_map(iz(ii),ix(ii)));
end
end

function y = patch_cs_stat(A, iz, ix, step_x, step_z, mode)
n = numel(iz); y = nan(n,1);
rx = max(1, round(step_x/2)); rz = max(1, round(step_z/2));
for ii = 1:n
    rr = max(1,iz(ii)-rz):min(size(A,1),iz(ii)+rz);
    cc = max(1,ix(ii)-rx):min(size(A,2),ix(ii)+rx);
    v = double(A(rr,cc)); v = v(isfinite(v));
    switch mode
        case "std", y(ii) = std(v(:));
        case "range", y(ii) = max(v(:))-min(v(:));
        case "iqr", y(ii) = prctile(v(:),75)-prctile(v(:),25);
    end
end
end

function y = sample_map(SIM, field, idx, fallback)
if isfield(SIM, field) && ~isempty(SIM.(field))
    A = SIM.(field);
    y = A(idx);
else
    y = fallback;
end
end

function roi = roi_labels(SIM, idx, n)
roi = repmat("none", n, 1);
if ~isfield(SIM, 'roi_masks') || isempty(SIM.roi_masks), return; end
names = string(fieldnames(SIM.roi_masks));
for name = names(:)'
    M = SIM.roi_masks.(name);
    if isequal(size(M), size(SIM.cs_map))
        roi(M(idx) > 0) = name;
    end
end
end

function bin = purity_bin(p)
bin = repmat("unknown", size(p));
bin(p < 0.50) = "[0,0.50)";
bin(p >= 0.50 & p < 0.75) = "[0.50,0.75)";
bin(p >= 0.75 & p < 0.90) = "[0.75,0.90)";
bin(p >= 0.90 & p < 0.95) = "[0.90,0.95)";
bin(p >= 0.95 & p < 0.99) = "[0.95,0.99)";
bin(p >= 0.99) = "[0.99,1.00]";
end

function bin = distance_window_bin(r)
bin = repmat("unknown", size(r));
bin(r < 0.25) = "0-0.25";
bin(r >= 0.25 & r < 0.5) = "0.25-0.5";
bin(r >= 0.5 & r < 1) = "0.5-1";
bin(r >= 1 & r < 2) = "1-2";
bin(r >= 2) = ">2";
end

function q = invert_mapping_to_q(mapping, k_target)
q = NaN;
if isempty(mapping) || ~isfinite(k_target), return; end
if isfield(mapping, 'k_cent')
    k = double(mapping.k_cent(:));
elseif isfield(mapping, 'k')
    k = double(mapping.k(:));
else
    return;
end
if isfield(mapping, 'Ecum')
    e = double(mapping.Ecum(:));
elseif isfield(mapping, 'cdf')
    e = double(mapping.cdf(:));
else
    return;
end
ok = isfinite(k) & isfinite(e);
k = k(ok); e = e(ok);
if numel(k) < 2, return; end
[ku, ia] = unique(k, 'stable');
q = clamp01(interp1(ku, e(ia), k_target, 'linear', 'extrap'));
end

function q = theory_q(SIM, M, cs_guess, regime, CFG)
feat = adaptive_req.config.default_feature_config('M', M, 'cs_guess', cs_guess, ...
    'gamma_win', CFG.REQ.Gamma, 'pad_factor', CFG.REQ.PadFactor);
if contains(lower(string(regime)), "diffuse")
    field_type = "Diffuse3D";
else
    field_type = "SingleWave";
end
try
    out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ(SIM.dx, SIM.dz, SIM.f0, ...
        feat.cs_guess_used, 'M', feat.M, 'Gamma', feat.gamma_win, ...
        'PadFactor', feat.pad_factor, 'Nbins', 'auto', ...
        'SmoothSigma', CFG.REQ.SmoothSigma, 'TheoryMode', 'S2D', ...
        'FieldType', field_type, 'Plot', false);
    q = out.q_th;
catch
    q = 0.5;
end
q = clamp01(q);
end

function sws = q_to_sws(mappings, q, f0)
sws = nan(numel(q),1);
for ii = 1:numel(q)
    if ~isfinite(q(ii)), continue; end
    try
        sws(ii) = adaptive_req.quantile.quantile_to_cs(mappings{ii}, q(ii), f0(ii));
    catch
        sws(ii) = NaN;
    end
end
end

function s = positive_class_score(model, score)
idx = size(score,2);
if isprop(model, 'ClassNames')
    classes = string(model.ClassNames);
    hit = find(classes == "true" | classes == "1", 1);
    if ~isempty(hit), idx = hit; end
end
s = clamp01(score(:,idx));
end

function x = clamp01(x)
x = min(max(double(x), 0), 1);
end

function y = omitmean(x)
y = mean(x, 'omitnan');
end

function y = omitmedian(x)
y = median(x, 'omitnan');
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_=-]+', '_');
end

function v = getenv_default(name, fallback)
v = getenv(name);
if isempty(v), v = fallback; end
end

function tf = env_true(name, fallback)
raw = getenv(name);
if isempty(raw), tf = fallback; return; end
tf = any(strcmpi(raw, {'1','true','yes','on'}));
end

function x = env_number(name, fallback)
raw = getenv(name);
if isempty(raw), x = fallback; return; end
x = str2double(raw);
if ~isfinite(x), x = fallback; end
end

function values = env_num_list(name, fallback)
raw = strtrim(getenv(name));
if isempty(raw), values = fallback; return; end
raw = erase(raw, '['); raw = erase(raw, ']'); raw = strrep(raw, ',', ' ');
values = sscanf(raw, '%f')';
if isempty(values), values = fallback; end
end
