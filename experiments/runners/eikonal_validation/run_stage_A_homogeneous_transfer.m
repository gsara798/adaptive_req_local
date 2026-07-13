%% run_stage_A_homogeneous_transfer.m
% Stage A: homogeneous Eikonal transfer validation for baseline_minimal_v1.
%
% This runner evaluates frozen q/SWS models on clean homogeneous Eikonal
% fields. It does not train or correct anything. The analysis companion is:
%   experiments/analysis/eikonal_validation/analyze_stage_A_homogeneous_transfer.m
%
% Runtime controls:
%   ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE = validate_only | quick | full
%   ADAPTIVE_REQ_EIKONAL_STAGE_A_FORCE_REBUILD = true | false

clear; clc; close all;
format compact;

%% Setup
this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(fileparts(this_file))));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = load_stage_config(root_dir);
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
assert(ismember(CFG.Mode, ["validate_only","quick","full"]), ...
    'ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE must be validate_only, quick, or full.');
CFG.ForceRebuild = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_A_FORCE_REBUILD', false);
CFG.CacheVersion = 3;

OUT = make_output_dirs(root_dir, CFG);
write_json_copy(CFG, fullfile(OUT.root_dir, 'stage_A_homogeneous_transfer_effective_config.json'));
waveSrc = fullfile(char(CFG.WaveSimProjectRoot), 'src');
addpath(waveSrc);
rehash;
assert(exist(fullfile(waveSrc, '+simcore', '+config', 'defaultConfig.m'), 'file') == 2, ...
    'Could not find simcore defaultConfig.m. Check WaveSimProjectRoot in config.');

fprintf('\nStage A: Homogeneous Eikonal Transfer Validation\n');
fprintf('Mode: %s | force rebuild: %d\n', CFG.Mode, CFG.ForceRebuild);
fprintf('Output root: %s\n', OUT.root_dir);
fprintf('Wave sim root: %s\n', CFG.WaveSimProjectRoot);
fprintf('Target REQ step: %.3f mm | M=%s | cs_guess=%.2f m/s\n', ...
    CFG.REQ.TargetStepM*1e3, mat2str(CFG.MList), CFG.REQ.cs_guess);

BUNDLE = load_baseline_bundle(root_dir, CFG);
assert_no_forbidden_predictors(BUNDLE.BASE_FEATURES);
fprintf('Loaded frozen models: %s\n', strjoin(string(BUNDLE.MODELS.model_names), ', '));
fprintf('Required base predictors: %d\n', numel(BUNDLE.BASE_FEATURES));

CONDS = build_condition_list(CFG);
fprintf('Stage A conditions selected: %d\n', numel(CONDS));
if CFG.Mode == "validate_only"
    fprintf('Validation-only runs one small condition end-to-end.\n');
end

%% Run conditions
patchParts = {};
sourceParts = {};
logParts = {};
mapCount = 0;
spectrumDone = containers.Map('KeyType','char','ValueType','logical');
totalTimer = tic;

for ci = 1:numel(CONDS)
    C = CONDS(ci);
    fprintf('\n[%d/%d] %s\n', ci, numel(CONDS), C.condition_id);
    cacheFile = fullfile(OUT.cache_dir, "stageA__" + sanitize(C.condition_id) + ".mat");
    conditionTimer = tic;
    try
        if ~CFG.ForceRebuild && exist(cacheFile, 'file') == 2
            S = load(cacheFile);
            if isfield(S, 'cond_signature')
                cached_signature = S.cond_signature;
            elseif isfield(S, 'condition_signature')
                cached_signature = S.condition_signature;
            else
                cached_signature = struct();
            end
            if isequaln(cached_signature, condition_signature(CFG, C))
                T_condition = S.T_condition;
                T_sources = S.T_sources;
                fprintf('  reused cache: %d patch/model rows.\n', height(T_condition));
                status = "reused";
                elapsed = toc(conditionTimer);
                patchParts{end+1,1} = T_condition; %#ok<SAGROW>
                sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
                logParts{end+1,1} = condition_log_row(C, status, height(T_condition), elapsed, ""); %#ok<SAGROW>
                continue;
            end
        end

        [sim, srcTable, simDiag] = run_eikonal_condition(C, CFG);
        quality_check_field(sim, C, CFG);
        [F, reqDiag] = extract_req_feature_table(sim, C, CFG, BUNDLE.BASE_FEATURES);
        T_condition = apply_frozen_models(F, BUNDLE, C, CFG);
        T_sources = srcTable;
        T_sources.condition_id = repmat(string(C.condition_id), height(T_sources), 1);
        T_sources = movevars(T_sources, 'condition_id', 'Before', 1);

        cond_signature = condition_signature(CFG, C); %#ok<NASGU>
        MAPS = compact_maps_for_cache(sim, T_condition, C, CFG); %#ok<NASGU>
        save(cacheFile, 'T_condition', 'T_sources', 'cond_signature', 'MAPS', 'simDiag', 'reqDiag', '-v7.3');

        if should_save_representative_map(C, CFG)
            plot_representative_maps(sim, T_condition, C, OUT);
            mapCount = mapCount + 1;
        end
        if CFG.Figures.SaveSourceGeometry
            plot_source_geometry(sim, srcTable, C, OUT);
        end
        if CFG.Figures.SaveCentralSpectrum && contains(C.field_regime, "diffuse_like")
            key = char(C.field_regime);
            if ~isKey(spectrumDone, key)
                plot_central_power_spectrum(sim, C, OUT);
                spectrumDone(key) = true;
            end
        end

        elapsed = toc(conditionTimer);
        fprintf('  completed: %d centers x %d models in %.1f s.\n', ...
            height(F), numel(BUNDLE.MODELS.model_names), elapsed);
        patchParts{end+1,1} = T_condition; %#ok<SAGROW>
        sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
        logParts{end+1,1} = condition_log_row(C, "completed", height(T_condition), elapsed, ""); %#ok<SAGROW>
    catch ME
        elapsed = toc(conditionTimer);
        warning('StageA:ConditionFailed', 'Condition failed (%s): %s', C.condition_id, ME.message);
        logParts{end+1,1} = condition_log_row(C, "failed", 0, elapsed, string(ME.message)); %#ok<SAGROW>
    end
end

if isempty(patchParts)
    error('No Stage A conditions completed successfully.');
end

T_patch = vertcat(patchParts{:});
T_source = vertcat(sourceParts{:});
T_log = vertcat(logParts{:});

%% Tables and summaries
writetable(T_patch, fullfile(OUT.table_dir, 'stage_A_patch_level_results.csv'));
writetable(T_source, fullfile(OUT.table_dir, 'stage_A_source_table.csv'));
writetable(T_log, fullfile(OUT.table_dir, 'stage_A_condition_log.csv'));

T_overall = summarize_metrics(T_patch, "model_name");
T_by_model = T_overall;
T_by_geometry = summarize_metrics(T_patch, ["model_name","geometry"]);
T_by_frequency = summarize_metrics(T_patch, ["model_name","f0"]);
T_by_field = summarize_metrics(T_patch, ["model_name","field_regime"]);
T_by_geom_freq = summarize_metrics(T_patch, ["model_name","geometry","f0"]);
T_by_geom_field = summarize_metrics(T_patch, ["model_name","geometry","field_regime"]);
T_worst = summarize_metrics(T_patch, ["model_name","condition_id","geometry","true_sws","f0","field_regime"]);
T_worst = sortrows(T_worst, 'MAPE_pct', 'descend');
T_worst = T_worst(1:min(30,height(T_worst)),:);

writetable(T_overall, fullfile(OUT.table_dir, 'stage_A_summary_overall.csv'));
writetable(T_by_model, fullfile(OUT.table_dir, 'stage_A_summary_by_model.csv'));
writetable(T_by_geometry, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'stage_A_summary_by_frequency.csv'));
writetable(T_by_field, fullfile(OUT.table_dir, 'stage_A_summary_by_field_regime.csv'));
writetable(T_by_geom_freq, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry_frequency.csv'));
writetable(T_by_geom_field, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry_field.csv'));
writetable(T_worst, fullfile(OUT.table_dir, 'stage_A_worst_conditions.csv'));
save(fullfile(OUT.data_dir, 'stage_A_homogeneous_transfer_results.mat'), ...
    'T_patch','T_source','T_log','T_overall','T_by_geometry','T_by_frequency', ...
    'T_by_field','T_by_geom_freq','T_by_geom_field','T_worst','CFG','-v7.3');

write_stage_readme(OUT, CFG, T_overall, T_by_geometry, T_by_frequency, T_by_field, T_worst, toc(totalTimer));
write_docs_readme(root_dir, CFG);

%% Console summary
nFailed = sum(string(T_log.status) == "failed");
fprintf('\nStage A runner complete.\n');
fprintf('Completed conditions: %d | failed: %d | representative maps: %d\n', ...
    sum(ismember(string(T_log.status), ["completed","reused"])), nFailed, mapCount);
print_summary_block(T_overall, 'Global summary');
print_summary_block(T_by_geometry, 'By true SWS / geometry');
print_summary_block(T_by_field, 'By field regime');
print_summary_block(T_by_frequency, 'By frequency');
if ~isempty(T_worst)
    X = T_worst(string(T_worst.model_name)=="q_spectrum_plus_composition",:);
    if ~isempty(X)
        fprintf('Worst q_spectrum_plus_composition condition: %s | MAPE %.2f%% | bias %.2f%%\n', ...
            X.condition_id(1), X.MAPE_pct(1), X.signed_bias_pct(1));
    end
end
fprintf('Tables: %s\nFigures: %s\nREADME: %s\n', OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir,'README_results.md'));

%% Local functions

function CFG = load_stage_config(root_dir)
configPath = fullfile(root_dir, 'configs', 'eikonal_validation', 'stage_A_homogeneous_transfer.json');
assert(exist(configPath, 'file') == 2, 'Missing Stage A config: %s', configPath);
CFG = jsondecode(fileread(configPath));
CFG.ConfigPath = string(configPath);
CFG.AdaptiveReqRoot = string(CFG.AdaptiveReqRoot);
CFG.WaveSimProjectRoot = string(CFG.WaveSimProjectRoot);
CFG.ModelBundlePath = string(CFG.ModelBundlePath);
CFG.OutputRoot = string(CFG.OutputRoot);
CFG.Geometries = normalize_geometry_array(CFG.Geometries);
CFG.FrequenciesHz = double(CFG.FrequenciesHz(:)).';
CFG.FieldRegimes = string(CFG.FieldRegimes(:)).';
CFG.MList = double(CFG.MList(:)).';
CFG.RealismLevels = string(CFG.RealismLevels(:)).';
CFG.REQ.Nbins = string(CFG.REQ.Nbins);
CFG.REQ.EdgeMode = string(CFG.REQ.EdgeMode);
CFG.Simulation.AmplitudeNormalization = string(CFG.Simulation.AmplitudeNormalization);
CFG.Figures.RepresentativeGeometries = string(CFG.Figures.RepresentativeGeometries(:)).';
CFG.Figures.RepresentativeFrequenciesHz = double(CFG.Figures.RepresentativeFrequenciesHz(:)).';
end

function G = normalize_geometry_array(Gin)
if isstruct(Gin)
    G = Gin(:).';
else
    error('Config Geometries must decode as a struct array.');
end
for i = 1:numel(G)
    G(i).id = string(G(i).id);
    G(i).true_sws = double(G(i).true_sws);
end
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, char(CFG.OutputRoot));
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.cache_dir = fullfile(OUT.data_dir, 'condition_cache');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.source_dir = fullfile(OUT.figure_dir, 'source_geometry');
OUT.spectrum_dir = fullfile(OUT.figure_dir, 'central_power_spectra');
ensure_dir(OUT.root_dir); ensure_dir(OUT.table_dir); ensure_dir(OUT.figure_dir);
ensure_dir(OUT.data_dir); ensure_dir(OUT.cache_dir); ensure_dir(OUT.map_dir);
ensure_dir(OUT.source_dir); ensure_dir(OUT.spectrum_dir);
end

function BUNDLE = load_baseline_bundle(root_dir, CFG)
path = fullfile(root_dir, char(CFG.ModelBundlePath));
assert(exist(path, 'file') == 2, 'Missing baseline_minimal_v1 model bundle: %s', path);
S = load(path, 'MODEL_BUNDLE');
BUNDLE = S.MODEL_BUNDLE;
assert(isfield(BUNDLE, 'MODELS') && isfield(BUNDLE, 'BASE_FEATURES'), 'Invalid model bundle.');
assert(all(ismember(["q_spectrum_only","q_spectrum_plus_composition"], string(BUNDLE.MODELS.model_names))), ...
    'Bundle does not contain both Stage A models.');
end

function assert_no_forbidden_predictors(features)
features = lower(string(features(:)));
% Keep this guard strict for oracle/evaluation variables, but do not block
% legitimate spectral names such as centroid/proxy features.
forbiddenExact = ["true_sws","sws_true","cs_true","k_true","q_oracle", ...
    "q_true","q_pred","sws_pred","signed_error_percent","abs_error_percent", ...
    "patch_purity","true_patch_purity","material_label","label_map", ...
    "distance_to_interface","distance_to_interface_m","roi","roi_label", ...
    "confidence","risk","req_mapping","map_ix","map_iz","cx","cz", ...
    "x_center_m","z_center_m","condition_id"];
forbiddenContains = ["oracle","signed_error","abs_error","high_error", ...
    "true_purity","material_side","interface_distance"];
bad = ismember(features, forbiddenExact);
for p = forbiddenContains
    bad = bad | contains(features, p);
end
if any(bad)
    error('Forbidden leakage predictor in frozen feature list: %s', strjoin(features(bad), ', '));
end
end

function CONDS = build_condition_list(CFG)
modeBlock = CFG.RunModes.(char(CFG.Mode));
geomIds = string(modeBlock.Geometries(:)).';
freqs = double(modeBlock.FrequenciesHz(:)).';
regimes = string(modeBlock.FieldRegimes(:)).';
Ms = double(modeBlock.MList(:)).';
CONDS = struct('geometry',{},'true_sws',{},'f0',{},'M',{},'field_regime',{}, ...
    'realism_level',{},'source_layout_id',{},'source_seed',{},'noise_seed',{},'condition_id',{});
idx = 0;
for gi = 1:numel(geomIds)
    G = find_geometry(CFG, geomIds(gi));
    for fi = 1:numel(freqs)
        for ri = 1:numel(regimes)
            for mi = 1:numel(Ms)
                idx = idx + 1;
                C = struct();
                C.geometry = G.id;
                C.true_sws = G.true_sws;
                C.f0 = freqs(fi);
                C.M = Ms(mi);
                C.field_regime = regimes(ri);
                C.realism_level = "clean";
                C.source_layout_id = regimes(ri);
                C.source_seed = source_seed_for(CFG, regimes(ri));
                C.noise_seed = 0;
                C.condition_id = sprintf('%s__f%d__clean__%s__M%d', ...
                    C.geometry, C.f0, C.field_regime, C.M);
                CONDS(idx) = C; %#ok<AGROW>
            end
        end
    end
end
end

function G = find_geometry(CFG, id)
ids = string({CFG.Geometries.id});
idx = find(ids == string(id), 1);
assert(~isempty(idx), 'Unknown Stage A geometry: %s', id);
G = CFG.Geometries(idx);
end

function seed = source_seed_for(CFG, regime)
S = CFG.SourceLayouts.(char(regime));
seed = double(S.source_seed);
end

function sig = condition_signature(CFG, C)
sig = struct();
sig.CacheVersion = CFG.CacheVersion;
sig.geometry = C.geometry;
sig.true_sws = C.true_sws;
sig.f0 = C.f0;
sig.M = C.M;
sig.field_regime = C.field_regime;
sig.realism_level = C.realism_level;
sig.source_seed = C.source_seed;
sig.dx = CFG.Simulation.dx_m;
sig.dz = CFG.Simulation.dz_m;
sig.TargetStepM = CFG.REQ.TargetStepM;
sig.REQ = CFG.REQ;
sig.SourceLayout = build_source_layout(C, CFG);
end

function [sim, Tsrc, D] = run_eikonal_condition(C, CFG)
cfg = simcore.config.defaultConfig();
cfg.Name = string(C.condition_id);
cfg.Lx = CFG.Simulation.Lx_m;
cfg.Lz = CFG.Simulation.Lz_m;
cfg.dx = CFG.Simulation.dx_m;
cfg.dz = CFG.Simulation.dz_m;
cfg.f0 = C.f0;
cfg.Geometry.Type = "homogeneous";
% Some simcore utilities still read legacy geometry parameters even when
% material structs are provided, so keep these explicit and harmless.
cfg.Geometry.Params = struct('cs_bg', C.true_sws, 'cs_inc', C.true_sws, ...
    'rho_bg', 1000, 'rho_inc', 1000);
mat = simcore.materials.makeMaterial("stageA_homogeneous", 'Model','elastic', ...
    'c', C.true_sws, 'alpha0', 0, 'USAlpha0', 0.5, 'BackscatterDB', 0);
cfg.Materials.Background = mat;
cfg.Materials.Inclusion = mat;
cfg.Source.Type = "point_stageA";
cfg.Source.PositionXYZ = [-10e-3 0 25e-3];
cfg.Source.ApertureMM = 0;
cfg.Source.NumPoints = 1;
cfg.Source.Axis = [0 0 1];
cfg.Source.MotionAxis = CFG.Simulation.MotionDirection(:).';
cfg.Source.Amplitude = 1;
cfg.Measurement.Axis = CFG.Simulation.MeasurementAxis(:).';
cfg.Shear.PhaseModel = string(CFG.Simulation.PhaseModel);
cfg.Shear.GeometricDecayPower = CFG.Simulation.GeometricDecayPower;
cfg.Shear.Attenuation.Enabled = false;
cfg.Shear.InterfaceShadow.Enabled = false;
cfg.Shear.InterfaceRT.Enabled = false;
cfg.Shear.Noise.Enabled = false;
cfg.Noise.Enabled = false;
cfg.Ultrasound.Enabled = false;

S = build_source_layout(C, CFG);
validate_source_layout(S, C);
cfg.SourceSet.Enabled = true;
cfg.SourceSet.PositionsXYZ = S.positions_xyz;
cfg.SourceSet.Amplitudes = S.amplitudes;
cfg.SourceSet.PhasesRad = S.phases_rad;
cfg.SourceSet.ComplexWeights = S.complex_weights;
cfg.SourceSet.MotionAxes = S.motion_axes;
cfg.SourceSet.UseParallel = false;

out = simcore.shear.simulateMultiSource2p5D(cfg);
U_pre = out.U;
preRms = sqrt(mean(abs(U_pre(:)).^2, 'omitnan'));
preMedian = median(abs(U_pre(:)), 'omitnan');
scale = 1;
if CFG.Simulation.NormalizeFieldAmplitude
    switch lower(char(CFG.Simulation.AmplitudeNormalization))
        case 'rms_to_one'
            scale = 1 / max(preRms, eps);
        case 'median_to_one'
            scale = 1 / max(preMedian, eps);
        otherwise
            error('Unsupported Stage A amplitude normalization: %s', CFG.Simulation.AmplitudeNormalization);
    end
end
out.U = U_pre .* scale;
out.U_clean = out.U_clean .* scale;
out.amplitude = abs(out.U);
out.phase = angle(out.U);
postRms = sqrt(mean(abs(out.U(:)).^2, 'omitnan'));
postMedian = median(abs(out.U(:)), 'omitnan');

sim = struct();
sim.Uxz = out.U;
sim.U = out.U;
sim.real_U = real(out.U);
sim.amplitude = abs(out.U);
sim.phase = angle(out.U);
sim.cs_map = out.medium.cs_map;
sim.label_map = out.medium.label_map;
sim.x = out.x;
sim.z = out.z;
sim.dx = cfg.dx;
sim.dz = cfg.dz;
sim.f0 = cfg.f0;
sim.source_positions_xyz = S.positions_xyz;
sim.source_phases_rad = S.phases_rad;
sim.source_amplitudes = S.amplitudes;
sim.source_motion_axes = S.motion_axes;
sim.source_layout_id = string(C.source_layout_id);
sim.amplitude_normalization_scale = scale;
sim.pre_norm_rms = preRms;
sim.pre_norm_median = preMedian;
sim.post_norm_rms = postRms;
sim.post_norm_median = postMedian;
sim.tracking_snr_db = nan(size(sim.amplitude));

Tsrc = source_table(S, C, CFG);
D = struct('pre_rms',preRms,'pre_median',preMedian,'post_rms',postRms, ...
    'post_median',postMedian,'normalization_scale',scale);
end

function S = build_source_layout(C, CFG)
Lx = CFG.Simulation.Lx_m;
Lz = CFG.Simulation.Lz_m;
cx = Lx/2; cz = Lz/2;
yOff = 18e-3;
regime = string(C.field_regime);
switch regime
    case "single_source_lateral"
        P = [-10e-3 0 cz];
        phaseSeed = 0;
        phases = 0;
    case "diffuse_like_8src_layout1"
        P = [
            -10e-3    0       cz
             Lx+10e-3 0       cz
             cx       0      -10e-3
             cx       0       Lz+10e-3
            -9e-3    -yOff   -8e-3
             Lx+9e-3  yOff    Lz+8e-3
            -9e-3     yOff    Lz+8e-3
             Lx+9e-3 -yOff   -8e-3
        ];
        phaseSeed = 6601;
        phases = deterministic_phases(size(P,1), phaseSeed);
    case "diffuse_like_8src_layout2"
        P = [
            -10e-3    0       10e-3
             Lx+10e-3 0       40e-3
             10e-3    0      -10e-3
             40e-3    0       Lz+10e-3
            -8e-3     yOff    42e-3
             Lx+8e-3 -yOff    8e-3
             6e-3    -yOff    Lz+8e-3
             44e-3    yOff   -8e-3
        ];
        phaseSeed = 6602;
        phases = deterministic_phases(size(P,1), phaseSeed);
    otherwise
        error('Unsupported Stage A field regime: %s', regime);
end
N = size(P,1);
S = struct();
S.layout_id = regime;
S.positions_xyz = P;
S.amplitudes = ones(N,1);
S.phases_rad = phases(:);
S.complex_weights = S.amplitudes .* exp(1i*S.phases_rad);
S.motion_axes = repmat(CFG.Simulation.MotionDirection(:).', N, 1);
S.source_seed = source_seed_for(CFG, regime);
S.phase_seed = phaseSeed;
end

function phases = deterministic_phases(N, seed)
rng(seed, 'twister');
phases = 2*pi*rand(N,1);
end

function validate_source_layout(S, C)
P = S.positions_xyz;
inPlane = abs(P(:,2)) < 1e-12;
assert(all(vecnorm(S.motion_axes - repmat([0 0 1], size(S.motion_axes,1), 1), 2, 2) < 1e-12), ...
    'All Stage A source motion vectors must be [0 0 1].');
if string(C.field_regime) == "single_source_lateral"
    assert(size(P,1) == 1, 'single_source_lateral must contain exactly one source.');
    assert(sum(inPlane) == 1, 'single_source_lateral must have exactly one in-plane source.');
else
    assert(size(P,1) == 8, '%s must contain exactly 8 sources.', C.field_regime);
    assert(sum(inPlane) == 4, '%s must have exactly 4 in-plane sources.', C.field_regime);
    assert(sum(~inPlane) == 4, '%s must have exactly 4 off-plane sources.', C.field_regime);
end
end

function T = source_table(S, C, CFG)
P = S.positions_xyz;
N = size(P,1);
center = [CFG.Simulation.Lx_m/2, 0, CFG.Simulation.Lz_m/2];
T = table();
T.source_index = (1:N)';
T.source_layout_id = repmat(string(C.source_layout_id), N, 1);
T.source_seed = repmat(C.source_seed, N, 1);
T.phase_seed = repmat(S.phase_seed, N, 1);
T.x_m = P(:,1); T.y_m = P(:,2); T.z_m = P(:,3);
T.x_mm = P(:,1)*1e3; T.y_mm = P(:,2)*1e3; T.z_mm = P(:,3)*1e3;
T.amplitude = S.amplitudes(:);
T.phase_rad = S.phases_rad(:);
T.motion_x = S.motion_axes(:,1);
T.motion_y = S.motion_axes(:,2);
T.motion_z = S.motion_axes(:,3);
T.is_in_plane = abs(P(:,2)) < 1e-12;
T.distance_to_domain_center_m = vecnorm(P - center, 2, 2);
T.distance_to_domain_center_mm = T.distance_to_domain_center_m * 1e3;
end

function quality_check_field(sim, C, CFG)
assert(all(isfinite(sim.U(:))), 'Eikonal field contains non-finite values.');
assert(all(isfinite(sim.cs_map(:))), 'True SWS map contains non-finite values.');
amp = sim.amplitude(:);
valid = isfinite(amp);
assert(any(valid), 'No finite amplitude values.');
lowFrac = mean(amp(valid) < 1e-8);
if lowFrac > 0.20
    warning('StageA:LowAmplitude', '%s has %.1f%% near-zero amplitude pixels.', C.condition_id, 100*lowFrac);
end
assert(abs(median(sim.cs_map(:), 'omitnan') - C.true_sws) < 1e-10, ...
    'Homogeneous true SWS map does not match condition true_sws.');
end

function [F, D] = extract_req_feature_table(sim, C, CFG, baseFeatures)
feat = adaptive_req.config.default_feature_config('M', C.M, ...
    'cs_guess', CFG.REQ.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
reqCfg = struct('f0', C.f0, 'dx', CFG.Simulation.dx_m, 'dz', CFG.Simulation.dz_m, ...
    'cs_bg', C.true_sws);
stepX = max(1, round(CFG.REQ.TargetStepM / CFG.Simulation.dx_m));
stepZ = max(1, round(CFG.REQ.TargetStepM / CFG.Simulation.dz_m));
O = adaptive_req.estimators.req_estimator_map(sim.Uxz, reqCfg, feat, ...
    'StepX', stepX, 'StepZ', stepZ, ...
    'EdgeMode', char(CFG.REQ.EdgeMode), 'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins', char(CFG.REQ.Nbins), ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, 'smooth_sigma', CFG.REQ.smooth_sigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', false, ...
    'StoreReqCurves', false, 'Verbose', false);
F = O.feature_table;
assert(height(F) > 0, 'REQ generated zero valid centers.');
F.condition_id = repmat(string(C.condition_id), height(F), 1);
F.geometry = repmat(string(C.geometry), height(F), 1);
F.true_sws = C.true_sws * ones(height(F),1);
F.SWS_true = F.true_sws;
F.f0 = C.f0 * ones(height(F),1);
F.SIM_f0 = F.f0;
F.M = C.M * ones(height(F),1);
F.REQ_M = F.M;
F.dx = CFG.Simulation.dx_m * ones(height(F),1);
F.dz = CFG.Simulation.dz_m * ones(height(F),1);
F.realism_level = repmat(string(C.realism_level), height(F), 1);
F.field_regime = repmat(string(C.field_regime), height(F), 1);
F.source_layout_id = repmat(string(C.source_layout_id), height(F), 1);
F.source_seed = C.source_seed * ones(height(F),1);
F.noise_seed = C.noise_seed * ones(height(F),1);
F.amplitude_normalization_scale = sim.amplitude_normalization_scale * ones(height(F),1);
F.pre_norm_rms = sim.pre_norm_rms * ones(height(F),1);
F.post_norm_rms = sim.post_norm_rms * ones(height(F),1);
F.REQ_StepX = stepX * ones(height(F),1);
F.REQ_StepZ = stepZ * ones(height(F),1);
F.TargetStepM = CFG.REQ.TargetStepM * ones(height(F),1);
F.local_amplitude = sample_map_at_centers(sim.amplitude, F);
F.tracking_snr_db = sample_map_at_centers(sim.tracking_snr_db, F);
F.source_to_center_distance_mm = source_to_center_distance(F, sim) * 1e3;
ensure_predictor_columns(F, baseFeatures);
D = struct('win_size', O.win_size, 'StepX', stepX, 'StepZ', stepZ, 'NumCenters', height(F));
end

function vals = sample_map_at_centers(M, F)
vals = nan(height(F),1);
if all(ismember(["cz","cx"], string(F.Properties.VariableNames)))
    idx = sub2ind(size(M), F.cz, F.cx);
    vals = M(idx);
elseif all(ismember(["map_iz","map_ix"], string(F.Properties.VariableNames)))
    idx = sub2ind(size(M), F.map_iz, F.map_ix);
    vals = M(idx);
end
end

function d = source_to_center_distance(F, sim)
if all(ismember(["x_center_m","z_center_m"], string(F.Properties.VariableNames)))
    x = F.x_center_m; z = F.z_center_m;
elseif all(ismember(["cx","cz"], string(F.Properties.VariableNames)))
    x = sim.x(F.cx).'; z = sim.z(F.cz).';
else
    x = nan(height(F),1); z = nan(height(F),1);
end
P = sim.source_positions_xyz;
d = nan(height(F),1);
for i = 1:height(F)
    if isfinite(x(i)) && isfinite(z(i))
        q = [x(i), 0, z(i)];
        d(i) = min(vecnorm(P - q, 2, 2));
    end
end
end

function ensure_predictor_columns(F, features)
vars = string(F.Properties.VariableNames);
missing = setdiff(string(features(:)), vars);
if ~isempty(missing)
    error('Required frozen baseline predictor missing from Stage A feature table: %s', strjoin(missing, ', '));
end
X = table2array(F(:, cellstr(features)));
if any(~isfinite(X), 'all')
    badCols = any(~isfinite(X), 1);
    error('Non-finite frozen predictor values for: %s', strjoin(string(features(badCols)), ', '));
end
end

function T = apply_frozen_models(F, BUNDLE, C, CFG)
MODELS = BUNDLE.MODELS;
baseFeatures = string(BUNDLE.BASE_FEATURES(:));
COMP = predict_composition(MODELS.composition, F, baseFeatures);
Faug = F;
Faug.predicted_patch_purity = COMP.predicted_patch_purity;
Faug.p_mixed = COMP.p_mixed;
Faug.p_strong_mixed = COMP.p_strong_mixed;
modelNames = ["q_spectrum_only", "q_spectrum_plus_composition"];
parts = cell(numel(modelNames),1);
for mi = 1:numel(modelNames)
    name = modelNames(mi);
    switch name
        case "q_spectrum_only"
            q = predict_q_model(MODELS.q.spectrum_only, Faug);
        case "q_spectrum_plus_composition"
            q = predict_q_model(MODELS.q.spectrum_plus_composition, Faug);
        otherwise
            error('Unsupported Stage A model: %s', name);
    end
    if any(~isfinite(q)), error('Non-finite q predictions for %s.', name); end
    if any(q < -0.05 | q > 1.05)
        warning('StageA:QRange', '%s produced q outside [-0.05, 1.05]; values are clamped for SWS conversion.', name);
    end
    q = clamp01(q);
    sws = q_to_sws(Faug.req_mapping, q, Faug.f0);
    if any(~isfinite(sws)), error('Non-finite SWS predictions for %s.', name); end
    R = table();
    R.condition_id = F.condition_id;
    R.geometry = F.geometry;
    R.true_sws = F.true_sws;
    R.f0 = F.f0;
    R.M = F.M;
    R.field_regime = F.field_regime;
    R.realism_level = F.realism_level;
    R.model_name = repmat(name, height(F), 1);
    R.q_pred = q;
    R.SWS_pred = sws;
    R.SWS_true = F.SWS_true;
    R.signed_error_percent = 100*(sws - F.SWS_true) ./ F.SWS_true;
    R.abs_error_percent = abs(R.signed_error_percent);
    R.high_error10 = R.abs_error_percent > 10;
    R.high_error20 = R.abs_error_percent > 20;
    R.local_amplitude = F.local_amplitude;
    R.tracking_snr_db = F.tracking_snr_db;
    R.source_layout_id = F.source_layout_id;
    R.source_seed = F.source_seed;
    R.noise_seed = F.noise_seed;
    R.amplitude_normalization_scale = F.amplitude_normalization_scale;
    R.pre_norm_rms = F.pre_norm_rms;
    R.post_norm_rms = F.post_norm_rms;
    R.REQ_StepX = F.REQ_StepX;
    R.REQ_StepZ = F.REQ_StepZ;
    R.TargetStepM = F.TargetStepM;
    R.REQ_Nbins_effective = F.REQ_Nbins_effective;
    R.predicted_patch_purity = Faug.predicted_patch_purity;
    R.p_mixed = Faug.p_mixed;
    R.p_strong_mixed = Faug.p_strong_mixed;
    R.source_to_center_distance_mm = F.source_to_center_distance_mm;
    if ismember("map_iz", string(F.Properties.VariableNames)), R.map_iz = F.map_iz; end
    if ismember("map_ix", string(F.Properties.VariableNames)), R.map_ix = F.map_ix; end
    if ismember("cx", string(F.Properties.VariableNames)), R.cx = F.cx; end
    if ismember("cz", string(F.Properties.VariableNames)), R.cz = F.cz; end
    if ismember("x_center_m", string(F.Properties.VariableNames)), R.x_center_m = F.x_center_m; end
    if ismember("z_center_m", string(F.Properties.VariableNames)), R.z_center_m = F.z_center_m; end
    parts{mi} = R;
end
T = vertcat(parts{:});
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

function q = predict_q_model(M, F)
X = F(:, cellstr(M.features));
q = predict(M.model, X);
end

function y = q_to_sws(mappings, q, f0)
y = nan(numel(q),1);
for i = 1:numel(q)
    if isscalar(f0), fi = f0; else, fi = f0(i); end
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), fi);
end
end

function T = summarize_metrics(T, groupVars)
if isstring(groupVars), groupVars = cellstr(groupVars); end
[G, groups] = findgroups(T(:, groupVars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G == gi,:);
    err = X.SWS_pred - X.SWS_true;
    rows{gi} = table(height(X), ...
        mean(X.abs_error_percent,'omitnan'), ...
        median(X.abs_error_percent,'omitnan'), ...
        mean(abs(err),'omitnan'), ...
        sqrt(mean(err.^2,'omitnan')), ...
        mean(X.signed_error_percent,'omitnan'), ...
        median(X.signed_error_percent,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), ...
        mean(X.SWS_pred,'omitnan'), ...
        mean(X.SWS_true,'omitnan'), ...
        mean(X.q_pred,'omitnan'), ...
        'VariableNames', {'N_valid_patches','MAPE_pct','median_abs_pct_error', ...
        'MAE_SWS','RMSE_SWS','signed_bias_pct','median_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','mean_predicted_SWS', ...
        'mean_true_SWS','mean_q_pred'});
end
T = [groups vertcat(rows{:})];
end

function row = condition_log_row(C, status, nRows, elapsed, message)
row = table(string(C.condition_id), string(C.geometry), C.true_sws, C.f0, C.M, ...
    string(C.field_regime), string(C.realism_level), string(status), nRows, elapsed, string(message), ...
    'VariableNames', {'condition_id','geometry','true_sws','f0','M','field_regime', ...
    'realism_level','status','N_patch_model_rows','elapsed_sec','message'});
end

function MAPS = compact_maps_for_cache(sim, T, C, CFG)
MAPS = struct();
MAPS.condition_id = string(C.condition_id);
MAPS.x = sim.x; MAPS.z = sim.z;
MAPS.true_sws = sim.cs_map;
MAPS.real_U = sim.real_U;
MAPS.amplitude = sim.amplitude;
MAPS.phase = sim.phase;
MAPS.sign_cos_phase = sign(cos(sim.phase));
MAPS.q_comp = rows_to_grid(T(T.model_name=="q_spectrum_plus_composition",:), 'q_pred');
MAPS.sws_comp = rows_to_grid(T(T.model_name=="q_spectrum_plus_composition",:), 'SWS_pred');
MAPS.signed_error_comp = rows_to_grid(T(T.model_name=="q_spectrum_plus_composition",:), 'signed_error_percent');
MAPS.abs_error_comp = rows_to_grid(T(T.model_name=="q_spectrum_plus_composition",:), 'abs_error_percent');
MAPS.cfg_note = 'Compact plotting cache for Stage A only.';
end

function tf = should_save_representative_map(C, CFG)
tf = CFG.Figures.SaveRepresentativeMaps && ...
    ismember(string(C.geometry), CFG.Figures.RepresentativeGeometries) && ...
    ismember(double(C.f0), CFG.Figures.RepresentativeFrequenciesHz);
end

function plot_representative_maps(sim, T, C, OUT)
X = T(T.model_name=="q_spectrum_plus_composition",:);
qMap = rows_to_grid(X, 'q_pred');
swsMap = rows_to_grid(X, 'SWS_pred');
signedMap = rows_to_grid(X, 'signed_error_percent');
absMap = rows_to_grid(X, 'abs_error_percent');
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13]);
tl = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), sim.cs_map, 'True SWS', 'm/s');
plot_map(nexttile(tl), real(sim.U), 'Real displacement', 'a.u.');
plot_map(nexttile(tl), sim.amplitude, 'Displacement amplitude', 'a.u.');
plot_map(nexttile(tl), sign(cos(sim.phase)), 'sign(cos phase)', 'sign');
plot_map(nexttile(tl), qMap, 'Predicted q', 'REQ quantile q');
plot_map(nexttile(tl), swsMap, 'Predicted SWS', 'm/s');
plot_map(nexttile(tl), signedMap, 'Signed SWS error', '%');
plot_map(nexttile(tl), absMap, 'Absolute SWS error', '%');
title(tl, pretty_condition(C), 'Interpreter','none', 'FontWeight','normal');
outDir = fullfile(OUT.map_dir, char(C.geometry), char(C.field_regime));
ensure_dir(outDir);
export_fig(fig, fullfile(outDir, "stage_A_map__" + sanitize(C.condition_id) + ".png"));
end

function plot_source_geometry(sim, Tsrc, C, OUT)
outDir = fullfile(OUT.source_dir, char(C.field_regime));
ensure_dir(outDir);
P = [Tsrc.x_m, Tsrc.y_m, Tsrc.z_m] * 1e3;
fig = figure('Color','w','Units','centimeters','Position',[1 1 25 10]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
hold(ax,'on'); grid(ax,'on'); view(ax, 35, 22);
box3d(ax, [0 50], [-20 20], [0 50]);
plane_x = [0 50 50 0]; plane_y = [0 0 0 0]; plane_z = [0 0 50 50];
patch(ax, plane_x, plane_y, plane_z, 'b', 'FaceAlpha',0.05, 'EdgeColor','b', 'LineWidth',1.2);
scatter3(ax, P(:,1), P(:,2), P(:,3), 45, 'r', 'filled');
quiver3(ax, P(:,1), P(:,2), P(:,3), Tsrc.motion_x*0, Tsrc.motion_y*0, Tsrc.motion_z*6, 0, 'k', 'LineWidth',1.3);
for i = 1:size(P,1)
    text(ax, P(i,1), P(i,2), P(i,3), sprintf(' S%d', i), 'FontSize',7);
end
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); zlabel(ax,'z (mm)');
title(ax,'3D source geometry', 'FontWeight','normal');
axis(ax,'equal');
ax = nexttile(tl);
imagesc(ax, sim.x*1e3, sim.z*1e3, sim.cs_map); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on');
scatter(ax, P(:,1), P(:,3), 45, double(Tsrc.is_in_plane), 'filled', 'MarkerEdgeColor','k');
quiver(ax, P(:,1), P(:,3), Tsrc.motion_x*0, Tsrc.motion_z*5, 0, 'k', 'LineWidth',1.2);
cb = colorbar(ax); ylabel(cb,'True SWS (m/s)'); xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
title(ax,'x-z source projection', 'FontWeight','normal');
sgtitle(fig, pretty_condition(C), 'Interpreter','none', 'FontWeight','normal');
export_fig(fig, fullfile(outDir, "stage_A_sources__" + sanitize(C.condition_id) + ".png"));
end

function box3d(ax, xr, yr, zr)
X = [xr(1) xr(2) xr(2) xr(1) xr(1)];
Y = [yr(1) yr(1) yr(2) yr(2) yr(1)];
plot3(ax, X, Y, zr(1)*ones(size(X)), 'k-');
plot3(ax, X, Y, zr(2)*ones(size(X)), 'k-');
for ix = xr
    for iy = yr
        plot3(ax, [ix ix], [iy iy], zr, 'k-');
    end
end
end

function plot_central_power_spectrum(sim, C, OUT)
Nz = size(sim.U,1); Nx = size(sim.U,2);
cz = round(Nz/2); cx = round(Nx/2); rad = min([50, cz-1, cx-1, Nz-cz, Nx-cx]);
patch = sim.U((cz-rad):(cz+rad), (cx-rad):(cx+rad));
wz = hann_local(size(patch,1)); wx = hann_local(size(patch,2));
Pw = patch .* (wz(:) * wx(:).');
P2 = abs(fftshift(fft2(Pw))).^2;
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 7]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), abs(patch), 'Central window amplitude', 'a.u.');
plot_map(nexttile(tl), log10(P2 + eps), '2D power spectrum', 'log10 power');
title(tl, "Central-window spectrum: " + pretty_condition(C), 'Interpreter','none', 'FontWeight','normal');
export_fig(fig, fullfile(OUT.spectrum_dir, "stage_A_central_spectrum__" + sanitize(C.condition_id) + ".png"));
end

function w = hann_local(N)
if N <= 1
    w = 1;
else
    n = (0:N-1)';
    w = 0.5 - 0.5*cos(2*pi*n/(N-1));
end
end

function Z = rows_to_grid(T, valueName)
if isempty(T), Z = nan; return; end
if all(ismember(["map_iz","map_ix"], string(T.Properties.VariableNames)))
    iz = T.map_iz; ix = T.map_ix;
elseif all(ismember(["cz","cx"], string(T.Properties.VariableNames)))
    [~,~,iz] = unique(T.cz); [~,~,ix] = unique(T.cx);
else
    Z = nan; return;
end
nz = max(iz); nx = max(ix);
Z = nan(nz,nx);
valid = isfinite(iz) & isfinite(ix) & iz>=1 & ix>=1;
Z(sub2ind([nz nx], iz(valid), ix(valid))) = T.(valueName)(valid);
end

function plot_map(ax, Z, ttl, cblabel)
imagesc(ax, Z); axis(ax,'image'); set(ax,'YDir','normal'); title(ax, ttl, 'FontWeight','normal');
cb = colorbar(ax); ylabel(cb, cblabel);
end

function print_summary_block(T, titleText)
fprintf('\n%s\n', titleText);
showVars = intersect(string(T.Properties.VariableNames), ...
    ["model_name","geometry","f0","field_regime","N_valid_patches","MAPE_pct","signed_bias_pct","high_error20_pct"], 'stable');
disp(T(:, cellstr(showVars)));
end

function write_stage_readme(OUT, CFG, T_overall, T_geom, T_freq, T_field, T_worst, elapsedSec)
path = fullfile(OUT.root_dir, 'README_results.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Stage A: Homogeneous Eikonal Transfer Validation\n\n');
fprintf(fid, '## Objective\n\n');
fprintf(fid, ['Validate zero-shot transfer of frozen `baseline_minimal_v1` q/SWS models ', ...
    'to clean homogeneous Eikonal simulations, without interfaces, material mixing, readout noise, ', ...
    'risk masks, corrections, or oracle q.\n\n']);
fprintf(fid, '## Scientific Questions\n\n');
fprintf(fid, ['1. Does the baseline predict SWS correctly in homogeneous Eikonal media?\n', ...
    '2. Is there a systematic clean-domain transfer bias?\n', ...
    '3. Does bias depend on true SWS: 2, 3, or 4 m/s?\n', ...
    '4. Does performance differ between directional and diffuse-like fields?\n', ...
    '5. Does performance depend on frequency?\n', ...
    '6. Does `q_spectrum_plus_composition` remain comparable or better than `q_spectrum_only`?\n\n']);
fprintf(fid, '## Simulation Design\n\n');
fprintf(fid, '- Geometries: homogeneous 2, 3, and 4 m/s.\n');
fprintf(fid, '- Frequencies: %s Hz.\n', strjoin(string(CFG.FrequenciesHz), ', '));
fprintf(fid, '- Realism: clean only. No readout noise, no shear attenuation, no interface modules.\n');
fprintf(fid, '- Measurement plane: `y = 0`; measured component: z displacement.\n');
fprintf(fid, '- Source motion is fixed to `[0 0 1]` for every source.\n');
fprintf(fid, ['- Field regimes: `single_source_lateral`, `diffuse_like_8src_layout1`, ', ...
    '`diffuse_like_8src_layout2`. Diffuse-like layouts use exactly 8 fixed sources: ', ...
    '4 in-plane and 4 off-plane.\n']);
fprintf(fid, '- Amplitude normalization: `%s`; scale and pre/post RMS are saved per patch.\n\n', CFG.Simulation.AmplitudeNormalization);
fprintf(fid, '## REQ And Model Settings\n\n');
fprintf(fid, '- M = %s.\n', mat2str(CFG.MList));
fprintf(fid, '- `cs_guess = %.2f m/s`.\n', CFG.REQ.cs_guess);
fprintf(fid, '- `dx = dz = %.1f mm`.\n', CFG.Simulation.dx_m*1e3);
fprintf(fid, '- `TargetStepM = %.1f mm`. This uses a 1 mm REQ center step to keep Stage A fast while preserving enough spatial samples for homogeneous validation.\n', CFG.REQ.TargetStepM*1e3);
fprintf(fid, '- `Nbins = auto`, `Nbins_min = %d`, `smooth_sigma = %.1f`, `Gamma = %.1f`, `PadFactor = %.1f`, `EdgeMode = %s`.\n', ...
    CFG.REQ.Nbins_min, CFG.REQ.smooth_sigma, CFG.REQ.Gamma, CFG.REQ.PadFactor, CFG.REQ.EdgeMode);
fprintf(fid, '- Frozen models: `q_spectrum_only`, `q_spectrum_plus_composition`.\n\n');
fprintf(fid, '## Metrics\n\n');
fprintf(fid, 'MAPE, signed bias, MAE/RMSE SWS, median absolute percentage error, high-error fractions >10%% and >20%%, mean predicted/true SWS, mean q, and valid patch count.\n\n');
fprintf(fid, '## Generated Figures\n\n');
fprintf(fid, ['Summary figures are generated by `experiments/analysis/eikonal_validation/analyze_stage_A_homogeneous_transfer.m`. ', ...
    'Representative maps and source geometry diagnostics are generated by the runner.\n\n']);
fprintf(fid, '## How To Run\n\n');
fprintf(fid, 'Validation-only:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''validate_only''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, 'Quick:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''quick''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, 'Full:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''full''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''full''); run(''experiments/analysis/eikonal_validation/analyze_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, '## Runtime\n\n');
fprintf(fid, 'Observed runner time for this run: %.1f minutes. Runtime scales approximately with number of conditions and valid REQ centers. Cached conditions are reused unless `ADAPTIVE_REQ_EIKONAL_STAGE_A_FORCE_REBUILD=true`.\n\n', elapsedSec/60);
fprintf(fid, '## Current Results\n\n');
fprintf(fid, '### Overall\n\n'); write_markdown_table(fid, T_overall);
fprintf(fid, '\n### By Geometry\n\n'); write_markdown_table(fid, T_geom);
fprintf(fid, '\n### By Frequency\n\n'); write_markdown_table(fid, T_freq);
fprintf(fid, '\n### By Field Regime\n\n'); write_markdown_table(fid, T_field);
fprintf(fid, '\n### Worst Conditions\n\n'); write_markdown_table(fid, T_worst(1:min(10,height(T_worst)),:));
fprintf(fid, '\n## Brief Interpretation\n\n');
fprintf(fid, '%s', interpretation_text(T_overall, T_geom, T_field, T_freq));
fclose(fid);
end

function txt = interpretation_text(T_overall, T_geom, T_field, T_freq)
txt = "";
if isempty(T_overall)
    txt = "No completed rows are available yet." + newline;
    return;
end
comp = T_overall(string(T_overall.model_name)=="q_spectrum_plus_composition",:);
only = T_overall(string(T_overall.model_name)=="q_spectrum_only",:);
if ~isempty(comp)
    txt = txt + sprintf('`q_spectrum_plus_composition` global MAPE is %.2f%% with signed bias %.2f%%. ', comp.MAPE_pct(1), comp.signed_bias_pct(1));
end
if ~isempty(comp) && ~isempty(only)
    delta = comp.MAPE_pct(1) - only.MAPE_pct(1);
    if delta <= 0
        txt = txt + sprintf('Composition is comparable or better than q-spectrum-only by %.2f points. ', abs(delta));
    else
        txt = txt + sprintf('Composition is worse than q-spectrum-only by %.2f points in this run. ', delta);
    end
end
if ~isempty(comp)
    txt = txt + "If homogeneous clean MAPE is low and bias is small across true SWS and field regimes, later Eikonal failures are more likely caused by interfaces, mixed windows, attenuation/noise, or domain complexity rather than a basic homogeneous transfer failure." + newline;
end
end

function write_docs_readme(root_dir, CFG)
docDir = fullfile(root_dir, 'docs', 'eikonal_validation');
ensure_dir(docDir);
path = fullfile(docDir, 'stage_A_homogeneous_transfer.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Stage A: Homogeneous Eikonal Transfer Validation\n\n');
fprintf(fid, 'This document mirrors the runnable Stage A configuration in `%s`.\n\n', CFG.ConfigPath);
fprintf(fid, 'Stage A evaluates the frozen `baseline_minimal_v1` models on clean homogeneous Eikonal fields only. It is intentionally separated from later interface/noise/readout stages.\n\n');
fprintf(fid, 'Run the runner from the adaptive REQ project root, then run the analysis script. The output README under `outputs/eikonal_validation/stage_A_homogeneous_transfer/README_results.md` is refreshed with observed metrics after each run.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T), fprintf(fid, '_No rows._\n'); return; end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
maxRows = min(height(T), 60);
for i = 1:maxRows
    vals = strings(1, numel(names));
    for j = 1:numel(names)
        v = T.(names(j))(i);
        if isnumeric(v) || islogical(v)
            vals(j) = sprintf('%.4g', v);
        else
            vals(j) = string(v);
        end
    end
    fprintf(fid, '| %s |\n', strjoin(vals, ' | '));
end
if height(T) > maxRows, fprintf(fid, '\n_Only first %d rows shown._\n', maxRows); end
end

function write_json_copy(CFG, path)
try
    txt = jsonencode(CFG, 'PrettyPrint', true);
catch
    txt = jsonencode(CFG);
end
fid = fopen(path, 'w'); assert(fid > 0); fwrite(fid, txt); fclose(fid);
end

function export_fig(fig, path)
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
close(fig);
end

function ensure_dir(d)
if exist(d, 'dir') ~= 7, mkdir(d); end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function s = pretty_condition(C)
s = sprintf('%s | f=%d Hz | %s | M=%d', C.geometry, C.f0, C.field_regime, C.M);
end

function y = clamp01(x)
y = min(max(x,0),1);
end

function tf = env_true(name, defaultVal)
s = lower(strtrim(string(getenv(name))));
if s == "", tf = defaultVal; return; end
tf = ismember(s, ["1","true","yes","y","on"]);
end
