%% run_stage_c_field_complexity.m
% Stage C: clean Eikonal field-complexity transfer validation.
%
% This runner evaluates frozen baseline_minimal_v1 q/SWS models on clean
% Eikonal fields while increasing source count and off-plane complexity. It
% does not train, correct, mask, or use oracle quantities as predictors.
%
% Runtime controls:
%   ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE = validate_only | quick | full
%   ADAPTIVE_REQ_EIKONAL_STAGE_C_FORCE_REBUILD = true | false
%   ADAPTIVE_REQ_EIKONAL_STAGE_C_USE_SOURCE_PARFOR = true | false
%   ADAPTIVE_REQ_EIKONAL_STAGE_C_USE_WINDOW_PARFOR = true | false

clear; clc; close all;
format compact;

%% Setup
this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(fileparts(this_file))));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.0);

CFG = load_stage_config(root_dir);
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
assert(ismember(CFG.Mode, ["validate_only","quick","full"]), ...
    'ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE must be validate_only, quick, or full.');
CFG.ForceRebuild = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_C_FORCE_REBUILD', false);
CFG.UseSourceParfor = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_C_USE_SOURCE_PARFOR', true);
CFG.UseWindowParfor = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_C_USE_WINDOW_PARFOR', false);
CFG.UseParallel = CFG.UseSourceParfor || CFG.UseWindowParfor;
CFG.CacheVersion = 1;

OUT = make_output_dirs(root_dir, CFG);
write_json_copy(CFG, fullfile(OUT.root_dir, 'stage_C_field_complexity_effective_config.json'));
waveSrc = fullfile(char(CFG.WaveSimProjectRoot), 'src');
addpath(waveSrc); rehash;
assert(exist(fullfile(waveSrc, '+simcore', '+eikonal', 'fromPointSource2p5D.m'), 'file') == 2, ...
    'Could not find simcore eikonal backend. Check WaveSimProjectRoot in config.');

fprintf('\nStage C: Field Complexity Transfer Validation\n');
fprintf('Mode: %s | force rebuild: %d\n', CFG.Mode, CFG.ForceRebuild);
fprintf('Output root: %s\n', OUT.root_dir);
fprintf('Domain: %.1f x %.1f mm | dx=dz=%.1f mm\n', ...
    CFG.Simulation.Lx_m*1e3, CFG.Simulation.Lz_m*1e3, CFG.Simulation.dx_m*1e3);
fprintf('Target REQ step: %.3f mm | M=%s | cs_guess=%.2f m/s\n', ...
    CFG.REQ.TargetStepM*1e3, mat2str(CFG.MList), CFG.REQ.cs_guess);
fprintf('Parallel: source parfor=%d | REQ window parfor=%d\n', CFG.UseSourceParfor, CFG.UseWindowParfor);
if CFG.UseParallel, ensure_parallel_pool(); end

BUNDLE = load_baseline_bundle(root_dir, CFG);
assert_no_forbidden_predictors(BUNDLE.BASE_FEATURES);
fprintf('Loaded frozen models: %s\n', strjoin(string(BUNDLE.MODELS.model_names), ', '));
fprintf('Required base predictors: %d\n', numel(BUNDLE.BASE_FEATURES));

CONDS = build_condition_list(CFG);
fprintf('Stage C conditions selected: %d\n', numel(CONDS));

%% Run conditions
patchParts = {};
sourceParts = {};
logParts = {};
mapCount = 0;
totalTimer = tic;

for ci = 1:numel(CONDS)
    C = CONDS(ci);
    fprintf('\n[%d/%d] %s\n', ci, numel(CONDS), C.condition_id);
    cacheFile = fullfile(OUT.cache_dir, "stageC__" + sanitize(C.condition_id) + ".mat");
    tCondition = tic;
    try
        sig = condition_signature(CFG, C);
        if ~CFG.ForceRebuild && exist(cacheFile, 'file') == 2
            S = load(cacheFile);
            if isfield(S, 'cond_signature') && isequaln(S.cond_signature, sig)
                T_condition = S.T_condition;
                T_sources = S.T_sources;
                elapsed = toc(tCondition);
                fprintf('  reused cache: %d patch/model rows.\n', height(T_condition));
                patchParts{end+1,1} = T_condition; %#ok<SAGROW>
                sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
                logParts{end+1,1} = condition_log_row(C, "reused", height(T_condition), elapsed, ""); %#ok<SAGROW>
                continue;
            end
        end

        [sim, Tsrc, simDiag] = run_eikonal_condition(C, CFG); %#ok<ASGLU>
        quality_check_field(sim, C);
        [F, reqDiag] = extract_req_feature_table(sim, C, CFG, BUNDLE.BASE_FEATURES); %#ok<ASGLU>
        T_condition = apply_frozen_models(F, BUNDLE);
        T_sources = Tsrc;
        T_sources.condition_id = repmat(string(C.condition_id), height(T_sources), 1);
        T_sources = movevars(T_sources, 'condition_id', 'Before', 1);

        cond_signature = sig; %#ok<NASGU>
        MAPS = compact_maps_for_cache(sim, T_condition, C, CFG); %#ok<NASGU>
        save(cacheFile, 'T_condition', 'T_sources', 'cond_signature', 'MAPS', 'simDiag', 'reqDiag', '-v7.3');

        if should_save_representative_map(C, CFG)
            plot_representative_maps(sim, T_condition, C, OUT);
            mapCount = mapCount + 1;
        end
        if CFG.Figures.SaveSourceGeometry
            plot_source_geometry(sim, Tsrc, C, OUT);
        end

        elapsed = toc(tCondition);
        fprintf('  completed: %d centers x %d models in %.1f s.\n', ...
            height(F), numel(BUNDLE.MODELS.model_names), elapsed);
        patchParts{end+1,1} = T_condition; %#ok<SAGROW>
        sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
        logParts{end+1,1} = condition_log_row(C, "completed", height(T_condition), elapsed, ""); %#ok<SAGROW>
    catch ME
        elapsed = toc(tCondition);
        warning('StageC:ConditionFailed', 'Condition failed (%s): %s', C.condition_id, ME.message);
        logParts{end+1,1} = condition_log_row(C, "failed", 0, elapsed, string(ME.message)); %#ok<SAGROW>
    end
end

if isempty(patchParts), error('No Stage C conditions completed successfully.'); end
T_patch = vertcat(patchParts{:});
T_source = vertcat(sourceParts{:});
T_log = vertcat(logParts{:});

%% Tables and summaries
writetable(T_patch, fullfile(OUT.table_dir, 'stage_c_patch_level_results.csv'));
writetable(T_source, fullfile(OUT.table_dir, 'stage_c_source_table.csv'));
writetable(T_log, fullfile(OUT.table_dir, 'stage_c_condition_log.csv'));
write_all_summaries(T_patch, OUT);
save(fullfile(OUT.data_dir, 'stage_C_field_complexity_results.mat'), 'T_patch','T_source','T_log','CFG','-v7.3');

write_stage_readme(OUT, CFG, T_log, toc(totalTimer));
write_docs_readme(root_dir, CFG);

%% Console summary
nFailed = sum(string(T_log.status) == "failed");
fprintf('\nStage C runner complete.\n');
fprintf('Completed conditions: %d | failed: %d | representative maps: %d\n', ...
    sum(ismember(string(T_log.status), ["completed","reused"])), nFailed, mapCount);
S_overall = summarize_metrics(T_patch, "model_name");
S_geom = summarize_metrics(T_patch, ["model_name","geometry"]);
S_field = summarize_metrics(T_patch, ["model_name","field_regime"]);
S_nsrc = summarize_metrics(T_patch, ["model_name","N_sources"]);
S_worst = sortrows(summarize_metrics(T_patch, ["model_name","condition_id","geometry","f0","field_regime"]), 'MAPE_pct', 'descend');
print_summary_block(S_overall, 'Global summary');
print_summary_block(S_geom, 'By geometry');
print_summary_block(S_field, 'By field regime');
print_summary_block(S_nsrc, 'By source count');
X = S_worst(string(S_worst.model_name)=="q_spectrum_plus_composition",:);
if ~isempty(X)
    fprintf('Worst q_spectrum_plus_composition condition: %s | MAPE %.2f%% | bias %.2f%%\n', ...
        X.condition_id(1), X.MAPE_pct(1), X.signed_bias_pct(1));
end
fprintf('Tables: %s\nFigures: %s\nREADME: %s\n', OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir,'README_results.md'));

%% Local functions
function ensure_parallel_pool()
try
    pool = gcp('nocreate');
    if isempty(pool)
        try, parpool('threads'); catch, parpool; end
    end
catch ME
    warning('StageC:ParallelPool', 'Could not start a parallel pool. Continuing serially: %s', ME.message);
end
end

function CFG = load_stage_config(root_dir)
configPath = fullfile(root_dir, 'configs', 'eikonal_validation', 'stage_c_field_complexity.json');
assert(exist(configPath, 'file') == 2, 'Missing Stage C config: %s', configPath);
CFG = jsondecode(fileread(configPath));
CFG.ConfigPath = string(configPath);
CFG.AdaptiveReqRoot = string(CFG.AdaptiveReqRoot);
CFG.WaveSimProjectRoot = string(CFG.WaveSimProjectRoot);
CFG.ModelBundlePath = string(CFG.ModelBundlePath);
CFG.OutputRoot = string(CFG.OutputRoot);
CFG.FrequenciesHz = double(CFG.FrequenciesHz(:)).';
CFG.FieldRegimes = string(CFG.FieldRegimes(:)).';
CFG.MList = double(CFG.MList(:)).';
CFG.RealismLevels = string(CFG.RealismLevels(:)).';
CFG.REQ.Nbins = string(CFG.REQ.Nbins);
CFG.REQ.EdgeMode = string(CFG.REQ.EdgeMode);
CFG.Simulation.AmplitudeNormalization = string(CFG.Simulation.AmplitudeNormalization);
CFG.Geometries = CFG.Geometries(:).';
for i = 1:numel(CFG.Geometries)
    CFG.Geometries(i).id = string(CFG.Geometries(i).id);
    CFG.Geometries(i).type = string(CFG.Geometries(i).type);
end
CFG.Figures.RepresentativeGeometries = string(CFG.Figures.RepresentativeGeometries(:)).';
CFG.Figures.RepresentativeFrequenciesHz = double(CFG.Figures.RepresentativeFrequenciesHz(:)).';
CFG.Figures.RepresentativeFieldRegimes = string(CFG.Figures.RepresentativeFieldRegimes(:)).';
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, char(CFG.OutputRoot));
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.cache_dir = fullfile(OUT.data_dir, 'condition_cache');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.source_dir = fullfile(OUT.figure_dir, 'source_geometry');
ensure_dir(OUT.root_dir); ensure_dir(OUT.table_dir); ensure_dir(OUT.figure_dir);
ensure_dir(OUT.data_dir); ensure_dir(OUT.cache_dir); ensure_dir(OUT.map_dir); ensure_dir(OUT.source_dir);
end

function BUNDLE = load_baseline_bundle(root_dir, CFG)
path = fullfile(root_dir, char(CFG.ModelBundlePath));
assert(exist(path, 'file') == 2, 'Missing baseline_minimal_v1 model bundle: %s', path);
S = load(path, 'MODEL_BUNDLE');
BUNDLE = S.MODEL_BUNDLE;
assert(isfield(BUNDLE, 'MODELS') && isfield(BUNDLE, 'BASE_FEATURES'), 'Invalid model bundle.');
assert(all(ismember(["q_spectrum_only","q_spectrum_plus_composition"], string(BUNDLE.MODELS.model_names))), ...
    'Bundle does not contain both Stage C models.');
end

function assert_no_forbidden_predictors(features)
features = lower(string(features(:)));
forbiddenExact = ["true_sws","sws_true","cs_true","k_true","q_oracle","q_true", ...
    "q_pred","sws_pred","signed_error_percent","abs_error_percent","patch_purity", ...
    "true_patch_purity","material_label","label_map","distance_to_interface", ...
    "distance_to_interface_m","roi","roi_label","confidence","risk","req_mapping", ...
    "map_ix","map_iz","cx","cz","x_center_m","z_center_m","condition_id"];
forbiddenContains = ["oracle","signed_error","abs_error","high_error","true_purity", ...
    "material_side","interface_distance"];
bad = ismember(features, forbiddenExact);
for p = forbiddenContains, bad = bad | contains(features, p); end
if any(bad), error('Forbidden leakage predictor in frozen feature list: %s', strjoin(features(bad), ', ')); end
end

function CONDS = build_condition_list(CFG)
modeBlock = CFG.RunModes.(char(CFG.Mode));
geomIds = string(modeBlock.Geometries(:)).';
freqs = double(modeBlock.FrequenciesHz(:)).';
regimes = string(modeBlock.FieldRegimes(:)).';
Ms = double(modeBlock.MList(:)).';
CONDS = struct('geometry',{},'geometry_type',{},'background_sws',{},'inclusion_sws',{}, ...
    'inclusion_diameter_mm',{},'f0',{},'M',{},'field_regime',{},'realism_level',{}, ...
    'source_layout_id',{},'source_seed',{},'noise_seed',{},'condition_id',{});
idx = 0;
for gi = 1:numel(geomIds)
    G = find_geometry(CFG, geomIds(gi));
    for fi = 1:numel(freqs)
        for ri = 1:numel(regimes)
            for mi = 1:numel(Ms)
                idx = idx + 1;
                C = struct();
                C.geometry = G.id;
                C.geometry_type = G.type;
                C.background_sws = double(G.background_sws);
                C.inclusion_sws = double(G.inclusion_sws);
                C.inclusion_diameter_mm = double(G.inclusion_diameter_mm);
                C.f0 = freqs(fi);
                C.M = Ms(mi);
                C.field_regime = regimes(ri);
                C.realism_level = "clean";
                C.source_layout_id = regimes(ri);
                C.source_seed = source_seed_for(CFG, regimes(ri));
                C.noise_seed = 0;
                C.condition_id = sprintf('%s__f%d__clean__%s__M%d', C.geometry, C.f0, C.field_regime, C.M);
                CONDS(idx) = C; %#ok<AGROW>
            end
        end
    end
end
end

function G = find_geometry(CFG, id)
ids = string({CFG.Geometries.id});
idx = find(ids == string(id), 1);
assert(~isempty(idx), 'Unknown Stage C geometry: %s', id);
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
sig.geometry_type = C.geometry_type;
sig.background_sws = C.background_sws;
sig.inclusion_sws = C.inclusion_sws;
sig.inclusion_diameter_mm = C.inclusion_diameter_mm;
sig.f0 = C.f0;
sig.M = C.M;
sig.field_regime = C.field_regime;
sig.realism_level = C.realism_level;
sig.source_seed = C.source_seed;
sig.dx = CFG.Simulation.dx_m;
sig.dz = CFG.Simulation.dz_m;
sig.Lx = CFG.Simulation.Lx_m;
sig.Lz = CFG.Simulation.Lz_m;
sig.TargetStepM = CFG.REQ.TargetStepM;
sig.REQ = CFG.REQ;
sig.SourceLayout = build_source_layout(C, CFG);
end

function [sim, Tsrc, D] = run_eikonal_condition(C, CFG)
phant = build_stage_c_geometry(C, CFG);
S = build_source_layout(C, CFG);
validate_source_layout(S, C);
U_pre = synthesize_clean_eikonal_field(phant, S, C, CFG);
validMask = central_valid_mask(phant, CFG);
preRms = sqrt(mean(abs(U_pre(validMask)).^2, 'omitnan'));
preMedian = median(abs(U_pre(validMask)), 'omitnan');
scale = 1;
if CFG.Simulation.NormalizeFieldAmplitude
    switch lower(char(CFG.Simulation.AmplitudeNormalization))
        case 'central_rms_to_one'
            scale = 1 / max(preRms, eps);
        case 'central_median_to_one'
            scale = 1 / max(preMedian, eps);
        otherwise
            error('Unsupported Stage C amplitude normalization: %s', CFG.Simulation.AmplitudeNormalization);
    end
end
U = U_pre .* scale;
postRms = sqrt(mean(abs(U(validMask)).^2, 'omitnan'));
postMedian = median(abs(U(validMask)), 'omitnan');

sim = struct();
sim.Uxz = U;
sim.U = U;
sim.real_U = real(U);
sim.amplitude = abs(U);
sim.phase = angle(U);
sim.cs_map = phant.cs_map;
sim.label_map = phant.label_map;
sim.signed_distance_m = phant.signed_distance_m;
sim.abs_distance_m = phant.abs_distance_m;
sim.x = phant.x;
sim.z = phant.z;
sim.dx = CFG.Simulation.dx_m;
sim.dz = CFG.Simulation.dz_m;
sim.f0 = C.f0;
sim.source_positions_xyz = S.positions_xyz;
sim.source_phases_rad = S.phases_rad;
sim.source_amplitudes = S.amplitudes;
sim.source_motion_axes = S.motion_axes;
sim.source_layout_id = string(C.source_layout_id);
sim.N_sources = size(S.positions_xyz,1);
sim.N_in_plane_sources = sum(abs(S.positions_xyz(:,2)) < 1e-12);
sim.layout_family = S.layout_family;
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

function phant = build_stage_c_geometry(C, CFG)
x = 0:CFG.Simulation.dx_m:CFG.Simulation.Lx_m;
z = 0:CFG.Simulation.dz_m:CFG.Simulation.Lz_m;
[X,Z] = meshgrid(x,z);
cs = C.background_sws * ones(size(X));
label = zeros(size(X));
center = [CFG.Simulation.Lx_m/2, CFG.Simulation.Lz_m/2];
infMap = inf(size(X));
signedD = infMap; absD = infMap;
if C.geometry_type == "homogeneous"
    roiType = zeros(size(X)); %#ok<NASGU>
elseif C.geometry_type == "inclusion"
    r = C.inclusion_diameter_mm/2 * 1e-3;
    distCenter = hypot(X-center(1), Z-center(2));
    inside = distCenter <= r;
    cs(inside) = C.inclusion_sws;
    label(inside) = 1;
    signedD = r - distCenter;
    absD = abs(signedD);
elseif C.geometry_type == "bilayer"
    interfaceX = CFG.Simulation.Lx_m/2;
    hard = X >= interfaceX;
    cs(hard) = C.inclusion_sws;
    label(hard) = 1;
    signedD = X - interfaceX;
    absD = abs(signedD);
else
    error('Unsupported Stage C geometry type: %s', C.geometry_type);
end
phant = struct('x',x,'z',z,'X',X,'Z',Z,'cs_map',cs,'label_map',label, ...
    'signed_distance_m',signedD,'abs_distance_m',absD);
end

function mask = central_valid_mask(phant, CFG)
margin = CFG.Simulation.ValidRegionMarginM;
mask = phant.X >= margin & phant.X <= max(phant.x)-margin & ...
       phant.Z >= margin & phant.Z <= max(phant.z)-margin;
end

function U = synthesize_clean_eikonal_field(phant, S, C, CFG)
omega = 2*pi*C.f0;
Nsrc = size(S.positions_xyz,1);
parts = cell(Nsrc,1);
usePar = CFG.UseSourceParfor && Nsrc > 1;
if usePar
    parfor si = 1:Nsrc
        parts{si} = synthesize_one_source(phant, S, C, CFG, omega, si);
    end
else
    for si = 1:Nsrc
        parts{si} = synthesize_one_source(phant, S, C, CFG, omega, si);
    end
end
U = zeros(size(phant.cs_map));
for si = 1:Nsrc
    U = U + parts{si};
end
end

function Usi = synthesize_one_source(phant, S, C, CFG, omega, si)
p = S.positions_xyz(si,:);
opts = struct();
opts.ReferenceSpeed = median(phant.cs_map(:));
opts.MaxIterations = 220;
opts.Tolerance = 1e-5;
opts.NormalizeBoundaryTime = true;
[T,~] = simcore.eikonal.fromPointSource2p5D(phant.cs_map, phant.x, phant.z, 0, p, opts);
R = sqrt((phant.X-p(1)).^2 + (0-p(2)).^2 + (phant.Z-p(3)).^2);
R = max(R, eps);
R0 = max(min(R(:)), CFG.Simulation.dx_m);
geomAmp = S.amplitudes(si) ./ max((R./R0).^CFG.Simulation.GeometricDecayPower, eps);
[~,~,pz] = transverse_polarization(phant, p, S.motion_axes(si,:));
Usi = geomAmp .* pz .* exp(1i*(omega*T + S.phases_rad(si)));
end

function [px,py,pz] = transverse_polarization(phant, p, motionAxis)
motionAxis = motionAxis(:).' ./ norm(motionAxis);
dxv = phant.X - p(1); dyv = zeros(size(phant.X)) - p(2); dzv = phant.Z - p(3);
nrm = max(sqrt(dxv.^2 + dyv.^2 + dzv.^2), eps);
nx = dxv./nrm; ny = dyv./nrm; nz = dzv./nrm;
dotMN = motionAxis(1).*nx + motionAxis(2).*ny + motionAxis(3).*nz;
px = motionAxis(1) - dotMN.*nx;
py = motionAxis(2) - dotMN.*ny;
pz = motionAxis(3) - dotMN.*nz;
pnorm = max(sqrt(px.^2 + py.^2 + pz.^2), eps);
px = px./pnorm; py = py./pnorm; pz = pz./pnorm;
end

function S = build_source_layout(C, CFG)
regime = string(C.field_regime);
Lx = CFG.Simulation.Lx_m; Lz = CFG.Simulation.Lz_m;
cz = Lz/2;
if regime == "single_source_lateral"
    P = [-12e-3 0 cz];
    phases = 0;
    phaseSeed = 0;
    layoutFamily = "single";
else
    L = CFG.SourceLayouts.(char(regime));
    N = double(L.N_sources);
    Nin = double(L.N_in_plane);
    layoutFamily = string(L.layout_family);
    P = shell_source_positions(N, Nin, layoutFamily, CFG);
    phaseSeed = double(L.phase_seed);
    phases = deterministic_phases(N, phaseSeed);
end
Nsrc = size(P,1);
S = struct();
S.layout_id = regime;
S.layout_family = layoutFamily;
S.positions_xyz = P;
S.amplitudes = ones(Nsrc,1);
S.phases_rad = phases(:);
S.complex_weights = S.amplitudes .* exp(1i*S.phases_rad);
S.motion_axes = repmat(CFG.Simulation.MotionDirection(:).', Nsrc, 1);
S.source_seed = source_seed_for(CFG, regime);
S.phase_seed = phaseSeed;
end

function P = shell_source_positions(N, Nin, family, CFG)
Lx = CFG.Simulation.Lx_m; Lz = CFG.Simulation.Lz_m;
cx = Lx/2; cz = Lz/2; pad = 13e-3; yOff = 26e-3;
inA = [-pad 0 cz; Lx+pad 0 cz; cx 0 -pad; cx 0 Lz+pad];
inB = [-pad 0 0.30*Lz; Lx+pad 0 0.70*Lz; 0.30*Lx 0 -pad; 0.70*Lx 0 Lz+pad];
if family == "A"
    inPool = inA;
    theta0 = pi/8;
else
    inPool = inB;
    theta0 = pi/4;
end
P = inPool(1:Nin,:);
Noff = N - Nin;
if Noff > 0
    off = zeros(Noff,3);
    R = 0.5*sqrt(Lx^2 + Lz^2) + pad;
    for i = 1:Noff
        th = theta0 + 2*pi*(i-1)/Noff;
        off(i,1) = cx + R*cos(th);
        off(i,3) = cz + R*sin(th);
        off(i,2) = yOff * (-1)^(i + (family=="B"));
    end
    P = [P; off];
end
end

function phases = deterministic_phases(N, seed)
rng(seed, 'twister');
phases = 2*pi*rand(N,1);
end

function validate_source_layout(S, C)
P = S.positions_xyz;
inPlane = abs(P(:,2)) < 1e-12;
assert(all(vecnorm(S.motion_axes - repmat([0 0 1], size(S.motion_axes,1), 1), 2, 2) < 1e-12), ...
    'All Stage C source motion vectors must be [0 0 1].');
if string(C.field_regime) == "single_source_lateral"
    assert(size(P,1) == 1 && sum(inPlane) == 1, 'single_source_lateral must have one in-plane source.');
else
    parts = split(string(C.field_regime), "_");
    nstr = erase(parts(3), "src");
    expectedN = str2double(nstr);
    assert(size(P,1) == expectedN, '%s must contain %d sources.', C.field_regime, expectedN);
    if endsWith(string(C.field_regime), "layoutA")
        expectedIn = 1;
    elseif contains(string(C.field_regime), "4src")
        expectedIn = 1;
    elseif contains(string(C.field_regime), "8src")
        expectedIn = 2;
    else
        expectedIn = 4;
    end
    assert(sum(inPlane) == expectedIn, '%s must have %d in-plane sources.', C.field_regime, expectedIn);
end
end

function T = source_table(S, C, CFG)
P = S.positions_xyz;
N = size(P,1);
center = [CFG.Simulation.Lx_m/2, 0, CFG.Simulation.Lz_m/2];
T = table();
T.source_index = (1:N)';
T.source_layout_id = repmat(string(C.source_layout_id), N, 1);
T.layout_family = repmat(string(S.layout_family), N, 1);
T.N_sources = repmat(N, N, 1);
T.N_in_plane_sources = repmat(sum(abs(P(:,2))<1e-12), N, 1);
T.source_seed = repmat(C.source_seed, N, 1);
T.phase_seed = repmat(S.phase_seed, N, 1);
T.x_m = P(:,1); T.y_m = P(:,2); T.z_m = P(:,3);
T.x_mm = P(:,1)*1e3; T.y_mm = P(:,2)*1e3; T.z_mm = P(:,3)*1e3;
T.amplitude = S.amplitudes(:);
T.phase_rad = S.phases_rad(:);
T.motion_x = S.motion_axes(:,1); T.motion_y = S.motion_axes(:,2); T.motion_z = S.motion_axes(:,3);
T.is_in_plane = abs(P(:,2)) < 1e-12;
T.distance_to_domain_center_m = vecnorm(P - center, 2, 2);
T.distance_to_domain_center_mm = T.distance_to_domain_center_m * 1e3;
if min(T.distance_to_domain_center_mm) < 12
    warning('StageC:SourceTooClose', '%s has a source close to the domain center.', C.condition_id);
end
end

function quality_check_field(sim, C)
assert(all(isfinite(sim.U(:))), 'Eikonal field contains non-finite values.');
assert(all(isfinite(sim.cs_map(:))), 'True SWS map contains non-finite values.');
amp = sim.amplitude(:);
lowFrac = mean(amp(isfinite(amp)) < 1e-8);
if lowFrac > 0.20
    warning('StageC:LowAmplitude', '%s has %.1f%% near-zero amplitude pixels.', C.condition_id, 100*lowFrac);
end
end

function [F, D] = extract_req_feature_table(sim, C, CFG, baseFeatures)
feat = adaptive_req.config.default_feature_config('M', C.M, ...
    'cs_guess', CFG.REQ.cs_guess, 'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
reqCfg = struct('f0', C.f0, 'dx', CFG.Simulation.dx_m, 'dz', CFG.Simulation.dz_m, ...
    'cs_bg', C.background_sws);
stepX = max(1, round(CFG.REQ.TargetStepM / CFG.Simulation.dx_m));
stepZ = max(1, round(CFG.REQ.TargetStepM / CFG.Simulation.dz_m));
O = adaptive_req.estimators.req_estimator_map(sim.Uxz, reqCfg, feat, ...
    'StepX', stepX, 'StepZ', stepZ, ...
    'EdgeMode', char(CFG.REQ.EdgeMode), 'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins', char(CFG.REQ.Nbins), ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, 'smooth_sigma', CFG.REQ.smooth_sigma}, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, ...
    'StoreReqCurves', false, 'Verbose', false, ...
    'UseWindowParfor', CFG.UseWindowParfor);
F = O.feature_table;
assert(height(F) > 0, 'REQ generated zero valid centers.');
F.condition_id = repmat(string(C.condition_id), height(F), 1);
F.geometry = repmat(string(C.geometry), height(F), 1);
F.geometry_type = repmat(string(C.geometry_type), height(F), 1);
F.background_sws = C.background_sws * ones(height(F),1);
F.inclusion_sws = C.inclusion_sws * ones(height(F),1);
F.inclusion_diameter_mm = C.inclusion_diameter_mm * ones(height(F),1);
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
F.N_sources = sim.N_sources * ones(height(F),1);
F.N_in_plane_sources = sim.N_in_plane_sources * ones(height(F),1);
F.layout_family = repmat(string(sim.layout_family), height(F), 1);
F.amplitude_normalization_scale = sim.amplitude_normalization_scale * ones(height(F),1);
F.pre_norm_rms = sim.pre_norm_rms * ones(height(F),1);
F.post_norm_rms = sim.post_norm_rms * ones(height(F),1);
F.REQ_StepX = stepX * ones(height(F),1);
F.REQ_StepZ = stepZ * ones(height(F),1);
F.TargetStepM = CFG.REQ.TargetStepM * ones(height(F),1);
F.local_amplitude = sample_map_at_centers(sim.amplitude, F);
F.tracking_snr_db = sample_map_at_centers(sim.tracking_snr_db, F);
F.SWS_true = sample_map_at_centers(sim.cs_map, F);
F.true_sws = F.SWS_true;
F.material_label = sample_map_at_centers(sim.label_map, F);
F.material_label = round(F.material_label);
F.distance_to_interface_m = sample_map_at_centers(sim.abs_distance_m, F);
F.distance_to_interface_m(~isfinite(F.distance_to_interface_m)) = inf;
winRadiusM = 0.5 * C.M * CFG.REQ.cs_guess / C.f0;
F.distance_to_interface_over_window_radius = F.distance_to_interface_m ./ max(winRadiusM, eps);
F.patch_purity = compute_patch_purity(sim.label_map, F, C, CFG);
[F.patch_cs_std, F.patch_cs_range, F.patch_cs_iqr] = compute_patch_cs_stats(sim.cs_map, F, C, CFG);
F.M_eff = C.M * CFG.REQ.cs_guess ./ F.SWS_true;
F.lambda_over_Lwin = 1 ./ F.M_eff;
F.wavelengths_per_window = F.M_eff;
F.roi_label = assign_roi_labels(F, C, CFG);
F.source_to_center_distance_mm = source_to_center_distance(F, sim) * 1e3;
ensure_predictor_columns(F, baseFeatures);
D = struct('win_size', O.win_size, 'StepX', stepX, 'StepZ', stepZ, 'NumCenters', height(F));
end

function vals = sample_map_at_centers(M, F)
vals = nan(height(F),1);
if all(ismember(["cz","cx"], string(F.Properties.VariableNames)))
    idx = sub2ind(size(M), F.cz, F.cx); vals = M(idx);
elseif all(ismember(["map_iz","map_ix"], string(F.Properties.VariableNames)))
    idx = sub2ind(size(M), F.map_iz, F.map_ix); vals = M(idx);
end
end

function purity = compute_patch_purity(labelMap, F, C, CFG)
winM = C.M * CFG.REQ.cs_guess / C.f0;
rx = max(1, round(0.5*winM / CFG.Simulation.dx_m));
rz = max(1, round(0.5*winM / CFG.Simulation.dz_m));
purity = nan(height(F),1);
for i = 1:height(F)
    [cz,cx] = center_indices(F, i);
    zidx = max(1,cz-rz):min(size(labelMap,1),cz+rz);
    xidx = max(1,cx-rx):min(size(labelMap,2),cx+rx);
    patch = labelMap(zidx,xidx);
    labs = unique(patch(:));
    counts = arrayfun(@(v) sum(patch(:)==v), labs);
    purity(i) = max(counts) / numel(patch);
end
end

function [csStd, csRange, csIqr] = compute_patch_cs_stats(csMap, F, C, CFG)
winM = C.M * CFG.REQ.cs_guess / C.f0;
rx = max(1, round(0.5*winM / CFG.Simulation.dx_m));
rz = max(1, round(0.5*winM / CFG.Simulation.dz_m));
csStd = nan(height(F),1); csRange = nan(height(F),1); csIqr = nan(height(F),1);
for i = 1:height(F)
    [cz,cx] = center_indices(F, i);
    zidx = max(1,cz-rz):min(size(csMap,1),cz+rz);
    xidx = max(1,cx-rx):min(size(csMap,2),cx+rx);
    vals = csMap(zidx,xidx); vals = vals(:);
    csStd(i) = std(vals, 'omitnan');
    csRange(i) = max(vals) - min(vals);
    csIqr(i) = iqr(vals);
end
end

function [cz,cx] = center_indices(F, i)
if all(ismember(["cz","cx"], string(F.Properties.VariableNames)))
    cz = F.cz(i); cx = F.cx(i);
else
    cz = F.map_iz(i); cx = F.map_ix(i);
end
end

function roi = assign_roi_labels(F, C, CFG)
roi = repmat("other", height(F), 1);
dmm = F.distance_to_interface_m * 1e3;
if C.geometry_type == "homogeneous"
    cx0 = CFG.Simulation.Lx_m/2; cz0 = CFG.Simulation.Lz_m/2;
    hw = CFG.ROI.HomogeneousCenterHalfWidthMm * 1e-3;
    [x,z] = centers_from_F(F);
    centerMask = abs(x-cx0) <= hw & abs(z-cz0) <= hw;
    roi(:) = "homogeneous_other";
    roi(centerMask) = "homogeneous_center";
elseif C.geometry_type == "inclusion"
    inside = F.SWS_true > F.background_sws + 0.1;
    farBg = ~inside & dmm >= CFG.ROI.BackgroundFarMinDistanceMm;
    roi(farBg) = "background_far";
    for bi = 1:size(CFG.ROI.InterfaceBandsMm,1)
        lo = CFG.ROI.InterfaceBandsMm(bi,1); hi = CFG.ROI.InterfaceBandsMm(bi,2);
        mask = dmm >= lo & dmm < hi;
        roi(mask) = sprintf('interface_%g_%gmm', lo, hi);
    end
    centerR = 0.5 * C.inclusion_diameter_mm/2; % central half-radius in mm
    centerMask = inside & dmm >= centerR;
    roi(centerMask) = "inclusion_center";
    core = inside & dmm >= CFG.ROI.CoreMinDistanceMm & F.patch_purity >= 0.95;
    roi(core) = "inclusion_core";
elseif C.geometry_type == "bilayer"
    soft = F.SWS_true <= F.background_sws + 0.1;
    hard = ~soft;
    roi(soft & dmm >= CFG.ROI.BackgroundFarMinDistanceMm) = "soft_layer_far";
    roi(hard & dmm >= CFG.ROI.BackgroundFarMinDistanceMm) = "hard_layer_far";
    for bi = 1:size(CFG.ROI.InterfaceBandsMm,1)
        lo = CFG.ROI.InterfaceBandsMm(bi,1); hi = CFG.ROI.InterfaceBandsMm(bi,2);
        mask = dmm >= lo & dmm < hi;
        roi(mask) = sprintf('interface_%g_%gmm', lo, hi);
    end
end
end

function [x,z] = centers_from_F(F)
if all(ismember(["x_center_m","z_center_m"], string(F.Properties.VariableNames)))
    x = F.x_center_m; z = F.z_center_m;
else
    x = nan(height(F),1); z = nan(height(F),1);
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
    error('Required frozen baseline predictor missing from Stage C feature table: %s', strjoin(missing, ', '));
end
X = table2array(F(:, cellstr(features)));
if any(~isfinite(X), 'all')
    badCols = any(~isfinite(X), 1);
    error('Non-finite frozen predictor values for: %s', strjoin(string(features(badCols)), ', '));
end
end

function T = apply_frozen_models(F, BUNDLE)
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
    end
    if any(~isfinite(q)), error('Non-finite q predictions for %s.', name); end
    if any(q < -0.05 | q > 1.05)
        warning('StageC:QRange', '%s produced q outside [-0.05, 1.05]; clamping for SWS conversion.', name);
    end
    q = clamp01(q);
    sws = q_to_sws(Faug.req_mapping, q, Faug.f0);
    if any(~isfinite(sws)), error('Non-finite SWS predictions for %s.', name); end
    R = table();
    keep = ["condition_id","geometry","geometry_type","background_sws","inclusion_sws", ...
        "inclusion_diameter_mm","f0","M","field_regime","realism_level", ...
        "source_layout_id","source_seed","noise_seed","N_sources","N_in_plane_sources", ...
        "layout_family","patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr", ...
        "distance_to_interface_m","distance_to_interface_over_window_radius", ...
        "M_eff","lambda_over_Lwin","wavelengths_per_window","roi_label", ...
        "local_amplitude","tracking_snr_db","source_to_center_distance_mm", ...
        "amplitude_normalization_scale","pre_norm_rms","post_norm_rms", ...
        "REQ_StepX","REQ_StepZ","TargetStepM","REQ_Nbins_effective"];
    for v = keep
        if ismember(v, string(F.Properties.VariableNames)), R.(v) = F.(v); end
    end
    R.model_name = repmat(name, height(F), 1);
    R.q_pred = q;
    R.SWS_pred = sws;
    R.k_pred = 2*pi*F.f0 ./ sws;
    R.SWS_true = F.SWS_true;
    R.signed_error_pct = 100*(sws - F.SWS_true) ./ F.SWS_true;
    R.abs_error_pct = abs(R.signed_error_pct);
    R.high_error10 = R.abs_error_pct > 10;
    R.high_error20 = R.abs_error_pct > 20;
    R.mean_true_sws = F.SWS_true;
    R.predicted_patch_purity = Faug.predicted_patch_purity;
    R.p_mixed = Faug.p_mixed;
    R.p_strong_mixed = Faug.p_strong_mixed;
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
    fi = f0(i);
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), fi);
end
end

function write_all_summaries(T, OUT)
writetable(summarize_metrics(T, "model_name"), fullfile(OUT.table_dir, 'stage_c_summary_overall.csv'));
writetable(summarize_metrics(T, ["model_name","geometry"]), fullfile(OUT.table_dir, 'stage_c_summary_by_geometry.csv'));
writetable(summarize_metrics(T, ["model_name","field_regime"]), fullfile(OUT.table_dir, 'stage_c_summary_by_field_regime.csv'));
writetable(summarize_metrics(T, ["model_name","N_sources"]), fullfile(OUT.table_dir, 'stage_c_summary_by_source_count.csv'));
writetable(summarize_metrics(T, ["model_name","N_in_plane_sources"]), fullfile(OUT.table_dir, 'stage_c_summary_by_in_plane_source_count.csv'));
writetable(summarize_metrics(T, ["model_name","f0"]), fullfile(OUT.table_dir, 'stage_c_summary_by_frequency.csv'));
writetable(summarize_metrics(T, ["model_name","geometry","field_regime"]), fullfile(OUT.table_dir, 'stage_c_summary_by_geometry_field.csv'));
writetable(summarize_metrics(T, ["model_name","geometry","f0"]), fullfile(OUT.table_dir, 'stage_c_summary_by_geometry_frequency.csv'));
writetable(summarize_metrics(T, ["model_name","roi_label"]), fullfile(OUT.table_dir, 'stage_c_summary_by_roi.csv'));
writetable(summarize_metrics(T, ["model_name","roi_label","field_regime"]), fullfile(OUT.table_dir, 'stage_c_summary_by_roi_field.csv'));
writetable(summarize_metrics(T, ["model_name","roi_label","N_sources"]), fullfile(OUT.table_dir, 'stage_c_summary_by_roi_source_count.csv'));
W = sortrows(summarize_metrics(T, ["model_name","condition_id","geometry","f0","field_regime","N_sources"]), 'MAPE_pct', 'descend');
W = W(1:min(40,height(W)),:);
writetable(W, fullfile(OUT.table_dir, 'stage_c_worst_conditions.csv'));
end

function S = summarize_metrics(T, groups)
groups = string(groups);
[G, keys] = findgroups(T(:, cellstr(groups)));
S = keys;
S.N = splitapply(@numel, T.abs_error_pct, G);
S.MAPE_pct = splitapply(@(x) mean(x,'omitnan'), T.abs_error_pct, G);
S.signed_bias_pct = splitapply(@(x) mean(x,'omitnan'), T.signed_error_pct, G);
S.high_error10_pct = splitapply(@(x) 100*mean(x,'omitnan'), double(T.high_error10), G);
S.high_error20_pct = splitapply(@(x) 100*mean(x,'omitnan'), double(T.high_error20), G);
S.median_abs_error_pct = splitapply(@(x) median(x,'omitnan'), T.abs_error_pct, G);
S.mean_predicted_SWS = splitapply(@(x) mean(x,'omitnan'), T.SWS_pred, G);
S.mean_true_SWS = splitapply(@(x) mean(x,'omitnan'), T.SWS_true, G);
S.mean_q_pred = splitapply(@(x) mean(x,'omitnan'), T.q_pred, G);
end

function row = condition_log_row(C, status, nRows, elapsed, msg)
row = table();
row.condition_id = string(C.condition_id);
row.geometry = string(C.geometry);
row.field_regime = string(C.field_regime);
row.f0 = C.f0;
row.M = C.M;
row.status = string(status);
row.patch_model_rows = nRows;
row.elapsed_sec = elapsed;
row.message = string(msg);
end

function tf = should_save_representative_map(C, CFG)
tf = CFG.Figures.SaveRepresentativeMaps && ...
    ismember(string(C.geometry), CFG.Figures.RepresentativeGeometries) && ...
    ismember(C.f0, CFG.Figures.RepresentativeFrequenciesHz) && ...
    ismember(string(C.field_regime), CFG.Figures.RepresentativeFieldRegimes);
end

function plot_representative_maps(sim, T, C, OUT)
Tc = T(string(T.model_name)=="q_spectrum_plus_composition",:);
swsMap = scatter_to_map(Tc, Tc.SWS_pred, size(sim.cs_map));
qMap = scatter_to_map(Tc, Tc.q_pred, size(sim.cs_map));
signErr = scatter_to_map(Tc, Tc.signed_error_pct, size(sim.cs_map));
absErr = scatter_to_map(Tc, Tc.abs_error_pct, size(sim.cs_map));
qGrad = gradient_magnitude(qMap);
fig = figure('Color','w','Position',[60 60 1500 820], 'Visible','off');
tl = tiledlayout(fig, 2, 5, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s | f=%d Hz | %s', pretty_id(C.geometry), C.f0, pretty_id(C.field_regime)), 'Interpreter','none', 'FontSize', 12);
plot_map(nexttile(tl), sim.cs_map, 'True SWS', 'm/s');
plot_map(nexttile(tl), abs(sim.U), 'Measured amplitude', 'a.u.');
plot_map(nexttile(tl), sign(cos(angle(sim.U))), 'sign(cos phase)', 'sign');
plot_map(nexttile(tl), real(sim.U), 'Real displacement', 'a.u.');
plot_map(nexttile(tl), qMap, 'Predicted q', 'REQ quantile');
plot_map(nexttile(tl), swsMap, 'Predicted SWS', 'm/s');
plot_map(nexttile(tl), signErr, 'Signed error', '%');
plot_map(nexttile(tl), absErr, 'Absolute error', '%');
plot_map(nexttile(tl), qGrad, '|grad q|', 'a.u.');
plot_source_projection(nexttile(tl), sim, C);
folder = fullfile(OUT.map_dir, char(C.geometry), char(C.field_regime)); ensure_dir(folder);
fn = fullfile(folder, "stageC_map__" + sanitize(C.condition_id) + ".png");
exportgraphics(fig, fn, 'Resolution', 180); close(fig);
end

function plot_source_geometry(sim, Tsrc, C, OUT)
folder = fullfile(OUT.source_dir, char(C.field_regime)); ensure_dir(folder);
fig = figure('Color','w','Position',[80 80 1350 480], 'Visible','off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Source geometry: %s', pretty_id(C.field_regime)), 'Interpreter','none', 'FontSize', 11);
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); view(ax, 3);
plot_box3(ax, max(sim.x), max(sim.z));
scatter3(ax, Tsrc.x_mm, Tsrc.y_mm, Tsrc.z_mm, 50, Tsrc.is_in_plane, 'filled');
quiver3(ax, Tsrc.x_mm, Tsrc.y_mm, Tsrc.z_mm, 0*Tsrc.x_mm, 0*Tsrc.x_mm, 6+0*Tsrc.x_mm, 'k','LineWidth',1.2);
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); zlabel(ax,'z (mm)'); title(ax,'3D positions');
axis(ax,'equal');
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
rectangle(ax,'Position',[0 -10 max(sim.x)*1e3 20],'EdgeColor',[0 0 0],'LineWidth',1.1);
yline(ax,0,'b-','Measurement plane y=0');
scatter(ax,Tsrc.x_mm,Tsrc.y_mm,50,Tsrc.is_in_plane,'filled');
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); title(ax,'x-y view'); axis(ax,'equal');
ax = nexttile(tl); imagesc(ax, sim.x*1e3, sim.z*1e3, sim.cs_map); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on');
scatter(ax,Tsrc.x_mm,Tsrc.z_mm,45,Tsrc.is_in_plane,'filled','MarkerEdgeColor','k');
quiver(ax,Tsrc.x_mm,Tsrc.z_mm,0*Tsrc.x_mm,5+0*Tsrc.x_mm,'k','LineWidth',1.1);
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)'); title(ax,'x-z projection'); cb=colorbar(ax); ylabel(cb,'SWS (m/s)');
fn = fullfile(folder, "stageC_sources__" + sanitize(C.field_regime) + ".png");
exportgraphics(fig, fn, 'Resolution', 180);
savefig(fig, fullfile(folder, "stageC_sources__" + sanitize(C.field_regime) + ".fig"));
close(fig);
end

function plot_box3(ax, Lx, Lz)
x = [0 Lx Lx 0 0]*1e3; z = [0 0 Lz Lz 0]*1e3;
for y = [-10 10]
    plot3(ax,x,y+zeros(size(x)),z,'k-','LineWidth',1);
end
for k = 1:4
    plot3(ax,[x(k) x(k)],[ -10 10],[z(k) z(k)],'k-','LineWidth',1);
end
plot3(ax,[0 max(x) max(x) 0 0], zeros(1,5), [0 0 max(z) max(z) 0], 'b-', 'LineWidth',1.2);
end

function plot_source_projection(ax, sim, C)
imagesc(ax, sim.x*1e3, sim.z*1e3, sim.cs_map); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on');
P = sim.source_positions_xyz;
scatter(ax, P(:,1)*1e3, P(:,3)*1e3, 35, abs(P(:,2))<1e-12, 'filled','MarkerEdgeColor','k');
quiver(ax, P(:,1)*1e3, P(:,3)*1e3, zeros(size(P,1),1), 5*ones(size(P,1),1), 'k','LineWidth',1);
title(ax,'Source projection'); xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)'); cb=colorbar(ax); ylabel(cb,'m/s');
end

function M = scatter_to_map(T, vals, sz)
M = nan(sz);
if ismember('cz', T.Properties.VariableNames) && ismember('cx', T.Properties.VariableNames)
    idx = sub2ind(sz, T.cz, T.cx);
elseif ismember('map_iz', T.Properties.VariableNames) && ismember('map_ix', T.Properties.VariableNames)
    idx = sub2ind(sz, T.map_iz, T.map_ix);
else
    return;
end
M(idx) = vals;
end

function G = gradient_magnitude(A)
A2 = A;
if all(isnan(A2(:))), G = A2; return; end
A2 = fillmissing(A2,'nearest',1); A2 = fillmissing(A2,'nearest',2);
[gx,gz] = gradient(A2);
G = hypot(gx,gz);
G(isnan(A)) = nan;
end

function plot_map(ax, M, ttl, cblabel)
imagesc(ax, M); axis(ax,'image'); set(ax,'YDir','normal'); title(ax, ttl, 'FontSize',9);
cb = colorbar(ax); ylabel(cb, cblabel, 'FontSize',8);
end

function MAPS = compact_maps_for_cache(sim, T, C, CFG)
MAPS = struct();
MAPS.cs_map = sim.cs_map;
MAPS.label_map = sim.label_map;
MAPS.amplitude = sim.amplitude;
MAPS.phase = sim.phase;
MAPS.x = sim.x; MAPS.z = sim.z;
Tmain = T(string(T.model_name)=="q_spectrum_plus_composition",:);
MAPS.sws_pred = scatter_to_map(Tmain, Tmain.SWS_pred, size(sim.cs_map));
MAPS.q_pred = scatter_to_map(Tmain, Tmain.q_pred, size(sim.cs_map));
MAPS.signed_error_pct = scatter_to_map(Tmain, Tmain.signed_error_pct, size(sim.cs_map));
MAPS.abs_error_pct = scatter_to_map(Tmain, Tmain.abs_error_pct, size(sim.cs_map));
end

function y = clamp01(x)
y = min(max(x,0),1);
end

function s = sanitize(x)
s = regexprep(char(string(x)), '[^A-Za-z0-9]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
end

function ensure_dir(d)
if exist(d,'dir') ~= 7, mkdir(d); end
end

function tf = env_true(name, defaultVal)
val = string(getenv(name));
if val == "", tf = defaultVal; return; end
tf = any(lower(strtrim(val)) == ["1","true","yes","on"]);
end

function write_json_copy(CFG, path)
fid = fopen(path,'w'); assert(fid>0);
fprintf(fid, '%s', jsonencode(CFG, PrettyPrint=true));
fclose(fid);
end

function name = pretty_id(x)
name = strrep(char(string(x)), '_', ' ');
end

function print_summary_block(T, label)
fprintf('\n%s\n', label);
disp(T(1:min(12,height(T)),:));
end

function write_stage_readme(OUT, CFG, T_log, elapsedSec)
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w'); assert(fid>0);
fprintf(fid, '# Stage C: Field Complexity Transfer Validation\n\n');
fprintf(fid, '## Objective\n\nEvaluate whether frozen `baseline_minimal_v1` q/SWS models transfer to clean Eikonal fields as source count and off-plane field complexity increase. No retraining, correction, risk mask, readout noise, or shear-amplitude SNR coupling is used.\n\n');
fprintf(fid, '## Design\n\n');
fprintf(fid, '- Geometries: %s.\n', strjoin(string({CFG.Geometries.id}), ', '));
fprintf(fid, '- Frequencies: %s Hz.\n', strjoin(string(CFG.FrequenciesHz), ', '));
fprintf(fid, '- Field regimes: %s.\n', strjoin(CFG.FieldRegimes, ', '));
fprintf(fid, '- M=%s, cs_guess=%.1f m/s, TargetStepM=%.1f mm.\n', mat2str(CFG.MList), CFG.REQ.cs_guess, CFG.REQ.TargetStepM*1e3);
fprintf(fid, '- All source motion vectors are `[0 0 1]`, parallel to the measured displacement component.\n');
fprintf(fid, '- The final summed complex field is globally normalized by central-region RMS. This isolates field complexity from global amplitude/SNR changes while preserving local interference.\n\n');
fprintf(fid, '## Source Layouts\n\n');
fprintf(fid, '- `single_source_lateral`: one in-plane lateral source.\n');
fprintf(fid, '- Layout A: one in-plane source for N=4, 8, 16; the rest are off-plane shell sources.\n');
fprintf(fid, '- Layout B: practical multi-actuator-like source counts with 1, 2, and 4 in-plane sources for N=4, 8, and 16 respectively.\n\n');
fprintf(fid, '## Outputs\n\nTables are under `tables/`; figures are under `figures/`. Representative maps include true SWS, amplitude, sign phase, real displacement, predicted SWS, signed error, absolute error, q, q-gradient, and source projection.\n\n');
fprintf(fid, '## Runtime\n\nObserved runner time for this run: %.1f minutes. Cached conditions are reused unless `ADAPTIVE_REQ_EIKONAL_STAGE_C_FORCE_REBUILD=true`. Source-level parallelism used: %d. REQ-window parallelism used: %d.\n\n', elapsedSec/60, CFG.UseSourceParfor, CFG.UseWindowParfor);
fprintf(fid, 'Completed conditions: %d. Failed conditions: %d.\n\n', sum(ismember(string(T_log.status), ["completed","reused"])), sum(string(T_log.status)=="failed"));
fprintf(fid, '## Commands\n\n');
fprintf(fid, 'Validation-only:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE'',''validate_only''); run(''experiments/runners/eikonal_validation/run_stage_c_field_complexity.m'')"\n```\n\n');
fprintf(fid, 'Quick:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE'',''quick''); run(''experiments/runners/eikonal_validation/run_stage_c_field_complexity.m'')"\n```\n\n');
fprintf(fid, 'Full and analysis:\n\n```bash\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE'',''full''); run(''experiments/runners/eikonal_validation/run_stage_c_field_complexity.m'')"\nmatlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE'',''full''); run(''experiments/analysis/eikonal_validation/analyze_stage_c_field_complexity.m'')"\n```\n');
fclose(fid);
end

function write_docs_readme(root_dir, CFG)
docDir = fullfile(root_dir, 'docs', 'eikonal_validation'); ensure_dir(docDir);
fid = fopen(fullfile(docDir, 'stage_C_field_complexity.md'), 'w'); assert(fid>0);
fprintf(fid, '# Stage C: Field Complexity Transfer Validation\n\n');
fprintf(fid, 'Stage C tests whether clean Eikonal field complexity alone changes zero-shot transfer of frozen `baseline_minimal_v1` q/SWS models. It uses no noise, no correction, and no estimability mask.\n\n');
fprintf(fid, 'Outputs are written to `%s`.\n', CFG.OutputRoot);
fclose(fid);
end
