%% run_stage_B_resolution_phantom.m
% Stage B: resolution phantom / geometry transfer validation.
%
% This runner evaluates frozen baseline_minimal_v1 q/SWS models on clean
% Eikonal multi-inclusion phantoms. It does not train, correct, or apply
% reliability masks. Truth maps, ROI labels, diameter, and interface distance
% are used only after prediction for evaluation.
%
% Runtime controls:
%   ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE = validate_only | quick | full
%   ADAPTIVE_REQ_EIKONAL_STAGE_B_FORCE_REBUILD = true | false
%   ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_SOURCE_PARFOR = true | false
%   ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_WINDOW_PARFOR = true | false

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
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
assert(ismember(CFG.Mode, ["validate_only","quick","full"]), ...
    'ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE must be validate_only, quick, or full.');
CFG.ForceRebuild = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_B_FORCE_REBUILD', false);
CFG.UseSourceParfor = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_SOURCE_PARFOR', true);
CFG.UseWindowParfor = env_true('ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_WINDOW_PARFOR', false);
CFG.UseParallel = CFG.UseSourceParfor || CFG.UseWindowParfor;
CFG.CacheVersion = 3;

OUT = make_output_dirs(root_dir, CFG);
write_json_copy(CFG, fullfile(OUT.root_dir, 'stage_B_resolution_phantom_effective_config.json'));
waveSrc = fullfile(char(CFG.WaveSimProjectRoot), 'src');
addpath(waveSrc);
rehash;
assert(exist(fullfile(waveSrc, '+simcore', '+config', 'defaultConfig.m'), 'file') == 2, ...
    'Could not find simcore. Check WaveSimProjectRoot in config.');

fprintf('\nStage B: Resolution Phantom / Geometry Transfer Validation\n');
fprintf('Mode: %s | force rebuild: %d\n', CFG.Mode, CFG.ForceRebuild);
fprintf('Output root: %s\n', OUT.root_dir);
fprintf('Wave sim root: %s\n', CFG.WaveSimProjectRoot);
fprintf('Domain: %.1f x %.1f mm | dx=dz=%.1f mm\n', ...
    CFG.Simulation.Lx_m*1e3, CFG.Simulation.Lz_m*1e3, CFG.Simulation.dx_m*1e3);
fprintf('Target REQ step: %.3f mm | M=%s | cs_guess=%.2f m/s\n', ...
    CFG.REQ.TargetStepM*1e3, mat2str(CFG.MList), CFG.REQ.cs_guess);
fprintf('Parallel: source parfor=%d | REQ window parfor=%d\n', CFG.UseSourceParfor, CFG.UseWindowParfor);
if CFG.UseParallel
    ensure_parallel_pool();
end

BUNDLE = load_baseline_bundle(root_dir, CFG);
assert_no_forbidden_predictors(BUNDLE.BASE_FEATURES);
fprintf('Loaded frozen models: %s\n', strjoin(string(BUNDLE.MODELS.model_names), ', '));
fprintf('Required base predictors: %d\n', numel(BUNDLE.BASE_FEATURES));

CONDS = build_condition_list(CFG);
fprintf('Stage B conditions selected: %d\n', numel(CONDS));

%% Run conditions
patchParts = {};
sourceParts = {};
inclusionParts = {};
logParts = {};
mapCount = 0;
spectrumDone = containers.Map('KeyType','char','ValueType','logical');
totalTimer = tic;

for ci = 1:numel(CONDS)
    C = CONDS(ci);
    fprintf('\n[%d/%d] %s\n', ci, numel(CONDS), C.condition_id);
    cacheFile = fullfile(OUT.cache_dir, "stageB__" + sanitize(C.condition_id) + ".mat");
    conditionTimer = tic;
    try
        sig = condition_signature(CFG, C);
        if ~CFG.ForceRebuild && exist(cacheFile, 'file') == 2
            S = load(cacheFile);
            if isfield(S, 'cond_signature') && isequaln(S.cond_signature, sig)
                T_condition = S.T_condition;
                T_sources = S.T_sources;
                T_inclusions = S.T_inclusions;
                fprintf('  reused cache: %d patch/model rows.\n', height(T_condition));
                elapsed = toc(conditionTimer);
                patchParts{end+1,1} = T_condition; %#ok<SAGROW>
                sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
                inclusionParts{end+1,1} = T_inclusions; %#ok<SAGROW>
                logParts{end+1,1} = condition_log_row(C, "reused", height(T_condition), elapsed, ""); %#ok<SAGROW>
                continue;
            end
        end

        [sim, Tsrc, Tinc, simDiag] = run_eikonal_condition(C, CFG); %#ok<ASGLU>
        quality_check_field(sim, C, CFG);
        [F, reqDiag] = extract_req_feature_table(sim, C, CFG, BUNDLE.BASE_FEATURES); %#ok<ASGLU>
        T_condition = apply_frozen_models(F, BUNDLE, C, CFG);
        T_sources = Tsrc;
        T_sources.condition_id = repmat(string(C.condition_id), height(T_sources), 1);
        T_sources = movevars(T_sources, 'condition_id', 'Before', 1);
        T_inclusions = Tinc;
        T_inclusions.condition_id = repmat(string(C.condition_id), height(T_inclusions), 1);
        T_inclusions = movevars(T_inclusions, 'condition_id', 'Before', 1);

        cond_signature = sig; %#ok<NASGU>
        MAPS = compact_maps_for_cache(sim, T_condition, C, CFG); %#ok<NASGU>
        save(cacheFile, 'T_condition', 'T_sources', 'T_inclusions', 'cond_signature', 'MAPS', 'simDiag', 'reqDiag', '-v7.3');

        if should_save_representative_map(C, CFG)
            plot_representative_maps(sim, T_condition, C, OUT);
            mapCount = mapCount + 1;
        end
        if CFG.Figures.SaveSourceGeometry
            plot_source_geometry(sim, Tsrc, C, OUT);
        end
        if CFG.Figures.SaveCentralSpectrum && ~isKey(spectrumDone, char(C.field_regime))
            plot_central_power_spectrum(sim, C, OUT);
            spectrumDone(char(C.field_regime)) = true;
        end

        elapsed = toc(conditionTimer);
        fprintf('  completed: %d centers x %d models in %.1f s.\n', ...
            height(F), numel(BUNDLE.MODELS.model_names), elapsed);
        patchParts{end+1,1} = T_condition; %#ok<SAGROW>
        sourceParts{end+1,1} = T_sources; %#ok<SAGROW>
        inclusionParts{end+1,1} = T_inclusions; %#ok<SAGROW>
        logParts{end+1,1} = condition_log_row(C, "completed", height(T_condition), elapsed, ""); %#ok<SAGROW>
    catch ME
        elapsed = toc(conditionTimer);
        warning('StageB:ConditionFailed', 'Condition failed (%s): %s', C.condition_id, ME.message);
        logParts{end+1,1} = condition_log_row(C, "failed", 0, elapsed, string(ME.message)); %#ok<SAGROW>
    end
end

if isempty(patchParts)
    error('No Stage B conditions completed successfully.');
end

T_patch = vertcat(patchParts{:});
T_source = vertcat(sourceParts{:});
T_inclusion = vertcat(inclusionParts{:});
T_log = vertcat(logParts{:});

%% Tables and summaries
writetable(T_patch, fullfile(OUT.table_dir, 'stage_B_patch_level_results.csv'));
writetable(T_source, fullfile(OUT.table_dir, 'stage_B_source_table.csv'));
writetable(T_inclusion, fullfile(OUT.table_dir, 'stage_B_inclusion_table.csv'));
writetable(T_log, fullfile(OUT.table_dir, 'stage_B_condition_log.csv'));

write_all_summaries(T_patch, OUT);
save(fullfile(OUT.data_dir, 'stage_B_resolution_phantom_results.mat'), ...
    'T_patch','T_source','T_inclusion','T_log','CFG','-v7.3');

write_stage_readme(OUT, CFG, T_log, toc(totalTimer));
write_docs_readme(root_dir, CFG);

%% Console summary
nFailed = sum(string(T_log.status) == "failed");
fprintf('\nStage B runner complete.\n');
fprintf('Completed conditions: %d | failed: %d | representative maps: %d\n', ...
    sum(ismember(string(T_log.status), ["completed","reused"])), nFailed, mapCount);
S_overall = summarize_metrics(T_patch, "model_name");
S_diam = summarize_metrics(T_patch, ["model_name","inclusion_sws","inclusion_diameter_mm"]);
S_field = summarize_metrics(T_patch, ["model_name","field_regime"]);
S_worst = sortrows(summarize_metrics(T_patch, ["model_name","condition_id","phantom_name","f0","field_regime"]), 'MAPE_pct', 'descend');
print_summary_block(S_overall, 'Global summary');
print_summary_block(S_diam, 'By inclusion diameter');
print_summary_block(S_field, 'By field regime');
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
        try
            parpool('threads');
        catch
            parpool;
        end
    end
catch ME
    warning('StageB:ParallelPool', 'Could not start a parallel pool. Continuing serially: %s', ME.message);
end
end

function CFG = load_stage_config(root_dir)
configPath = fullfile(root_dir, 'configs', 'eikonal_validation', 'stage_B_resolution_phantom.json');
assert(exist(configPath, 'file') == 2, 'Missing Stage B config: %s', configPath);
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
CFG.Phantoms = CFG.Phantoms(:).';
for i = 1:numel(CFG.Phantoms)
    CFG.Phantoms(i).id = string(CFG.Phantoms(i).id);
    CFG.Phantoms(i).diameters_mm = double(CFG.Phantoms(i).diameters_mm(:)).';
end
CFG.Figures.RepresentativePhantoms = string(CFG.Figures.RepresentativePhantoms(:)).';
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
    'Bundle does not contain both Stage B models.');
end

function assert_no_forbidden_predictors(features)
features = lower(string(features(:)));
forbiddenExact = ["true_sws","sws_true","cs_true","k_true","q_oracle", ...
    "q_true","q_pred","sws_pred","signed_error_percent","abs_error_percent", ...
    "patch_purity","true_patch_purity","material_label","label_map", ...
    "distance_to_interface","distance_to_interface_m","roi","roi_label", ...
    "confidence","risk","req_mapping","map_ix","map_iz","cx","cz", ...
    "x_center_m","z_center_m","condition_id","inclusion_diameter_mm", ...
    "inclusion_id","d_over_lambda","d_over_lwin","m_eff"];
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
phantomIds = string(modeBlock.Phantoms(:)).';
freqs = double(modeBlock.FrequenciesHz(:)).';
regimes = string(modeBlock.FieldRegimes(:)).';
Ms = double(modeBlock.MList(:)).';
CONDS = struct('phantom_name',{},'background_sws',{},'inclusion_sws',{},'diameters_mm',{}, ...
    'f0',{},'M',{},'field_regime',{},'realism_level',{},'source_layout_id',{}, ...
    'source_seed',{},'noise_seed',{},'condition_id',{});
idx = 0;
for pi = 1:numel(phantomIds)
    P = find_phantom(CFG, phantomIds(pi));
    for fi = 1:numel(freqs)
        for ri = 1:numel(regimes)
            for mi = 1:numel(Ms)
                idx = idx + 1;
                C = struct();
                C.phantom_name = P.id;
                C.background_sws = double(P.background_sws);
                C.inclusion_sws = double(P.inclusion_sws);
                C.diameters_mm = double(P.diameters_mm(:)).';
                C.f0 = freqs(fi);
                C.M = Ms(mi);
                C.field_regime = regimes(ri);
                C.realism_level = "clean";
                C.source_layout_id = regimes(ri);
                C.source_seed = source_seed_for(CFG, regimes(ri));
                C.noise_seed = 0;
                C.condition_id = sprintf('%s__f%d__clean__%s__M%d', ...
                    C.phantom_name, C.f0, C.field_regime, C.M);
                CONDS(idx) = C; %#ok<AGROW>
            end
        end
    end
end
end

function P = find_phantom(CFG, id)
ids = string({CFG.Phantoms.id});
idx = find(ids == string(id), 1);
assert(~isempty(idx), 'Unknown Stage B phantom: %s', id);
P = CFG.Phantoms(idx);
end

function seed = source_seed_for(CFG, regime)
S = CFG.SourceLayouts.(char(regime));
seed = double(S.source_seed);
end

function sig = condition_signature(CFG, C)
sig = struct();
sig.CacheVersion = CFG.CacheVersion;
sig.phantom_name = C.phantom_name;
sig.background_sws = C.background_sws;
sig.inclusion_sws = C.inclusion_sws;
sig.diameters_mm = C.diameters_mm;
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

function [sim, Tsrc, Tinc, D] = run_eikonal_condition(C, CFG)
phant = build_multi_inclusion_resolution_phantom(C, CFG);
S = build_source_layout(C, CFG);
validate_source_layout(S, C);
U_pre = synthesize_clean_eikonal_field(phant, S, C, CFG);
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
            error('Unsupported Stage B amplitude normalization: %s', CFG.Simulation.AmplitudeNormalization);
    end
end
U = U_pre .* scale;
postRms = sqrt(mean(abs(U(:)).^2, 'omitnan'));
postMedian = median(abs(U(:)), 'omitnan');

sim = struct();
sim.Uxz = U;
sim.U = U;
sim.real_U = real(U);
sim.amplitude = abs(U);
sim.phase = angle(U);
sim.cs_map = phant.cs_map;
sim.label_map = phant.label_map;
sim.inclusion_id_map = phant.inclusion_id_map;
sim.signed_distance_m = phant.signed_distance_m;
sim.abs_distance_m = phant.abs_distance_m;
sim.nearest_diameter_mm_map = phant.nearest_diameter_mm_map;
sim.nearest_inclusion_id_map = phant.nearest_inclusion_id_map;
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
sim.amplitude_normalization_scale = scale;
sim.pre_norm_rms = preRms;
sim.pre_norm_median = preMedian;
sim.post_norm_rms = postRms;
sim.post_norm_median = postMedian;
sim.tracking_snr_db = nan(size(sim.amplitude));
Tsrc = source_table(S, C, CFG);
Tinc = phant.inclusion_table;
D = struct('pre_rms',preRms,'pre_median',preMedian,'post_rms',postRms, ...
    'post_median',postMedian,'normalization_scale',scale);
end

function phant = build_multi_inclusion_resolution_phantom(C, CFG)
x = 0:CFG.Simulation.dx_m:CFG.Simulation.Lx_m;
z = 0:CFG.Simulation.dz_m:CFG.Simulation.Lz_m;
[X,Z] = meshgrid(x,z);
cs = C.background_sws * ones(size(X));
label = zeros(size(X));
incIdMap = zeros(size(X));
Dmm = double(C.diameters_mm(:));
radii = Dmm(:)'/2 * 1e-3;
Lx = CFG.Simulation.Lx_m; Lz = CFG.Simulation.Lz_m;
centers = [0.037 0.037; 0.099 0.037; 0.037 0.099; 0.099 0.099];
assert(numel(radii) == 4, 'Stage B expects four inclusion diameters.');
for ii = 1:4
    r = hypot(X-centers(ii,1), Z-centers(ii,2));
    mask = r <= radii(ii);
    cs(mask) = C.inclusion_sws;
    label(mask) = ii;
    incIdMap(mask) = ii;
end
% Distance to nearest inclusion boundary. Positive inside any inclusion,
% negative in background. The nearest inclusion metadata is used for D/lambda.
signedAll = nan([size(X), 4]);
for ii = 1:4
    signedAll(:,:,ii) = radii(ii) - hypot(X-centers(ii,1), Z-centers(ii,2));
end
absAll = abs(signedAll);
[absD, nearestId] = min(absAll, [], 3);
signedD = -absD;
inside = label > 0;
signedD(inside) = absD(inside);
nearestDmm = Dmm(nearestId);
% For pixels inside an inclusion, force nearest ID to the actual inclusion.
nearestId(inside) = label(inside);
nearestDmm(inside) = Dmm(label(inside));

Tinc = table();
Tinc.inclusion_id = (1:4)';
Tinc.inclusion_diameter_mm = Dmm(:);
Tinc.inclusion_radius_mm = Dmm(:)/2;
Tinc.center_x_m = centers(:,1);
Tinc.center_z_m = centers(:,2);
Tinc.center_x_mm = centers(:,1)*1e3;
Tinc.center_z_mm = centers(:,2)*1e3;
Tinc.inclusion_sws = C.inclusion_sws * ones(4,1);
Tinc.background_sws = C.background_sws * ones(4,1);
% Layout QA.
for ii = 1:4
    border = min([centers(ii,1)-radii(ii), Lx-centers(ii,1)-radii(ii), centers(ii,2)-radii(ii), Lz-centers(ii,2)-radii(ii)]);
    if border < 0.020 - 1e-12
        warning('StageB:BorderMargin', 'Inclusion %d border margin is %.1f mm (<20 mm).', ii, border*1e3);
    end
end
for ii = 1:4
    for jj = ii+1:4
        edgeSep = hypot(centers(ii,1)-centers(jj,1), centers(ii,2)-centers(jj,2)) - radii(ii) - radii(jj);
        if edgeSep < 0.025 - 1e-12
            warning('StageB:InclusionSpacing', 'Inclusion %d-%d edge separation is %.1f mm (<25 mm).', ii, jj, edgeSep*1e3);
        end
    end
end
phant = struct('x',x,'z',z,'X',X,'Z',Z,'cs_map',cs,'label_map',label, ...
    'inclusion_id_map',incIdMap,'signed_distance_m',signedD, ...
    'abs_distance_m',absD,'nearest_inclusion_id_map',nearestId, ...
    'nearest_diameter_mm_map',nearestDmm,'inclusion_table',Tinc);
end

function U = synthesize_clean_eikonal_field(phant, S, C, CFG)
omega = 2*pi*C.f0;
Nsrc = size(S.positions_xyz,1);
usePar = CFG.UseSourceParfor && Nsrc > 1;
parts = cell(Nsrc,1);
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
opts.MaxIterations = 200;
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
Lx = CFG.Simulation.Lx_m;
Lz = CFG.Simulation.Lz_m;
cx = Lx/2; cz = Lz/2;
yOff = 24e-3;
regime = string(C.field_regime);
switch regime
    case "single_source_lateral"
        P = [-12e-3 0 cz];
        phaseSeed = 0;
        phases = 0;
    case "diffuse_like_8src_layout1"
        P = [
            -12e-3    0       cz
             Lx+12e-3 0       cz
             cx       0      -12e-3
             cx       0       Lz+12e-3
            -10e-3   -yOff   -10e-3
             Lx+10e-3 yOff    Lz+10e-3
            -10e-3    yOff    Lz+10e-3
             Lx+10e-3 -yOff  -10e-3
        ];
        phaseSeed = 6601;
        phases = deterministic_phases(size(P,1), phaseSeed);
    case "diffuse_like_8src_layout2"
        P = [
            -12e-3    0       0.25*Lz
             Lx+12e-3 0       0.75*Lz
             0.25*Lx  0      -12e-3
             0.75*Lx  0       Lz+12e-3
            -10e-3    yOff    0.80*Lz
             Lx+10e-3 -yOff   0.20*Lz
             0.15*Lx -yOff    Lz+10e-3
             0.85*Lx  yOff   -10e-3
        ];
        phaseSeed = 6602;
        phases = deterministic_phases(size(P,1), phaseSeed);
    otherwise
        error('Unsupported Stage B field regime: %s', regime);
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
    'All Stage B source motion vectors must be [0 0 1].');
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
assert(numel(unique(sim.inclusion_id_map(sim.inclusion_id_map>0))) == 4, 'Expected four inclusions.');
amp = sim.amplitude(:);
valid = isfinite(amp);
assert(any(valid), 'No finite amplitude values.');
lowFrac = mean(amp(valid) < 1e-8);
if lowFrac > 0.20
    warning('StageB:LowAmplitude', '%s has %.1f%% near-zero amplitude pixels.', C.condition_id, 100*lowFrac);
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
F.phantom_name = repmat(string(C.phantom_name), height(F), 1);
F.geometry = F.phantom_name;
F.background_sws = C.background_sws * ones(height(F),1);
F.inclusion_sws = C.inclusion_sws * ones(height(F),1);
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
F.SWS_true = sample_map_at_centers(sim.cs_map, F);
F.true_sws = F.SWS_true;
F.inclusion_id = sample_map_at_centers(sim.inclusion_id_map, F);
F.inclusion_id = round(F.inclusion_id);
F.nearest_inclusion_id = round(sample_map_at_centers(sim.nearest_inclusion_id_map, F));
F.inclusion_diameter_mm = sample_map_at_centers(sim.nearest_diameter_mm_map, F);
F.distance_to_interface_m = sample_map_at_centers(sim.abs_distance_m, F);
winRadiusM = 0.5 * C.M * CFG.REQ.cs_guess / C.f0;
F.distance_to_interface_over_window_radius = F.distance_to_interface_m ./ max(winRadiusM, eps);
F.patch_purity = compute_patch_purity(sim.label_map, F, C, CFG);
[F.patch_cs_std, F.patch_cs_range, F.patch_cs_iqr] = compute_patch_cs_stats(sim.cs_map, F, C, CFG);
F.D_over_lambda = (F.inclusion_diameter_mm*1e-3) ./ (C.inclusion_sws / C.f0);
F.D_over_Lwin = (F.inclusion_diameter_mm*1e-3) ./ (C.M * CFG.REQ.cs_guess / C.f0);
F.M_eff = C.M * CFG.REQ.cs_guess ./ F.SWS_true;
F.roi_label = assign_roi_labels(F, CFG);
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

function roi = assign_roi_labels(F, CFG)
roi = repmat("background_other", height(F), 1);
dmm = F.distance_to_interface_m * 1e3;
inside = F.SWS_true > F.background_sws + 0.1;
farBg = ~inside & dmm >= CFG.ROI.BackgroundFarMinDistanceMm;
roi(farBg) = "background_far";
for bi = 1:size(CFG.ROI.InterfaceBandsMm,1)
    lo = CFG.ROI.InterfaceBandsMm(bi,1); hi = CFG.ROI.InterfaceBandsMm(bi,2);
    mask = dmm >= lo & dmm < hi;
    roi(mask) = sprintf('interface_%g_%gmm', lo, hi);
end
% Core has priority over broad interface bands. This keeps small-but-resolvable
% inclusions visible in core summaries when they have clean central patches.
core = inside & dmm >= CFG.ROI.CoreMinDistanceMm & F.patch_purity >= 0.95;
roi(core) = "inclusion_core";
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
    error('Required frozen baseline predictor missing from Stage B feature table: %s', strjoin(missing, ', '));
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
            error('Unsupported Stage B model: %s', name);
    end
    if any(~isfinite(q)), error('Non-finite q predictions for %s.', name); end
    if any(q < -0.05 | q > 1.05)
        warning('StageB:QRange', '%s produced q outside [-0.05, 1.05]; values are clamped for SWS conversion.', name);
    end
    q = clamp01(q);
    sws = q_to_sws(Faug.req_mapping, q, Faug.f0);
    if any(~isfinite(sws)), error('Non-finite SWS predictions for %s.', name); end
    R = table();
    R.condition_id = F.condition_id;
    R.phantom_name = F.phantom_name;
    R.geometry = F.geometry;
    R.background_sws = F.background_sws;
    R.inclusion_sws = F.inclusion_sws;
    R.inclusion_id = F.inclusion_id;
    R.nearest_inclusion_id = F.nearest_inclusion_id;
    R.inclusion_diameter_mm = F.inclusion_diameter_mm;
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
    R.patch_purity = F.patch_purity;
    R.patch_cs_std = F.patch_cs_std;
    R.patch_cs_range = F.patch_cs_range;
    R.patch_cs_iqr = F.patch_cs_iqr;
    R.distance_to_interface_m = F.distance_to_interface_m;
    R.distance_to_interface_over_window_radius = F.distance_to_interface_over_window_radius;
    R.D_over_lambda = F.D_over_lambda;
    R.D_over_Lwin = F.D_over_Lwin;
    R.M_eff = F.M_eff;
    R.roi_label = F.roi_label;
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

function write_all_summaries(T_patch, OUT)
T_overall = summarize_metrics(T_patch, "model_name");
T_by_model = T_overall;
T_by_phantom = summarize_metrics(T_patch, ["model_name","phantom_name"]);
T_by_diam = summarize_metrics(T_patch, ["model_name","inclusion_sws","inclusion_diameter_mm"]);
T_by_sws = summarize_metrics(T_patch, ["model_name","inclusion_sws"]);
T_by_freq = summarize_metrics(T_patch, ["model_name","f0"]);
T_by_field = summarize_metrics(T_patch, ["model_name","field_regime"]);
T_by_Dlambda = summarize_metrics(add_bin(T_patch,"D_over_lambda","D_over_lambda_bin",[0 1 2 3 4 6 8 12 inf]), ["model_name","D_over_lambda_bin"]);
T_by_DLwin = summarize_metrics(add_bin(T_patch,"D_over_Lwin","D_over_Lwin_bin",[0 0.5 1 1.5 2 3 4 6 inf]), ["model_name","D_over_Lwin_bin"]);
T_by_Meff = summarize_metrics(add_bin(T_patch,"M_eff","M_eff_bin",[0 1.5 2 2.5 3 4 inf]), ["model_name","M_eff_bin"]);
T_by_dist = summarize_metrics(add_bin(T_patch,"distance_to_interface_over_window_radius","distance_over_window_bin",[0 0.25 0.5 1 2 inf]), ["model_name","distance_over_window_bin"]);
T_by_roi = summarize_metrics(T_patch, ["model_name","roi_label"]);
T_by_roi_freq = summarize_metrics(T_patch, ["model_name","roi_label","f0"]);
T_by_diam_freq = summarize_metrics(T_patch, ["model_name","inclusion_sws","inclusion_diameter_mm","f0"]);
T_by_diam_field = summarize_metrics(T_patch, ["model_name","inclusion_sws","inclusion_diameter_mm","field_regime"]);
T_worst = sortrows(summarize_metrics(T_patch, ["model_name","condition_id","phantom_name","f0","field_regime"]), 'MAPE_pct', 'descend');
T_worst = T_worst(1:min(40,height(T_worst)),:);

writetable(T_overall, fullfile(OUT.table_dir, 'stage_B_summary_overall.csv'));
writetable(T_by_model, fullfile(OUT.table_dir, 'stage_B_summary_by_model.csv'));
writetable(T_by_phantom, fullfile(OUT.table_dir, 'stage_B_summary_by_phantom.csv'));
writetable(T_by_diam, fullfile(OUT.table_dir, 'stage_B_summary_by_inclusion_diameter.csv'));
writetable(T_by_sws, fullfile(OUT.table_dir, 'stage_B_summary_by_inclusion_sws.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_frequency.csv'));
writetable(T_by_field, fullfile(OUT.table_dir, 'stage_B_summary_by_field_regime.csv'));
writetable(T_by_Dlambda, fullfile(OUT.table_dir, 'stage_B_summary_by_D_over_lambda_bin.csv'));
writetable(T_by_DLwin, fullfile(OUT.table_dir, 'stage_B_summary_by_D_over_Lwin_bin.csv'));
writetable(T_by_Meff, fullfile(OUT.table_dir, 'stage_B_summary_by_Meff_bin.csv'));
writetable(T_by_dist, fullfile(OUT.table_dir, 'stage_B_summary_by_distance_over_window_bin.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir, 'stage_B_summary_by_roi.csv'));
writetable(T_by_roi_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_frequency.csv'));
writetable(T_by_diam_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_diameter_frequency.csv'));
writetable(T_by_diam_field, fullfile(OUT.table_dir, 'stage_B_summary_by_diameter_field.csv'));
writetable(T_worst, fullfile(OUT.table_dir, 'stage_B_worst_conditions.csv'));
end

function Tout = add_bin(T, varName, binName, edges)
Tout = T;
x = T.(varName);
labels = strings(numel(x),1);
for i = 1:numel(edges)-1
    lo = edges(i); hi = edges(i+1);
    if isinf(hi)
        lab = sprintf('>%g', lo);
        mask = x >= lo;
    else
        lab = sprintf('%g-%g', lo, hi);
        mask = x >= lo & x < hi;
    end
    labels(mask) = lab;
end
labels(labels=="") = "unknown";
Tout.(binName) = labels;
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
        mean(X.patch_purity,'omitnan'), ...
        mean(X.distance_to_interface_over_window_radius,'omitnan'), ...
        'VariableNames', {'N_valid_patches','MAPE_pct','median_abs_pct_error', ...
        'MAE_SWS','RMSE_SWS','signed_bias_pct','median_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','mean_predicted_SWS', ...
        'mean_true_SWS','mean_q_pred','mean_patch_purity','mean_distance_over_window_radius'});
end
T = [groups vertcat(rows{:})];
end

function row = condition_log_row(C, status, nRows, elapsed, message)
row = table(string(C.condition_id), string(C.phantom_name), C.background_sws, C.inclusion_sws, C.f0, C.M, ...
    string(C.field_regime), string(C.realism_level), string(status), nRows, elapsed, string(message), ...
    'VariableNames', {'condition_id','phantom_name','background_sws','inclusion_sws','f0','M','field_regime', ...
    'realism_level','status','N_patch_model_rows','elapsed_sec','message'});
end

function MAPS = compact_maps_for_cache(sim, T, C, CFG)
MAPS = struct();
MAPS.condition_id = string(C.condition_id);
MAPS.x = sim.x; MAPS.z = sim.z;
MAPS.true_sws = sim.cs_map;
MAPS.label_map = sim.label_map;
MAPS.real_U = sim.real_U;
MAPS.amplitude = sim.amplitude;
MAPS.phase = sim.phase;
MAPS.sign_cos_phase = sign(cos(sim.phase));
X = T(T.model_name=="q_spectrum_plus_composition",:);
MAPS.q_comp = rows_to_grid(X, 'q_pred');
MAPS.sws_comp = rows_to_grid(X, 'SWS_pred');
MAPS.signed_error_comp = rows_to_grid(X, 'signed_error_percent');
MAPS.abs_error_comp = rows_to_grid(X, 'abs_error_percent');
MAPS.patch_purity = rows_to_grid(X, 'patch_purity');
MAPS.distance_over_window = rows_to_grid(X, 'distance_to_interface_over_window_radius');
end

function tf = should_save_representative_map(C, CFG)
tf = CFG.Figures.SaveRepresentativeMaps && ...
    ismember(string(C.phantom_name), CFG.Figures.RepresentativePhantoms) && ...
    ismember(double(C.f0), CFG.Figures.RepresentativeFrequenciesHz) && ...
    ismember(string(C.field_regime), CFG.Figures.RepresentativeFieldRegimes);
end

function plot_representative_maps(sim, T, C, OUT)
X = T(T.model_name=="q_spectrum_plus_composition",:);
qMap = rows_to_grid(X, 'q_pred');
swsMap = rows_to_grid(X, 'SWS_pred');
signedMap = rows_to_grid(X, 'signed_error_percent');
absMap = rows_to_grid(X, 'abs_error_percent');
purityMap = rows_to_grid(X, 'patch_purity');
distMap = rows_to_grid(X, 'distance_to_interface_over_window_radius');
roiMap = categorical_rows_to_grid(X, 'roi_label');
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 18]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), sim.cs_map, 'True SWS', 'm/s');
plot_map(nexttile(tl), sim.label_map, 'Inclusion labels', 'label');
plot_map(nexttile(tl), real(sim.U), 'Real displacement', 'a.u.');
plot_map(nexttile(tl), sim.amplitude, 'Displacement amplitude', 'a.u.');
plot_map(nexttile(tl), sign(cos(sim.phase)), 'sign(cos phase)', 'sign');
plot_map(nexttile(tl), qMap, 'Predicted q', 'REQ quantile q');
plot_map(nexttile(tl), swsMap, 'Predicted SWS', 'm/s');
plot_map(nexttile(tl), signedMap, 'Signed SWS error', '%');
plot_map(nexttile(tl), absMap, 'Absolute SWS error', '%');
plot_map(nexttile(tl), purityMap, 'Patch purity', 'fraction');
plot_map(nexttile(tl), distMap, 'Distance/window radius', 'ratio');
plot_map(nexttile(tl), roiMap, 'ROI labels', 'ROI code');
title(tl, pretty_condition(C), 'Interpreter','none', 'FontWeight','normal', 'FontSize',10);
outDir = fullfile(OUT.map_dir, char(C.phantom_name), char(C.field_regime));
ensure_dir(outDir);
export_fig(fig, fullfile(outDir, "stage_B_map__" + sanitize(C.condition_id) + ".png"));
end

function plot_source_geometry(sim, Tsrc, C, OUT)
outDir = fullfile(OUT.source_dir, char(C.field_regime));
ensure_dir(outDir);
P = [Tsrc.x_m, Tsrc.y_m, Tsrc.z_m] * 1e3;
Lx = max(sim.x)*1e3; Lz = max(sim.z)*1e3;
fig = figure('Color','w','Units','centimeters','Position',[1 1 25 10]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); view(ax, 35, 22);
box3d(ax, [0 Lx], [-30 30], [0 Lz]);
patch(ax, [0 Lx Lx 0], [0 0 0 0], [0 0 Lz Lz], 'b', 'FaceAlpha',0.05, 'EdgeColor','b', 'LineWidth',1.2);
scatter3(ax, P(:,1), P(:,2), P(:,3), 45, 'r', 'filled');
quiver3(ax, P(:,1), P(:,2), P(:,3), Tsrc.motion_x*0, Tsrc.motion_y*0, Tsrc.motion_z*8, 0, 'k', 'LineWidth',1.3);
for i = 1:size(P,1), text(ax, P(i,1), P(i,2), P(i,3), sprintf(' S%d', i), 'FontSize',7); end
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); zlabel(ax,'z (mm)'); title(ax,'3D source geometry','FontWeight','normal'); axis(ax,'equal');
ax = nexttile(tl);
imagesc(ax, sim.x*1e3, sim.z*1e3, sim.cs_map); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on');
scatter(ax, P(:,1), P(:,3), 45, double(Tsrc.is_in_plane), 'filled', 'MarkerEdgeColor','k');
quiver(ax, P(:,1), P(:,3), Tsrc.motion_x*0, Tsrc.motion_z*8, 0, 'k', 'LineWidth',1.2);
cb = colorbar(ax); ylabel(cb,'True SWS (m/s)'); xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)'); title(ax,'x-z source projection','FontWeight','normal');
sgtitle(fig, pretty_condition(C), 'Interpreter','none', 'FontWeight','normal', 'FontSize',10);
export_fig(fig, fullfile(outDir, "stage_B_sources__" + sanitize(C.condition_id) + ".png"));
end

function box3d(ax, xr, yr, zr)
X = [xr(1) xr(2) xr(2) xr(1) xr(1)]; Y = [yr(1) yr(1) yr(2) yr(2) yr(1)];
plot3(ax, X, Y, zr(1)*ones(size(X)), 'k-'); plot3(ax, X, Y, zr(2)*ones(size(X)), 'k-');
for ix = xr, for iy = yr, plot3(ax, [ix ix], [iy iy], zr, 'k-'); end, end
end

function plot_central_power_spectrum(sim, C, OUT)
Nz = size(sim.U,1); Nx = size(sim.U,2); cz = round(Nz/2); cx = round(Nx/2);
rad = min([80, cz-1, cx-1, Nz-cz, Nx-cx]);
patch = sim.U((cz-rad):(cz+rad), (cx-rad):(cx+rad));
wz = hann_local(size(patch,1)); wx = hann_local(size(patch,2));
P2 = abs(fftshift(fft2(patch .* (wz(:) * wx(:).')))).^2;
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 7]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), abs(patch), 'Central window amplitude', 'a.u.');
plot_map(nexttile(tl), log10(P2 + eps), '2D power spectrum', 'log10 power');
title(tl, "Central-window spectrum: " + pretty_condition(C), 'Interpreter','none', 'FontWeight','normal', 'FontSize',10);
export_fig(fig, fullfile(OUT.spectrum_dir, "stage_B_central_spectrum__" + sanitize(C.condition_id) + ".png"));
end

function w = hann_local(N)
if N <= 1, w = 1; else, n = (0:N-1)'; w = 0.5 - 0.5*cos(2*pi*n/(N-1)); end
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
nz = max(iz); nx = max(ix); Z = nan(nz,nx);
valid = isfinite(iz) & isfinite(ix) & iz>=1 & ix>=1;
Z(sub2ind([nz nx], iz(valid), ix(valid))) = T.(valueName)(valid);
end

function Z = categorical_rows_to_grid(T, valueName)
vals = double(categorical(T.(valueName)));
Tmp = T; Tmp.tmp_stageB_category_code = vals;
Z = rows_to_grid(Tmp, 'tmp_stageB_category_code');
end

function plot_map(ax, Z, ttl, cblabel)
imagesc(ax, Z); axis(ax,'image'); set(ax,'YDir','normal'); title(ax, ttl, 'FontWeight','normal', 'FontSize',8);
cb = colorbar(ax); ylabel(cb, cblabel);
end

function print_summary_block(T, titleText)
fprintf('\n%s\n', titleText);
showVars = intersect(string(T.Properties.VariableNames), ...
    ["model_name","phantom_name","inclusion_sws","inclusion_diameter_mm","f0","field_regime", ...
    "N_valid_patches","MAPE_pct","signed_bias_pct","high_error20_pct"], 'stable');
disp(T(:, cellstr(showVars)));
end

function write_stage_readme(OUT, CFG, T_log, elapsedSec)
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w'); assert(fid > 0);
fprintf(fid, '# Stage B: Resolution Phantom / Geometry Transfer Validation\n\n');
fprintf(fid, '## Objective\n\nEvaluate whether frozen `baseline_minimal_v1` models reconstruct heterogeneous inclusions when diameter is controlled relative to wavelength and REQ window size. No retraining, correction, risk mask, or oracle q is used.\n\n');
fprintf(fid, '## Simulation Design\n\n');
fprintf(fid, '- Phantoms: B1 background 2 m/s with inclusions 3 m/s; B2 background 2 m/s with inclusions 4 m/s.\n');
fprintf(fid, '- Diameters: 10, 14, 24, and 32 mm, arranged in a 2x2 layout inside one large %.0f x %.0f mm field.\n', CFG.Simulation.Lx_m*1e3, CFG.Simulation.Lz_m*1e3);
fprintf(fid, '- Frequencies: %s Hz. Realism: clean only. M=%s, cs_guess=%.1f m/s.\n', strjoin(string(CFG.FrequenciesHz), ', '), mat2str(CFG.MList), CFG.REQ.cs_guess);
fprintf(fid, '- TargetStepM = %.1f mm, matching Stage A to control runtime and keep comparisons aligned.\n', CFG.REQ.TargetStepM*1e3);
fprintf(fid, '- Field regimes: %s. Sources are fixed across phantoms/frequencies. Motion direction is [0 0 1].\n\n', strjoin(CFG.FieldRegimes, ', '));
fprintf(fid, '## Physical Ratios\n\n');
fprintf(fid, '- `lambda_inc = cs_inc / f0`: wavelength inside the inclusion.\n');
fprintf(fid, '- `D_over_lambda = D / lambda_inc`: inclusion diameter in inclusion wavelengths.\n');
fprintf(fid, '- `Lwin = M * cs_guess / f0`: nominal REQ window length.\n');
fprintf(fid, '- `D_over_Lwin = D / Lwin`: inclusion diameter in REQ window lengths.\n');
fprintf(fid, '- `M_eff = M * cs_guess / SWS_true`: effective wavelength count in local material.\n\n');
fprintf(fid, '## Outputs\n\nTables are under `tables/`, figures under `figures/`, and condition caches under `data/condition_cache/`.\n\n');
fprintf(fid, '## Runtime\n\nObserved runner time for this run: %.1f minutes. Cached conditions are reused unless `ADAPTIVE_REQ_EIKONAL_STAGE_B_FORCE_REBUILD=true`.\n\n', elapsedSec/60);
fprintf(fid, 'Parallel controls: `ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_SOURCE_PARFOR=true` parallelizes independent Eikonal sources inside diffuse-like fields. `ADAPTIVE_REQ_EIKONAL_STAGE_B_USE_WINDOW_PARFOR=true` also parallelizes REQ windows, but uses more memory and is optional. This run used source parfor=%d and REQ window parfor=%d.\n\n', CFG.UseSourceParfor, CFG.UseWindowParfor);
fprintf(fid, '## Completed Conditions\n\n');
write_markdown_table(fid, T_log);
fprintf(fid, '\n## Interpretation Template\n\nLarge inclusions should have lower core error than small inclusions. Error should generally decrease as `D/lambda_inc` and `D/Lwin` increase. Interface bands are expected to remain harder than inclusion cores. B2 hard inclusions may retain extra negative bias relative to B1; compare this against the homogeneous 4 m/s Stage A bias.\n');
fclose(fid);
end

function write_docs_readme(root_dir, CFG)
docDir = fullfile(root_dir, 'docs', 'eikonal_validation'); ensure_dir(docDir);
fid = fopen(fullfile(docDir, 'stage_B_resolution_phantom.md'), 'w'); assert(fid > 0);
fprintf(fid, '# Stage B: Resolution Phantom / Geometry Transfer Validation\n\n');
fprintf(fid, 'Stage B evaluates frozen `baseline_minimal_v1` models on clean Eikonal multi-inclusion phantoms. It is designed to separate inclusion core performance from interface-band behavior as a function of inclusion diameter, `D/lambda_inc`, and `D/Lwin`.\n\n');
fprintf(fid, 'Configuration file: `%s`\n\n', CFG.ConfigPath);
fprintf(fid, 'Run commands are listed in `outputs/eikonal_validation/stage_B_resolution_phantom/README_results.md` after the runner is executed.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T), fprintf(fid, '_No rows._\n'); return; end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
maxRows = min(height(T), 80);
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
try, txt = jsonencode(CFG, 'PrettyPrint', true); catch, txt = jsonencode(CFG); end
fid = fopen(path, 'w'); assert(fid > 0); fwrite(fid, txt); fclose(fid);
end

function export_fig(fig, path)
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
savefig(fig, replace(path, '.png', '.fig'));
close(fig);
end

function ensure_dir(d)
if exist(d, 'dir') ~= 7, mkdir(d); end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function s = pretty_condition(C)
s = sprintf('%s | f=%d Hz | %s | M=%d', C.phantom_name, C.f0, C.field_regime, C.M);
end

function y = clamp01(x)
y = min(max(x,0),1);
end

function tf = env_true(name, defaultVal)
s = lower(strtrim(string(getenv(name))));
if s == "", tf = defaultVal; return; end
tf = ismember(s, ["1","true","yes","y","on"]);
end
