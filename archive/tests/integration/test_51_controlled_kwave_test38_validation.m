%% test_51_controlled_kwave_test38_validation.m
% Apply frozen Test38 q models to controlled k-Wave fields from Test 50.
%
% This test does not train any q model. It consumes controlled k-Wave fields,
% extracts REQ/features with the Test38 profile, applies the frozen Test38
% models, and writes patch-level metrics plus map diagnostics.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST51_MODE = validate | quick | full
%   ADAPTIVE_REQ_TEST51_TEST50_SOURCE = quick | full | /path/to/test50/output
%   ADAPTIVE_REQ_TEST51_MODEL_SOURCE = quick | medium | full | /path/to/bundle.mat
%   ADAPTIVE_REQ_TEST51_TARGET_STEP_M = 0.001
%   ADAPTIVE_REQ_TEST51_USE_PARFOR = true | false
%   ADAPTIVE_REQ_TEST51_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST51_MODELS = comma-separated Test38 model names
%   ADAPTIVE_REQ_TEST51_SOURCE_SIDE = all | left | right | top | bottom
%   ADAPTIVE_REQ_TEST51_SOURCE_POLARIZATION = all | axial | lateral | radial | transverse
%   ADAPTIVE_REQ_TEST51_ANALYSIS_ROI = all | exclude_source_buffer | full

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));
run(fullfile(root_dir, 'setup_adaptive_req.m'));

CFG = struct();
CFG.Mode = lower(env_string('ADAPTIVE_REQ_TEST51_MODE', "quick"));
CFG.Test50Source = env_string('ADAPTIVE_REQ_TEST51_TEST50_SOURCE', "full");
CFG.ModelSource = env_string('ADAPTIVE_REQ_TEST51_MODEL_SOURCE', "full");
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST51_TARGET_STEP_M', 1.0e-3);
CFG.UseWindowParfor = env_true('ADAPTIVE_REQ_TEST51_USE_PARFOR', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST51_SAVE_ALL_MAPS', true);
CFG.SourceSideFilter = lower(env_string('ADAPTIVE_REQ_TEST51_SOURCE_SIDE', "all"));
CFG.SourcePolarizationFilter = lower(env_string('ADAPTIVE_REQ_TEST51_SOURCE_POLARIZATION', "all"));
CFG.AnalysisROIFilter = lower(env_string('ADAPTIVE_REQ_TEST51_ANALYSIS_ROI', "all"));
CFG.ModelsToEvaluate = env_string_list('ADAPTIVE_REQ_TEST51_MODELS', ...
    ["q_spectrum_only","q_spectrum_plus_composition","q_spectrum_plus_theory_composition"]);
CFG.cs_guess = 3.0;
CFG.REQ = struct('Gamma',1,'PadFactor',1,'Nbins',"auto", ...
    'Nbins_auto_oversample',1,'Nbins_min',16,'SmoothSigma',1,'EdgeMode',"valid");

OUT = make_output_dirs(root_dir, CFG);
fprintf('\nTest 51: controlled k-Wave validation with frozen Test38 models\n');
fprintf('Mode: %s | Test50 source: %s | model source: %s\n', ...
    CFG.Mode, CFG.Test50Source, CFG.ModelSource);
fprintf('Target REQ step %.3f mm | parfor windows: %d\n', ...
    1e3*CFG.TargetStepM, CFG.UseWindowParfor);
fprintf('No training. True SWS/material maps are evaluation-only.\n');

BUNDLE = load_test38_bundle(root_dir, CFG.ModelSource);
BASE_FEATURES = string(BUNDLE.BASE_FEATURES(:));
assert_no_forbidden_predictors(BASE_FEATURES);
fprintf('Loaded Test38 bundle: %s\n', BUNDLE.bundle_file);
fprintf('Models requested: %s\n', strjoin(CFG.ModelsToEvaluate, ', '));

field_files = list_test50_fields(root_dir, CFG);
fprintf('Matched %d controlled k-Wave fields.\n', numel(field_files));
assert(~isempty(field_files), 'No Test50 field files found.');

if CFG.Mode == "validate"
    validate_one_field(field_files(1), BUNDLE, BASE_FEATURES, CFG, OUT);
    fprintf('Test 51 validation-only checks passed.\n');
    return;
elseif CFG.Mode == "quick" && numel(field_files) > 2
    field_files = quick_subset(field_files);
    fprintf('Quick mode subset: %d fields.\n', numel(field_files));
end

T_all = table();
for fi = 1:numel(field_files)
    [S, key] = load_test50_field(field_files(fi));
    fprintf('\n[%d/%d] %s | component: %s\n', ...
        fi, numel(field_files), key, string_or(S, 'velocity_component', "unknown"));

    M_list = get_M_list(CFG.Mode);
    for mi = 1:numel(M_list)
        M = M_list(mi);
        fprintf('  M=%g\n', M);
        [F, cache_file] = get_or_extract_features(S, key, M, CFG, OUT);
        fprintf('  features: %d windows (%s)\n', height(F), cache_file);
        T = apply_test38_models(F, BUNDLE.MODELS, BASE_FEATURES, CFG);
        T_all = vertcat_compatible(T_all, T);
        if CFG.SaveAllMaps && ~isempty(T)
            plot_condition_maps(T, F, S, key, M, OUT);
        end
    end
end

T_all = add_analysis_regions(T_all);
T_overall = summarize_long(T_all, "strategy_name");
T_by_geometry = summarize_long(T_all, ["strategy_name","geometry"]);
T_by_source = summarize_long(T_all, ["strategy_name","source_mode"]);
T_by_M = summarize_long(T_all, ["strategy_name","M"]);
T_by_side = summarize_long(T_all, ["strategy_name","material_side"]);
T_by_region = summarize_long(T_all, ["strategy_name","analysis_region"]);
T_by_source_distance = summarize_long(T_all, ["strategy_name","source_region"]);

writetable(T_all, fullfile(OUT.table_dir, 'test51_patch_level_predictions.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test51_summary_overall.csv'));
writetable(T_by_geometry, fullfile(OUT.table_dir, 'test51_summary_by_geometry.csv'));
writetable(T_by_source, fullfile(OUT.table_dir, 'test51_summary_by_source_mode.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test51_summary_by_M.csv'));
writetable(T_by_side, fullfile(OUT.table_dir, 'test51_summary_by_material_side.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test51_summary_by_analysis_region.csv'));
writetable(T_by_source_distance, fullfile(OUT.table_dir, 'test51_summary_by_source_distance.csv'));
save(fullfile(OUT.data_dir, 'test51_controlled_kwave_test38_validation.mat'), ...
    'T_all','T_overall','T_by_geometry','T_by_source','T_by_M','T_by_side', ...
    'T_by_region','T_by_source_distance','CFG','-v7.3');

plot_overall_ranking(T_overall, OUT);
plot_group_bars(T_by_geometry, "geometry", OUT, 'test51_mape_by_geometry.png');
plot_group_bars(T_by_source, "source_mode", OUT, 'test51_mape_by_source_mode.png');
plot_group_bars(T_by_M, "M", OUT, 'test51_mape_by_M.png');
plot_group_bars(T_by_region, "analysis_region", OUT, 'test51_mape_by_analysis_region.png');

print_summary(T_overall, T_by_side, T_by_region);
fprintf('\nTables: %s\nFigures: %s\nTest 51 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Main helpers

function OUT = make_output_dirs(root_dir, CFG)
mode_dir = char(CFG.Mode);
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_51_controlled_kwave_test38_validation', mode_dir);
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.cache_dir = fullfile(OUT.data_dir, 'feature_cache');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
dirs = {OUT.root_dir, OUT.data_dir, OUT.cache_dir, OUT.table_dir, OUT.figure_dir, OUT.map_dir};
for i = 1:numel(dirs)
    if ~exist(dirs{i}, 'dir'), mkdir(dirs{i}); end
end
end

function B = load_test38_bundle(root_dir, source)
src = string(source);
if exist(src, 'file') == 2
    file = char(src);
else
    switch lower(src)
        case "full"
            file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
                'models', 'test38_velocity_field_diverse_q_models.mat');
        case "medium"
            file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
                'medium', 'models', 'test38_velocity_field_diverse_q_models.mat');
        case "quick"
            file = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
                'quick', 'models', 'test38_velocity_field_diverse_q_models.mat');
        otherwise
            error('Unknown ADAPTIVE_REQ_TEST51_MODEL_SOURCE=%s.', src);
    end
end
assert(exist(file,'file') == 2, 'Missing Test38 bundle: %s', file);
S = load(file, 'MODEL_BUNDLE');
B = S.MODEL_BUNDLE;
B.bundle_file = file;
end

function files = list_test50_fields(root_dir, CFG)
src = string(CFG.Test50Source);
if exist(src, 'dir') == 7
    base_dir = char(src);
else
    switch lower(src)
        case "full"
            base_dir = fullfile(root_dir, 'outputs', 'test_50_controlled_kwave_phase1', 'full');
        case "quick"
            base_dir = fullfile(root_dir, 'outputs', 'test_50_controlled_kwave_phase1', 'quick');
        otherwise
            error('Unknown ADAPTIVE_REQ_TEST51_TEST50_SOURCE=%s.', src);
    end
end
data_dir = fullfile(base_dir, 'data');
listing = dir(fullfile(data_dir, '*_field.mat'));
files = fullfile(string({listing.folder}), string({listing.name}));
files = sort(files(:));
if CFG.SourceSideFilter ~= "all"
    files = files(contains(files, "__" + CFG.SourceSideFilter + "_"));
end
if CFG.SourcePolarizationFilter ~= "all"
    files = files(contains(files, "_" + CFG.SourcePolarizationFilter + "__"));
end
if CFG.AnalysisROIFilter ~= "all"
    files = files(contains(files, "__" + CFG.AnalysisROIFilter + "__"));
end
end

function files = quick_subset(files)
names = lower(string(files));
keep = contains(names, "homogeneous_cs2__single_sine") | ...
    contains(names, "inclusion_2_3__single_sine");
if any(keep)
    files = files(keep);
else
    files = files(1:min(2,numel(files)));
end
end

function [S, key] = load_test50_field(file)
L = load(char(file));
if isfield(L, 'S')
    S = L.S;
elseif isfield(L, 'Slim')
    S = L.Slim;
else
    error('File does not contain S or Slim: %s', file);
end
[~, key] = fileparts(char(file));
key = erase(string(key), "_field");
end

function M_list = get_M_list(mode)
if string(mode) == "quick"
    M_list = 2;
else
    M_list = [2 3];
end
end

function validate_one_field(file, BUNDLE, BASE_FEATURES, CFG, OUT)
[S, key] = load_test50_field(file);
M = 2;
[F, ~] = get_or_extract_features(S, key, M, CFG, OUT);
F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
P = predict_composition(BUNDLE.MODELS.composition, F, BASE_FEATURES); %#ok<NASGU>
assert(height(F) > 0);
assert(all(isfinite(F.true_SWS)));
assert(all(ismember(BASE_FEATURES, string(F.Properties.VariableNames))));
end

function [F, cache_file] = get_or_extract_features(S, key, M, CFG, OUT)
cache_file = fullfile(OUT.cache_dir, sprintf('%s__M%g__step%gum.mat', ...
    sanitize(key), M, round(1e6*CFG.TargetStepM)));
cache_sig = feature_cache_signature(S, M, CFG); %#ok<NASGU>
if exist(cache_file, 'file') == 2
    C = load(cache_file, 'F', 'cache_sig');
    if isfield(C, 'cache_sig') && isequaln(C.cache_sig, cache_sig)
        F = C.F;
        return;
    end
end

sim_cfg = S.CFG;
cfg_req = struct();
cfg_req.dx = sim_cfg.dx;
cfg_req.dz = sim_cfg.dz;
cfg_req.f0 = sim_cfg.f0;
cfg_req.cs_bg = sim_cfg.cs_soft;
cfg_req.WaveModel = char(string_or(S, 'source_mode', "unknown"));

feat = adaptive_req.config.default_feature_config( ...
    'M', M, 'cs_guess', CFG.cs_guess, ...
    'gamma_win', CFG.REQ.Gamma, 'pad_factor', CFG.REQ.PadFactor);
step_x = max(1, round(CFG.TargetStepM / sim_cfg.dx));
step_z = max(1, round(CFG.TargetStepM / sim_cfg.dz));
req_options = {'Nbins', CFG.REQ.Nbins, ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, 'smooth_sigma', CFG.REQ.SmoothSigma};

t_req = tic;
O = adaptive_req.estimators.req_estimator_map(S.Uxz, cfg_req, feat, ...
    'StepX', step_x, 'StepZ', step_z, ...
    'EdgeMode', CFG.REQ.EdgeMode, 'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, 'UseWindowParfor', CFG.UseWindowParfor, ...
    'StoreReqCurves', false, 'Verbose', false);
F = O.feature_table;
fprintf('    REQ/features extracted in %.1f s | StepX=%d StepZ=%d\n', ...
    toc(t_req), step_x, step_z);

F = attach_metadata(F, S, key, M, feat, step_x, step_z, CFG, O.win_size);
F = ensure_predictor_columns(F, string([]), CFG);
save(cache_file, 'F', 'cache_sig', '-v7.3');
end

function sig = feature_cache_signature(S, M, CFG)
sig = struct();
sig.version = 3;
sig.M = M;
sig.TargetStepM = CFG.TargetStepM;
sig.REQ = CFG.REQ;
sig.field_size = size(S.Uxz);
sig.geometry = string(S.CFG.Geometry);
sig.source_mode = string(S.CFG.SourceMode);
sig.velocity_component = string_or(S, 'velocity_component', "");
sig.f0 = S.CFG.f0;
sig.dx = S.CFG.dx;
sig.dz = S.CFG.dz;
end

function F = attach_metadata(F, S, key, M, feat, step_x, step_z, CFG, win_size)
n = height(F);
sim_cfg = S.CFG;
geometry = string(sim_cfg.Geometry);
source_mode = string(sim_cfg.SourceMode);
origin = analysis_origin(S);
F.x_center_roi_m = F.x_center_m;
F.z_center_roi_m = F.z_center_m;
F.x_center_m = F.x_center_roi_m + origin(1);
F.z_center_m = F.z_center_roi_m + origin(2);
F.condition_key = repmat(string(key), n, 1);
F.geometry = repmat(geometry, n, 1);
F.geometry_type = repmat(geometry_type(geometry), n, 1);
F.source_mode = repmat(source_mode, n, 1);
F.source_side = repmat(string_or(sim_cfg, 'SourceSide', "unknown"), n, 1);
F.source_polarization = repmat(string_or(sim_cfg, 'SourcePolarization', "unknown"), n, 1);
F.field_regime = repmat(source_to_regime(source_mode), n, 1);
F.wave_model = F.field_regime;
F.velocity_component = repmat(string_or(S, 'velocity_component', "unknown"), n, 1);
F.f0 = sim_cfg.f0 * ones(n, 1);
F.SIM_f0 = F.f0;
F.dx = sim_cfg.dx * ones(n, 1);
F.dz = sim_cfg.dz * ones(n, 1);
F.M = M * ones(n, 1);
F.REQ_M = F.M;
F.REQ_cs_guess = feature_cs_guess(feat, CFG.cs_guess) * ones(n, 1);
F.REQ_StepX = step_x * ones(n, 1);
F.REQ_StepZ = step_z * ones(n, 1);
F.TargetStepM = CFG.TargetStepM * ones(n, 1);
F.analysis_origin_x_m = origin(1) * ones(n, 1);
F.analysis_origin_z_m = origin(2) * ones(n, 1);
F.true_SWS = sample_map(S.cs_map, F.x_center_roi_m, F.z_center_roi_m, sim_cfg);
F.cs_true = F.true_SWS;
F.k_true = 2*pi*F.f0 ./ F.true_SWS;
F.material_side = material_side(F.true_SWS, sim_cfg);
F.source_distance_mm = source_distance(F.x_center_m, F.z_center_m, S);
F.source_region = source_region(F.source_distance_mm);
F.distance_to_interface_mm = distance_to_interface(F.x_center_m, F.z_center_m, geometry, sim_cfg);
F.distance_bin = distance_bin(F.distance_to_interface_mm);
F.roi_name = roi_name(F.x_center_m, F.z_center_m, geometry, sim_cfg);
F.patch_purity = patch_purity_from_material(S.material_id, F.cx, F.cz, win_size, sim_cfg);
F.purity_bin = purity_bin(F.patch_purity, geometry);
F.q_theory_prior = repmat(theory_q_for_condition(F, source_mode), n, 1);
F.sws_theory = q_to_sws(F.req_mapping, F.q_theory_prior, F.f0);
F.REQ_Nbins_effective = ensure_nbins_effective(F);
end

function F = ensure_predictor_columns(F, features, CFG)
if nargin < 2, features = string([]); end
features = string(features(:));
for f = features.'
    if ismember(f, string(F.Properties.VariableNames)), continue; end
    switch f
        case "dx", F.dx = F.dx;
        case "dz", F.dz = F.dz;
        case "f0", F.f0 = F.f0;
        case "M", F.M = F.M;
        case "REQ_M", F.REQ_M = F.M;
        case "SIM_f0", F.SIM_f0 = F.f0;
        case "REQ_StepX", F.REQ_StepX = max(1, round(CFG.TargetStepM ./ F.dx));
        case "REQ_StepZ", F.REQ_StepZ = max(1, round(CFG.TargetStepM ./ F.dz));
        case "TargetStepM", F.TargetStepM = CFG.TargetStepM * ones(height(F),1);
        otherwise
            error('Required Test38 predictor missing in Test51 feature table: %s', f);
    end
end
end

function T_out = apply_test38_models(F, MODELS, BASE_FEATURES, CFG)
F = ensure_predictor_columns(F, BASE_FEATURES, CFG);
P = predict_composition(MODELS.composition, F, BASE_FEATURES);
F.predicted_patch_purity = P.predicted_patch_purity;
F.p_mixed = P.p_mixed;
F.p_strong_mixed = P.p_strong_mixed;

parts = {};
for name = string(CFG.ModelsToEvaluate(:)).'
    if name == "theory_discrete"
        q_pred = F.q_theory_prior;
        sws_pred = F.sws_theory;
    elseif isfield(MODELS.q, model_field_name(name))
        Mdl = MODELS.q.(model_field_name(name));
        q_pred = predict_q_model(Mdl, F);
        sws_pred = q_to_sws(F.req_mapping, q_pred, F.f0);
    else
        warning('Skipping unavailable Test38 model: %s', name);
        continue;
    end
    R = base_prediction_table(F, name);
    R.q_pred = q_pred;
    R.sws_pred = sws_pred;
    R.sws_signed_error_pct = 100 * (sws_pred - F.true_SWS) ./ F.true_SWS;
    R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
    R.high_error10 = R.sws_abs_error_pct > 10;
    R.high_error20 = R.sws_abs_error_pct > 20;
    parts{end+1,1} = R; %#ok<AGROW>
end
T_out = vertcat(parts{:});
end

function fld = model_field_name(name)
s = string(name);
switch s
    case "q_spectrum_only"
        fld = 'spectrum_only';
    case "q_spectrum_plus_theory"
        fld = 'spectrum_plus_theory';
    case "q_spectrum_plus_composition"
        fld = 'spectrum_plus_composition';
    case "q_spectrum_plus_theory_composition"
        fld = 'spectrum_plus_theory_composition';
    case "delta_q_theory_composition"
        fld = 'delta_q_theory_composition';
    case "delta_logk_theory_composition"
        fld = 'delta_logk_theory_composition';
    otherwise
        fld = char(s);
end
end

function R = base_prediction_table(F, strategy_name)
keep = ["condition_key","geometry","geometry_type","source_mode","field_regime", ...
    "wave_model","velocity_component","M","REQ_M","f0","dx","dz","REQ_StepX", ...
    "REQ_StepZ","TargetStepM","source_side","source_polarization", ...
    "map_iz","map_ix","x_center_m","z_center_m","x_center_roi_m","z_center_roi_m", ...
    "true_SWS","cs_true","roi_name","material_side","distance_to_interface_mm", ...
    "distance_bin","source_distance_mm","source_region","purity_bin","patch_purity","predicted_patch_purity", ...
    "p_mixed","p_strong_mixed","q_theory_prior","sws_theory"];
keep = keep(ismember(keep, string(F.Properties.VariableNames)));
R = F(:, cellstr(keep));
R.x = F.x_center_m;
R.z = F.z_center_m;
R.model_family = repmat("Test38_frozen_controlled_kWave", height(F), 1);
R.strategy_name = repmat(string(strategy_name), height(F), 1);
end

%% Metrics and plotting

function T = add_analysis_regions(T)
T.analysis_region = repmat("other", height(T), 1);
is_hom = T.geometry_type == "homogeneous";
T.analysis_region(is_hom) = "homogeneous";
T.analysis_region(~is_hom & T.distance_to_interface_mm <= 1) = "interface_0_1mm";
T.analysis_region(~is_hom & T.distance_to_interface_mm > 1 & T.distance_to_interface_mm <= 2) = "interface_1_2mm";
T.analysis_region(~is_hom & T.distance_to_interface_mm > 2 & T.distance_to_interface_mm <= 4) = "near_interface_2_4mm";
T.analysis_region(~is_hom & T.distance_to_interface_mm > 4 & T.material_side == "soft") = "soft_core_gt4mm";
T.analysis_region(~is_hom & T.distance_to_interface_mm > 4 & T.material_side == "hard") = "hard_core_gt4mm";
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
under = splitapply(@(x) 100*mean(x < 0,'omitnan'), T.sws_signed_error_pct, G);
over = splitapply(@(x) 100*mean(x > 0,'omitnan'), T.sws_signed_error_pct, G);
he10 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error10, G);
he20 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
mean_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
S = [groups table(n,mape,medae,signed,under,over,he10,he20,mean_sws, ...
    'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
    'underestimate_pct','overestimate_pct','high_error10_pct','high_error20_pct', ...
    'mean_sws_pred'})];
end

function plot_condition_maps(T, F, S, key, M, OUT)
strategies = unique(T.strategy_name, 'stable');
main = strategies(1:min(numel(strategies), 3));
Tref = T(T.strategy_name == strategies(1), :);
ncols = 4;
nrows = 4;
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 28], 'Visible','off');
tiledlayout(fig, nrows, ncols, 'TileSpacing','compact','Padding','compact');

plot_map(nexttile, S.cs_map, 'True SWS', 'm/s');
plot_map(nexttile, real(S.Uxz), 'Real axial harmonic field', 'a.u.');
plot_map(nexttile, abs(S.Uxz), 'Axial harmonic amplitude', 'a.u.');
plot_map(nexttile, angle(S.Uxz), 'Axial harmonic phase', 'rad');

for i = 1:numel(main)
    Ti = T(T.strategy_name == main(i), :);
    plot_grid_map(nexttile, Ti, Ti.sws_pred, pretty_strategy(main(i)) + " SWS", 'm/s');
end
if numel(main) < 3
    for i = (numel(main)+1):3, axis(nexttile, 'off'); end
end

for i = 1:numel(main)
    Ti = T(T.strategy_name == main(i), :);
    plot_grid_map(nexttile, Ti, Ti.sws_abs_error_pct, pretty_strategy(main(i)) + " abs error", '%');
end
if numel(main) < 3
    for i = (numel(main)+1):3, axis(nexttile, 'off'); end
end
plot_grid_map(nexttile, Tref, Tref.p_mixed, 'Predicted mixedness', 'probability');
plot_grid_map(nexttile, F, F.patch_purity, 'True patch purity', 'fraction');
plot_grid_map(nexttile, Tref, Tref.predicted_patch_purity, 'Predicted patch purity', 'probability');
plot_grid_map(nexttile, F, F.source_distance_mm, 'Distance to source', 'mm');

sgtitle(sprintf('Test 51 controlled k-Wave | %s | M=%g', strrep(key,'_','\_'), M), ...
    'FontWeight','normal');
out_dir = fullfile(OUT.map_dir, sanitize(key));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
export_fig(fig, fullfile(out_dir, sprintf('test51__%s__M%g.png', sanitize(key), M)));
end

function plot_map(ax, Z, ttl, cb_label)
imagesc(ax, Z); axis(ax, 'image'); set(ax, 'YDir', 'normal');
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
cb = colorbar(ax); ylabel(cb, cb_label);
end

function plot_grid_map(ax, T, values, ttl, cb_label)
nz = max(T.map_iz);
nx = max(T.map_ix);
Z = nan(nz, nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
plot_map(ax, Z, ttl, cb_label);
end

function plot_overall_ranking(T, OUT)
T = sortrows(T, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 12], 'Visible','off');
bar(categorical(T.strategy_name), T.MAPE_pct); grid on; xtickangle(30);
ylabel('MAPE (%)');
title('Test 51 controlled k-Wave: frozen Test38 model ranking', 'FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, 'test51_overall_model_ranking.png'));
end

function plot_group_bars(T, group_var, OUT, filename)
groups = unique(string(T.(group_var)), 'stable');
strategies = unique(T.strategy_name, 'stable');
Y = nan(numel(groups), numel(strategies));
for i = 1:numel(groups)
    for j = 1:numel(strategies)
        idx = string(T.(group_var)) == groups(i) & T.strategy_name == strategies(j);
        if any(idx), Y(i,j) = mean(T.MAPE_pct(idx), 'omitnan'); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 13], 'Visible','off');
bar(categorical(groups), Y); grid on; xtickangle(30);
ylabel('MAPE (%)'); legend(strategies, 'Interpreter','none', 'Location','bestoutside');
title("Test 51 MAPE by " + group_var, 'Interpreter','none','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir, filename));
end

function print_summary(T_overall, T_side, T_region)
fprintf('\nSummary:\n');
T = sortrows(T_overall, 'MAPE_pct', 'ascend');
disp(T(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
if ~isempty(T_side)
    fprintf('\nBy material side:\n');
    disp(T_side(:, {'strategy_name','material_side','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
end
if ~isempty(T_region)
    fprintf('\nBy analysis region:\n');
    disp(T_region(:, {'strategy_name','analysis_region','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
end
end

%% Numeric helpers

function P = predict_composition(MIX, T, base_features)
X = T(:, cellstr(base_features));
P = struct();
P.predicted_patch_purity = min(max(predict(MIX.purity, X), 0), 1);
[~, score] = predict(MIX.mixed, X);
P.p_mixed = positive_score(MIX.mixed, score);
[~, score] = predict(MIX.strong_mixed, X);
P.p_strong_mixed = positive_score(MIX.strong_mixed, score);
end

function y = predict_q_model(M, T, do_clamp)
if nargin < 3, do_clamp = true; end
X = T(:, cellstr(M.features));
y = predict(M.model, X);
if do_clamp, y = clamp01(y); end
end

function y = q_to_sws(mappings, q, f0)
n = numel(q);
y = nan(n,1);
for i = 1:n
    if isscalar(f0), f = f0; else, f = f0(i); end
    if isfinite(q(i)) && ~isempty(mappings{i})
        y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), f);
    end
end
end

function q = theory_q_for_condition(F, source_mode)
switch string(source_mode)
    case {"single_sine","single_square"}
        type = "SingleWave";
    case {"sources8_sine","sources128_sine"}
        type = "Diffuse2D";
    otherwise
        type = "Diffuse2D";
end
o = adaptive_req.theory.q_theory_REQ_discrete_shearUZ(F.dx(1), F.dz(1), ...
    F.f0(1), F.REQ_cs_guess(1), 'M', F.REQ_M(1), 'Gamma', 1, ...
    'PadFactor', 1, 'Nbins', 'auto', 'SmoothSigma', 1, ...
    'TheoryMode', 'S2D', 'FieldType', type, 'Plot', false);
q = clamp01(o.q_th);
end

function p = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score, 2); end
p = score(:, idx);
end

function y = clamp01(y)
y = min(max(y, 0), 1);
end

function cs = sample_map(cs_map_zx, x_m, z_m, cfg)
ix = min(max(round(x_m ./ cfg.dx) + 1, 1), size(cs_map_zx, 2));
iz = min(max(round(z_m ./ cfg.dz) + 1, 1), size(cs_map_zx, 1));
cs = cs_map_zx(sub2ind(size(cs_map_zx), iz, ix));
end

function side = material_side(true_sws, cfg)
side = repmat("soft", numel(true_sws), 1);
thr = 0.5 * (cfg.cs_soft + cfg.cs_hard);
side(true_sws > thr) = "hard";
if all(abs(true_sws - true_sws(1)) < eps)
    side(:) = "homogeneous";
end
end

function d = distance_to_interface(x, z, geometry, cfg)
geometry = string(geometry);
d = nan(numel(x), 1);
switch geometry
    case {"inclusion_2_3","circular_inclusion_2_3"}
        d = abs(hypot(x - cfg.inclusion_center_m(1), z - cfg.inclusion_center_m(2)) - ...
            cfg.inclusion_radius_m) * 1e3;
    case "bilayer_2_3"
        d = abs(x - median(cfg.x_m)) * 1e3;
end
end

function d = source_distance(x, z, S)
d = nan(numel(x), 1);
if ~isfield(S, 'source_points_xz') || isempty(S.source_points_xz)
    return;
end
P = S.source_points_xz;
if ~all(ismember(["x_m","z_m"], string(P.Properties.VariableNames)))
    return;
end
px = P.x_m(:).';
pz = P.z_m(:).';
for i = 1:numel(x)
    d(i) = min(hypot(x(i) - px, z(i) - pz)) * 1e3;
end
end

function origin = analysis_origin(S)
if isfield(S, 'analysis_origin_m') && numel(S.analysis_origin_m) >= 2
    origin = double(S.analysis_origin_m(1:2));
else
    origin = [0, 0];
end
end

function r = source_region(d)
r = repmat("source_far", numel(d), 1);
r(isfinite(d) & d <= 10) = "source_5_10mm";
r(isfinite(d) & d <= 5) = "source_0_5mm";
r(~isfinite(d)) = "source_unknown";
end

function r = roi_name(x, z, geometry, cfg)
r = repmat("outside_roi", numel(x), 1);
geometry = string(geometry);
if geometry_type(geometry) == "homogeneous"
    cx = median(cfg.x_m); cz = median(cfg.z_m); half = 4e-3;
    r(abs(x-cx) <= half & abs(z-cz) <= half) = "center_roi";
elseif contains(geometry, "inclusion")
    rr = hypot(x - cfg.inclusion_center_m(1), z - cfg.inclusion_center_m(2));
    r(rr <= max(0, cfg.inclusion_radius_m - 3e-3)) = "inclusion_roi";
    r(rr >= cfg.inclusion_radius_m + 5e-3) = "background_roi";
elseif contains(geometry, "bilayer")
    r(x < median(cfg.x_m) - 4e-3) = "soft_roi";
    r(x > median(cfg.x_m) + 4e-3) = "hard_roi";
end
end

function p = patch_purity_from_material(material_id, cx, cz, win_size, cfg)
p = nan(numel(cx), 1);
if ~isfinite(win_size)
    return;
end
half = floor(win_size/2);
for i = 1:numel(cx)
    xidx = max(1,cx(i)-half):min(size(material_id,2),cx(i)+half);
    zidx = max(1,cz(i)-half):min(size(material_id,1),cz(i)+half);
    patch = material_id(zidx, xidx);
    center_label = material_id(min(max(cz(i),1),size(material_id,1)), ...
        min(max(cx(i),1),size(material_id,2)));
    p(i) = mean(patch(:) == center_label);
end
end

function b = purity_bin(p, geometry)
b = repmat("unknown", numel(p), 1);
if geometry_type(geometry) == "homogeneous"
    b(:) = "homogeneous";
    return;
end
b(p >= 0.99) = "pure";
b(p >= 0.90 & p < 0.99) = "near_pure";
b(p >= 0.70 & p < 0.90) = "moderately_mixed";
b(p < 0.70) = "strongly_mixed";
end

function b = distance_bin(d)
b = repmat("homogeneous_or_unknown", numel(d), 1);
b(isfinite(d) & d <= 8) = "4_8mm";
b(isfinite(d) & d <= 4) = "2_4mm";
b(isfinite(d) & d <= 2) = "1_2mm";
b(isfinite(d) & d <= 1) = "0_1mm";
b(isfinite(d) & d > 8) = "gt_8mm";
end

function t = geometry_type(g)
g = string(g);
if startsWith(g, "homogeneous")
    t = "homogeneous";
elseif contains(g, "inclusion")
    t = "inclusion";
elseif contains(g, "bilayer")
    t = "bilayer";
else
    t = "unknown";
end
end

function r = source_to_regime(s)
s = string(s);
if contains(s, "single")
    r = "directional_2D";
elseif contains(s, "sources")
    r = "diffuse_2D";
else
    r = "unknown";
end
end

function n = ensure_nbins_effective(F)
if ismember("REQ_Nbins_effective", string(F.Properties.VariableNames))
    n = F.REQ_Nbins_effective;
else
    n = nan(height(F),1);
end
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

function val = string_or(S, field, default_val)
if isstruct(S) && isfield(S, field)
    val = string(S.(field));
else
    val = string(default_val);
end
end

%% Boilerplate helpers

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

function s = sanitize(x)
s = regexprep(char(string(x)), '[^A-Za-z0-9_=-]+', '_');
end

function name = pretty_strategy(s)
switch string(s)
    case "q_spectrum_only", name = "q spectrum only";
    case "q_spectrum_plus_composition", name = "q spectrum + composition";
    case "q_spectrum_plus_theory_composition", name = "q spectrum + theory + composition";
    case "theory_discrete", name = "Theory discrete";
    otherwise, name = string(s);
end
end

function export_fig(fig, path)
drawnow;
try
    exportgraphics(fig, path, 'Resolution', 200);
catch
    saveas(fig, path);
end
close(fig);
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
string_vars = ["condition_key","geometry","geometry_type","source_mode","field_regime", ...
    "wave_model","velocity_component","source_side","source_polarization", ...
    "roi_name","material_side","distance_bin", ...
    "source_region","purity_bin","model_family","strategy_name","analysis_region"];
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

function out = env_string(name, default_val)
val = string(getenv(name));
if strlength(val) == 0
    out = string(default_val);
else
    out = val;
end
end

function out = env_number(name, default_val)
val = string(getenv(name));
if strlength(val) == 0
    out = default_val;
else
    out = str2double(val);
    if ~isfinite(out)
        error('Environment variable %s must be numeric.', name);
    end
end
end

function out = env_true(name, default_val)
val = lower(string(getenv(name)));
if strlength(val) == 0
    out = logical(default_val);
else
    out = any(val == ["1", "true", "yes", "y", "on"]);
end
end

function out = env_string_list(name, default_val)
val = string(getenv(name));
if strlength(val) == 0
    out = string(default_val);
else
    parts = split(val, ',');
    out = strip(parts(:)).';
    out = out(strlength(out) > 0);
end
end
