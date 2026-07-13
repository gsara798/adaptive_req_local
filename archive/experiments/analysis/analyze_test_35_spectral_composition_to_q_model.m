%% analyze_test_35_spectral_composition_to_q_model.m
% Test 35: spectral composition -> direct q/SWS model.
%
% This is a clean-model experiment. It deliberately does not use frozen
% Local/Hybrid q models, confidence detectors, Test 30 structures, or any
% previous correction maps as predictors.
%
% Operational inputs:
%   - primitive local spectral / Ecum features from the REQ patch,
%   - physical acquisition parameters (M, f0, dx, dz),
%   - optional analytic TheoryQDiscrete q as a physics prior.
%
% Training labels only:
%   - patch_purity = fraction of the dominant material inside the REQ window,
%     which remains well-defined for 2, 3, or more materials,
%   - mixed labels derived from patch_purity,
%   - q_oracle, obtained by inverting each patch's Ecum(k) at k_true.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST35_MODE          = quick | full
%   ADAPTIVE_REQ_TEST35_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST35_USE_KWAVE     = true | false
%   ADAPTIVE_REQ_TEST35_SAVE_MAPS     = true | false

clear; clc; close all;
format compact;

%% Setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.12, ...
    'defaultAxesLabelFontSizeMultiplier',1.04);

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir,'test35_configuration.json'));

fprintf('\nTest 35: spectral composition -> direct q/SWS model\n');
fprintf('Mode: %s | Validate only: %d | kWave external: %d\n', ...
    CFG.RunMode, CFG.ValidateOnly, CFG.UseKWave);
fprintf('No previous q/confidence/correction model is used as a predictor.\n');
fprintf('Composition target: scalar dominant-material patch purity.\n');

if CFG.ValidateOnly
    validate_test35(root_dir, CFG);
    fprintf('Test 35 validation-only checks passed.\n');
    return;
end

if CFG.KWaveOnly
    run_kwave_only(root_dir, CFG, OUT);
    return;
end

%% Build/load clean synthetic dataset

dataset_file = fullfile(OUT.data_dir, 'test35_synthetic_feature_dataset.mat');
if exist(dataset_file,'file') == 2
    S = load(dataset_file, 'T_syn', 'CFG_saved');
    if isfield(S,'CFG_saved') && cache_signature(S.CFG_saved) == cache_signature(CFG)
        T_syn = S.T_syn;
        fprintf('Reused synthetic feature dataset: %d patches.\n', height(T_syn));
    else
        fprintf('Ignoring incompatible synthetic dataset cache.\n');
        T_syn = build_synthetic_dataset(root_dir, CFG, OUT);
        CFG_saved = CFG; %#ok<NASGU>
        save(dataset_file, 'T_syn', 'CFG_saved', '-v7.3');
    end
else
    T_syn = build_synthetic_dataset(root_dir, CFG, OUT);
    CFG_saved = CFG; %#ok<NASGU>
    save(dataset_file, 'T_syn', 'CFG_saved', '-v7.3');
end

assert(~isempty(T_syn), 'Synthetic dataset is empty.');
assert_policy_no_forbidden_predictors(T_syn);

%% Train composition and q models

[train_mask, test_mask] = condition_split(T_syn, CFG);
fprintf('Train patches: %d | held-out synthetic patches: %d | conditions: %d\n', ...
    sum(train_mask), sum(test_mask), numel(unique(T_syn.condition_key)));

BASE_FEATURES = select_base_predictors(T_syn);
assert(~isempty(BASE_FEATURES), 'No primitive predictor variables were found.');
T_syn.is_train_row = train_mask;
T_syn.is_heldout_row = test_mask;

MODELS = train_models(T_syn, train_mask, BASE_FEATURES, CFG);
T_syn_pred = apply_models(T_syn, MODELS, BASE_FEATURES, "synthetic");
heldout_pred = T_syn_pred(T_syn_pred.is_heldout_row,:);

%% Synthetic summaries

T_comp = composition_metrics(heldout_pred);
T_q_overall = summarize_q_predictions(heldout_pred, "model_name");
T_q_geom = summarize_q_predictions(heldout_pred, ["model_name","geometry"]);
T_q_regime = summarize_q_predictions(heldout_pred, ["model_name","field_regime"]);
T_q_M = summarize_q_predictions(heldout_pred, ["model_name","M"]);
T_q_purity = summarize_q_predictions(heldout_pred, ["model_name","purity_bin"]);

writetable(remove_cell_columns(T_syn_pred), fullfile(OUT.table_dir, ...
    'test35_synthetic_patch_level_predictions.csv'));
writetable(T_comp, fullfile(OUT.table_dir, 'test35_composition_metrics.csv'));
writetable(T_q_overall, fullfile(OUT.table_dir, ...
    'test35_q_model_summary_synthetic_heldout.csv'));
writetable(T_q_geom, fullfile(OUT.table_dir, ...
    'test35_q_model_summary_by_geometry.csv'));
writetable(T_q_regime, fullfile(OUT.table_dir, ...
    'test35_q_model_summary_by_regime.csv'));
writetable(T_q_M, fullfile(OUT.table_dir, ...
    'test35_q_model_summary_by_M.csv'));
writetable(T_q_purity, fullfile(OUT.table_dir, ...
    'test35_q_model_summary_by_purity_bin.csv'));

MODEL_BUNDLE = struct('MODELS', MODELS, 'BASE_FEATURES', BASE_FEATURES, ...
    'CFG', CFG, 'policy', "clean spectral composition model; no frozen model predictors");
save(fullfile(OUT.model_dir, 'test35_spectral_composition_q_models.mat'), ...
    'MODEL_BUNDLE', '-v7.3');

%% Optional k-Wave external transfer

T_kw_pred = table();
if CFG.UseKWave
    try
        T_kw = build_kwave_dataset(root_dir, CFG, OUT);
        if ~isempty(T_kw)
            T_kw_pred = apply_models(T_kw, MODELS, BASE_FEATURES, "kwave");
            T_kw_summary = summarize_q_predictions(T_kw_pred, ["model_name","case_name","M"]);
            writetable(remove_cell_columns(T_kw_pred), fullfile(OUT.table_dir, ...
                'test35_kwave_patch_level_predictions.csv'));
            writetable(T_kw_summary, fullfile(OUT.table_dir, ...
                'test35_kwave_external_summary.csv'));
        end
    catch ME
        warning('kWave external transfer skipped: %s', ME.message);
    end
end

%% Figures and interpretation

safe_plot(@() plot_composition_diagnostics(heldout_pred, OUT), ...
    'composition diagnostics');
safe_plot(@() plot_q_model_summary(T_q_overall, T_q_purity, OUT), ...
    'q model summary');
safe_plot(@() plot_predicted_vs_true_q(heldout_pred, OUT), ...
    'q predicted-vs-oracle');
if CFG.SaveMaps
    safe_plot(@() plot_representative_maps(heldout_pred, OUT, CFG), ...
        'representative maps');
end
if ~isempty(T_kw_pred)
    safe_plot(@() plot_kwave_summary(T_kw_pred, OUT), ...
        'kWave summary');
end

print_summary(T_comp, T_q_overall, T_q_purity, T_kw_pred);

fprintf('\nTables: %s\nFigures: %s\nModels: %s\n', ...
    OUT.table_dir, OUT.figure_dir, OUT.model_dir);
fprintf('Test 35 complete.\n');

%% Local functions

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST35_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), ...
    'ADAPTIVE_REQ_TEST35_MODE must be quick or full.');

CFG = struct();
CFG.RunMode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST35_VALIDATE_ONLY', false);
CFG.UseKWave = env_true('ADAPTIVE_REQ_TEST35_USE_KWAVE', true);
CFG.SaveMaps = env_true('ADAPTIVE_REQ_TEST35_SAVE_MAPS', true);
CFG.KWaveOnly = env_true('ADAPTIVE_REQ_TEST35_KWAVE_ONLY', false);
CFG.Version = 3;
CFG.RandomSeed = 35001;
CFG.Dx = 0.2e-3;
CFG.TargetStepM = 1.0e-3;
CFG.MaxPatchesPerCondition = 900;
CFG.MaxKWavePatchesPerCondition = env_double( ...
    'ADAPTIVE_REQ_TEST35_MAX_KWAVE_PATCHES', 1200);
CFG.TrainFraction = 0.70;
CFG.PurityMixedThreshold = 0.95;
CFG.PurityStrongMixedThreshold = 0.75;
CFG.TreeLearners = ternary(CFG.QuickMode, 60, 160);
CFG.MinLeafSize = ternary(CFG.QuickMode, 12, 8);
CFG.UseParfor = false;
CFG.Geometries = ["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"];
if CFG.QuickMode
    CFG.Frequencies = 500;
    CFG.M = [2 3];
    CFG.Regimes = ["directional_2D","diffuse_3D"];
else
    CFG.Frequencies = [300 400 500 600];
    CFG.M = [2 3 4];
    CFG.Regimes = ["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
end

CFG.FieldCacheDir = fullfile(root_dir,'outputs', ...
    'test_22_confidence_external_validation','analysis','data','field_cache');
CFG.KWaveFile = fullfile(root_dir,'outputs','test_12_cs_guess_window_sweep', ...
    'test_12_cs_guess_window_sweep_2026-06-10_122903','analysis', ...
    'level_06_kwave_transfer','data','level12_level06_kwave_transfer.mat');
end

function run_kwave_only(root_dir, CFG, OUT)
assert(CFG.UseKWave, 'Set ADAPTIVE_REQ_TEST35_USE_KWAVE=true for kWave-only mode.');
model_file = fullfile(OUT.model_dir, 'test35_spectral_composition_q_models.mat');
assert(exist(model_file,'file') == 2, ...
    'Missing saved Test 35 model bundle: %s', model_file);
S = load(model_file, 'MODEL_BUNDLE');
MODELS = S.MODEL_BUNDLE.MODELS;
BASE_FEATURES = S.MODEL_BUNDLE.BASE_FEATURES;
fprintf('\nTest 35 kWave-only transfer from saved models.\n');
fprintf('Max kWave patches per condition: %g\n', CFG.MaxKWavePatchesPerCondition);
T_kw = build_kwave_dataset(root_dir, CFG, OUT);
assert(~isempty(T_kw), 'No kWave rows were built.');
T_kw_pred = apply_models(T_kw, MODELS, BASE_FEATURES, "kwave");
T_kw_summary = summarize_q_predictions(T_kw_pred, ["model_name","case_name","M"]);
writetable(remove_cell_columns(T_kw_pred), fullfile(OUT.table_dir, ...
    'test35_kwave_patch_level_predictions.csv'));
writetable(T_kw_summary, fullfile(OUT.table_dir, ...
    'test35_kwave_external_summary.csv'));
save(fullfile(OUT.data_dir,'test35_kwave_transfer_results.mat'), ...
    'T_kw_pred','T_kw_summary','CFG','-v7.3');
safe_plot(@() plot_kwave_summary(T_kw_pred, OUT), 'kWave summary');
S_overall = summarize_q_predictions(T_kw_pred, "model_name");
fprintf('\nkWave-only summary:\n');
disp(S_overall(:, {'model_name','N','MAPE_pct','mean_signed_error_pct', ...
    'high_error20_pct'}));
fprintf('\nUpdated kWave tables:\n%s\n', OUT.table_dir);
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_35_spectral_composition_to_q_model');
if CFG.QuickMode
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
end
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.model_dir = fullfile(OUT.root_dir, 'models');
OUT.condition_dir = fullfile(OUT.data_dir, 'condition_features');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if exist(dirs(i),'dir') ~= 7, mkdir(dirs(i)); end
end
end

function validate_test35(root_dir, CFG)
assert(exist(CFG.FieldCacheDir,'dir') == 7, ...
    'Missing Test 22 field cache: %s', CFG.FieldCacheDir);
files = dir(fullfile(CFG.FieldCacheDir, 'field__*.mat'));
assert(~isempty(files), 'No field caches found in %s', CFG.FieldCacheDir);
sample = fullfile(files(1).folder, files(1).name);
S = load(sample, 'sim', 'cfg_sim');
assert(isfield(S.sim,'Uxz') && isfield(S.sim,'cs_map'), ...
    'Field cache must contain sim.Uxz and sim.cs_map.');
feat_cfg = adaptive_req.config.default_feature_config('M', 2, ...
    'cs_guess', 3, 'gamma_win', 1, 'pad_factor', 1);
OUT = adaptive_req.estimators.req_estimator_map(S.sim.Uxz, S.cfg_sim, feat_cfg, ...
    'StepX', max(1,round(CFG.TargetStepM/S.cfg_sim.dx)), ...
    'StepZ', max(1,round(CFG.TargetStepM/S.cfg_sim.dz)), ...
    'StoreReqCurves', false, 'ReturnFeatures', false, ...
    'ReturnFeatureTable', true, 'ReuseReqSpectrumForFeatures', true, ...
    'UseWindowParfor', false);
T = OUT.feature_table;
assert(ismember('req_mapping', T.Properties.VariableNames), ...
    'REQ feature table did not include req_mapping.');
q = invert_mapping_to_q(T.req_mapping{1}, 2*pi*S.cfg_sim.f0/2);
assert(isfinite(q), 'q_oracle inversion failed on validation sample.');
forbidden = ["confidence","patch_purity","true_SWS","sws_error", ...
    "material_side","distance_to_interface","q_local_req","q_pred","cs_pred"];
fprintf('Validation checked %s with %d feature rows.\n', files(1).name, height(T));
fprintf('Forbidden previous-model/oracle predictors are excluded by name: %s\n', ...
    strjoin(forbidden, ', '));
end

function T_all = build_synthetic_dataset(root_dir, CFG, OUT)
field_dir = CFG.FieldCacheDir;
assert(exist(field_dir,'dir') == 7, 'Missing field cache: %s', field_dir);
parts = {};
part_idx = 0;
total_conditions = numel(CFG.Geometries) * numel(CFG.Frequencies) * ...
    numel(CFG.Regimes) * numel(CFG.M);
ci = 0;
for gi = 1:numel(CFG.Geometries)
    for fi = 1:numel(CFG.Frequencies)
        for ri = 1:numel(CFG.Regimes)
            sim_key = sprintf('%s__f%g__%s__dx%gum', CFG.Geometries(gi), ...
                CFG.Frequencies(fi), lower(CFG.Regimes(ri)), round(1e6*CFG.Dx));
            field_file = fullfile(field_dir, "field__" + sanitize(sim_key) + ".mat");
            if exist(field_file,'file') ~= 2
                warning('Missing field cache, skipping: %s', field_file);
                continue;
            end
            S = load(field_file, 'sim', 'cfg_sim');
            for mi = 1:numel(CFG.M)
                ci = ci + 1;
                M = CFG.M(mi);
                condition_key = string(sim_key) + "__M" + string(M);
                checkpoint = fullfile(OUT.condition_dir, "features__" + ...
                    sanitize(condition_key) + ".mat");
                if exist(checkpoint,'file') == 2
                    C = load(checkpoint, 'T_condition');
                    T_condition = C.T_condition;
                    has_good_purity = ismember('patch_purity', T_condition.Properties.VariableNames) && ...
                        any(isfinite(T_condition.patch_purity));
                    has_good_theory = ismember('q_theory_prior', T_condition.Properties.VariableNames) && ...
                        any(abs(T_condition.q_theory_prior - 0.5) > 1e-6);
                    if has_good_purity && has_good_theory
                        fprintf('[%d/%d] Reused %s (%d patches).\n', ...
                            ci, total_conditions, condition_key, height(T_condition));
                    else
                        fprintf('[%d/%d] Regenerating stale %s (purity/theory cache).\n', ...
                            ci, total_conditions, condition_key);
                        T_condition = extract_condition_features(S.sim, S.cfg_sim, ...
                            CFG.Geometries(gi), CFG.Regimes(ri), M, CFG, condition_key);
                        save(checkpoint, 'T_condition', '-v7.3');
                    end
                else
                    timer = tic;
                    T_condition = extract_condition_features(S.sim, S.cfg_sim, ...
                        CFG.Geometries(gi), CFG.Regimes(ri), M, CFG, condition_key);
                    save(checkpoint, 'T_condition', '-v7.3');
                    fprintf('[%d/%d] Built %s: %d sampled patches in %.1f s.\n', ...
                        ci, total_conditions, condition_key, height(T_condition), toc(timer));
                end
                part_idx = part_idx + 1;
                parts{part_idx,1} = T_condition; %#ok<AGROW>
            end
        end
    end
end
assert(~isempty(parts), 'No synthetic conditions were built.');
T_all = vertcat(parts{:});
T_all.dataset = repmat("synthetic", height(T_all), 1);
T_all = movevars(T_all, 'dataset', 'Before', 1);
end

function T = extract_condition_features(sim, cfg_sim, geometry, regime, M, CFG, condition_key)
feat_cfg = adaptive_req.config.default_feature_config('M', M, ...
    'cs_guess', 3, 'gamma_win', 1, 'pad_factor', 1);
step_x = max(1, round(CFG.TargetStepM / cfg_sim.dx));
step_z = max(1, round(CFG.TargetStepM / cfg_sim.dz));
OUT = adaptive_req.estimators.req_estimator_map(sim.Uxz, cfg_sim, feat_cfg, ...
    'StepX', step_x, 'StepZ', step_z, 'EdgeMode', 'valid', ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', {'Nbins','auto','Nbins_auto_oversample',1, ...
    'Nbins_min',16,'smooth_sigma',1}, ...
    'StoreReqCurves', false, 'ReturnFeatures', false, ...
    'ReturnFeatureTable', true, 'ReuseReqSpectrumForFeatures', true, ...
    'UseWindowParfor', CFG.UseParfor);
T = OUT.feature_table;
if height(T) > CFG.MaxPatchesPerCondition
    rng(CFG.RandomSeed + sum(double(char(condition_key))), 'twister');
    T = T(sort(randperm(height(T), CFG.MaxPatchesPerCondition)), :);
end

half_win = floor(OUT.win_size / 2);
n = height(T);
true_sws = nan(n,1);
purity = nan(n,1);
for i = 1:n
    cz = T.cz(i); cx = T.cx(i);
    true_sws(i) = sim.cs_map(cz, cx);
    z_idx = (cz-half_win):(cz+half_win);
    x_idx = (cx-half_win):(cx+half_win);
    purity(i) = dominant_material_fraction(sim.cs_map(z_idx, x_idx));
end
q_theory = theory_q_for_condition(cfg_sim, feat_cfg, regime);
q_oracle = nan(n,1);
for i = 1:n
    k_true = 2*pi*cfg_sim.f0 / true_sws(i);
    q_oracle(i) = invert_mapping_to_q(T.req_mapping{i}, k_true);
end
T.condition_key = repmat(string(condition_key), n, 1);
T.geometry = repmat(string(geometry), n, 1);
T.field_regime = repmat(string(regime), n, 1);
T.f0 = cfg_sim.f0 * ones(n,1);
T.dx = cfg_sim.dx * ones(n,1);
T.dz = cfg_sim.dz * ones(n,1);
T.M = M * ones(n,1);
T.true_SWS = true_sws;
T.k_true = 2*pi*cfg_sim.f0 ./ true_sws;
T.patch_purity = purity;
T.is_mixed = purity < CFG.PurityMixedThreshold;
T.is_strong_mixed = purity < CFG.PurityStrongMixedThreshold;
T.purity_bin = purity_bin(purity);
T.q_theory_prior = repmat(q_theory, n, 1);
T.q_oracle = q_oracle;
T.sws_theory = q_to_sws(T.req_mapping, T.q_theory_prior, cfg_sim.f0);
T.sws_oracle = q_to_sws(T.req_mapping, T.q_oracle, cfg_sim.f0);
T = movevars(T, {'condition_key','geometry','field_regime','f0','M','dx','dz', ...
    'true_SWS','patch_purity','purity_bin','is_mixed','is_strong_mixed', ...
    'q_theory_prior','q_oracle','sws_theory','sws_oracle'}, 'Before', 1);
end

function T_kw = build_kwave_dataset(root_dir, CFG, OUT)
if exist(CFG.KWaveFile,'file') ~= 2
    warning('Missing kWave cache: %s', CFG.KWaveFile);
    T_kw = table();
    return;
end
S = load(CFG.KWaveFile, 'CASE_OUT');
parts = {};
for ci = 1:numel(S.CASE_OUT)
    C = S.CASE_OUT(ci);
    T = C.T_feat;
    if ~ismember(C.REQ_M, CFG.M), continue; end
    if height(T) > CFG.MaxKWavePatchesPerCondition
        rng(CFG.RandomSeed + ci, 'twister');
        T = T(sort(randperm(height(T), CFG.MaxKWavePatchesPerCondition)), :);
    end
    n = height(T);
    f0 = C.cfg.f0;
    true_sws = infer_kwave_true_sws(T, C);
    q_theory = theory_q_for_condition(C.cfg, ...
        adaptive_req.config.default_feature_config('M', C.REQ_M, ...
        'cs_guess', 3, 'gamma_win', 1, 'pad_factor', 1), "diffuse_3D");
    q_oracle = nan(n,1);
    for i = 1:n
        q_oracle(i) = invert_mapping_to_q(T.req_mapping{i}, 2*pi*f0/true_sws(i));
    end
    T.condition_key = repmat("kwave__" + string(C.case.case_name) + ...
        "__M" + string(C.REQ_M), n, 1);
    T.case_name = repmat(string(C.case.case_name), n, 1);
    T.geometry = repmat("kwave_circular_inclusion_2_3", n, 1);
    T.field_regime = repmat(canonical_kwave_regime(C.case.case_name), n, 1);
    T.f0 = f0 * ones(n,1);
    T.dx = C.cfg.dx * ones(n,1);
    T.dz = C.cfg.dz * ones(n,1);
    T.M = C.REQ_M * ones(n,1);
    T.true_SWS = true_sws;
    T.k_true = 2*pi*f0 ./ true_sws;
    T.q_theory_prior = repmat(q_theory, n, 1);
    T.q_oracle = q_oracle;
    T.sws_theory = q_to_sws(T.req_mapping, T.q_theory_prior, f0);
    T.sws_oracle = q_to_sws(T.req_mapping, T.q_oracle, f0);
    T.patch_purity = nan(n,1);
    T.purity_bin = repmat("unknown", n, 1);
    T.is_mixed = false(n,1);
    T.is_strong_mixed = false(n,1);
    parts{end+1,1} = T; %#ok<AGROW>
end
if isempty(parts)
    T_kw = table();
else
    T_kw = vertcat(parts{:});
    T_kw.dataset = repmat("kwave", height(T_kw), 1);
    T_kw = movevars(T_kw, 'dataset', 'Before', 1);
    save(fullfile(OUT.data_dir,'test35_kwave_feature_dataset.mat'), 'T_kw', '-v7.3');
end
end

function MODELS = train_models(T, train_mask, base_features, CFG)
Fbase = T(:, cellstr(base_features));
Ypurity = T.patch_purity;
Ymix = logical(T.is_mixed);
Ystrong = logical(T.is_strong_mixed);

template = templateTree('MinLeafSize', CFG.MinLeafSize);
fprintf('Training composition models with %d primitive predictors.\n', numel(base_features));
MIX = struct();
MIX.base_features = base_features;
MIX.purity = fitrensemble(Fbase(train_mask,:), Ypurity(train_mask), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template);
MIX.mixed = fitcensemble(Fbase(train_mask,:), Ymix(train_mask), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template);
MIX.strong_mixed = fitcensemble(Fbase(train_mask,:), Ystrong(train_mask), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template);

P = predict_composition(MIX, T, base_features);
Ttrain = T;
Ttrain.predicted_patch_purity = P.predicted_patch_purity;
Ttrain.p_mixed = P.p_mixed;
Ttrain.p_strong_mixed = P.p_strong_mixed;

Q = struct();
Q.spectrum_only = fit_q_regressor(Ttrain, train_mask, base_features, ...
    T.q_oracle, CFG);
Q.spectrum_plus_theory = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "q_theory_prior"], T.q_oracle, CFG);
Q.spectrum_plus_composition = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "predicted_patch_purity"; "p_mixed"; "p_strong_mixed"], ...
    T.q_oracle, CFG);
Q.spectrum_plus_theory_composition = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "q_theory_prior"; "predicted_patch_purity"; ...
    "p_mixed"; "p_strong_mixed"], T.q_oracle, CFG);
Q.delta_q_theory_composition = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "q_theory_prior"; "predicted_patch_purity"; ...
    "p_mixed"; "p_strong_mixed"], T.q_oracle - T.q_theory_prior, CFG);
Q.delta_logk_theory_composition = fit_q_regressor(Ttrain, train_mask, ...
    [base_features; "q_theory_prior"; "predicted_patch_purity"; ...
    "p_mixed"; "p_strong_mixed"], log(T.k_true) - log(2*pi*T.f0./T.sws_theory), CFG);

MODELS = struct('composition', MIX, 'q', Q, ...
    'model_names', ["theory_discrete","q_spectrum_only", ...
    "q_spectrum_plus_theory","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition","delta_q_theory_composition", ...
    "delta_logk_theory_composition"]);
end

function M = fit_q_regressor(T, train_mask, features, y, CFG)
valid = train_mask & all(isfinite(table2array(T(:, cellstr(features)))),2) & isfinite(y);
template = templateTree('MinLeafSize', CFG.MinLeafSize);
M = struct();
M.features = features(:);
M.model = fitrensemble(T(valid, cellstr(features)), y(valid), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template);
M.n_train = sum(valid);
end

function T_out = apply_models(T, MODELS, base_features, dataset_name)
P = predict_composition(MODELS.composition, T, base_features);
T_aug = T;
T_aug.predicted_patch_purity = P.predicted_patch_purity;
T_aug.p_mixed = P.p_mixed;
T_aug.p_strong_mixed = P.p_strong_mixed;

parts = cell(numel(MODELS.model_names),1);
for mi = 1:numel(MODELS.model_names)
    name = MODELS.model_names(mi);
    R = T_aug(:, intersect(string(T_aug.Properties.VariableNames), ...
        ["dataset","condition_key","case_name","geometry","field_regime", ...
        "f0","M","dx","dz","map_iz","map_ix","cx","cz","x_center_m", ...
        "z_center_m","true_SWS","k_true","patch_purity","purity_bin", ...
        "is_mixed","is_strong_mixed","q_oracle","q_theory_prior", ...
        "sws_theory","is_train_row","is_heldout_row", ...
        "predicted_patch_purity","p_mixed","p_strong_mixed"], ...
        'stable'));
    R.model_name = repmat(name, height(R), 1);
    R.dataset = repmat(dataset_name, height(R), 1);
    switch name
        case "theory_discrete"
            q_pred = T_aug.q_theory_prior;
            sws_pred = T_aug.sws_theory;
        case "q_spectrum_only"
            q_pred = predict_q_model(MODELS.q.spectrum_only, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "q_spectrum_plus_theory_composition"
            q_pred = predict_q_model(MODELS.q.spectrum_plus_theory_composition, T_aug);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "delta_q_theory_composition"
            delta = predict_q_model(MODELS.q.delta_q_theory_composition, T_aug);
            q_pred = clamp01(T_aug.q_theory_prior + delta);
            sws_pred = q_to_sws(T_aug.req_mapping, q_pred, T_aug.f0);
        case "delta_logk_theory_composition"
            delta = predict_q_model(MODELS.q.delta_logk_theory_composition, T_aug);
            k_theory = 2*pi*T_aug.f0 ./ T_aug.sws_theory;
            k_pred = k_theory .* exp(delta);
            sws_pred = 2*pi*T_aug.f0 ./ k_pred;
            q_pred = arrayfun(@(i) invert_mapping_to_q(T_aug.req_mapping{i}, k_pred(i)), ...
                (1:height(T_aug))');
        otherwise
            error('Unknown model: %s', name);
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
T_out = vertcat(parts{:});
end

function P = predict_composition(MIX, T, features)
F = T(:, cellstr(features));
P = struct();
P.predicted_patch_purity = min(max(predict(MIX.purity, F), 0), 1);
[~, score] = predict(MIX.mixed, F);
P.p_mixed = positive_score(MIX.mixed, score);
[~, score] = predict(MIX.strong_mixed, F);
P.p_strong_mixed = positive_score(MIX.strong_mixed, score);
end

function y = predict_q_model(M, T)
X = T(:, cellstr(M.features));
y = predict(M.model, X);
y = clamp01(y);
end

function features = select_base_predictors(T)
vars = string(T.Properties.VariableNames);
numeric = false(size(vars));
for i = 1:numel(vars)
    numeric(i) = isnumeric(T.(vars(i))) || islogical(T.(vars(i)));
end
forbidden_patterns = ["true","oracle","purity","mixed","confidence", ...
    "error","pred","sws","cs_","k_true","q_local","q_pred", ...
    "q_theory","req_mapping","patch_idx","map_ix","map_iz","cx","cz", ...
    "x_center","z_center","condition"];
allowed = numeric;
for p = forbidden_patterns
    allowed = allowed & ~contains(lower(vars), lower(p));
end
keep_names = ["REQ_M","M","SIM_f0","f0","dx","dz","REQ_Nbins_effective"];
allowed = allowed | ismember(vars, keep_names);
features = vars(allowed);
features = features(:);
end

function [train_mask, test_mask] = condition_split(T, CFG)
rng(CFG.RandomSeed, 'twister');
cond = unique(T(:, {'condition_key','geometry','field_regime','f0','M'}), 'rows', 'stable');
train_cond = false(height(cond),1);
groups = unique(cond(:, {'geometry','field_regime'}), 'rows', 'stable');
for gi = 1:height(groups)
    idx = find(cond.geometry == groups.geometry(gi) & ...
        cond.field_regime == groups.field_regime(gi));
    idx = idx(randperm(numel(idx)));
    ntrain = max(1, round(CFG.TrainFraction * numel(idx)));
    train_cond(idx(1:ntrain)) = true;
end
train_keys = cond.condition_key(train_cond);
train_mask = ismember(T.condition_key, train_keys);
test_mask = ~train_mask;
if ~any(test_mask)
    warning('Condition split produced no test rows; using a 70/30 row split fallback.');
    idx = randperm(height(T));
    train_mask = false(height(T),1);
    train_mask(idx(1:round(CFG.TrainFraction*height(T)))) = true;
    test_mask = ~train_mask;
end
end

function T = composition_metrics(T)
y = logical(T.is_mixed);
ys = logical(T.is_strong_mixed);
T = table( ...
    height(T), ...
    rmse(T.patch_purity, T.predicted_patch_purity), ...
    mean(abs(T.patch_purity - T.predicted_patch_purity), 'omitnan'), ...
    auc_binary(y, T.p_mixed, "roc"), ...
    auc_binary(y, T.p_mixed, "pr"), ...
    auc_binary(ys, T.p_strong_mixed, "roc"), ...
    auc_binary(ys, T.p_strong_mixed, "pr"), ...
    'VariableNames', {'N','purity_RMSE','purity_MAE', ...
    'mixed_ROC_AUC','mixed_PR_AUC','strong_mixed_ROC_AUC', ...
    'strong_mixed_PR_AUC'});
end

function S = summarize_q_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(group_vars)
    G = ones(height(T),1); groups = table();
else
    [G, groups] = findgroups(T(:, group_vars));
end
nG = max(G);
rows = cell(nG,1);
for gi = 1:nG
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
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','median_signed_error_pct', ...
        'underestimate_pct','high_error10_pct','high_error20_pct', ...
        'mean_abs_q_error'});
end
S = [groups vertcat(rows{:})];
end

function plot_composition_diagnostics(T, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact');
ax = nexttile(tl);
scatter(ax, T.patch_purity, T.predicted_patch_purity, 8, T.p_mixed, ...
    'filled', 'MarkerFaceAlpha', 0.20); axis(ax,'square'); grid(ax,'on');
xlabel(ax,'True patch purity'); ylabel(ax,'Predicted patch purity');
title(ax,'Composition regressor','FontWeight','normal'); colorbar(ax);
ax = nexttile(tl);
edges = linspace(0,1,21);
histogram(ax, T.predicted_patch_purity(T.is_mixed), edges, ...
    'Normalization','probability'); hold(ax,'on');
histogram(ax, T.predicted_patch_purity(~T.is_mixed), edges, ...
    'Normalization','probability');
legend(ax, {'mixed true','pure true'}, 'Location','northwest');
xlabel(ax,'Predicted purity'); ylabel(ax,'Probability'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test35_composition_diagnostics.png'));
end

function plot_q_model_summary(T_overall, T_purity, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 21 9]);
ax = axes(fig);
bar(ax, categorical(T_overall.model_name), T_overall.MAPE_pct);
xtickangle(ax,30); ylabel(ax,'Held-out MAPE (%)'); grid(ax,'on');
title(ax,'Test 35 q/SWS model comparison','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir,'test35_q_model_mape_summary.png'));

fig = figure('Color','w','Units','centimeters','Position',[2 2 22 10]);
models = unique(T_purity.model_name,'stable');
bins = unique(T_purity.purity_bin,'stable');
C = nan(numel(models), numel(bins));
for mi = 1:numel(models)
    for bi = 1:numel(bins)
        idx = T_purity.model_name == models(mi) & T_purity.purity_bin == bins(bi);
        if any(idx), C(mi,bi) = T_purity.MAPE_pct(idx); end
    end
end
heatmap(string(bins), string(models), C, 'ColorbarVisible','on');
title('MAPE by purity bin');
export_fig(fig, fullfile(OUT.figure_dir,'test35_mape_by_purity_bin.png'));
end

function plot_predicted_vs_true_q(T, OUT)
models = unique(T.model_name,'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 14]);
tl = tiledlayout(fig,2,ceil(numel(models)/2),'TileSpacing','compact');
for mi = 1:numel(models)
    ax = nexttile(tl);
    X = T(T.model_name == models(mi),:);
    scatter(ax, X.q_oracle, X.q_pred, 6, X.patch_purity, ...
        'filled', 'MarkerFaceAlpha', 0.18);
    hold(ax,'on'); plot(ax,[0 1],[0 1],'k:'); axis(ax,'square');
    xlabel(ax,'q oracle'); ylabel(ax,'q predicted'); grid(ax,'on');
    title(ax, models(mi), 'Interpreter','none', 'FontWeight','normal');
end
export_fig(fig, fullfile(OUT.figure_dir,'test35_q_pred_vs_oracle.png'));
end

function plot_representative_maps(T, OUT, CFG)
main = "q_spectrum_plus_theory_composition";
geoms = unique(T.geometry,'stable');
for gi = 1:min(numel(geoms), 4)
    X = T(T.geometry == geoms(gi) & T.model_name == main,:);
    if isempty(X), continue; end
    key = X.condition_key(1);
    X = X(X.condition_key == key,:);
    [Ztrue,nz,nx] = vector_to_grid(X.true_SWS, X.map_iz, X.map_ix);
    Zpred = vector_to_grid(X.sws_pred, X.map_iz, X.map_ix, nz, nx);
    Zerr = abs(100*(Zpred-Ztrue)./Ztrue);
    Zpur = vector_to_grid(X.predicted_patch_purity, X.map_iz, X.map_ix, nz, nx);
    fig = figure('Color','w','Units','centimeters','Position',[2 2 24 8]);
    tl = tiledlayout(fig,1,4,'TileSpacing','compact');
    maps = {Ztrue,Zpred,Zerr,Zpur};
    titles = {'True SWS','Predicted SWS','Abs error %','Predicted purity'};
    for i = 1:4
        ax = nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); colorbar(ax);
        title(ax,titles{i},'FontWeight','normal');
    end
    title(tl, key, 'Interpreter','none');
    export_fig(fig, fullfile(OUT.figure_dir, "test35_map__" + sanitize(key) + ".png"));
end
end

function plot_kwave_summary(T, OUT)
S = summarize_q_predictions(T, ["model_name","case_name","M"]);
fig = figure('Color','w','Units','centimeters','Position',[2 2 22 10]);
main = S(S.model_name == "q_spectrum_plus_theory_composition",:);
if isempty(main), main = S; end
cases = unique(main.case_name, 'stable');
Mvals = unique(main.M, 'stable');
C = nan(numel(cases), numel(Mvals));
for ci = 1:numel(cases)
    for mi = 1:numel(Mvals)
        idx = main.case_name == cases(ci) & main.M == Mvals(mi);
        if any(idx)
            C(ci, mi) = mean(main.MAPE_pct(idx), 'omitnan');
        end
    end
end
bar(categorical(cases), C); grid on;
ylabel('kWave MAPE (%)');
legend("M=" + string(Mvals), 'Location','best');
title('Test 35 kWave external summary','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir,'test35_kwave_external_summary.png'));
end

function safe_plot(fun, label)
try
    fun();
catch ME
    warning('Skipping %s figure: %s', label, ME.message);
end
end

function print_summary(T_comp, T_q, T_purity, T_kw)
fprintf('\n================ Test 35 preliminary summary ================\n');
fprintf('Composition: purity MAE %.3f, mixed PR AUC %.3f, strong-mixed PR AUC %.3f.\n', ...
    T_comp.purity_MAE(1), T_comp.mixed_PR_AUC(1), T_comp.strong_mixed_PR_AUC(1));
[best_mape, idx] = min(T_q.MAPE_pct);
fprintf('Best held-out synthetic SWS MAPE: %s = %.2f%%.\n', ...
    T_q.model_name(idx), best_mape);
for b = unique(T_purity.purity_bin,'stable')'
    X = T_purity(T_purity.purity_bin == b,:);
    [v,j] = min(X.MAPE_pct);
    fprintf('  Best in %s: %s = %.2f%%.\n', b, X.model_name(j), v);
end
if ~isempty(T_kw)
    S = summarize_q_predictions(T_kw, "model_name");
    [v,j] = min(S.MAPE_pct);
    fprintf('kWave external best overall: %s = %.2f%%.\n', S.model_name(j), v);
end
fprintf('==============================================================\n');
end

function q = invert_mapping_to_q(mapping, k_target)
k = double(mapping.k_cent(:));
E = double(mapping.Ecum(:));
valid = isfinite(k) & isfinite(E);
k = k(valid); E = E(valid);
if numel(k) < 2 || ~isfinite(k_target)
    q = NaN; return;
end
[kuniq, ia] = unique(k, 'stable');
Euniq = E(ia);
if k_target <= min(kuniq)
    q = Euniq(1);
elseif k_target >= max(kuniq)
    q = Euniq(end);
else
    q = interp1(kuniq, Euniq, k_target, 'linear');
end
q = clamp01(q);
end

function q = theory_q_for_condition(cfg, feat_cfg, regime)
field_type = regime_to_theory_field(regime);
try
    if field_type == "Partial3D"
        q = 0.5 * (theory_one(cfg, feat_cfg, "Diffuse2D") + ...
            theory_one(cfg, feat_cfg, "Diffuse3D"));
    else
        q = theory_one(cfg, feat_cfg, field_type);
    end
catch
    q = 0.5;
end
q = clamp01(q);
end

function field_type = regime_to_theory_field(regime)
switch string(regime)
    case "directional_2D"
        field_type = "SingleWave";
    case "diffuse_2D"
        field_type = "Diffuse2D";
    case "partial_3D"
        field_type = "Partial3D";
    otherwise
        field_type = "Diffuse3D";
end
end

function q = theory_one(cfg, feat_cfg, field_type)
out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    cfg.dx, cfg.dz, cfg.f0, feat_cfg.cs_guess_used, ...
    'M', feat_cfg.M, 'Gamma', feat_cfg.gamma_win, ...
    'PadFactor', feat_cfg.pad_factor, 'Nbins', 'auto', ...
    'SmoothSigma', 1, 'TheoryMode', 'S2D', ...
    'FieldType', field_type, 'Plot', false);
q = out.q_th;
end

function p = dominant_material_fraction(cs_patch)
v = round(double(cs_patch(:))*1e6)/1e6;
v = v(isfinite(v));
if isempty(v), p = NaN; return; end
[u,~,ic] = unique(v);
counts = accumarray(ic, 1, [numel(u),1]);
p = max(counts) / numel(v);
end

function y = q_to_sws(mappings, q, f0)
n = numel(q); y = nan(n,1);
for i = 1:n
    if isscalar(f0), f = f0; else, f = f0(i); end
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), f);
end
end

function y = infer_kwave_true_sws(T, C)
x = T.x_center_m;
z = T.z_center_m;
if ismember('x', T.Properties.VariableNames), x = T.x; end
if ismember('z', T.Properties.VariableNames), z = T.z; end
center = [0.025 0.025]; radius = 0.010;
r = hypot(x - center(1), z - center(2));
y = 2 * ones(height(T),1);
y(r <= radius) = 3;
end

function label = canonical_kwave_regime(case_name)
switch string(case_name)
    case "Directional 2D", label = "directional_2D";
    case "Diffuse 2D", label = "diffuse_2D";
    case "Projected diffuse 3D", label = "diffuse_3D";
    otherwise, label = "partial_3D";
end
end

function bin = purity_bin(p)
bin = repmat("strongly_mixed", size(p));
bin(p >= 0.75) = "moderately_mixed";
bin(p >= 0.95) = "near_pure";
bin(p >= 0.99) = "pure";
bin(~isfinite(p)) = "unknown";
end

function auc = auc_binary(y, score, mode)
y = logical(y); score = double(score);
valid = isfinite(score) & ~isnan(y);
y = y(valid); score = score(valid);
if numel(unique(y)) < 2
    auc = NaN; return;
end
[score, ord] = sort(score, 'descend');
y = y(ord);
tp = cumsum(y); fp = cumsum(~y);
P = sum(y); N = sum(~y);
if mode == "roc"
    x = [0; fp/N; 1]; z = [0; tp/P; 1];
else
    x = [0; tp/P]; z = [1; tp ./ max(tp+fp,1)];
end
auc = trapz(x, z);
end

function s = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
s = score(:, idx);
end

function y = clamp01(x)
y = min(max(x, 0), 1);
end

function e = rmse(a,b)
e = sqrt(mean((a-b).^2,'omitnan'));
end

function [Z,nz,nx] = vector_to_grid(v, iz, ix, nz, nx)
if nargin < 4
    nz = max(iz); nx = max(ix);
end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], iz, ix)) = v;
end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars)
    drop(i) = iscell(T.(vars(i)));
end
T(:, cellstr(vars(drop))) = [];
end

function assert_policy_no_forbidden_predictors(T)
vars = string(T.Properties.VariableNames);
assert(~any(vars == "confidence"), ...
    'Unexpected confidence variable found in Test 35 dataset.');
end

function sig = cache_signature(CFG)
sig = string(jsonencode(struct('version',CFG.Version,'mode',CFG.RunMode, ...
    'freq',CFG.Frequencies,'M',CFG.M,'regimes',CFG.Regimes, ...
    'geometries',CFG.Geometries,'dx',CFG.Dx, ...
    'max_patches',CFG.MaxPatchesPerCondition,'target_step',CFG.TargetStepM)));
end

function write_config_json(CFG, path)
fid = fopen(path,'w');
assert(fid > 0, 'Could not write %s', path);
fprintf(fid, '%s', jsonencode(CFG, PrettyPrint=true));
fclose(fid);
end

function tf = env_true(name, default_value)
v = strtrim(lower(string(getenv(name))));
if v == "", tf = logical(default_value); return; end
tf = ismember(v, ["1","true","yes","y","on"]);
end

function x = env_double(name, default_value)
v = strtrim(lower(string(getenv(name))));
if v == ""
    x = default_value;
else
    x = str2double(v);
    if isnan(x)
        error('%s must be numeric or inf. Received: %s', name, v);
    end
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
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
