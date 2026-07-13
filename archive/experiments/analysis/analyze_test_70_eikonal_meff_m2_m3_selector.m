%% analyze_test_70_eikonal_meff_m2_m3_selector.m
% Test 70: high-reliability M-effective M2/M3 selector on dense Eikonal maps.
%
% Purpose:
%   Reduce hard-core underestimation in dense Test66 Eikonal/readout maps
%   without retraining q, composition, or risk models. This test compares
%   fixed M=2, fixed M=3, and a conservative operational M2->M3 switch.
%
% Operational decision rule:
%   Switch from M=2 to M=3 only when the M=2 patch is high reliability,
%   predicted pure/non-mixed, hard-like, and M=3 increases SWS smoothly.
%   True SWS, true patch purity, ROI labels, and distance-to-interface are
%   used only for evaluation and plotting.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST70_MODE          = quick | full
%   ADAPTIVE_REQ_TEST70_TEST66_SOURCE = full_a | /path/to/test66/run
%   ADAPTIVE_REQ_TEST70_RISK_BUNDLE   = /path/to/test68_estimability_mask_compact.mat
%   ADAPTIVE_REQ_TEST70_MODEL         = q_spectrum_plus_composition
%   ADAPTIVE_REQ_TEST70_RISK_THRESHOLD = high-reliability max risk, default 0.20
%   ADAPTIVE_REQ_TEST70_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST70_MAX_MAPS      = max representative maps
%   ADAPTIVE_REQ_TEST70_GEOMETRIES    = optional comma-separated filter
%   ADAPTIVE_REQ_TEST70_FREQUENCIES   = optional comma-separated filter
%   ADAPTIVE_REQ_TEST70_REALISM_LEVELS = optional comma-separated filter
%   ADAPTIVE_REQ_TEST70_FIELD_REGIMES = optional comma-separated filter

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test70_configuration.json'));

fprintf('\nTest 70: high-reliability M-eff M2/M3 selector\n');
fprintf('Source Test66 run: %s\n', CFG.Test66RunDir);
fprintf('Base q model: %s | risk threshold %.2f | save maps: %d\n', ...
    CFG.BaseModel, CFG.RiskThreshold, CFG.SaveAllMaps);
fprintf('No q, composition, risk, or correction model is retrained.\n');

T = load_and_prepare_test66(CFG);
assert(any(abs(T.M-2)<1e-8), 'Test70 requires M=2 rows in Test66 source.');
assert(any(abs(T.M-3)<1e-8), ['Test70 requires M=3 rows in Test66 source. ' ...
    'Generate a filtered Test66 run with ADAPTIVE_REQ_TEST66_M_VALUES=2,3.']);

W = build_m2_m3_wide(T, CFG);
assert(~isempty(W), 'No aligned M=2/M=3 rows found.');

fprintf('Loading frozen Test68 risk ensemble...\n');
S = load(CFG.RiskBundle, 'MODELS');
RISK_MODELS = select_test68_models(S.MODELS, CFG);
W = attach_m2_risk(W, RISK_MODELS);

R = apply_selector_strategies(W, CFG);
R_save = remove_cell_columns(R);
writetable(R_save, fullfile(OUT.table_dir, 'test70_patch_level_strategy_results.csv'));

SUM = write_summaries(R, OUT);
plot_summary_figures(R, SUM, OUT, CFG);
if CFG.SaveAllMaps
    plot_representative_maps(R, OUT, CFG);
end
save(fullfile(OUT.data_dir, 'test70_compact_results.mat'), 'CFG','SUM','-v7.3');
write_readme(R, SUM, OUT, CFG);
print_console_summary(SUM.overall, SUM.by_roi);

fprintf('\nTables: %s\nFigures: %s\nTest 70 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST70_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST70_MODE must be quick or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST70_MODEL', "q_spectrum_plus_composition");
src = env_string('ADAPTIVE_REQ_TEST70_TEST66_SOURCE', "full_a");
if contains(src, filesep) || startsWith(src, "/")
    CFG.Test66RunDir = char(src);
else
    CFG.Test66RunDir = fullfile(root_dir, 'outputs', ...
        'test_66_eikonal_realistic_transfer_validation', char(src));
end
CFG.Test66PatchCsv = fullfile(CFG.Test66RunDir, 'tables', 'test66_patch_level_results.csv');
bundle = env_string('ADAPTIVE_REQ_TEST70_RISK_BUNDLE', "");
if bundle == ""
    bundle = fullfile(root_dir, 'outputs', 'test_68_test66_estimability_mask', ...
        'data', 'test68_estimability_mask_compact.mat');
end
CFG.RiskBundle = char(bundle);

CFG.CsGuess = env_number('ADAPTIVE_REQ_TEST70_CS_GUESS', 3.0);
CFG.RiskThreshold = env_number('ADAPTIVE_REQ_TEST70_RISK_THRESHOLD', 0.20);
CFG.PurityThreshold = env_number('ADAPTIVE_REQ_TEST70_PURITY_THRESHOLD', 0.90);
CFG.MixedThreshold = env_number('ADAPTIVE_REQ_TEST70_MIXED_THRESHOLD', 0.25);
CFG.HardSwsRatio = env_number('ADAPTIVE_REQ_TEST70_HARD_SWS_RATIO', 1.05);
CFG.MaxM3RelDiff = env_number('ADAPTIVE_REQ_TEST70_MAX_M3_REL_DIFF', 0.25);
CFG.MaxAlignToleranceM = env_number('ADAPTIVE_REQ_TEST70_ALIGN_TOLERANCE_M', 0.30e-3);
CFG.M2MeffMaxForSwitch = env_number('ADAPTIVE_REQ_TEST70_M2_MEFF_MAX_FOR_SWITCH', 1.90);
CFG.MinM3GainPct = env_number('ADAPTIVE_REQ_TEST70_MIN_M3_GAIN_PCT', 0.0);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST70_SAVE_ALL_MAPS', true);
CFG.MaxMaps = env_number('ADAPTIVE_REQ_TEST70_MAX_MAPS', ternary(CFG.QuickMode, 12, 40));

if CFG.QuickMode
    CFG.Geometries = ["homogeneous_cs4","inclusion_2_4","bilayer_2_4","bilayer_inclusion_2_3_4"];
    CFG.Frequencies = [300 500];
    CFG.RealismLevels = ["clean","readout_medium"];
    CFG.FieldRegimes = ["single_source_lateral","diffuse_like_8src"];
else
    CFG.Geometries = strings(0,1);
    CFG.Frequencies = [];
    CFG.RealismLevels = strings(0,1);
    CFG.FieldRegimes = strings(0,1);
end
CFG.Geometries = env_string_list('ADAPTIVE_REQ_TEST70_GEOMETRIES', CFG.Geometries);
CFG.Frequencies = env_number_list('ADAPTIVE_REQ_TEST70_FREQUENCIES', CFG.Frequencies);
CFG.RealismLevels = env_string_list('ADAPTIVE_REQ_TEST70_REALISM_LEVELS', CFG.RealismLevels);
CFG.FieldRegimes = env_string_list('ADAPTIVE_REQ_TEST70_FIELD_REGIMES', CFG.FieldRegimes);
CFG.SourceSeeds = env_number_list('ADAPTIVE_REQ_TEST70_SOURCE_SEEDS', []);
CFG.NoiseSeeds = env_number_list('ADAPTIVE_REQ_TEST70_NOISE_SEEDS', []);
end

function OUT = make_output_dirs(root_dir, CFG)
base = fullfile(root_dir, 'outputs', 'test_70_eikonal_meff_m2_m3_selector');
OUT.root_dir = fullfile(base, char(CFG.Mode));
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
dirs = {OUT.root_dir, OUT.table_dir, OUT.figure_dir, OUT.map_dir, OUT.data_dir};
for i = 1:numel(dirs)
    if exist(dirs{i}, 'dir') ~= 7, mkdir(dirs{i}); end
end
end

%% Loading and alignment

function T = load_and_prepare_test66(CFG)
assert(exist(CFG.Test66PatchCsv,'file') == 2, 'Missing Test66 patch table: %s', CFG.Test66PatchCsv);
fprintf('Loading Test66 patch table...\n');
opts = detectImportOptions(CFG.Test66PatchCsv, 'TextType','string');
T = readtable(CFG.Test66PatchCsv, opts);
T = normalize_loaded_types(T);
T = T(T.model_name == CFG.BaseModel & ismember(T.M, [2 3]), :);
T.condition_base_key = strip_m_suffix(T.condition_key);

if ~isempty(CFG.Geometries), T = T(ismember(T.geometry, CFG.Geometries), :); end
if ~isempty(CFG.Frequencies), T = T(ismember(T.f0, CFG.Frequencies), :); end
if ~isempty(CFG.RealismLevels), T = T(ismember(T.realism_level, CFG.RealismLevels), :); end
if ~isempty(CFG.FieldRegimes), T = T(ismember(T.field_regime, CFG.FieldRegimes), :); end
if ~isempty(CFG.SourceSeeds), T = T(ismember(T.source_seed, CFG.SourceSeeds), :); end
if ~isempty(CFG.NoiseSeeds), T = T(ismember(T.noise_seed, CFG.NoiseSeeds), :); end
fprintf('Selected %d patch rows after filters.\n', height(T));

T = add_operational_features(T, CFG);
end

function T = normalize_loaded_types(T)
names = string(T.Properties.VariableNames);
for v = ["condition_key","geometry","geometry_family","case_id","case_family", ...
        "realism_level","field_regime","field_regime_ood","model_name", ...
        "purity_bin","distance_bin","distance_over_window_bin","roi_region", ...
        "analysis_region","snr_bin","amplitude_bin","depth_bin"]
    if ismember(v, names)
        T.(v) = string(T.(v));
    end
end
end

function T = add_operational_features(T, CFG)
T.SWS_pred = first_existing_numeric(T, ["SWS_pred","sws_pred"]);
T.SWS_pred = clamp(T.SWS_pred, 0.2, 20);
T.q_pred = clamp(T.q_pred, 0, 1);
if ~ismember("k_pred", string(T.Properties.VariableNames))
    T.k_pred = 2*pi*T.f0 ./ T.SWS_pred;
else
    bad = ~isfinite(T.k_pred) | T.k_pred <= 0;
    T.k_pred(bad) = 2*pi*T.f0(bad) ./ T.SWS_pred(bad);
end
T.log_sws_pred = log(T.SWS_pred);
T.lambda_pred_m = T.SWS_pred ./ T.f0;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.M_eff_pred = T.window_length_m ./ max(T.lambda_pred_m, eps);
T.q_theory_prior = ensure_numeric_var(T, "q_theory_prior", nan(height(T),1));
T.sws_theory = ensure_numeric_var(T, "sws_theory", nan(height(T),1));
T.abs_q_minus_theory = abs(T.q_pred - T.q_theory_prior);
T.rel_sws_minus_theory = abs(T.SWS_pred - T.sws_theory) ./ max(abs(T.SWS_pred), eps);
T.predicted_patch_purity = clamp(ensure_numeric_var(T, "predicted_patch_purity", ones(height(T),1)), 0, 1);
T.p_mixed = clamp(ensure_numeric_var(T, "p_mixed", zeros(height(T),1)), 0, 1);
T.p_strong_mixed = clamp(ensure_numeric_var(T, "p_strong_mixed", zeros(height(T),1)), 0, 1);
T.tracking_snr_db = ensure_numeric_var(T, "tracking_snr_db", ensure_numeric_var(T, "snr_proxy_db", nan(height(T),1)));
T.acoustic_readout_amplitude_factor = ensure_numeric_var(T, "acoustic_readout_amplitude_factor", nan(height(T),1));
T.local_amplitude = ensure_numeric_var(T, "local_amplitude", nan(height(T),1));
T.source_to_patch_distance_m = ensure_numeric_var(T, "source_to_patch_distance_m", nan(height(T),1));
T.depth_m = ensure_numeric_var(T, "depth_m", T.z_center_m);
T = add_map_neighborhood_features(T);
end

function T = add_map_neighborhood_features(T)
T.grad_q_pred = nan(height(T),1);
T.grad_sws_pred = nan(height(T),1);
T.local_std_q_pred = nan(height(T),1);
T.local_std_sws_pred = nan(height(T),1);
conds = unique(T.condition_key, 'stable');
for ci = 1:numel(conds)
    % Neighborhood features are only used by the frozen risk model. Compute
    % them separately for each condition/M so maps with different valid
    % borders do not contaminate each other.
    mids = unique(T.M(T.condition_key == conds(ci)));
    for mi = 1:numel(mids)
        idx = find(T.condition_key == conds(ci) & T.M == mids(mi));
        [Q, rows, cols] = vector_to_map(T.map_iz(idx), T.map_ix(idx), T.q_pred(idx));
        S = vector_to_map(T.map_iz(idx), T.map_ix(idx), T.SWS_pred(idx));
        Gq = gradient_magnitude(Q);
        Gs = gradient_magnitude(S);
        Lq = local_std_nan(Q, 3);
        Ls = local_std_nan(S, 3);
        lin = sub2ind(size(Q), rows, cols);
        T.grad_q_pred(idx) = Gq(lin);
        T.grad_sws_pred(idx) = Gs(lin);
        T.local_std_q_pred(idx) = Lq(lin);
        T.local_std_sws_pred(idx) = Ls(lin);
    end
end
end

function W = build_m2_m3_wide(T, CFG)
T2 = T(abs(T.M-2)<1e-8, :);
T3 = T(abs(T.M-3)<1e-8, :);
conds = intersect(unique(T2.condition_base_key, 'stable'), unique(T3.condition_base_key, 'stable'), 'stable');
parts = cell(numel(conds),1);
tol_um = round(1e6*CFG.MaxAlignToleranceM);
for ci = 1:numel(conds)
    A = sortrows(T2(T2.condition_base_key == conds(ci), :), {'x_center_m','z_center_m'});
    B = sortrows(T3(T3.condition_base_key == conds(ci), :), {'x_center_m','z_center_m'});
    if isempty(A) || isempty(B), continue; end
    A.x_um = round(1e6*A.x_center_m);
    A.z_um = round(1e6*A.z_center_m);
    B.x_um = round(1e6*B.x_center_m);
    B.z_um = round(1e6*B.z_center_m);
    keyA = string(A.x_um) + "_" + string(A.z_um);
    keyB = string(B.x_um) + "_" + string(B.z_um);
    [tf, loc] = ismember(keyB, keyA);
    if ~all(tf)
        [loc, dist_um] = nearest_center(A.x_um, A.z_um, B.x_um, B.z_um);
        tf = dist_um <= tol_um;
    end
    if any(tf)
        parts{ci} = assemble_wide(A(loc(tf),:), B(tf,:));
    end
end
parts = parts(~cellfun(@isempty, parts));
W = vertcat(parts{:});
fprintf('Aligned %d M=3 centers to M=2 centers across %d conditions.\n', height(W), numel(parts));
end

function [loc, dist_um] = nearest_center(xA, zA, xB, zB)
loc = zeros(numel(xB),1);
dist_um = inf(numel(xB),1);
for i = 1:numel(xB)
    d = hypot(double(xA) - double(xB(i)), double(zA) - double(zB(i)));
    [dist_um(i), loc(i)] = min(d);
end
end

function W = assemble_wide(T2, T3)
meta = ["condition_key","geometry","geometry_family","case_id","case_family", ...
    "condition_base_key", ...
    "realism_level","field_regime","field_regime_ood","source_seed","noise_seed", ...
    "f0","dx","dz","map_iz","map_ix","x_center_m","z_center_m","x_um","z_um", ...
    "true_SWS","k_true","patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr", ...
    "distance_to_interface_m","distance_to_interface_mm", ...
    "distance_to_interface_over_window_radius","purity_bin","distance_bin", ...
    "distance_over_window_bin","depth_bin","snr_bin","amplitude_bin","roi_region", ...
    "analysis_region","tracking_snr_db","local_amplitude","depth_m", ...
    "shear_attenuation_factor","acoustic_readout_amplitude_factor"];
meta = intersect(meta, string(T3.Properties.VariableNames), 'stable');
W = T3(:, cellstr(meta));
W.condition_key_M2 = T2.condition_key;
W.condition_key_M3 = T3.condition_key;
W.match_distance_m = hypot(T2.x_center_m - T3.x_center_m, T2.z_center_m - T3.z_center_m);

W.M2_sws = T2.SWS_pred;
W.M2_q = T2.q_pred;
W.M2_k = T2.k_pred;
W.M2_M_eff = T2.M_eff_pred;
W.M2_predicted_patch_purity = T2.predicted_patch_purity;
W.M2_p_mixed = T2.p_mixed;
W.M2_p_strong_mixed = T2.p_strong_mixed;
W.M2_grad_q_pred = T2.grad_q_pred;
W.M2_grad_sws_pred = T2.grad_sws_pred;
W.M2_local_std_q_pred = T2.local_std_q_pred;
W.M2_local_std_sws_pred = T2.local_std_sws_pred;
W.M2_tracking_snr_db = T2.tracking_snr_db;
W.M2_local_amplitude = T2.local_amplitude;

W.M3_sws = T3.SWS_pred;
W.M3_q = T3.q_pred;
W.M3_k = T3.k_pred;
W.M3_M_eff = T3.M_eff_pred;
W.M3_predicted_patch_purity = T3.predicted_patch_purity;
W.M3_p_mixed = T3.p_mixed;
W.M3_p_strong_mixed = T3.p_strong_mixed;
W.sws_M3_minus_M2 = W.M3_sws - W.M2_sws;
W.rel_sws_M3_minus_M2 = W.sws_M3_minus_M2 ./ max(W.M2_sws, eps);
end

%% Risk and strategies

function RISK_MODELS = select_test68_models(MODELS, CFG)
RISK_MODELS = struct('split_name', {}, 'fold_name', {}, 'detectors', {});
for i = 1:numel(MODELS)
    if ~isfield(MODELS(i), 'detectors'), continue; end
    D = MODELS(i).detectors;
    if ~isfield(D, 'high20') || ~isfield(D.high20, 'bagged_trees'), continue; end
    if CFG.QuickMode || string(MODELS(i).split_name) == "grouped_condition_repeated"
        RISK_MODELS(end+1).split_name = string(MODELS(i).split_name); %#ok<AGROW>
        RISK_MODELS(end).fold_name = string(MODELS(i).fold_name);
        RISK_MODELS(end).detectors = D;
    end
end
assert(~isempty(RISK_MODELS), 'No Test68 high20 bagged-tree risk models found.');
end

function W = attach_m2_risk(W, RISK_MODELS)
T = table();
T.condition_key = W.condition_key;
T.field_regime = W.field_regime;
T.q_pred = W.M2_q;
T.SWS_pred = W.M2_sws;
T.sws_pred = W.M2_sws;
T.k_pred = W.M2_k;
T.f0 = W.f0;
T.M = 2*ones(height(W),1);
T.dx = W.dx;
T.dz = W.dz;
T.REQ_StepX = nan(height(W),1);
T.REQ_StepZ = nan(height(W),1);
T.TargetStepM = nan(height(W),1);
T.predicted_patch_purity = W.M2_predicted_patch_purity;
T.p_mixed = W.M2_p_mixed;
T.p_strong_mixed = W.M2_p_strong_mixed;
T.tracking_snr_db = W.M2_tracking_snr_db;
T.local_amplitude = W.M2_local_amplitude;
T.acoustic_readout_amplitude_factor = ensure_wide_numeric(W, "acoustic_readout_amplitude_factor", nan(height(W),1));
T.source_to_patch_distance_m = nan(height(W),1);
T.depth_m = ensure_wide_numeric(W, "depth_m", W.z_center_m);
T.q_theory_prior = nan(height(W),1);
T.sws_theory = nan(height(W),1);
T.grad_q_pred = W.M2_grad_q_pred;
T.grad_sws_pred = W.M2_grad_sws_pred;
T.local_std_q_pred = W.M2_local_std_q_pred;
T.local_std_sws_pred = W.M2_local_std_sws_pred;
T.log_sws_pred = log(clamp(T.SWS_pred,0.2,20));
T.M_eff_pred = W.M2_M_eff;
T.abs_q_minus_theory = nan(height(W),1);
T.rel_sws_minus_theory = nan(height(W),1);

R = nan(height(W), numel(RISK_MODELS));
for i = 1:numel(RISK_MODELS)
    D = RISK_MODELS(i).detectors;
    X = design_matrix_for_inference(T, D.design);
    R(:,i) = predict_positive(D.high20.bagged_trees, X);
end
W.risk_M2_high20 = mean(R, 2, 'omitnan');
W.risk_M2_std = std(R, 0, 2, 'omitnan');
W.reliability_M2 = 1 - W.risk_M2_high20;
fprintf('Attached Test68 risk to %d aligned rows using %d detectors.\n', height(W), numel(RISK_MODELS));
end

function X = design_matrix_for_inference(T, DESIGN)
features = string(DESIGN.numeric_features);
Xnum = zeros(height(T), numel(features));
for i = 1:numel(features)
    v = features(i);
    if ismember(v, string(T.Properties.VariableNames))
        a = double(T.(v));
    else
        a = nan(height(T),1);
    end
    med = DESIGN.numeric_median(i);
    a(~isfinite(a)) = med;
    Xnum(:,i) = a;
end
mu = DESIGN.numeric_mu;
sd = DESIGN.numeric_sd;
sd(~isfinite(sd) | sd <= 0) = 1;
Xnum = (Xnum - mu) ./ sd;
Xcat = onehot(string(T.field_regime), string(DESIGN.field_regime_categories));
X = [Xnum Xcat];
end

function risk = predict_positive(model, X)
[~, score] = predict(model, X);
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
risk = clamp(score(:,idx), 0, 1);
end

function R = apply_selector_strategies(W, CFG)
hard_like = W.M2_sws >= CFG.HardSwsRatio * CFG.CsGuess | W.M2_M_eff <= CFG.M2MeffMaxForSwitch;
high_rel = W.risk_M2_high20 < CFG.RiskThreshold;
pred_pure = W.M2_predicted_patch_purity >= CFG.PurityThreshold & W.M2_p_mixed <= CFG.MixedThreshold;
m3_increases = W.M3_sws >= W.M2_sws .* (1 + CFG.MinM3GainPct/100);
m3_stable = abs(W.rel_sws_M3_minus_M2) <= CFG.MaxM3RelDiff;
switch_mask = high_rel & pred_pure & hard_like & m3_increases & m3_stable;

R = [
    make_strategy(W, "fixed_M2", W.M2_sws, W.M2_q, 2*ones(height(W),1), false(height(W),1), "none");
    make_strategy(W, "fixed_M3", W.M3_sws, W.M3_q, 3*ones(height(W),1), true(height(W),1), "all_M3");
    make_strategy(W, "reliable_hard_M2_to_M3_switch", choose(W.M2_sws,W.M3_sws,switch_mask), ...
        choose(W.M2_q,W.M3_q,switch_mask), choose(2*ones(height(W),1),3*ones(height(W),1),switch_mask), ...
        switch_mask, "reliable_hard_M3");
    make_oracle(W)
    ];

fprintf('Rule switched %.2f%% of aligned pixels to M=3.\n', 100*mean(switch_mask,'omitnan'));
end

function R = make_strategy(W, name, sws, q, source_M, was_corrected, mechanism)
meta = ["condition_key","geometry","geometry_family","case_id","case_family", ...
    "realism_level","field_regime","field_regime_ood","source_seed","noise_seed", ...
    "f0","dx","dz","map_iz","map_ix","x_center_m","z_center_m","x_um","z_um", ...
    "true_SWS","k_true","patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr", ...
    "distance_to_interface_m","distance_to_interface_mm", ...
    "distance_to_interface_over_window_radius","purity_bin","distance_bin", ...
    "distance_over_window_bin","depth_bin","snr_bin","amplitude_bin","roi_region", ...
    "analysis_region","tracking_snr_db","local_amplitude","depth_m", ...
    "shear_attenuation_factor","acoustic_readout_amplitude_factor", ...
    "M2_sws","M3_sws","M2_q","M3_q","M2_M_eff","M3_M_eff", ...
    "M2_predicted_patch_purity","M2_p_mixed","M2_p_strong_mixed", ...
    "risk_M2_high20","risk_M2_std","reliability_M2","rel_sws_M3_minus_M2"];
meta = intersect(meta, string(W.Properties.VariableNames), 'stable');
R = W(:, cellstr(meta));
R.strategy_name = repmat(string(name), height(W), 1);
R.SWS_pred_strategy = sws;
R.q_pred_strategy = q;
R.source_M = source_M;
R.was_corrected = was_corrected;
if isscalar(mechanism), mechanism = repmat(string(mechanism), height(W), 1); end
R.correction_mechanism = string(mechanism);
R.correction_pct = 100*(sws - W.M2_sws)./max(W.M2_sws, eps);
R.sws_signed_error_pct = 100*(sws - W.true_SWS)./W.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.base_abs_error_pct = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
R.improved_vs_M2 = R.sws_abs_error_pct + 0.5 < R.base_abs_error_pct;
R.harmed_vs_M2 = R.sws_abs_error_pct > R.base_abs_error_pct + 0.5;
R.diagnostic_only = repmat(name == "oracle_best_M2_M3", height(W), 1);
end

function R = make_oracle(W)
err2 = abs(100*(W.M2_sws - W.true_SWS)./W.true_SWS);
err3 = abs(100*(W.M3_sws - W.true_SWS)./W.true_SWS);
use3 = err3 < err2;
sws = choose(W.M2_sws, W.M3_sws, use3);
q = choose(W.M2_q, W.M3_q, use3);
R = make_strategy(W, "oracle_best_M2_M3", sws, q, choose(2*ones(height(W),1),3*ones(height(W),1),use3), ...
    use3, "oracle_best");
end

%% Summaries

function SUM = write_summaries(R, OUT)
SUM = struct();
SUM.overall = summarize_group(R, ["strategy_name"]);
SUM.by_geometry = summarize_group(R, ["strategy_name","geometry"]);
SUM.by_frequency = summarize_group(R, ["strategy_name","f0"]);
SUM.by_geometry_frequency = summarize_group(R, ["strategy_name","geometry","f0"]);
SUM.by_realism = summarize_group(R, ["strategy_name","realism_level"]);
SUM.by_field_regime = summarize_group(R, ["strategy_name","field_regime"]);
SUM.by_roi = summarize_group(R, ["strategy_name","roi_region"]);
SUM.by_true_sws = summarize_group(R, ["strategy_name","true_SWS"]);
SUM.by_distance = summarize_group(R, ["strategy_name","distance_over_window_bin"]);

writetable(SUM.overall, fullfile(OUT.table_dir,'test70_summary_overall.csv'));
writetable(SUM.by_geometry, fullfile(OUT.table_dir,'test70_summary_by_geometry.csv'));
writetable(SUM.by_frequency, fullfile(OUT.table_dir,'test70_summary_by_frequency.csv'));
writetable(SUM.by_geometry_frequency, fullfile(OUT.table_dir,'test70_summary_by_geometry_frequency.csv'));
writetable(SUM.by_realism, fullfile(OUT.table_dir,'test70_summary_by_realism.csv'));
writetable(SUM.by_field_regime, fullfile(OUT.table_dir,'test70_summary_by_field_regime.csv'));
writetable(SUM.by_roi, fullfile(OUT.table_dir,'test70_summary_by_roi.csv'));
writetable(SUM.by_true_sws, fullfile(OUT.table_dir,'test70_summary_by_true_sws.csv'));
writetable(SUM.by_distance, fullfile(OUT.table_dir,'test70_summary_by_distance_over_window.csv'));
end

function S = summarize_group(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, keys] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    corrected = X.was_corrected;
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), 100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.sws_signed_error_pct > 0,'omitnan'), 100*mean(corrected,'omitnan'), ...
        100*mean(X.improved_vs_M2 & corrected,'omitnan'), 100*mean(X.harmed_vs_M2 & corrected,'omitnan'), ...
        mean(X.risk_M2_high20,'omitnan'), mean(X.M2_M_eff,'omitnan'), mean(X.M3_M_eff,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_APE_pct','bias_pct','median_bias_pct', ...
        'high10_pct','high20_pct','underestimate_pct','overestimate_pct', ...
        'corrected_pct','improved_corrected_pct','harmed_corrected_pct', ...
        'mean_risk_M2','mean_M2_Meff','mean_M3_Meff'});
end
S = [keys vertcat(rows{:})];
end

%% Figures

function plot_summary_figures(R, SUM, OUT, CFG)
fig = figure('Color','w','Position',[80 80 1300 750]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl, 'Test70: M-eff hard-core M2/M3 selector', 'FontWeight','bold');
plot_strategy_bar(nexttile, SUM.overall, 'MAPE_pct', 'Overall MAPE (%)');
plot_strategy_bar(nexttile, SUM.overall, 'bias_pct', 'Overall bias (%)');
plot_strategy_bar(nexttile, SUM.overall, 'high20_pct', 'High error >20% (%)');
H = SUM.by_roi(ismember(SUM.by_roi.roi_region, ["hard_core","hard_inclusion_core","soft_core","background_far","interface_0_1mm"]),:);
plot_grouped_lines(nexttile, H, "roi_region", "MAPE_pct", "ROI", "MAPE (%)");
F = SUM.by_frequency;
plot_grouped_lines(nexttile, F, "f0", "MAPE_pct", "Frequency (Hz)", "MAPE (%)");
plot_grouped_lines(nexttile, SUM.by_geometry, "geometry", "MAPE_pct", "Geometry", "MAPE (%)");
export_fig(fig, fullfile(OUT.figure_dir,'test70_strategy_summary.png'));

fig = figure('Color','w','Position',[100 100 1100 500]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(tl, 'Hard-core bias and correction behavior', 'FontWeight','bold');
HR = SUM.by_roi(ismember(SUM.by_roi.roi_region, ["hard_core","hard_inclusion_core"]),:);
plot_grouped_lines(nexttile, HR, "roi_region", "bias_pct", "ROI", "Mean signed error (%)");
plot_grouped_lines(nexttile, HR, "roi_region", "corrected_pct", "ROI", "Pixels switched/corrected (%)");
export_fig(fig, fullfile(OUT.figure_dir,'test70_hard_core_bias_and_switching.png'));
end

function plot_strategy_bar(ax, S, metric, ylab)
strategies = ["fixed_M2","fixed_M3","reliable_hard_M2_to_M3_switch","oracle_best_M2_M3"];
S = S(ismember(S.strategy_name, strategies),:);
[~,ord] = ismember(S.strategy_name, strategies);
S = S(sortidx(ord),:);
bar(ax, categorical(pretty_strategy(S.strategy_name), pretty_strategy(strategies)), S.(metric));
ylabel(ax, ylab); grid(ax,'on'); xtickangle(ax,25);
end

function plot_grouped_lines(ax, S, xvar, metric, xlab, ylab)
strategies = ["fixed_M2","fixed_M3","reliable_hard_M2_to_M3_switch","oracle_best_M2_M3"];
hold(ax,'on');
for si = 1:numel(strategies)
    X = S(S.strategy_name == strategies(si), :);
    if isempty(X), continue; end
    if isnumeric(X.(xvar))
        [x,ord] = sort(X.(xvar));
        y = X.(metric)(ord);
        plot(ax, x, y, '-o', 'DisplayName', pretty_strategy(strategies(si)));
    else
        cats = unique(string(X.(xvar)), 'stable');
        y = nan(numel(cats),1);
        for ci = 1:numel(cats)
            y(ci) = mean(X.(metric)(string(X.(xvar))==cats(ci)), 'omitnan');
        end
        plot(ax, 1:numel(cats), y, '-o', 'DisplayName', pretty_strategy(strategies(si)));
        xticks(ax, 1:numel(cats)); xticklabels(ax, pretty_category(cats)); xtickangle(ax,30);
    end
end
ylabel(ax, ylab); xlabel(ax, xlab); grid(ax,'on');
legend(ax,'Location','bestoutside');
end

function plot_representative_maps(R, OUT, CFG)
X = R(R.strategy_name=="reliable_hard_M2_to_M3_switch", :);
if isempty(X), return; end
[G, groups] = findgroups(X(:, {'condition_key','geometry','realism_level','field_regime','f0'}));
case_mape = splitapply(@(x) mean(x,'omitnan'), X.sws_abs_error_pct, G);
case_hard = splitapply(@(x) mean(x,'omitnan'), X.was_corrected, G);
CASE = [groups table(case_mape, case_hard, 'VariableNames', {'MAPE_pct','switch_fraction'})];
CASE = sortrows(CASE, {'switch_fraction','MAPE_pct'}, {'descend','descend'});
CASE = CASE(1:min(height(CASE), CFG.MaxMaps), :);
for i = 1:height(CASE)
    idx = X.condition_key == CASE.condition_key(i) & X.f0 == CASE.f0(i);
    Xi = X(idx,:);
    out_dir = fullfile(OUT.map_dir, sanitize(CASE.geometry(i)), sanitize(CASE.realism_level(i)), sanitize(CASE.field_regime(i)));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    plot_one_map(Xi, fullfile(out_dir, "test70_map__" + sanitize(CASE.condition_key(i)) + ".png"));
end
end

function plot_one_map(X, path)
fig = figure('Color','w','Position',[80 80 1450 920]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s | f=%g Hz | %s | %s', X.geometry(1), X.f0(1), ...
    X.realism_level(1), X.field_regime(1)), 'Interpreter','none', 'FontWeight','bold');
vars = ["true_SWS","M2_sws","M3_sws","SWS_pred_strategy", ...
    "risk_M2_high20","was_corrected","sws_signed_error_pct","base_abs_error_pct", ...
    "sws_abs_error_pct","M2_M_eff","M3_M_eff","patch_purity"];
titles = ["True SWS","M=2 SWS","M=3 SWS","Selected SWS", ...
    "M=2 high-error risk","Switch mask","Selected signed error","M=2 abs error", ...
    "Selected abs error","M=2 M_{eff}","M=3 M_{eff}","True patch purity"];
labels = ["m/s","m/s","m/s","m/s", ...
    "probability","0/1","%","%", ...
    "%","wavelengths","wavelengths","fraction"];
for pi = 1:numel(vars)
    ax = nexttile(tl);
    [A, x, z] = vector_to_grid(X.x_center_m, X.z_center_m, double(X.(vars(pi))));
    imagesc(ax, 1e3*x, 1e3*z, A); axis(ax,'image'); set(ax,'YDir','normal');
    title(ax, titles(pi), 'FontWeight','normal');
    xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
    cb = colorbar(ax); ylabel(cb, labels(pi));
end
export_fig(fig, path);
end

%% README and console

function write_readme(R, SUM, OUT, CFG)
lines = strings(0,1);
lines(end+1) = "# Test70 M-eff M2/M3 Selector";
lines(end+1) = "";
lines(end+1) = sprintf("- Source Test66 run: `%s`", CFG.Test66RunDir);
lines(end+1) = sprintf("- Base model: `%s`", CFG.BaseModel);
lines(end+1) = sprintf("- Switch risk threshold: %.2f", CFG.RiskThreshold);
lines(end+1) = "- No q, composition, risk, or correction model was retrained.";
lines(end+1) = "- Operational switch features: predicted SWS/q, predicted purity/mixedness, Test68 risk, M_eff proxies, and M2/M3 disagreement.";
lines(end+1) = "- Evaluation-only variables: true SWS, true patch purity, ROI, distance-to-interface, and errors.";
lines(end+1) = "";
lines(end+1) = "## Overall";
S = SUM.overall;
for i = 1:height(S)
    lines(end+1) = sprintf("- %s: MAPE %.2f%%, bias %.2f%%, high20 %.2f%%, switched %.2f%%.", ...
        pretty_strategy(S.strategy_name(i)), S.MAPE_pct(i), S.bias_pct(i), ...
        S.high20_pct(i), S.corrected_pct(i));
end
lines(end+1) = "";
lines(end+1) = "## Hard-Core";
H = SUM.by_roi(ismember(SUM.by_roi.roi_region, ["hard_core","hard_inclusion_core"]),:);
for i = 1:height(H)
    lines(end+1) = sprintf("- %s / %s: MAPE %.2f%%, bias %.2f%%, switched %.2f%%.", ...
        pretty_region(H.roi_region(i)), pretty_strategy(H.strategy_name(i)), ...
        H.MAPE_pct(i), H.bias_pct(i), H.corrected_pct(i));
end
lines(end+1) = "";
lines(end+1) = "## Outputs";
lines(end+1) = "- `tables/test70_patch_level_strategy_results.csv`";
lines(end+1) = "- `figures/test70_strategy_summary.png`";
lines(end+1) = "- `figures/test70_hard_core_bias_and_switching.png`";
lines(end+1) = "- `figures/maps_by_condition/`";
fid = fopen(fullfile(OUT.root_dir,'README_results.md'), 'w');
fprintf(fid, '%s\n', lines);
fclose(fid);
end

function print_console_summary(S, Sroi)
fprintf('\nInterpretive summary:\n');
for i = 1:height(S)
    fprintf('  %s: MAPE %.2f%% | bias %.2f%% | high20 %.2f%% | switched %.2f%%\n', ...
        pretty_strategy(S.strategy_name(i)), S.MAPE_pct(i), S.bias_pct(i), ...
        S.high20_pct(i), S.corrected_pct(i));
end
H = Sroi(ismember(Sroi.roi_region, ["hard_core","hard_inclusion_core"]),:);
if ~isempty(H)
    fprintf('  Hard-core rows written in test70_summary_by_roi.csv.\n');
end
end

%% Small utilities

function [A, rows, cols] = vector_to_map(rowv, colv, values)
ur = unique(rowv, 'stable');
uc = unique(colv, 'stable');
[~, rows] = ismember(rowv, ur);
[~, cols] = ismember(colv, uc);
A = nan(numel(ur), numel(uc));
lin = sub2ind(size(A), rows, cols);
A(lin) = values;
end

function [A, ux, uz] = vector_to_grid(x, z, values)
ux = unique(x, 'stable');
uz = unique(z, 'stable');
[~, cols] = ismember(x, ux);
[~, rows] = ismember(z, uz);
A = nan(numel(uz), numel(ux));
A(sub2ind(size(A), rows, cols)) = values;
end

function G = gradient_magnitude(A)
[gx, gz] = gradient(fillmissing2(A));
G = hypot(gx, gz);
G(~isfinite(A)) = NaN;
end

function L = local_std_nan(A, w)
pad = floor(w/2);
% Avoid depending on Image Processing Toolbox just for NaN padding.
Ap = nan(size(A,1) + 2*pad, size(A,2) + 2*pad);
Ap(1+pad:end-pad, 1+pad:end-pad) = A;
L = nan(size(A));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        B = Ap(i:i+w-1, j:j+w-1);
        L(i,j) = std(B(:), 'omitnan');
    end
end
end

function B = fillmissing2(A)
B = A;
med = median(A(:),'omitnan');
if ~isfinite(med), med = 0; end
B(~isfinite(B)) = med;
end

function X = onehot(x, cats)
X = zeros(numel(x), numel(cats));
for i = 1:numel(cats)
    X(:,i) = x == cats(i);
end
end

function v = ensure_numeric_var(T, name, default)
if ismember(name, string(T.Properties.VariableNames))
    v = double(T.(name));
else
    v = default;
end
if isscalar(v), v = repmat(v, height(T), 1); end
end

function v = ensure_wide_numeric(T, name, default)
if ismember(name, string(T.Properties.VariableNames))
    v = double(T.(name));
else
    v = default;
end
if isscalar(v), v = repmat(v, height(T), 1); end
end

function v = first_existing_numeric(T, names)
for n = string(names)
    if ismember(n, string(T.Properties.VariableNames))
        v = double(T.(n));
        return;
    end
end
error('None of the requested variables exists: %s', strjoin(string(names), ', '));
end

function y = choose(a, b, mask)
y = a;
y(mask) = b(mask);
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function out = sortidx(x)
[~, out] = sort(x);
end

function s = pretty_strategy(x)
s = string(x);
s(s=="fixed_M2") = "Fixed M=2";
s(s=="fixed_M3") = "Fixed M=3";
s(s=="reliable_hard_M2_to_M3_switch") = "Reliable hard M2→M3";
s(s=="oracle_best_M2_M3") = "Oracle best M2/M3";
end

function s = pretty_region(x)
s = string(x);
s = strrep(s, "_", " ");
s = regexprep(s, '\<hard\>', 'hard');
end

function s = pretty_category(x)
s = string(x);
s = strrep(s, "_", " ");
s = strrep(s, "cs", "c_s ");
s = strrep(s, "diffuse like 8src", "diffuse-like 8 sources");
s = strrep(s, "single source lateral", "single lateral source");
end

function T = remove_cell_columns(T)
names = string(T.Properties.VariableNames);
drop = false(size(names));
for i = 1:numel(names)
    drop(i) = iscell(T.(names(i)));
end
T(:, drop) = [];
end

function s = sanitize(x)
s = regexprep(string(x), '[^A-Za-z0-9_=-]+', '_');
s = char(s);
end

function y = strip_m_suffix(x)
% Test66 condition keys include the evaluated M and step. Remove that suffix
% so M=2 and M=3 predictions from the same physical field can be aligned.
y = regexprep(string(x), '__M[0-9p\.]+__step[0-9]+um$', '');
end

function export_fig(fig, path)
folder = fileparts(path);
if exist(folder,'dir') ~= 7, mkdir(folder); end
exportgraphics(fig, path, 'Resolution', 220);
close(fig);
end

function write_config_json(CFG, path)
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

function v = env_string(name, default)
raw = string(getenv(name));
if strlength(raw) == 0, v = string(default); else, v = raw; end
end

function v = env_number(name, default)
raw = string(getenv(name));
if strlength(raw) == 0
    v = default;
else
    v = str2double(raw);
    if ~isfinite(v), v = default; end
end
end

function v = env_true(name, default)
raw = lower(strtrim(string(getenv(name))));
if raw == "", v = default; else, v = ismember(raw, ["1","true","yes","on"]); end
end

function v = env_string_list(name, default)
raw = strtrim(string(getenv(name)));
if raw == "", v = string(default); return; end
parts = regexp(char(raw), '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
v = string(parts(:))';
end

function v = env_number_list(name, default)
raw = strtrim(string(getenv(name)));
if raw == "", v = default; return; end
parts = regexp(char(raw), '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
v = str2double(string(parts(:)))';
v = v(isfinite(v));
end

function y = ternary(tf, a, b)
if tf, y = a; else, y = b; end
end
