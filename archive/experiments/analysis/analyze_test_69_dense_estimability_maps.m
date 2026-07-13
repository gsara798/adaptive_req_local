%% analyze_test_69_dense_estimability_maps.m
% Test 69: apply the frozen Test68 estimability mask to dense Test66 maps.
%
% This script does not retrain anything. It loads the frozen Test68 bagged-tree
% high-error detectors, applies them to selected dense Test66 patch maps, and
% exports map panels showing SWS, error, risk, reliability, and masked SWS.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST69_MODE = quick | full
%   ADAPTIVE_REQ_TEST69_TEST66_SOURCE = full_a | quick | /path/to/Test66/run
%   ADAPTIVE_REQ_TEST69_TEST68_BUNDLE = /path/to/test68_estimability_mask_compact.mat
%   ADAPTIVE_REQ_TEST69_RUN_LABEL = optional output suffix, e.g. readout_medium
%   ADAPTIVE_REQ_TEST69_MODEL = q_spectrum_plus_composition
%   ADAPTIVE_REQ_TEST69_M = 2
%   ADAPTIVE_REQ_TEST69_RISK_THRESHOLD = 0.5
%   ADAPTIVE_REQ_TEST69_MAX_CONDITIONS = integer
%   ADAPTIVE_REQ_TEST69_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST69_GEOMETRIES = comma-separated geometry names
%   ADAPTIVE_REQ_TEST69_FREQUENCIES = comma-separated frequency list
%   ADAPTIVE_REQ_TEST69_REALISM_LEVELS = comma-separated realism names
%   ADAPTIVE_REQ_TEST69_FIELD_REGIMES = comma-separated field-regime names
%   ADAPTIVE_REQ_TEST69_SOURCE_SEEDS = comma-separated source seeds
%   ADAPTIVE_REQ_TEST69_NOISE_SEEDS = comma-separated noise seeds

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
write_config_json(CFG, fullfile(OUT.root_dir, 'test69_configuration.json'));

fprintf('\nTest 69: dense Test66 estimability maps using frozen Test68 risk models\n');
fprintf('Mode: %s | Test66 source: %s\n', CFG.Mode, CFG.Test66RunDir);
if strlength(CFG.RunLabel) > 0
    fprintf('Run label: %s\n', CFG.RunLabel);
end
fprintf('Frozen Test68 bundle: %s\n', CFG.Test68Bundle);
fprintf('Base Test66 q model: %s | M=%g | risk threshold %.2f\n', ...
    CFG.BaseModel, CFG.M, CFG.RiskThreshold);
fprintf('Filters: geometry=[%s] | f0=[%s] | realism=[%s] | regime=[%s]\n', ...
    join_or_all(CFG.Geometries), join_or_all(CFG.Frequencies), ...
    join_or_all(CFG.RealismLevels), join_or_all(CFG.FieldRegimes));
fprintf('Seed filters: source=[%s] | noise=[%s]\n', ...
    join_or_all(CFG.SourceSeeds), join_or_all(CFG.NoiseSeeds));
fprintf('No q, composition, correction, or risk model is retrained.\n');

assert(exist(CFG.Test66PatchCsv,'file') == 2, 'Missing Test66 patch CSV: %s', CFG.Test66PatchCsv);
assert(exist(CFG.Test66ConditionCsv,'file') == 2, 'Missing Test66 condition CSV: %s', CFG.Test66ConditionCsv);
assert(exist(CFG.Test68Bundle,'file') == 2, 'Missing Test68 bundle: %s', CFG.Test68Bundle);

COND = select_conditions(CFG);
assert(~isempty(COND), 'No Test66 conditions matched Test69 filters.');
fprintf('Selected %d dense Test66 conditions.\n', height(COND));

fprintf('Loading selected dense patch rows by chunks...\n');
T = load_selected_patch_rows(CFG, COND.condition_key);
assert(~isempty(T), 'No patch rows loaded for the selected conditions.');
fprintf('Loaded %d patch rows for Test69 overlays.\n', height(T));

T = add_operational_features(T, CFG);
assert_no_forbidden_predictors([CFG.NumericFeatureVars(:); CFG.CategoricalFeatureVars(:)]);

fprintf('Loading frozen Test68 risk ensemble...\n');
S = load(CFG.Test68Bundle, 'MODELS');
RISK_MODELS = select_test68_models(S.MODELS, CFG);
fprintf('Applying %d frozen Test68 high20 bagged-tree detectors.\n', numel(RISK_MODELS));

[risk_mean, risk_std, risk_n] = predict_test68_ensemble(T, RISK_MODELS);
T.predicted_risk = risk_mean;
T.predicted_risk_std = risk_std;
T.predicted_risk_n_models = risk_n;
T.predicted_reliability = 1 - T.predicted_risk;
T.nonestimable_mask = T.predicted_risk >= CFG.RiskThreshold;
T.estimable_mask = ~T.nonestimable_mask;
T.SWS_estimable_masked = T.SWS_pred;
T.SWS_estimable_masked(T.nonestimable_mask) = NaN;

fprintf('Writing Test69 tables...\n');
COND_SUM = summarize_group(T, ["condition_key","geometry","realism_level","field_regime","f0","M"], CFG);
ROI_SUM = summarize_group(T, ["geometry","realism_level","field_regime","f0","roi_region"], CFG);
GEOM_SUM = summarize_group(T, ["geometry"], CFG);
REAL_SUM = summarize_group(T, ["realism_level"], CFG);
FREQ_SUM = summarize_group(T, ["f0"], CFG);
REGIME_SUM = summarize_group(T, ["field_regime"], CFG);
RISK_BIN_SUM = summarize_risk_bins(T, CFG);

writetable(select_patch_columns(T), fullfile(OUT.table_dir, 'test69_dense_patch_risk_overlay.csv'));
writetable(COND_SUM, fullfile(OUT.table_dir, 'test69_dense_condition_summary.csv'));
writetable(ROI_SUM, fullfile(OUT.table_dir, 'test69_dense_roi_summary.csv'));
writetable(GEOM_SUM, fullfile(OUT.table_dir, 'test69_summary_by_geometry.csv'));
writetable(REAL_SUM, fullfile(OUT.table_dir, 'test69_summary_by_realism_level.csv'));
writetable(FREQ_SUM, fullfile(OUT.table_dir, 'test69_summary_by_frequency.csv'));
writetable(REGIME_SUM, fullfile(OUT.table_dir, 'test69_summary_by_field_regime.csv'));
writetable(RISK_BIN_SUM, fullfile(OUT.table_dir, 'test69_summary_by_risk_bin.csv'));

fprintf('Writing Test69 figures...\n');
safe_plot(@() plot_summary_figures(GEOM_SUM, REAL_SUM, FREQ_SUM, REGIME_SUM, RISK_BIN_SUM, OUT, CFG), ...
    'summary figures');
if CFG.SaveAllMaps
    safe_plot(@() plot_dense_condition_maps(T, OUT, CFG), 'dense condition maps');
end
write_readme(T, COND_SUM, ROI_SUM, OUT, CFG);
save(fullfile(OUT.data_dir, 'test69_dense_estimability_overlay_compact.mat'), ...
    'CFG','COND_SUM','ROI_SUM','GEOM_SUM','REAL_SUM','FREQ_SUM','REGIME_SUM','RISK_BIN_SUM','-v7.3');

print_console_summary(COND_SUM, ROI_SUM, OUT, CFG);
fprintf('\nTest 69 complete.\nTables: %s\nFigures: %s\nREADME: %s\n', ...
    OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST69_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST69_MODE must be quick or full.');

CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.RunLabel = env_string('ADAPTIVE_REQ_TEST69_RUN_LABEL', "");
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST69_MODEL', "q_spectrum_plus_composition");
CFG.M = env_number('ADAPTIVE_REQ_TEST69_M', 2);
CFG.CsGuess = 3.0;
CFG.RiskThreshold = env_number('ADAPTIVE_REQ_TEST69_RISK_THRESHOLD', 0.5);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST69_SAVE_ALL_MAPS', true);
CFG.MaxConditions = env_number('ADAPTIVE_REQ_TEST69_MAX_CONDITIONS', ternary(CFG.QuickMode, 10, 60));
CFG.MaxPatchRowsToSave = env_number('ADAPTIVE_REQ_TEST69_MAX_PATCH_ROWS_TO_SAVE', ternary(CFG.QuickMode, 120000, 600000));

src = strtrim(string(getenv('ADAPTIVE_REQ_TEST69_TEST66_SOURCE')));
if src == "", src = "full_a"; end
if contains(src, filesep) || startsWith(src, "/")
    CFG.Test66RunDir = char(src);
else
    CFG.Test66RunDir = fullfile(root_dir, 'outputs', 'test_66_eikonal_realistic_transfer_validation', char(src));
end
CFG.Test66PatchCsv = fullfile(CFG.Test66RunDir, 'tables', 'test66_patch_level_results.csv');
CFG.Test66ConditionCsv = fullfile(CFG.Test66RunDir, 'tables', 'test66_condition_summary.csv');

bundle = strtrim(string(getenv('ADAPTIVE_REQ_TEST69_TEST68_BUNDLE')));
if bundle == ""
    bundle = fullfile(root_dir, 'outputs', 'test_68_test66_estimability_mask', ...
        'data', 'test68_estimability_mask_compact.mat');
end
CFG.Test68Bundle = char(bundle);

if CFG.QuickMode
    CFG.Geometries = ["homogeneous_cs2","homogeneous_cs4","inclusion_2_4", ...
        "bilayer_2_3","bilayer_inclusion_2_3_4"];
    CFG.Frequencies = [300 500];
    CFG.RealismLevels = ["clean","readout_medium"];
    CFG.FieldRegimes = ["single_source_lateral","diffuse_like_8src"];
else
    CFG.Geometries = strings(0,1);
    CFG.Frequencies = [];
    CFG.RealismLevels = strings(0,1);
    CFG.FieldRegimes = strings(0,1);
end
CFG.Geometries = env_string_list('ADAPTIVE_REQ_TEST69_GEOMETRIES', CFG.Geometries);
CFG.Frequencies = env_number_list('ADAPTIVE_REQ_TEST69_FREQUENCIES', CFG.Frequencies);
CFG.RealismLevels = env_string_list('ADAPTIVE_REQ_TEST69_REALISM_LEVELS', CFG.RealismLevels);
CFG.FieldRegimes = env_string_list('ADAPTIVE_REQ_TEST69_FIELD_REGIMES', CFG.FieldRegimes);
CFG.SourceSeeds = env_number_list('ADAPTIVE_REQ_TEST69_SOURCE_SEEDS', []);
CFG.NoiseSeeds = env_number_list('ADAPTIVE_REQ_TEST69_NOISE_SEEDS', []);

CFG.NumericFeatureVars = ["q_pred","SWS_pred","k_pred","log_sws_pred","M_eff_pred", ...
    "q_theory_prior","sws_theory","abs_q_minus_theory","rel_sws_minus_theory", ...
    "predicted_patch_purity","p_mixed","p_strong_mixed","f0","M","dx","dz", ...
    "REQ_StepX","REQ_StepZ","TargetStepM","local_amplitude","tracking_snr_db", ...
    "acoustic_readout_amplitude_factor","source_to_patch_distance_m","depth_m", ...
    "grad_q_pred","grad_sws_pred","local_std_q_pred","local_std_sws_pred"];
CFG.CategoricalFeatureVars = "field_regime";
end

function OUT = make_output_dirs(root_dir, CFG)
base = fullfile(root_dir, 'outputs', 'test_69_dense_estimability_maps');
if CFG.QuickMode
    run_name = "quick";
else
    run_name = "full";
end
if strlength(CFG.RunLabel) > 0
    run_name = run_name + "_" + sanitize(CFG.RunLabel);
end
root = fullfile(base, char(run_name));
OUT = struct();
OUT.root_dir = root;
OUT.table_dir = fullfile(root, 'tables');
OUT.figure_dir = fullfile(root, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(root, 'data');
dirs = {OUT.root_dir, OUT.table_dir, OUT.figure_dir, OUT.map_dir, OUT.data_dir};
for i = 1:numel(dirs)
    if exist(dirs{i}, 'dir') ~= 7, mkdir(dirs{i}); end
end
end

%% Loading selected Test66 rows

function COND = select_conditions(CFG)
opts = detectImportOptions(CFG.Test66ConditionCsv, 'TextType','string');
COND = readtable(CFG.Test66ConditionCsv, opts);
COND = normalize_loaded_types(COND);
COND = COND(COND.model_name == CFG.BaseModel & COND.M == CFG.M, :);
if ~isempty(CFG.Geometries)
    COND = COND(ismember(COND.geometry, CFG.Geometries), :);
end
if ~isempty(CFG.Frequencies)
    COND = COND(ismember(COND.f0, CFG.Frequencies), :);
end
if ~isempty(CFG.RealismLevels)
    COND = COND(ismember(COND.realism_level, CFG.RealismLevels), :);
end
if ~isempty(CFG.FieldRegimes)
    COND = COND(ismember(COND.field_regime, CFG.FieldRegimes), :);
end
if ~isempty(CFG.SourceSeeds) && ismember("source_seed", string(COND.Properties.VariableNames))
    COND = COND(ismember(COND.source_seed, CFG.SourceSeeds), :);
end
if ~isempty(CFG.NoiseSeeds) && ismember("noise_seed", string(COND.Properties.VariableNames))
    COND = COND(ismember(COND.noise_seed, CFG.NoiseSeeds), :);
end

% Clean conditions have noise_seed 0; readout conditions may have multiple
% noise seeds. Keep the first few per physics/source combination for compact
% map output.
COND = sortrows(COND, {'geometry','realism_level','field_regime','f0','source_seed','noise_seed'});
[~, ia] = unique(COND.condition_key, 'stable');
COND = COND(ia,:);
if CFG.MaxConditions > 0 && height(COND) > CFG.MaxConditions
    COND = choose_representative_conditions(COND, CFG.MaxConditions);
end
end

function COND = choose_representative_conditions(COND, max_conditions)
% Round-robin by geometry first. This avoids alphabetic ordering selecting
% only bilayer-like cases in quick mode.
take = false(height(COND),1);
geoms = unique(COND.geometry, 'stable');
while sum(take) < max_conditions
    changed = false;
    for g = geoms'
        idx = find(COND.geometry == g & ~take, 1, 'first');
        if ~isempty(idx)
            take(idx) = true;
            changed = true;
            if sum(take) >= max_conditions
                break;
            end
        end
    end
    if ~changed
        break;
    end
end
COND = COND(take,:);
end

function T = load_selected_patch_rows(CFG, keys)
keys = string(keys(:));
vars = ["condition_key","geometry","geometry_family","case_id","case_family", ...
    "realism_level","field_regime","field_regime_ood","source_seed","noise_seed", ...
    "f0","M","REQ_M","map_iz","map_ix","x_center_m","z_center_m", ...
    "REQ_StepX","REQ_StepZ","TargetStepM","dx","dz","true_SWS","k_true", ...
    "patch_purity","patch_cs_std","patch_cs_range","patch_cs_iqr","q_oracle", ...
    "distance_to_interface_m","distance_to_interface_mm","depth_m","local_amplitude", ...
    "tracking_snr_db","snr_proxy_db","shear_attenuation_factor", ...
    "acoustic_readout_amplitude_factor","source_to_patch_distance_m", ...
    "distance_to_interface_over_window_radius","purity_bin","distance_bin", ...
    "distance_over_window_bin","depth_bin","snr_bin","amplitude_bin","roi_region", ...
    "analysis_region","q_theory_prior","sws_theory","predicted_patch_purity", ...
    "p_mixed","p_strong_mixed","model_name","q_pred","k_pred","SWS_pred","sws_pred", ...
    "q_error","k_error","sws_signed_error_pct","sws_abs_error_pct","high_error10","high_error20"];

ds = tabularTextDatastore(CFG.Test66PatchCsv, 'TextType','string');
try
    ds.SelectedVariableNames = cellstr(vars);
catch
    warning('Test69:SelectedVariableNames', 'Could not set selected variables; reading full CSV chunks.');
end
parts = {};
chunk_count = 0;
while hasdata(ds)
    C = read(ds);
    C = normalize_loaded_types(C);
    mask = C.model_name == CFG.BaseModel & C.M == CFG.M & ismember(C.condition_key, keys);
    if any(mask)
        parts{end+1,1} = C(mask,:); %#ok<AGROW>
    end
    chunk_count = chunk_count + 1;
    if mod(chunk_count, 25) == 0
        fprintf('  read %d datastore chunks, retained %d chunks.\n', chunk_count, numel(parts));
    end
end
if isempty(parts)
    T = table();
else
    T = vertcat(parts{:});
end
end

function T = normalize_loaded_types(T)
names = string(T.Properties.VariableNames);
string_vars = intersect(names, ["dataset","condition_key","geometry","geometry_family", ...
    "case_id","case_family","realism_level","field_regime","field_regime_ood", ...
    "purity_bin","distance_bin","distance_over_window_bin","depth_bin", ...
    "snr_bin","amplitude_bin","roi_region","analysis_region","model_name"], 'stable');
for v = string_vars
    T.(v) = string(T.(v));
end
numeric_vars = setdiff(names, string_vars, 'stable');
for v = numeric_vars
    if ~isnumeric(T.(v)) && ~islogical(T.(v))
        T.(v) = str2double(string(T.(v)));
    end
end
if ismember("high_error10", names), T.high_error10 = logical(T.high_error10); end
if ismember("high_error20", names), T.high_error20 = logical(T.high_error20); end
end

%% Operational Test68 feature reconstruction

function T = add_operational_features(T, CFG)
T.SWS_pred = first_existing_numeric(T, ["SWS_pred","sws_pred"]);
T.q_pred = clamp(T.q_pred, 0, 1);
T.SWS_pred = clamp(T.SWS_pred, 0.2, 20);
T.k_pred = first_existing_numeric(T, ["k_pred"]);
bad_k = ~isfinite(T.k_pred) | T.k_pred <= 0;
T.k_pred(bad_k) = 2*pi*T.f0(bad_k) ./ T.SWS_pred(bad_k);

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
for v = CFG.NumericFeatureVars(:)'
    if ~ismember(v, string(T.Properties.VariableNames))
        T.(v) = nan(height(T),1);
    end
end
end

function T = add_map_neighborhood_features(T)
T.grad_q_pred = nan(height(T),1);
T.grad_sws_pred = nan(height(T),1);
T.local_std_q_pred = nan(height(T),1);
T.local_std_sws_pred = nan(height(T),1);
conds = unique(T.condition_key, 'stable');
for ci = 1:numel(conds)
    idx = find(T.condition_key == conds(ci));
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
if isempty(RISK_MODELS)
    error('No usable Test68 high20 bagged-tree detectors were found.');
end
end

function [risk_mean, risk_std, risk_n] = predict_test68_ensemble(T, RISK_MODELS)
R = nan(height(T), numel(RISK_MODELS));
for i = 1:numel(RISK_MODELS)
    D = RISK_MODELS(i).detectors;
    X = design_matrix_for_inference(T, D.design);
    R(:,i) = predict_positive(D.high20.bagged_trees, X);
end
risk_mean = mean(R, 2, 'omitnan');
risk_std = std(R, 0, 2, 'omitnan');
risk_n = sum(isfinite(R), 2);
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

function X = onehot(x, cats)
X = zeros(numel(x), numel(cats));
for i = 1:numel(cats)
    X(:,i) = x == cats(i);
end
end

%% Summaries and plots

function P = select_patch_columns(T)
vars = ["condition_key","geometry","realism_level","field_regime","source_seed","noise_seed", ...
    "f0","M","map_iz","map_ix","x_center_m","z_center_m","true_SWS","SWS_pred", ...
    "SWS_estimable_masked","q_pred","sws_abs_error_pct","sws_signed_error_pct", ...
    "high_error10","high_error20","predicted_risk","predicted_risk_std", ...
    "predicted_reliability","estimable_mask","nonestimable_mask", ...
    "tracking_snr_db","local_amplitude","predicted_patch_purity","p_mixed", ...
    "patch_purity","distance_to_interface_over_window_radius","roi_region"];
vars = intersect(vars, string(T.Properties.VariableNames), 'stable');
P = T(:, cellstr(vars));
if height(P) > 0 && isfinite(max(T.predicted_risk_n_models))
    P.predicted_risk_n_models = T.predicted_risk_n_models;
end
end

function S = summarize_group(T, group_vars, CFG)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = {};
for gi = 1:max(G)
    X = T(G == gi,:);
    keep = X.predicted_risk < CFG.RiskThreshold;
    reject = ~keep;
    rows{end+1,1} = table(height(X), ...
        mean(X.sws_abs_error_pct,'omitnan'), median(X.sws_abs_error_pct,'omitnan'), ...
        mean(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), 100*mean(reject,'omitnan'), ...
        mean(X.sws_abs_error_pct(keep),'omitnan'), mean(X.sws_signed_error_pct(keep),'omitnan'), ...
        100*mean(X.high_error20(keep),'omitnan'), ...
        mean(X.sws_abs_error_pct(reject),'omitnan'), mean(X.sws_signed_error_pct(reject),'omitnan'), ...
        100*mean(X.high_error20(reject),'omitnan'), ...
        mean(X.predicted_risk,'omitnan'), mean(X.predicted_risk_std,'omitnan'), ...
        mean(X.tracking_snr_db,'omitnan'), mean(X.local_amplitude,'omitnan'), ...
        'VariableNames', {'N','MAPE_all_pct','median_APE_all_pct','bias_all_pct', ...
        'high10_all_pct','high20_all_pct','nonestimable_pct', ...
        'MAPE_estimable_pct','bias_estimable_pct','high20_estimable_pct', ...
        'MAPE_nonestimable_pct','bias_nonestimable_pct','high20_nonestimable_pct', ...
        'mean_risk','mean_risk_std','mean_tracking_snr_db','mean_local_amplitude'}); %#ok<AGROW>
end
S = [groups vertcat(rows{:})];
end

function S = summarize_risk_bins(T, CFG)
edges = 0:0.1:1;
labels = strings(numel(edges)-1,1);
bins = strings(height(T),1);
for i = 1:numel(edges)-1
    lo = edges(i); hi = edges(i+1);
    labels(i) = sprintf('%.1f-%.1f', lo, hi);
    idx = T.predicted_risk >= lo & T.predicted_risk < hi;
    if i == numel(edges)-1
        idx = T.predicted_risk >= lo & T.predicted_risk <= hi;
    end
    bins(idx) = labels(i);
end
T.risk_bin = categorical(bins, labels, 'Ordinal', true);
S = summarize_group(T, ["risk_bin"], CFG);
end

function plot_summary_figures(GEOM_SUM, REAL_SUM, FREQ_SUM, REGIME_SUM, RISK_BIN_SUM, OUT, CFG)
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 20]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Test69 dense estimability overlay, risk threshold %.2f', CFG.RiskThreshold), ...
    'FontWeight','bold');

ax = nexttile(tl);
bar(ax, categorical(string(GEOM_SUM.geometry)), [GEOM_SUM.MAPE_all_pct GEOM_SUM.MAPE_estimable_pct]);
xtickangle(ax,25); ylabel(ax,'MAPE (%)'); title(ax,'Geometry','FontWeight','normal');
legend(ax, {'all patches','estimable patches'}, 'Location','best'); grid(ax,'on');

ax = nexttile(tl);
bar(ax, categorical(pretty_realism(REAL_SUM.realism_level)), [REAL_SUM.MAPE_all_pct REAL_SUM.MAPE_estimable_pct]);
xtickangle(ax,20); ylabel(ax,'MAPE (%)'); title(ax,'Realism level','FontWeight','normal');
legend(ax, {'all patches','estimable patches'}, 'Location','best'); grid(ax,'on');

ax = nexttile(tl);
plot(ax, FREQ_SUM.f0, FREQ_SUM.MAPE_all_pct, '-o', FREQ_SUM.f0, FREQ_SUM.MAPE_estimable_pct, '-o');
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'Frequency','FontWeight','normal');
legend(ax, {'all patches','estimable patches'}, 'Location','best'); grid(ax,'on');

ax = nexttile(tl);
bar(ax, categorical(pretty_regime(REGIME_SUM.field_regime)), [REGIME_SUM.MAPE_all_pct REGIME_SUM.MAPE_estimable_pct]);
xtickangle(ax,20); ylabel(ax,'MAPE (%)'); title(ax,'Field regime','FontWeight','normal');
legend(ax, {'all patches','estimable patches'}, 'Location','best'); grid(ax,'on');

ax = nexttile(tl);
bar(ax, categorical(string(GEOM_SUM.geometry)), GEOM_SUM.nonestimable_pct);
xtickangle(ax,25); ylabel(ax,'Rejected pixels (%)'); title(ax,'Non-estimable fraction','FontWeight','normal'); grid(ax,'on');

ax = nexttile(tl);
plot(ax, 1:height(RISK_BIN_SUM), RISK_BIN_SUM.high20_all_pct, '-o', ...
    1:height(RISK_BIN_SUM), RISK_BIN_SUM.MAPE_all_pct, '-o');
xticks(ax, 1:height(RISK_BIN_SUM)); xticklabels(ax, string(RISK_BIN_SUM.risk_bin)); xtickangle(ax,35);
xlabel(ax,'Predicted risk bin'); ylabel(ax,'Metric'); title(ax,'Observed error versus risk','FontWeight','normal');
legend(ax, {'High-error >20% (%)','MAPE (%)'}, 'Location','best'); grid(ax,'on');

export_fig(fig, fullfile(OUT.figure_dir, 'test69_estimability_summary.png'));
end

function plot_dense_condition_maps(T, OUT, CFG)
conds = unique(T.condition_key, 'stable');
for i = 1:numel(conds)
    X = T(T.condition_key == conds(i), :);
    out_dir = fullfile(OUT.map_dir, sanitize(X.geometry(1)), sanitize(X.realism_level(1)), sanitize(X.field_regime(1)));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    plot_one_condition_map(X, fullfile(out_dir, "test69_estimability__" + sanitize(conds(i)) + ".png"), CFG);
end
end

function plot_one_condition_map(X, path, CFG)
panels = ["true_SWS","SWS_pred","SWS_estimable_masked","sws_abs_error_pct", ...
    "sws_signed_error_pct","predicted_risk","predicted_reliability","nonestimable_mask", ...
    "tracking_snr_db","local_amplitude","predicted_patch_purity","p_mixed", ...
    "patch_purity","distance_to_interface_over_window_radius","q_pred"];
titles = ["True SWS","Predicted SWS","SWS kept after mask","Absolute SWS error", ...
    "Signed SWS error","Predicted risk","Reliability","Non-estimable mask", ...
    "Tracking SNR","Local amplitude","Predicted patch purity","Predicted mixedness", ...
    "True patch purity","Distance/window radius","Predicted q"];
units = ["m/s","m/s","m/s","%","%","probability","probability","0/1", ...
    "dB","a.u.","probability","probability","fraction","ratio","REQ quantile"];

fig = figure('Color','w','Units','centimeters','Position',[1 1 40 24]);
tl = tiledlayout(fig,3,5,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s | f=%g Hz | M=%g | %s | %s | risk threshold %.2f', ...
    X.geometry(1), X.f0(1), X.M(1), pretty_realism(X.realism_level(1)), ...
    pretty_regime(X.field_regime(1)), CFG.RiskThreshold), ...
    'Interpreter','none','FontWeight','bold');
for i = 1:numel(panels)
    ax = nexttile(tl);
    A = map_from_table(X, panels(i));
    imagesc_nan(ax, A);
    axis(ax,'image'); axis(ax,'off');
    title(ax, titles(i), 'FontWeight','normal');
    cb = colorbar(ax); ylabel(cb, units(i));
end
export_fig(fig, path);
end

function imagesc_nan(ax, A)
finite = isfinite(A);
h = imagesc(ax, A);
set(h, 'AlphaData', finite);
set(ax, 'Color', [0.88 0.88 0.88]);
if any(finite(:))
    colormap(ax, parula);
end
end

function A = map_from_table(T, var)
[A, ~, ~] = vector_to_map(T.map_iz, T.map_ix, double(T.(var)));
end

%% README and console

function write_readme(T, COND_SUM, ROI_SUM, OUT, CFG)
keep = T.predicted_risk < CFG.RiskThreshold;
lines = strings(0,1);
lines(end+1) = "# Test 69 Results: Dense Test66 Estimability Maps";
lines(end+1) = "";
lines(end+1) = sprintf("- Mode: `%s`", CFG.Mode);
if strlength(CFG.RunLabel) > 0
    lines(end+1) = sprintf("- Run label: `%s`", CFG.RunLabel);
end
lines(end+1) = sprintf("- Test66 source: `%s`", CFG.Test66RunDir);
lines(end+1) = sprintf("- Test68 frozen bundle: `%s`", CFG.Test68Bundle);
lines(end+1) = sprintf("- Base q model: `%s`, M=%g", CFG.BaseModel, CFG.M);
lines(end+1) = sprintf("- Risk threshold for masking: %.2f", CFG.RiskThreshold);
lines(end+1) = sprintf("- Geometry filter: `%s`", join_or_all(CFG.Geometries));
lines(end+1) = sprintf("- Frequency filter: `%s`", join_or_all(CFG.Frequencies));
lines(end+1) = sprintf("- Realism filter: `%s`", join_or_all(CFG.RealismLevels));
lines(end+1) = sprintf("- Field-regime filter: `%s`", join_or_all(CFG.FieldRegimes));
lines(end+1) = sprintf("- Source-seed filter: `%s`", join_or_all(CFG.SourceSeeds));
lines(end+1) = sprintf("- Noise-seed filter: `%s`", join_or_all(CFG.NoiseSeeds));
lines(end+1) = "";
lines(end+1) = "No q model, composition model, correction model, or risk model was retrained. The risk map is the mean prediction of frozen Test68 high-error >20% bagged-tree detectors.";
lines(end+1) = "";
lines(end+1) = "## Overall Mask Effect";
lines(end+1) = sprintf("- Loaded dense patches: %d", height(T));
lines(end+1) = sprintf("- Rejected as non-estimable: %.1f%%", 100*mean(~keep,'omitnan'));
lines(end+1) = sprintf("- MAPE all patches: %.2f%%", mean(T.sws_abs_error_pct,'omitnan'));
lines(end+1) = sprintf("- MAPE estimable patches: %.2f%%", mean(T.sws_abs_error_pct(keep),'omitnan'));
lines(end+1) = sprintf("- MAPE non-estimable patches: %.2f%%", mean(T.sws_abs_error_pct(~keep),'omitnan'));
lines(end+1) = sprintf("- High-error >20%% all patches: %.1f%%", 100*mean(T.high_error20,'omitnan'));
lines(end+1) = sprintf("- High-error >20%% estimable patches: %.1f%%", 100*mean(T.high_error20(keep),'omitnan'));
lines(end+1) = sprintf("- High-error >20%% non-estimable patches: %.1f%%", 100*mean(T.high_error20(~keep),'omitnan'));
lines(end+1) = "";
lines(end+1) = "## Interpretation";
lines(end+1) = "- If the kept SWS map removes most high-error interface/noisy pixels while preserving core regions, Test68 is useful as a reporting mask.";
lines(end+1) = "- If core hard regions are kept but remain biased, that is a q/SWS estimator bias problem rather than a mask problem.";
lines(end+1) = "- If readout-noisy zones are rejected, the mask is responding to SNR/amplitude diagnostics as intended.";
lines(end+1) = "";
lines(end+1) = "## Main Tables";
lines(end+1) = "- `test69_dense_condition_summary.csv`";
lines(end+1) = "- `test69_dense_roi_summary.csv`";
lines(end+1) = "- `test69_dense_patch_risk_overlay.csv`";
lines(end+1) = "";
lines(end+1) = "## Map Figures";
lines(end+1) = "- Per-condition maps are under `figures/maps_by_condition/<geometry>/<realism>/<field_regime>/`.";
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w');
fprintf(fid, '%s\n', lines);
fclose(fid);

% Store compact human-readable extracts too.
writetable(COND_SUM(1:min(height(COND_SUM), 30),:), fullfile(OUT.root_dir, 'README_condition_summary_excerpt.csv'));
writetable(ROI_SUM(1:min(height(ROI_SUM), 50),:), fullfile(OUT.root_dir, 'README_roi_summary_excerpt.csv'));
end

function print_console_summary(COND_SUM, ROI_SUM, OUT, CFG)
fprintf('\n================ Test 69 summary ================\n');
fprintf('Risk threshold %.2f. Outputs: %s\n', CFG.RiskThreshold, OUT.root_dir);
disp(COND_SUM(:, {'geometry','realism_level','field_regime','f0','N', ...
    'MAPE_all_pct','MAPE_estimable_pct','MAPE_nonestimable_pct','nonestimable_pct'}));
fprintf('\nROI summary excerpt:\n');
disp(ROI_SUM(1:min(height(ROI_SUM), 12), {'geometry','roi_region','N','MAPE_all_pct', ...
    'MAPE_estimable_pct','nonestimable_pct'}));
end

%% Shared utilities

function [A, rr, cc] = vector_to_map(r, c, v)
r = double(r); c = double(c);
ur = unique(r, 'stable'); uc = unique(c, 'stable');
[~, rr] = ismember(r, ur);
[~, cc] = ismember(c, uc);
A = nan(numel(ur), numel(uc));
A(sub2ind(size(A), rr, cc)) = v;
end

function G = gradient_magnitude(A)
Af = fillmissing2(A);
[gz, gx] = gradient(Af);
G = sqrt(gx.^2 + gz.^2);
G(~isfinite(A)) = NaN;
end

function S = local_std_nan(A, win)
Af = fillmissing2(A);
k = ones(win, win);
n = conv2(double(isfinite(A)), k, 'same');
s1 = conv2(Af, k, 'same');
s2 = conv2(Af.^2, k, 'same');
mu = s1 ./ max(n, 1);
S = sqrt(max(s2 ./ max(n, 1) - mu.^2, 0));
S(n <= 1) = 0;
S(~isfinite(A)) = NaN;
end

function B = fillmissing2(A)
B = A;
if all(~isfinite(B(:)))
    B(:) = 0;
    return;
end
med = median(B(isfinite(B)), 'omitnan');
B(~isfinite(B)) = med;
end

function y = first_existing_numeric(T, names)
for n = names
    if ismember(n, string(T.Properties.VariableNames))
        y = double(T.(n));
        return;
    end
end
error('Missing all candidate variables: %s', strjoin(names, ', '));
end

function y = ensure_numeric_var(T, name, default)
if ismember(name, string(T.Properties.VariableNames))
    y = double(T.(name));
else
    y = default;
end
end

function assert_no_forbidden_predictors(features)
low = lower(string(features));
bad_exact = ["patch_purity","true_sws","k_true","q_oracle","roi_region", ...
    "analysis_region","distance_to_interface_m","distance_to_interface_mm", ...
    "distance_to_interface_over_window_radius","sws_abs_error_pct", ...
    "sws_signed_error_pct","high_error10","high_error20"];
bad_patterns = ["oracle","true","error","distance_to_interface","roi","purity_bin", ...
    "distance_bin","analysis_region"];
hit = low(ismember(low, bad_exact));
assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
for p = bad_patterns
    hit = low(contains(low, p));
    assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
end
end

function y = pretty_regime(x)
x = string(x);
y = replace(x, "_", " ");
y = replace(y, "single source lateral", "single source lateral");
y = replace(y, "diffuse like 8src", "diffuse-like 8 sources");
end

function y = pretty_realism(x)
y = replace(string(x), "_", " ");
end

function s = sanitize(x)
s = regexprep(char(string(x)), '[^A-Za-z0-9_=-]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'unnamed'; end
end

function export_fig(fig, path)
folder = fileparts(path);
if exist(folder,'dir') ~= 7, mkdir(folder); end
exportgraphics(fig, path, 'Resolution', 220);
close(fig);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Test69:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
end
end

function write_config_json(CFG, path)
try
    txt = jsonencode(CFG, 'PrettyPrint', true);
catch
    txt = jsonencode(CFG);
end
fid = fopen(path, 'w');
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function y = ternary(tf, a, b)
if tf
    y = a;
else
    y = b;
end
end

function v = env_string(name, default)
raw = string(getenv(name));
if strlength(raw) == 0
    v = string(default);
else
    v = raw;
end
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

function v = env_string_list(name, default)
raw = strtrim(string(getenv(name)));
if strlength(raw) == 0
    v = string(default);
    return;
end
parts = regexp(char(raw), '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
v = string(parts(:));
end

function v = env_number_list(name, default)
raw = strtrim(string(getenv(name)));
if strlength(raw) == 0
    v = default;
    return;
end
parts = regexp(char(raw), '[,;\s]+', 'split');
parts = parts(~cellfun(@isempty, parts));
v = str2double(string(parts(:)))';
v = v(isfinite(v));
end

function txt = join_or_all(v)
if isempty(v)
    txt = 'all';
    return;
end
if isnumeric(v)
    txt = char(strjoin(string(v(:)'), ', '));
else
    txt = char(strjoin(string(v(:)'), ', '));
end
end

function tf = env_true(name, default)
raw = lower(strtrim(string(getenv(name))));
if raw == ""
    tf = default;
else
    tf = ismember(raw, ["1","true","yes","y","on"]);
end
end
