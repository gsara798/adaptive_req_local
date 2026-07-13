%% analyze_test_17_synthetic_inclusion_kWave_like.m
% Test 17: synthetic inclusion, k-Wave-like geometry.
%
% Goal
%   Run the project-native simulator with a circular inclusion matching the
%   k-Wave geometry and compare:
%     1) deployed HybridLocalGlobal adaptive-q
%     2) discrete theory-q matched to the nominal wavefield type
%
% This test does not train models. It is designed to separate
% heterogeneity/interface error from possible k-Wave domain shift.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLegendFontSize', 11);

%% Project setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = adaptive_req.config.load_profile_config( ...
    'test_17_synthetic_inclusion_kWave_like', ...
    'RootDir', root_dir);

%% Output folders

OUT = make_output_dirs(root_dir);
checkpoint_file = fullfile(OUT.data_dir, 'level17_checkpoint.mat');

%% Load deployed model

PATHS12 = locate_latest_test12_paths(root_dir);
deployment_dir = fullfile(PATHS12.analysis_dir, ...
    'level_01_model_comparison', 'models', 'deployment');
[MODEL_DEPLOY, MODEL_INFO, model_file] = ...
    adaptive_req.analysis.load_q_model_deployment( ...
    deployment_dir, ...
    'ModelName', CFG.MODEL.ModelName, ...
    'FeatureSet', CFG.MODEL.FeatureSet, ...
    'ModelType', CFG.MODEL.ModelType);

fprintf('\nLoaded Hybrid deployment model:\n%s\n', model_file);
disp(struct2table(MODEL_INFO));

%% Main loop

T_all_pred = table();
T_region = table();
T_case_summary = table();
T_theory = table();
CASE_OUT = struct([]);
completed = strings(0, 1);

if exist(checkpoint_file, 'file') == 2
    fprintf('\nLoading Test 17 checkpoint:\n%s\n', checkpoint_file);
    S = load(checkpoint_file, 'T_all_pred', 'T_region', ...
        'T_case_summary', 'T_theory', 'CASE_OUT', 'completed');
    T_all_pred = S.T_all_pred;
    T_region = S.T_region;
    T_case_summary = S.T_case_summary;
    T_theory = S.T_theory;
    CASE_OUT = S.CASE_OUT;
    completed = S.completed;
    fprintf('Completed case/M combinations: %d\n', numel(completed));
end

for ci = 1:numel(CFG.CASES)
    case_i = CFG.CASES(ci);
    fprintf('\n=== Test 17 case %d / %d: %s ===\n', ...
        ci, numel(CFG.CASES), case_i.case_name);

    cfg_sim = build_sim_cfg(CFG, case_i, ci);
    fprintf('Running synthetic inclusion simulation: %s | Nwaves=%d | Is2D=%d\n', ...
        case_i.case_name, cfg_sim.Nwaves, cfg_sim.Is2D);
    t_sim = tic;
    sim = adaptive_req.simulate.run_single_simulation(cfg_sim);
    fprintf('Simulation completed in %.2f s. Field size: %s\n', ...
        toc(t_sim), mat2str(size(sim.Uxz)));
    print_in_plane_diagnostic(sim, case_i);

    for mi = 1:numel(CFG.REQ.M_list)
        M_i = CFG.REQ.M_list(mi);
        combo_key = sprintf('%s__M%g', case_i.wave_label, M_i);
        if any(completed == combo_key)
            fprintf('Skipping completed %s.\n', combo_key);
            continue;
        end

        fprintf('\n--- Extracting REQ/features: %s | M=%g ---\n', ...
            case_i.case_name, M_i);
        [feat_cfg, req_options] = build_req_settings(CFG, M_i);
        condition_id = (ci - 1) * numel(CFG.REQ.M_list) + mi;

        t_req = tic;
        [T_feat, req_out, global_req] = extract_feature_table( ...
            sim, cfg_sim, feat_cfg, req_options, case_i, CFG, condition_id);
        fprintf('REQ/features completed in %.2f s for %d windows.\n', ...
            toc(t_req), height(T_feat));

        q_theory = compute_theory_q(CFG, cfg_sim, feat_cfg, req_options, case_i);
        T_theory_i = make_theory_row(case_i, CFG, cfg_sim, feat_cfg, q_theory);
        T_theory = concat_tables(T_theory, T_theory_i);

        fprintf('Applying HybridLocalGlobal and theory-q maps...\n');
        T_pred_hybrid = predict_hybrid_map(MODEL_DEPLOY, MODEL_INFO, T_feat, cfg_sim);
        T_pred_theory = predict_theory_map(T_feat, cfg_sim, q_theory, case_i);
        T_pred_cond = concat_tables(T_pred_hybrid, T_pred_theory);
        T_pred_cond = add_map_diagnostics(T_pred_cond, CFG);
        T_pred_cond = add_error_metrics(T_pred_cond);

        T_region_i = summarize_region_metrics(T_pred_cond);
        T_case_i = summarize_case_metrics(T_pred_cond, sim, case_i, CFG, M_i);

        T_all_pred = concat_tables(T_all_pred, T_pred_cond);
        T_region = concat_tables(T_region, T_region_i);
        T_case_summary = concat_tables(T_case_summary, T_case_i);

        out_idx = numel(CASE_OUT) + 1;
        CASE_OUT(out_idx).case = case_i;
        CASE_OUT(out_idx).REQ_M = M_i;
        CASE_OUT(out_idx).cfg = cfg_sim;
        CASE_OUT(out_idx).sim = sim;
        CASE_OUT(out_idx).req_out = req_out;
        CASE_OUT(out_idx).global_req = global_req;
        CASE_OUT(out_idx).T_feat = T_feat;
        CASE_OUT(out_idx).T_pred = T_pred_cond;
        CASE_OUT(out_idx).q_theory = q_theory;

        save(fullfile(OUT.condition_dir, sprintf( ...
            'level17_%s_M%g.mat', sanitize_filename(case_i.wave_label), M_i)), ...
            'sim', 'cfg_sim', 'case_i', 'M_i', 'T_feat', 'T_pred_cond', ...
            'T_region_i', 'T_case_i', 'req_out', 'global_req', ...
            'q_theory', '-v7.3');

        completed(end + 1, 1) = string(combo_key); %#ok<SAGROW>
        save(checkpoint_file, 'CFG', 'MODEL_INFO', 'model_file', ...
            'T_all_pred', 'T_region', 'T_case_summary', 'T_theory', ...
            'CASE_OUT', 'completed', '-v7.3');
    end
end

assert(~isempty(T_all_pred), 'No Test 17 predictions were generated.');

%% Save outputs

writetable(remove_cell_columns(T_all_pred), fullfile(OUT.table_dir, ...
    'level17_predictions_by_pixel.csv'));
writetable(T_region, fullfile(OUT.table_dir, ...
    'level17_region_metrics.csv'));
writetable(T_case_summary, fullfile(OUT.table_dir, ...
    'level17_case_summary.csv'));
writetable(make_model_comparison(T_region), fullfile(OUT.table_dir, ...
    'level17_model_comparison.csv'));
writetable(T_theory, fullfile(OUT.table_dir, ...
    'level17_theory_q_values.csv'));

save(fullfile(OUT.data_dir, 'level17_synthetic_inclusion_kWave_like.mat'), ...
    'CFG', 'MODEL_INFO', 'model_file', 'T_all_pred', 'T_region', ...
    'T_case_summary', 'T_theory', 'CASE_OUT', '-v7.3');

%% Figures

plot_particle_velocity_maps(CASE_OUT, CFG, OUT.fig_dir);
plot_power_spectra(CASE_OUT, CFG, OUT.fig_dir);
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "HybridLocalGlobal", ...
    "cs_pred", 'level17_hybrid_sws_maps.png', 'Hybrid SWS maps', 'c_s (m/s)');
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "TheoryQDiscrete", ...
    "cs_pred", 'level17_theory_sws_maps.png', 'Theory-q SWS maps', 'c_s (m/s)');
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "HybridLocalGlobal", ...
    "cs_error_pct", 'level17_signed_error_maps_hybrid.png', ...
    'Hybrid signed SWS error maps', 'error (%)');
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "TheoryQDiscrete", ...
    "cs_error_pct", 'level17_signed_error_maps_theory.png', ...
    'Theory-q signed SWS error maps', 'error (%)');
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "HybridLocalGlobal", ...
    "q_pred", 'level17_q_maps_hybrid.png', 'Hybrid q maps', 'q');
plot_method_maps(T_all_pred, CFG, OUT.fig_dir, "TheoryQDiscrete", ...
    "q_theory_used", 'level17_q_maps_theory.png', 'Theory-q maps used', 'q');
plot_region_mape(T_region, OUT.fig_dir);
plot_error_vs_distance(T_all_pred, OUT.fig_dir);
plot_error_vs_q_gradient(T_all_pred, OUT.fig_dir);

%% Console summary

fprintf('\nTest 17 synthetic inclusion k-Wave-like analysis complete.\n');
fprintf('Analysis folder:\n%s\n', OUT.analysis_dir);
fprintf('\nRegion metrics preview:\n');
disp(sortrows(T_region(:, {'case_name','REQ_M','method_name','region_name', ...
    'MAPE_pct','bias_pct','CoV_pct','HighError_gt20_pct'}), ...
    {'case_name','REQ_M','region_name','MAPE_pct'}));

print_interpretation(T_region);

%% Local functions

function OUT = make_output_dirs(root_dir)
OUT.analysis_dir = fullfile(root_dir, 'outputs', ...
    'test_17_synthetic_inclusion_kWave_like', 'analysis');
OUT.fig_dir = fullfile(OUT.analysis_dir, 'figures');
OUT.table_dir = fullfile(OUT.analysis_dir, 'tables');
OUT.data_dir = fullfile(OUT.analysis_dir, 'data');
OUT.condition_dir = fullfile(OUT.data_dir, 'conditions');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if ~exist(dirs(i), 'dir')
        mkdir(dirs(i));
    end
end
end

function PATHS = locate_latest_test12_paths(root_dir)
out_root = fullfile(root_dir, 'outputs', 'test_12_cs_guess_window_sweep');
runs = dir(fullfile(out_root, 'test_12_cs_guess_window_sweep_*'));
runs = runs([runs.isdir]);
assert(~isempty(runs), 'No Test 12 run folders found in %s', out_root);
[~, idx] = max([runs.datenum]);
run_dir = fullfile(runs(idx).folder, runs(idx).name);
PATHS.run_dir = run_dir;
PATHS.analysis_dir = fullfile(run_dir, 'analysis');
end

function cfg = build_sim_cfg(CFG, case_i, case_idx)
cfg = struct();
cfg.Lx = CFG.SIM.Lx;
cfg.Lz = CFG.SIM.Lz;
cfg.dx = CFG.SIM.dx;
cfg.dz = CFG.SIM.dz;
cfg.f0 = CFG.SIM.f0;
cfg.cs_bg = CFG.SIM.cs_bg;
cfg.cs_inc = CFG.SIM.cs_inc;
cfg.Nwaves = case_i.Nwaves;
cfg.Is2D = logical(case_i.Is2D);
cfg.WaveModel = char(case_i.WaveModel);
cfg.AngularSamplingMethod = char(case_i.AngularSamplingMethod);
cfg.ForceInPlaneWave = logical(case_i.ForceInPlaneWave);
cfg.SNR = CFG.SIM.SNR;
cfg.AmpJitter = CFG.SIM.AmpJitter;
cfg.DecayAlpha = CFG.SIM.DecayAlpha;
cfg.Seed = CFG.SIM.Seed + 101 * case_idx;
cfg.UseParfor = CFG.SIM.UseParfor;
cfg.PhiRange = [0, 2*pi];
cfg.ThetaRange = [0, pi];
cfg.AngleRange2D = [0, 2*pi];
cfg.SourceSampling = 'ranges';
cfg.MaskConfig = struct( ...
    'cs_bg', CFG.SIM.cs_bg, ...
    'CombineMode', 'overlay', ...
    'Masks', {{struct('Type', 'circle', ...
        'cs_inc', CFG.SIM.cs_inc, ...
        'Params', struct('Center', CFG.INCLUSION.Center, ...
        'Radius', CFG.INCLUSION.Radius, ...
        'SigmaEdge', CFG.INCLUSION.SigmaEdge))}});
end

function [feat_cfg, req_options] = build_req_settings(CFG, M)
feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', M, ...
    'cs_guess', CFG.REQ.cs_guess, ...
    'gamma_win', CFG.REQ.Gamma, ...
    'pad_factor', CFG.REQ.PadFactor);
req_options = { ...
    'Nbins', CFG.REQ.Nbins, ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, ...
    'smooth_sigma', CFG.REQ.SmoothSigma};
end

function [T, OUT, global_req] = extract_feature_table( ...
    sim, cfg, feat_cfg, req_options, case_i, CFG, condition_id)

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz, cfg, req_cfg, feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features(global_curve);

global_req = struct();
global_req.q = q_global;
global_req.mapping = adaptive_req.quantile.make_req_mapping(global_curve);
global_req.features = global_features;
global_req.shape_features = global_shape;

OUT = adaptive_req.estimators.req_estimator_map( ...
    sim.Uxz, cfg, feat_cfg, ...
    'StepX', CFG.REQ.StepX, ...
    'StepZ', CFG.REQ.StepZ, ...
    'EdgeMode', CFG.REQ.EdgeMode, ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', true, ...
    'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, ...
    'UseWindowParfor', true, ...
    'StoreReqCurves', false, ...
    'Verbose', false);

T = OUT.feature_table;
T.condition_id = condition_id * ones(height(T), 1);
T.condition_label = repmat(case_i.case_name, height(T), 1);
T.case_name = repmat(case_i.case_name, height(T), 1);
T.wave_label = repmat(string(case_i.wave_label), height(T), 1);
T.SIM_WaveModel = repmat(string(case_i.case_name), height(T), 1);
T.step_idx = ones(height(T), 1);
T.realization_idx = ones(height(T), 1);
T.frequency_hz = cfg.f0 * ones(height(T), 1);
T.SIM_f0 = cfg.f0 * ones(height(T), 1);
T.SIM_cs_bg = CFG.SIM.cs_bg * ones(height(T), 1);
T.SIM_cs_inc = CFG.SIM.cs_inc * ones(height(T), 1);
T.REQ_M = feat_cfg.M * ones(height(T), 1);
cs_guess = get_cs_guess(feat_cfg);
T.REQ_cs_guess = cs_guess * ones(height(T), 1);
T.M_eff_guess = feat_cfg.M * ones(height(T), 1);
T.M_eff_true_diag = feat_cfg.M * cs_guess ./ ...
    true_cs_at_points(T.x_center_m, T.z_center_m, CFG);
T.lambda_guess = cs_guess / cfg.f0 * ones(height(T), 1);
T.cs_true = true_cs_at_points(T.x_center_m, T.z_center_m, CFG);
T.region_name = classify_region(T.x_center_m, T.z_center_m, CFG);
T.distance_to_interface_m = abs(distance_signed_to_interface( ...
    T.x_center_m, T.z_center_m, CFG));
T.global_req_mapping = repmat({global_req.mapping}, height(T), 1);
T.q_global_req = q_global * ones(height(T), 1);
T.global_REQ_Nbins_effective = global_curve.Nbins_effective * ones(height(T), 1);
T = assign_global_feature_columns(T, global_features, "global_");
T = assign_global_feature_columns(T, global_shape, "global_");
end

function q = compute_theory_q(CFG, cfg, feat_cfg, req_options, case_i)
if string(case_i.theory_field_type) == "Partial3D"
    q2 = theory_one(CFG, cfg, feat_cfg, req_options, "Diffuse2D");
    q3 = theory_one(CFG, cfg, feat_cfg, req_options, "Diffuse3D");
    q = 0.5 * (q2 + q3);
else
    q = theory_one(CFG, cfg, feat_cfg, req_options, string(case_i.theory_field_type));
end
end

function q = theory_one(~, cfg, feat_cfg, req_options, field_type)
[req_cfg, ~] = adaptive_req.config.default_req_config(cfg, feat_cfg, req_options{:});
out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    cfg.dx, cfg.dz, cfg.f0, get_cs_guess(feat_cfg), ...
    'M', feat_cfg.M, ...
    'Gamma', feat_cfg.gamma_win, ...
    'PadFactor', feat_cfg.pad_factor, ...
    'Nbins', req_cfg.Nbins, ...
    'SmoothSigma', req_cfg.smooth_sigma, ...
    'TheoryMode', 'S2D', ...
    'FieldType', field_type, ...
    'Plot', false);
q = out.q_th;
end

function cs_guess = get_cs_guess(feat_cfg)
if isfield(feat_cfg, 'cs_guess_used')
    cs_guess = feat_cfg.cs_guess_used;
elseif isfield(feat_cfg, 'cs_guess')
    cs_guess = feat_cfg.cs_guess;
else
    error('Feature config does not contain cs_guess or cs_guess_used.');
end
end

function T = make_theory_row(case_i, CFG, cfg, feat_cfg, q)
T = table();
T.case_name = string(case_i.case_name);
T.wave_label = string(case_i.wave_label);
T.theory_label = string(case_i.theory_label);
T.theory_field_type = string(case_i.theory_field_type);
T.SIM_f0 = cfg.f0;
T.SIM_dx = cfg.dx;
T.SIM_dz = cfg.dz;
T.REQ_M = feat_cfg.M;
T.REQ_cs_guess = get_cs_guess(feat_cfg);
T.q_theory_discrete = q;
if string(case_i.theory_field_type) == "Partial3D"
    T.notes = "0.5 * (diffuse 2D + projected diffuse 3D)";
else
    T.notes = "direct discrete theory quantile";
end
T.cs_bg = CFG.SIM.cs_bg;
T.cs_inc = CFG.SIM.cs_inc;
end

function T_pred = predict_hybrid_map(MODEL, MODEL_INFO, T_feat, cfg)
T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL, T_feat, ...
    'ModelType', MODEL_INFO.model_type, ...
    'ModelName', MODEL_INFO.model_name);

T_pred = base_prediction_table(T_feat);
T_pred.method_name = repmat("HybridLocalGlobal", height(T_feat), 1);
T_pred.model_name = repmat(string(MODEL_INFO.model_name), height(T_feat), 1);
T_pred.feature_set = repmat(string(MODEL_INFO.feature_set), height(T_feat), 1);
T_pred.model_type = repmat(string(MODEL_INFO.model_type), height(T_feat), 1);
T_pred.q_pred_raw = T_q.q_pred_raw;
T_pred.q_pred = min(max(T_q.q_pred, 0.001), 0.999);
T_pred.q_theory_used = nan(height(T_feat), 1);
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, cfg.f0);
end

function T_pred = predict_theory_map(T_feat, cfg, q_theory, case_i)
T_pred = base_prediction_table(T_feat);
T_pred.method_name = repmat("TheoryQDiscrete", height(T_feat), 1);
T_pred.model_name = repmat("TheoryQDiscrete", height(T_feat), 1);
T_pred.feature_set = repmat(string(case_i.theory_label), height(T_feat), 1);
T_pred.model_type = repmat("theory_no_ml", height(T_feat), 1);
T_pred.q_pred_raw = q_theory * ones(height(T_feat), 1);
T_pred.q_pred = T_pred.q_pred_raw;
T_pred.q_theory_used = T_pred.q_pred;
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, cfg.f0);
end

function T = base_prediction_table(T_feat)
vars = ["condition_id", "condition_label", "case_name", "wave_label", ...
    "step_idx", "realization_idx", "frequency_hz", "SIM_f0", ...
    "SIM_cs_bg", "SIM_cs_inc", "SIM_WaveModel", "REQ_M", ...
    "REQ_cs_guess", "M_eff_guess", "M_eff_true_diag", ...
    "patch_idx", "map_iz", "map_ix", "x_center_m", "z_center_m", ...
    "region_name", "distance_to_interface_m", "cs_true"];
vars = vars(ismember(vars, string(T_feat.Properties.VariableNames)));
T = T_feat(:, cellstr(vars));
end

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)
cs_pred = nan(numel(q_pred), 1);
for i = 1:numel(q_pred)
    if isempty(mappings{i}) || ~isfinite(q_pred(i))
        continue;
    end
    cs_pred(i) = adaptive_req.quantile.quantile_to_cs( ...
        mappings{i}, q_pred(i), f0);
end
end

function T = add_error_metrics(T)
T.cs_error = T.cs_pred - T.cs_true;
T.cs_abs_error = abs(T.cs_error);
T.cs_error_pct = 100 * T.cs_error ./ T.cs_true;
T.cs_abs_error_pct = abs(T.cs_error_pct);
end

function T = add_map_diagnostics(T, CFG)
groups = unique(T(:, {'case_name','REQ_M','method_name'}), 'rows');
T.q_gradient_mag = nan(height(T), 1);
T.cs_gradient_mag = nan(height(T), 1);
for gi = 1:height(groups)
    idx = T.case_name == groups.case_name(gi) & ...
        T.REQ_M == groups.REQ_M(gi) & ...
        T.method_name == groups.method_name(gi);
    Ti = T(idx, :);
    [Q, ~, ~] = map_from_table(Ti, 'q_pred');
    [C, ~, ~] = map_from_table(Ti, 'cs_pred');
    [dqz, dqx] = gradient(Q, CFG.SIM.dz, CFG.SIM.dx);
    [dcz, dcx] = gradient(C, CFG.SIM.dz, CFG.SIM.dx);
    Gq = sqrt(dqx.^2 + dqz.^2);
    Gc = sqrt(dcx.^2 + dcz.^2);
    T.q_gradient_mag(idx) = values_from_map(Ti, Gq);
    T.cs_gradient_mag(idx) = values_from_map(Ti, Gc);
end
end

function vals = values_from_map(T, A)
vals = nan(height(T), 1);
idx = sub2ind(size(A), T.map_iz, T.map_ix);
vals(:) = A(idx);
end

function cs = true_cs_at_points(x, z, CFG)
r = hypot(x - CFG.INCLUSION.Center(1), z - CFG.INCLUSION.Center(2));
cs = CFG.SIM.cs_bg * ones(numel(x), 1);
cs(r <= CFG.INCLUSION.Radius) = CFG.SIM.cs_inc;
end

function d = distance_signed_to_interface(x, z, CFG)
r = hypot(x - CFG.INCLUSION.Center(1), z - CFG.INCLUSION.Center(2));
d = r - CFG.INCLUSION.Radius;
end

function region = classify_region(x, z, CFG)
d = distance_signed_to_interface(x, z, CFG);
region = repmat("whole_valid_map", numel(x), 1);
core = d <= -CFG.INCLUSION.CoreMargin;
interface = abs(d) <= CFG.INCLUSION.InterfaceHalfWidth;
bg_far = d >= CFG.INCLUSION.BackgroundFarMargin;
region(core) = "inclusion_core";
region(interface) = "interface_band";
region(bg_far) = "background_far";
region(~(core | interface | bg_far)) = "transition_excluded";
end

function T_region = summarize_region_metrics(T)
T_keep = T(T.region_name ~= "transition_excluded", :);
T_whole = T;
T_whole.region_name = repmat("whole_valid_map", height(T_whole), 1);
T_sum_input = concat_tables(T_keep, T_whole);
[G, T_region] = findgroups(T_sum_input(:, {'case_name', 'wave_label', ...
    'REQ_M', 'method_name', 'model_name', 'feature_set', 'region_name'}));
T_region.N = splitapply(@numel, T_sum_input.cs_error_pct, G);
T_region.true_cs_median = splitapply(@(x) median(x, 'omitnan'), T_sum_input.cs_true, G);
T_region.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), T_sum_input.cs_error_pct, G);
T_region.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T_sum_input.cs_error_pct, G);
T_region.CoV_pct = splitapply(@cov_pct, T_sum_input.cs_pred, G);
T_region.median_error_pct = splitapply(@(x) median(x, 'omitnan'), T_sum_input.cs_error_pct, G);
T_region.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T_sum_input.cs_error_pct, G);
T_region.P95_abs_error_pct = splitapply(@(x) prctile(abs(x), 95), T_sum_input.cs_error_pct, G);
T_region.HighError_gt10_pct = splitapply(@(x) 100 * mean(abs(x) > 10, 'omitnan'), T_sum_input.cs_error_pct, G);
T_region.HighError_gt20_pct = splitapply(@(x) 100 * mean(abs(x) > 20, 'omitnan'), T_sum_input.cs_error_pct, G);
T_region.cs_pred_mean = splitapply(@(x) mean(x, 'omitnan'), T_sum_input.cs_pred, G);
T_region.cs_pred_std = splitapply(@(x) std(x, 'omitnan'), T_sum_input.cs_pred, G);
T_region = sortrows(T_region, {'case_name','REQ_M','region_name','MAPE_pct'});
end

function T = summarize_case_metrics(T_pred_cond, sim, case_i, CFG, M)
T = table();
T.case_name = string(case_i.case_name);
T.wave_label = string(case_i.wave_label);
T.REQ_M = M;
T.in_plane_fraction = in_plane_fraction(sim);
T.min_abs_uy = min(abs(sim.diag.waveDirs.uy));
T.Nwaves = case_i.Nwaves;
T.Is2D = logical(case_i.Is2D);
T.field_abs_mean = mean(abs(sim.Uxz(:)), 'omitnan');
T.field_abs_std = std(abs(sim.Uxz(:)), 0, 'omitnan');
Th = T_pred_cond(T_pred_cond.method_name == "HybridLocalGlobal" & ...
    T_pred_cond.region_name == "interface_band", :);
T.hybrid_interface_MAPE_pct = mean(Th.cs_abs_error_pct, 'omitnan');
T.hybrid_interface_HighError_gt20_pct = 100 * mean(Th.cs_abs_error_pct > 20, 'omitnan');
end

function y = cov_pct(x)
x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end
end

function plot_particle_velocity_maps(CASE_OUT, CFG, fig_dir)
cases = unique(string(arrayfun(@(s) string(s.case.case_name), CASE_OUT)), 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 7.5*numel(cases) 7]);
tl = tiledlayout(1, numel(cases), 'TileSpacing', 'compact', 'Padding', 'compact');
for ci = 1:numel(cases)
    idx = find(arrayfun(@(s) string(s.case.case_name) == cases(ci), CASE_OUT), 1);
    sim = CASE_OUT(idx).sim;
    ax = nexttile(tl);
    imagesc(ax, sim.x * 100, sim.z * 100, real(sim.Uxz));
    axis(ax, 'image'); set(ax, 'YDir', 'normal', 'FontSize', 9);
    colormap(ax, parula); colorbar(ax);
    draw_inclusion(ax, CFG);
    title(ax, cases(ci), 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 10);
    xlabel(ax, 'x (cm)'); ylabel(ax, 'z (cm)');
end
title(tl, 'Particle velocity maps, real component used by REQ', ...
    'FontWeight', 'normal');
export_clean(fig, fullfile(fig_dir, 'level17_particle_velocity_maps.png'), CFG);
end

function plot_power_spectra(CASE_OUT, CFG, fig_dir)
cases = unique(string(arrayfun(@(s) string(s.case.case_name), CASE_OUT)), 'stable');
M_values = CFG.REQ.M_list(:)';
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 6.7*numel(cases) 5.8*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        idx = find(arrayfun(@(s) string(s.case.case_name) == cases(ci) && ...
            s.REQ_M == M_values(mi), CASE_OUT), 1);
        sim = CASE_OUT(idx).sim;
        win = CASE_OUT(idx).req_out.win_size;
        patch = center_patch(sim.Uxz, win);
        S = log1p(abs(fftshift(fft2(patch))).^2);
        ax = nexttile(tl);
        imagesc(ax, S);
        axis(ax, 'image'); set(ax, 'YDir', 'normal', 'FontSize', 8);
        colormap(ax, turbo); colorbar(ax);
        title(ax, sprintf('%s | M=%g', cases(ci), M_values(mi)), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
        xlabel(ax, 'k_x bin'); ylabel(ax, 'k_z bin');
    end
end
title(tl, 'Representative central-window power spectra log(1+|FFT|^2)', ...
    'FontWeight', 'normal');
export_clean(fig, fullfile(fig_dir, 'level17_power_spectra.png'), CFG);
end

function patch = center_patch(U, win)
[nz, nx] = size(U);
cz = round((nz + 1) / 2);
cx = round((nx + 1) / 2);
half = floor(win / 2);
patch = U((cz-half):(cz+half), (cx-half):(cx+half));
end

function plot_method_maps(T, CFG, fig_dir, method_name, value_var, file_name, title_text, cbar_label)
T = T(T.method_name == method_name, :);
cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 7.2*numel(cases) 5.8*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Tc = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile(tl);
        [A, x_cm, z_cm] = map_from_table(Tc, value_var);
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image'); set(ax, 'YDir', 'normal', 'FontSize', 8);
        if contains(value_var, "error")
            colormap(ax, parula);
        else
            colormap(ax, turbo);
        end
        cb = colorbar(ax);
        ylabel(cb, cbar_label);
        draw_inclusion(ax, CFG);
        mape = mean(Tc.cs_abs_error_pct, 'omitnan');
        title(ax, sprintf('%s | M=%g | MAPE %.2f%%', ...
            cases(ci), M_values(mi), mape), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
        xlabel(ax, 'x (cm)'); ylabel(ax, 'z (cm)');
    end
end
title(tl, title_text, 'FontWeight', 'normal');
export_clean(fig, fullfile(fig_dir, file_name), CFG);
end

function plot_region_mape(T_region, fig_dir)
T = T_region(T_region.region_name ~= "whole_valid_map", :);
regions = ["background_far", "inclusion_core", "interface_band"];
cases = unique(T.case_name, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 6.3*numel(cases) 12.5]);
tl = tiledlayout(numel(regions), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for ri = 1:numel(regions)
    for ci = 1:numel(cases)
        ax = nexttile(tl);
        Tc = T(T.case_name == cases(ci) & T.region_name == regions(ri), :);
        plot_grouped_lines(ax, Tc, "REQ_M", "method_name", "MAPE_pct");
        title(ax, sprintf('%s | %s', cases(ci), regions(ri)), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 7.5);
        ylabel(ax, 'MAPE (%)', 'FontSize', 8);
        xlabel(ax, 'REQ M', 'FontSize', 8);
        set(ax, 'FontSize', 8);
        grid(ax, 'on');
    end
end
lg = legend(tl.Children(end), {'HybridLocalGlobal', 'TheoryQDiscrete'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'Interpreter', 'none', 'FontSize', 8);
lg.Layout.Tile = 'north';
title(tl, 'Regional MAPE: Hybrid vs Theory-q', ...
    'FontWeight', 'normal', 'FontSize', 11);
export_clean(fig, fullfile(fig_dir, 'level17_region_mape_hybrid_vs_theory.png'), ...
    struct('FIG', struct('Resolution', 260)));
end

function plot_error_vs_distance(T, fig_dir)
T = T(T.region_name ~= "transition_excluded", :);
methods = unique(T.method_name, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8*numel(methods) 10]);
tl = tiledlayout(1, numel(methods), 'TileSpacing', 'compact', 'Padding', 'compact');
rng(1701);
for mi = 1:numel(methods)
    ax = nexttile(tl);
    idx = find(T.method_name == methods(mi));
    idx = idx(isfinite(T.distance_to_interface_m(idx)) & isfinite(T.cs_abs_error_pct(idx)));
    if numel(idx) > 40000
        idx = idx(randperm(numel(idx), 40000));
    end
    scatter(ax, T.distance_to_interface_m(idx) * 1000, ...
        T.cs_abs_error_pct(idx), 5, 'filled', 'MarkerFaceAlpha', 0.12);
    yline(ax, 20, 'k--');
    grid(ax, 'on');
    title(ax, methods(mi), 'Interpreter', 'none', 'FontWeight', 'normal');
    xlabel(ax, 'distance to interface (mm)');
    ylabel(ax, '|SWS error| (%)');
end
title(tl, 'Error vs distance to inclusion interface', 'FontWeight', 'normal');
export_clean(fig, fullfile(fig_dir, 'level17_error_vs_distance_to_interface.png'), ...
    struct('FIG', struct('Resolution', 260)));
end

function plot_error_vs_q_gradient(T, fig_dir)
T = T(T.method_name == "HybridLocalGlobal", :);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 11]);
idx = find(isfinite(T.q_gradient_mag) & isfinite(T.cs_abs_error_pct));
rng(1702);
if numel(idx) > 50000
    idx = idx(randperm(numel(idx), 50000));
end
scatter(T.q_gradient_mag(idx), T.cs_abs_error_pct(idx), 5, ...
    'filled', 'MarkerFaceAlpha', 0.12);
yline(20, 'k--');
grid on;
xlabel('|grad q|');
ylabel('|SWS error| (%)');
title('Hybrid error vs q-map gradient', 'FontWeight', 'normal');
export_clean(fig, fullfile(fig_dir, 'level17_error_vs_q_gradient.png'), ...
    struct('FIG', struct('Resolution', 260)));
end

function plot_grouped_lines(ax, T, x_var, series_var, metric_var)
x_values = unique(T.(char(x_var)), 'stable');
series_values = ["HybridLocalGlobal", "TheoryQDiscrete"];
Y = nan(numel(x_values), numel(series_values));
for i = 1:numel(x_values)
    for j = 1:numel(series_values)
        idx = T.(char(x_var)) == x_values(i) & ...
            T.(char(series_var)) == series_values(j);
        if any(idx)
            Y(i, j) = mean(T.(char(metric_var))(idx), 'omitnan');
        end
    end
end
hold(ax, 'on');
colors = lines(numel(series_values));
for j = 1:numel(series_values)
    plot(ax, x_values, Y(:, j), '-o', 'LineWidth', 1.2, ...
        'MarkerSize', 4, 'Color', colors(j, :));
end
hold(ax, 'off');
xticks(ax, x_values);
xlim(ax, [min(x_values)-0.15, max(x_values)+0.15]);
end

function [A, x_cm, z_cm] = map_from_table(T, value_var)
nz = max(T.map_iz);
nx = max(T.map_ix);
A = nan(nz, nx);
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
A(idx) = T.(char(value_var));
x_cm = unique(T.x_center_m, 'stable') * 100;
z_cm = unique(T.z_center_m, 'stable') * 100;
end

function draw_inclusion(ax, CFG)
hold(ax, 'on');
th = linspace(0, 2*pi, 300);
xc = CFG.INCLUSION.Center(1) * 100;
zc = CFG.INCLUSION.Center(2) * 100;
r = CFG.INCLUSION.Radius * 100;
plot(ax, xc + r*cos(th), zc + r*sin(th), 'w--', 'LineWidth', 1.1);
hold(ax, 'off');
end

function T = make_model_comparison(T_region)
T = T_region(T_region.region_name == "whole_valid_map", :);
T = sortrows(T, {'case_name','REQ_M','MAPE_pct'});
end

function T = assign_global_feature_columns(T, values, prefix)
if ~isstruct(values)
    return;
end
names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if isnumeric(value_i) && isscalar(value_i)
        T.(char(string(prefix) + string(names{i}))) = ...
            double(value_i) * ones(height(T), 1);
    elseif islogical(value_i) && isscalar(value_i)
        T.(char(string(prefix) + string(names{i}))) = ...
            double(value_i) * ones(height(T), 1);
    end
end
end

function T = concat_tables(A, B)
if isempty(A)
    T = B;
    return;
end
if isempty(B)
    T = A;
    return;
end
vars_all = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');
A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);
T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];
end

function T = add_missing_columns(T, vars_all)
for i = 1:numel(vars_all)
    v = char(vars_all(i));
    if ismember(v, T.Properties.VariableNames)
        continue;
    end
    T.(v) = missing_column(height(T));
end
end

function x = missing_column(n)
x = nan(n, 1);
end

function T = remove_cell_columns(T)
vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];
end

function name = sanitize_filename(x)
name = regexprep(char(string(x)), '[^A-Za-z0-9_=-]+', '_');
end

function export_clean(fig, file_path, CFG)
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end
drawnow;
res = 260;
if isfield(CFG, 'FIG') && isfield(CFG.FIG, 'Resolution')
    res = CFG.FIG.Resolution;
end
exportgraphics(fig, file_path, 'Resolution', res, 'BackgroundColor', 'white');
close(fig);
end

function print_in_plane_diagnostic(sim, case_i)
covg = sim.diag.inPlaneCoverage;
fprintf('In-plane source fraction: %.3f | min |uy| = %.3g\n', ...
    in_plane_fraction(sim), min(abs(sim.diag.waveDirs.uy)));
if string(case_i.wave_label) == "partial_diffuse_3d"
    assert(min(abs(sim.diag.waveDirs.uy)) < 1e-7, ...
        'Partially diffuse 3D case does not contain an in-plane source.');
end
end

function frac = in_plane_fraction(sim)
covg = sim.diag.inPlaneCoverage;
if isfield(covg, 'fraction_in_plane')
    frac = covg.fraction_in_plane;
elseif isfield(covg, 'count_in_plane') && isfield(covg, 'n_waves') && covg.n_waves > 0
    frac = covg.count_in_plane / covg.n_waves;
else
    frac = NaN;
end
end

function print_interpretation(T_region)
Hybrid = T_region(T_region.method_name == "HybridLocalGlobal", :);
Hwhole = Hybrid(Hybrid.region_name == "whole_valid_map", :);
Hint = Hybrid(Hybrid.region_name == "interface_band", :);
Hcore = Hybrid(Hybrid.region_name == "inclusion_core", :);
Hbg = Hybrid(Hybrid.region_name == "background_far", :);

fprintf('\nAutomatic interpretation guide:\n');
if isempty(Hint) || isempty(Hwhole)
    fprintf('Insufficient regional data for interpretation.\n');
    return;
end

interface_mape = mean(Hint.MAPE_pct, 'omitnan');
core_mape = mean(Hcore.MAPE_pct, 'omitnan');
bg_mape = mean(Hbg.MAPE_pct, 'omitnan');
whole_mape = mean(Hwhole.MAPE_pct, 'omitnan');

fprintf('Mean Hybrid MAPE, whole valid map: %.2f%%\n', whole_mape);
fprintf('Mean Hybrid MAPE, background far: %.2f%%\n', bg_mape);
fprintf('Mean Hybrid MAPE, inclusion core: %.2f%%\n', core_mape);
fprintf('Mean Hybrid MAPE, interface band: %.2f%%\n', interface_mape);

if interface_mape > 1.5 * max([bg_mape, core_mape], [], 'omitnan')
    fprintf(['Interpretation: error is concentrated at the synthetic interface. ', ...
        'This supports heterogeneity/material mixing inside the REQ window as ', ...
        'a major failure mode.\n']);
elseif whole_mape < 5
    fprintf(['Interpretation: synthetic heterogeneous maps remain relatively accurate. ', ...
        'If k-Wave remains much worse, k-Wave domain shift is likely important.\n']);
else
    fprintf(['Interpretation: errors are not isolated to the interface. ', ...
        'Inspect wavefield-specific maps and spectra before assigning the cause.\n']);
end
end
