%% analyze_test_17_model_comparison_heterogeneous_cases.m
% Test 17 extension: compare existing q models across homogeneous and
% heterogeneous synthetic maps.
%
% This script does not train any model. LocalOnly, GlobalOnly, and
% HybridLocalGlobal are loaded from the model registry. TheoryQDiscrete is
% a diagnostic no-ML baseline. Each simulation and each local REQ spectrum
% extraction is shared by all four methods.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultTextFontSize', 11);
set(groot, 'defaultLegendFontSize', 10);

%% Project setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = adaptive_req.config.load_profile_config( ...
    'test_17_synthetic_inclusion_kWave_like', ...
    'RootDir', root_dir);

% The comparison uses the exact Test 17 numerical settings.
CFG.REQ.StepX = 1;
CFG.REQ.StepZ = 1;
CFG.REQ.EdgeMode = 'valid';
CFG.FIG.Resolution = 220;

GEOMETRIES = build_geometry_specs(CFG);
OUT = make_output_dirs(root_dir);
checkpoint_file = fullfile(OUT.data_dir, 'level17_model_comparison_checkpoint.mat');

%% Load existing operational models

MODEL_SPECS = load_registry_models(root_dir);
validate_operational_model_predictors(MODEL_SPECS);

fprintf('\nLoaded operational models from the registry:\n');
disp(struct2table(rmfield(MODEL_SPECS, 'model')));
fprintf(['TheoryQDiscrete is a diagnostic no-ML baseline. It is not an ', ...
    'operational learned model.\n']);

%% Resume state

T_pixel = table();
T_region = table();
T_case = table();
T_theory = table();
CASE_OUT = struct([]);
completed = strings(0, 1);

if exist(checkpoint_file, 'file') == 2
    fprintf('\nLoading checkpoint:\n%s\n', checkpoint_file);
    S = load(checkpoint_file, 'T_pixel', 'T_region', 'T_case', ...
        'T_theory', 'CASE_OUT', 'completed');
    T_pixel = S.T_pixel;
    T_region = S.T_region;
    T_case = S.T_case;
    T_theory = S.T_theory;
    CASE_OUT = S.CASE_OUT;
    completed = S.completed;
    fprintf('Completed geometry/wavefield/M combinations: %d / %d\n', ...
        numel(completed), numel(GEOMETRIES) * numel(CFG.CASES) * numel(CFG.REQ.M_list));
end

%% Main experiment

condition_id = 0;
if ~isempty(T_pixel) && ismember('condition_id', T_pixel.Properties.VariableNames)
    condition_id = max(T_pixel.condition_id, [], 'omitnan');
end
for gi = 1:numel(GEOMETRIES)
    geometry = GEOMETRIES(gi);

    for wi = 1:numel(CFG.CASES)
        wave_case = CFG.CASES(wi);
        sim_key = sprintf('%s__%s', geometry.geometry_id, wave_case.wave_label);

        pending = false;
        for mi = 1:numel(CFG.REQ.M_list)
            key_i = sprintf('%s__M%g', sim_key, CFG.REQ.M_list(mi));
            pending = pending || ~any(completed == string(key_i));
        end
        if ~pending
            fprintf('\nSkipping completed simulation group %s.\n', sim_key);
            continue;
        end

        cfg_sim = build_sim_cfg(CFG, geometry, wave_case, gi, wi);
        fprintf('\n============================================================\n');
        fprintf('Geometry: %s | Wavefield: %s\n', ...
            geometry.geometry_name, wave_case.case_name);

        prior_file = prior_inclusion_condition_file( ...
            root_dir, geometry, wave_case, CFG.REQ.M_list(1));
        if strlength(prior_file) > 0
            S_prior = load(prior_file, 'sim');
            sim = S_prior.sim;
            fprintf('Reusing Test 17 inclusion simulation from:\n%s\n', prior_file);
        else
            t_sim = tic;
            sim = adaptive_req.simulate.run_single_simulation(cfg_sim);
            fprintf('Simulation completed in %.2f s | field size %s\n', ...
                toc(t_sim), mat2str(size(sim.Uxz)));
        end
        print_in_plane_diagnostic(sim, wave_case);

        for mi = 1:numel(CFG.REQ.M_list)
            M_i = CFG.REQ.M_list(mi);
            combo_key = sprintf('%s__M%g', sim_key, M_i);
            if any(completed == string(combo_key))
                fprintf('Skipping completed %s.\n', combo_key);
                continue;
            end

            condition_id = condition_id + 1;
            fprintf('\n--- M=%g: shared REQ/features extraction ---\n', M_i);
            [feat_cfg, req_options] = build_req_settings(CFG, M_i);

            prior_file = prior_inclusion_condition_file( ...
                root_dir, geometry, wave_case, M_i);
            if strlength(prior_file) > 0
                S_prior = load(prior_file, 'T_feat', 'req_out', 'global_req');
                T_feat = attach_feature_metadata(S_prior.T_feat, ...
                    cfg_sim, feat_cfg, geometry, wave_case, condition_id);
                req_out = S_prior.req_out;
                global_req = S_prior.global_req;
                fprintf('Reused Test 17 inclusion REQ/features for %d windows.\n', ...
                    height(T_feat));
            else
                t_req = tic;
                [T_feat, req_out, global_req] = extract_feature_table( ...
                    sim, cfg_sim, feat_cfg, req_options, geometry, ...
                    wave_case, CFG, condition_id);
                fprintf('REQ/features completed in %.2f s for %d windows.\n', ...
                    toc(t_req), height(T_feat));
            end

            q_theory = compute_theory_q(cfg_sim, feat_cfg, req_options, wave_case);
            T_theory = concat_tables(T_theory, make_theory_row( ...
                geometry, wave_case, cfg_sim, feat_cfg, q_theory));

            fprintf('Applying LocalOnly, GlobalOnly, HybridLocalGlobal, and theory q...\n');
            T_pred = table();
            for model_idx = 1:numel(MODEL_SPECS)
                T_pred = concat_tables(T_pred, predict_operational_model( ...
                    MODEL_SPECS(model_idx), T_feat, cfg_sim));
            end
            T_pred = concat_tables(T_pred, predict_theory_model( ...
                T_feat, cfg_sim, q_theory, wave_case));
            T_pred = add_error_metrics(T_pred);
            T_pred = add_map_gradients(T_pred, CFG);

            T_region_i = summarize_region_metrics(T_pred);
            T_case_i = summarize_case_metrics(T_pred, sim, geometry, wave_case, M_i);

            T_pixel = concat_tables(T_pixel, T_pred);
            T_region = concat_tables(T_region, T_region_i);
            T_case = concat_tables(T_case, T_case_i);

            out_idx = numel(CASE_OUT) + 1;
            CASE_OUT(out_idx).geometry = geometry;
            CASE_OUT(out_idx).wave_case = wave_case;
            CASE_OUT(out_idx).REQ_M = M_i;
            CASE_OUT(out_idx).cfg = cfg_sim;
            CASE_OUT(out_idx).sim = sim;
            CASE_OUT(out_idx).req_out = req_out;
            CASE_OUT(out_idx).global_req = global_req;
            CASE_OUT(out_idx).T_pred = T_pred;
            CASE_OUT(out_idx).q_theory = q_theory;

            condition_file = fullfile(OUT.condition_dir, sprintf( ...
                'level17_compare_%s_%s_M%g.mat', geometry.geometry_id, ...
                sanitize_filename(wave_case.wave_label), M_i));
            save(condition_file, 'geometry', 'wave_case', 'M_i', ...
                'cfg_sim', 'sim', 'req_out', 'global_req', 'T_feat', ...
                'T_pred', 'T_region_i', 'T_case_i', 'q_theory', '-v7.3');

            completed(end + 1, 1) = string(combo_key); %#ok<SAGROW>
            save(checkpoint_file, 'CFG', 'MODEL_SPECS', 'T_pixel', ...
                'T_region', 'T_case', 'T_theory', 'CASE_OUT', ...
                'completed', '-v7.3');
        end
    end
end

assert(~isempty(T_pixel), 'No model-comparison predictions are available.');

%% Tables and MAT output

T_model_comparison = make_model_comparison(T_region);
T_distance = summarize_error_vs_distance(T_pixel);

writetable(remove_cell_columns(T_pixel), fullfile(OUT.table_dir, ...
    'level17_model_comparison_predictions_by_pixel.csv'));
writetable(T_region, fullfile(OUT.table_dir, ...
    'level17_model_comparison_region_metrics.csv'));
writetable(T_model_comparison, fullfile(OUT.table_dir, ...
    'level17_model_comparison_summary.csv'));
writetable(T_case, fullfile(OUT.table_dir, ...
    'level17_model_comparison_case_summary.csv'));
writetable(T_theory, fullfile(OUT.table_dir, ...
    'level17_model_comparison_theory_q_values.csv'));
writetable(T_distance, fullfile(OUT.table_dir, ...
    'level17_model_comparison_error_vs_interface_distance.csv'));

save(fullfile(OUT.data_dir, 'level17_model_comparison_heterogeneous_cases.mat'), ...
    'CFG', 'GEOMETRIES', 'MODEL_SPECS', 'T_pixel', 'T_region', ...
    'T_model_comparison', 'T_case', 'T_theory', 'T_distance', ...
    'CASE_OUT', '-v7.3');

%% Figures

plot_particle_velocity_overview(CASE_OUT, GEOMETRIES, CFG, OUT.fig_dir);
plot_power_spectrum_overview(CASE_OUT, GEOMETRIES, CFG, OUT.fig_dir);
plot_all_comparison_maps(T_pixel, GEOMETRIES, CFG, OUT.fig_dir);
plot_metric_bars(T_region, "MAPE_pct", 'MAPE (%)', ...
    'level17_compare_mape_by_model_case_wavefield_M', CFG, OUT.fig_dir);
plot_metric_bars(T_region, "HighError_gt20_pct", 'high-error >20% (%)', ...
    'level17_compare_high_error_by_model_case_wavefield_M', CFG, OUT.fig_dir);
plot_error_vs_distance(T_distance, GEOMETRIES, CFG, OUT.fig_dir);
plot_error_vs_q_gradient(T_pixel, CFG, OUT.fig_dir);

%% Console conclusions

fprintf('\nTest 17 model comparison complete.\n');
fprintf('Analysis folder:\n%s\n', OUT.analysis_dir);
print_summary(T_region);

%% Local functions

function OUT = make_output_dirs(root_dir)
OUT.analysis_dir = fullfile(root_dir, 'outputs', ...
    'test_17_model_comparison_heterogeneous_cases', 'analysis');
OUT.fig_dir = fullfile(OUT.analysis_dir, 'figures');
OUT.table_dir = fullfile(OUT.analysis_dir, 'tables');
OUT.data_dir = fullfile(OUT.analysis_dir, 'data');
OUT.condition_dir = fullfile(OUT.data_dir, 'conditions');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if exist(dirs(i), 'dir') ~= 7
        mkdir(dirs(i));
    end
end
end

function specs = build_geometry_specs(CFG)
specs = [
    struct('geometry_id', "homogeneous_cs2", ...
        'geometry_name', "Homogeneous c_s=2 m/s", ...
        'geometry_type', "homogeneous", ...
        'cs_low', 2.0, 'cs_high', 2.0, ...
        'center', CFG.INCLUSION.Center, 'radius', CFG.INCLUSION.Radius, ...
        'interface_half_width', CFG.INCLUSION.InterfaceHalfWidth, ...
        'far_margin', CFG.INCLUSION.BackgroundFarMargin)
    struct('geometry_id', "homogeneous_cs3", ...
        'geometry_name', "Homogeneous c_s=3 m/s", ...
        'geometry_type', "homogeneous", ...
        'cs_low', 3.0, 'cs_high', 3.0, ...
        'center', CFG.INCLUSION.Center, 'radius', CFG.INCLUSION.Radius, ...
        'interface_half_width', CFG.INCLUSION.InterfaceHalfWidth, ...
        'far_margin', CFG.INCLUSION.BackgroundFarMargin)
    struct('geometry_id', "bilayer_2_3", ...
        'geometry_name', "Bilayer 2/3 m/s", ...
        'geometry_type', "bilayer", ...
        'cs_low', 2.0, 'cs_high', 3.0, ...
        'center', [0.025 0.025], 'radius', NaN, ...
        'interface_half_width', CFG.INCLUSION.InterfaceHalfWidth, ...
        'far_margin', CFG.INCLUSION.BackgroundFarMargin)
    struct('geometry_id', "inclusion_2_3", ...
        'geometry_name', "Circular inclusion 2/3 m/s", ...
        'geometry_type', "inclusion", ...
        'cs_low', 2.0, 'cs_high', 3.0, ...
        'center', CFG.INCLUSION.Center, 'radius', CFG.INCLUSION.Radius, ...
        'interface_half_width', CFG.INCLUSION.InterfaceHalfWidth, ...
        'far_margin', CFG.INCLUSION.BackgroundFarMargin)
    ];
end

function file_path = prior_inclusion_condition_file(root_dir, geometry, wave_case, M)
file_path = "";
if string(geometry.geometry_type) ~= "inclusion"
    return;
end
candidate = fullfile(root_dir, 'outputs', ...
    'test_17_synthetic_inclusion_kWave_like', 'analysis', 'data', ...
    'conditions', sprintf('level17_%s_M%g.mat', ...
    sanitize_filename(wave_case.wave_label), M));
if exist(candidate, 'file') == 2
    file_path = string(candidate);
end
end

function MODEL_SPECS = load_registry_models(root_dir)
model_dir = fullfile(root_dir, 'outputs', 'model_registry', ...
    'test12_hybrid_baseline');
requests = [
    struct('name', "LocalOnly", 'feature_set', "NoCsGuess")
    struct('name', "GlobalOnly", 'feature_set', "NoCsGuess")
    struct('name', "HybridLocalGlobal", 'feature_set', "WithCsGuess")
    ];

MODEL_SPECS = struct([]);
for i = 1:numel(requests)
    [model_i, info_i, file_i] = adaptive_req.analysis.load_q_model_deployment( ...
        model_dir, ...
        'ModelName', requests(i).name, ...
        'FeatureSet', requests(i).feature_set, ...
        'ModelType', 'bagged_trees');
    MODEL_SPECS(i).model_name = string(info_i.model_name);
    MODEL_SPECS(i).feature_set = string(info_i.feature_set);
    MODEL_SPECS(i).model_type = string(info_i.model_type);
    MODEL_SPECS(i).model_role = "operational";
    MODEL_SPECS(i).model_file = string(file_i);
    MODEL_SPECS(i).model = model_i;
end
end

function validate_operational_model_predictors(MODEL_SPECS)
forbidden = lower([
    "q_theory"; "q_true"; "q_global_theory"; "q_local_minus_global"
    "abs_q_local_minus_global"; "M_eff_true_diag"; "cs_true"; "cs_pred"
    "sws_error"; "abs_sws_error"; "sws_error_pct"; "abs_sws_error_pct"
    "residual"; "abs_error"; "aperture_weight"; "solid_angle_weight"
    "true_aperture_weight"]);
for i = 1:numel(MODEL_SPECS)
    predictors = lower(string({MODEL_SPECS(i).model.encoder.entries.name}));
    leaked = intersect(predictors, forbidden);
    assert(isempty(leaked), ...
        'Operational model %s contains forbidden predictors: %s', ...
        MODEL_SPECS(i).model_name, strjoin(leaked, ', '));
end
end

function cfg = build_sim_cfg(CFG, geometry, wave_case, geometry_idx, wave_idx)
cfg = struct();
cfg.Lx = CFG.SIM.Lx;
cfg.Lz = CFG.SIM.Lz;
cfg.dx = CFG.SIM.dx;
cfg.dz = CFG.SIM.dz;
cfg.f0 = CFG.SIM.f0;
cfg.cs_bg = geometry.cs_low;
cfg.cs_inc = geometry.cs_high;
cfg.Nwaves = wave_case.Nwaves;
cfg.Is2D = logical(wave_case.Is2D);
cfg.WaveModel = char(wave_case.WaveModel);
cfg.AngularSamplingMethod = char(wave_case.AngularSamplingMethod);
cfg.ForceInPlaneWave = logical(wave_case.ForceInPlaneWave);
cfg.SNR = CFG.SIM.SNR;
cfg.AmpJitter = CFG.SIM.AmpJitter;
cfg.DecayAlpha = CFG.SIM.DecayAlpha;
cfg.Seed = CFG.SIM.Seed + 1000 * geometry_idx + 101 * wave_idx;
cfg.UseParfor = CFG.SIM.UseParfor;
cfg.PhiRange = [0, 2*pi];
cfg.ThetaRange = [0, pi];
cfg.AngleRange2D = [0, 2*pi];
cfg.SourceSampling = 'ranges';

switch geometry.geometry_type
    case "homogeneous"
        masks = {};
    case "inclusion"
        masks = {struct('Type', 'circle', 'cs_inc', geometry.cs_high, ...
            'Params', struct('Center', geometry.center, ...
            'Radius', geometry.radius, 'SigmaEdge', 1e-6))};
    case "bilayer"
        masks = {struct('Type', 'bilayer', 'cs_inc', geometry.cs_high, ...
            'Params', struct('Bi_Angle', 0, ...
            'Bi_Offset', geometry.center(1), 'SigmaEdge', 1e-6))};
    otherwise
        error('Unsupported geometry type: %s', geometry.geometry_type);
end
cfg.MaskConfig = struct('cs_bg', geometry.cs_low, ...
    'CombineMode', 'overlay', 'Masks', {masks});
end

function [feat_cfg, req_options] = build_req_settings(CFG, M)
feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', M, 'cs_guess', CFG.REQ.cs_guess, ...
    'gamma_win', CFG.REQ.Gamma, 'pad_factor', CFG.REQ.PadFactor);
req_options = { ...
    'Nbins', CFG.REQ.Nbins, ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, ...
    'smooth_sigma', CFG.REQ.SmoothSigma};
end

function [T, OUT, global_req] = extract_feature_table( ...
    sim, cfg, feat_cfg, req_options, geometry, wave_case, CFG, condition_id)

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

cfg_req = cfg;
cfg_req.UseParfor = false;
[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz, cfg_req, req_cfg, feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features(global_curve);

global_req = struct('q', q_global, ...
    'mapping', adaptive_req.quantile.make_req_mapping(global_curve), ...
    'features', global_features, 'shape_features', global_shape);

OUT = adaptive_req.estimators.req_estimator_map( ...
    sim.Uxz, cfg_req, feat_cfg, ...
    'StepX', CFG.REQ.StepX, 'StepZ', CFG.REQ.StepZ, ...
    'EdgeMode', CFG.REQ.EdgeMode, ...
    'QuantileMode', 'local_req', ...
    'ReqOptions', req_options, ...
    'ReturnFeatures', true, ...
    'ReturnFeatureTable', true, ...
    'ReuseReqSpectrumForFeatures', true, ...
    'UseWindowParfor', true, ...
    'StoreReqCurves', false, ...
    'Verbose', false);

T = attach_feature_metadata(OUT.feature_table, cfg, feat_cfg, ...
    geometry, wave_case, condition_id);
n = height(T);
T.global_req_mapping = repmat({global_req.mapping}, n, 1);
T.q_global_req = q_global * ones(n, 1);
T.global_REQ_Nbins_effective = global_curve.Nbins_effective * ones(n, 1);
T = assign_global_feature_columns(T, global_features, "global_");
T = assign_global_feature_columns(T, global_shape, "global_");
end

function T = attach_feature_metadata(T, cfg, feat_cfg, geometry, wave_case, condition_id)
n = height(T);
cs_guess = get_cs_guess(feat_cfg);
[cs_true, region_name, material_region, distance_m] = geometry_at_points( ...
    T.x_center_m, T.z_center_m, geometry);

T.condition_id = condition_id * ones(n, 1);
T.condition_label = repmat(geometry.geometry_id + "__" + ...
    string(wave_case.wave_label) + "__M" + string(feat_cfg.M), n, 1);
T.geometry_id = repmat(string(geometry.geometry_id), n, 1);
T.geometry_name = repmat(string(geometry.geometry_name), n, 1);
T.geometry_type = repmat(string(geometry.geometry_type), n, 1);
T.wavefield_name = repmat(string(wave_case.case_name), n, 1);
T.wave_label = repmat(string(wave_case.wave_label), n, 1);
T.SIM_WaveModel = repmat(string(wave_case.case_name), n, 1);
T.step_idx = ones(n, 1);
T.realization_idx = ones(n, 1);
T.SIM_f0 = cfg.f0 * ones(n, 1);
T.SIM_cs_bg = geometry.cs_low * ones(n, 1);
T.SIM_cs_inc = geometry.cs_high * ones(n, 1);
T.SIM_Nwaves = cfg.Nwaves * ones(n, 1);
T.REQ_M = feat_cfg.M * ones(n, 1);
T.REQ_cs_guess = cs_guess * ones(n, 1);
T.M_eff_guess = feat_cfg.M * ones(n, 1);
T.lambda_guess = cs_guess / cfg.f0 * ones(n, 1);
T.cs_true = cs_true;
T.region_name = region_name;
T.material_region = material_region;
T.distance_to_interface_m = distance_m;
end

function [cs, region, material, distance_m] = geometry_at_points(x, z, geometry)
n = numel(x);
cs = geometry.cs_low * ones(n, 1);
region = strings(n, 1);
material = strings(n, 1);

switch geometry.geometry_type
    case "homogeneous"
        distance_m = nan(n, 1);
        region(:) = "homogeneous";
        material(:) = "homogeneous";

    case "inclusion"
        signed_d = hypot(x - geometry.center(1), z - geometry.center(2)) - geometry.radius;
        distance_m = abs(signed_d);
        inside = signed_d < 0;
        cs(inside) = geometry.cs_high;
        material(~inside) = "background";
        material(inside) = "inclusion";
        region(:) = "transition_excluded";
        region(signed_d <= -geometry.far_margin) = "inclusion_core";
        region(abs(signed_d) <= geometry.interface_half_width) = "interface_band";
        region(signed_d >= geometry.far_margin) = "background_far";

    case "bilayer"
        signed_d = x - geometry.center(1);
        distance_m = abs(signed_d);
        high = signed_d >= 0;
        cs(high) = geometry.cs_high;
        material(~high) = "layer_1";
        material(high) = "layer_2";
        region(:) = "transition_excluded";
        region(signed_d <= -geometry.far_margin) = "layer_1_far";
        region(signed_d >= geometry.far_margin) = "layer_2_far";
        region(abs(signed_d) <= geometry.interface_half_width) = "interface_band";
end
end

function q = compute_theory_q(cfg, feat_cfg, req_options, wave_case)
if string(wave_case.theory_field_type) == "Partial3D"
    q2 = theory_one(cfg, feat_cfg, req_options, "Diffuse2D");
    q3 = theory_one(cfg, feat_cfg, req_options, "Diffuse3D");
    q = 0.5 * (q2 + q3);
else
    q = theory_one(cfg, feat_cfg, req_options, ...
        string(wave_case.theory_field_type));
end
end

function q = theory_one(cfg, feat_cfg, req_options, field_type)
[req_cfg, ~] = adaptive_req.config.default_req_config(cfg, feat_cfg, req_options{:});
out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    cfg.dx, cfg.dz, cfg.f0, get_cs_guess(feat_cfg), ...
    'M', feat_cfg.M, 'Gamma', feat_cfg.gamma_win, ...
    'PadFactor', feat_cfg.pad_factor, 'Nbins', req_cfg.Nbins, ...
    'SmoothSigma', req_cfg.smooth_sigma, 'TheoryMode', 'S2D', ...
    'FieldType', field_type, 'Plot', false);
q = out.q_th;
end

function T = predict_operational_model(spec, T_feat, cfg)
Tq = adaptive_req.analysis.predict_q_model_from_table( ...
    spec.model, T_feat, 'ModelType', spec.model_type, ...
    'ModelName', spec.model_name);
T = prediction_base(T_feat);
T.method_name = repmat(spec.model_name, height(T), 1);
T.model_name = T.method_name;
T.feature_set = repmat(spec.feature_set, height(T), 1);
T.model_type = repmat(spec.model_type, height(T), 1);
T.model_role = repmat("operational", height(T), 1);
T.q_pred_raw = Tq.q_pred_raw;
T.q_pred = min(max(Tq.q_pred, 0.001), 0.999);
T.q_theory_used = nan(height(T), 1);
T.cs_pred = q_to_cs(T.q_pred, T_feat.req_mapping, cfg.f0);
end

function T = predict_theory_model(T_feat, cfg, q_theory, wave_case)
T = prediction_base(T_feat);
T.method_name = repmat("TheoryQDiscrete", height(T), 1);
T.model_name = T.method_name;
T.feature_set = repmat(string(wave_case.theory_label), height(T), 1);
T.model_type = repmat("theory_no_ml", height(T), 1);
T.model_role = repmat("diagnostic_only", height(T), 1);
T.q_pred_raw = q_theory * ones(height(T), 1);
T.q_pred = T.q_pred_raw;
T.q_theory_used = T.q_pred;
T.cs_pred = q_to_cs(T.q_pred, T_feat.req_mapping, cfg.f0);
end

function T = prediction_base(T_feat)
vars = ["condition_id", "condition_label", "geometry_id", ...
    "geometry_name", "geometry_type", "wavefield_name", "wave_label", ...
    "SIM_WaveModel", "SIM_f0", "SIM_cs_bg", "SIM_cs_inc", ...
    "SIM_Nwaves", "REQ_M", "REQ_cs_guess", "M_eff_guess", ...
    "patch_idx", "map_iz", "map_ix", "x_center_m", "z_center_m", ...
    "region_name", "material_region", "distance_to_interface_m", "cs_true"];
vars = vars(ismember(vars, string(T_feat.Properties.VariableNames)));
T = T_feat(:, cellstr(vars));
end

function cs = q_to_cs(q, mappings, f0)
cs = nan(numel(q), 1);
for i = 1:numel(q)
    if ~isempty(mappings{i}) && isfinite(q(i))
        cs(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), f0);
    end
end
end

function T = add_error_metrics(T)
T.cs_error = T.cs_pred - T.cs_true;
T.cs_error_pct = 100 * T.cs_error ./ T.cs_true;
T.cs_abs_error_pct = abs(T.cs_error_pct);
T.is_high_error_10 = T.cs_abs_error_pct > 10;
T.is_high_error_20 = T.cs_abs_error_pct > 20;
end

function T = add_map_gradients(T, CFG)
groups = unique(T(:, {'condition_id', 'method_name'}), 'rows');
T.q_gradient_mag = nan(height(T), 1);
T.cs_gradient_mag = nan(height(T), 1);
for i = 1:height(groups)
    idx = T.condition_id == groups.condition_id(i) & ...
        T.method_name == groups.method_name(i);
    Ti = T(idx, :);
    [Q, ~, ~] = map_from_table(Ti, 'q_pred');
    [C, ~, ~] = map_from_table(Ti, 'cs_pred');
    [dqz, dqx] = gradient(Q, CFG.SIM.dz, CFG.SIM.dx);
    [dcz, dcx] = gradient(C, CFG.SIM.dz, CFG.SIM.dx);
    T.q_gradient_mag(idx) = values_from_map(Ti, hypot(dqx, dqz));
    T.cs_gradient_mag(idx) = values_from_map(Ti, hypot(dcx, dcz));
end
end

function T_out = summarize_region_metrics(T)
T_named = T(T.region_name ~= "transition_excluded", :);
T_whole = T;
T_whole.region_name = repmat("whole_valid_map", height(T), 1);
T_in = concat_tables(T_named, T_whole);
group_vars = {'geometry_id', 'geometry_name', 'geometry_type', ...
    'wavefield_name', 'wave_label', 'REQ_M', 'method_name', ...
    'feature_set', 'model_type', 'model_role', 'region_name'};
[G, T_out] = findgroups(T_in(:, group_vars));
T_out.N = splitapply(@numel, T_in.cs_error_pct, G);
T_out.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), T_in.cs_error_pct, G);
T_out.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T_in.cs_error_pct, G);
T_out.CoV_pct = splitapply(@cov_pct, T_in.cs_pred, G);
T_out.median_error_pct = splitapply(@(x) median(x, 'omitnan'), T_in.cs_error_pct, G);
T_out.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T_in.cs_error_pct, G);
T_out.HighError_gt10_pct = splitapply(@(x) 100*mean(abs(x)>10, 'omitnan'), T_in.cs_error_pct, G);
T_out.HighError_gt20_pct = splitapply(@(x) 100*mean(abs(x)>20, 'omitnan'), T_in.cs_error_pct, G);
T_out.cs_pred_mean = splitapply(@(x) mean(x, 'omitnan'), T_in.cs_pred, G);
T_out.cs_pred_std = splitapply(@(x) std(x, 'omitnan'), T_in.cs_pred, G);
end

function T = summarize_case_metrics(T_pred, sim, geometry, wave_case, M)
T = table();
T.geometry_id = string(geometry.geometry_id);
T.geometry_name = string(geometry.geometry_name);
T.geometry_type = string(geometry.geometry_type);
T.wavefield_name = string(wave_case.case_name);
T.wave_label = string(wave_case.wave_label);
T.REQ_M = M;
T.Nwaves = wave_case.Nwaves;
T.Is2D = logical(wave_case.Is2D);
T.in_plane_fraction = in_plane_fraction(sim);
T.min_abs_uy = min(abs(sim.diag.waveDirs.uy));
T.field_abs_mean = mean(abs(sim.Uxz(:)), 'omitnan');
T.field_abs_std = std(abs(sim.Uxz(:)), 0, 'omitnan');
T.N_valid_pixels_per_method = height(T_pred) / numel(unique(T_pred.method_name));
end

function T = make_theory_row(geometry, wave_case, cfg, feat_cfg, q)
T = table(string(geometry.geometry_id), string(wave_case.case_name), ...
    string(wave_case.theory_label), string(wave_case.theory_field_type), ...
    cfg.f0, cfg.dx, cfg.dz, feat_cfg.M, get_cs_guess(feat_cfg), q, ...
    'VariableNames', {'geometry_id','wavefield_name','theory_label', ...
    'theory_field_type','SIM_f0','SIM_dx','SIM_dz','REQ_M', ...
    'REQ_cs_guess','q_theory_discrete'});
if string(wave_case.theory_field_type) == "Partial3D"
    T.notes = "diagnostic approximation: 0.5*(Diffuse2D + Diffuse3D)";
else
    T.notes = "direct discrete theory quantile";
end
end

function T = make_model_comparison(T_region)
T = T_region(T_region.region_name == "whole_valid_map", :);
T = sortrows(T, {'geometry_id','wavefield_name','REQ_M','MAPE_pct'});
end

function T_out = summarize_error_vs_distance(T)
T = T(isfinite(T.distance_to_interface_m), :);
if isempty(T)
    T_out = table();
    return;
end
edges_mm = [0 0.5 1 2 3 4 5 7.5 10 15 Inf];
bin_labels = strings(numel(edges_mm)-1, 1);
for i = 1:numel(bin_labels)
    bin_labels(i) = sprintf('[%g,%g)', edges_mm(i), edges_mm(i+1));
end
bin_idx = discretize(T.distance_to_interface_m * 1000, edges_mm);
T.distance_bin = repmat("outside", height(T), 1);
valid_bin = isfinite(bin_idx);
T.distance_bin(valid_bin) = bin_labels(bin_idx(valid_bin));
[G, T_out] = findgroups(T(:, {'geometry_id','wavefield_name','REQ_M', ...
    'method_name','model_role','distance_bin'}));
T_out.N = splitapply(@numel, T.cs_abs_error_pct, G);
T_out.distance_mm_median = splitapply(@(x) median(x, 'omitnan'), ...
    T.distance_to_interface_m * 1000, G);
T_out.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), T.cs_abs_error_pct, G);
T_out.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.cs_error_pct, G);
T_out.HighError_gt20_pct = splitapply(@(x) 100*mean(x>20, 'omitnan'), ...
    T.cs_abs_error_pct, G);
end

function plot_particle_velocity_overview(CASE_OUT, GEOMETRIES, CFG, fig_dir)
wave_names = string({CFG.CASES.case_name});
fig = figure('Color','w','Units','centimeters', ...
    'Position',[1 1 6.2*numel(wave_names) 5.5*numel(GEOMETRIES)]);
tl = tiledlayout(numel(GEOMETRIES), numel(wave_names), ...
    'TileSpacing','compact','Padding','compact');
for gi = 1:numel(GEOMETRIES)
    for wi = 1:numel(wave_names)
        idx = find_case(CASE_OUT, GEOMETRIES(gi).geometry_id, ...
            CFG.CASES(wi).wave_label, CFG.REQ.M_list(1));
        ax = nexttile(tl);
        sim = CASE_OUT(idx).sim;
        imagesc(ax, sim.x*100, sim.z*100, real(sim.Uxz));
        axis(ax,'image'); set(ax,'YDir','normal','FontSize',8);
        colorbar(ax); colormap(ax, parula);
        draw_interface(ax, GEOMETRIES(gi));
        title(ax, sprintf('%s | %s', GEOMETRIES(gi).geometry_name, wave_names(wi)), ...
            'Interpreter','none','FontWeight','normal','FontSize',8);
        xlabel(ax,'x (cm)'); ylabel(ax,'z (cm)');
    end
end
title(tl,'Particle velocity: real component used by REQ','FontWeight','normal');
export_clean(fig, fullfile(fig_dir,'level17_compare_particle_velocity_maps.png'), CFG);
end

function plot_power_spectrum_overview(CASE_OUT, GEOMETRIES, CFG, fig_dir)
wave_names = string({CFG.CASES.case_name});
M_ref = 3;
fig = figure('Color','w','Units','centimeters', ...
    'Position',[1 1 6.2*numel(wave_names) 5.5*numel(GEOMETRIES)]);
tl = tiledlayout(numel(GEOMETRIES), numel(wave_names), ...
    'TileSpacing','compact','Padding','compact');
for gi = 1:numel(GEOMETRIES)
    for wi = 1:numel(wave_names)
        idx = find_case(CASE_OUT, GEOMETRIES(gi).geometry_id, ...
            CFG.CASES(wi).wave_label, M_ref);
        sim = CASE_OUT(idx).sim;
        patch = center_patch(sim.Uxz, CASE_OUT(idx).req_out.win_size);
        S = log1p(abs(fftshift(fft2(patch))).^2);
        ax = nexttile(tl);
        imagesc(ax, S); axis(ax,'image'); set(ax,'YDir','normal','FontSize',8);
        colorbar(ax); colormap(ax, turbo);
        title(ax, sprintf('%s | %s', GEOMETRIES(gi).geometry_name, wave_names(wi)), ...
            'Interpreter','none','FontWeight','normal','FontSize',8);
        xlabel(ax,'k_x bin'); ylabel(ax,'k_z bin');
    end
end
title(tl,'Representative central-window spectra, M=3','FontWeight','normal');
export_clean(fig, fullfile(fig_dir,'level17_compare_power_spectra.png'), CFG);
end

function plot_all_comparison_maps(T, GEOMETRIES, CFG, fig_dir)
value_specs = [
    struct('var',"cs_pred",'label',"c_s (m/s)",'suffix',"sws",'title',"SWS")
    struct('var',"cs_error_pct",'label',"signed error (%)",'suffix',"signed_error",'title',"Signed SWS error")
    struct('var',"cs_abs_error_pct",'label',"absolute error (%)",'suffix',"absolute_error",'title',"Absolute SWS error")
    struct('var',"q_pred",'label',"q",'suffix',"q",'title',"Predicted/used q")
    ];
for gi = 1:numel(GEOMETRIES)
    for wi = 1:numel(CFG.CASES)
        Tc = T(T.geometry_id == GEOMETRIES(gi).geometry_id & ...
            T.wave_label == string(CFG.CASES(wi).wave_label), :);
        for vi = 1:numel(value_specs)
            file_name = sprintf('level17_compare_%s_%s_%s_maps.png', ...
                GEOMETRIES(gi).geometry_id, ...
                sanitize_filename(CFG.CASES(wi).wave_label), ...
                value_specs(vi).suffix);
            plot_one_map_grid(Tc, GEOMETRIES(gi), value_specs(vi), ...
                CFG, fullfile(fig_dir, file_name));
        end
    end
end
end

function plot_one_map_grid(T, geometry, value_spec, CFG, file_path)
methods = ["LocalOnly","GlobalOnly","HybridLocalGlobal","TheoryQDiscrete"];
M_values = CFG.REQ.M_list;
fig = figure('Color','w','Units','centimeters', ...
    'Position',[1 1 7*numel(methods) 5.8*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(methods), ...
    'TileSpacing','compact','Padding','compact');
for mi = 1:numel(M_values)
    for model_i = 1:numel(methods)
        Ti = T(T.REQ_M == M_values(mi) & T.method_name == methods(model_i), :);
        ax = nexttile(tl);
        [A,x_cm,z_cm] = map_from_table(Ti, value_spec.var);
        imagesc(ax,x_cm,z_cm,A); axis(ax,'image');
        set(ax,'YDir','normal','FontSize',8);
        if contains(value_spec.var,"error")
            colormap(ax, parula);
        else
            colormap(ax, turbo);
        end
        cb = colorbar(ax); ylabel(cb,value_spec.label,'FontSize',8);
        draw_interface(ax,geometry);
        title(ax,sprintf('%s | M=%g | MAPE %.2f%%', ...
            methods(model_i),M_values(mi),mean(Ti.cs_abs_error_pct,'omitnan')), ...
            'Interpreter','none','FontWeight','normal','FontSize',8);
        xlabel(ax,'x (cm)'); ylabel(ax,'z (cm)');
    end
end
title(tl,sprintf('%s | %s | %s',geometry.geometry_name, ...
    T.wavefield_name(1),value_spec.title), ...
    'Interpreter','none','FontWeight','normal','FontSize',11);
export_clean(fig,file_path,CFG);
end

function plot_metric_bars(T_region, metric_var, y_label, file_prefix, CFG, fig_dir)
T = T_region(T_region.region_name == "whole_valid_map", :);
geometries = unique(T.geometry_id,'stable');
methods = ["LocalOnly","GlobalOnly","HybridLocalGlobal","TheoryQDiscrete"];
for gi = 1:numel(geometries)
    Tg = T(T.geometry_id == geometries(gi), :);
    waves = unique(Tg.wavefield_name,'stable');
    M_values = unique(Tg.REQ_M,'stable');
    fig = figure('Color','w','Units','centimeters', ...
        'Position',[1 1 8*numel(waves) 11]);
    tl = tiledlayout(1,numel(waves),'TileSpacing','compact','Padding','compact');
    for wi = 1:numel(waves)
        ax = nexttile(tl);
        Y = nan(numel(M_values),numel(methods));
        for mi = 1:numel(M_values)
            for model_i = 1:numel(methods)
                idx = Tg.wavefield_name == waves(wi) & ...
                    Tg.REQ_M == M_values(mi) & Tg.method_name == methods(model_i);
                if any(idx)
                    Y(mi,model_i) = Tg.(char(metric_var))(find(idx,1));
                end
            end
        end
        bar(ax,M_values,Y,'grouped'); grid(ax,'on');
        title(ax,waves(wi),'Interpreter','none','FontWeight','normal');
        xlabel(ax,'REQ M'); ylabel(ax,y_label); xticks(ax,M_values);
    end
    lg = legend(tl.Children(end),methods,'Location','northoutside', ...
        'Orientation','horizontal','Interpreter','none');
    lg.Layout.Tile = 'north';
    title(tl,sprintf('%s | %s',Tg.geometry_name(1),y_label), ...
        'Interpreter','none','FontWeight','normal');
    export_clean(fig,fullfile(fig_dir,sprintf('%s_%s.png', ...
        file_prefix,geometries(gi))),CFG);
end
end

function plot_error_vs_distance(T_distance, GEOMETRIES, CFG, fig_dir)
for gi = 1:numel(GEOMETRIES)
    if GEOMETRIES(gi).geometry_type == "homogeneous"
        continue;
    end
    T = T_distance(T_distance.geometry_id == GEOMETRIES(gi).geometry_id, :);
    waves = unique(T.wavefield_name,'stable');
    methods = ["LocalOnly","GlobalOnly","HybridLocalGlobal","TheoryQDiscrete"];
    fig = figure('Color','w','Units','centimeters', ...
        'Position',[1 1 8*numel(waves) 12]);
    tl = tiledlayout(1,numel(waves),'TileSpacing','compact','Padding','compact');
    for wi = 1:numel(waves)
        ax = nexttile(tl); hold(ax,'on');
        colors = lines(numel(methods));
        for model_i = 1:numel(methods)
            Ti = T(T.wavefield_name == waves(wi) & ...
                T.method_name == methods(model_i), :);
            [x,ord] = sort(Ti.distance_mm_median);
            y = Ti.MAPE_pct(ord);
            plot(ax,x,y,'-o','LineWidth',1.2,'MarkerSize',3, ...
                'Color',colors(model_i,:));
        end
        yline(ax,20,'k--'); grid(ax,'on');
        title(ax,waves(wi),'Interpreter','none','FontWeight','normal');
        xlabel(ax,'distance to interface (mm)'); ylabel(ax,'MAPE (%)');
    end
    lg = legend(tl.Children(end),methods,'Location','northoutside', ...
        'Orientation','horizontal','Interpreter','none');
    lg.Layout.Tile = 'north';
    title(tl,GEOMETRIES(gi).geometry_name + " | error vs interface distance", ...
        'Interpreter','none','FontWeight','normal');
    export_clean(fig,fullfile(fig_dir,sprintf( ...
        'level17_compare_error_vs_distance_%s.png', ...
        GEOMETRIES(gi).geometry_id)),CFG);
end
end

function plot_error_vs_q_gradient(T, CFG, fig_dir)
T = T(T.model_role == "operational" & isfinite(T.q_gradient_mag), :);
methods = ["LocalOnly","GlobalOnly","HybridLocalGlobal"];
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 9]);
tl = tiledlayout(1,numel(methods),'TileSpacing','compact','Padding','compact');
rng(1717);
for i = 1:numel(methods)
    ax = nexttile(tl);
    idx = find(T.method_name == methods(i) & isfinite(T.cs_abs_error_pct));
    if numel(idx) > 30000
        idx = idx(randperm(numel(idx),30000));
    end
    scatter(ax,T.q_gradient_mag(idx),T.cs_abs_error_pct(idx),5,'filled', ...
        'MarkerFaceAlpha',0.12);
    yline(ax,20,'k--'); grid(ax,'on');
    title(ax,methods(i),'Interpreter','none','FontWeight','normal');
    xlabel(ax,'|grad q|'); ylabel(ax,'|SWS error| (%)');
end
title(tl,'Operational-model error vs q-map gradient','FontWeight','normal');
export_clean(fig,fullfile(fig_dir, ...
    'level17_compare_error_vs_q_gradient.png'),CFG);
end

function print_summary(T_region)
T = T_region(T_region.region_name == "whole_valid_map", :);
[G,Tmean] = findgroups(T(:,{'geometry_type','method_name','model_role'}));
Tmean.MAPE_pct = splitapply(@(x) mean(x,'omitnan'),T.MAPE_pct,G);
Tmean.HighError_gt20_pct = splitapply(@(x) mean(x,'omitnan'), ...
    T.HighError_gt20_pct,G);

fprintf('\nMean whole-map performance:\n');
disp(sortrows(Tmean,{'geometry_type','MAPE_pct'}));

geometry_types = ["homogeneous","bilayer","inclusion"];
fprintf('\nBest method by geometry (whole-map mean MAPE):\n');
for i = 1:numel(geometry_types)
    Ti = Tmean(Tmean.geometry_type == geometry_types(i),:);
    [~,idx] = min(Ti.MAPE_pct);
    fprintf('  %-12s: %-22s MAPE %.2f%% | high-error >20 %.2f%%\n', ...
        geometry_types(i),Ti.method_name(idx),Ti.MAPE_pct(idx), ...
        Ti.HighError_gt20_pct(idx));
end

Top = T_region(T_region.region_name == "interface_band",:);
Far = T_region(contains(T_region.region_name,"far"),:);
fprintf('\nInterface diagnostic:\n');
for method = unique(Top.method_name,'stable')'
    a = mean(Top.MAPE_pct(Top.method_name == method),'omitnan');
    b = mean(Far.MAPE_pct(Far.method_name == method),'omitnan');
    fprintf('  %-22s interface %.2f%% vs far-region %.2f%% (ratio %.2f)\n', ...
        method,a,b,a/max(b,eps));
end

Operational = Tmean(Tmean.model_role == "operational",:);
fprintf('\nOperational-model interpretation:\n');
for geometry_type = geometry_types
    Ti = Operational(Operational.geometry_type == geometry_type,:);
    [~,idx] = min(Ti.MAPE_pct);
    fprintf('  Best for %s: %s.\n',geometry_type,Ti.method_name(idx));
end

Tint = mean(Top.MAPE_pct,'omitnan');
Tfar = mean(Far.MAPE_pct,'omitnan');
if Tint > 1.5*Tfar
    fprintf(['Conclusion: errors increase strongly near interfaces, so ', ...
        'heterogeneity/window mixing is an important failure mode.\n']);
else
    fprintf(['Conclusion: errors are not strongly isolated to interfaces; ', ...
        'wavefield type and model transfer should be inspected jointly.\n']);
end

[G2,Twave] = findgroups(T(:,{'wavefield_name','method_name'}));
Twave.MAPE_pct = splitapply(@(x) mean(x,'omitnan'),T.MAPE_pct,G2);
Twave_summary = groupsummary(Twave,'wavefield_name','mean','MAPE_pct');
spread_by_wave = max(Twave_summary.mean_MAPE_pct) - ...
    min(Twave_summary.mean_MAPE_pct);
fprintf('Mean performance spread across wavefield types: %.2f percentage points.\n', ...
    spread_by_wave);
fprintf(['TheoryQDiscrete is diagnostic_only. No q target, true SWS, or ', ...
    'oracle effective-window variable was used by operational predictors.\n']);
end

function idx = find_case(CASE_OUT, geometry_id, wave_label, M)
idx = find(arrayfun(@(s) string(s.geometry.geometry_id) == string(geometry_id) && ...
    string(s.wave_case.wave_label) == string(wave_label) && s.REQ_M == M, ...
    CASE_OUT),1);
assert(~isempty(idx),'Missing CASE_OUT entry for %s | %s | M=%g.', ...
    geometry_id,wave_label,M);
end

function patch = center_patch(U, win)
[nz,nx] = size(U);
cz = round((nz+1)/2); cx = round((nx+1)/2);
half = floor(win/2);
patch = U((cz-half):(cz+half),(cx-half):(cx+half));
end

function draw_interface(ax, geometry)
hold(ax,'on');
switch geometry.geometry_type
    case "inclusion"
        th = linspace(0,2*pi,300);
        plot(ax,(geometry.center(1)+geometry.radius*cos(th))*100, ...
            (geometry.center(2)+geometry.radius*sin(th))*100, ...
            'w--','LineWidth',1);
    case "bilayer"
        xline(ax,geometry.center(1)*100,'w--','LineWidth',1);
end
hold(ax,'off');
end

function [A,x_cm,z_cm] = map_from_table(T,var_name)
nz = max(T.map_iz); nx = max(T.map_ix);
A = nan(nz,nx);
A(sub2ind([nz nx],T.map_iz,T.map_ix)) = T.(char(var_name));
x_cm = unique(T.x_center_m,'stable')*100;
z_cm = unique(T.z_center_m,'stable')*100;
end

function vals = values_from_map(T,A)
vals = A(sub2ind(size(A),T.map_iz,T.map_ix));
end

function y = cov_pct(x)
x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100*std(x)/mean(x);
end
end

function cs_guess = get_cs_guess(feat_cfg)
if isfield(feat_cfg,'cs_guess_used')
    cs_guess = feat_cfg.cs_guess_used;
elseif isfield(feat_cfg,'cs_guess')
    cs_guess = feat_cfg.cs_guess;
else
    error('Feature config does not contain cs_guess.');
end
end

function T = assign_global_feature_columns(T,values,prefix)
if ~isstruct(values), return; end
names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if (isnumeric(value_i) || islogical(value_i)) && isscalar(value_i)
        T.(char(string(prefix)+string(names{i}))) = double(value_i)*ones(height(T),1);
    end
end
end

function T = concat_tables(A,B)
if isempty(A), T = B; return; end
if isempty(B), T = A; return; end
vars = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)],'stable');
A = add_missing_columns(A,vars);
B = add_missing_columns(B,vars);
T = [A(:,cellstr(vars)); B(:,cellstr(vars))];
end

function T = add_missing_columns(T,vars)
for i = 1:numel(vars)
    if ~ismember(vars(i),string(T.Properties.VariableNames))
        T.(char(vars(i))) = nan(height(T),1);
    end
end
end

function T = remove_cell_columns(T)
remove = false(1,width(T));
for i = 1:width(T)
    remove(i) = iscell(T.(T.Properties.VariableNames{i}));
end
T(:,remove) = [];
end

function print_in_plane_diagnostic(sim,wave_case)
fprintf('In-plane fraction %.3f | min |u_y| %.3g\n', ...
    in_plane_fraction(sim),min(abs(sim.diag.waveDirs.uy)));
if string(wave_case.wave_label) == "partial_diffuse_3d"
    assert(min(abs(sim.diag.waveDirs.uy)) < 1e-7, ...
        'Partially diffuse 3D has no explicitly in-plane source.');
end
end

function frac = in_plane_fraction(sim)
covg = sim.diag.inPlaneCoverage;
if isfield(covg,'fraction_in_plane')
    frac = covg.fraction_in_plane;
elseif isfield(covg,'count_in_plane') && isfield(covg,'n_waves') && covg.n_waves > 0
    frac = covg.count_in_plane/covg.n_waves;
else
    frac = NaN;
end
end

function name = sanitize_filename(x)
name = regexprep(char(string(x)),'[^A-Za-z0-9_=-]+','_');
end

function export_clean(fig,file_path,CFG)
axs = findall(fig,'Type','axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end
drawnow;
exportgraphics(fig,file_path,'Resolution',CFG.FIG.Resolution, ...
    'BackgroundColor','white');
close(fig);
end
