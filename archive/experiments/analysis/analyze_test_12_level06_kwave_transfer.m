%% analyze_test_12_level06_kwave_transfer.m
% Test 12 Level 06: transfer the trained adaptive-q REQ model to k-Wave data.
%
% The k-Wave data are 500 Hz inclusion phantoms with a 2 m/s background and
% a 3 m/s circular inclusion. This script loads the complex harmonic fields,
% extracts local REQ mappings/features, predicts q with the deployed Test 12
% HybridLocalGlobal + WithCsGuess model, converts q_pred to SWS, and plots
% predicted SWS/error/q maps across k-Wave cases.

clear; clc; close all;
format compact;

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLegendFontSize', 11);

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

%% Analysis settings

KW = struct();
KW.f0 = 500;
KW.cs_bg = 2.0;
KW.cs_inc = 3.0;
KW.inclusion_center_m = [0.025, 0.025];
KW.inclusion_radius_m = 0.010;

KW.REQ.M_list = [2 3 4];
KW.REQ.M = KW.REQ.M_list(1);
KW.REQ.cs_guess = 3.0;
KW.REQ.StepX = 1;
KW.REQ.StepZ = 1;
KW.REQ.Gamma = 1;
KW.REQ.PadFactor = 1;
KW.REQ.Nbins = 'auto';
KW.REQ.Nbins_auto_oversample = 1;
KW.REQ.Nbins_min = 16;
KW.REQ.SmoothSigma = 1;
KW.REQ.EdgeMode = 'valid';

KW.MODEL.ModelName = "HybridLocalGlobal";
KW.MODEL.FeatureSet = "WithCsGuess";
KW.MODEL.ModelType = "bagged_trees";

KW.ROI(1).name = "inclusion_roi";
KW.ROI(1).true_cs = KW.cs_inc;
KW.ROI(1).xlim_m = [0.021, 0.029];
KW.ROI(1).zlim_m = [0.021, 0.029];

KW.ROI(2).name = "background_roi";
KW.ROI(2).true_cs = KW.cs_bg;
KW.ROI(2).xlim_m = [0.006, 0.014];
KW.ROI(2).zlim_m = [0.021, 0.029];

KW.cases = [
    struct('case_name', "Directional 2D", 'folder', "2D-SS", ...
        'wave_model', "SingleWave")
    struct('case_name', "Diffuse 2D", 'folder', "2D-diffuse", ...
        'wave_model', "Diffuse2D")
    struct('case_name', "Projected diffuse 3D", 'folder', "3D-diffuse", ...
        'wave_model', "Diffuse3D")
    struct('case_name', "3D-rev", 'folder', "3D-rev", ...
        'wave_model', "Rev3D")
    ];

%% Locate Test 12 run and output folders

PATHS12 = locate_latest_test12_paths(root_dir);
fprintf('\nUsing Test 12 run folder:\n%s\n', PATHS12.run_dir);

analysis_dir = fullfile(PATHS12.analysis_dir, 'level_06_kwave_transfer');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
data_dir = fullfile(analysis_dir, 'data');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(data_dir);

%% Load deployed model

deployment_dir = fullfile(PATHS12.analysis_dir, ...
    'level_01_model_comparison', 'models', 'deployment');
[MODEL_DEPLOY, MODEL_INFO, model_file] = ...
    adaptive_req.analysis.load_q_model_deployment( ...
    deployment_dir, ...
    'ModelName', KW.MODEL.ModelName, ...
    'FeatureSet', KW.MODEL.FeatureSet, ...
    'ModelType', KW.MODEL.ModelType);

fprintf('\nLoaded deployment model:\n%s\n', model_file);
disp(struct2table(MODEL_INFO));

%% Main k-Wave loop

T_all_pred = table();
T_roi = table();
CASE_OUT = struct([]);
out_idx = 0;

for ci = 1:numel(KW.cases)
    case_i = KW.cases(ci);
    data_file = fullfile(root_dir, 'data', 'k-wave', ...
        char(case_i.folder), sprintf('data_%dHz.mat', KW.f0));

    fprintf('\n=== k-Wave case %d / %d: %s ===\n', ...
        ci, numel(KW.cases), case_i.case_name);
    fprintf('Loading %s\n', data_file);

    if ~exist(data_file, 'file')
        warning('Skipping missing k-Wave data file:\n%s', data_file);
        continue;
    end

    S = load(data_file, 'Vz_mg_ph', 'dinf');
    field = S.Vz_mg_ph;
    dinf = S.dinf;

    cfg = struct();
    cfg.f0 = KW.f0;
    cfg.cs_bg = KW.cs_bg;
    cfg.dx = dinf.dx;
    cfg.dz = dinf.dy;
    cfg.UseParfor = false;

    fprintf('Field size: %s | dx=%.4g m | dz=%.4g m\n', ...
        mat2str(size(field)), cfg.dx, cfg.dz);

    for mi = 1:numel(KW.REQ.M_list)
        KW_i = KW;
        KW_i.REQ.M = KW.REQ.M_list(mi);
        condition_id = (ci - 1) * numel(KW.REQ.M_list) + mi;

        [feat_cfg, req_options] = build_req_settings(KW_i);

        fprintf('Computing local REQ mappings/features: M=%g, cs_guess=%g, Step=(%d,%d), EdgeMode=%s...\n', ...
            KW_i.REQ.M, KW_i.REQ.cs_guess, KW_i.REQ.StepX, KW_i.REQ.StepZ, KW_i.REQ.EdgeMode);
        t_req = tic;
        [T_feat, req_out, global_req] = extract_kwave_feature_table( ...
            field, cfg, feat_cfg, req_options, case_i, KW_i, condition_id);
        fprintf('REQ/features completed in %.2f s for %d windows.\n', ...
            toc(t_req), height(T_feat));

        fprintf('Applying adaptive-q model...\n');
        t_pred = tic;
        T_pred_i = predict_kwave_map(MODEL_DEPLOY, MODEL_INFO, T_feat, cfg);
        T_pred_i = add_error_metrics(T_pred_i);
        fprintf('Prediction completed in %.2f s.\n', toc(t_pred));

        T_roi_i = summarize_roi_metrics(T_pred_i);

        T_all_pred = concat_tables(T_all_pred, T_pred_i);
        T_roi = concat_tables(T_roi, T_roi_i);

        out_idx = out_idx + 1;
        CASE_OUT(out_idx).case = case_i;
        CASE_OUT(out_idx).REQ_M = KW_i.REQ.M;
        CASE_OUT(out_idx).cfg = cfg;
        CASE_OUT(out_idx).dinf = dinf;
        CASE_OUT(out_idx).req_out = req_out;
        CASE_OUT(out_idx).global_req = global_req;
        CASE_OUT(out_idx).T_feat = T_feat;
        CASE_OUT(out_idx).T_pred = T_pred_i;
        CASE_OUT(out_idx).T_roi = T_roi_i;
    end
end

assert(~isempty(T_all_pred), 'No k-Wave predictions were generated.');

%% Save tables/data

writetable(remove_cell_columns(T_all_pred), fullfile(table_dir, ...
    'level12_level06_kwave_predictions.csv'));
writetable(T_roi, fullfile(table_dir, ...
    'level12_level06_kwave_roi_metrics.csv'));

save(fullfile(data_dir, 'level12_level06_kwave_transfer.mat'), ...
    'KW', 'MODEL_INFO', 'model_file', 'CASE_OUT', ...
    'T_all_pred', 'T_roi', '-v7.3');

%% Figures

plot_kwave_sws_panel(T_all_pred, KW, fig_dir);
plot_kwave_error_panel(T_all_pred, KW, fig_dir);
plot_kwave_q_panel(T_all_pred, KW, fig_dir);
plot_kwave_roi_summary(T_roi, fig_dir);

fprintf('\nTest 12 Level 06 k-Wave transfer complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);
fprintf('\nROI metrics:\n');
disp(T_roi(:, {'case_name', 'REQ_M', 'roi_name', 'true_cs', ...
    'MAPE_pct', 'CoV_pct', 'cs_pred_mean', 'cs_pred_std', 'n'}));

%% Local functions

function PATHS = locate_latest_test12_paths(root_dir)

out_root = fullfile(root_dir, 'outputs', 'test_12_cs_guess_window_sweep');
runs = dir(fullfile(out_root, 'test_12_cs_guess_window_sweep_*'));
runs = runs([runs.isdir]);
assert(~isempty(runs), 'No Test 12 run folders found in %s', out_root);

[~, idx] = max([runs.datenum]);
run_dir = fullfile(runs(idx).folder, runs(idx).name);

PATHS = struct();
PATHS.run_dir = run_dir;
PATHS.analysis_dir = fullfile(run_dir, 'analysis');

end

function [feat_cfg, req_options] = build_req_settings(KW)

R = KW.REQ;
feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', R.M, ...
    'cs_guess', R.cs_guess, ...
    'gamma_win', R.Gamma, ...
    'pad_factor', R.PadFactor);

req_options = { ...
    'Nbins', R.Nbins, ...
    'Nbins_auto_oversample', R.Nbins_auto_oversample, ...
    'Nbins_min', R.Nbins_min, ...
    'smooth_sigma', R.SmoothSigma};

end

function [T, OUT, global_req] = extract_kwave_feature_table( ...
    field, cfg, feat_cfg, req_options, case_i, KW, condition_id)

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
    field, cfg, req_cfg, feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features(global_curve);

global_req = struct();
global_req.q = q_global;
global_req.mapping = adaptive_req.quantile.make_req_mapping(global_curve);
global_req.features = global_features;
global_req.shape_features = global_shape;

OUT = adaptive_req.estimators.req_estimator_map( ...
    field, cfg, feat_cfg, ...
    'StepX', KW.REQ.StepX, ...
    'StepZ', KW.REQ.StepZ, ...
    'EdgeMode', KW.REQ.EdgeMode, ...
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
T.SIM_WaveModel = repmat(case_i.wave_model, height(T), 1);
T.step_idx = ones(height(T), 1);
T.realization_idx = ones(height(T), 1);
T.frequency_hz = cfg.f0 * ones(height(T), 1);
T.SIM_f0 = cfg.f0 * ones(height(T), 1);
T.SIM_cs_bg = KW.cs_bg * ones(height(T), 1);
T.SIM_cs_inc = KW.cs_inc * ones(height(T), 1);
T.SIM_SNR_dB = Inf(height(T), 1);
T.SNR_label = repmat("kWave", height(T), 1);
T.REQ_M = feat_cfg.M * ones(height(T), 1);
T.REQ_cs_guess = feat_cfg.cs_guess_used * ones(height(T), 1);
T.M_eff_guess = feat_cfg.M * ones(height(T), 1);
T.M_eff_true_diag = feat_cfg.M * feat_cfg.cs_guess_used ./ ...
    KW.cs_bg * ones(height(T), 1);
T.lambda_guess = feat_cfg.cs_guess_used / cfg.f0 * ones(height(T), 1);
T.lambda_true_bg = KW.cs_bg / cfg.f0 * ones(height(T), 1);
T.lambda_true_inc = KW.cs_inc / cfg.f0 * ones(height(T), 1);
T.aperture_label = repmat("kWave", height(T), 1);
T.ConeHalfAngleDeg = NaN(height(T), 1);
T.Omega_sr = NaN(height(T), 1);
T.omega_mean = T.Omega_sr;

T.cs_true = true_cs_at_points(T.x_center_m, T.z_center_m, KW);
T.roi_name = classify_kwave_roi(T.x_center_m, T.z_center_m, KW.ROI);
T.global_req_mapping = repmat({global_req.mapping}, height(T), 1);
T.q_global_req = q_global * ones(height(T), 1);
T.global_REQ_Nbins_effective = global_curve.Nbins_effective * ones(height(T), 1);

T = assign_global_feature_columns(T, global_features, "global_");
T = assign_global_feature_columns(T, global_shape, "global_");

end

function cs = true_cs_at_points(x, z, KW)

r = hypot(x - KW.inclusion_center_m(1), z - KW.inclusion_center_m(2));
cs = KW.cs_bg * ones(numel(x), 1);
cs(r <= KW.inclusion_radius_m) = KW.cs_inc;

end

function roi = classify_kwave_roi(x, z, rois)

roi = repmat("outside_roi", numel(x), 1);
for ri = 1:numel(rois)
    in_roi = x >= rois(ri).xlim_m(1) & x <= rois(ri).xlim_m(2) & ...
        z >= rois(ri).zlim_m(1) & z <= rois(ri).zlim_m(2);
    roi(in_roi) = string(rois(ri).name);
end

end

function T_pred = predict_kwave_map(MODEL, MODEL_INFO, T_feat, cfg)

T_q = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL, T_feat, ...
    'ModelType', MODEL_INFO.model_type, ...
    'ModelName', MODEL_INFO.model_name);

vars = ["condition_id", "condition_label", "case_name", ...
    "step_idx", "realization_idx", "frequency_hz", "SIM_f0", ...
    "SIM_cs_bg", "SIM_cs_inc", "SIM_WaveModel", "SIM_SNR_dB", ...
    "SNR_label", "REQ_M", "REQ_cs_guess", "M_eff_guess", ...
    "M_eff_true_diag", "aperture_label", "ConeHalfAngleDeg", ...
    "Omega_sr", "patch_idx", "map_iz", "map_ix", ...
    "x_center_m", "z_center_m", "roi_name", "cs_true"];
vars = vars(ismember(vars, string(T_feat.Properties.VariableNames)));

T_pred = T_feat(:, cellstr(vars));
T_pred.model_name = repmat(string(MODEL_INFO.model_name), height(T_feat), 1);
T_pred.feature_set = repmat(string(MODEL_INFO.feature_set), height(T_feat), 1);
T_pred.model_type = repmat(string(MODEL_INFO.model_type), height(T_feat), 1);
T_pred.model_role = repmat("kwave_transfer", height(T_feat), 1);
T_pred.q_pred_raw = T_q.q_pred_raw;
T_pred.q_pred = min(max(T_q.q_pred, 0.001), 0.999);
T_pred.cs_pred = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, cfg.f0);

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

function T_roi = summarize_roi_metrics(T)

T = T(T.roi_name ~= "outside_roi", :);
if isempty(T)
    T_roi = table();
    return;
end

[G, T_roi] = findgroups(T(:, {'case_name', 'model_name', ...
    'feature_set', 'model_type', 'roi_name', 'SIM_WaveModel', ...
    'SIM_f0', 'REQ_M', 'REQ_cs_guess'}));
T_roi.n = splitapply(@numel, T.cs_error_pct, G);
T_roi.true_cs = splitapply(@(x) median(x, 'omitnan'), T.cs_true, G);
T_roi.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), ...
    T.cs_error_pct, G);
T_roi.bias_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.cs_error_pct, G);
T_roi.CoV_pct = splitapply(@cov_pct, T.cs_pred, G);
T_roi.cs_pred_mean = splitapply(@(x) mean(x, 'omitnan'), T.cs_pred, G);
T_roi.cs_pred_median = splitapply(@(x) median(x, 'omitnan'), T.cs_pred, G);
T_roi.cs_pred_std = splitapply(@(x) std(x, 'omitnan'), T.cs_pred, G);

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function plot_kwave_sws_panel(T, KW, fig_dir)

cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.2*numel(cases) 6.2*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Tc = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile(tl);
        [A, x_cm, z_cm] = map_from_table(Tc, 'cs_pred');
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontSize', 9);
        colormap(ax, turbo);
        colorbar(ax);
        draw_geometry(ax, KW);
        mape = mean(Tc.cs_abs_error_pct, 'omitnan');
        title(ax, sprintf('%s | M=%g | MAPE %.2f%%', ...
            cases(ci), M_values(mi), mape), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
        xlabel(ax, 'x (cm)');
        ylabel(ax, 'z (cm)');
    end
end

title(tl, sprintf('k-Wave adaptive-q REQ SWS maps | c_{s,guess}=%g m/s', ...
    KW.REQ.cs_guess), 'FontWeight', 'normal', 'FontSize', 12);
exportgraphics(fig, fullfile(fig_dir, ...
    'level12_level06_kwave_sws_maps.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function plot_kwave_error_panel(T, KW, fig_dir)

cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.2*numel(cases) 6.2*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Tc = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile(tl);
        [A, x_cm, z_cm] = map_from_table(Tc, 'cs_error_pct');
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontSize', 9);
        colormap(ax, parula);
        colorbar(ax);
        draw_geometry(ax, KW);
        title(ax, sprintf('%s | M=%g', cases(ci), M_values(mi)), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
        xlabel(ax, 'x (cm)');
        ylabel(ax, 'z (cm)');
    end
end

title(tl, 'k-Wave signed SWS error maps (%)', ...
    'FontWeight', 'normal', 'FontSize', 12);
exportgraphics(fig, fullfile(fig_dir, ...
    'level12_level06_kwave_error_maps.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function plot_kwave_q_panel(T, KW, fig_dir)

cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.2*numel(cases) 6.2*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Tc = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile(tl);
        [A, x_cm, z_cm] = map_from_table(Tc, 'q_pred');
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontSize', 9);
        colormap(ax, turbo);
        colorbar(ax);
        draw_geometry(ax, KW);
        title(ax, sprintf('%s | M=%g | median q %.3f', ...
            cases(ci), M_values(mi), median(Tc.q_pred, 'omitnan')), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
        xlabel(ax, 'x (cm)');
        ylabel(ax, 'z (cm)');
    end
end

title(tl, 'k-Wave predicted q maps', 'FontWeight', 'normal', 'FontSize', 12);
exportgraphics(fig, fullfile(fig_dir, ...
    'level12_level06_kwave_q_maps.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function plot_kwave_roi_summary(T_roi, fig_dir)

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8*numel(unique(T_roi.case_name, 'stable')) 12]);
tl = tiledlayout(2, numel(unique(T_roi.case_name, 'stable')), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

cases = unique(T_roi.case_name, 'stable');
M_values = unique(T_roi.REQ_M, 'stable');
roi_order = ["background_roi", "inclusion_roi"];
x = categorical(string(M_values));
x = reordercats(x, cellstr(string(M_values)));

for ci = 1:numel(cases)
    mape = nan(numel(M_values), numel(roi_order));
    covv = nan(numel(M_values), numel(roi_order));

    for mi = 1:numel(M_values)
        for ri = 1:numel(roi_order)
            idx = T_roi.case_name == cases(ci) & ...
                T_roi.REQ_M == M_values(mi) & ...
                T_roi.roi_name == roi_order(ri);
            if any(idx)
                mape(mi, ri) = T_roi.MAPE_pct(find(idx, 1, 'first'));
                covv(mi, ri) = T_roi.CoV_pct(find(idx, 1, 'first'));
            end
        end
    end

    ax = nexttile(tl, ci);
    bar(ax, x, mape);
    grid(ax, 'on');
    ylabel(ax, 'MAPE (%)');
    title(ax, sprintf('%s | accuracy', cases(ci)), ...
        'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 10);
    legend(ax, {'Background', 'Inclusion'}, 'Location', 'northwest');
    xlabel(ax, 'REQ M');
    set(ax, 'FontSize', 9);

    ax = nexttile(tl, numel(cases) + ci);
    bar(ax, x, covv);
    grid(ax, 'on');
    ylabel(ax, 'CoV (%)');
    title(ax, sprintf('%s | precision', cases(ci)), ...
        'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 10);
    legend(ax, {'Background', 'Inclusion'}, 'Location', 'northwest');
    xlabel(ax, 'REQ M');
    set(ax, 'FontSize', 9);
end

exportgraphics(fig, fullfile(fig_dir, ...
    'level12_level06_kwave_roi_mape_cov.png'), ...
    'Resolution', 260, 'BackgroundColor', 'white');
close(fig);

end

function [A, x_cm, z_cm] = map_from_table(T, value_var)

nz = max(T.map_iz);
nx = max(T.map_ix);
A = nan(nz, nx);
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
A(idx) = T.(value_var);
x_cm = unique(T.x_center_m, 'stable') * 100;
z_cm = unique(T.z_center_m, 'stable') * 100;

end

function draw_geometry(ax, KW)

hold(ax, 'on');
th = linspace(0, 2*pi, 300);
xc = KW.inclusion_center_m(1) * 100;
zc = KW.inclusion_center_m(2) * 100;
r = KW.inclusion_radius_m * 100;
plot(ax, xc + r*cos(th), zc + r*sin(th), 'w--', 'LineWidth', 1.2);
for ri = 1:numel(KW.ROI)
    roi = KW.ROI(ri);
    x0 = roi.xlim_m(1) * 100;
    z0 = roi.zlim_m(1) * 100;
    w = diff(roi.xlim_m) * 100;
    h = diff(roi.zlim_m) * 100;
    rectangle(ax, 'Position', [x0 z0 w h], ...
        'EdgeColor', 'w', 'LineWidth', 1.1);
end
hold(ax, 'off');

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

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

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

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    if ismember(vars_all(i), vars)
        continue;
    end
    name_i = char(vars_all(i));
    string_like = any(endsWith(vars_all(i), ...
        ["name", "label", "type", "role", "set", "model"]));
    if string_like
        T.(name_i) = strings(height(T), 1);
    else
        T.(name_i) = nan(height(T), 1);
    end
end

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end
