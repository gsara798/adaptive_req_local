%% run_test_10_heterogeneous_maps.m
% Test 10: simple heterogeneous phantoms with current residual-corrected model.
%
% This is a qualitative/diagnostic test, separate from Test 09 robustness.
% It simulates simple bilayer and inclusion phantoms, predicts patch-level SWS,
% summarizes ROIs, and plots maps/statistics.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
experiments_dir = fileparts(this_file);
root_dir = fileparts(experiments_dir);

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Output folder

out_root = fullfile(root_dir, 'outputs', 'test_10_heterogeneous_maps');
run_name = "test_10_heterogeneous_maps_" + ...
    string(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
output_dir = fullfile(out_root, char(run_name));
figure_dir = fullfile(output_dir, 'figures');
table_dir = fullfile(output_dir, 'tables');
data_dir = fullfile(output_dir, 'data');

make_dir(output_dir);
make_dir(figure_dir);
make_dir(table_dir);
make_dir(data_dir);

%% Train current Level 16 residual-corrected model on Test 08

fprintf('\nTraining current residual-corrected model from Test 08...\n');
[MODEL_base, MODEL_residual, predictors] = train_current_residual_model(root_dir);

%% Heterogeneous test settings

cfg_base = adaptive_req.config.default_sim_config( ...
    'WaveModel', 'spherical', ...
    'f0', 500, ...
    'cs_bg', 3.0, ...
    'Nwaves', 2000, ...
    'SNR', Inf, ...
    'AmpJitter', 0, ...
    'Lx', 0.06, ...
    'Lz', 0.06, ...
    'dx', 2.5e-4, ...
    'dz', 2.5e-4, ...
    'SourceSampling', 'cone', ...
    'AngularSamplingMethod', 'fibonacci', ...
    'ConeAxis', [1 0 0], ...
    'ConeHalfAngleDeg', 70, ...
    'ForceInPlaneWave', true, ...
    'UseParfor', true);

feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', 3, ...
    'cs_guess', 3.0, ...
    'gamma_win', 1.0, ...
    'pad_factor', 1.0);

phantoms = build_phantom_specs();

T_all = table();
phantom_outputs = struct([]);

for pidx = 1:numel(phantoms)
    spec = phantoms(pidx);
    fprintf('\nRunning heterogeneous phantom: %s\n', spec.name);

    cfg = cfg_base;
    cfg.Seed = 10000 + pidx;
    cfg.MaskType = spec.mask_type;
    cfg.cs_inc = spec.cs_inc;
    cfg.MaskParams = spec.mask_params;

    [T_phantom, sim, patch_pack] = predict_heterogeneous_phantom( ...
        cfg, feat_cfg, MODEL_base, MODEL_residual);

    T_phantom.phantom_name = repmat(string(spec.name), height(T_phantom), 1);
    T_phantom.phantom_type = repmat(string(spec.mask_type), height(T_phantom), 1);

    T_all = [T_all; T_phantom]; %#ok<AGROW>

    phantom_outputs(pidx).name = spec.name;
    phantom_outputs(pidx).sim = sim;
    phantom_outputs(pidx).patch_pack = patch_pack;
    phantom_outputs(pidx).T = T_phantom;

    plot_heterogeneous_maps(sim, patch_pack, T_phantom, spec, figure_dir);
end

%% ROI statistics

T_roi = summarize_roi_metrics(T_all);
T_overall = summarize_roi_metrics(add_all_roi(T_all));

writetable(T_all, fullfile(table_dir, 'test10_patch_predictions.csv'));
writetable(T_roi, fullfile(table_dir, 'test10_roi_metrics.csv'));
writetable(T_overall, fullfile(table_dir, 'test10_overall_metrics.csv'));

plot_roi_statistics(T_roi, figure_dir);
plot_patch_error_distributions(T_all, figure_dir);

save(fullfile(data_dir, 'test_10_heterogeneous_maps.mat'), ...
    'T_all', 'T_roi', 'T_overall', 'phantom_outputs', ...
    'MODEL_base', 'MODEL_residual', 'predictors', '-v7.3');

fprintf('\n============================================================\n');
fprintf('TEST 10 HETEROGENEOUS ROI METRICS\n');
fprintf('============================================================\n');
disp(T_roi(:, {'phantom_name', 'roi_name', 'n', 'MAPE_pct', ...
    'RMSE_pct', 'bias_pct', 'CoV_pct'}));
fprintf('\nTest 10 complete. Results saved to:\n%s\n', output_dir);

%% Local helpers

function [MODEL_base, MODEL_residual, predictors] = train_current_residual_model(root_dir)

[T_mc, ~, ~] = adaptive_req.analysis.load_mc_results( ...
    'test_08_advanced_angular_features', ...
    'RootDir', root_dir, ...
    'Verbose', false);

COL = adaptive_req.analysis.detect_mc_columns(T_mc);
T_feat = add_ecum_features_from_mapping(T_mc);
T_feat = adaptive_req.analysis.add_effective_window_metrics( ...
    T_feat, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', 3.0);

ecum_features = string(fieldnames( ...
    adaptive_req.quantile.extract_ecum_shape_features( ...
        T_feat.req_mapping{find_first_mapping(T_feat)})));

base_predictors = [
    "ang_entropy"
    "radial_entropy"
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"
    "width_75_25_rel"
    "width_90_50_rel"
    "width_90_10_rel"
    "lowk_frac_rel"
    "midband_frac_rel"
    "highk_frac_rel"
    "circ_var"
    "dom_dir_frac"
    "window_max_frac"
    "window_cf"
];

predictors = keep_existing_predictors(T_feat, [base_predictors; ecum_features]);

[MODEL_base, T_base_pred, ~] = ...
    adaptive_req.analysis.train_q_model_from_predictors( ...
        T_feat, ...
        predictors, ...
        'QVar', COL.q, ...
        'ModelName', "current_base_ecum_srad_proxy", ...
        'SplitMode', 'condition', ...
        'ConditionVar', COL.condition, ...
        'TrainFraction', 0.70, ...
        'RandomSeed', 1701, ...
        'ModelTypes', "bagged_trees", ...
        'NumLearningCycles', 300, ...
        'MinLeafSize', 8, ...
        'ClipPredictions', true, ...
        'Verbose', false);

T_residual_data = add_base_predictions_to_feature_table(T_feat, T_base_pred);
T_residual_data.q_residual_target = ...
    T_residual_data.q_base_true - T_residual_data.q_base_pred;
T_residual_data.q_base_pred_feature = T_residual_data.q_base_pred;
train_mask = T_residual_data.split == "train";

residual_predictors = keep_existing_predictors(T_residual_data, ...
    [predictors; "q_base_pred_feature"]);

[MODEL_residual, ~, ~] = ...
    adaptive_req.analysis.train_q_model_from_predictors( ...
        T_residual_data, ...
        residual_predictors, ...
        'QVar', 'q_residual_target', ...
        'ModelName', "current_residual_corrector", ...
        'SplitMode', 'condition', ...
        'ConditionVar', COL.condition, ...
        'TrainFraction', 0.70, ...
        'TrainMask', train_mask, ...
        'RandomSeed', 2701, ...
        'ModelTypes', "bagged_trees", ...
        'NumLearningCycles', 220, ...
        'MinLeafSize', 12, ...
        'ClipPredictions', false, ...
        'Verbose', false);

end

function specs = build_phantom_specs()

specs = struct([]);

specs(1).name = "bilayer_slow_fast";
specs(1).mask_type = 'bilayer';
specs(1).cs_inc = 4.0;
specs(1).mask_params = struct( ...
    'Bi_Angle', pi/2, ...
    'Bi_Offset', 0.030, ...
    'SigmaEdge', 0.0015);

specs(2).name = "circular_fast_inclusion";
specs(2).mask_type = 'circle';
specs(2).cs_inc = 4.0;
specs(2).mask_params = struct( ...
    'Center', [0.030, 0.030], ...
    'Radius', 0.012, ...
    'SigmaEdge', 0.0015);

end

function [T_pred, sim, patch_pack] = predict_heterogeneous_phantom( ...
    cfg, feat_cfg, MODEL_base, MODEL_residual)

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, ...
    'Nbins', 'auto', ...
    'Nbins_auto_oversample', 1, ...
    'Nbins_min', 16, ...
    'smooth_sigma', 1);

sim = adaptive_req.simulate.run_single_simulation(cfg);
patch_pack = adaptive_req.simulate.build_patch_windows( ...
    cfg, feat_cfg, ...
    'Pattern', 'grid', ...
    'GridSize', [9 9], ...
    'CoverageFraction', 0.78);

T_feat = extract_patch_feature_table(sim, cfg, feat_cfg, req_cfg, patch_pack);

T_base = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_base, T_feat, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "test10_base");

T_resid_input = T_feat;
T_resid_input.q_base_pred_feature = T_base.q_pred;

T_resid = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_residual, T_resid_input, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "test10_residual");

T_pred = T_feat(:, {'patch_idx', 'patch_label', 'cx', 'cz', ...
    'x_center_m', 'z_center_m', 'roi_name', ...
    'cs_true_patch_mean', 'cs_true_patch_median', ...
    'cs_true_patch_std', 'SIM_f0', 'SIM_cs_bg', 'REQ_M', ...
    'Omega_sr', 'SIM_WaveModel'});
T_pred.q_base_pred = T_base.q_pred;
T_pred.q_residual_pred = T_resid.q_pred;
T_pred.q_pred = min(max(T_pred.q_base_pred + T_pred.q_residual_pred, 0), 1);

T_pred.cs_pred_from_q = q_to_cs_for_table(T_pred.q_pred, T_feat.req_mapping, cfg.f0);
T_pred.cs_true = T_pred.cs_true_patch_median;
T_pred.cs_error = T_pred.cs_pred_from_q - T_pred.cs_true;
T_pred.cs_abs_error = abs(T_pred.cs_error);
T_pred.cs_error_pct = 100 * T_pred.cs_error ./ T_pred.cs_true;
T_pred.cs_abs_error_pct = abs(T_pred.cs_error_pct);

end

function T = extract_patch_feature_table(sim, cfg, feat_cfg, req_cfg, patch_pack)

n = patch_pack.n_patches;
rows(n) = struct();

for pidx = 1:n
    xi = patch_pack.x_idx_list{pidx};
    zi = patch_pack.z_idx_list{pidx};
    V_patch = sim.Uxz(zi, xi);
    cs_patch = sim.cs_map(zi, xi);

    feat = adaptive_req.features.req_extract_patch_features( ...
        V_patch, cfg.dx, cfg.dz, cfg.f0, cfg.cs_bg, feat_cfg);
    [~, req_curve] = adaptive_req.quantile.compute_quantile_from_patch( ...
        V_patch, cfg, req_cfg);
    req_shape = adaptive_req.quantile.extract_ecum_shape_features(req_curve);

    rows(pidx).patch_idx = pidx;
    rows(pidx).patch_label = patch_pack.patch_labels(pidx);
    rows(pidx).cx = patch_pack.cx_list(pidx);
    rows(pidx).cz = patch_pack.cz_list(pidx);
    rows(pidx).x_center_m = sim.x(patch_pack.cx_list(pidx));
    rows(pidx).z_center_m = sim.z(patch_pack.cz_list(pidx));
    rows(pidx).cs_true_patch_mean = mean(cs_patch(:), 'omitnan');
    rows(pidx).cs_true_patch_median = median(cs_patch(:), 'omitnan');
    rows(pidx).cs_true_patch_std = std(cs_patch(:), 'omitnan');
    rows(pidx).roi_name = classify_roi(rows(pidx).cs_true_patch_median, cfg);

    rows(pidx).SIM_WaveModel = string(cfg.WaveModel);
    rows(pidx).SIM_f0 = cfg.f0;
    rows(pidx).SIM_cs_bg = cfg.cs_bg;
    rows(pidx).REQ_M = feat_cfg.M;
    rows(pidx).Omega_sr = cone_omega_sr(cfg.ConeHalfAngleDeg);
    rows(pidx).REQ_Nbins_effective = req_curve.Nbins_effective;
    rows(pidx).req_mapping = {adaptive_req.quantile.make_req_mapping(req_curve)};

    feature_names = fieldnames(feat.scalar);
    for j = 1:numel(feature_names)
        rows(pidx).(feature_names{j}) = feat.scalar.(feature_names{j});
    end

    req_names = fieldnames(req_shape);
    for j = 1:numel(req_names)
        rows(pidx).(req_names{j}) = req_shape.(req_names{j});
    end
end

T = struct2table(rows);

end

function roi = classify_roi(cs_value, cfg)

if cs_value > cfg.cs_bg + 0.35
    roi = "fast";
elseif cs_value < cfg.cs_bg - 0.35
    roi = "slow";
else
    roi = "background_or_boundary";
end

end

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)

n = numel(q_pred);
cs_pred = NaN(n, 1);
for i = 1:n
    cs_pred(i) = adaptive_req.quantile.quantile_to_cs( ...
        mappings{i}, q_pred(i), f0);
end

end

function T_roi = summarize_roi_metrics(T)

[G, T_roi] = findgroups(T(:, {'phantom_name', 'roi_name'}));
T_roi.n = splitapply(@numel, T.cs_error_pct, G);
T_roi.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), T.cs_error_pct, G);
T_roi.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.cs_error_pct, G);
T_roi.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.cs_error_pct, G);
T_roi.CoV_pct = splitapply(@cov_pct, T.cs_pred_from_q, G);
T_roi.cs_true_median = splitapply(@(x) median(x, 'omitnan'), T.cs_true, G);
T_roi.cs_pred_median = splitapply(@(x) median(x, 'omitnan'), T.cs_pred_from_q, G);

end

function T = add_all_roi(T)

T.roi_name = repmat("all", height(T), 1);

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || mean(x) == 0
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function plot_heterogeneous_maps(sim, patch_pack, T, spec, figure_dir)

fig = figure('Color', 'w', 'Position', [100 100 1450 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

x_cm = sim.x * 100;
z_cm = sim.z * 100;

ax = nexttile(tl);
imagesc(ax, x_cm, z_cm, sim.cs_map);
axis(ax, 'image');
set(ax, 'YDir', 'normal');
colorbar(ax);
title(ax, 'True c_s map', 'Interpreter', 'tex');
xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');
hold(ax, 'on');
plot_patch_centers(ax, sim, patch_pack);

ax = nexttile(tl);
scatter(ax, 100*T.x_center_m, 100*T.z_center_m, 55, ...
    T.cs_pred_from_q, 'filled');
axis(ax, 'image');
set(ax, 'YDir', 'normal');
xlim(ax, [min(x_cm), max(x_cm)]);
ylim(ax, [min(z_cm), max(z_cm)]);
colorbar(ax);
title(ax, 'Predicted c_s at patch centers', 'Interpreter', 'tex');
xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');

ax = nexttile(tl);
scatter(ax, 100*T.x_center_m, 100*T.z_center_m, 55, ...
    T.cs_error_pct, 'filled');
axis(ax, 'image');
set(ax, 'YDir', 'normal');
xlim(ax, [min(x_cm), max(x_cm)]);
ylim(ax, [min(z_cm), max(z_cm)]);
colorbar(ax);
title(ax, 'Signed SWS error (%)');
xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');

ax = nexttile(tl);
boxchart(ax, categorical(T.roi_name), T.cs_abs_error_pct);
ylabel(ax, 'Absolute SWS error (%)');
title(ax, 'Patch error by ROI');
grid(ax, 'on');

title(tl, "Test 10 " + string(spec.name), ...
    'Interpreter', 'none', 'FontWeight', 'bold');
save_fig(fig, figure_dir, "test10_maps_" + string(spec.name));

end

function plot_patch_centers(ax, sim, patch_pack)

x = sim.x(patch_pack.cx_list) * 100;
z = sim.z(patch_pack.cz_list) * 100;
scatter(ax, x, z, 16, 'w', 'filled', 'MarkerEdgeColor', 'k');

end

function plot_roi_statistics(T_roi, figure_dir)

fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
bar(categorical(T_roi.phantom_name + " / " + T_roi.roi_name), ...
    T_roi.MAPE_pct);
ylabel(ax, 'MAPE (%)');
title(ax, 'ROI accuracy');
grid(ax, 'on');
xtickangle(ax, 35);

ax = nexttile(tl);
bar(categorical(T_roi.phantom_name + " / " + T_roi.roi_name), ...
    T_roi.CoV_pct);
ylabel(ax, 'CoV (%)');
title(ax, 'ROI precision');
grid(ax, 'on');
xtickangle(ax, 35);

save_fig(fig, figure_dir, 'test10_roi_statistics');

end

function plot_patch_error_distributions(T, figure_dir)

fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
boxchart(categorical(T.phantom_name + " / " + T.roi_name), ...
    T.cs_abs_error_pct);
ylabel('Absolute SWS error (%)');
title('Test 10 patch-level error distributions');
grid on;
xtickangle(35);
save_fig(fig, figure_dir, 'test10_patch_error_distributions');

end

function omega = cone_omega_sr(half_angle_deg)

theta = deg2rad(half_angle_deg);
omega = 2*pi*(1 - cos(theta));

end

function idx = find_first_mapping(T)

idx = find(~cellfun(@isempty, T.req_mapping), 1, 'first');
if isempty(idx)
    error('No non-empty req_mapping entries were found.');
end

end

function T = add_ecum_features_from_mapping(T)

n = height(T);
feat0 = adaptive_req.quantile.extract_ecum_shape_features( ...
    T.req_mapping{find_first_mapping(T)});
names = fieldnames(feat0);

for j = 1:numel(names)
    T.(names{j}) = NaN(n, 1);
end

for i = 1:n
    mapping_i = T.req_mapping{i};
    if isempty(mapping_i)
        continue;
    end
    feat_i = adaptive_req.quantile.extract_ecum_shape_features(mapping_i);
    for j = 1:numel(names)
        T.(names{j})(i) = feat_i.(names{j});
    end
end

end

function predictors_out = keep_existing_predictors(T, predictors_in)

vars = string(T.Properties.VariableNames);
predictors_in = string(predictors_in);
predictors_out = strings(0, 1);
for i = 1:numel(predictors_in)
    p_i = predictors_in(i);
    candidates = [p_i; p_i + "_mean"];
    if any(ismember(vars, candidates))
        predictors_out(end + 1, 1) = p_i; %#ok<AGROW>
    end
end

end

function T_out = add_base_predictions_to_feature_table(T_feat, T_base_pred)

key_vars = {'condition_id', 'step_idx', 'realization_idx', 'patch_idx'};
T_base = T_base_pred(:, [key_vars, {'split', 'q_true', 'q_pred'}]);
T_base.Properties.VariableNames{'q_true'} = 'q_base_true';
T_base.Properties.VariableNames{'q_pred'} = 'q_base_pred';
[~, ia] = unique(T_base(:, key_vars), 'rows', 'stable');
T_base = T_base(ia, :);
T_out = innerjoin(T_feat, T_base, 'Keys', key_vars);

end

function make_dir(path_in)

if ~exist(path_in, 'dir')
    mkdir(path_in);
end

end

function save_fig(fig, output_dir, base_name)

exportgraphics(fig, fullfile(output_dir, base_name + ".png"), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, base_name + ".pdf"), ...
    'ContentType', 'vector');
close(fig);

end
