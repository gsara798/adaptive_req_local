%% analyze_test_09_level01_stress_robustness.m
% Analyze Test 09 stress-screening dataset with current residual model.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load Test 09 and Test 08

[T09, ~, PATHS09] = adaptive_req.analysis.load_mc_results( ...
    'test_09_stress_robustness', ...
    'RootDir', root_dir, ...
    'Verbose', true);

analysis_dir = fullfile(PATHS09.analysis_dir, ...
    'level_01_stress_robustness_model_eval');
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

fprintf('\nTraining current residual model from Test 08...\n');
[MODEL_base, MODEL_residual] = train_current_residual_model(root_dir);

%% Add Ecum/Srad-proxy features to Test 09

T09_feat = add_ecum_features_from_mapping(T09);
T09_feat = adaptive_req.analysis.add_effective_window_metrics( ...
    T09_feat, ...
    'MVar', 'REQ_M', ...
    'F0Var', 'SIM_f0', ...
    'CsTrueVar', 'SIM_cs_bg', ...
    'CsGuess', 3.0);

%% Predict q and SWS

T_base = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_base, T09_feat, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "stress_base");

T_resid_input = T09_feat;
T_resid_input.q_base_pred_feature = T_base.q_pred;

T_resid = adaptive_req.analysis.predict_q_model_from_table( ...
    MODEL_residual, T_resid_input, ...
    'ModelType', 'bagged_trees', ...
    'ModelName', "stress_residual");

T_pred = T_base;
T_pred.model_name = repmat("stress_residual_corrected", height(T_pred), 1);
T_pred.model_role = repmat("stress_eval", height(T_pred), 1);
T_pred.q_base_pred = T_base.q_pred;
T_pred.q_residual_pred = T_resid.q_pred;
T_pred.q_pred = min(max(T_pred.q_base_pred + T_pred.q_residual_pred, 0), 1);
T_pred.q_true = T09_feat.q_theory;
T_pred.q_error = T_pred.q_pred - T_pred.q_true;
T_pred.q_abs_error = abs(T_pred.q_error);

T_pred.cs_pred_from_q = q_to_cs_for_table( ...
    T_pred.q_pred, T09_feat.req_mapping, T09_feat.SIM_f0);
T_pred.cs_true = T09_feat.SIM_cs_bg;
T_pred.cs_error = T_pred.cs_pred_from_q - T_pred.cs_true;
T_pred.cs_abs_error = abs(T_pred.cs_error);
T_pred.cs_error_pct = 100 * T_pred.cs_error ./ T_pred.cs_true;
T_pred.cs_abs_error_pct = abs(T_pred.cs_error_pct);

copy_vars = ["SIM_Nwaves", "SIM_SNR", "SIM_ForceInPlaneWave", ...
    "SIM_WaveModel", "REQ_M", "Omega_sr", "step_idx"];
for i = 1:numel(copy_vars)
    name = copy_vars(i);
    if ismember(name, string(T09_feat.Properties.VariableNames))
        T_pred.(char(name)) = T09_feat.(char(name));
    end
end

%% Summaries

T_overall = summarize_metrics(T_pred, "model_name");
T_by_snr = summarize_metrics(T_pred, ["SIM_SNR", "REQ_M"]);
T_by_nwaves = summarize_metrics(T_pred, ["SIM_Nwaves", "REQ_M"]);
T_by_force = summarize_metrics(T_pred, ["SIM_ForceInPlaneWave", "REQ_M"]);
T_by_wave = summarize_metrics(T_pred, ["SIM_WaveModel", "REQ_M"]);
T_by_aperture = summarize_metrics(T_pred, ["Omega_sr", "REQ_M"]);

writetable(T_pred, fullfile(analysis_dir, 'test09_stress_predictions.csv'));
writetable(T_overall, fullfile(analysis_dir, 'test09_overall_metrics.csv'));
writetable(T_by_snr, fullfile(analysis_dir, 'test09_metrics_by_SNR_M.csv'));
writetable(T_by_nwaves, fullfile(analysis_dir, 'test09_metrics_by_Nwaves_M.csv'));
writetable(T_by_force, fullfile(analysis_dir, 'test09_metrics_by_force_in_plane_M.csv'));
writetable(T_by_wave, fullfile(analysis_dir, 'test09_metrics_by_wave_model_M.csv'));
writetable(T_by_aperture, fullfile(analysis_dir, 'test09_metrics_by_aperture_M.csv'));

plot_stress_summary(T_by_snr, T_by_nwaves, T_by_force, T_by_wave, analysis_dir);
plot_stress_aperture(T_by_aperture, analysis_dir);
plot_stress_error_distributions(T_pred, analysis_dir);

save(fullfile(analysis_dir, 'test09_level01_stress_robustness.mat'), ...
    'T_pred', 'T_overall', 'T_by_snr', 'T_by_nwaves', ...
    'T_by_force', 'T_by_wave', 'T_by_aperture', '-v7.3');

fprintf('\n============================================================\n');
fprintf('TEST 09 STRESS OVERALL METRICS\n');
fprintf('============================================================\n');
disp(T_overall);
fprintf('\nStress analysis complete. Results saved to:\n%s\n', analysis_dir);

%% Local helpers

function [MODEL_base, MODEL_residual] = train_current_residual_model(root_dir)

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
        T_feat, predictors, ...
        'QVar', COL.q, ...
        'ModelName', "stress_train_base", ...
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
        T_residual_data, residual_predictors, ...
        'QVar', 'q_residual_target', ...
        'ModelName', "stress_train_residual", ...
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

function T = add_ecum_features_from_mapping(T)

n = height(T);
feat0 = adaptive_req.quantile.extract_ecum_shape_features( ...
    T.req_mapping{find_first_mapping(T)});
names = fieldnames(feat0);
for j = 1:numel(names)
    T.(names{j}) = NaN(n, 1);
end
for i = 1:n
    feat_i = adaptive_req.quantile.extract_ecum_shape_features(T.req_mapping{i});
    for j = 1:numel(names)
        T.(names{j})(i) = feat_i.(names{j});
    end
end

end

function idx = find_first_mapping(T)
idx = find(~cellfun(@isempty, T.req_mapping), 1, 'first');
if isempty(idx), error('No req_mapping found.'); end
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

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)
n = numel(q_pred);
cs_pred = NaN(n, 1);
for i = 1:n
    cs_pred(i) = adaptive_req.quantile.quantile_to_cs( ...
        mappings{i}, q_pred(i), f0(i));
end
end

function T_metrics = summarize_metrics(T, group_vars)
[G, T_metrics] = findgroups(T(:, cellstr(group_vars)));
T_metrics.n = splitapply(@numel, T.cs_error_pct, G);
T_metrics.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), T.cs_error_pct, G);
T_metrics.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.cs_error_pct, G);
T_metrics.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.cs_error_pct, G);
T_metrics.p95_APE_pct = splitapply(@(x) prctile(abs(x), 95), T.cs_error_pct, G);
T_metrics.high_error_pct = splitapply(@(x) 100*mean(abs(x) > 20), T.cs_error_pct, G);
end

function plot_stress_summary(T_snr, T_nwaves, T_force, T_wave, output_dir)
fig = figure('Color', 'w', 'Position', [100 100 1450 800]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_grouped(nexttile(tl), T_snr, "SIM_SNR", "REQ_M", "RMSE_pct", 'SNR (dB)');
plot_grouped(nexttile(tl), T_nwaves, "SIM_Nwaves", "REQ_M", "RMSE_pct", 'N waves');
plot_grouped(nexttile(tl), T_force, "SIM_ForceInPlaneWave", "REQ_M", "RMSE_pct", 'Force in-plane');
plot_grouped(nexttile(tl), T_wave, "SIM_WaveModel", "REQ_M", "RMSE_pct", 'Wave model');
title(tl, 'Test 09 stress robustness: SWS RMSE', 'FontWeight', 'bold');
save_fig(fig, output_dir, 'test09_stress_summary');
end

function plot_stress_aperture(T_ap, output_dir)
fig = figure('Color', 'w', 'Position', [100 100 900 560]);
hold on;
Mvals = unique(T_ap.REQ_M, 'sorted');
for i = 1:numel(Mvals)
    Ti = sortrows(T_ap(T_ap.REQ_M == Mvals(i), :), 'Omega_sr');
    plot(Ti.Omega_sr, Ti.RMSE_pct, '-o', 'DisplayName', "M=" + string(Mvals(i)));
end
xlabel('\Omega (sr)', 'Interpreter', 'tex');
ylabel('SWS RMSE (%)');
title('Test 09 stress error vs aperture');
legend('Location', 'best');
grid on;
save_fig(fig, output_dir, 'test09_stress_aperture');
end

function plot_stress_error_distributions(T, output_dir)
fig = figure('Color', 'w', 'Position', [100 100 1200 520]);
boxchart(categorical("SNR " + string(T.SIM_SNR) + " / M " + string(T.REQ_M)), ...
    T.cs_abs_error_pct);
ylabel('Absolute SWS error (%)');
title('Test 09 patch errors by SNR and M');
grid on;
xtickangle(45);
save_fig(fig, output_dir, 'test09_stress_error_distributions');
end

function plot_grouped(ax, T, xvar, groupvar, yvar, xlabel_text)
xvals = unique(T.(char(xvar)), 'stable');
gvals = unique(T.(char(groupvar)), 'sorted');
Y = NaN(numel(xvals), numel(gvals));
for i = 1:numel(xvals)
    for j = 1:numel(gvals)
        mask = T.(char(xvar)) == xvals(i) & T.(char(groupvar)) == gvals(j);
        if any(mask), Y(i,j) = T.(char(yvar))(find(mask, 1)); end
    end
end
bar(ax, categorical(string(xvals)), Y);
xlabel(ax, xlabel_text);
ylabel(ax, 'SWS RMSE (%)');
legend(ax, "M=" + string(gvals), 'Location', 'best');
grid(ax, 'on');
end

function save_fig(fig, output_dir, base_name)
exportgraphics(fig, fullfile(output_dir, base_name + ".png"), 'Resolution', 300);
exportgraphics(fig, fullfile(output_dir, base_name + ".pdf"), 'ContentType', 'vector');
close(fig);
end
