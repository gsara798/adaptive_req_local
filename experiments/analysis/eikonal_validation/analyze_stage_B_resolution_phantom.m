%% analyze_stage_B_resolution_phantom.m
% Analysis for Stage B resolution phantom / geometry transfer validation.
%
% This script reads Stage B runner outputs only. It does not rerun Eikonal
% simulations, retrain models, apply corrections, or use oracle quantities as
% predictors. Truth, ROI labels, diameter, purity, and interface distance are
% evaluation-only diagnostics.

clear; clc; close all;
format compact;

%% Setup
this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(fileparts(this_file))));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.0);

CFG = load_stage_config(root_dir);
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
OUT = make_output_dirs(root_dir, CFG);

fprintf('\nStage B analysis: resolution phantom / geometry transfer\n');
fprintf('Mode label: %s\nOutput root: %s\n', CFG.Mode, OUT.root_dir);

patchCsv = fullfile(OUT.table_dir, 'stage_B_patch_level_results.csv');
assert(exist(patchCsv, 'file') == 2, ['Missing patch-level table. Run the Stage B runner first: ' patchCsv]);
T = readtable(patchCsv, 'TextType','string');
T = normalize_loaded_table(T);
assert(height(T) > 0, 'Patch-level table is empty.');
assert(all(isfinite(T.SWS_pred)) && all(isfinite(T.q_pred)), 'Non-finite predictions found.');

% Older Stage B runner outputs labeled the broad 4-8 mm distance band after
% inclusion core, which hid D=10/14 mm core points. Rebuild analysis ROI labels
% from saved distance/purity diagnostics so figures reflect the configured core
% definition without rerunning expensive Eikonal/REQ conditions.
T = repair_roi_labels_for_analysis(T, CFG);
T = add_inclusion_center_roi(T);
T = add_analysis_bins(T);

%% Summaries
S_overall = summarize_metrics(T, "model_name");
S_model = S_overall;
S_phantom = summarize_metrics(T, ["model_name","phantom_name"]);
S_diam = summarize_metrics(T, ["model_name","inclusion_sws","inclusion_diameter_mm"]);
S_incSws = summarize_metrics(T, ["model_name","inclusion_sws"]);
S_freq = summarize_metrics(T, ["model_name","f0"]);
S_field = summarize_metrics(T, ["model_name","field_regime"]);
S_Dlambda = summarize_metrics(T, ["model_name","D_over_lambda_bin"]);
S_DLwin = summarize_metrics(T, ["model_name","D_over_Lwin_bin"]);
S_Meff = summarize_metrics(T, ["model_name","M_eff_bin"]);
S_dist = summarize_metrics(T, ["model_name","distance_over_window_bin"]);
S_roi = summarize_metrics(T, ["model_name","roi_label"]);
S_roi_freq = summarize_metrics(T, ["model_name","roi_label","f0"]);
S_diam_freq = summarize_metrics(T, ["model_name","inclusion_sws","inclusion_diameter_mm","f0"]);
S_diam_field = summarize_metrics(T, ["model_name","inclusion_sws","inclusion_diameter_mm","field_regime"]);
S_condition = summarize_metrics(T, ["model_name","condition_id","phantom_name","inclusion_sws","f0","field_regime"]);
S_worst = sortrows(S_condition, 'MAPE_pct', 'descend');
S_worst = S_worst(1:min(40,height(S_worst)),:);
S_roi_diam_freq = summarize_metrics(T, ["model_name","roi_group","inclusion_sws","inclusion_diameter_mm","f0"]);
S_roi_diam_field = summarize_metrics(T, ["model_name","roi_group","inclusion_sws","inclusion_diameter_mm","field_regime"]);
S_roi_physics = summarize_metrics(T, ["model_name","roi_group","D_over_lambda_bin","D_over_Lwin_bin"]);
S_roi_resolution_axes = summarize_metrics(T, ["model_name","roi_group","roi_label","inclusion_sws","inclusion_diameter_mm","f0","field_regime"]);
writetable(S_roi_resolution_axes, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_resolution_axes.csv'));
S_window_lambda = summarize_metrics(T, ["model_name","roi_group","roi_label","true_sws_level","field_regime"]);
writetable(S_window_lambda, fullfile(OUT.table_dir, 'stage_B_summary_by_window_lambda_roi.csv'));
S_condition_avg_diam = summarize_condition_averaged(T, ["model_name","inclusion_sws","inclusion_diameter_mm"]);
S_condition_avg_roi = summarize_condition_averaged(T, ["model_name","roi_group"]);

write_summary_tables(OUT, S_overall, S_model, S_phantom, S_diam, S_incSws, S_freq, S_field, ...
    S_Dlambda, S_DLwin, S_Meff, S_dist, S_roi, S_roi_freq, S_diam_freq, S_diam_field, ...
    S_worst, S_roi_diam_freq, S_roi_diam_field, S_roi_physics, S_condition_avg_diam, S_condition_avg_roi);

%% Figures
make_analysis_figures(T, OUT);
make_extra_diagnostics(T, OUT);

%% README refresh
write_stage_readme(OUT, CFG, T, S_overall, S_diam, S_roi, S_worst);
write_docs_readme(root_dir, CFG);

%% Console summary
fprintf('\nStage B analysis complete.\n');
fprintf('Patch/model rows: %d\n', height(T));
print_summary_block(S_overall, 'Global summary');
print_summary_block(S_diam, 'By inclusion diameter');
print_summary_block(S_roi, 'By ROI');
X = S_worst(string(S_worst.model_name)=="q_spectrum_plus_composition",:);
if ~isempty(X)
    fprintf('Worst q_spectrum_plus_composition condition: %s | MAPE %.2f%% | bias %.2f%%\n', ...
        X.condition_id(1), X.MAPE_pct(1), X.signed_bias_pct(1));
end
fprintf('Analysis figures: %s\nREADME: %s\n', OUT.analysis_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Local functions
function CFG = load_stage_config(root_dir)
configPath = fullfile(root_dir, 'configs', 'eikonal_validation', 'stage_B_resolution_phantom.json');
assert(exist(configPath, 'file') == 2, 'Missing Stage B config: %s', configPath);
CFG = jsondecode(fileread(configPath));
CFG.ConfigPath = string(configPath);
CFG.OutputRoot = string(CFG.OutputRoot);
CFG.FrequenciesHz = double(CFG.FrequenciesHz(:)).';
CFG.FieldRegimes = string(CFG.FieldRegimes(:)).';
CFG.MList = double(CFG.MList(:)).';
CFG.RealismLevels = string(CFG.RealismLevels(:)).';
CFG.REQ.Nbins = string(CFG.REQ.Nbins);
CFG.REQ.EdgeMode = string(CFG.REQ.EdgeMode);
end

function OUT = make_output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, char(CFG.OutputRoot));
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.analysis_dir = fullfile(OUT.figure_dir, 'analysis_summary');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
ensure_dir(OUT.root_dir); ensure_dir(OUT.table_dir); ensure_dir(OUT.figure_dir);
ensure_dir(OUT.analysis_dir); ensure_dir(OUT.data_dir);
end

function T = normalize_loaded_table(T)
stringVars = ["condition_id","phantom_name","geometry","field_regime","realism_level", ...
    "model_name","roi_label","source_layout_id"];
for v = stringVars
    if ismember(v, string(T.Properties.VariableNames))
        T.(v) = string(T.(v));
    end
end
logicalVars = ["high_error10","high_error20"];
for v = logicalVars
    if ismember(v, string(T.Properties.VariableNames)) && ~islogical(T.(v))
        T.(v) = T.(v) ~= 0;
    end
end
end

function T = repair_roi_labels_for_analysis(T, CFG)
required = ["distance_to_interface_m","SWS_true","background_sws","patch_purity","roi_label"];
if ~all(ismember(required, string(T.Properties.VariableNames)))
    warning('StageB:ROILabelRepairSkipped', 'Skipping ROI label repair because required columns are missing.');
    return;
end

T.roi_label_original = T.roi_label;
dmm = T.distance_to_interface_m * 1e3;
inside = T.SWS_true > T.background_sws + 0.1;
roi = repmat("background_other", height(T), 1);

farBg = ~inside & dmm >= CFG.ROI.BackgroundFarMinDistanceMm;
roi(farBg) = "background_far";
for bi = 1:size(CFG.ROI.InterfaceBandsMm,1)
    lo = CFG.ROI.InterfaceBandsMm(bi,1); hi = CFG.ROI.InterfaceBandsMm(bi,2);
    mask = dmm >= lo & dmm < hi;
    roi(mask) = sprintf('interface_%g_%gmm', lo, hi);
end
% Core labels intentionally override the broad 4-8 mm interface band.
core = inside & dmm >= CFG.ROI.CoreMinDistanceMm & T.patch_purity >= 0.95;
roi(core) = "inclusion_core";
T.roi_label = roi;
end

function T = add_inclusion_center_roi(T)
% Add a geometry-centered inclusion ROI for small inclusions. Unlike
% inclusion_core, this ROI does not require patch purity >= 0.95; it answers
% the practical question: what does the estimator report near the center of a
% small inclusion, even if the REQ window is partially mixed?
required = ["distance_to_interface_m","SWS_true","background_sws","inclusion_diameter_mm","roi_label"];
if ~all(ismember(required, string(T.Properties.VariableNames)))
    warning('StageB:CenterROISkipped', 'Skipping center ROI because required columns are missing.');
    return;
end
if ~ismember("roi_label_original", string(T.Properties.VariableNames))
    T.roi_label_original = T.roi_label;
end

dmm = T.distance_to_interface_m * 1e3;
radiusMm = T.inclusion_diameter_mm / 2;
inside = T.SWS_true > T.background_sws + 0.1;
% Central disk: inner half-radius of each inclusion. This is intentionally
% smaller than the inclusion but available even when no patch is pure enough
% to satisfy the conservative core definition.
centerMask = inside & dmm >= 0.5 * radiusMm;
T.roi_label(centerMask) = "inclusion_center_roi";
end

function T = add_analysis_bins(T)
T.lambda_over_Lwin = 1 ./ T.M_eff;
T.wavelengths_per_window = T.M_eff;
T.true_sws_level = strings(height(T),1);
T.true_sws_level(abs(T.SWS_true-2)<0.05) = "2 m/s background";
T.true_sws_level(abs(T.SWS_true-3)<0.05) = "3 m/s inclusion";
T.true_sws_level(abs(T.SWS_true-4)<0.05) = "4 m/s inclusion";
T.true_sws_level(T.true_sws_level=="") = "other";
T = add_bin(T, 'D_over_lambda', 'D_over_lambda_bin', [0 1 2 3 4 6 8 12 inf], ...
    ["0-1","1-2","2-3","3-4","4-6","6-8","8-12",">12"]);
T = add_bin(T, 'D_over_Lwin', 'D_over_Lwin_bin', [0 0.5 1 1.5 2 3 4 6 inf], ...
    ["0-0.5","0.5-1","1-1.5","1.5-2","2-3","3-4","4-6",">6"]);
T = add_bin(T, 'M_eff', 'M_eff_bin', [0 1.5 2 2.5 3 4 inf], ...
    ["0-1.5","1.5-2","2-2.5","2.5-3","3-4",">4"]);
T = add_bin(T, 'distance_to_interface_over_window_radius', 'distance_over_window_bin', [0 0.25 0.5 1 2 inf], ...
    ["0-0.25","0.25-0.5","0.5-1","1-2",">2"]);
T.roi_group = strings(height(T),1);
T.roi_group(startsWith(T.roi_label, "interface_")) = "interface";
T.roi_group(T.roi_label == "inclusion_core") = "inclusion_core";
T.roi_group(T.roi_label == "inclusion_center_roi") = "inclusion_center_roi";
T.roi_group(T.roi_label == "background_far") = "background_far";
T.roi_group(T.roi_group == "") = "other";
end

function T = add_bin(T, varName, binName, edges, labels)
x = T.(varName);
b = strings(numel(x),1);
for i = 1:numel(edges)-1
    lo = edges(i); hi = edges(i+1);
    if isinf(hi)
        mask = x >= lo;
    else
        mask = x >= lo & x < hi;
    end
    b(mask) = labels(i);
end
b(b=="") = "unknown";
T.(binName) = b;
end

function Tsum = summarize_metrics(T, groupVars)
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
        std(X.SWS_pred,'omitnan'), ...
        mean(X.SWS_true,'omitnan'), ...
        mean(X.q_pred,'omitnan'), ...
        std(X.q_pred,'omitnan'), ...
        mean(X.patch_purity,'omitnan'), ...
        mean(X.distance_to_interface_over_window_radius,'omitnan'), ...
        mean(X.D_over_lambda,'omitnan'), ...
        mean(X.D_over_Lwin,'omitnan'), ...
        mean(X.M_eff,'omitnan'), ...
        mean(X.lambda_over_Lwin,'omitnan'), ...
        mean(X.wavelengths_per_window,'omitnan'), ...
        'VariableNames', {'N_valid_patches','MAPE_pct','median_abs_pct_error', ...
        'MAE_SWS','RMSE_SWS','signed_bias_pct','median_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','mean_predicted_SWS','std_predicted_SWS', ...
        'mean_true_SWS','mean_q_pred','std_q_pred','mean_patch_purity', ...
        'mean_distance_over_window_radius','mean_D_over_lambda','mean_D_over_Lwin','mean_M_eff', ...
        'mean_lambda_over_Lwin','mean_wavelengths_per_window'});
end
Tsum = [groups vertcat(rows{:})];
end

function Tout = summarize_condition_averaged(T, groupVars)
if isstring(groupVars), groupVars = cellstr(groupVars); end
C = summarize_metrics(T, ["model_name","condition_id", string(setdiff(groupVars, {'model_name'}, 'stable'))]);
[G, groups] = findgroups(C(:, groupVars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = C(G == gi,:);
    rows{gi} = table(height(X), mean(X.MAPE_pct,'omitnan'), std(X.MAPE_pct,'omitnan'), ...
        mean(X.signed_bias_pct,'omitnan'), std(X.signed_bias_pct,'omitnan'), ...
        mean(X.high_error20_pct,'omitnan'), std(X.high_error20_pct,'omitnan'), ...
        'VariableNames', {'N_conditions','condition_mean_MAPE_pct','condition_std_MAPE_pct', ...
        'condition_mean_bias_pct','condition_std_bias_pct','condition_mean_high20_pct','condition_std_high20_pct'});
end
Tout = [groups vertcat(rows{:})];
end

function write_summary_tables(OUT, S_overall, S_model, S_phantom, S_diam, S_incSws, S_freq, S_field, ...
    S_Dlambda, S_DLwin, S_Meff, S_dist, S_roi, S_roi_freq, S_diam_freq, S_diam_field, ...
    S_worst, S_roi_diam_freq, S_roi_diam_field, S_roi_physics, S_condition_avg_diam, S_condition_avg_roi)
writetable(S_overall, fullfile(OUT.table_dir, 'stage_B_summary_overall.csv'));
writetable(S_model, fullfile(OUT.table_dir, 'stage_B_summary_by_model.csv'));
writetable(S_phantom, fullfile(OUT.table_dir, 'stage_B_summary_by_phantom.csv'));
writetable(S_diam, fullfile(OUT.table_dir, 'stage_B_summary_by_inclusion_diameter.csv'));
writetable(S_incSws, fullfile(OUT.table_dir, 'stage_B_summary_by_inclusion_sws.csv'));
writetable(S_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_frequency.csv'));
writetable(S_field, fullfile(OUT.table_dir, 'stage_B_summary_by_field_regime.csv'));
writetable(S_Dlambda, fullfile(OUT.table_dir, 'stage_B_summary_by_D_over_lambda_bin.csv'));
writetable(S_DLwin, fullfile(OUT.table_dir, 'stage_B_summary_by_D_over_Lwin_bin.csv'));
writetable(S_Meff, fullfile(OUT.table_dir, 'stage_B_summary_by_Meff_bin.csv'));
writetable(S_dist, fullfile(OUT.table_dir, 'stage_B_summary_by_distance_over_window_bin.csv'));
writetable(S_roi, fullfile(OUT.table_dir, 'stage_B_summary_by_roi.csv'));
writetable(S_roi_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_frequency.csv'));
writetable(S_diam_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_diameter_frequency.csv'));
writetable(S_diam_field, fullfile(OUT.table_dir, 'stage_B_summary_by_diameter_field.csv'));
writetable(S_worst, fullfile(OUT.table_dir, 'stage_B_worst_conditions.csv'));
writetable(S_roi_diam_freq, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_diameter_frequency.csv'));
writetable(S_roi_diam_field, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_diameter_field.csv'));
writetable(S_roi_physics, fullfile(OUT.table_dir, 'stage_B_summary_by_roi_physics_bins.csv'));
writetable(S_condition_avg_diam, fullfile(OUT.table_dir, 'stage_B_condition_averaged_by_diameter.csv'));
writetable(S_condition_avg_roi, fullfile(OUT.table_dir, 'stage_B_condition_averaged_by_roi.csv'));
end

function make_analysis_figures(T, OUT)
models = unique(string(T.model_name), 'stable');
comp = "q_spectrum_plus_composition";

% 1. MAPE vs inclusion diameter.
S = summarize_metrics(T(ismember(T.roi_group,["inclusion_center_roi","interface"]),:), ...
    ["model_name","roi_group","inclusion_sws","inclusion_diameter_mm","field_regime"]);
fig = newfig([22 12]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
roiList = ["inclusion_center_roi","interface"];
diameterTicks = [10 14 24 32];
for ri = 1:numel(roiList)
    for si = 1:2
        ax = nexttile(tl); hold(ax,'on');
        incSws = [3 4];
        if si > numel(incSws), continue; end
        sws = incSws(si);
        Y = S(S.roi_group==roiList(ri) & S.inclusion_sws==sws & S.model_name==comp,:);
        fields = unique(Y.field_regime, 'stable');
        for fi = 1:numel(fields)
            Z = sortrows(Y(Y.field_regime==fields(fi),:), 'inclusion_diameter_mm');
            plot(ax, Z.inclusion_diameter_mm, Z.MAPE_pct, '-o', 'LineWidth',1.4, 'DisplayName', pretty_regime(fields(fi)));
        end
        title(ax, sprintf('%s, inclusion %.0f m/s', pretty_roi_group(roiList(ri)), sws), 'FontWeight','normal');
        xlabel(ax,'Inclusion diameter D (mm)'); ylabel(ax,'MAPE (%)'); grid(ax,'on');
        xticks(ax, diameterTicks); xlim(ax, [9 33]);
        if ri == 1 && si == 2, legend(ax,'Location','bestoutside'); end
    end
end
title(tl, 'MAPE vs inclusion diameter (q spectrum + composition)', 'FontWeight','normal', 'FontSize',11);
savefigs(fig, OUT.analysis_dir, 'stage_B_mape_vs_inclusion_diameter');
% 1b. Signed bias vs inclusion diameter with the same ROI layout.
fig = newfig([22 12]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for ri = 1:numel(roiList)
    for si = 1:2
        ax = nexttile(tl); hold(ax,'on');
        incSws = [3 4];
        sws = incSws(si);
        Y = S(S.roi_group==roiList(ri) & S.inclusion_sws==sws & S.model_name==comp,:);
        fields = unique(Y.field_regime, 'stable');
        for fi = 1:numel(fields)
            Z = sortrows(Y(Y.field_regime==fields(fi),:), 'inclusion_diameter_mm');
            plot(ax, Z.inclusion_diameter_mm, Z.signed_bias_pct, '-o', 'LineWidth',1.4, 'DisplayName', pretty_regime(fields(fi)));
        end
        yline(ax, 0, 'k-', 'LineWidth',0.8, 'DisplayName','zero bias');
        title(ax, sprintf('%s, inclusion %.0f m/s', pretty_roi_group(roiList(ri)), sws), 'FontWeight','normal');
        xlabel(ax,'Inclusion diameter D (mm)'); ylabel(ax,'Signed bias (%)'); grid(ax,'on');
        xticks(ax, diameterTicks); xlim(ax, [9 33]);
        if ri == 1 && si == 2, legend(ax,'Location','bestoutside'); end
    end
end
title(tl, 'Signed bias vs inclusion diameter (q spectrum + composition)', 'FontWeight','normal', 'FontSize',11);
savefigs(fig, OUT.analysis_dir, 'stage_B_bias_vs_inclusion_diameter_center_roi');

% 2. Error and signed bias vs physical resolution axes.
resolutionRois = ["inclusion_center_roi","inclusion_core","interface"];
Sphys = summarize_metrics(T(ismember(T.roi_group,resolutionRois),:), ...
    ["model_name","roi_group","inclusion_sws","inclusion_diameter_mm","f0","field_regime"]);
plot_metric_vs_resolution_axis(Sphys, OUT, comp, resolutionRois, 'mean_D_over_lambda', ...
    'D / \lambda_{inc}', 'MAPE_pct', 'MAPE (%)', ...
    'MAPE vs D / \lambda_{inc}', 'stage_B_mape_vs_D_over_lambda');
plot_metric_vs_resolution_axis(Sphys, OUT, comp, resolutionRois, 'mean_D_over_lambda', ...
    'D / \lambda_{inc}', 'signed_bias_pct', 'Signed bias (%)', ...
    'Signed bias vs D / \lambda_{inc}', 'stage_B_bias_vs_D_over_lambda');
plot_metric_vs_resolution_axis(Sphys, OUT, comp, resolutionRois, 'mean_D_over_Lwin', ...
    'D / L_{win}', 'MAPE_pct', 'MAPE (%)', ...
    'MAPE vs D / L_{win}', 'stage_B_mape_vs_D_over_Lwin');
plot_metric_vs_resolution_axis(Sphys, OUT, comp, resolutionRois, 'mean_D_over_Lwin', ...
    'D / L_{win}', 'signed_bias_pct', 'Signed bias (%)', ...
    'Signed bias vs D / L_{win}', 'stage_B_bias_vs_D_over_Lwin');

% 3b. Error and bias vs wavelengths per REQ window.
windowRois = ["background_far","inclusion_center_roi","inclusion_core","interface"];
Sw = summarize_metrics(T(ismember(T.roi_group,windowRois),:), ...
    ["model_name","roi_group","true_sws_level","field_regime"]);
plot_metric_vs_window_lambda(Sw, OUT, comp, windowRois, 'mean_lambda_over_Lwin', ...
    '\lambda / L_{win}', 'MAPE_pct', 'MAPE (%)', ...
    'MAPE vs \lambda / L_{win}', 'stage_B_mape_vs_lambda_over_Lwin_by_roi');
plot_metric_vs_window_lambda(Sw, OUT, comp, windowRois, 'mean_lambda_over_Lwin', ...
    '\lambda / L_{win}', 'signed_bias_pct', 'Signed bias (%)', ...
    'Signed bias vs \lambda / L_{win}', 'stage_B_bias_vs_lambda_over_Lwin_by_roi');
plot_metric_vs_window_lambda(Sw, OUT, comp, windowRois, 'mean_wavelengths_per_window', ...
    'L_{win} / \lambda', 'MAPE_pct', 'MAPE (%)', ...
    'MAPE vs wavelengths per window', 'stage_B_mape_vs_wavelengths_per_window_by_roi');
plot_metric_vs_window_lambda(Sw, OUT, comp, windowRois, 'mean_wavelengths_per_window', ...
    'L_{win} / \lambda', 'signed_bias_pct', 'Signed bias (%)', ...
    'Signed bias vs wavelengths per window', 'stage_B_bias_vs_wavelengths_per_window_by_roi');

% 4. Bias vs inclusion diameter.
S = summarize_metrics(T(T.roi_group=="inclusion_core",:), ["model_name","inclusion_sws","inclusion_diameter_mm","field_regime"]);
fig = newfig([17 9]); hold on;
for sws = unique(S.inclusion_sws).'
    for fi = 1:numel(unique(S.field_regime))
        fields = unique(S.field_regime,'stable');
        X = sortrows(S(S.model_name==comp & S.inclusion_sws==sws & S.field_regime==fields(fi),:), 'inclusion_diameter_mm');
        plot(X.inclusion_diameter_mm, X.signed_bias_pct, '-o', 'LineWidth',1.3, ...
            'DisplayName', sprintf('inc %.0f m/s, %s', sws, pretty_regime(fields(fi))));
    end
end
yline(0,'k-','Zero bias'); xlabel('Inclusion diameter D (mm)'); ylabel('Signed bias (%)'); xticks([10 14 24 32]); xlim([9 33]);
title('Core bias vs inclusion diameter', 'FontWeight','normal'); grid on; legend('Location','bestoutside');
savefigs(fig, OUT.analysis_dir, 'stage_B_bias_vs_inclusion_diameter');

% 5. Core vs interface error.
S = summarize_metrics(T, ["model_name","roi_label","inclusion_sws"]);
fig = newfig([18 9]);
Y = S(S.model_name==comp & ismember(S.roi_label, ["inclusion_center_roi","inclusion_core","background_far","interface_0_1mm","interface_1_2mm","interface_2_4mm","interface_4_8mm"]),:);
plot_grouped_bar(Y, 'roi_label', 'inclusion_sws', 'MAPE_pct', 'Core and interface-band MAPE', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_B_core_vs_interface_error');

% 6. Diameter x frequency heatmaps for composition.
for sws = unique(T.inclusion_sws).'
    X = T(T.model_name==comp & T.roi_group=="inclusion_core" & T.inclusion_sws==sws,:);
    if isempty(X), X = T(T.model_name==comp & T.inclusion_sws==sws,:); end
    Sx = summarize_metrics(X, ["inclusion_diameter_mm","f0"]);
    fig = newfig([11 8]);
    plot_heatmap(Sx, 'inclusion_diameter_mm', 'f0', 'MAPE_pct', sprintf('Core MAPE: inclusion %.0f m/s', sws), 'MAPE (%)');
    savefigs(fig, OUT.analysis_dir, sprintf('stage_B_heatmap_diameter_frequency_inc%.0f', sws));
end

% 7. Field regime comparison.
S = summarize_metrics(T, ["model_name","field_regime","roi_group"]);
fig = newfig([16 9]);
Y = S(S.model_name==comp & ismember(S.roi_group,["inclusion_center_roi","inclusion_core","interface","background_far"]),:);
plot_grouped_bar(Y, 'field_regime', 'roi_group', 'MAPE_pct', 'MAPE by field regime and ROI', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_B_field_regime_comparison');

% 8. Error vs distance/window radius.
S = summarize_metrics(T, ["model_name","distance_over_window_bin","inclusion_sws"]);
fig = newfig([16 9]); hold on;
for sws = unique(S.inclusion_sws).'
    X = S(S.model_name==comp & S.inclusion_sws==sws,:);
    [xv, ord] = ordered_bin_positions(X.distance_over_window_bin, ["0-0.25","0.25-0.5","0.5-1","1-2",">2"]);
    plot(xv, X.MAPE_pct(ord), '-o', 'LineWidth',1.4, 'DisplayName', sprintf('inclusion %.0f m/s', sws));
end
xticks(1:5); xticklabels(["0-0.25","0.25-0.5","0.5-1","1-2",">2"]);
xlabel('Distance to interface / window radius'); ylabel('MAPE (%)');
title('Error vs normalized interface distance', 'FontWeight','normal'); grid on; legend('Location','bestoutside');
savefigs(fig, OUT.analysis_dir, 'stage_B_error_vs_distance_over_window');

% 9. Error vs M_eff.
S = summarize_metrics(T, ["model_name","M_eff_bin","roi_group"]);
fig = newfig([16 9]);
Y = S(S.model_name==comp & ismember(S.roi_group,["inclusion_center_roi","inclusion_core","interface","background_far"]),:);
plot_grouped_bar(Y, 'M_eff_bin', 'roi_group', 'MAPE_pct', 'MAPE by M_{eff} bin and ROI', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_B_error_vs_Meff');

% 10. Predicted vs true ROI mean SWS.
S = summarize_metrics(T, ["model_name","condition_id","roi_label","inclusion_sws","inclusion_diameter_mm","f0","field_regime"]);
fig = newfig([12 10]); hold on;
colors = lines(numel(models));
for mi = 1:numel(models)
    X = S(string(S.model_name)==models(mi) & ismember(S.roi_label,["inclusion_center_roi","inclusion_core","background_far"]),:);
    scatter(X.mean_true_SWS, X.mean_predicted_SWS, 28, colors(mi,:), 'filled', ...
        'DisplayName', pretty_model(models(mi)), 'MarkerFaceAlpha',0.65);
end
lims = [min([S.mean_true_SWS; S.mean_predicted_SWS])-0.15, max([S.mean_true_SWS; S.mean_predicted_SWS])+0.15];
plot(lims, lims, 'k--', 'DisplayName','Ideal'); xlim(lims); ylim(lims); axis square;
xlabel('True ROI mean SWS (m/s)'); ylabel('Predicted ROI mean SWS (m/s)');
title('Predicted vs true ROI mean SWS', 'FontWeight','normal'); grid on; legend('Location','bestoutside');
savefigs(fig, OUT.analysis_dir, 'stage_B_predicted_vs_true_roi_mean_sws');
end

function make_extra_diagnostics(T, OUT)
comp = "q_spectrum_plus_composition";
X = T(T.model_name==comp,:);

% High-error rate by diameter and ROI.
S = summarize_metrics(X, ["roi_group","inclusion_sws","inclusion_diameter_mm"]);
fig = newfig([18 9]);
Y = S(ismember(S.roi_group,["inclusion_center_roi","inclusion_core","interface","background_far"]),:);
plot_grouped_bar(Y, 'inclusion_diameter_mm', 'roi_group', 'high_error20_pct', 'High-error rate by diameter and ROI', 'Pixels with error >20% (%)');
savefigs(fig, OUT.analysis_dir, 'stage_B_high20_by_diameter_roi');

% q prediction by frequency and inclusion SWS.
S = summarize_metrics(X, ["roi_group","inclusion_sws","f0"]);
fig = newfig([17 9]); hold on;
for rg = ["inclusion_core","background_far","interface"]
    for sws = unique(S.inclusion_sws).'
        Y = sortrows(S(S.roi_group==rg & S.inclusion_sws==sws,:), 'f0');
        if isempty(Y), continue; end
        plot(Y.f0, Y.mean_q_pred, '-o', 'LineWidth',1.2, ...
            'DisplayName', sprintf('%s, inc %.0f', pretty_roi_group(rg), sws));
    end
end
xlabel('Frequency (Hz)'); ylabel('Mean predicted q'); title('Predicted q vs frequency', 'FontWeight','normal');
grid on; legend('Location','bestoutside');
savefigs(fig, OUT.analysis_dir, 'stage_B_q_pred_vs_frequency_roi');

% Bias by ROI and true SWS level.
T2 = X;
T2.true_sws_level = strings(height(T2),1);
T2.true_sws_level(abs(T2.SWS_true-2)<0.05) = "2 m/s background";
T2.true_sws_level(abs(T2.SWS_true-3)<0.05) = "3 m/s inclusion";
T2.true_sws_level(abs(T2.SWS_true-4)<0.05) = "4 m/s inclusion";
T2.true_sws_level(T2.true_sws_level=="") = "other";
S = summarize_metrics(T2, ["roi_group","true_sws_level"]);
fig = newfig([15 8]);
Y = S(ismember(S.roi_group,["inclusion_center_roi","inclusion_core","interface","background_far"]),:);
plot_grouped_bar(Y, 'true_sws_level', 'roi_group', 'signed_bias_pct', 'Signed bias by ROI and true SWS', 'Signed bias (%)');
yline(0,'k-');
savefigs(fig, OUT.analysis_dir, 'stage_B_bias_by_roi_true_sws');
end

function plot_metric_vs_resolution_axis(S, OUT, modelName, roiList, xVar, xLabelText, yVar, yLabelText, titleText, fileStem)
S = S(S.model_name == modelName,:);
incSwsList = [3 4];
fig = newfig([24 16]);
tl = tiledlayout(fig, numel(roiList), numel(incSwsList), 'TileSpacing','compact', 'Padding','compact');
for ri = 1:numel(roiList)
    for si = 1:numel(incSwsList)
        ax = nexttile(tl); hold(ax,'on');
        sws = incSwsList(si);
        Y = S(S.roi_group == roiList(ri) & S.inclusion_sws == sws,:);
        fields = unique(Y.field_regime, 'stable');
        for fi = 1:numel(fields)
            Z = sortrows(Y(Y.field_regime == fields(fi),:), xVar);
            if isempty(Z), continue; end
            plot(ax, Z.(xVar), Z.(yVar), '-o', 'LineWidth',1.2, ...
                'MarkerSize',4.5, 'DisplayName', pretty_regime(fields(fi)));
        end
        if contains(string(yVar), "bias")
            yline(ax, 0, 'k-', 'LineWidth',0.7, 'HandleVisibility','off');
        end
        title(ax, sprintf('%s, inclusion %.0f m/s', pretty_roi_group(roiList(ri)), sws), ...
            'FontWeight','normal', 'FontSize',9);
        xlabel(ax, xLabelText); ylabel(ax, yLabelText); grid(ax,'on');
        xvals = unique(round(Y.(xVar), 2));
        if numel(xvals) <= 10 && ~isempty(xvals)
            xticks(ax, xvals);
        end
        if ri == 1 && si == numel(incSwsList)
            legend(ax, 'Location','bestoutside');
        end
    end
end
title(tl, titleText + " (q spectrum + composition)", 'FontWeight','normal', 'FontSize',11);
savefigs(fig, OUT.analysis_dir, fileStem + "_by_roi");
end

function plot_metric_vs_window_lambda(S, OUT, modelName, roiList, xVar, xLabelText, yVar, yLabelText, titleText, fileStem)
S = S(S.model_name == modelName,:);
fig = newfig([22 14]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
for ri = 1:numel(roiList)
    ax = nexttile(tl); hold(ax,'on');
    Y = S(S.roi_group == roiList(ri),:);
    fields = unique(Y.field_regime, 'stable');
    for fi = 1:numel(fields)
        Z = sortrows(Y(Y.field_regime == fields(fi),:), xVar);
        if isempty(Z), continue; end
        plot(ax, Z.(xVar), Z.(yVar), '-o', 'LineWidth',1.25, 'MarkerSize',5, ...
            'DisplayName', pretty_regime(fields(fi)));
        for zi = 1:height(Z)
            txt = erase(string(Z.true_sws_level(zi)), " m/s");
            txt = replace(txt, " background", ""); txt = replace(txt, " inclusion", "");
            text(ax, Z.(xVar)(zi), Z.(yVar)(zi), " " + txt, 'FontSize',7, ...
                'VerticalAlignment','middle', 'HandleVisibility','off');
        end
    end
    if contains(string(yVar), "bias")
        yline(ax, 0, 'k-', 'LineWidth',0.7, 'HandleVisibility','off');
    end
    title(ax, pretty_roi_group(roiList(ri)), 'FontWeight','normal', 'FontSize',9);
    xlabel(ax, xLabelText); ylabel(ax, yLabelText); grid(ax,'on');
    if ri == numel(roiList)
        legend(ax, 'Location','bestoutside');
    end
end
title(tl, titleText + " (q spectrum + composition)", 'FontWeight','normal', 'FontSize',11);
savefigs(fig, OUT.analysis_dir, fileStem);
end

function plot_grouped_bar(T, xVar, groupVar, yVar, ttl, ylab)
if isempty(T)
    text(0.5,0.5,'No data available','HorizontalAlignment','center'); axis off; return;
end
xCats = unique(string(T.(xVar)), 'stable');
gCats = unique(string(T.(groupVar)), 'stable');
Y = nan(numel(xCats), numel(gCats));
for xi = 1:numel(xCats)
    for gi = 1:numel(gCats)
        idx = string(T.(xVar))==xCats(xi) & string(T.(groupVar))==gCats(gi);
        if any(idx), Y(xi,gi) = T.(yVar)(find(idx,1)); end
    end
end
bar(categorical(arrayfun(@pretty_category, xCats)), Y); grid on;
xlabel(pretty_axis(xVar)); ylabel(ylab); title(ttl, 'FontWeight','normal');
legend(arrayfun(@pretty_category, gCats), 'Location','bestoutside'); xtickangle(20);
end

function plot_heatmap(T, rowVar, colVar, valueVar, ttl, cblabel)
if isempty(T)
    text(0.5,0.5,'No data available','HorizontalAlignment','center'); axis off; return;
end
rows = unique(T.(rowVar), 'stable');
cols = unique(T.(colVar), 'stable');
Z = nan(numel(rows), numel(cols));
for i = 1:numel(rows)
    for j = 1:numel(cols)
        idx = T.(rowVar)==rows(i) & T.(colVar)==cols(j);
        if any(idx), Z(i,j) = T.(valueVar)(find(idx,1)); end
    end
end
imagesc(Z); axis tight; grid off;
set(gca,'XTick',1:numel(cols),'XTickLabel',arrayfun(@pretty_category,string(cols)), ...
    'YTick',1:numel(rows),'YTickLabel',arrayfun(@pretty_category,string(rows)));
xtickangle(0); xlabel(pretty_axis(colVar)); ylabel(pretty_axis(rowVar));
title(ttl, 'FontWeight','normal'); cb = colorbar; ylabel(cb,cblabel);
textStrings = compose('%.1f', Z);
[xg,yg] = meshgrid(1:numel(cols),1:numel(rows));
text(xg(:), yg(:), textStrings(:), 'HorizontalAlignment','center', 'FontSize',7, 'Color','w');
end

function [xv, ord] = ordered_bin_positions(vals, order)
vals = string(vals);
ord = zeros(numel(order),1);
for i = 1:numel(order)
    idx = find(vals == order(i), 1);
    if isempty(idx), ord(i) = NaN; else, ord(i) = idx; end
end
keep = isfinite(ord); ord = ord(keep); xv = find(keep);
end

function fig = newfig(posCm)
fig = figure('Color','w','Units','centimeters','Position',[2 2 posCm]);
end

function savefigs(fig, outDir, stem)
ensure_dir(outDir);
set(findall(fig,'Type','axes'), 'FontName','Times New Roman');
try
    exportgraphics(fig, fullfile(outDir, stem + ".png"), 'Resolution', 220);
catch
    saveas(fig, fullfile(outDir, stem + ".png"));
end
savefig(fig, fullfile(outDir, stem + ".fig"));
close(fig);
end

function write_stage_readme(OUT, CFG, T, S_overall, S_diam, S_roi, S_worst)
fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w'); assert(fid > 0);
fprintf(fid, '# Stage B: Resolution Phantom / Geometry Transfer Validation\n\n');
fprintf(fid, '## Objective\n\n');
fprintf(fid, 'Evaluate whether frozen `baseline_minimal_v1` models reconstruct clean Eikonal heterogeneous phantoms when inclusion diameter is controlled relative to wavelength and REQ window size. No retraining, correction, reliability mask, or oracle q is used.\n\n');
fprintf(fid, '## Simulation Design\n\n');
fprintf(fid, '- Phantoms: B1 has 2 m/s background and 3 m/s inclusions; B2 has 2 m/s background and 4 m/s inclusions.\n');
fprintf(fid, '- Diameters: 10, 14, 24, and 32 mm in one 2x2 phantom. The layout is large enough to preserve valid-window ROI pixels after edge cropping.\n');
fprintf(fid, '- ROI note: `inclusion_center_roi` uses the central half-radius of each inclusion and is reported even when the REQ patch is mixed; `inclusion_core` remains the stricter high-purity core diagnostic.\n');
fprintf(fid, '- Frequencies: %s Hz. Realism: clean only. REQ: M=%s, cs_guess=%.1f m/s, TargetStepM=%.1f mm.\n', ...
    strjoin(string(CFG.FrequenciesHz), ', '), mat2str(CFG.MList), CFG.REQ.cs_guess, CFG.REQ.TargetStepM*1e3);
fprintf(fid, '- Field regimes: %s. Source motion/polarization is [0 0 1].\n\n', strjoin(CFG.FieldRegimes, ', '));
fprintf(fid, '## Physical Ratios\n\n');
fprintf(fid, '- `D/lambda_inc`: inclusion diameter divided by wavelength inside the inclusion.\n');
fprintf(fid, '- `D/Lwin`: inclusion diameter divided by the nominal REQ window length, where `Lwin = M*cs_guess/f0`.\n');
fprintf(fid, '- `M_eff = M*cs_guess/SWS_true`: local effective wavelength count. Hard inclusions have smaller M_eff when cs_guess remains fixed at 3 m/s.\n\n');
fprintf(fid, '## ROI Notes\n\n');
fprintf(fid, 'Core ROIs require enough distance from the interface and patch purity >= 0.95. Small inclusions may have few or no clean core patches at low frequency because valid REQ windows and physical mixing consume the interior.\n\n');
fprintf(fid, '## Generated Figures\n\n');
fprintf(fid, '- MAPE vs inclusion diameter, D/lambda_inc, D/Lwin, and M_eff.\n');
fprintf(fid, '- Core versus interface-band error.\n');
fprintf(fid, '- Diameter x frequency heatmaps.\n');
fprintf(fid, '- Field-regime comparison and predicted-vs-true ROI mean SWS.\n');
fprintf(fid, '- Representative condition maps from the runner.\n\n');
fprintf(fid, '## Summary: Overall\n\n');
write_markdown_table(fid, S_overall);
fprintf(fid, '\n## Summary: Diameter\n\n');
write_markdown_table(fid, S_diam);
fprintf(fid, '\n## Summary: ROI\n\n');
write_markdown_table(fid, S_roi);
fprintf(fid, '\n## Worst Conditions\n\n');
write_markdown_table(fid, S_worst(1:min(10,height(S_worst)),:));
fprintf(fid, '\n## Automatic Interpretation\n\n');
write_interpretation(fid, T);
fprintf(fid, '\n## How to Run\n\n');
fprintf(fid, 'Validation-only:\n\n```matlab\nsetenv(''ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'',''validate_only''); run(''experiments/runners/eikonal_validation/run_stage_B_resolution_phantom.m'')\n```\n\n');
fprintf(fid, 'Quick:\n\n```matlab\nsetenv(''ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'',''quick''); run(''experiments/runners/eikonal_validation/run_stage_B_resolution_phantom.m'')\n```\n\n');
fprintf(fid, 'Full:\n\n```matlab\nsetenv(''ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'',''full''); run(''experiments/runners/eikonal_validation/run_stage_B_resolution_phantom.m'')\n```\n\n');
fprintf(fid, 'Analysis:\n\n```matlab\nsetenv(''ADAPTIVE_REQ_EIKONAL_STAGE_B_MODE'',''full''); run(''experiments/analysis/eikonal_validation/analyze_stage_B_resolution_phantom.m'')\n```\n');
fclose(fid);
end

function write_interpretation(fid, T)
comp = T(T.model_name=="q_spectrum_plus_composition",:);
if isempty(comp)
    fprintf(fid, '- No composition-model rows were available.\n'); return;
end
core = comp(comp.roi_group=="inclusion_core",:);
iface = comp(comp.roi_group=="interface",:);
if ~isempty(core)
    S = summarize_metrics(core, ["inclusion_sws","inclusion_diameter_mm"]);
    for sws = unique(S.inclusion_sws).'
        X = sortrows(S(S.inclusion_sws==sws,:), 'inclusion_diameter_mm');
        if height(X) >= 2
            fprintf(fid, '- Inclusion %.0f m/s core MAPE changes from %.2f%% at D=%.0f mm to %.2f%% at D=%.0f mm.\n', ...
                sws, X.MAPE_pct(1), X.inclusion_diameter_mm(1), X.MAPE_pct(end), X.inclusion_diameter_mm(end));
        end
    end
end
if ~isempty(core) && ~isempty(iface)
    Sc = summarize_metrics(core, "model_name"); Si = summarize_metrics(iface, "model_name");
    fprintf(fid, '- Core MAPE is %.2f%%; interface-band MAPE is %.2f%% for q_spectrum_plus_composition.\n', Sc.MAPE_pct(1), Si.MAPE_pct(1));
end
Sbias = summarize_metrics(comp(comp.roi_group=="inclusion_core",:), ["inclusion_sws"]);
for i = 1:height(Sbias)
    fprintf(fid, '- Core signed bias for %.0f m/s inclusions: %.2f%%.\n', Sbias.inclusion_sws(i), Sbias.signed_bias_pct(i));
end
fprintf(fid, '- Interpret interface bands cautiously: they include true physical mixing, finite REQ windows, and multi-source field complexity.\n');
end

function write_docs_readme(root_dir, CFG)
docDir = fullfile(root_dir, 'docs', 'eikonal_validation'); ensure_dir(docDir);
fid = fopen(fullfile(docDir, 'stage_B_resolution_phantom.md'), 'w'); assert(fid > 0);
fprintf(fid, '# Stage B: Resolution Phantom / Geometry Transfer Validation\n\n');
fprintf(fid, 'Stage B evaluates the frozen `baseline_minimal_v1` models on clean Eikonal multi-inclusion phantoms with controlled inclusion diameter. It is designed to quantify core-versus-interface behavior as a function of `D/lambda_inc`, `D/Lwin`, and `M_eff`.\n\n');
fprintf(fid, 'Runner: `experiments/runners/eikonal_validation/run_stage_B_resolution_phantom.m`\n\n');
fprintf(fid, 'Analysis: `experiments/analysis/eikonal_validation/analyze_stage_B_resolution_phantom.m`\n\n');
fprintf(fid, 'Config: `%s`\n\n', CFG.ConfigPath);
fprintf(fid, 'Outputs: `outputs/eikonal_validation/stage_B_resolution_phantom/`\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T), fprintf(fid, '_No rows._\n'); return; end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
maxRows = min(height(T), 30);
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

function print_summary_block(T, titleText)
fprintf('\n%s\n', titleText);
showVars = intersect(string(T.Properties.VariableNames), ...
    ["model_name","inclusion_sws","inclusion_diameter_mm","roi_label", ...
    "N_valid_patches","MAPE_pct","signed_bias_pct","high_error20_pct"], 'stable');
disp(T(:, cellstr(showVars)));
end

function label = pretty_model(name)
switch string(name)
    case "q_spectrum_only", label = "q spectrum only";
    case "q_spectrum_plus_composition", label = "q spectrum + composition";
    otherwise, label = pretty_category(name);
end
end

function label = pretty_regime(name)
switch string(name)
    case "single_source_lateral", label = "single source";
    case "diffuse_like_8src_layout1", label = "diffuse-like 8 sources";
    case "diffuse_like_8src_layout2", label = "diffuse-like layout 2";
    otherwise, label = pretty_category(name);
end
end

function label = pretty_roi_group(name)
switch string(name)
    case "inclusion_core", label = "clean inclusion core";
    case "inclusion_center_roi", label = "inclusion center ROI";
    case "background_far", label = "background far";
    case "interface", label = "interface bands";
    case "other", label = "other";
    otherwise, label = pretty_category(name);
end
end

function label = pretty_axis(name)
switch string(name)
    case "inclusion_diameter_mm", label = "Inclusion diameter D (mm)";
    case "inclusion_sws", label = "Inclusion SWS (m/s)";
    case "field_regime", label = "Field regime";
    case "roi_label", label = "ROI";
    case "roi_group", label = "ROI group";
    case "f0", label = "Frequency (Hz)";
    case "D_over_lambda_bin", label = "D / lambda_inc bin";
    case "D_over_Lwin_bin", label = "D / Lwin bin";
    case "M_eff_bin", label = "M_eff bin";
    case "distance_over_window_bin", label = "Distance/window bin";
    case "true_sws_level", label = "True SWS level";
    otherwise, label = pretty_category(name);
end
end

function label = pretty_category(name)
label = string(name);
label = replace(label, "q_spectrum_plus_composition", "q spectrum + composition");
label = replace(label, "q_spectrum_only", "q spectrum only");
label = replace(label, "phantom_B1_bg2_inc3", "B1: 2 to 3 m/s");
label = replace(label, "phantom_B2_bg2_inc4", "B2: 2 to 4 m/s");
label = replace(label, "single_source_lateral", "single source");
label = replace(label, "diffuse_like_8src_layout1", "diffuse-like 8 sources");
label = replace(label, "inclusion_center_roi", "inclusion center ROI");
label = replace(label, "inclusion_core", "clean inclusion core");
label = replace(label, "background_far", "background far");
label = replace(label, "interface_0_1mm", "interface 0-1 mm");
label = replace(label, "interface_1_2mm", "interface 1-2 mm");
label = replace(label, "interface_2_4mm", "interface 2-4 mm");
label = replace(label, "interface_4_8mm", "interface 4-8 mm");
label = replace(label, "_", " ");
end

function ensure_dir(d)
if exist(d, 'dir') ~= 7, mkdir(d); end
end
