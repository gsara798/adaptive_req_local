%% analyze_test_47_test38_full_results_diagnostics.m
% Test 47: diagnostics for the saved Test 38 full result table.
%
% This script does not run REQ and does not train any model. It reads the
% held-out prediction CSV produced by Test 38 and makes diagnostic summaries
% for the frozen q-spectrum model family.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = default_config();
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir,'test47_configuration.json'));

fprintf('\nTest 47: Test 38 full results diagnostics\n');
fprintf('Mode: %s | source: %s\n', CFG.Mode, CFG.SourceCsv);
fprintf('No REQ recalculation and no training in this script.\n');

assert(exist(CFG.SourceCsv,'file') == 2, 'Missing Test 38 CSV: %s', CFG.SourceCsv);

T = load_test38_predictions(CFG);
T = add_test47_labels(T, CFG);

main_models = CFG.MainModels;
diagnostic_models = [CFG.MainModels "theory_discrete"];
T_main = T(ismember(T.model_name, main_models),:);
T_diag = T(ismember(T.model_name, diagnostic_models),:);

fprintf('Loaded %d prediction rows (%d main-model rows, %d conditions).\n', ...
    height(T), height(T_main), numel(unique(T.condition_key)));

S_overall = summarize_predictions(T_diag, "model_name");
S_main_overall = S_overall(ismember(S_overall.model_name, main_models),:);
S_by_freq = summarize_predictions(T_main, ["model_name","f0"]);
S_by_M = summarize_predictions(T_main, ["model_name","M"]);
S_by_dx = summarize_predictions(T_main, ["model_name","dx_mm"]);
S_by_regime = summarize_predictions(T_main, ["model_name","field_regime_ood"]);
S_by_family = summarize_predictions(T_main, ["model_name","case_family"]);
S_by_case = summarize_predictions(T_main, ["model_name","case_id"]);
S_by_purity = summarize_predictions(T_main, ["model_name","purity_bin"]);
S_by_distance = summarize_predictions(T_main, ["model_name","distance_bin_test47"]);
S_by_roi = summarize_predictions(T_main, ["model_name","case_family","roi_region"]);
S_by_condition = summarize_predictions(T_main, ["model_name","condition_key","case_id", ...
    "case_family","field_regime_ood","f0","M","dx_mm"]);
S_worst = top_rows(S_by_condition, "MAPE_pct", 40);
S_best = bottom_rows(S_by_condition, "MAPE_pct", 40);
S_roi_sws = summarize_roi_sws(T_main);

writetable(S_overall, fullfile(OUT.table_dir,'test47_model_summary_with_theory_diagnostic.csv'));
writetable(S_main_overall, fullfile(OUT.table_dir,'test47_model_summary_main_no_theory.csv'));
writetable(S_by_freq, fullfile(OUT.table_dir,'test47_summary_by_frequency.csv'));
writetable(S_by_M, fullfile(OUT.table_dir,'test47_summary_by_M.csv'));
writetable(S_by_dx, fullfile(OUT.table_dir,'test47_summary_by_dx.csv'));
writetable(S_by_regime, fullfile(OUT.table_dir,'test47_summary_by_regime.csv'));
writetable(S_by_family, fullfile(OUT.table_dir,'test47_summary_by_family.csv'));
writetable(S_by_case, fullfile(OUT.table_dir,'test47_summary_by_case.csv'));
writetable(S_by_purity, fullfile(OUT.table_dir,'test47_summary_by_purity_bin.csv'));
writetable(S_by_distance, fullfile(OUT.table_dir,'test47_summary_by_distance_bin.csv'));
writetable(S_by_roi, fullfile(OUT.table_dir,'test47_summary_by_roi_region.csv'));
writetable(S_roi_sws, fullfile(OUT.table_dir,'test47_roi_sws_statistics.csv'));
writetable(S_worst, fullfile(OUT.table_dir,'test47_worst_conditions.csv'));
writetable(S_best, fullfile(OUT.table_dir,'test47_best_conditions.csv'));

safe_plot(@() plot_model_ranking(S_main_overall, S_overall, OUT), 'model ranking');
safe_plot(@() plot_mape_vs_frequency_M(S_by_freq, S_by_M, OUT), 'MAPE vs frequency/M');
safe_plot(@() plot_mape_vs_dx(S_by_dx, OUT), 'MAPE vs dx');
safe_plot(@() plot_family_case_diagnostics(S_by_family, S_by_case, OUT), 'family/case diagnostics');
safe_plot(@() plot_distance_roi_diagnostics(S_by_distance, S_by_roi, OUT), 'distance/ROI diagnostics');
safe_plot(@() plot_roi_sws_stats(S_roi_sws, OUT), 'ROI SWS statistics');
safe_plot(@() plot_roi_overlay_examples(T_main, OUT, CFG), 'ROI overlay examples');

save(fullfile(OUT.data_dir,'test47_compact_results.mat'), ...
    'CFG','S_overall','S_main_overall','S_by_freq','S_by_M','S_by_dx', ...
    'S_by_regime','S_by_family','S_by_case','S_by_purity','S_by_distance', ...
    'S_by_roi','S_roi_sws','S_worst','S_best','-v7.3');

print_summary(S_main_overall, S_by_freq, S_by_M, S_by_family, S_by_purity, S_worst, OUT);

fprintf('\nTables: %s\nFigures: %s\nTest 47 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config()
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST47_MODE'))));
if mode == "", mode = "full"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST47_MODE must be quick or full.');
root_dir = setup_adaptive_req();
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.SourceCsv = fullfile(root_dir, 'outputs', 'test_38_velocity_field_diverse_q_training', ...
    'tables', 'test38_patch_level_predictions.csv');
source_csv = strtrim(string(getenv('ADAPTIVE_REQ_TEST47_SOURCE_CSV')));
if source_csv ~= ""
    CFG.SourceCsv = char(source_csv);
end
CFG.MaxRowsQuick = env_num('ADAPTIVE_REQ_TEST47_MAX_ROWS', 300000);
CFG.MainModels = ["q_spectrum_plus_composition","q_spectrum_only", ...
    "q_spectrum_plus_theory_composition"];
CFG.PrimaryModel = "q_spectrum_plus_composition";
CFG.CenterRoiSizeMm = 8;
CFG.SoftCoreMinDistanceMm = 8;
CFG.HardCoreMinDistanceMm = 4;
CFG.DistanceEdgesMm = [0 0.5 1 2 4 8 Inf];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_47_test38_full_results_diagnostics','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_47_test38_full_results_diagnostics');
end
OUT.table_dir = fullfile(OUT.root_dir,'tables');
OUT.figure_dir = fullfile(OUT.root_dir,'figures');
OUT.roi_dir = fullfile(OUT.figure_dir,'roi_overlays');
OUT.data_dir = fullfile(OUT.root_dir,'data');
dirs = string(struct2cell(OUT));
for d = dirs(:)'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function T = load_test38_predictions(CFG)
opts = detectImportOptions(CFG.SourceCsv, 'TextType','string');
want = ["dataset","map_iz","map_ix","x_center_m","z_center_m","M", ...
    "condition_key","case_id","case_family","field_regime","field_regime_ood", ...
    "f0","dx","dz","true_SWS","patch_purity","q_oracle","purity_bin", ...
    "distance_to_boundary_mm","distance_bin","predicted_patch_purity","p_mixed", ...
    "p_strong_mixed","model_name","q_pred","sws_pred","q_error", ...
    "sws_signed_error_pct","sws_abs_error_pct","high_error10","high_error20"];
opts.SelectedVariableNames = cellstr(intersect(want, string(opts.VariableNames), 'stable'));
for v = ["dataset","condition_key","case_id","case_family","field_regime", ...
        "field_regime_ood","purity_bin","distance_bin","model_name"]
    if ismember(v, string(opts.VariableNames))
        opts = setvartype(opts, char(v), 'string');
    end
end
T = readtable(CFG.SourceCsv, opts);
if CFG.QuickMode && height(T) > CFG.MaxRowsQuick
    rng(47001,'twister');
    idx = sort(randperm(height(T), CFG.MaxRowsQuick));
    T = T(idx,:);
end
T.dx_mm = 1e3*T.dx;
T.dz_mm = 1e3*T.dz;
T.high_error10 = logical(T.high_error10);
T.high_error20 = logical(T.high_error20);
end

function T = add_test47_labels(T, CFG)
T.material_region = repmat("unknown", height(T), 1);
T.roi_region = repmat("other", height(T), 1);
T.distance_bin_test47 = repmat("no_interface", height(T), 1);
T.is_center_roi = false(height(T),1);

[G, keys] = findgroups(T.condition_key);
[Gs, ord] = sort(G);
starts = find([true; diff(Gs) ~= 0]);
stops = [starts(2:end)-1; numel(Gs)];
for si = 1:numel(starts)
    gi = Gs(starts(si));
    idx = ord(starts(si):stops(si));
    X = T(idx,:);
    fam = string(X.case_family(1));
    cs = X.true_SWS;
    cmin = min(cs, [], 'omitnan');
    cmax = max(cs, [], 'omitnan');
    crange = cmax - cmin;
    mat = repmat("unknown", height(X), 1);
    if ~isfinite(crange) || crange < 0.05 || fam == "homogeneous"
        mat(:) = "homogeneous";
    elseif contains(fam, ["inclusion","ellipse","two_inclusions"])
        mat(:) = "soft_background";
        mat(cs >= cmin + 0.5*crange) = "hard_inclusion";
    elseif fam == "three_material"
        mat(:) = "soft";
        mat(cs >= cmin + crange/3 & cs < cmin + 2*crange/3) = "mid";
        mat(cs >= cmin + 2*crange/3) = "hard";
    else
        mat(:) = "soft";
        mat(cs >= cmin + 0.5*crange) = "hard";
    end
    d = X.distance_to_boundary_mm;
    db = bin_distance(d, CFG.DistanceEdgesMm);
    roi = repmat("other", height(X), 1);

    x0 = median(X.x_center_m, 'omitnan');
    z0 = median(X.z_center_m, 'omitnan');
    half = 0.5*CFG.CenterRoiSizeMm/1e3;
    center_roi = abs(X.x_center_m - x0) <= half & abs(X.z_center_m - z0) <= half;
    if fam == "homogeneous"
        roi(center_roi) = "homogeneous_center_8mm";
    else
        roi(d <= 0.5) = "interface_0_0p5mm";
        roi(d > 0.5 & d <= 1) = "interface_0p5_1mm";
        roi(d > 1 & d <= 2) = "interface_1_2mm";
        roi(d > 2 & d <= 4) = "interface_2_4mm";
        soft_like = ismember(mat, ["soft","soft_background"]);
        hard_like = ismember(mat, ["hard","hard_inclusion"]);
        roi(soft_like & d > CFG.SoftCoreMinDistanceMm) = "soft_core_gt8mm";
        roi(hard_like & d > CFG.HardCoreMinDistanceMm) = "hard_core_gt4mm";
        roi(mat == "mid" & d > CFG.HardCoreMinDistanceMm) = "mid_core_gt4mm";
    end
    T.material_region(idx) = mat;
    T.roi_region(idx) = roi;
    T.distance_bin_test47(idx) = db;
    T.is_center_roi(idx) = center_roi;
end
end

function b = bin_distance(d, edges)
b = repmat("no_interface", numel(d), 1);
labels = ["0-0.5 mm","0.5-1 mm","1-2 mm","2-4 mm","4-8 mm",">8 mm"];
for i = 1:numel(labels)
    hit = d >= edges(i) & d < edges(i+1);
    if i == numel(labels), hit = d >= edges(i); end
    b(hit) = labels(i);
end
b(~isfinite(d)) = "no_interface";
end

%% Summaries

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T)
    S = table();
    return;
end
[G, groups] = findgroups(T(:, group_vars));
S = groups;
S.N = splitapply(@numel, T.sws_abs_error_pct, G);
S.MAPE_pct = splitapply(@(x) mean(x,'omitnan'), T.sws_abs_error_pct, G);
S.median_abs_error_pct = splitapply(@(x) median(x,'omitnan'), T.sws_abs_error_pct, G);
S.mean_signed_error_pct = splitapply(@(x) mean(x,'omitnan'), T.sws_signed_error_pct, G);
S.median_signed_error_pct = splitapply(@(x) median(x,'omitnan'), T.sws_signed_error_pct, G);
S.underestimate_pct = 100*splitapply(@(x) mean(x < 0,'omitnan'), T.sws_signed_error_pct, G);
S.overestimate_pct = 100*splitapply(@(x) mean(x > 0,'omitnan'), T.sws_signed_error_pct, G);
S.high_error10_pct = 100*splitapply(@(x) mean(x,'omitnan'), double(T.high_error10), G);
S.high_error20_pct = 100*splitapply(@(x) mean(x,'omitnan'), double(T.high_error20), G);
S.mean_abs_q_error = splitapply(@(x) mean(abs(x),'omitnan'), T.q_error, G);
S.mean_pred_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
S.std_pred_sws = splitapply(@std_omitnan, T.sws_pred, G);
S.mean_true_sws = splitapply(@(x) mean(x,'omitnan'), T.true_SWS, G);
end

function S = summarize_roi_sws(T)
T = T(T.roi_region ~= "other",:);
S = summarize_predictions(T, ["model_name","case_family","case_id","f0","M","dx_mm", ...
    "field_regime_ood","material_region","roi_region"]);
end

function T = top_rows(T, metric, n)
T = sortrows(T, char(metric), 'descend');
T = T(1:min(n,height(T)),:);
end

function T = bottom_rows(T, metric, n)
T = sortrows(T, char(metric), 'ascend');
T = T(1:min(n,height(T)),:);
end

%% Plots

function plot_model_ranking(S_main, S_all, OUT)
S_main = sortrows(S_main,'MAPE_pct','ascend');
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 10]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile;
barh(ax, categorical(pretty_model(S_main.model_name)), S_main.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'Model');
title(ax,'Main model ranking (Theory excluded)','FontWeight','normal');
grid(ax,'on');
ax = nexttile;
S_all = sortrows(S_all,'MAPE_pct','ascend');
barh(ax, categorical(pretty_model(S_all.model_name)), S_all.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'Model');
title(ax,'Diagnostic ranking including Theory','FontWeight','normal');
grid(ax,'on');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_model_ranking.png'), 'Resolution', 220);
close(fig);
end

function plot_mape_vs_frequency_M(S_freq, S_M, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 28 11]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on');
for m = unique(S_freq.model_name,'stable')'
    X = sortrows(S_freq(S_freq.model_name == m,:), 'f0');
    plot(ax, X.f0, X.MAPE_pct, '-o', 'LineWidth',1.2, ...
        'DisplayName', pretty_model(m));
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)');
title(ax,'Error versus frequency','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','best','Interpreter','none');
ax = nexttile; hold(ax,'on');
for m = unique(S_M.model_name,'stable')'
    X = sortrows(S_M(S_M.model_name == m,:), 'M');
    plot(ax, X.M, X.MAPE_pct, '-o', 'LineWidth',1.2, ...
        'DisplayName', pretty_model(m));
end
xlabel(ax,'REQ M'); ylabel(ax,'MAPE (%)');
title(ax,'Error versus REQ window M','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','best','Interpreter','none');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_mape_vs_frequency_M.png'), 'Resolution', 220);
close(fig);
end

function plot_mape_vs_dx(S_dx, OUT)
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 10]);
ax = axes(fig); hold(ax,'on');
for m = unique(S_dx.model_name,'stable')'
    X = sortrows(S_dx(S_dx.model_name == m,:), 'dx_mm');
    plot(ax, X.dx_mm, X.MAPE_pct, '-o', 'LineWidth',1.2, ...
        'DisplayName', pretty_model(m));
end
xlabel(ax,'Spatial resolution dx (mm)'); ylabel(ax,'MAPE (%)');
title(ax,'Error versus spatial resolution','FontWeight','normal');
grid(ax,'on'); legend(ax,'Location','best','Interpreter','none');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_mape_vs_dx.png'), 'Resolution', 220);
close(fig);
end

function plot_family_case_diagnostics(S_family, S_case, OUT)
primary = "q_spectrum_plus_composition";
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile;
X = sortrows(S_family(S_family.model_name == primary,:), 'MAPE_pct', 'ascend');
barh(ax, categorical(pretty_region(X.case_family)), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'Geometry family');
title(ax,'Primary model error by family','FontWeight','normal');
grid(ax,'on');
ax = nexttile;
X = sortrows(S_case(S_case.model_name == primary,:), 'MAPE_pct', 'descend');
X = X(1:min(18,height(X)),:);
barh(ax, categorical(pretty_region(X.case_id)), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'Case');
title(ax,'Worst cases for primary model','FontWeight','normal');
grid(ax,'on');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_family_case_diagnostics.png'), 'Resolution', 220);
close(fig);
end

function plot_distance_roi_diagnostics(S_dist, S_roi, OUT)
primary = "q_spectrum_plus_composition";
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on');
order = ["0-0.5 mm","0.5-1 mm","1-2 mm","2-4 mm","4-8 mm",">8 mm","no_interface"];
for m = unique(S_dist.model_name,'stable')'
    X = S_dist(S_dist.model_name == m,:);
    [~,loc] = ismember(X.distance_bin_test47, order);
    [~,ord] = sort(loc);
    X = X(ord,:);
    plot(ax, categorical(X.distance_bin_test47, order, 'Ordinal', true), ...
        X.MAPE_pct, '-o', 'LineWidth',1.2, 'DisplayName', pretty_model(m));
end
xlabel(ax,'Distance to interface'); ylabel(ax,'MAPE (%)');
title(ax,'Error versus distance to interface','FontWeight','normal');
grid(ax,'on'); legend(ax,'Location','best','Interpreter','none');
ax = nexttile;
X = S_roi(S_roi.model_name == primary & S_roi.roi_region ~= "other",:);
X = aggregate_for_plot(X, "roi_region");
X = sortrows(X, 'MAPE_pct', 'descend');
barh(ax, categorical(pretty_region(X.roi_region)), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'ROI / zone');
title(ax,'Primary model error by ROI/core/interface','FontWeight','normal');
grid(ax,'on');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_distance_roi_diagnostics.png'), 'Resolution', 220);
close(fig);
end

function plot_roi_sws_stats(S_roi, OUT)
keep = ismember(S_roi.model_name, ["q_spectrum_plus_composition","q_spectrum_only"]);
X = S_roi(keep & ismember(S_roi.roi_region, ...
    ["homogeneous_center_8mm","soft_core_gt8mm","hard_core_gt4mm","mid_core_gt4mm"]),:);
if isempty(X), return; end
X = aggregate_for_plot(X, ["model_name","roi_region"]);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile;
G = categorical(pretty_region(X.roi_region) + " | " + pretty_model(X.model_name));
barh(ax, G, X.mean_pred_sws);
xlabel(ax,'Mean predicted SWS (m/s)'); ylabel(ax,'ROI and model');
title(ax,'Predicted SWS levels in core ROIs','FontWeight','normal');
grid(ax,'on');
ax = nexttile;
barh(ax, G, X.MAPE_pct);
xlabel(ax,'MAPE (%)'); ylabel(ax,'ROI and model');
title(ax,'Core ROI error','FontWeight','normal');
grid(ax,'on');
exportgraphics(fig, fullfile(OUT.figure_dir,'test47_roi_sws_statistics.png'), 'Resolution', 220);
close(fig);
end

function plot_roi_overlay_examples(T, OUT, CFG)
primary = T(T.model_name == CFG.PrimaryModel,:);
if isempty(primary), return; end
families = ["homogeneous","bilayer","inclusion","ellipse","thin_layer","two_inclusions","three_material"];
for fam = families
    Xfam = primary(primary.case_family == fam,:);
    if isempty(Xfam), continue; end
    key = Xfam.condition_key(1);
    X = primary(primary.condition_key == key,:);
    plot_single_roi_overlay(X, OUT, CFG);
end
end

function plot_single_roi_overlay(X, OUT, CFG)
[true_map,nz,nx] = rows_to_grid(X, X.true_SWS);
roi_map = rows_to_grid(X, double(categorical(X.roi_region)), nz, nx);
dist_map = rows_to_grid(X, X.distance_to_boundary_mm, nz, nx);
x_mm = 1e3*rows_to_grid(X, X.x_center_m, nz, nx);
z_mm = 1e3*rows_to_grid(X, X.z_center_m, nz, nx);
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 10]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile;
imagesc(ax, true_map); axis(ax,'image'); colorbar(ax); colormap(ax, turbo);
title(ax, 'True SWS with ROI overlay', 'FontWeight','normal');
hold(ax,'on');
if string(X.case_family(1)) == "homogeneous"
    cx = median(X.x_center_m,'omitnan'); cz = median(X.z_center_m,'omitnan');
    dx = median(diff(unique(sort(X.x_center_m))),'omitnan');
    dz = median(diff(unique(sort(X.z_center_m))),'omitnan');
    if ~isfinite(dx) || dx <= 0, dx = median(X.dx,'omitnan'); end
    if ~isfinite(dz) || dz <= 0, dz = median(X.dz,'omitnan'); end
    ix0 = median(X.map_ix,'omitnan'); iz0 = median(X.map_iz,'omitnan');
    w = CFG.CenterRoiSizeMm/1e3/dx; h = CFG.CenterRoiSizeMm/1e3/dz;
    rectangle(ax, 'Position', [ix0-w/2 iz0-h/2 w h], 'EdgeColor','w', ...
        'LineWidth',1.8, 'LineStyle','--');
    text(ax, ix0, iz0, sprintf('8 x 8 mm center ROI'), 'Color','w', ...
        'HorizontalAlignment','center','FontWeight','bold');
else
    contour(ax, dist_map <= 0.5, [1 1], 'w', 'LineWidth',1.4);
    contour(ax, dist_map > 4, [1 1], 'k', 'LineWidth',1.0);
    text(ax, 3, 5, 'white: 0-0.5 mm interface | black: core >4 mm', ...
        'Color','w','FontWeight','bold','BackgroundColor',[0 0 0 0.35]);
end
ax = nexttile;
imagesc(ax, roi_map); axis(ax,'image'); colorbar(ax);
title(ax, 'ROI labels used for summaries', 'FontWeight','normal');
subtitle(ax, 'homogeneous center 8x8 mm; soft core >8 mm; hard core >4 mm');
sgtitle(fig, sprintf('ROI diagnostic: %s', X.condition_key(1)), 'Interpreter','none');
out = fullfile(OUT.roi_dir, "test47_roi_overlay__" + sanitize(X.condition_key(1)) + ".png");
exportgraphics(fig, out, 'Resolution', 220);
close(fig);
end

%% Utilities

function [M,nz,nx] = rows_to_grid(T, v, nz, nx)
if nargin < 3
    nz = max(T.map_iz);
    nx = max(T.map_ix);
end
M = nan(nz,nx);
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
M(idx) = v;
end

function s = pretty_model(x)
x = string(x);
s = x;
s(x=="q_spectrum_plus_composition") = "q spectrum + composition";
s(x=="q_spectrum_only") = "q spectrum only";
s(x=="q_spectrum_plus_theory_composition") = "q spectrum + theory + composition";
s(x=="theory_discrete") = "Theory discrete";
s = strrep(s, "_", " ");
end

function s = pretty_region(x)
s = string(x);
s = strrep(s, "_", " ");
s = strrep(s, "0p5", "0.5");
s = strrep(s, "gt", ">");
end

function s = sanitize(x)
s = regexprep(string(x), '[^A-Za-z0-9_=-]+', '_');
s = char(s);
end

function y = std_omitnan(x)
x = x(isfinite(x));
if numel(x) <= 1, y = NaN; else, y = std(x); end
end

function safe_plot(fun, label)
try
    fun();
catch ME
    warning('Test47:PlotFailed', 'Skipping %s: %s', label, ME.message);
end
end

function S = aggregate_for_plot(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
S = groups;
S.N = splitapply(@(x) sum(x,'omitnan'), T.N, G);
S.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), T.MAPE_pct, T.N, G);
S.mean_signed_error_pct = splitapply(@(x,w) weighted_mean(x,w), T.mean_signed_error_pct, T.N, G);
S.high_error20_pct = splitapply(@(x,w) weighted_mean(x,w), T.high_error20_pct, T.N, G);
if ismember("mean_pred_sws", string(T.Properties.VariableNames))
    S.mean_pred_sws = splitapply(@(x,w) weighted_mean(x,w), T.mean_pred_sws, T.N, G);
else
    S.mean_pred_sws = nan(height(S),1);
end
end

function y = weighted_mean(x,w)
ok = isfinite(x) & isfinite(w) & w > 0;
if ~any(ok), y = NaN; else, y = sum(x(ok).*w(ok))/sum(w(ok)); end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w');
fprintf(fid,'%s\n',txt);
fclose(fid);
end

function x = env_num(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x), x = default; end
end

function print_summary(S_overall, S_freq, S_M, S_family, S_purity, S_worst, OUT)
primary = "q_spectrum_plus_composition";
P = S_overall(S_overall.model_name == primary,:);
fprintf('\n================ Test 47 summary ================\n');
if ~isempty(P)
    fprintf('Primary model (%s): MAPE %.2f%% | signed %.2f%% | high>20 %.2f%%.\n', ...
        primary, P.MAPE_pct, P.mean_signed_error_pct, P.high_error20_pct);
end
best = sortrows(S_overall, 'MAPE_pct', 'ascend');
if ~isempty(best)
    fprintf('Best main model: %s (MAPE %.2f%%).\n', ...
        pretty_model(best.model_name(1)), best.MAPE_pct(1));
end
show_minmax(S_freq(S_freq.model_name==primary,:), "f0", "frequency");
show_minmax(S_M(S_M.model_name==primary,:), "M", "M");
show_minmax(S_family(S_family.model_name==primary,:), "case_family", "geometry family");
X = S_purity(S_purity.model_name==primary,:);
if ~isempty(X)
    X = sortrows(X,'MAPE_pct','descend');
    fprintf('Worst purity bin: %s (MAPE %.2f%%).\n', pretty_region(X.purity_bin(1)), X.MAPE_pct(1));
end
if ~isempty(S_worst)
    W = S_worst(S_worst.model_name==primary,:);
    if ~isempty(W)
        fprintf('Worst primary condition: %s (MAPE %.2f%%, high>20 %.2f%%).\n', ...
            W.condition_key(1), W.MAPE_pct(1), W.high_error20_pct(1));
    end
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

function show_minmax(T, col, label)
if isempty(T), return; end
T = sortrows(T, 'MAPE_pct', 'ascend');
fprintf('Best %s: %s (MAPE %.2f%%). Worst %s: %s (MAPE %.2f%%).\n', ...
    label, string(T.(col)(1)), T.MAPE_pct(1), label, ...
    string(T.(col)(end)), T.MAPE_pct(end));
end
