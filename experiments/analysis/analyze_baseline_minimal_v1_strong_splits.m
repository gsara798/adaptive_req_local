%% analyze_baseline_minimal_v1_strong_splits.m
% Analysis-only companion for baseline_minimal_v1 strong splits.
%
% This script never trains models and never recomputes REQ. It reads tables
% produced by experiments/runners/run_baseline_minimal_v1_strong_splits.m,
% creates paper/debug figures, and refreshes README_results.md.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = analysis_config(root_dir);
OUT = output_dirs(root_dir, CFG);
assert(exist(OUT.table_dir, 'dir') == 7, 'Missing table directory: %s', OUT.table_dir);
fig_dir = fullfile(OUT.figure_dir, 'analysis');
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end

fprintf('\nBaseline minimal v1 strong-splits analysis\n');
fprintf('Reading tables from: %s\n', OUT.table_dir);

T_split = readtable(fullfile(OUT.table_dir, 'summary_metrics_by_split.csv'));
T_fold = readtable(fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'));
T_freq = readtable(fullfile(OUT.table_dir, 'metrics_by_frequency.csv'));
T_geomfam = readtable(fullfile(OUT.table_dir, 'metrics_by_geometry_family.csv'));
T_geomfreq = readtable(fullfile(OUT.table_dir, 'metrics_by_geometry_frequency.csv'));
T_fieldfam = readtable(fullfile(OUT.table_dir, 'metrics_by_field_family.csv'));
T_purity = readtable(fullfile(OUT.table_dir, 'metrics_by_patch_purity_bin.csv'));
T_purity_group = readtable(fullfile(OUT.table_dir, 'metrics_by_purity_group.csv'));
T_roi = readtable(fullfile(OUT.table_dir, 'metrics_by_roi.csv'));
T_roi_freq = readtable(fullfile(OUT.table_dir, 'metrics_by_roi_frequency.csv'));
q_scatter_path = fullfile(OUT.table_dir, 'q_scatter_sample.csv');
if exist(q_scatter_path, 'file') == 2
    T_q = readtable(q_scatter_path);
else
    T_q = table();
end

safe_plot(@() plot_split_ranking(T_split, fig_dir), 'split ranking');
safe_plot(@() plot_fold_boxplots(T_fold, fig_dir), 'fold boxplots');
safe_plot(@() plot_leave_one_frequency(T_fold, fig_dir), 'leave frequency');
safe_plot(@() plot_geometry_family_ood(T_fold, fig_dir), 'geometry family OOD');
safe_plot(@() plot_field_family_ood(T_fold, fig_dir), 'field family OOD');
safe_plot(@() plot_frequency_curves(T_freq, fig_dir), 'frequency curves');
safe_plot(@() plot_purity_curves(T_purity, fig_dir), 'purity curves');
safe_plot(@() plot_purity_group_box(T_purity_group, fig_dir), 'purity group');
safe_plot(@() plot_roi_curves(T_roi, fig_dir), 'ROI bars');
safe_plot(@() plot_roi_frequency(T_roi_freq, fig_dir), 'ROI frequency');
safe_plot(@() plot_geometry_frequency_heatmaps(T_geomfreq, fig_dir), 'geometry-frequency heatmaps');
safe_plot(@() plot_model_delta(T_fold, fig_dir), 'model delta');
if ~isempty(T_q)
    safe_plot(@() plot_q_scatter(T_q, fig_dir), 'q scatter');
end

write_analysis_readme(OUT, CFG, T_split, T_fold, T_purity_group, T_roi);

fprintf('Analysis figures: %s\n', fig_dir);
fprintf('README: %s\n', fullfile(OUT.root_dir, 'README_results.md'));
fprintf('Strong-splits analysis complete.\n');

function CFG = analysis_config(root_dir)
CFG = struct();
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_STRONG_SPLITS_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
CFG.OutputName = lower(strtrim(string(getenv('ADAPTIVE_REQ_STRONG_SPLITS_OUTPUT_NAME'))));
if CFG.OutputName == "", CFG.OutputName = "baseline_minimal_v1_strong_splits"; end
CFG.ConfigPath = fullfile(root_dir, 'configs', 'final_training', 'baseline_minimal_v1_strong_splits.json');
end

function OUT = output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', char(CFG.OutputName));
if CFG.Mode == "quick"
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
end

function plot_split_ranking(T, fig_dir)
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 10]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
metrics = ["MAPE_pct_mean","high_error20_pct_mean"];
ylabs = ["Mean fold MAPE (%)","Mean fold high-error >20% (%)"];
for pi = 1:2
    ax = nexttile(tl);
    X = T;
    cats = categorical(pretty_split(X.split_name) + " | " + pretty_model(X.model_name));
    bar(ax, cats, X.(metrics(pi)));
    ylabel(ax, ylabs(pi)); grid(ax,'on'); xtickangle(ax, 35);
end
title(tl, 'Strong-split model ranking', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_strategy_ranking.png'));
end

function plot_fold_boxplots(T, fig_dir)
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 10]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for pi = 1:2
    ax = nexttile(tl);
    if pi == 1, y = T.MAPE_pct; ylabel_text = 'MAPE (%)'; ttl = 'Fold MAPE distribution';
    else, y = T.high_error20_pct; ylabel_text = 'High-error >20% (%)'; ttl = 'Fold high-error distribution'; end
    g = categorical(pretty_split(T.split_name) + " | " + pretty_model(T.model_name));
    boxchart(ax, g, y); ylabel(ax, ylabel_text); title(ax, ttl, 'FontWeight','normal');
    grid(ax,'on'); xtickangle(ax, 35);
end
export_fig(fig, fullfile(fig_dir, 'strong_splits_fold_boxplots.png'));
end

function plot_leave_one_frequency(T, fig_dir)
X = T(string(T.split_name)=="leave_one_frequency_out",:);
if isempty(X), return; end
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
plot_fold_lines(nexttile(tl), X, 'MAPE_pct', 'MAPE (%)');
plot_fold_lines(nexttile(tl), X, 'mean_signed_error_pct', 'Mean signed error (%)');
title(tl, 'Leave-one-frequency-out', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_leave_one_frequency.png'));
end

function plot_geometry_family_ood(T, fig_dir)
X = T(string(T.split_name)=="leave_one_geometry_family_out",:);
if isempty(X), return; end
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
plot_fold_lines(nexttile(tl), X, 'MAPE_pct', 'MAPE (%)');
plot_fold_lines(nexttile(tl), X, 'high_error20_pct', 'High-error >20% (%)');
title(tl, 'Leave-one-geometry-family-out', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_leave_geometry_family.png'));
end

function plot_field_family_ood(T, fig_dir)
X = T(string(T.split_name)=="leave_field_family_out",:);
if isempty(X), return; end
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
plot_fold_lines(nexttile(tl), X, 'MAPE_pct', 'MAPE (%)');
plot_fold_lines(nexttile(tl), X, 'high_error20_pct', 'High-error >20% (%)');
title(tl, 'Leave-field-family-out', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_leave_field_family.png'));
end

function plot_fold_lines(ax, X, metric, ylabel_text)
models = unique(string(X.model_name), 'stable'); hold(ax,'on');
for mi = 1:numel(models)
    M = X(string(X.model_name)==models(mi),:);
    plot(ax, categorical(string(M.fold_key)), M.(metric), '-o', 'LineWidth',1.2, ...
        'DisplayName', pretty_model(models(mi)));
end
ylabel(ax, ylabel_text); grid(ax,'on'); legend(ax,'Location','best'); xtickangle(ax,25);
end

function plot_frequency_curves(T, fig_dir)
models = unique(string(T.model_name), 'stable');
splits = unique(string(T.split_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 21 13]);
tl = tiledlayout(fig,ceil(numel(splits)/2),2,'TileSpacing','compact','Padding','compact');
for si = 1:numel(splits)
    ax = nexttile(tl); hold(ax,'on');
    for mi = 1:numel(models)
        X = sortrows(T(string(T.split_name)==splits(si) & string(T.model_name)==models(mi),:), 'f0');
        if isempty(X), continue; end
        [G, f] = findgroups(X.f0);
        y = splitapply(@(v) mean(v,'omitnan'), X.MAPE_pct, G);
        plot(ax, f, y, '-o', 'LineWidth',1.1, 'DisplayName', pretty_model(models(mi)));
    end
    title(ax, pretty_split(splits(si)), 'FontWeight','normal'); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); grid(ax,'on');
end
legend(nexttile(tl,1), 'Location','best');
export_fig(fig, fullfile(fig_dir, 'strong_splits_mape_vs_frequency_by_split.png'));
end

function plot_purity_curves(T, fig_dir)
models = unique(string(T.model_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
ax = axes(fig); hold(ax,'on');
bins = unique(string(T.purity_bin), 'stable');
for mi = 1:numel(models)
    y = nan(size(bins));
    for bi = 1:numel(bins)
        idx = string(T.model_name)==models(mi) & string(T.purity_bin)==bins(bi);
        y(bi) = mean(T.MAPE_pct(idx), 'omitnan');
    end
    plot(ax, 1:numel(bins), y, '-o', 'LineWidth',1.2, 'DisplayName', pretty_model(models(mi)));
end
set(ax,'XTick',1:numel(bins),'XTickLabel',bins); xtickangle(ax,25);
ylabel(ax,'MAPE (%)'); xlabel(ax,'Patch purity bin'); grid(ax,'on'); legend(ax,'Location','best');
title(ax,'Error versus patch purity across strong splits','FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_error_vs_patch_purity.png'));
end

function plot_purity_group_box(T, fig_dir)
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 9]);
g = categorical(pretty_roi(T.purity_group) + " | " + pretty_model(T.model_name));
boxchart(g, T.MAPE_pct); grid on; xtickangle(30);
ylabel('MAPE (%)'); title('Pure versus mixed patches across folds', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_purity_group_boxplot.png'));
end

function plot_roi_curves(T, fig_dir)
X = T(ismember(string(T.roi_label), ["homogeneous_center","soft_core","intermediate_core","hard_core","interface_0_1mm","interface_1_2mm","interface_2_4mm"]),:);
if isempty(X), return; end
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 10]);
g = categorical(pretty_roi(X.roi_label) + " | " + pretty_model(X.model_name));
boxchart(g, X.MAPE_pct); grid on; xtickangle(35);
ylabel('MAPE (%)'); title('ROI error across strong splits', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_roi_mape_boxplot.png'));
end

function plot_roi_frequency(T, fig_dir)
models = unique(string(T.model_name), 'stable');
rois = ["homogeneous_center","soft_core","intermediate_core","hard_core","interface_0_1mm","interface_1_2mm","interface_2_4mm"];
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi) & ismember(string(T.roi_label), rois),:);
    if isempty(X), continue; end
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11]);
    ax = axes(fig); hold(ax,'on');
    for ri = 1:numel(rois)
        R = X(string(X.roi_label)==rois(ri),:);
        if isempty(R), continue; end
        [G,f] = findgroups(R.f0);
        y = splitapply(@(v) mean(v,'omitnan'), R.MAPE_pct, G);
        plot(ax, f, y, '-o', 'LineWidth',1.1, 'DisplayName', pretty_roi(rois(ri)));
    end
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); grid(ax,'on'); legend(ax,'Location','eastoutside');
    title(ax, "ROI x frequency: " + pretty_model(models(mi)), 'FontWeight','normal', 'Interpreter','none');
    export_fig(fig, fullfile(fig_dir, "strong_splits_roi_frequency__" + sanitize(models(mi)) + ".png"));
end
end

function plot_geometry_frequency_heatmaps(T, fig_dir)
models = unique(string(T.model_name), 'stable');
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi),:);
    cases = unique(string(X.case_id), 'stable'); freqs = unique(X.f0);
    Z = nan(numel(cases), numel(freqs));
    for ci = 1:numel(cases)
        for fi = 1:numel(freqs)
            idx = string(X.case_id)==cases(ci) & X.f0==freqs(fi);
            Z(ci,fi) = mean(X.MAPE_pct(idx), 'omitnan');
        end
    end
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 15]);
    imagesc(Z); axis tight; colormap(parula); cb=colorbar; ylabel(cb,'MAPE (%)');
    set(gca,'XTick',1:numel(freqs),'XTickLabel',string(freqs), 'YTick',1:numel(cases),'YTickLabel',pretty_case(cases));
    xlabel('Frequency (Hz)'); ylabel('Geometry'); title("Geometry x frequency MAPE: " + pretty_model(models(mi)), 'FontWeight','normal');
    export_fig(fig, fullfile(fig_dir, "strong_splits_geometry_frequency__" + sanitize(models(mi)) + ".png"));
end
end

function plot_model_delta(T, fig_dir)
only = T(string(T.model_name)=="q_spectrum_only",:);
comp = T(string(T.model_name)=="q_spectrum_plus_composition",:);
if isempty(only) || isempty(comp), return; end

keys_only = string(only.split_name) + "||" + string(only.fold_key);
keys_comp = string(comp.split_name) + "||" + string(comp.fold_key);
[tf, loc] = ismember(keys_only, keys_comp);
if ~any(tf), return; end

labels = pretty_split(only.split_name(tf)) + " | " + string(only.fold_key(tf));
mape_gain = only.MAPE_pct(tf) - comp.MAPE_pct(loc(tf));
high20_gain = only.high_error20_pct(tf) - comp.high_error20_pct(loc(tf));

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 9]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); bar(ax, categorical(labels), mape_gain);
ylabel(ax,'MAPE gain from composition (points)'); grid(ax,'on'); xtickangle(ax,35); yline(ax,0,'k-');
ax = nexttile(tl); bar(ax, categorical(labels), high20_gain);
ylabel(ax,'High-error >20% gain (points)'); grid(ax,'on'); xtickangle(ax,35); yline(ax,0,'k-');
title(tl,'Composition model gain over q-spectrum-only', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'strong_splits_composition_gain.png'));
end

function plot_q_scatter(T, fig_dir)
models = unique(string(T.model_name), 'stable');
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi),:);
    if isempty(X), continue; end
    if height(X) > 50000
        rng(1); X = X(randperm(height(X), 50000),:);
    end
    fig = figure('Color','w','Units','centimeters','Position',[2 2 14 12]);
    scatter(X.q_oracle, X.q_pred, 8, X.sws_abs_error_pct, 'filled', 'MarkerFaceAlpha',0.25);
    axis square; grid on; xlim([0 1]); ylim([0 1]); hold on; plot([0 1],[0 1],'k--');
    cb=colorbar; ylabel(cb,'SWS absolute error (%)'); xlabel('Oracle REQ quantile q'); ylabel('Predicted REQ quantile q');
    title("q prediction scatter: " + pretty_model(models(mi)), 'FontWeight','normal');
    export_fig(fig, fullfile(fig_dir, "strong_splits_q_scatter__" + sanitize(models(mi)) + ".png"));
end
end

function write_analysis_readme(OUT, CFG, T_split, T_fold, T_purity_group, T_roi)
path = fullfile(OUT.root_dir, 'README_results.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Baseline minimal v1 strong-splits results\n\n');
fprintf(fid, 'Generated/refreshed by `experiments/analysis/analyze_baseline_minimal_v1_strong_splits.m`.\n\n');
fprintf(fid, '## What was run\n\n');
fprintf(fid, '- Mode: `%s`\n', CFG.Mode);
fprintf(fid, '- Runner output folder: `%s`\n', OUT.root_dir);
fprintf(fid, '- Models: `q_spectrum_only`, `q_spectrum_plus_composition`.\n');
fprintf(fid, '- Composition auxiliary models are trained fold-internally.\n\n');
fprintf(fid, '## Summary by split\n\n');
write_markdown_table(fid, T_split);
fprintf(fid, '\n## Summary by fold\n\n');
write_markdown_table(fid, T_fold(:, intersect(string(T_fold.Properties.VariableNames), ...
    ["split_name","fold_key","model_name","N","MAPE_pct","mean_signed_error_pct","high_error20_pct","mean_abs_q_error"], 'stable')));
fprintf(fid, '\n## Purity groups\n\n');
write_markdown_table(fid, T_purity_group(:, intersect(string(T_purity_group.Properties.VariableNames), ...
    ["split_name","fold_key","model_name","purity_group","N","MAPE_pct","mean_signed_error_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## ROI summary\n\n');
write_markdown_table(fid, T_roi(:, intersect(string(T_roi.Properties.VariableNames), ...
    ["split_name","fold_key","model_name","roi_label","N","MAPE_pct","mean_signed_error_pct","high_error20_pct","mean_true_SWS","mean_pred_SWS"], 'stable')));
fprintf(fid, '\n## Figures\n\nAnalysis figures are saved under `figures/analysis/`.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T), fprintf(fid, '_No rows available._\n'); return; end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
max_rows = min(height(T), 80);
for i = 1:max_rows
    vals = strings(1, numel(names));
    for j = 1:numel(names)
        v = T.(names(j))(i);
        if isnumeric(v), vals(j) = sprintf('%.4g', v); else, vals(j) = string(v); end
    end
    fprintf(fid, '| %s |\n', strjoin(vals, ' | '));
end
if height(T) > max_rows, fprintf(fid, '\n_Only first %d rows shown._\n', max_rows); end
end

function plot_fold_lines_dummy %#ok<DEFNU>
end

function s = pretty_model(s)
s = string(s);
s = replace(s, "q_spectrum_only", "q spectrum only");
s = replace(s, "q_spectrum_plus_composition", "q spectrum + composition");
end

function s = pretty_split(s)
s = replace(string(s), "grouped_condition_repeated", "Grouped condition");
s = replace(s, "leave_one_frequency_out", "Leave frequency out");
s = replace(s, "leave_one_geometry_family_out", "Leave geometry family out");
s = replace(s, "field_regime_within_family_ood", "Field regime OOD");
s = replace(s, "leave_field_family_out", "Leave field family out");
s = replace(s, "leave_one_M_out", "Leave M out");
s = replace(s, "_", " ");
end

function s = pretty_roi(s)
s = replace(string(s), "_", " ");
s = replace(s, "interface 0 1mm", "interface 0-1 mm");
s = replace(s, "interface 1 2mm", "interface 1-2 mm");
s = replace(s, "interface 2 4mm", "interface 2-4 mm");
end

function s = pretty_case(s)
s = replace(string(s), "_", " ");
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function export_fig(fig, path)
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
close(fig);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('StrongSplitsAnalysis:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
end
end
