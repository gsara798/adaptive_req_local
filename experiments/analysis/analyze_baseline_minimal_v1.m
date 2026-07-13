%% analyze_baseline_minimal_v1.m
% Analysis-only companion for baseline_minimal_v1.
%
% This script never trains models and never recomputes REQ. It reads the
% tables produced by experiments/runners/run_baseline_minimal_v1.m, creates
% paper/debug figures, and writes a README_results.md file.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = analysis_config(root_dir);
OUT = output_dirs(root_dir, CFG);
assert(exist(OUT.table_dir, 'dir') == 7, 'Missing table directory: %s', OUT.table_dir);

fprintf('\nBaseline minimal v1 analysis\n');
fprintf('Output root: %s\n', OUT.root_dir);

T_overall = readtable(fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_overall.csv'));
T_freq = readtable(fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_frequency.csv'));
T_geom_freq = readtable(fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_geometry_frequency.csv'));
T_roi_freq = readtable(fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_roi_frequency.csv'));
T_purity = readtable(fullfile(OUT.table_dir, 'baseline_minimal_v1_summary_by_purity_bin.csv'));

fig_dir = fullfile(OUT.figure_dir, 'analysis');
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end

plot_overall_ranking(T_overall, fig_dir);
plot_frequency_curves(T_freq, fig_dir);
plot_geometry_frequency_heatmaps(T_geom_freq, fig_dir);
plot_roi_frequency_curves(T_roi_freq, fig_dir);
plot_purity_curves(T_purity, fig_dir);

write_results_readme(OUT, T_overall, T_freq, T_geom_freq, T_roi_freq, T_purity);

fprintf('Analysis figures: %s\n', fig_dir);
fprintf('README: %s\n', fullfile(OUT.root_dir, 'README_results.md'));
fprintf('Baseline minimal v1 analysis complete.\n');

function CFG = analysis_config(root_dir)
CFG = struct();
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_MODE'))));
if CFG.Mode == "", CFG.Mode = "full"; end
CFG.OutputName = lower(strtrim(string(getenv('ADAPTIVE_REQ_BASELINE_OUTPUT_NAME'))));
if CFG.OutputName == "", CFG.OutputName = "baseline_minimal_v1"; end
CFG.ConfigPath = fullfile(root_dir, 'configs', 'final_training', ...
    'baseline_minimal_v1.json');
end

function OUT = output_dirs(root_dir, CFG)
OUT.root_dir = fullfile(root_dir, 'outputs', char(CFG.OutputName));
if CFG.Mode == "quick"
    OUT.root_dir = fullfile(OUT.root_dir, 'quick');
elseif CFG.Mode == "medium"
    OUT.root_dir = fullfile(OUT.root_dir, 'medium');
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
end

function plot_overall_ranking(T, fig_dir)
T = sortrows(T, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 9]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
bar(ax, categorical(string(T.model_name)), T.MAPE_pct);
ylabel(ax, 'MAPE (%)'); title(ax, 'Overall MAPE', 'FontWeight','normal');
grid(ax, 'on'); xtickangle(ax, 25);
ax = nexttile(tl);
bar(ax, categorical(string(T.model_name)), T.high_error20_pct);
ylabel(ax, 'Pixels with error >20% (%)');
title(ax, 'High-error fraction', 'FontWeight','normal');
grid(ax, 'on'); xtickangle(ax, 25);
title(tl, 'Baseline minimal model ranking', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'baseline_minimal_v1_overall_ranking.png'));
end

function plot_frequency_curves(T, fig_dir)
models = unique(string(T.model_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    metrics = ["MAPE_pct","mean_signed_error_pct"];
ylabs = ["MAPE (%)","Mean signed error (%)"];
for pi = 1:2
    ax = nexttile(tl);
    hold(ax, 'on');
    for mi = 1:numel(models)
        X = sortrows(T(string(T.model_name)==models(mi),:), 'f0');
        plot(ax, X.f0, X.(metrics(pi)), '-o', 'LineWidth',1.2, ...
            'DisplayName', pretty_model(models(mi)));
    end
    xlabel(ax, 'Frequency (Hz)'); ylabel(ax, ylabs(pi));
    grid(ax, 'on'); legend(ax, 'Location','best');
end
title(tl, 'Frequency dependence', 'FontWeight','normal');
export_fig(fig, fullfile(fig_dir, 'baseline_minimal_v1_frequency_curves.png'));
end

function plot_geometry_frequency_heatmaps(T, fig_dir)
models = unique(string(T.model_name), 'stable');
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi),:);
    cases = unique(string(X.case_id), 'stable');
    freqs = unique(X.f0);
    Z = nan(numel(cases), numel(freqs));
    for ci = 1:numel(cases)
        for fi = 1:numel(freqs)
            idx = string(X.case_id)==cases(ci) & X.f0==freqs(fi);
            if any(idx), Z(ci,fi) = X.MAPE_pct(find(idx,1)); end
        end
    end
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 14]);
    imagesc(Z); axis tight;
    colormap(parula); cb = colorbar; ylabel(cb, 'MAPE (%)');
    set(gca, 'XTick',1:numel(freqs), 'XTickLabel', string(freqs), ...
        'YTick',1:numel(cases), 'YTickLabel', pretty_case(cases));
    xlabel('Frequency (Hz)'); ylabel('Geometry');
    title(sprintf('Geometry x frequency MAPE: %s', pretty_model(models(mi))), ...
        'FontWeight','normal', 'Interpreter','none');
    export_fig(fig, fullfile(fig_dir, ...
        "baseline_minimal_v1_geometry_frequency__" + sanitize(models(mi)) + ".png"));
end
end

function plot_roi_frequency_curves(T, fig_dir)
models = unique(string(T.model_name), 'stable');
rois = ["homogeneous_center","soft_core","intermediate_core","hard_core", ...
    "interface_0_1mm","interface_1_2mm","interface_2_4mm"];
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi) & ismember(string(T.roi_label), rois),:);
    if isempty(X), continue; end
    fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11]);
    tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    for pi = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        for ri = 1:numel(rois)
            R = sortrows(X(string(X.roi_label)==rois(ri),:), 'f0');
            if isempty(R), continue; end
            y = R.MAPE_pct;
        if pi == 2, y = R.mean_signed_error_pct; end
            plot(ax, R.f0, y, '-o', 'LineWidth',1.1, ...
                'DisplayName', pretty_roi(rois(ri)));
        end
        xlabel(ax, 'Frequency (Hz)');
        if pi == 1, ylabel(ax, 'MAPE (%)'); else, ylabel(ax, 'Mean signed error (%)'); end
        grid(ax, 'on'); legend(ax, 'Location','eastoutside');
    end
    title(tl, "ROI x frequency: " + pretty_model(models(mi)), ...
        'FontWeight','normal', 'Interpreter','none');
    export_fig(fig, fullfile(fig_dir, ...
        "baseline_minimal_v1_roi_frequency__" + sanitize(models(mi)) + ".png"));
end
end

function plot_purity_curves(T, fig_dir)
models = unique(string(T.model_name), 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 10]);
ax = axes(fig); hold(ax, 'on');
bins = unique(string(T.purity_bin), 'stable');
for mi = 1:numel(models)
    y = nan(size(bins));
    for bi = 1:numel(bins)
        idx = string(T.model_name)==models(mi) & string(T.purity_bin)==bins(bi);
        if any(idx), y(bi) = T.MAPE_pct(find(idx,1)); end
    end
    plot(ax, 1:numel(bins), y, '-o', 'LineWidth',1.2, ...
        'DisplayName', pretty_model(models(mi)));
end
set(ax, 'XTick',1:numel(bins), 'XTickLabel', pretty_roi(bins));
xtickangle(ax, 25);
ylabel(ax, 'MAPE (%)'); xlabel(ax, 'Patch purity bin');
title(ax, 'Error versus patch purity', 'FontWeight','normal');
grid(ax, 'on'); legend(ax, 'Location','best');
export_fig(fig, fullfile(fig_dir, 'baseline_minimal_v1_error_vs_patch_purity.png'));
end

function write_results_readme(OUT, T_overall, T_freq, T_geom_freq, T_roi_freq, T_purity)
path = fullfile(OUT.root_dir, 'README_results.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Baseline minimal v1 results\n\n');
fprintf(fid, 'This file was generated by `experiments/analysis/analyze_baseline_minimal_v1.m`.\n\n');
fprintf(fid, '## Models\n\n');
fprintf(fid, '- `q_spectrum_only`: bagged-tree regressor from clean operational spectral features to `q_oracle`.\n');
fprintf(fid, '- `q_spectrum_plus_composition`: same base features plus internally predicted `predicted_patch_purity`, `p_mixed`, and `p_strong_mixed`.\n');
fprintf(fid, '- The auxiliary composition models are trained only on the experiment train split and saved separately in `models/composition_auxiliary_models.mat`.\n\n');
fprintf(fid, '## Overall metrics\n\n');
write_markdown_table(fid, T_overall(:, intersect(string(T_overall.Properties.VariableNames), ...
    ["model_name","N","MAPE_pct","median_abs_error_pct","mean_signed_error_pct","high_error10_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## Frequency summary\n\n');
write_markdown_table(fid, T_freq(:, intersect(string(T_freq.Properties.VariableNames), ...
    ["model_name","f0","N","MAPE_pct","mean_signed_error_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## ROI x frequency\n\n');
write_markdown_table(fid, T_roi_freq(:, intersect(string(T_roi_freq.Properties.VariableNames), ...
    ["model_name","roi_label","f0","N","MAPE_pct","mean_signed_error_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## Geometry x frequency\n\n');
write_markdown_table(fid, T_geom_freq(:, intersect(string(T_geom_freq.Properties.VariableNames), ...
    ["model_name","case_id","f0","N","MAPE_pct","mean_signed_error_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## Patch purity\n\n');
write_markdown_table(fid, T_purity(:, intersect(string(T_purity.Properties.VariableNames), ...
    ["model_name","purity_bin","N","MAPE_pct","mean_signed_error_pct","high_error20_pct"], 'stable')));
fprintf(fid, '\n## Figures\n\n');
fprintf(fid, 'Analysis figures are saved under `figures/analysis/`.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T)
    fprintf(fid, '_No rows available._\n');
    return;
end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
max_rows = min(height(T), 40);
for i = 1:max_rows
    vals = strings(1, numel(names));
    for j = 1:numel(names)
        v = T.(names(j))(i);
        if isnumeric(v)
            vals(j) = sprintf('%.4g', v);
        elseif islogical(v)
            vals(j) = string(v);
        else
            vals(j) = string(v);
        end
    end
    fprintf(fid, '| %s |\n', strjoin(vals, ' | '));
end
if height(T) > max_rows
    fprintf(fid, '\n_Only first %d rows shown._\n', max_rows);
end
end

function s = pretty_model(s)
s = string(s);
s = replace(s, "q_spectrum_only", "q spectrum only");
s = replace(s, "q_spectrum_plus_composition", "q spectrum + composition");
end

function s = pretty_case(s)
s = replace(string(s), "_", " ");
end

function s = pretty_roi(s)
s = replace(string(s), "_", " ");
s = replace(s, "interface 0 1mm", "interface 0-1 mm");
s = replace(s, "interface 1 2mm", "interface 1-2 mm");
s = replace(s, "interface 2 4mm", "interface 2-4 mm");
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end

function export_fig(fig, path)
drawnow;
try
    exportgraphics(fig, path, 'Resolution', 220);
catch
    saveas(fig, path);
end
close(fig);
end
