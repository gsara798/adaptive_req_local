%% refresh_test_53_strong_splits_figures.m
% Rebuild paper/debug figures for Test 53 strong-split validation.
%
% This script reads the CSV outputs produced by
% analyze_test_53_strong_splits_q_training.m and creates cleaner figures
% without retraining any model.
%
% Runtime controls
%   ADAPTIVE_REQ_TEST53_STRONG_MODE = quick | full
%
% Outputs
%   outputs/test53_strong_splits/<mode>/figures_paper/

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST53_STRONG_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), ...
    'ADAPTIVE_REQ_TEST53_STRONG_MODE must be quick or full.');

OUT.root_dir = fullfile(root_dir, 'outputs', 'test53_strong_splits', char(mode));
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures_paper');
if exist(OUT.figure_dir, 'dir') ~= 7
    mkdir(OUT.figure_dir);
end

fprintf('Refreshing Test 53 strong-split figures from:\n%s\n', OUT.table_dir);
assert(exist(fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'), 'file') == 2, ...
    'Missing summary_metrics_by_fold.csv. Run Test 53 strong-splits first.');

T.summary = readtable(fullfile(OUT.table_dir, 'summary_metrics_by_fold.csv'));
T.split = readtable(fullfile(OUT.table_dir, 'summary_metrics_by_split.csv'));
T.frequency = readtable(fullfile(OUT.table_dir, 'metrics_by_frequency.csv'));
T.M = readtable(fullfile(OUT.table_dir, 'metrics_by_M.csv'));
T.geometry_family = readtable(fullfile(OUT.table_dir, 'metrics_by_geometry_family.csv'));
T.field_family = readtable(fullfile(OUT.table_dir, 'metrics_by_field_family.csv'));
T.purity_bin = readtable(fullfile(OUT.table_dir, 'metrics_by_patch_purity_bin.csv'));
T.purity_group = readtable(fullfile(OUT.table_dir, 'metrics_by_purity_group.csv'));
roi_file = fullfile(OUT.table_dir, 'metrics_by_roi.csv');
if exist(roi_file, 'file') == 2
    T.roi = readtable(roi_file);
else
    T.roi = table();
end

plot_strategy_ranking(T.summary, OUT, "MAPE_pct", ...
    "Strong-split MAPE, learned models", "MAPE (%)", ...
    'paper_test53_mape_by_split_learned.png', false);
plot_strategy_ranking(T.summary, OUT, "high_error20_pct", ...
    "Strong-split high-error rate, learned models", ...
    "Pixels with |error| >20% (%)", ...
    'paper_test53_high20_by_split_learned.png', false);
plot_strategy_ranking(T.summary, OUT, "MAPE_pct", ...
    "Strong-split MAPE with diagnostic baselines", "MAPE (%)", ...
    'paper_test53_mape_by_split_with_baselines.png', true);

plot_composition_gain(T.summary, OUT);
plot_metric_by_group(T.frequency, OUT, "frequency", "MAPE_pct", ...
    "Frequency dependence across strong splits", "Frequency (Hz)", ...
    "MAPE (%)", 'paper_test53_mape_by_frequency.png');
plot_metric_by_group(T.M, OUT, "M", "MAPE_pct", ...
    "Window-size dependence across strong splits", "REQ M", ...
    "MAPE (%)", 'paper_test53_mape_by_M.png');
plot_metric_by_group(T.geometry_family, OUT, "geometry family", "MAPE_pct", ...
    "Geometry-family error", "Geometry family", ...
    "MAPE (%)", 'paper_test53_mape_by_geometry_family.png');
plot_metric_by_group(T.field_family, OUT, "field family", "MAPE_pct", ...
    "Field-family error", "Field family", ...
    "MAPE (%)", 'paper_test53_mape_by_field_family.png');

plot_metric_by_group(T.purity_bin, OUT, "patch purity bin", "MAPE_pct", ...
    "Error versus true patch purity", "True patch-purity bin", ...
    "MAPE (%)", 'paper_test53_mape_by_patch_purity.png');
plot_metric_by_group(T.purity_bin, OUT, "patch purity bin", "high_error20_pct", ...
    "High-error rate versus true patch purity", "True patch-purity bin", ...
    "Pixels with |error| >20% (%)", 'paper_test53_high20_by_patch_purity.png');
plot_purity_group_summary(T.purity_group, OUT);

if ~isempty(T.roi)
    plot_roi_summary(T.roi, OUT);
    plot_roi_pred_vs_true(T.roi, OUT);
    plot_roi_bias_by_true_sws(T.roi, OUT);
end

write_figure_readme(OUT);
fprintf('Paper/debug figures written to:\n%s\n', OUT.figure_dir);

%% Local plotting helpers

function plot_strategy_ranking(T, OUT, metric, ttl, ylab, fname, include_baselines)
T = T(string(T.group_type) == "overall", :);
if include_baselines
    models = ["fixed_q_train_median","q_spectrum_only", ...
        "q_spectrum_plus_composition","q_spectrum_plus_theory_composition"];
else
    models = ["q_spectrum_only","q_spectrum_plus_composition", ...
        "q_spectrum_plus_theory_composition"];
end
T = T(ismember(string(T.model_name), models), :);
A = aggregate_metric(T, ["split_name","model_name"], metric);
splits = ordered_split_names(unique(string(A.split_name), 'stable'));
Y = nan(numel(splits), numel(models));
E = nan(numel(splits), numel(models));
for si = 1:numel(splits)
    for mi = 1:numel(models)
        idx = string(A.split_name) == splits(si) & string(A.model_name) == models(mi);
        if any(idx)
            Y(si, mi) = A.mean_value(find(idx, 1));
            E(si, mi) = A.std_value(find(idx, 1));
        end
    end
end

fig = figure('Color','w','Position',[80 80 1250 620]);
b = bar(categorical(pretty_split(splits)), Y, 'grouped');
hold on;
for mi = 1:numel(models)
    x = b(mi).XEndPoints;
    errorbar(x, Y(:,mi), E(:,mi), 'k.', 'LineWidth', 1.0, 'CapSize', 8);
end
grid on;
ylabel(ylab);
title(ttl);
legend(pretty_model(models), 'Location','northoutside', ...
    'Orientation','horizontal', 'Interpreter','none');
set(gca, 'FontSize', 12);
xtickangle(20);
save_figure(fig, OUT, fname);
end

function plot_composition_gain(T, OUT)
T = T(string(T.group_type) == "overall" & ...
    ismember(string(T.model_name), ["q_spectrum_only","q_spectrum_plus_composition"]), :);
wide = unstack(T(:, {'split_name','fold_key','model_name','MAPE_pct'}), ...
    'MAPE_pct', 'model_name');
if ~all(ismember(["q_spectrum_only","q_spectrum_plus_composition"], string(wide.Properties.VariableNames)))
    return;
end
wide.delta_mape = wide.q_spectrum_only - wide.q_spectrum_plus_composition;
A = aggregate_metric(renamevars(wide, 'delta_mape', 'MAPE_gain_pct'), ...
    "split_name", "MAPE_gain_pct");
splits = ordered_split_names(unique(string(A.split_name), 'stable'));
[~, ord] = ismember(splits, string(A.split_name));
A = A(ord,:);
fig = figure('Color','w','Position',[80 80 1050 520]);
bar(categorical(pretty_split(splits)), A.mean_value, 'FaceColor',[0.20 0.45 0.75]);
hold on;
errorbar(categorical(pretty_split(splits)), A.mean_value, A.std_value, ...
    'k.', 'LineWidth', 1.0, 'CapSize', 8);
yline(0, 'k-', 'LineWidth', 1.0);
grid on;
ylabel('MAPE gain versus q spectrum only (percentage points)');
title('Does composition help? Positive means lower MAPE');
xtickangle(20);
save_figure(fig, OUT, 'paper_test53_composition_gain_by_split.png');
end

function plot_metric_by_group(T, OUT, group_label, metric, ttl, xlab, ylab, fname)
models = ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"];
T = T(ismember(string(T.model_name), models), :);
A = aggregate_metric(T, ["model_name","group_value"], metric);
groups = sort_group_values(unique(string(A.group_value), 'stable'));
fig = figure('Color','w','Position',[80 80 1100 560]);
hold on;
for mi = 1:numel(models)
    X = A(string(A.model_name) == models(mi), :);
    [~, ord] = ismember(groups, string(X.group_value));
    ok = ord > 0;
    vals = nan(size(groups));
    errs = nan(size(groups));
    vals(ok) = X.mean_value(ord(ok));
    errs(ok) = X.std_value(ord(ok));
    errorbar(1:numel(groups), vals, errs, '-o', ...
        'LineWidth', 1.6, 'MarkerSize', 6, ...
        'DisplayName', pretty_model(models(mi)));
end
grid on;
xlim([0.5 numel(groups)+0.5]);
set(gca, 'XTick', 1:numel(groups), 'XTickLabel', pretty_group(groups));
xtickangle(25);
xlabel(xlab);
ylabel(ylab);
title(ttl);
legend('Location','northoutside', 'Orientation','horizontal', 'Interpreter','none');
text(0.01, 0.02, "Grouped over folds and available conditions", ...
    'Units','normalized', 'FontSize', 10, 'Color',[0.3 0.3 0.3]);
save_figure(fig, OUT, fname);
end

function plot_purity_group_summary(T, OUT)
models = ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"];
T = T(ismember(string(T.model_name), models), :);
A = aggregate_metric(T, ["model_name","group_value"], "MAPE_pct");
groups = ["pure_ge_0p95","mixed_lt_0p95","strong_mixed_lt_0p75"];
Y = nan(numel(groups), numel(models));
E = nan(numel(groups), numel(models));
for gi = 1:numel(groups)
    for mi = 1:numel(models)
        idx = string(A.group_value) == groups(gi) & string(A.model_name) == models(mi);
        if any(idx)
            Y(gi,mi) = A.mean_value(find(idx,1));
            E(gi,mi) = A.std_value(find(idx,1));
        end
    end
end
fig = figure('Color','w','Position',[80 80 1050 560]);
b = bar(categorical(pretty_group(groups)), Y, 'grouped');
hold on;
for mi = 1:numel(models)
    errorbar(b(mi).XEndPoints, Y(:,mi), E(:,mi), 'k.', 'LineWidth', 1.0, 'CapSize', 8);
end
grid on;
ylabel('MAPE (%)');
title('Pure versus mixed patch performance');
legend(pretty_model(models), 'Location','northoutside', ...
    'Orientation','horizontal', 'Interpreter','none');
save_figure(fig, OUT, 'paper_test53_pure_vs_mixed_mape.png');
end

function plot_roi_summary(T, OUT)
models = ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"];
T = T(ismember(string(T.model_name), models), :);
A = aggregate_metric(T, ["model_name","roi_type"], "MAPE_pct");
rois = unique(string(A.roi_type), 'stable');
Y = nan(numel(rois), numel(models));
E = nan(numel(rois), numel(models));
for ri = 1:numel(rois)
    for mi = 1:numel(models)
        idx = string(A.roi_type) == rois(ri) & string(A.model_name) == models(mi);
        if any(idx)
            Y(ri,mi) = A.mean_value(find(idx,1));
            row = find(idx, 1);
            E(ri,mi) = A.std_value(row) ./ sqrt(max(A.n_rows(row), 1));
        end
    end
end
fig = figure('Color','w','Position',[80 80 1150 560]);
b = bar(categorical(pretty_group(rois)), Y, 'grouped');
hold on;
for mi = 1:numel(models)
    errorbar(b(mi).XEndPoints, Y(:,mi), E(:,mi), 'k.', 'LineWidth', 1.0, 'CapSize', 8);
end
grid on;
ylabel('MAPE (%)');
title('Core and interface ROI performance');
yline(0, 'k-', 'LineWidth', 0.8);
legend(pretty_model(models), 'Location','northoutside', ...
    'Orientation','horizontal', 'Interpreter','none');
xtickangle(20);
save_figure(fig, OUT, 'paper_test53_roi_mape_grouped.png');
end

function plot_roi_pred_vs_true(T, OUT)
models = ["q_spectrum_only","q_spectrum_plus_composition"];
T = T(ismember(string(T.model_name), models), :);
fig = figure('Color','w','Position',[80 80 860 720]);
hold on;
colors = lines(numel(models));
for mi = 1:numel(models)
    X = T(string(T.model_name) == models(mi), :);
    scatter(X.true_SWS_mean, X.pred_SWS_mean, 18, colors(mi,:), ...
        'filled', 'MarkerFaceAlpha', 0.25, ...
        'DisplayName', pretty_model(models(mi)));
end
lims = [min([T.true_SWS_mean; T.pred_SWS_mean], [], 'omitnan'), ...
    max([T.true_SWS_mean; T.pred_SWS_mean], [], 'omitnan')];
plot(lims, lims, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Ideal');
axis equal;
xlim(lims); ylim(lims);
grid on;
xlabel('True mean SWS in ROI (m/s)');
ylabel('Predicted mean SWS in ROI (m/s)');
title('ROI mean SWS: predicted versus true');
legend('Location','northwest', 'Interpreter','none');
save_figure(fig, OUT, 'paper_test53_roi_predicted_vs_true_sws.png');
end

function plot_roi_bias_by_true_sws(T, OUT)
models = ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition"];
T = T(ismember(string(T.model_name), models), :);
T.true_sws_bin = round(T.true_SWS_mean * 4) / 4;
A = aggregate_metric(T, ["model_name","true_sws_bin"], "mean_bias_pct");
bins = sort(unique(A.true_sws_bin));
fig = figure('Color','w','Position',[80 80 1050 560]);
hold on;
for mi = 1:numel(models)
    X = A(string(A.model_name) == models(mi), :);
    [~, ord] = ismember(bins, X.true_sws_bin);
    ok = ord > 0;
    vals = nan(size(bins));
    errs = nan(size(bins));
    vals(ok) = X.mean_value(ord(ok));
    errs(ok) = X.std_value(ord(ok)) ./ sqrt(max(X.n_rows(ord(ok)), 1));
    errorbar(bins, vals, errs, '-o', 'LineWidth', 1.6, ...
        'MarkerSize', 6, 'DisplayName', pretty_model(models(mi)));
end
h0 = yline(0, 'k--', 'LineWidth', 1.0);
h0.Annotation.LegendInformation.IconDisplayStyle = 'off';
grid on;
xlabel('True mean SWS in ROI (m/s)');
ylabel('Mean signed SWS error in ROI (%)');
title('ROI bias by true SWS level');
legend('Location','northoutside', 'Orientation','horizontal', 'Interpreter','none');
save_figure(fig, OUT, 'paper_test53_roi_bias_by_true_sws.png');
end

function A = aggregate_metric(T, group_vars, metric)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
mean_value = splitapply(@(x) mean(x, 'omitnan'), T.(metric), G);
std_value = splitapply(@(x) std(x, 'omitnan'), T.(metric), G);
n_rows = splitapply(@numel, T.(metric), G);
A = [groups table(mean_value, std_value, n_rows)];
end

function save_figure(fig, OUT, fname)
set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 13);
exportgraphics(fig, fullfile(OUT.figure_dir, fname), 'Resolution', 220);
savefig(fig, fullfile(OUT.figure_dir, replace(fname, ".png", ".fig")));
close(fig);
end

function s = ordered_split_names(s)
order = ["grouped_condition_repeated","leave_one_frequency_out", ...
    "leave_one_geometry_family_out","field_regime_within_family_ood", ...
    "leave_field_family_out","leave_one_M_out"];
s = string(s);
[~, ia] = ismember(s, order);
[~, ord] = sortrows([ia(:)==0 ia(:)]);
s = s(ord);
end

function s = sort_group_values(s)
s = string(s);
x = str2double(s);
if all(isfinite(x))
    [~, ord] = sort(x);
    s = s(ord);
else
    s = s(:);
end
end

function out = pretty_split(s)
s = string(s);
out = s;
out(s == "grouped_condition_repeated") = "Grouped random";
out(s == "leave_one_frequency_out") = "Leave frequency";
out(s == "leave_one_geometry_family_out") = "Leave geometry";
out(s == "field_regime_within_family_ood") = "Regime OOD";
out(s == "leave_field_family_out") = "Leave field family";
out(s == "leave_one_M_out") = "Leave M";
end

function out = pretty_model(s)
s = string(s);
out = s;
out(s == "fixed_q_train_median") = "Fixed q";
out(s == "q_spectrum_only") = "Spectrum only";
out(s == "q_spectrum_plus_composition") = "Spectrum + composition";
out(s == "q_spectrum_plus_theory_composition") = "Spectrum + theory + composition";
out(s == "theory_discrete") = "Theory discrete";
end

function out = pretty_group(s)
s = string(s);
out = s;
out = replace(out, "_", " ");
out = replace(out, "pure ge 0p95", "Pure >=0.95");
out = replace(out, "mixed lt 0p95", "Mixed <0.95");
out = replace(out, "strong mixed lt 0p75", "Strong mixed <0.75");
out = replace(out, "homogeneous center", "Homogeneous center");
out = replace(out, "soft core", "Soft core");
out = replace(out, "hard core", "Hard core");
out = replace(out, "mid core", "Mid core");
out = replace(out, "interface 0 1mm", "Interface 0-1 mm");
end

function write_figure_readme(OUT)
fid = fopen(fullfile(OUT.figure_dir, 'README_figures.md'), 'w');
fprintf(fid, '# Test 53 Strong-Split Paper Figures\n\n');
fprintf(fid, 'These figures are regenerated from CSV tables only. No model was retrained.\n\n');
fprintf(fid, 'Key files:\n\n');
fprintf(fid, '- `paper_test53_mape_by_split_learned.png`\n');
fprintf(fid, '- `paper_test53_high20_by_split_learned.png`\n');
fprintf(fid, '- `paper_test53_composition_gain_by_split.png`\n');
fprintf(fid, '- `paper_test53_mape_by_frequency.png`\n');
fprintf(fid, '- `paper_test53_mape_by_patch_purity.png`\n');
fprintf(fid, '- `paper_test53_pure_vs_mixed_mape.png`\n');
fprintf(fid, '- `paper_test53_roi_mape_grouped.png`\n');
fprintf(fid, '- `paper_test53_roi_predicted_vs_true_sws.png`\n');
fprintf(fid, '- `paper_test53_roi_bias_by_true_sws.png`\n');
fclose(fid);
end
