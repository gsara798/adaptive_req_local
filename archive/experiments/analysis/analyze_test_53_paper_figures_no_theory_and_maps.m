%% analyze_test_53_paper_figures_no_theory_and_maps.m
% Post-hoc paper figures and geometry-organized maps for Test 53.
%
% This script does not train and does not run REQ. It reads the saved Test53
% prediction table from test38_results.mat, excludes theory_discrete from the
% main figures, and exports condition maps under:
%
%   outputs/test_53_paper_final_clean_q_training/figures_no_theory/maps_by_geometry/<case_id>/
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST53_POSTHOC_MAX_MAPS = Inf by default
%   ADAPTIVE_REQ_TEST53_POSTHOC_MAP_STYLE = interp | patch

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = struct();
CFG.SourceRoot = fullfile(root_dir, 'outputs', 'test_53_paper_final_clean_q_training');
CFG.ResultsFile = fullfile(CFG.SourceRoot, 'data', 'test38_results.mat');
CFG.MapStyle = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST53_POSTHOC_MAP_STYLE'))));
if CFG.MapStyle == "", CFG.MapStyle = "interp"; end
assert(ismember(CFG.MapStyle, ["interp","patch"]), ...
    'ADAPTIVE_REQ_TEST53_POSTHOC_MAP_STYLE must be interp or patch.');
CFG.MaxMaps = env_number('ADAPTIVE_REQ_TEST53_POSTHOC_MAX_MAPS', Inf);
CFG.MainModels = ["q_spectrum_only","q_spectrum_plus_composition", ...
    "q_spectrum_plus_theory_composition","delta_q_theory_composition"];
CFG.PrimaryModel = "q_spectrum_plus_composition";

OUT = struct();
OUT.root_dir = fullfile(CFG.SourceRoot, 'figures_no_theory');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'summary');
OUT.map_dir = fullfile(OUT.root_dir, 'maps_by_geometry');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end

fprintf('\nTest 53 post-hoc paper figures without theory_discrete\n');
fprintf('Source: %s\n', CFG.ResultsFile);
fprintf('Map style: %s | max maps: %s\n', CFG.MapStyle, string(CFG.MaxMaps));
assert(exist(CFG.ResultsFile,'file') == 2, 'Missing Test53 results file: %s', CFG.ResultsFile);

S = load(CFG.ResultsFile, 'T_held');
T = S.T_held;
clear S;
T.model_name = string(T.model_name);
T.case_id = string(T.case_id);
T.case_family = string(T.case_family);
T.field_regime_ood = string(T.field_regime_ood);
T.condition_key = string(T.condition_key);

Tmain = T(T.model_name ~= "theory_discrete", :);
Tmain = Tmain(ismember(Tmain.model_name, CFG.MainModels), :);
fprintf('Loaded %d long rows, using %d rows across %d non-theory models.\n', ...
    height(T), height(Tmain), numel(unique(Tmain.model_name)));

%% Summaries

T_overall = summarize_predictions(Tmain, "model_name");
T_by_M = summarize_predictions(Tmain, ["model_name","M"]);
T_by_frequency = summarize_predictions(Tmain, ["model_name","f0"]);
T_by_family = summarize_predictions(Tmain, ["model_name","case_family"]);
T_by_case = summarize_predictions(Tmain, ["model_name","case_id"]);
T_by_regime = summarize_predictions(Tmain, ["model_name","field_regime_ood"]);
T_by_purity = summarize_predictions(Tmain, ["model_name","purity_bin"]);
T_worst = top_failure_conditions(summarize_predictions(Tmain, ...
    ["model_name","condition_key","case_id","case_family","field_regime_ood","f0","M"]), 50);

writetable(T_overall, fullfile(OUT.table_dir, 'test53_no_theory_summary_overall.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_M.csv'));
writetable(T_by_frequency, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_frequency.csv'));
writetable(T_by_family, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_family.csv'));
writetable(T_by_case, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_case.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_regime.csv'));
writetable(T_by_purity, fullfile(OUT.table_dir, 'test53_no_theory_summary_by_purity.csv'));
writetable(T_worst, fullfile(OUT.table_dir, 'test53_no_theory_worst_conditions.csv'));

%% Figures

plot_model_ranking(T_overall, OUT);
plot_metric_lines(T_by_frequency, "f0", 'Frequency (Hz)', ...
    'MAPE versus frequency', fullfile(OUT.figure_dir, 'test53_no_theory_mape_vs_frequency.png'));
plot_metric_lines(T_by_M, "M", 'REQ M', ...
    'MAPE versus REQ M', fullfile(OUT.figure_dir, 'test53_no_theory_mape_vs_M.png'));
plot_grouped_bars(T_by_family, "case_family", ...
    'MAPE by geometry family', fullfile(OUT.figure_dir, 'test53_no_theory_mape_by_family.png'));
plot_grouped_bars(T_by_purity, "purity_bin", ...
    'MAPE by patch purity bin', fullfile(OUT.figure_dir, 'test53_no_theory_mape_by_purity.png'));
plot_heatmap(T_by_case, "case_id", fullfile(OUT.figure_dir, 'test53_no_theory_mape_by_case.png'));
plot_heatmap(T_by_regime, "field_regime_ood", fullfile(OUT.figure_dir, 'test53_no_theory_mape_by_regime.png'));

%% Maps

keys = unique(Tmain.condition_key, 'stable');
if isfinite(CFG.MaxMaps)
    keys = keys(1:min(numel(keys), CFG.MaxMaps));
end
fprintf('Exporting %d condition maps grouped by geometry...\n', numel(keys));
for ki = 1:numel(keys)
    Xall = Tmain(Tmain.condition_key == keys(ki), :);
    if isempty(Xall), continue; end
    try
        export_condition_map(Xall, keys(ki), OUT, CFG);
    catch ME
        warning('Test53Posthoc:MapFailed', 'Failed map %s: %s', keys(ki), ME.message);
    end
    if mod(ki,25) == 0 || ki == numel(keys)
        fprintf('  maps %d/%d\n', ki, numel(keys));
    end
end

print_summary(T_overall, T_by_frequency, T_by_M, T_by_family, T_worst, OUT);
fprintf('\nTables: %s\nFigures/maps: %s\nDone.\n', OUT.table_dir, OUT.root_dir);

%% Local functions

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G == gi,:);
    rows{gi} = table(height(X), ...
        mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), ...
        mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), ...
        mean(abs(X.q_error),'omitnan'), ...
        mean(X.predicted_patch_purity,'omitnan'), ...
        mean(X.p_mixed,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','median_signed_error_pct', ...
        'underestimate_pct','high_error10_pct','high_error20_pct', ...
        'mean_abs_q_error','mean_predicted_patch_purity','mean_p_mixed'});
end
S = [groups vertcat(rows{:})];
end

function T = top_failure_conditions(T_condition_summary, n)
T = sortrows(T_condition_summary, 'MAPE_pct', 'descend');
T = T(1:min(n,height(T)),:);
end

function plot_model_ranking(T_overall, OUT)
T_overall = sortrows(T_overall, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[2 2 22 10]);
yyaxis left
bar(1:height(T_overall), T_overall.MAPE_pct);
ylabel('MAPE (%)');
yyaxis right
plot(1:height(T_overall), T_overall.high_error20_pct, 'ko-', ...
    'LineWidth',1.0,'MarkerFaceColor','k');
ylabel('High-error >20% (%)');
set(gca,'XTick',1:height(T_overall),'XTickLabel',cellstr(pretty_model(T_overall.model_name)));
xtickangle(30); grid on;
title('Test 53 model ranking (Theory discrete excluded)','FontWeight','normal');
export_fig(fig, fullfile(OUT.figure_dir,'test53_no_theory_model_ranking.png'));
end

function plot_metric_lines(T, xvar, xlab, ttl, file)
models = unique(T.model_name, 'stable');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 11]);
hold on;
for m = models(:)'
    X = sortrows(T(T.model_name == m,:), xvar);
    plot(X.(xvar), X.MAPE_pct, '-o', 'LineWidth', 1.2, 'DisplayName', pretty_model(m));
end
grid on; xlabel(xlab); ylabel('MAPE (%)'); title(ttl,'FontWeight','normal');
legend('Location','bestoutside');
export_fig(fig, file);
end

function plot_grouped_bars(T, group_var, ttl, file)
models = unique(T.model_name, 'stable');
groups = unique(string(T.(group_var)), 'stable');
Y = nan(numel(groups), numel(models));
for i = 1:numel(groups)
    for j = 1:numel(models)
        idx = T.model_name == models(j) & string(T.(group_var)) == groups(i);
        if any(idx), Y(i,j) = T.MAPE_pct(find(idx,1)); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 12]);
bar(Y); grid on;
set(gca,'XTick',1:numel(groups),'XTickLabel',cellstr(pretty_label(groups)));
xtickangle(35); ylabel('MAPE (%)'); title(ttl,'FontWeight','normal');
legend(cellstr(pretty_model(models)), 'Location','bestoutside');
export_fig(fig, file);
end

function plot_heatmap(T, group_var, file)
models = unique(T.model_name, 'stable');
groups = unique(string(T.(group_var)), 'stable');
Z = nan(numel(models), numel(groups));
for i = 1:numel(models)
    for j = 1:numel(groups)
        idx = T.model_name == models(i) & string(T.(group_var)) == groups(j);
        if any(idx), Z(i,j) = T.MAPE_pct(find(idx,1)); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 11]);
imagesc(Z); axis tight;
set(gca,'YTick',1:numel(models),'YTickLabel',cellstr(pretty_model(models)), ...
    'XTick',1:numel(groups),'XTickLabel',cellstr(pretty_label(groups)));
xtickangle(45); colorbar; ylabel(colorbar,'MAPE (%)');
title(sprintf('MAPE by %s (Theory discrete excluded)', strrep(group_var,'_',' ')), ...
    'FontWeight','normal');
export_fig(fig, file);
end

function export_condition_map(Xall, key, OUT, CFG)
models = CFG.MainModels(ismember(CFG.MainModels, unique(Xall.model_name)));
if isempty(models), return; end
X0 = Xall(Xall.model_name == CFG.PrimaryModel,:);
if isempty(X0), X0 = Xall(Xall.model_name == models(1),:); end

fig = figure('Color','w','Units','centimeters','Position',[1 1 34 26]);
tl = tiledlayout(fig, 4, numel(models), 'TileSpacing','compact','Padding','compact');

plot_map(nexttile(tl), rows_to_grid(X0, X0.true_SWS, CFG.MapStyle), 'True SWS', 'm/s');
plot_map(nexttile(tl), rows_to_grid(X0, X0.predicted_patch_purity, CFG.MapStyle), ...
    'Predicted purity', 'probability');
plot_map(nexttile(tl), rows_to_grid(X0, X0.p_mixed, CFG.MapStyle), ...
    'Predicted mixedness', 'probability');
plot_map(nexttile(tl), rows_to_grid(X0, X0.patch_purity, CFG.MapStyle), ...
    'True patch purity', 'fraction');

for mi = 1:numel(models)
    X = Xall(Xall.model_name == models(mi),:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, CFG.MapStyle), ...
        pretty_model(models(mi)) + " SWS", 'm/s');
end
for mi = 1:numel(models)
    X = Xall(Xall.model_name == models(mi),:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, CFG.MapStyle), ...
        pretty_model(models(mi)) + " abs error", '%');
end
for mi = 1:numel(models)
    X = Xall(Xall.model_name == models(mi),:);
    plot_map(nexttile(tl), rows_to_grid(X, X.sws_signed_error_pct, CFG.MapStyle), ...
        pretty_model(models(mi)) + " signed error", '%');
end
title(tl, strrep(string(key),'_','\_'), 'FontWeight','normal');

case_dir = fullfile(OUT.map_dir, sanitize(X0.case_id(1)));
if exist(case_dir,'dir') ~= 7, mkdir(case_dir); end
fname = sprintf('test53_%s__%s.png', CFG.MapStyle, sanitize(key));
export_fig(fig, fullfile(case_dir, fname));
end

function G = rows_to_grid(T, values, style)
x = double(T.x_center_m(:));
z = double(T.z_center_m(:));
v = double(values(:));
ok = isfinite(x) & isfinite(z) & isfinite(v);
x = x(ok); z = z(ok); v = v(ok);
if isempty(v), G = NaN; return; end
xu = unique(x); zu = unique(z);
[Xq,Zq] = meshgrid(xu, zu);
switch string(style)
    case "patch"
        [~,ix] = ismember(x, xu);
        [~,iz] = ismember(z, zu);
        G = accumarray([iz ix], v, [numel(zu) numel(xu)], @mean, NaN);
    otherwise
        if numel(unique(x)) < 2 || numel(unique(z)) < 2
            G = NaN(numel(zu), numel(xu));
            return;
        end
        F = scatteredInterpolant(x, z, v, 'natural', 'nearest');
        G = F(Xq, Zq);
end
end

function plot_map(ax, A, ttl, cb)
imagesc(ax, A); axis(ax,'image'); axis(ax,'off');
title(ax, ttl, 'Interpreter','none','FontWeight','normal');
c = colorbar(ax); ylabel(c, cb);
end

function s = pretty_model(name)
name = string(name);
s = name;
s(name == "q_spectrum_only") = "q spectrum only";
s(name == "q_spectrum_plus_composition") = "q spectrum + composition";
s(name == "q_spectrum_plus_theory_composition") = "q spectrum + theory + composition";
s(name == "delta_q_theory_composition") = "delta q + theory + composition";
s(name == "delta_logk_theory_composition") = "delta log-k + theory + composition";
s(name == "q_spectrum_plus_theory") = "q spectrum + theory";
end

function s = pretty_label(name)
s = strrep(string(name), "_", " ");
s = strrep(s, "p", ".");
end

function print_summary(T_overall, T_by_frequency, T_by_M, T_by_family, T_worst, OUT)
T = sortrows(T_overall, 'MAPE_pct', 'ascend');
fprintf('\nInterpretive summary without theory_discrete:\n');
fprintf('  Best global model: %s, MAPE %.2f%%, high>20 %.2f%%.\n', ...
    T.model_name(1), T.MAPE_pct(1), T.high_error20_pct(1));
F = sortrows(T_by_frequency(T_by_frequency.model_name == T.model_name(1),:), 'MAPE_pct', 'descend');
if ~isempty(F)
    fprintf('  Hardest frequency for best model: %.0f Hz, MAPE %.2f%%.\n', ...
        F.f0(1), F.MAPE_pct(1));
end
M = sortrows(T_by_M(T_by_M.model_name == T.model_name(1),:), 'MAPE_pct', 'ascend');
if ~isempty(M)
    fprintf('  Best M for best model: M=%.1f, MAPE %.2f%%.\n', M.M(1), M.MAPE_pct(1));
end
Fam = sortrows(T_by_family(T_by_family.model_name == T.model_name(1),:), 'MAPE_pct', 'descend');
if ~isempty(Fam)
    fprintf('  Hardest family for best model: %s, MAPE %.2f%%.\n', ...
        Fam.case_family(1), Fam.MAPE_pct(1));
end
fprintf('  Worst condition overall: %s / %s / f%.0f / M%.1f, MAPE %.2f%%.\n', ...
    T_worst.case_id(1), T_worst.field_regime_ood(1), T_worst.f0(1), ...
    T_worst.M(1), T_worst.MAPE_pct(1));
fprintf('  Outputs: %s\n', OUT.root_dir);
end

function export_fig(fig, file)
drawnow;
try
    exportgraphics(fig, file, 'Resolution', 180);
catch
    saveas(fig, file);
end
close(fig);
end

function x = env_number(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x) && isfinite(default), x = default; end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end
