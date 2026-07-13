%% analyze_stage_A_homogeneous_transfer.m
% Analysis for Stage A homogeneous Eikonal transfer validation.
%
% This script reads runner outputs only. It does not rerun simulations, does
% not retrain any model, and does not apply corrections or reliability masks.

clear; clc; close all;
format compact;

%% Setup
this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(fileparts(this_file))));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = load_stage_config(root_dir);
CFG.Mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'))));
if CFG.Mode == "", CFG.Mode = "quick"; end
OUT = make_output_dirs(root_dir, CFG);

fprintf('\nStage A analysis: homogeneous Eikonal transfer\n');
fprintf('Mode label: %s\nOutput root: %s\n', CFG.Mode, OUT.root_dir);

patchCsv = fullfile(OUT.table_dir, 'stage_A_patch_level_results.csv');
assert(exist(patchCsv, 'file') == 2, ['Missing patch-level table. Run the Stage A runner first: ' patchCsv]);
T = readtable(patchCsv, 'TextType','string');
T = normalize_loaded_table(T);
assert(height(T) > 0, 'Patch-level table is empty.');
assert(all(isfinite(T.SWS_pred)) && all(isfinite(T.q_pred)), 'Non-finite predictions found.');

%% Summaries
S_overall = summarize_metrics(T, "model_name");
S_model = S_overall;
S_geom = summarize_metrics(T, ["model_name","geometry"]);
S_freq = summarize_metrics(T, ["model_name","f0"]);
S_field = summarize_metrics(T, ["model_name","field_regime"]);
S_geom_freq = summarize_metrics(T, ["model_name","geometry","f0"]);
S_geom_field = summarize_metrics(T, ["model_name","geometry","field_regime"]);
S_condition = summarize_metrics(T, ["model_name","condition_id","geometry","true_sws","f0","field_regime"]);
S_worst = sortrows(S_condition, 'MAPE_pct', 'descend');
S_worst = S_worst(1:min(30,height(S_worst)),:);

writetable(S_overall, fullfile(OUT.table_dir, 'stage_A_summary_overall.csv'));
writetable(S_model, fullfile(OUT.table_dir, 'stage_A_summary_by_model.csv'));
writetable(S_geom, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry.csv'));
writetable(S_freq, fullfile(OUT.table_dir, 'stage_A_summary_by_frequency.csv'));
writetable(S_field, fullfile(OUT.table_dir, 'stage_A_summary_by_field_regime.csv'));
writetable(S_geom_freq, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry_frequency.csv'));
writetable(S_geom_field, fullfile(OUT.table_dir, 'stage_A_summary_by_geometry_field.csv'));
writetable(S_worst, fullfile(OUT.table_dir, 'stage_A_worst_conditions.csv'));

%% Figures
make_summary_figures(T, S_overall, S_geom, S_freq, S_field, S_geom_freq, S_condition, OUT);
make_additional_diagnostics(T, S_condition, OUT);
% Regenerate the analysis plots with publication-friendly labels, smaller
% titles, and larger margins. This intentionally overwrites the quick
% diagnostic versions above without changing any metric tables.
make_clean_analysis_figures(T, OUT);

%% README refresh
write_stage_readme(OUT, CFG, S_overall, S_geom, S_freq, S_field, S_worst);
write_docs_readme(root_dir, CFG);

%% Console summary
fprintf('\nStage A analysis complete.\n');
fprintf('Patch/model rows: %d\n', height(T));
print_summary_block(S_overall, 'Global summary');
print_summary_block(S_geom, 'By true SWS / geometry');
print_summary_block(S_field, 'By field regime');
print_summary_block(S_freq, 'By frequency');
X = S_worst(string(S_worst.model_name)=="q_spectrum_plus_composition",:);
if ~isempty(X)
    fprintf('Worst q_spectrum_plus_composition condition: %s | MAPE %.2f%% | bias %.2f%%\n', ...
        X.condition_id(1), X.MAPE_pct(1), X.signed_bias_pct(1));
end
fprintf('Figures: %s\nREADME: %s\n', OUT.analysis_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Local functions
function CFG = load_stage_config(root_dir)
configPath = fullfile(root_dir, 'configs', 'eikonal_validation', 'stage_A_homogeneous_transfer.json');
assert(exist(configPath, 'file') == 2, 'Missing Stage A config: %s', configPath);
CFG = jsondecode(fileread(configPath));
CFG.ConfigPath = string(configPath);
CFG.OutputRoot = string(CFG.OutputRoot);
CFG.FrequenciesHz = double(CFG.FrequenciesHz(:)).';
CFG.FieldRegimes = string(CFG.FieldRegimes(:)).';
CFG.MList = double(CFG.MList(:)).';
CFG.RealismLevels = string(CFG.RealismLevels(:)).';
CFG.REQ.Nbins = string(CFG.REQ.Nbins);
CFG.REQ.EdgeMode = string(CFG.REQ.EdgeMode);
CFG.Simulation.AmplitudeNormalization = string(CFG.Simulation.AmplitudeNormalization);
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
stringVars = ["condition_id","geometry","field_regime","realism_level","model_name","source_layout_id"];
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
        'VariableNames', {'N_valid_patches','MAPE_pct','median_abs_pct_error', ...
        'MAE_SWS','RMSE_SWS','signed_bias_pct','median_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','mean_predicted_SWS','std_predicted_SWS', ...
        'mean_true_SWS','mean_q_pred','std_q_pred'});
end
Tsum = [groups vertcat(rows{:})];
end

function make_summary_figures(T, S_overall, S_geom, S_freq, S_field, S_geom_freq, S_condition, OUT)
models = unique(string(T.model_name), 'stable');
compName = "q_spectrum_plus_composition";

% 1. Global MAPE by model.
fig = newfig([12 8]);
bar(categorical(S_overall.model_name), S_overall.MAPE_pct);
ylabel('MAPE (%)'); xlabel('Model'); title('Global MAPE by model'); grid on;
xtickangle(20);
savefigs(fig, OUT.analysis_dir, 'stage_A_global_mape_by_model');

% 2. Bias by true SWS / geometry.
fig = newfig([14 8]); hold on;
for mi = 1:numel(models)
    X = sortrows(S_geom(string(S_geom.model_name)==models(mi),:), 'mean_true_SWS');
    plot(X.mean_true_SWS, X.signed_bias_pct, '-o', 'LineWidth',1.4, 'DisplayName', pretty_model(models(mi)));
end
yline(0,'k-','Zero bias'); xlabel('True SWS (m/s)'); ylabel('Signed bias (%)');
title('Bias by true homogeneous SWS'); legend('Location','bestoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_bias_by_true_sws');

% 3. MAPE vs frequency.
fig = newfig([14 8]); hold on;
for mi = 1:numel(models)
    X = sortrows(S_freq(string(S_freq.model_name)==models(mi),:), 'f0');
    plot(X.f0, X.MAPE_pct, '-o', 'LineWidth',1.4, 'DisplayName', pretty_model(models(mi)));
end
xlabel('Frequency (Hz)'); ylabel('MAPE (%)'); title('MAPE vs frequency');
legend('Location','bestoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_mape_vs_frequency');

% 4. MAPE by field regime.
fig = newfig([16 8]);
X = S_field(:, {'model_name','field_regime','MAPE_pct'});
g = categorical(X.field_regime, unique(X.field_regime,'stable'));
fieldCats = categories(g);
barData = nan(numel(fieldCats), numel(models));
for ri = 1:numel(fieldCats)
    for mi = 1:numel(models)
        idx = string(X.field_regime)==string(fieldCats{ri}) & string(X.model_name)==models(mi);
        if any(idx), barData(ri,mi) = X.MAPE_pct(idx); end
    end
end
bar(categorical(fieldCats), barData); ylabel('MAPE (%)'); xlabel('Field regime');
title('MAPE by field regime'); legend(arrayfun(@pretty_model, models), 'Location','bestoutside'); grid on; xtickangle(20);
savefigs(fig, OUT.analysis_dir, 'stage_A_mape_by_field_regime');

% 5. Geometry x frequency heatmap for composition.
X = S_geom_freq(string(S_geom_freq.model_name)==compName,:);
fig = newfig([12 8]);
plot_heatmap(X, 'geometry', 'f0', 'MAPE_pct', 'Geometry x frequency MAPE', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_A_heatmap_geometry_frequency_composition');

% 6. Field regime x frequency heatmap for composition.
X = summarize_metrics(T(string(T.model_name)==compName,:), ["field_regime","f0"]);
fig = newfig([13 8]);
plot_heatmap(X, 'field_regime', 'f0', 'MAPE_pct', 'Field regime x frequency MAPE', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_A_heatmap_field_frequency_composition');

% 7. q_pred vs frequency grouped by geometry and field regime.
fig = newfig([18 10]);
X = T(string(T.model_name)==compName,:);
[G, Gtab] = findgroups(X(:, {'geometry','field_regime','f0'}));
Q = splitapply(@(q) mean(q,'omitnan'), X.q_pred, G);
Qtab = [Gtab table(Q, 'VariableNames', {'mean_q_pred'})];
geomList = unique(Qtab.geometry, 'stable');
tiledlayout(fig, 1, numel(geomList), 'TileSpacing','compact','Padding','compact');
for gi = 1:numel(geomList)
    ax = nexttile; hold(ax,'on');
    Y = Qtab(Qtab.geometry==geomList(gi),:);
    regimes = unique(Y.field_regime, 'stable');
    for ri = 1:numel(regimes)
        Z = sortrows(Y(Y.field_regime==regimes(ri),:), 'f0');
        plot(ax, Z.f0, Z.mean_q_pred, '-o', 'DisplayName', pretty_regime(regimes(ri)), 'LineWidth',1.2);
    end
    title(ax, pretty_geometry(geomList(gi))); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Mean predicted q'); grid(ax,'on');
    if gi == numel(geomList), legend(ax,'Location','bestoutside'); end
end
sgtitle(fig, 'Mean predicted REQ quantile q vs frequency');
savefigs(fig, OUT.analysis_dir, 'stage_A_q_pred_vs_frequency_by_geometry_field');

% 8. Predicted vs true SWS at condition-mean level.
fig = newfig([11 9]); hold on;
X = summarize_metrics(T, ["model_name","condition_id","geometry","true_sws","f0","field_regime"]);
for mi = 1:numel(models)
    Y = X(string(X.model_name)==models(mi),:);
    scatter(Y.mean_true_SWS, Y.mean_predicted_SWS, 38, 'filled', 'DisplayName', pretty_model(models(mi)), 'MarkerFaceAlpha',0.75);
end
lims = [min([X.mean_true_SWS; X.mean_predicted_SWS])-0.1, max([X.mean_true_SWS; X.mean_predicted_SWS])+0.1];
plot(lims, lims, 'k--', 'DisplayName','Ideal'); xlim(lims); ylim(lims); axis square;
xlabel('True SWS (m/s)'); ylabel('Predicted SWS, condition mean (m/s)');
title('Predicted vs true homogeneous SWS'); legend('Location','bestoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_predicted_vs_true_sws_condition_mean');
end

function make_additional_diagnostics(T, S_condition, OUT)
models = unique(string(T.model_name), 'stable');
compName = "q_spectrum_plus_composition";

% Bias vs frequency.
S_freq = summarize_metrics(T, ["model_name","f0"]);
fig = newfig([14 8]); hold on;
for mi = 1:numel(models)
    X = sortrows(S_freq(string(S_freq.model_name)==models(mi),:), 'f0');
    plot(X.f0, X.signed_bias_pct, '-o', 'LineWidth',1.4, 'DisplayName', pretty_model(models(mi)));
end
yline(0,'k-'); xlabel('Frequency (Hz)'); ylabel('Signed bias (%)');
title('Signed bias vs frequency'); legend('Location','bestoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_bias_vs_frequency');

% Condition ranking for composition.
X = sortrows(S_condition(string(S_condition.model_name)==compName,:), 'MAPE_pct', 'descend');
if ~isempty(X)
    n = min(20, height(X)); X = X(1:n,:);
    fig = newfig([18 8]);
    labels = categorical(short_condition_label(X));
    labels = reordercats(labels, cellstr(short_condition_label(X)));
    bar(labels, X.MAPE_pct); ylabel('MAPE (%)'); xlabel('Condition');
    title('Worst Stage A conditions for q spectrum + composition'); grid on; xtickangle(45);
    savefigs(fig, OUT.analysis_dir, 'stage_A_worst_conditions_composition');
end

% Error distributions.
fig = newfig([14 8]); hold on;
edges = linspace(0, max(1, prctile(T.abs_error_percent,99)), 36);
for mi = 1:numel(models)
    X = T(string(T.model_name)==models(mi),:);
    histogram(X.abs_error_percent, edges, 'Normalization','probability', 'DisplayName', pretty_model(models(mi)), 'FaceAlpha',0.5);
end
xlabel('Absolute SWS error (%)'); ylabel('Patch fraction'); title('Patch-level error distribution');
legend('Location','bestoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_abs_error_distribution');

% Composition diagnostics in a homogeneous-only transfer test.
if all(ismember(["predicted_patch_purity","p_mixed","p_strong_mixed"], string(T.Properties.VariableNames)))
    X = T(string(T.model_name)==compName,:);
    S = summarize_extra(X, ["geometry","field_regime"]);
    writetable(S, fullfile(OUT.table_dir, 'stage_A_composition_diagnostics.csv'));
    fig = newfig([16 8]);
    tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile; grouped_bar(ax, S, 'geometry', 'mean_predicted_patch_purity', 'Mean predicted patch purity');
    ax = nexttile; grouped_bar(ax, S, 'geometry', 'mean_p_mixed', 'Mean predicted mixed probability');
    sgtitle(fig, 'Composition auxiliary diagnostics in homogeneous Eikonal fields');
    savefigs(fig, OUT.analysis_dir, 'stage_A_composition_auxiliary_diagnostics');
end

% Error vs local amplitude, useful for detecting hidden near-zero-field issues.
if ismember('local_amplitude', string(T.Properties.VariableNames))
    fig = newfig([12 8]);
    X = T(string(T.model_name)==compName,:);
    scatter(X.local_amplitude, X.abs_error_percent, 8, double(X.f0), 'filled', 'MarkerFaceAlpha',0.35);
    xlabel('Local displacement amplitude (a.u.)'); ylabel('Absolute SWS error (%)');
    title('Error vs local amplitude'); cb = colorbar; ylabel(cb,'Frequency (Hz)'); grid on;
    savefigs(fig, OUT.analysis_dir, 'stage_A_error_vs_local_amplitude');
end
end

function S = summarize_extra(T, groupVars)
if isstring(groupVars), groupVars = cellstr(groupVars); end
[G, groups] = findgroups(T(:, groupVars));
meanPur = splitapply(@(x) mean(x,'omitnan'), T.predicted_patch_purity, G);
meanMix = splitapply(@(x) mean(x,'omitnan'), T.p_mixed, G);
meanStrong = splitapply(@(x) mean(x,'omitnan'), T.p_strong_mixed, G);
S = [groups table(meanPur, meanMix, meanStrong, 'VariableNames', {'mean_predicted_patch_purity','mean_p_mixed','mean_p_strong_mixed'})];
end

function grouped_bar(ax, S, xVar, yVar, ttl)
X = S(:, {xVar,'field_regime',yVar});
xs = unique(string(X.(xVar)), 'stable');
regs = unique(string(X.field_regime), 'stable');
Y = nan(numel(xs), numel(regs));
for i = 1:numel(xs)
    for j = 1:numel(regs)
        idx = string(X.(xVar))==xs(i) & string(X.field_regime)==regs(j);
        if any(idx), Y(i,j)=X.(yVar)(idx); end
    end
end
bar(ax, categorical(xs), Y); grid(ax,'on'); title(ax, ttl); xlabel(ax, strrep(xVar,'_',' ')); ylabel(ax, ttl);
legend(ax, arrayfun(@pretty_regime, regs), 'Location','bestoutside'); xtickangle(ax,20);
end

function labels = short_condition_label(T)
labels = string(T.geometry) + " f" + string(T.f0) + " " + erase(string(T.field_regime), "diffuse_like_8src_");
end

function plot_heatmap(T, rowVar, colVar, valueVar, ttl, cbLabel)
rows = unique(string(T.(rowVar)), 'stable');
cols = unique(T.(colVar), 'stable');
Z = nan(numel(rows), numel(cols));
for i = 1:numel(rows)
    for j = 1:numel(cols)
        idx = string(T.(rowVar))==rows(i) & T.(colVar)==cols(j);
        if any(idx), Z(i,j) = T.(valueVar)(find(idx,1)); end
    end
end
imagesc(Z); axis image; colormap(parula); cb = colorbar; ylabel(cb, cbLabel);
set(gca,'XTick',1:numel(cols),'XTickLabel',string(cols), 'YTick',1:numel(rows),'YTickLabel',arrayfun(@pretty_generic, rows));
xlabel(strrep(colVar,'_',' ')); ylabel(strrep(rowVar,'_',' ')); title(ttl, 'FontWeight','normal');
end

function fig = newfig(sizeCm)
fig = figure('Color','w','Units','centimeters','Position',[2 2 sizeCm]);
end

function savefigs(fig, outDir, baseName)
ensure_dir(outDir);
png = fullfile(outDir, baseName + ".png");
figp = fullfile(outDir, baseName + ".fig");
try
    exportgraphics(fig, png, 'Resolution', 220);
catch
    saveas(fig, png);
end
try
    savefig(fig, figp);
catch
end
close(fig);
end


function make_clean_analysis_figures(T, OUT)
% Publication/debug friendly figures for Stage A. These plots use readable
% labels everywhere so MATLAB does not render raw underscores as subscripts.
models = unique(string(T.model_name), 'stable');
compName = "q_spectrum_plus_composition";
S_overall = summarize_metrics(T, "model_name");
S_geom = summarize_metrics(T, ["model_name","geometry"]);
S_freq = summarize_metrics(T, ["model_name","f0"]);
S_field = summarize_metrics(T, ["model_name","field_regime"]);
S_geom_freq = summarize_metrics(T, ["model_name","geometry","f0"]);
S_condition = summarize_metrics(T, ["model_name","condition_id","geometry","true_sws","f0","field_regime"]);

% 1. Global MAPE by model.
fig = cleanfig([16 10]);
bar(categorical(arrayfun(@pretty_model, string(S_overall.model_name))), S_overall.MAPE_pct, 0.65);
ylabel('MAPE (%)'); xlabel('Model'); title('Global SWS Error by Model', 'FontSize', 11);
grid on; set(gca,'TickLabelInterpreter','none'); ylim padded;
savefigs(fig, OUT.analysis_dir, 'stage_A_global_mape_by_model');

% 2. Bias by true SWS.
fig = cleanfig([17 10]); hold on;
for mi = 1:numel(models)
    X = sortrows(S_geom(string(S_geom.model_name)==models(mi),:), 'mean_true_SWS');
    plot(X.mean_true_SWS, X.signed_bias_pct, '-o', 'LineWidth',1.7, 'MarkerSize',5, 'DisplayName', pretty_model(models(mi)));
end
yline(0,'k-','Zero bias', 'LabelHorizontalAlignment','left');
xlabel('True SWS (m/s)'); ylabel('Signed bias (%)'); title('Signed Bias by Homogeneous SWS', 'FontSize', 11);
legend('Location','eastoutside'); grid on; xlim([1.8 4.2]);
savefigs(fig, OUT.analysis_dir, 'stage_A_bias_by_true_sws');

% 3. MAPE vs frequency.
fig = cleanfig([17 10]); hold on;
for mi = 1:numel(models)
    X = sortrows(S_freq(string(S_freq.model_name)==models(mi),:), 'f0');
    plot(X.f0, X.MAPE_pct, '-o', 'LineWidth',1.7, 'MarkerSize',5, 'DisplayName', pretty_model(models(mi)));
end
xlabel('Frequency (Hz)'); ylabel('MAPE (%)'); title('SWS Error vs Frequency', 'FontSize', 11);
legend('Location','eastoutside'); grid on; xticks(unique(T.f0));
savefigs(fig, OUT.analysis_dir, 'stage_A_mape_vs_frequency');

% 4. MAPE by field regime.
fig = cleanfig([20 10]);
regs = unique(string(S_field.field_regime), 'stable');
Y = nan(numel(regs), numel(models));
for ri = 1:numel(regs)
    for mi = 1:numel(models)
        idx = string(S_field.field_regime)==regs(ri) & string(S_field.model_name)==models(mi);
        if any(idx), Y(ri,mi) = S_field.MAPE_pct(idx); end
    end
end
bar(categorical(arrayfun(@pretty_regime, regs)), Y, 0.72);
ylabel('MAPE (%)'); xlabel('Field regime'); title('SWS Error by Field Regime', 'FontSize', 11);
legend(arrayfun(@pretty_model, models), 'Location','eastoutside'); grid on; xtickangle(12); set(gca,'TickLabelInterpreter','none');
savefigs(fig, OUT.analysis_dir, 'stage_A_mape_by_field_regime');

% 5. Geometry x frequency heatmap for composition.
X = S_geom_freq(string(S_geom_freq.model_name)==compName,:);
fig = cleanfig([16 10]);
clean_heatmap(X, 'geometry', 'f0', 'MAPE_pct', 'q spectrum + composition: MAPE by SWS and Frequency', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_A_heatmap_geometry_frequency_composition');

% 6. Field regime x frequency heatmap for composition.
X = summarize_metrics(T(string(T.model_name)==compName,:), ["field_regime","f0"]);
fig = cleanfig([18 10]);
clean_heatmap(X, 'field_regime', 'f0', 'MAPE_pct', 'q spectrum + composition: MAPE by Field and Frequency', 'MAPE (%)');
savefigs(fig, OUT.analysis_dir, 'stage_A_heatmap_field_frequency_composition');

% 7. Mean q vs frequency by geometry and field regime.
fig = cleanfig([22 10]);
X = T(string(T.model_name)==compName,:);
[G, Gtab] = findgroups(X(:, {'geometry','field_regime','f0'}));
Q = splitapply(@(q) mean(q,'omitnan'), X.q_pred, G);
Qtab = [Gtab table(Q, 'VariableNames', {'mean_q_pred'})];
geomList = unique(Qtab.geometry, 'stable');
tl = tiledlayout(fig, 1, numel(geomList), 'TileSpacing','compact','Padding','compact');
for gi = 1:numel(geomList)
    ax = nexttile; hold(ax,'on');
    Yg = Qtab(Qtab.geometry==geomList(gi),:);
    regimes = unique(Yg.field_regime, 'stable');
    for ri = 1:numel(regimes)
        Z = sortrows(Yg(Yg.field_regime==regimes(ri),:), 'f0');
        plot(ax, Z.f0, Z.mean_q_pred, '-o', 'DisplayName', pretty_regime(regimes(ri)), 'LineWidth',1.4, 'MarkerSize',4);
    end
    title(ax, pretty_geometry(geomList(gi)), 'FontSize',10); xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Mean REQ quantile q'); grid(ax,'on'); xticks(ax, unique(T.f0));
    if gi == numel(geomList), legend(ax,'Location','eastoutside'); end
end
title(tl, 'Mean Predicted REQ Quantile vs Frequency', 'FontSize',12, 'FontWeight','bold');
savefigs(fig, OUT.analysis_dir, 'stage_A_q_pred_vs_frequency_by_geometry_field');

% 8. Predicted vs true SWS at condition-mean level.
fig = cleanfig([14 12]); hold on;
X = summarize_metrics(T, ["model_name","condition_id","geometry","true_sws","f0","field_regime"]);
for mi = 1:numel(models)
    Yc = X(string(X.model_name)==models(mi),:);
    scatter(Yc.mean_true_SWS, Yc.mean_predicted_SWS, 46, 'filled', 'DisplayName', pretty_model(models(mi)), 'MarkerFaceAlpha',0.70);
end
lims = [min([X.mean_true_SWS; X.mean_predicted_SWS])-0.08, max([X.mean_true_SWS; X.mean_predicted_SWS])+0.08];
plot(lims, lims, 'k--', 'DisplayName','Ideal'); xlim(lims); ylim(lims); axis square;
xlabel('True SWS (m/s)'); ylabel('Predicted SWS, condition mean (m/s)');
title('Condition-Mean Predicted vs True SWS', 'FontSize', 11); legend('Location','eastoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_predicted_vs_true_sws_condition_mean');

% 9. Bias vs frequency.
fig = cleanfig([17 10]); hold on;
for mi = 1:numel(models)
    Xf = sortrows(S_freq(string(S_freq.model_name)==models(mi),:), 'f0');
    plot(Xf.f0, Xf.signed_bias_pct, '-o', 'LineWidth',1.7, 'MarkerSize',5, 'DisplayName', pretty_model(models(mi)));
end
yline(0,'k-'); xlabel('Frequency (Hz)'); ylabel('Signed bias (%)'); title('Signed Bias vs Frequency', 'FontSize', 11);
legend('Location','eastoutside'); grid on; xticks(unique(T.f0));
savefigs(fig, OUT.analysis_dir, 'stage_A_bias_vs_frequency');

% 10. Worst conditions for composition.
Xw = sortrows(S_condition(string(S_condition.model_name)==compName,:), 'MAPE_pct', 'descend');
if ~isempty(Xw)
    Xw = Xw(1:min(18,height(Xw)),:);
    fig = cleanfig([24 11]);
    labels = categorical(clean_condition_label(Xw));
    labels = reordercats(labels, cellstr(clean_condition_label(Xw)));
    bar(labels, Xw.MAPE_pct, 0.75); ylabel('MAPE (%)'); xlabel('Condition');
    title('Worst Homogeneous Eikonal Conditions: q Spectrum + Composition', 'FontSize', 11);
    grid on; xtickangle(35); set(gca,'TickLabelInterpreter','none');
    savefigs(fig, OUT.analysis_dir, 'stage_A_worst_conditions_composition');
end

% 11. Error distributions.
fig = cleanfig([16 10]); hold on;
edges = linspace(0, max(0.5, prctile(T.abs_error_percent,99.5)), 42);
for mi = 1:numel(models)
    Xm = T(string(T.model_name)==models(mi),:);
    histogram(Xm.abs_error_percent, edges, 'Normalization','probability', 'DisplayName', pretty_model(models(mi)), 'FaceAlpha',0.50);
end
xlabel('Absolute SWS error (%)'); ylabel('Patch fraction'); title('Patch-Level Absolute Error Distribution', 'FontSize', 11);
legend('Location','eastoutside'); grid on;
savefigs(fig, OUT.analysis_dir, 'stage_A_abs_error_distribution');

% 12. Composition auxiliary diagnostics.
if all(ismember(["predicted_patch_purity","p_mixed","p_strong_mixed"], string(T.Properties.VariableNames)))
    Xc = T(string(T.model_name)==compName,:);
    S = summarize_extra(Xc, ["geometry","field_regime"]);
    fig = cleanfig([20 10]);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile; clean_grouped_bar(ax, S, 'geometry', 'mean_predicted_patch_purity', 'Predicted patch purity', 'Mean value');
    ax = nexttile; clean_grouped_bar(ax, S, 'geometry', 'mean_p_mixed', 'Predicted mixed probability', 'Mean probability');
    title(tl, 'Composition Auxiliary Outputs in Homogeneous Fields', 'FontSize',12, 'FontWeight','bold');
    savefigs(fig, OUT.analysis_dir, 'stage_A_composition_auxiliary_diagnostics');
end

% 13. Error vs local amplitude.
if ismember('local_amplitude', string(T.Properties.VariableNames))
    fig = cleanfig([15 11]);
    Xa = T(string(T.model_name)==compName,:);
    scatter(Xa.local_amplitude, Xa.abs_error_percent, 9, double(Xa.f0), 'filled', 'MarkerFaceAlpha',0.30);
    xlabel('Local displacement amplitude (a.u.)'); ylabel('Absolute SWS error (%)');
    title('Error vs Local Field Amplitude: q Spectrum + Composition', 'FontSize', 11);
    cb = colorbar; ylabel(cb,'Frequency (Hz)'); grid on;
    savefigs(fig, OUT.analysis_dir, 'stage_A_error_vs_local_amplitude');
end

% 14. Extra: hard-speed bias by field regime.
fig = cleanfig([18 10]);
X = summarize_metrics(T(string(T.model_name)==compName,:), ["geometry","field_regime"]);
clean_heatmap(X, 'geometry', 'field_regime', 'signed_bias_pct', 'q spectrum + composition: Signed Bias by SWS and Field', 'Signed bias (%)');
savefigs(fig, OUT.analysis_dir, 'stage_A_bias_heatmap_geometry_field_composition');
end

function fig = cleanfig(sizeCm)
fig = figure('Color','w','Units','centimeters','Position',[2 2 sizeCm], 'Renderer','painters');
set(fig, 'PaperPositionMode','auto');
end

function clean_heatmap(T, rowVar, colVar, valueVar, ttl, cbLabel)
rows = unique(string(T.(rowVar)), 'stable');
colsRaw = unique(T.(colVar), 'stable');
Z = nan(numel(rows), numel(colsRaw));
for i = 1:numel(rows)
    for j = 1:numel(colsRaw)
        idx = string(T.(rowVar))==rows(i) & T.(colVar)==colsRaw(j);
        if any(idx), Z(i,j) = T.(valueVar)(find(idx,1)); end
    end
end
imagesc(Z); axis normal; colormap(parula); cb = colorbar; ylabel(cb, cbLabel, 'FontSize',9);
if isnumeric(colsRaw)
    colLabels = string(colsRaw);
else
    colLabels = arrayfun(@pretty_generic, string(colsRaw));
end
rowLabels = arrayfun(@pretty_generic, rows);
set(gca,'XTick',1:numel(colsRaw),'XTickLabel',colLabels, 'YTick',1:numel(rows),'YTickLabel',rowLabels, 'TickLabelInterpreter','none', 'FontSize',9);
xlabel(clean_axis_label(colVar)); ylabel(clean_axis_label(rowVar)); title(ttl, 'FontSize',11, 'FontWeight','bold', 'Interpreter','none');
for i = 1:numel(rows)
    for j = 1:numel(colsRaw)
        if isfinite(Z(i,j))
            text(j, i, sprintf('%.2f', Z(i,j)), 'HorizontalAlignment','center', 'FontSize',8, 'Color','w', 'FontWeight','bold');
        end
    end
end
end

function clean_grouped_bar(ax, S, xVar, yVar, ttl, ylab)
X = S(:, {xVar,'field_regime',yVar});
xs = unique(string(X.(xVar)), 'stable');
regs = unique(string(X.field_regime), 'stable');
Y = nan(numel(xs), numel(regs));
for i = 1:numel(xs)
    for j = 1:numel(regs)
        idx = string(X.(xVar))==xs(i) & string(X.field_regime)==regs(j);
        if any(idx), Y(i,j)=X.(yVar)(idx); end
    end
end
bar(ax, categorical(arrayfun(@pretty_generic, xs)), Y, 0.72); grid(ax,'on');
title(ax, ttl, 'FontSize',10); xlabel(ax, clean_axis_label(xVar)); ylabel(ax, ylab);
legend(ax, arrayfun(@pretty_regime, regs), 'Location','eastoutside'); xtickangle(ax,15); set(ax,'TickLabelInterpreter','none');
end

function labels = clean_condition_label(T)
labels = arrayfun(@pretty_geometry, string(T.geometry)) + " | " + string(T.f0) + " Hz | " + arrayfun(@short_regime, string(T.field_regime));
end

function name = short_regime(name)
name = string(name);
switch name
    case "single_source_lateral"
        name = "single";
    case "diffuse_like_8src_layout1"
        name = "diffuse L1";
    case "diffuse_like_8src_layout2"
        name = "diffuse L2";
    otherwise
        name = strrep(name, '_', ' ');
end
end

function label = clean_axis_label(varName)
varName = string(varName);
switch varName
    case "geometry"
        label = "Homogeneous SWS";
    case "f0"
        label = "Frequency (Hz)";
    case "field_regime"
        label = "Field regime";
    otherwise
        label = strrep(varName, '_', ' ');
end
end

function write_stage_readme(OUT, CFG, T_overall, T_geom, T_freq, T_field, T_worst)
path = fullfile(OUT.root_dir, 'README_results.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Stage A: Homogeneous Eikonal Transfer Validation\n\n');
fprintf(fid, '## Objective\n\n');
fprintf(fid, ['Validate zero-shot transfer of frozen `baseline_minimal_v1` q/SWS models to clean homogeneous Eikonal simulations. ', ...
    'This stage intentionally removes interfaces, material mixing, readout noise, attenuation, correction, oracle q, and risk/reliability masks.\n\n']);
fprintf(fid, '## Scientific Questions\n\n');
fprintf(fid, ['1. Does the baseline predict SWS correctly in homogeneous Eikonal media?\n', ...
    '2. Is there a systematic transfer bias from ideal clean training simulations to Eikonal simulations?\n', ...
    '3. Does bias depend on true SWS: 2, 3, or 4 m/s?\n', ...
    '4. Does performance differ between directional and diffuse-like fields?\n', ...
    '5. Does performance depend on frequency?\n', ...
    '6. Does `q_spectrum_plus_composition` remain comparable or better than `q_spectrum_only`?\n\n']);
fprintf(fid, '## Simulation Design\n\n');
fprintf(fid, '- Geometry cases: `homogeneous_cs2`, `homogeneous_cs3`, `homogeneous_cs4`.\n');
fprintf(fid, '- Frequencies: %s Hz.\n', strjoin(string(CFG.FrequenciesHz), ', '));
fprintf(fid, '- Field regimes: `single_source_lateral`, `diffuse_like_8src_layout1`, `diffuse_like_8src_layout2`.\n');
fprintf(fid, '- Diffuse-like fields use deterministic fixed source layouts: exactly 8 sources, 4 in-plane and 4 off-plane, equal source amplitudes, deterministic phases.\n');
fprintf(fid, '- Realism level: clean only. No shear attenuation, no ultrasound readout noise, and no interface modules.\n');
fprintf(fid, '- Measurement plane: `y = 0`; measured displacement component: z; all source motion directions: `[0 0 1]`.\n');
fprintf(fid, '- Amplitude normalization: `%s`; pre/post RMS and normalization scale are saved.\n\n', CFG.Simulation.AmplitudeNormalization);
fprintf(fid, '## REQ And Model Settings\n\n');
fprintf(fid, '- `M = %s`, `cs_guess = %.2f m/s`.\n', mat2str(CFG.MList), CFG.REQ.cs_guess);
fprintf(fid, '- `dx = dz = %.1f mm`, `TargetStepM = %.1f mm`. The 1.0 mm step is chosen for Stage A because the maps are homogeneous and the goal is transfer bias, not edge-detail rendering.\n', CFG.Simulation.dx_m*1e3, CFG.REQ.TargetStepM*1e3);
fprintf(fid, '- `EdgeMode = %s`, `Nbins = %s`, `Nbins_min = %d`, `smooth_sigma = %.1f`, `Gamma = %.1f`, `PadFactor = %.1f`.\n', ...
    CFG.REQ.EdgeMode, CFG.REQ.Nbins, CFG.REQ.Nbins_min, CFG.REQ.smooth_sigma, CFG.REQ.Gamma, CFG.REQ.PadFactor);
fprintf(fid, '- Evaluated frozen models: `q_spectrum_only`, `q_spectrum_plus_composition`.\n\n');
fprintf(fid, '## Metrics\n\n');
fprintf(fid, 'MAPE, signed bias, MAE/RMSE SWS, median absolute percentage error, high-error fractions >10%% and >20%%, mean predicted SWS, mean true SWS, mean q, and valid patch count.\n\n');
fprintf(fid, '## Generated Figures\n\n');
fprintf(fid, '- Summary plots under `figures/analysis_summary/`.\n');
fprintf(fid, '- Representative field and prediction maps under `figures/maps_by_condition/`.\n');
fprintf(fid, '- Source geometry diagnostics under `figures/source_geometry/`.\n');
fprintf(fid, '- Central power spectra for diffuse-like layouts under `figures/central_power_spectra/`.\n\n');
fprintf(fid, '## How To Run\n\n');
fprintf(fid, 'Validation-only:\n\n```bash\ncd /Users/sara/local/adaptive_req_local\n/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''validate_only''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, 'Quick:\n\n```bash\ncd /Users/sara/local/adaptive_req_local\n/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''quick''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, 'Full runner and analysis:\n\n```bash\ncd /Users/sara/local/adaptive_req_local\n/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''full''); run(''experiments/runners/eikonal_validation/run_stage_A_homogeneous_transfer.m'')"\n/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv(''ADAPTIVE_REQ_EIKONAL_STAGE_A_MODE'',''full''); run(''experiments/analysis/eikonal_validation/analyze_stage_A_homogeneous_transfer.m'')"\n```\n\n');
fprintf(fid, '## Approximate Runtime\n\n');
fprintf(fid, 'Runtime is printed by the runner. Cached conditions are reused unless `ADAPTIVE_REQ_EIKONAL_STAGE_A_FORCE_REBUILD=true`. Full Stage A has 36 simulation conditions, each evaluated by two frozen models.\n\n');
fprintf(fid, '## Results Summary\n\n');
fprintf(fid, '### Overall\n\n'); write_markdown_table(fid, T_overall);
fprintf(fid, '\n### By Geometry\n\n'); write_markdown_table(fid, T_geom);
fprintf(fid, '\n### By Frequency\n\n'); write_markdown_table(fid, T_freq);
fprintf(fid, '\n### By Field Regime\n\n'); write_markdown_table(fid, T_field);
fprintf(fid, '\n### Worst Conditions\n\n'); write_markdown_table(fid, T_worst(1:min(10,height(T_worst)),:));
fprintf(fid, '\n## Brief Interpretation\n\n');
fprintf(fid, '%s', interpretation_text(T_overall, T_geom, T_field, T_freq));
fclose(fid);
end

function txt = interpretation_text(T_overall, T_geom, T_field, T_freq)
txt = "";
comp = T_overall(string(T_overall.model_name)=="q_spectrum_plus_composition",:);
only = T_overall(string(T_overall.model_name)=="q_spectrum_only",:);
if isempty(comp)
    txt = "No `q_spectrum_plus_composition` rows were found." + newline;
    return;
end
txt = txt + sprintf('`q_spectrum_plus_composition` global MAPE is %.2f%% with signed bias %.2f%%. ', comp.MAPE_pct(1), comp.signed_bias_pct(1));
if ~isempty(only)
    delta = comp.MAPE_pct(1) - only.MAPE_pct(1);
    if delta <= 0
        txt = txt + sprintf('It is comparable or better than `q_spectrum_only` by %.2f MAPE points. ', abs(delta));
    else
        txt = txt + sprintf('It is worse than `q_spectrum_only` by %.2f MAPE points in this run. ', delta);
    end
end
% Simple warnings that help decide whether to pause later stages.
if abs(comp.signed_bias_pct(1)) > 5 || comp.MAPE_pct(1) > 8
    txt = txt + 'This indicates a non-trivial homogeneous Eikonal transfer issue; later interface/noise stages should be interpreted cautiously. ';
else
    txt = txt + 'This supports moving to later interface/noise stages, because clean homogeneous transfer does not appear to be the dominant failure mode. ';
end
if ~isempty(T_geom)
    G = T_geom(string(T_geom.model_name)=="q_spectrum_plus_composition",:);
    if ~isempty(G)
        [~,ix] = max(G.MAPE_pct);
        txt = txt + sprintf('The hardest homogeneous SWS level in this run is `%s` with MAPE %.2f%%. ', G.geometry(ix), G.MAPE_pct(ix));
    end
end
txt = txt + newline;
end

function write_docs_readme(root_dir, CFG)
docDir = fullfile(root_dir, 'docs', 'eikonal_validation');
ensure_dir(docDir);
path = fullfile(docDir, 'stage_A_homogeneous_transfer.md');
fid = fopen(path, 'w'); assert(fid > 0);
fprintf(fid, '# Stage A: Homogeneous Eikonal Transfer Validation\n\n');
fprintf(fid, 'Stage A evaluates frozen `baseline_minimal_v1` models on clean homogeneous Eikonal simulations. It is the first gate before adding interfaces, attenuation, readout noise, or confidence masks.\n\n');
fprintf(fid, 'Configuration: `%s`\n\n', CFG.ConfigPath);
fprintf(fid, 'Primary outputs: `outputs/eikonal_validation/stage_A_homogeneous_transfer/`.\n');
fclose(fid);
end

function write_markdown_table(fid, T)
if isempty(T), fprintf(fid, '_No rows._\n'); return; end
names = string(T.Properties.VariableNames);
fprintf(fid, '| %s |\n', strjoin(names, ' | '));
fprintf(fid, '| %s |\n', strjoin(repmat("---", size(names)), ' | '));
maxRows = min(height(T), 80);
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
    ["model_name","geometry","f0","field_regime","N_valid_patches","MAPE_pct","signed_bias_pct","high_error20_pct"], 'stable');
disp(T(:, cellstr(showVars)));
end

function name = pretty_model(name)
name = string(name);
switch name
    case "q_spectrum_only"
        name = "q spectrum only";
    case "q_spectrum_plus_composition"
        name = "q spectrum + composition";
    otherwise
        name = strrep(name, '_', ' ');
end
end

function name = pretty_regime(name)
name = string(name);
switch name
    case "single_source_lateral"
        name = "single source lateral";
    case "diffuse_like_8src_layout1"
        name = "diffuse-like 8 src layout 1";
    case "diffuse_like_8src_layout2"
        name = "diffuse-like 8 src layout 2";
    otherwise
        name = strrep(name, '_', ' ');
end
end

function name = pretty_geometry(name)
name = string(name);
switch name
    case "homogeneous_cs2"
        name = "homogeneous 2 m/s";
    case "homogeneous_cs3"
        name = "homogeneous 3 m/s";
    case "homogeneous_cs4"
        name = "homogeneous 4 m/s";
    otherwise
        name = strrep(name, '_', ' ');
end
end

function name = pretty_generic(name)
name = string(name);
if startsWith(name, "homogeneous")
    name = pretty_geometry(name);
elseif contains(name, "source") || contains(name, "diffuse")
    name = pretty_regime(name);
else
    name = strrep(name, '_', ' ');
end
end

function ensure_dir(d)
if exist(d, 'dir') ~= 7, mkdir(d); end
end
