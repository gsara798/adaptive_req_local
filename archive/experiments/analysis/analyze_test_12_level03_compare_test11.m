%% analyze_test_12_level03_compare_test11.m
% Test 12 Level 03: compare Test 12 against Test 11.
%
% This analysis does not regenerate either dataset. It compares global
% grouped metrics and, when prediction tables are available, matched subsets.

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

%% Locate inputs

[T12_feat, MC12, PATHS12] = adaptive_req.analysis.load_mc_results( ...
    'test_12_cs_guess_window_sweep', 'RootDir', root_dir, 'Verbose', true);
[T11_feat, MC11, PATHS11] = adaptive_req.analysis.load_mc_results( ...
    'test_11_global_req_features', 'RootDir', root_dir, 'Verbose', true);

T12_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T12_feat);
T11_feat = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T11_feat);
T12_feat = adaptive_req.analysis.Test12Analysis.addBins(T12_feat);
T11_feat = adaptive_req.analysis.Test12Analysis.addBins(T11_feat);

OUT = adaptive_req.analysis.Test12Analysis.makeOutputDirs( ...
    PATHS12, 'level_03_compare_test11');

test11_l04 = fullfile(PATHS11.analysis_dir, ...
    'level_04_grouped_generalization', 'tables');
test12_l01 = fullfile(PATHS12.analysis_dir, ...
    'level_01_model_comparison', 'tables');
test12_l02 = fullfile(PATHS12.analysis_dir, ...
    'level_02_grouped_generalization', 'tables');

require_file(fullfile(test11_l04, 'level11_level04_sws_metrics_by_model.csv'));
require_file(fullfile(test12_l02, 'level12_level02_sws_metrics_by_feature_set.csv'));
require_file(fullfile(test12_l01, 'level12_level01_delta_mape_vs_NoCsGuess.csv'));
require_file(fullfile(test12_l02, 'level12_level02_delta_mape_vs_NoCsGuess.csv'));

T11_by_model = readtable(fullfile(test11_l04, ...
    'level11_level04_sws_metrics_by_model.csv'), 'TextType', 'string');
T12_by_feature = readtable(fullfile(test12_l02, ...
    'level12_level02_sws_metrics_by_feature_set.csv'), 'TextType', 'string');
T12_delta_l01 = readtable(fullfile(test12_l01, ...
    'level12_level01_delta_mape_vs_NoCsGuess.csv'), 'TextType', 'string');
T12_delta_l02 = readtable(fullfile(test12_l02, ...
    'level12_level02_delta_mape_vs_NoCsGuess.csv'), 'TextType', 'string');

%% A. Global Test 11 vs Test 12 comparison

T_global = build_global_comparison(T11_by_model, T12_by_feature);
writetable(T_global, fullfile(OUT.table_dir, ...
    'level12_level03_test11_vs_test12_global_metrics.csv'));

%% B. Matched/subset comparison

T_matched = table();
matched_filters = strings(0, 1);

pred11_file = fullfile(test11_l04, 'level11_level04_grouped_predictions.csv');
pred12_file = fullfile(test12_l02, 'level12_level02_grouped_predictions.csv');

if exist(pred11_file, 'file') && exist(pred12_file, 'file')
    fprintf('\nReading prediction tables for matched comparison. This can take a while.\n');
    P11 = readtable(pred11_file, 'TextType', 'string');
    P12 = readtable(pred12_file, 'TextType', 'string');

    P11 = join_prediction_metadata(P11, T11_feat);
    P12 = join_prediction_metadata(P12, T12_feat);

    [P11m, P12m, matched_filters] = apply_matched_filters(P11, P12);
    T_matched = summarize_matched(P11m, P12m);
else
    warning('Prediction CSV files were not found. Matched comparison will be TODO.');
end

writetable(T_matched, fullfile(OUT.table_dir, ...
    'level12_level03_test11_vs_test12_matched_metrics.csv'));

%% C. Test 12 cs_guess / M_eff_guess ablation

T_ablation = summarize_ablation(T12_delta_l01, T12_delta_l02);
writetable(T_ablation, fullfile(OUT.table_dir, ...
    'level12_level03_cs_guess_ablation_summary.csv'));

%% D. Outlier comparison

if exist('P11m', 'var') && exist('P12m', 'var')
    T_outliers = build_outlier_comparison(P11m, P12m);
else
    T_outliers = build_outlier_comparison_from_metrics(T11_by_model, T12_by_feature);
end
writetable(T_outliers, fullfile(OUT.table_dir, ...
    'level12_level03_test11_vs_test12_outlier_comparison.csv'));

%% Figures

plot_global_bars(T_global, 'MAPE_pct', ...
    fullfile(OUT.fig_dir, 'level12_level03_test11_vs_test12_mape.png'));
plot_global_bars(T_global, 'HighError_gt20_pct', ...
    fullfile(OUT.fig_dir, 'level12_level03_test11_vs_test12_high_error.png'));

if ~isempty(T_matched)
    plot_global_bars(T_matched, 'MAPE_pct', ...
        fullfile(OUT.fig_dir, 'level12_level03_matched_test11_vs_test12_mape.png'));
end

plot_ablation(T_ablation, 'Delta_MAPE_vs_NoCsGuess', ...
    fullfile(OUT.fig_dir, 'level12_level03_cs_guess_ablation_delta_mape.png'));
plot_ablation(T_ablation, 'Delta_HighError_gt20_vs_NoCsGuess', ...
    fullfile(OUT.fig_dir, 'level12_level03_cs_guess_ablation_delta_high_error.png'));

plot_outlier_bins(T_outliers, 'M_eff_true_diag_bin', ...
    fullfile(OUT.fig_dir, 'level12_level03_outliers_by_M_eff_true_diag_test11_vs_test12.png'));
plot_outlier_bins(T_outliers, 'M_eff_guess_bin', ...
    fullfile(OUT.fig_dir, 'level12_level03_outliers_by_M_eff_guess_test11_vs_test12.png'));

%% Markdown summary

write_markdown_summary(fullfile(root_dir, 'docs', ...
    'test_12_cs_guess_window_sweep_analysis_summary.md'), ...
    T_global, T_matched, T_ablation, matched_filters);

%% Console report

fprintf('\nTest 12 Level 03 comparison complete.\n');
fprintf('Analysis folder:\n%s\n', OUT.analysis_dir);
fprintf('\nGlobal comparison.\n');
disp(T_global);

if ~isempty(T_matched)
    fprintf('\nMatched comparison filters:\n');
    disp(matched_filters);
    fprintf('\nMatched comparison.\n');
    disp(T_matched);
else
    fprintf('\nMatched comparison was not computed because prediction tables are missing.\n');
end

fprintf('\nAblation summary. Negative deltas mean improvement over NoCsGuess.\n');
disp(T_ablation);

%% Local functions

function require_file(file_path)
if ~exist(file_path, 'file')
    error(['Required file not found:\n%s\nRun the prerequisite analysis ', ...
        'before Level 03.'], file_path);
end
end

function T_global = build_global_comparison(T11, T12)

T11 = T11(T11.model_type == "bagged_trees", :);
T12 = T12(T12.model_type == "bagged_trees" & ...
    T12.model_role == "operational", :);

models = ["LocalOnly"; "GlobalOnly"; "HybridLocalGlobal"];
T_global = table();

for i = 1:numel(models)
    m = models(i);

    A = T11(T11.model_name == m, :);
    if ~isempty(A)
        row = average_metric_row(A);
        row.dataset = "Test11";
        row.comparison_type = "global";
        row.feature_set = "Test11Baseline";
        T_global = adaptive_req.analysis.Test12Analysis.concatTables(T_global, row);
    end

    B = T12(T12.model_name == m, :);
    if ~isempty(B)
        Bavg = average_by_feature(B);
        Bavg = sortrows(Bavg, 'MAPE_pct', 'ascend');
        row = Bavg(1, :);
        row.dataset = "Test12";
        row.comparison_type = "global";
        T_global = adaptive_req.analysis.Test12Analysis.concatTables(T_global, row);
    end
end

T_global = movevars(T_global, {'dataset', 'comparison_type', 'model_name', ...
    'feature_set'}, 'Before', 1);

end

function T = average_by_feature(A)
[G, T] = findgroups(A(:, {'model_name', 'feature_set', 'model_type', 'model_role'}));
T.N = splitapply(@(x) sum(x, 'omitnan'), A.N, G);
T.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), A.MAPE_pct, G);
T.RMSE_pct = splitapply(@(x) mean(x, 'omitnan'), A.RMSE_pct, G);
T.HighError_gt20_pct = splitapply(@(x) mean(x, 'omitnan'), A.HighError_gt20_pct, G);
T.P95_abs_error_pct = splitapply(@(x) mean(x, 'omitnan'), A.P95_abs_error_pct, G);
end

function row = average_metric_row(A)
row = A(1, {'model_name', 'model_type'});
if ismember('model_role', string(A.Properties.VariableNames))
    row.model_role = A.model_role(1);
else
    row.model_role = "operational";
end
row.N = sum(A.N, 'omitnan');
row.MAPE_pct = mean(A.MAPE_pct, 'omitnan');
row.RMSE_pct = mean(A.RMSE_pct, 'omitnan');
row.HighError_gt20_pct = mean(A.HighError_gt20_pct, 'omitnan');
row.P95_abs_error_pct = mean(A.P95_abs_error_pct, 'omitnan');
end

function P = join_prediction_metadata(P, T_feat)

if ~ismember('row_key', string(P.Properties.VariableNames))
    error('Prediction table does not contain row_key.');
end

keep = ["row_key"; "SIM_f0"; "SIM_cs_bg"; "REQ_M"; "REQ_cs_guess"; ...
    "M_eff_guess"; "M_eff_true_diag"; "SIM_WaveModel"; "step_idx"];
keep = keep(ismember(keep, string(T_feat.Properties.VariableNames)));
M = T_feat(:, cellstr(keep));

P = outerjoin(P, M, 'Keys', 'row_key', 'MergeKeys', true, 'Type', 'left');
P = adaptive_req.analysis.Test12Analysis.addBins(P);

if ~ismember('feature_set', string(P.Properties.VariableNames))
    P.feature_set = repmat("Test11Baseline", height(P), 1);
end
if ~ismember('model_role', string(P.Properties.VariableNames))
    P.model_role = repmat("operational", height(P), 1);
end

end

function [P11, P12, filters] = apply_matched_filters(P11, P12)

filters = strings(0, 1);
P11.dataset = repmat("Test11", height(P11), 1);
P12.dataset = repmat("Test12", height(P12), 1);

common_vars = ["SIM_f0"; "SIM_cs_bg"; "REQ_M"; "SIM_WaveModel"];
for i = 1:numel(common_vars)
    v = common_vars(i);
    if ismember(v, string(P11.Properties.VariableNames)) && ...
            ismember(v, string(P12.Properties.VariableNames))
        c = intersect(string(unique(P11.(char(v)))), string(unique(P12.(char(v)))));
        if ~isempty(c)
            P11 = P11(ismember(string(P11.(char(v))), c), :);
            P12 = P12(ismember(string(P12.(char(v))), c), :);
            filters(end + 1, 1) = sprintf('%s in [%s]', v, strjoin(c, ', ')); %#ok<AGROW>
        end
    end
end

if ismember('REQ_cs_guess', string(P12.Properties.VariableNames)) && ...
        any(abs(P12.REQ_cs_guess - 3.0) < eps)
    P12 = P12(abs(P12.REQ_cs_guess - 3.0) < eps, :);
    filters(end + 1, 1) = "Test12 restricted to REQ_cs_guess == 3.0";
end

if ismember('SIM_WaveModel', string(P12.Properties.VariableNames)) && ...
        isscalar(unique(string(P12.SIM_WaveModel))) && ...
        ismember('SIM_WaveModel', string(P11.Properties.VariableNames))
    wm = unique(string(P12.SIM_WaveModel));
    P11 = P11(string(P11.SIM_WaveModel) == wm, :);
    filters(end + 1, 1) = "Test11 restricted to Test12 wave model: " + wm;
end

P11 = P11(P11.model_type == "bagged_trees" & P11.split == "test", :);
P12 = P12(P12.model_type == "bagged_trees" & P12.split == "test" & ...
    P12.model_role == "operational", :);

end

function T = summarize_matched(P11, P12)

T = table();

if ~isempty(P11)
    T11 = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
        P11, ["dataset"; "model_name"; "feature_set"; "model_type"]);
    T11.comparison_type = repmat("matched", height(T11), 1);
    T = adaptive_req.analysis.Test12Analysis.concatTables(T, T11);
end

if ~isempty(P12)
    T12 = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
        P12, ["dataset"; "model_name"; "feature_set"; "model_type"; "model_role"]);
    T12.comparison_type = repmat("matched", height(T12), 1);
    T = adaptive_req.analysis.Test12Analysis.concatTables(T, T12);
end

T = movevars(T, {'dataset', 'comparison_type', 'model_name', 'feature_set'}, 'Before', 1);

end

function T = summarize_ablation(T_l01, T_l02)

T_l01.analysis_level = repmat("Level01", height(T_l01), 1);
T_l02.analysis_level = repmat("Level02", height(T_l02), 1);
T = adaptive_req.analysis.Test12Analysis.concatTables(T_l01, T_l02);
T = T(T.model_type == "bagged_trees" & T.model_role == "operational", :);
keep = ["analysis_level"; "generalization_test"; "heldout_value"; ...
    "model_name"; "feature_set"; "Delta_MAPE_vs_NoCsGuess"; ...
    "Delta_HighError_gt20_vs_NoCsGuess"; "Delta_RMSE_vs_NoCsGuess"];
keep = keep(ismember(keep, string(T.Properties.VariableNames)));
T = T(:, cellstr(keep));
T = sortrows(T, 'Delta_MAPE_vs_NoCsGuess', 'ascend');

end

function T = build_outlier_comparison(P11, P12)

P11.dataset = repmat("Test11", height(P11), 1);
P12.dataset = repmat("Test12", height(P12), 1);
P = adaptive_req.analysis.Test12Analysis.concatTables(P11, P12);
P = adaptive_req.analysis.Test12Analysis.addBins(P);

T = table();
bin_vars = ["REQ_M"; "M_eff_guess_bin"; "M_eff_true_diag_bin"];
for i = 1:numel(bin_vars)
    v = bin_vars(i);
    if ~ismember(v, string(P.Properties.VariableNames))
        continue;
    end
    Ti = adaptive_req.analysis.Test12Analysis.summarizeSws( ...
        P, ["dataset"; "model_name"; "feature_set"; "model_type"; v]);
    Ti.bin_variable = repmat(v, height(Ti), 1);
    Ti.bin_value = string(Ti.(char(v)));
    T = adaptive_req.analysis.Test12Analysis.concatTables(T, Ti);
end

end

function T = build_outlier_comparison_from_metrics(T11, T12)
T11.dataset = repmat("Test11", height(T11), 1);
T11.feature_set = repmat("Test11Baseline", height(T11), 1);
T12.dataset = repmat("Test12", height(T12), 1);
T = adaptive_req.analysis.Test12Analysis.concatTables(T11, T12);
end

function plot_global_bars(T, metric_var, file_path)

if isempty(T) || ~ismember(metric_var, string(T.Properties.VariableNames))
    return;
end
T.series = T.dataset + " | " + T.feature_set;
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "model_name", "series", metric_var, ...
    sprintf('Test 11 vs Test 12: %s', metric_var), file_path);

end

function plot_ablation(T, metric_var, file_path)

if isempty(T) || ~ismember(metric_var, string(T.Properties.VariableNames))
    return;
end
T = T(T.analysis_level == "Level02", :);
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "feature_set", "model_name", metric_var, ...
    sprintf('Test 12 ablation: %s', metric_var), file_path);

end

function plot_outlier_bins(T, bin_var, file_path)

if isempty(T) || ~ismember('bin_variable', string(T.Properties.VariableNames))
    return;
end
T = T(T.bin_variable == bin_var & T.model_type == "bagged_trees", :);
if isempty(T)
    return;
end
T.series = T.dataset + " | " + T.model_name;
adaptive_req.analysis.Test12Analysis.plotMetricByTwoGroups( ...
    T, "bin_value", "series", "HighError_gt20_pct", ...
    sprintf('Outliers by %s', bin_var), file_path);

end

function write_markdown_summary(file_path, T_global, T_matched, T_ablation, matched_filters)

fid = fopen(file_path, 'w');
if fid < 0
    error('Could not write markdown summary: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Test 12: análisis del barrido de `cs_guess` y ventana efectiva\n\n');
fprintf(fid, 'Este documento fue generado por `analyze_test_12_level03_compare_test11.m`.\n\n');
fprintf(fid, '## Qué se hizo\n\n');
fprintf(fid, ['Test 12 expandió el diseño de Test 11 variando `REQ_cs_guess`, ', ...
    '`REQ_M`, `SIM_cs_bg`, `SIM_f0` y el modelo de onda. El objetivo fue ', ...
    'desacoplar el tamaño nominal de ventana del tamaño efectivo real.\n\n']);
fprintf(fid, '## Modelos y feature sets\n\n');
fprintf(fid, ['Se evaluaron `LocalOnly`, `GlobalOnly` y `HybridLocalGlobal` con ', ...
    '`NoCsGuess`, `WithCsGuess`, `WithMeffGuess`, ', ...
    '`WithCsGuessAndMeffGuess` y un set diagnóstico `DiagnosticWithMeffTrue`.\n\n']);
fprintf(fid, '`M_eff_true_diag` se considera diagnóstico/oracle y no debe usarse como predictor operacional.\n\n');

fprintf(fid, '## Comparación global Test 11 vs Test 12\n\n');
write_table_preview(fid, T_global);

fprintf(fid, '\n## Comparación matched/subset\n\n');
if isempty(T_matched)
    fprintf(fid, 'TODO: correr Level 01/02 y verificar que existan las predicciones para calcular la comparación matched.\n\n');
else
    fprintf(fid, 'Filtros usados:\n\n');
    for i = 1:numel(matched_filters)
        fprintf(fid, '- %s\n', matched_filters(i));
    end
    fprintf(fid, '\n');
    write_table_preview(fid, T_matched);
end

fprintf(fid, '\n## Ablation de `REQ_cs_guess` y `M_eff_guess`\n\n');
fprintf(fid, 'Delta negativo significa mejora frente a `NoCsGuess`.\n\n');
write_table_preview(fid, T_ablation);

fprintf(fid, '\n## Interpretación pendiente\n\n');
fprintf(fid, ['No se deben inventar conclusiones antes de correr Level 01, Level 02 y Level 03. ', ...
    'Usar las tablas en `outputs/test_12_cs_guess_window_sweep/.../analysis/` para cerrar si ', ...
    '`REQ_cs_guess`, `M_eff_guess` o su combinación mejoran realmente el modelo.\n\n']);
fprintf(fid, '## Próximos pasos sugeridos\n\n');
fprintf(fid, '1. Crear Test 13 para SNR.\n');
fprintf(fid, '2. Agregar SNR dependiente de profundidad.\n');
fprintf(fid, '3. Agregar atenuación.\n');
fprintf(fid, '4. Volver a bilayer.\n');
fprintf(fid, '5. Validar con k-Wave.\n');

end

function write_table_preview(fid, T)

if isempty(T)
    fprintf(fid, 'TODO: tabla pendiente de generar.\n\n');
    return;
end

T_preview = T(1:min(20, height(T)), :);
fprintf(fid, '```text\n');
disp_text = formattedDisplayText(T_preview);
fprintf(fid, '%s', disp_text);
fprintf(fid, '```\n');

end
