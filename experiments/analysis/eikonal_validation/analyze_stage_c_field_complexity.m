%% Analyze Stage C: Eikonal field-complexity transfer validation
% Purpose:
%   Summarize and visualize clean Eikonal transfer results for the frozen
%   baseline_minimal_v1 q/SWS models. Stage C isolates field complexity by
%   varying the source count/layout while keeping realism clean and M fixed.

clear; clc;

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if isempty(ROOT) || ~isfolder(fullfile(ROOT, 'configs'))
    ROOT = pwd;
end
addpath(fullfile(ROOT, 'src'));

CFG = load_stage_c_config(ROOT);
OUT = fullfile(ROOT, CFG.OutputRoot);
TBL_DIR = fullfile(OUT, 'tables');
FIG_DIR = fullfile(OUT, 'figures');
DOC_DIR = fullfile(ROOT, 'docs', 'eikonal_validation');
ensure_dir(TBL_DIR); ensure_dir(FIG_DIR); ensure_dir(DOC_DIR);

patch_csv = fullfile(TBL_DIR, 'stage_c_patch_level_results.csv');
assert(exist(patch_csv, 'file') == 2, 'Missing Stage C patch-level table: %s', patch_csv);

T = readtable(patch_csv, 'TextType', 'string');
if isempty(T)
    error('Stage C patch-level table is empty: %s', patch_csv);
end
T = normalize_stage_c_table(T);

primary_model = "q_spectrum_plus_composition";
models = unique(T.model_name, 'stable');
fprintf('\nStage C analysis: loaded %d patch/model rows across %d conditions.\n', height(T), numel(unique(T.condition_id)));
fprintf('Models: %s\n', strjoin(models, ', '));

% Recompute compact summary tables so the analysis can be rerun independently
% from the runner and still refresh paper-facing diagnostics.
write_summary_tables(T, TBL_DIR);

% Figures: small titles, plain readable labels, and no raw variable names in
% titles/legends. Use PNG plus MATLAB FIG for key figures when useful.
make_global_model_figure(T, FIG_DIR);
make_source_count_figure(T, primary_model, FIG_DIR);
make_field_regime_figure(T, primary_model, FIG_DIR);
make_geometry_field_heatmaps(T, primary_model, FIG_DIR);
make_roi_field_heatmaps(T, primary_model, FIG_DIR);
make_frequency_source_count_figure(T, primary_model, FIG_DIR);
make_model_comparison_figure(T, FIG_DIR);
make_q_source_count_figure(T, primary_model, FIG_DIR);
make_hard_bias_focus_figure(T, primary_model, FIG_DIR);
make_worst_condition_figure(T, primary_model, FIG_DIR);

write_stage_c_readme(T, OUT, DOC_DIR, CFG, primary_model);
print_console_summary(T, primary_model);

fprintf('\nStage C analysis complete.\nTables: %s\nFigures: %s\nREADME: %s\n', TBL_DIR, FIG_DIR, fullfile(OUT, 'README_results.md'));

%% Local functions
function CFG = load_stage_c_config(ROOT)
    cfg_path = fullfile(ROOT, 'configs', 'eikonal_validation', 'stage_c_field_complexity.json');
    assert(exist(cfg_path, 'file') == 2, 'Missing Stage C config: %s', cfg_path);
    CFG = jsondecode(fileread(cfg_path));
end

function T = normalize_stage_c_table(T)
    if ~ismember('ROI', T.Properties.VariableNames) && ismember('roi_label', T.Properties.VariableNames)
        T.ROI = T.roi_label;
    end
    if ~ismember('layout_id', T.Properties.VariableNames) && ismember('source_layout_id', T.Properties.VariableNames)
        T.layout_id = T.source_layout_id;
    end
    string_vars = ["condition_id","geometry","field_regime","model_name","ROI","layout_id","layout_family","realism_level"];
    for v = string_vars
        if ismember(v, T.Properties.VariableNames)
            T.(v) = string(T.(v));
        end
    end
    if ismember('source_seed', T.Properties.VariableNames)
        T.source_seed = double(T.source_seed);
    end
    required = ["condition_id","geometry","field_regime","model_name","f0","N_sources","N_in_plane_sources", ...
        "SWS_true","SWS_pred","q_pred","signed_error_pct","abs_error_pct","ROI"];
    missing = setdiff(required, string(T.Properties.VariableNames));
    assert(isempty(missing), 'Missing required Stage C columns: %s', strjoin(missing, ', '));

    T.geometry_pretty = arrayfun(@pretty_geometry, T.geometry);
    T.field_pretty = arrayfun(@pretty_field, T.field_regime);
    T.roi_pretty = arrayfun(@pretty_roi, T.ROI);
    T.model_pretty = arrayfun(@pretty_model, T.model_name);
    T.source_count_label = string(T.N_sources) + " sources";
    T.field_order = field_order(T.field_regime);
    T.geometry_order = arrayfun(@geometry_order, T.geometry);
    T.roi_order = roi_order(T.ROI);
    T.is_primary = T.model_name == "q_spectrum_plus_composition";
end

function write_summary_tables(T, tbl_dir)
    write_table(summarize_metrics(T, ["model_name"]), tbl_dir, 'stage_c_summary_by_model.csv');
    write_table(summarize_metrics(T, ["model_name","geometry"]), tbl_dir, 'stage_c_summary_by_geometry.csv');
    write_table(summarize_metrics(T, ["model_name","field_regime","N_sources","N_in_plane_sources"]), tbl_dir, 'stage_c_summary_by_field_regime.csv');
    write_table(summarize_metrics(T, ["model_name","N_sources"]), tbl_dir, 'stage_c_summary_by_source_count.csv');
    write_table(summarize_metrics(T, ["model_name","N_in_plane_sources"]), tbl_dir, 'stage_c_summary_by_in_plane_source_count.csv');
    write_table(summarize_metrics(T, ["model_name","f0"]), tbl_dir, 'stage_c_summary_by_frequency.csv');
    write_table(summarize_metrics(T, ["model_name","geometry","field_regime"]), tbl_dir, 'stage_c_summary_by_geometry_field.csv');
    write_table(summarize_metrics(T, ["model_name","geometry","f0"]), tbl_dir, 'stage_c_summary_by_geometry_frequency.csv');
    write_table(summarize_metrics(T, ["model_name","ROI"]), tbl_dir, 'stage_c_summary_by_roi.csv');
    write_table(summarize_metrics(T, ["model_name","ROI","field_regime"]), tbl_dir, 'stage_c_summary_by_roi_field.csv');
    write_table(summarize_metrics(T, ["model_name","ROI","N_sources"]), tbl_dir, 'stage_c_summary_by_roi_source_count.csv');
    write_table(summarize_metrics(T, ["model_name","geometry","ROI","field_regime"]), tbl_dir, 'stage_c_summary_by_geometry_roi_field.csv');

    cond = summarize_metrics(T, ["model_name","condition_id","geometry","field_regime","f0","N_sources","N_in_plane_sources"]);
    cond = sortrows(cond, 'MAPE_pct', 'descend');
    write_table(cond(1:min(30,height(cond)),:), tbl_dir, 'stage_c_worst_conditions.csv');
end

function S = summarize_metrics(T, groups)
    groups = string(groups);
    [G, key] = findgroups(T(:, cellstr(groups)));
    S = key;
    S.N_valid_patches = splitapply(@numel, T.SWS_pred, G);
    S.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), T.abs_error_pct, G);
    S.signed_bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.signed_error_pct, G);
    S.median_abs_error_pct = splitapply(@(x) median(x, 'omitnan'), T.abs_error_pct, G);
    S.high_error_10_pct = splitapply(@(x) 100*mean(x > 10, 'omitnan'), T.abs_error_pct, G);
    S.high_error_20_pct = splitapply(@(x) 100*mean(x > 20, 'omitnan'), T.abs_error_pct, G);
    S.mean_predicted_SWS = splitapply(@(x) mean(x, 'omitnan'), T.SWS_pred, G);
    S.mean_true_SWS = splitapply(@(x) mean(x, 'omitnan'), T.SWS_true, G);
    S.mean_q_pred = splitapply(@(x) mean(x, 'omitnan'), T.q_pred, G);
end

function make_global_model_figure(T, fig_dir)
    S = summarize_metrics(T, ["model_name"]);
    S.model_pretty = arrayfun(@pretty_model, string(S.model_name));
    [~,ord] = sort(string(S.model_name)); S = S(ord,:);
    fig = newfig([950 430]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; bar(categorical(S.model_pretty), S.MAPE_pct); grid on;
    ylabel('MAPE (%)'); title('Global MAPE by model'); xtickangle(20);
    nexttile; bar(categorical(S.model_pretty), S.signed_bias_pct); grid on; yline(0,'k-', 'HandleVisibility','off');
    ylabel('Signed bias (%)'); title('Global signed bias by model'); xtickangle(20);
    savefigures(fig, fig_dir, 'stage_C_global_model_summary');
end

function make_source_count_figure(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["geometry","N_sources"]);
    S.geometry_pretty = arrayfun(@pretty_geometry, string(S.geometry));
    geoms = sort_by_pretty(unique(string(S.geometry)), @geometry_order);
    fig = newfig([1050 760]);
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; grid on;
    for g = geoms'
        R = S(string(S.geometry)==g,:); R = sortrows(R,'N_sources');
        plot(R.N_sources, R.MAPE_pct, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', pretty_geometry(g));
    end
    xlabel('Number of sources'); ylabel('MAPE (%)'); title('MAPE versus field source count'); legend('Location','eastoutside'); xticks(unique(S.N_sources));
    nexttile; hold on; grid on;
    for g = geoms'
        R = S(string(S.geometry)==g,:); R = sortrows(R,'N_sources');
        plot(R.N_sources, R.signed_bias_pct, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', pretty_geometry(g));
    end
    yline(0,'k-', 'HandleVisibility','off'); xlabel('Number of sources'); ylabel('Signed bias (%)'); title('Signed bias versus field source count'); legend('Location','eastoutside'); xticks(unique(S.N_sources));
    savefigures(fig, fig_dir, 'stage_C_mape_bias_vs_source_count');
end

function make_field_regime_figure(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["field_regime","N_sources","N_in_plane_sources"]);
    S.field_pretty = arrayfun(@pretty_field, string(S.field_regime));
    S.field_order = field_order(string(S.field_regime));
    S = sortrows(S, 'field_order');
    fig = newfig([1200 520]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; bar(categorical(S.field_pretty, S.field_pretty), S.MAPE_pct); grid on;
    ylabel('MAPE (%)'); title('MAPE by field regime'); xtickangle(25);
    nexttile; bar(categorical(S.field_pretty, S.field_pretty), S.signed_bias_pct); grid on; yline(0,'k-', 'HandleVisibility','off');
    ylabel('Signed bias (%)'); title('Signed bias by field regime'); xtickangle(25);
    savefigures(fig, fig_dir, 'stage_C_field_regime_mape_bias');
end

function make_geometry_field_heatmaps(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["geometry","field_regime"]);
    rows = sort_by_pretty(unique(string(S.geometry)), @geometry_order);
    cols = sort_by_pretty(unique(string(S.field_regime)), @field_order);
    plot_heat(S, rows, cols, 'geometry', 'field_regime', 'MAPE_pct', 'MAPE (%)', 'Geometry x field MAPE', fullfile(fig_dir, 'stage_C_heatmap_geometry_field_mape'));
    plot_heat(S, rows, cols, 'geometry', 'field_regime', 'signed_bias_pct', 'Signed bias (%)', 'Geometry x field signed bias', fullfile(fig_dir, 'stage_C_heatmap_geometry_field_bias'));
    plot_heat(S, rows, cols, 'geometry', 'field_regime', 'high_error_20_pct', 'Pixels >20% error (%)', 'Geometry x field high-error rate', fullfile(fig_dir, 'stage_C_heatmap_geometry_field_high20'));
end

function make_roi_field_heatmaps(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["ROI","field_regime"]);
    rows = sort_by_pretty(unique(string(S.ROI)), @roi_order);
    cols = sort_by_pretty(unique(string(S.field_regime)), @field_order);
    plot_heat(S, rows, cols, 'ROI', 'field_regime', 'MAPE_pct', 'MAPE (%)', 'ROI x field MAPE', fullfile(fig_dir, 'stage_C_heatmap_roi_field_mape'), @pretty_roi, @pretty_field);
    plot_heat(S, rows, cols, 'ROI', 'field_regime', 'signed_bias_pct', 'Signed bias (%)', 'ROI x field signed bias', fullfile(fig_dir, 'stage_C_heatmap_roi_field_bias'), @pretty_roi, @pretty_field);
end

function make_frequency_source_count_figure(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["f0","N_sources"]);
    counts = unique(S.N_sources)';
    fig = newfig([1000 720]);
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
    nexttile; hold on; grid on;
    for n = counts
        R = sortrows(S(S.N_sources==n,:), 'f0');
        plot(R.f0, R.MAPE_pct, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', sprintf('%d sources', n));
    end
    xlabel('Frequency (Hz)'); ylabel('MAPE (%)'); title('MAPE versus frequency by source count'); legend('Location','eastoutside'); xticks(unique(S.f0));
    nexttile; hold on; grid on;
    for n = counts
        R = sortrows(S(S.N_sources==n,:), 'f0');
        plot(R.f0, R.signed_bias_pct, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', sprintf('%d sources', n));
    end
    yline(0,'k-', 'HandleVisibility','off'); xlabel('Frequency (Hz)'); ylabel('Signed bias (%)'); title('Signed bias versus frequency by source count'); legend('Location','eastoutside'); xticks(unique(S.f0));
    savefigures(fig, fig_dir, 'stage_C_frequency_dependence_by_source_count');
end

function make_model_comparison_figure(T, fig_dir)
    S = summarize_metrics(T, ["model_name","geometry"]);
    S.model_pretty = arrayfun(@pretty_model, string(S.model_name));
    S.geometry_pretty = arrayfun(@pretty_geometry, string(S.geometry));
    geoms = sort_by_pretty(unique(string(S.geometry)), @geometry_order);
    mods = unique(string(S.model_name), 'stable');
    X = categorical(arrayfun(@pretty_geometry, geoms));
    Y = nan(numel(geoms), numel(mods));
    B = nan(numel(geoms), numel(mods));
    for i=1:numel(geoms)
        for j=1:numel(mods)
            r = S(string(S.geometry)==geoms(i) & string(S.model_name)==mods(j),:);
            if ~isempty(r), Y(i,j)=r.MAPE_pct(1); B(i,j)=r.signed_bias_pct(1); end
        end
    end
    fig = newfig([1120 520]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; bar(X,Y); grid on; ylabel('MAPE (%)'); title('Model comparison by geometry'); xtickangle(20); legend(arrayfun(@pretty_model, mods),'Location','eastoutside');
    nexttile; bar(X,B); grid on; yline(0,'k-', 'HandleVisibility','off'); ylabel('Signed bias (%)'); title('Bias comparison by geometry'); xtickangle(20); legend(arrayfun(@pretty_model, mods),'Location','eastoutside');
    savefigures(fig, fig_dir, 'stage_C_model_comparison_by_geometry');
end

function make_q_source_count_figure(T, model, fig_dir)
    P = T(T.model_name == model,:);
    S = summarize_metrics(P, ["geometry","N_sources"]);
    geoms = sort_by_pretty(unique(string(S.geometry)), @geometry_order);
    fig = newfig([980 460]); hold on; grid on;
    for g = geoms'
        R = sortrows(S(string(S.geometry)==g,:), 'N_sources');
        plot(R.N_sources, R.mean_q_pred, '-o', 'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', pretty_geometry(g));
    end
    xlabel('Number of sources'); ylabel('Mean predicted REQ quantile q');
    title('Predicted q versus field source count'); legend('Location','eastoutside'); xticks(unique(S.N_sources));
    savefigures(fig, fig_dir, 'stage_C_q_pred_vs_source_count');
end

function make_hard_bias_focus_figure(T, model, fig_dir)
    P = T(T.model_name == model & ismember(T.geometry, ["homogeneous_cs4","inclusion_2_4_D24"]),:);
    if isempty(P), return; end
    S = summarize_metrics(P, ["geometry","ROI","N_sources"]);
    keep = string(S.ROI)=="homogeneous_center" | string(S.ROI)=="inclusion_core" | string(S.ROI)=="inclusion_center";
    S = S(keep,:);
    if isempty(S), return; end
    S.series = arrayfun(@(g,r) pretty_geometry(g) + " | " + pretty_roi(r), string(S.geometry), string(S.ROI));
    series = unique(S.series,'stable');
    fig = newfig([980 500]); hold on; grid on;
    for s = series'
        R = sortrows(S(S.series==s,:), 'N_sources');
        plot(R.N_sources, R.signed_bias_pct, '-o', 'LineWidth', 1.7, 'MarkerSize', 6, 'DisplayName', s);
    end
    yline(0,'k-', 'HandleVisibility','off'); xlabel('Number of sources'); ylabel('Signed bias (%)');
    title('Hard-speed bias focus: homogeneous 4 m/s vs inclusion core');
    legend('Location','eastoutside'); xticks(unique(S.N_sources));
    savefigures(fig, fig_dir, 'stage_C_hard_bias_focus');
end

function make_worst_condition_figure(T, model, fig_dir)
    S = summarize_metrics(T(T.model_name==model,:), ["condition_id","geometry","field_regime","f0","N_sources"]);
    S = sortrows(S, 'MAPE_pct', 'descend');
    S = S(1:min(12,height(S)),:);
    labels = arrayfun(@(g,f,n) sprintf('%s | %d Hz | %d src', pretty_geometry(g), f, n), string(S.geometry), S.f0, S.N_sources, 'UniformOutput', false);
    fig = newfig([1150 560]);
    barh(categorical(labels, flip(labels)), flip(S.MAPE_pct)); grid on;
    xlabel('MAPE (%)'); title('Worst Stage C conditions for q spectrum + composition');
    savefigures(fig, fig_dir, 'stage_C_worst_conditions');
end

function plot_heat(S, rows, cols, rowvar, colvar, metric, cb_label, title_text, outbase, rowpretty, colpretty)
    if nargin < 10, rowpretty = @(x) pretty_generic(x); end
    if nargin < 11, colpretty = @(x) pretty_generic(x); end
    Z = nan(numel(rows), numel(cols));
    for i=1:numel(rows)
        for j=1:numel(cols)
            r = S(string(S.(rowvar))==rows(i) & string(S.(colvar))==cols(j),:);
            if ~isempty(r), Z(i,j) = r.(metric)(1); end
        end
    end
    fig = newfig([1150 620]);
    imagesc(Z); axis tight; grid off;
    colormap(parula); cb = colorbar; cb.Label.String = cb_label;
    xticks(1:numel(cols)); xticklabels(arrayfun(colpretty, cols)); xtickangle(25);
    yticks(1:numel(rows)); yticklabels(arrayfun(rowpretty, rows));
    title(title_text, 'FontSize', 14);
    set(gca, 'FontSize', 10);
    for i=1:numel(rows)
        for j=1:numel(cols)
            if isfinite(Z(i,j))
                text(j, i, sprintf('%.2f', Z(i,j)), 'HorizontalAlignment','center', 'Color','w', 'FontWeight','bold', 'FontSize', 8);
            end
        end
    end
    savefigures(fig, fileparts(outbase), string(get_filename(outbase)));
end

function write_stage_c_readme(T, out_dir, doc_dir, CFG, primary_model)
    S_model = summarize_metrics(T, ["model_name"]);
    S_geom = summarize_metrics(T(T.model_name==primary_model,:), ["geometry"]);
    S_field = summarize_metrics(T(T.model_name==primary_model,:), ["field_regime","N_sources","N_in_plane_sources"]);
    S_freq = summarize_metrics(T(T.model_name==primary_model,:), ["f0"]);
    S_worst = summarize_metrics(T(T.model_name==primary_model,:), ["condition_id","geometry","field_regime","f0","N_sources"]);
    S_worst = sortrows(S_worst, 'MAPE_pct', 'descend');
    if isempty(S_worst)
        worst_line = 'No worst-condition summary was available.';
    else
        w = S_worst(1,:);
        worst_line = sprintf('Worst primary-model condition: `%s`, %s, %d Hz, %s, MAPE %.2f%%, bias %.2f%%.', ...
            w.condition_id, pretty_geometry(w.geometry), w.f0, pretty_field(w.field_regime), w.MAPE_pct, w.signed_bias_pct);
    end

    txt = strings(0,1);
    txt(end+1) = "# Stage C: Field Complexity Transfer Validation";
    txt(end+1) = "";
    txt(end+1) = "Stage C evaluates whether the frozen `baseline_minimal_v1` q/SWS models remain stable when clean Eikonal fields become less directional and more diffuse-like. No model is retrained, and no readout noise, risk mask, correction, or reliability layer is used.";
    txt(end+1) = "";
    txt(end+1) = "## Design";
    geom_ids = string({CFG.Geometries.id});
    freq_vals = double(CFG.FrequenciesHz);
    field_vals = string(CFG.FieldRegimes);
    req_M = double(CFG.MList(1));
    txt(end+1) = sprintf('- Geometries: %s.', strjoin(geom_ids, ', '));
    txt(end+1) = sprintf('- Frequencies: %s Hz.', strjoin(string(freq_vals), ', '));
    txt(end+1) = sprintf('- Field regimes: %s.', strjoin(field_vals, ', '));
    txt(end+1) = sprintf('- REQ: M=%g, cs_guess=%.2f m/s, target step %.2f mm, valid windows only.', req_M, CFG.REQ.cs_guess, 1e3*CFG.REQ.TargetStepM);
    txt(end+1) = "- Clean-only fields are globally RMS-normalized over the central evaluation region so source-count effects are not confounded with global amplitude scaling. Local interference patterns are preserved.";
    txt(end+1) = "";
    txt(end+1) = "## Models";
    txt(end+1) = "- Primary: `q_spectrum_plus_composition`.";
    txt(end+1) = "- Diagnostic: `q_spectrum_only`.";
    txt(end+1) = "";
    txt(end+1) = "## Main Results";
    txt(end+1) = table_to_markdown(S_model, ["model_name","N_valid_patches","MAPE_pct","signed_bias_pct","high_error_20_pct"]);
    txt(end+1) = "";
    txt(end+1) = "### Primary Model By Geometry";
    txt(end+1) = table_to_markdown(S_geom, ["geometry","MAPE_pct","signed_bias_pct","high_error_20_pct"]);
    txt(end+1) = "";
    txt(end+1) = "### Primary Model By Field Regime";
    txt(end+1) = table_to_markdown(S_field, ["field_regime","N_sources","N_in_plane_sources","MAPE_pct","signed_bias_pct","high_error_20_pct"]);
    txt(end+1) = "";
    txt(end+1) = "### Primary Model By Frequency";
    txt(end+1) = table_to_markdown(S_freq, ["f0","MAPE_pct","signed_bias_pct","high_error_20_pct"]);
    txt(end+1) = "";
    txt(end+1) = "## Interpretation Guide";
    txt(end+1) = "- If homogeneous 3 m/s stays accurate across all source counts, the feature extraction and frozen model are stable for clean Eikonal fields without hard-speed bias.";
    txt(end+1) = "- If homogeneous 4 m/s becomes increasingly negative with source count, the issue is hard-speed plus field-regime shift, not interface mixing.";
    txt(end+1) = "- If inclusion 2/4 degrades more than homogeneous 4 m/s, field complexity is interacting with the hard interface.";
    txt(end+1) = "- If bilayer error is concentrated in interface ROIs, this supports a window-mixing interpretation.";
    txt(end+1) = "";
    txt(end+1) = "## Worst Condition";
    txt(end+1) = worst_line;
    txt(end+1) = "";
    txt(end+1) = "## Figures";
    txt(end+1) = "Key analysis figures are under `figures/`, including source-count trends, geometry-field heatmaps, ROI-field heatmaps, frequency trends, model comparisons, and representative maps from the runner.";
    txt(end+1) = "";
    txt(end+1) = "## How To Run";
    txt(end+1) = "```matlab";
    txt(end+1) = "setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','validate_only'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')";
    txt(end+1) = "setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','quick'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')";
    txt(end+1) = "setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','full'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')";
    txt(end+1) = "run('experiments/analysis/eikonal_validation/analyze_stage_c_field_complexity.m')";
    txt(end+1) = "```";

    outfile = fullfile(out_dir, 'README_results.md');
    fid = fopen(outfile, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
    docfile = fullfile(doc_dir, 'stage_C_field_complexity.md');
    fid = fopen(docfile, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end

function print_console_summary(T, primary_model)
    S_model = summarize_metrics(T, ["model_name"]);
    disp('Overall model summary:'); disp(S_model(:, intersect(string(S_model.Properties.VariableNames), ["model_name","N_valid_patches","MAPE_pct","signed_bias_pct","high_error_20_pct"], 'stable')));
    S_geom = summarize_metrics(T(T.model_name==primary_model,:), ["geometry"]);
    disp('Primary model by geometry:'); disp(S_geom(:, ["geometry","MAPE_pct","signed_bias_pct","high_error_20_pct"]));
    S_field = summarize_metrics(T(T.model_name==primary_model,:), ["field_regime","N_sources"]);
    disp('Primary model by field/source count:'); disp(S_field(:, ["field_regime","N_sources","MAPE_pct","signed_bias_pct","high_error_20_pct"]));
end

function write_table(T, dirpath, name)
    ensure_dir(dirpath);
    writetable(T, fullfile(dirpath, name));
end

function ensure_dir(p)
    if ~exist(p, 'dir'), mkdir(p); end
end

function fig = newfig(pos)
    fig = figure('Visible','off','Color','w','Position',[100 100 pos]);
    set(fig, 'DefaultAxesFontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman');
    set(fig, 'DefaultAxesFontSize', 11);
    set(fig, 'DefaultTextFontSize', 11);
    set(fig, 'DefaultAxesTitleFontSizeMultiplier', 1.05);
    set(fig, 'DefaultAxesLabelFontSizeMultiplier', 1.0);
end

function savefigures(fig, fig_dir, stem)
    ensure_dir(fig_dir);
    stem = char(stem);
    exportgraphics(fig, fullfile(fig_dir, [stem '.png']), 'Resolution', 220);
    try
        savefig(fig, fullfile(fig_dir, [stem '.fig']));
    catch
    end
    close(fig);
end

function name = get_filename(p)
    [~,name] = fileparts(char(p));
end

function y = sort_by_pretty(x, orderfun)
    [~,ord] = sort(arrayfun(orderfun, x));
    y = x(ord);
end

function o = geometry_order(g)
    g = string(g);
    names = ["homogeneous_cs3","homogeneous_cs4","inclusion_2_3_D24","inclusion_2_4_D24","bilayer_2_3"];
    idx = find(names==g, 1);
    if isempty(idx), o = 100; else, o = idx; end
end

function o = field_order(f)
    f = string(f);
    names = ["single_source_lateral", "diffuse_like_4src_layoutA", "diffuse_like_8src_layoutA", "diffuse_like_16src_layoutA", ...
             "diffuse_like_4src_layoutB", "diffuse_like_8src_layoutB", "diffuse_like_16src_layoutB"];
    o = zeros(size(f));
    for k=1:numel(f)
        idx = find(names==f(k), 1);
        if isempty(idx), o(k)=100; else, o(k)=idx; end
    end
end

function o = roi_order(r)
    r = string(r);
    names = ["homogeneous_center","background_far","soft_layer_far","hard_layer_far","inclusion_center","inclusion_core", ...
             "interface_0_1mm","interface_1_2mm","interface_2_4mm","interface_4_8mm","unassigned"];
    o = zeros(size(r));
    for k=1:numel(r)
        idx = find(names==r(k), 1);
        if isempty(idx), o(k)=100; else, o(k)=idx; end
    end
end

function s = pretty_model(name)
    name = string(name);
    switch name
        case "q_spectrum_plus_composition"
            s = "q spectrum + composition";
        case "q_spectrum_only"
            s = "q spectrum only";
        otherwise
            s = strrep(name, '_', ' ');
    end
end

function s = pretty_geometry(g)
    g = string(g);
    switch g
        case "homogeneous_cs3", s = "homogeneous 3 m/s";
        case "homogeneous_cs4", s = "homogeneous 4 m/s";
        case "inclusion_2_3_D24", s = "inclusion 2/3, D=24 mm";
        case "inclusion_2_4_D24", s = "inclusion 2/4, D=24 mm";
        case "bilayer_2_3", s = "bilayer 2/3";
        otherwise, s = strrep(g, '_', ' ');
    end
end

function s = pretty_field(f)
    f = string(f);
    switch f
        case "single_source_lateral", s = "single lateral";
        case "diffuse_like_4src_layoutA", s = "4 sources, layout A";
        case "diffuse_like_8src_layoutA", s = "8 sources, layout A";
        case "diffuse_like_16src_layoutA", s = "16 sources, layout A";
        case "diffuse_like_4src_layoutB", s = "4 sources, layout B";
        case "diffuse_like_8src_layoutB", s = "8 sources, layout B";
        case "diffuse_like_16src_layoutB", s = "16 sources, layout B";
        otherwise, s = strrep(f, '_', ' ');
    end
end

function s = pretty_roi(r)
    r = string(r);
    switch r
        case "homogeneous_center", s = "homogeneous center";
        case "background_far", s = "background far";
        case "soft_layer_far", s = "soft layer far";
        case "hard_layer_far", s = "hard layer far";
        case "inclusion_center", s = "inclusion center";
        case "inclusion_core", s = "inclusion core";
        case "interface_0_1mm", s = "interface 0-1 mm";
        case "interface_1_2mm", s = "interface 1-2 mm";
        case "interface_2_4mm", s = "interface 2-4 mm";
        case "interface_4_8mm", s = "interface 4-8 mm";
        case "unassigned", s = "unassigned";
        otherwise, s = strrep(r, '_', ' ');
    end
end

function s = pretty_generic(x)
    x = string(x);
    if startsWith(x, "diffuse_like")
        s = pretty_field(x);
    elseif contains(x, "inclusion") || contains(x, "homogeneous") || contains(x, "bilayer")
        s = pretty_geometry(x);
    else
        s = strrep(x, '_', ' ');
    end
end

function md = table_to_markdown(T, cols)
    cols = string(cols);
    cols = cols(ismember(cols, string(T.Properties.VariableNames)));
    if isempty(T) || isempty(cols)
        md = "No rows available.";
        return;
    end
    U = T(:, cellstr(cols));
    maxRows = min(height(U), 20);
    U = U(1:maxRows,:);
    lines = strings(maxRows+2,1);
    header = "| " + strjoin(cols, " | ") + " |";
    sep = "| " + strjoin(repmat("---", size(cols)), " | ") + " |";
    lines(1) = header; lines(2) = sep;
    for i=1:maxRows
        vals = strings(1,numel(cols));
        for j=1:numel(cols)
            v = U.(cols(j))(i);
            if isnumeric(v)
                vals(j) = sprintf('%.3g', v);
            else
                vals(j) = string(v);
            end
        end
        lines(i+2) = "| " + strjoin(vals, " | ") + " |";
    end
    md = strjoin(lines, newline);
end
