%% analyze_test_12_level07_kwave_roi_diagnostics.m
% Test 12 Level 07: k-Wave ROI diagnostics.
%
% This script does not retrain models and does not modify Level 06. It reads
% the Level 06 k-Wave transfer outputs, compares adaptive-q against the
% discrete theoretical quantile and fixed quantiles, and diagnoses region
% errors and q-gradient/error relationships.

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
adaptive_req.templates.setup_style();

%% Locate Level 06 outputs

PATHS12 = locate_latest_test12_paths(root_dir);
level06_dir = fullfile(PATHS12.analysis_dir, 'level_06_kwave_transfer');
level06_file = fullfile(level06_dir, 'data', ...
    'level12_level06_kwave_transfer.mat');
assert(exist(level06_file, 'file') == 2, ...
    'Level 06 output not found:\n%s', level06_file);

analysis_dir = fullfile(PATHS12.analysis_dir, ...
    'level_07_kwave_roi_diagnostics');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
data_dir = fullfile(analysis_dir, 'data');
make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(data_dir);

S = load(level06_file, 'KW', 'CASE_OUT', 'T_all_pred', 'MODEL_INFO');
KW = S.KW;
CASE_OUT = S.CASE_OUT;
T_adaptive = S.T_all_pred;

fprintf('\nLoaded Level 06 k-Wave outputs:\n%s\n', level06_file);
fprintf('Available Level 06 case/M outputs: %d\n', numel(CASE_OUT));

available_cases = string(arrayfun(@(c) c.case.case_name, ...
    CASE_OUT(:), 'UniformOutput', false));
if ~any(contains(lower(available_cases), "partially")) && ...
        ~any(contains(lower(available_cases), "rev"))
    warning(['3D-rev was requested, but it is not present in the Level 06 ', ...
        'outputs. Level 07 will analyze the available Level 06 cases only.']);
end

%% Build method-level prediction table

fixed_q_values = [0.50 0.60 0.70 0.80 0.90];
T_methods = table();
T_grad = table();
theory_rows = table();

for oi = 1:numel(CASE_OUT)
    out_i = CASE_OUT(oi);
    T_feat = out_i.T_feat;
    T_adapt_i = out_i.T_pred;
    cfg_i = out_i.cfg;
    case_name = string(out_i.case.case_name);
    wave_model = string(out_i.case.wave_model);
    req_M = out_i.REQ_M;

    q_theory = compute_case_theory_q(cfg_i, KW, case_name, wave_model, req_M);
    theory_rows = [theory_rows; table(case_name, wave_model, req_M, q_theory, ...
        'VariableNames', {'case_name', 'SIM_WaveModel', 'REQ_M', 'q_theory_discrete'})]; %#ok<AGROW>

    T_base = make_eval_base(T_feat, KW);

    T_methods = concat_tables(T_methods, make_adaptive_table(T_base, T_adapt_i));
    T_methods = concat_tables(T_methods, make_q_method_table( ...
        T_base, T_feat.req_mapping, q_theory, "TheoryQDiscrete", cfg_i.f0));

    for qi = 1:numel(fixed_q_values)
        q_fixed = fixed_q_values(qi);
        method_name = "FixedQ_" + strrep(sprintf('%.2f', q_fixed), '.', 'p');
        T_methods = concat_tables(T_methods, make_q_method_table( ...
            T_base, T_feat.req_mapping, q_fixed, method_name, cfg_i.f0));
    end

    T_grad = concat_tables(T_grad, compute_gradient_diagnostics(T_adapt_i, KW, case_name, req_M));
end

T_methods = add_error_metrics(T_methods);

%% Metrics

T_metrics_global = summarize_metrics(T_methods, ...
    {'method_name', 'case_name', 'SIM_WaveModel', 'REQ_M'});
T_metrics_region = summarize_metrics(T_methods, ...
    {'method_name', 'case_name', 'SIM_WaveModel', 'REQ_M', 'region_name'});
T_region_adaptive = T_metrics_region(T_metrics_region.method_name == ...
    "AdaptiveQ_HybridLocalGlobal", :);
T_corr = summarize_gradient_correlations(T_grad);

writetable(T_metrics_global, fullfile(table_dir, ...
    'level12_level07_kwave_adaptive_vs_theory_vs_fixedq_metrics.csv'));
writetable(T_metrics_region, fullfile(table_dir, ...
    'level12_level07_kwave_adaptive_vs_theory_vs_fixedq_by_region.csv'));
writetable(T_region_adaptive, fullfile(table_dir, ...
    'level12_level07_kwave_region_metrics.csv'));
writetable(T_corr, fullfile(table_dir, ...
    'level12_level07_kwave_error_q_gradient_correlations.csv'));
writetable(theory_rows, fullfile(table_dir, ...
    'level12_level07_kwave_theory_quantiles.csv'));

save(fullfile(data_dir, 'level12_level07_kwave_roi_diagnostics.mat'), ...
    'KW', 'T_methods', 'T_metrics_global', 'T_metrics_region', ...
    'T_region_adaptive', 'T_grad', 'T_corr', 'theory_rows', '-v7.3');

%% Figures

plot_method_metric(T_metrics_global, 'MAPE_pct', fig_dir, ...
    'level12_level07_adaptive_vs_theory_vs_fixedq_mape.png', ...
    'Global MAPE by method');
plot_method_metric(T_metrics_global, 'bias_pct', fig_dir, ...
    'level12_level07_adaptive_vs_theory_vs_fixedq_bias.png', ...
    'Global bias by method');
plot_method_metric(T_metrics_global, 'CoV_pct', fig_dir, ...
    'level12_level07_adaptive_vs_theory_vs_fixedq_cov.png', ...
    'Global CoV by method');

plot_region_metric(T_region_adaptive, 'MAPE_pct', fig_dir, ...
    'level12_level07_region_mape_by_M.png', ...
    'AdaptiveQ region MAPE by M');
plot_region_metric(T_region_adaptive, 'bias_pct', fig_dir, ...
    'level12_level07_region_bias_by_M.png', ...
    'AdaptiveQ region bias by M');
plot_region_metric(T_region_adaptive, 'CoV_pct', fig_dir, ...
    'level12_level07_region_cov_by_M.png', ...
    'AdaptiveQ region CoV by M');

plot_adaptive_vs_theory_region(T_metrics_region, 'MAPE_pct', fig_dir, ...
    'level12_level07_region_mape_adaptive_vs_theory.png', ...
    'Region MAPE: AdaptiveQ vs TheoryQDiscrete');
plot_adaptive_vs_theory_region(T_metrics_region, 'bias_pct', fig_dir, ...
    'level12_level07_region_bias_adaptive_vs_theory.png', ...
    'Region bias: AdaptiveQ vs TheoryQDiscrete');

plot_method_map_panels(T_methods, KW, fig_dir);
plot_error_vs_distance(T_grad, fig_dir);
plot_error_vs_q_gradient(T_grad, fig_dir);
plot_q_gradient_maps(T_grad, KW, fig_dir);

%% Console summary

print_console_summary(T_metrics_global, T_metrics_region, T_corr);
fprintf('\nLevel 07 k-Wave ROI diagnostics complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);

%% Local functions

function PATHS = locate_latest_test12_paths(root_dir)

out_root = fullfile(root_dir, 'outputs', 'test_12_cs_guess_window_sweep');
runs = dir(fullfile(out_root, 'test_12_cs_guess_window_sweep_*'));
runs = runs([runs.isdir]);
assert(~isempty(runs), 'No Test 12 run folders found in %s', out_root);
[~, idx] = max([runs.datenum]);
PATHS.run_dir = fullfile(runs(idx).folder, runs(idx).name);
PATHS.analysis_dir = fullfile(PATHS.run_dir, 'analysis');

end

function q = compute_case_theory_q(cfg, KW, case_name, wave_model, req_M)

field_type = theory_field_type(case_name, wave_model);
common = {'M', req_M, 'Gamma', KW.REQ.Gamma, ...
    'PadFactor', KW.REQ.PadFactor, 'Nbins', KW.REQ.Nbins, ...
    'SmoothSigma', KW.REQ.SmoothSigma, 'TheoryMode', 'S2D', ...
    'Plot', false, 'UseDonutFilter', false};

if field_type == "Rev3D"
    q2 = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, KW.REQ.cs_guess, common{:}, ...
        'FieldType', 'Diffuse2D');
    q3 = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, KW.REQ.cs_guess, common{:}, ...
        'FieldType', 'Diffuse3D');
    q = 0.5 * (q2.q_th + q3.q_th);
else
    res = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, KW.REQ.cs_guess, common{:}, ...
        'FieldType', char(field_type));
    q = res.q_th;
end

end

function field_type = theory_field_type(case_name, wave_model)

s = lower(case_name + " " + wave_model);
if contains(s, "directional") || contains(s, "single")
    field_type = "SingleWave";
elseif contains(s, "2d") && contains(s, "diffuse") && ~contains(s, "3d")
    field_type = "Diffuse2D";
elseif contains(s, "rev") || contains(s, "partial") || contains(s, "semi")
    field_type = "Rev3D";
else
    field_type = "Diffuse3D";
end

end

function T = make_eval_base(T_feat, KW)

vars = ["condition_id", "condition_label", "case_name", "SIM_WaveModel", ...
    "SIM_f0", "SIM_cs_bg", "SIM_cs_inc", "REQ_M", "REQ_cs_guess", ...
    "map_iz", "map_ix", "patch_idx", "x_center_m", "z_center_m", "cs_true"];
vars = vars(ismember(vars, string(T_feat.Properties.VariableNames)));
T = T_feat(:, cellstr(vars));
T.region_name = classify_regions(T.x_center_m, T.z_center_m, KW);
T.distance_to_interface_m = abs(hypot(T.x_center_m - KW.inclusion_center_m(1), ...
    T.z_center_m - KW.inclusion_center_m(2)) - KW.inclusion_radius_m);

end

function region = classify_regions(x, z, KW)

r = hypot(x - KW.inclusion_center_m(1), z - KW.inclusion_center_m(2));
band = 0.0015;
core_radius = max(KW.inclusion_radius_m - band, 0);
region = repmat("background_far", numel(x), 1);
region(r <= core_radius) = "inclusion_core";
region(abs(r - KW.inclusion_radius_m) <= band) = "interface_band";

end

function T = make_adaptive_table(T_base, T_adapt)

T = T_base;
T.method_name = repmat("AdaptiveQ_HybridLocalGlobal", height(T), 1);
T.method_family = repmat("adaptive_q", height(T), 1);
T.q_used = T_adapt.q_pred;
T.cs_pred = T_adapt.cs_pred;

end

function T = make_q_method_table(T_base, mappings, q_value, method_name, f0)

T = T_base;
T.method_name = repmat(method_name, height(T), 1);
if startsWith(method_name, "FixedQ")
    T.method_family = repmat("fixed_q", height(T), 1);
else
    T.method_family = repmat("theory_q", height(T), 1);
end
T.q_used = q_value * ones(height(T), 1);
T.cs_pred = q_to_cs_for_table(T.q_used, mappings, f0);

end

function cs_pred = q_to_cs_for_table(q_pred, mappings, f0)

cs_pred = nan(numel(q_pred), 1);
for i = 1:numel(q_pred)
    if isempty(mappings{i}) || ~isfinite(q_pred(i))
        continue;
    end
    cs_pred(i) = adaptive_req.quantile.quantile_to_cs( ...
        mappings{i}, q_pred(i), f0);
end

end

function T = add_error_metrics(T)

T.cs_error = T.cs_pred - T.cs_true;
T.cs_abs_error = abs(T.cs_error);
T.cs_error_pct = 100 * T.cs_error ./ T.cs_true;
T.cs_abs_error_pct = abs(T.cs_error_pct);

end

function Tm = summarize_metrics(T, group_vars)

[G, Tm] = findgroups(T(:, group_vars));
Tm.N = splitapply(@numel, T.cs_error_pct, G);
Tm.MAPE_pct = splitapply(@(x) mean(abs(x), 'omitnan'), T.cs_error_pct, G);
Tm.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.cs_error_pct, G);
Tm.median_error_pct = splitapply(@(x) median(x, 'omitnan'), T.cs_error_pct, G);
Tm.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.cs_error_pct, G);
Tm.P95_abs_error_pct = splitapply(@(x) prctile_finite(abs(x), 95), T.cs_error_pct, G);
Tm.CoV_pct = splitapply(@cov_pct, T.cs_pred, G);
Tm.HighError_gt10_pct = splitapply(@(x) 100 * mean(abs(x) > 10, 'omitnan'), T.cs_error_pct, G);
Tm.HighError_gt20_pct = splitapply(@(x) 100 * mean(abs(x) > 20, 'omitnan'), T.cs_error_pct, G);

end

function T_grad = compute_gradient_diagnostics(T_adapt, KW, case_name, req_M)

T = T_adapt;
[q_map, x_cm, z_cm] = map_from_table(T, 'q_pred');
[cs_map, ~, ~] = map_from_table(T, 'cs_pred');
[err_map, ~, ~] = map_from_table(T, 'cs_abs_error_pct');

dx = median(diff(x_cm)) / 100;
dz = median(diff(z_cm)) / 100;
[dq_dz, dq_dx] = gradient(q_map, dz, dx);
[dc_dz, dc_dx] = gradient(cs_map, dz, dx);
grad_q_mag = hypot(dq_dx, dq_dz);
grad_sws_pred_mag = hypot(dc_dx, dc_dz);

[X, Z] = meshgrid(x_cm / 100, z_cm / 100);
dist = abs(hypot(X - KW.inclusion_center_m(1), ...
    Z - KW.inclusion_center_m(2)) - KW.inclusion_radius_m);

T_grad = table();
T_grad.case_name = repmat(case_name, numel(err_map), 1);
T_grad.REQ_M = req_M * ones(numel(err_map), 1);
T_grad.x_center_m = X(:);
T_grad.z_center_m = Z(:);
T_grad.abs_error_pct = err_map(:);
T_grad.grad_q_mag = grad_q_mag(:);
T_grad.grad_sws_pred_mag = grad_sws_pred_mag(:);
T_grad.distance_to_interface_m = dist(:);
T_grad.q_pred = q_map(:);
T_grad.cs_pred = cs_map(:);
T_grad.region_name = classify_regions(T_grad.x_center_m, T_grad.z_center_m, KW);

end

function Tc = summarize_gradient_correlations(T)

[G, Tc] = findgroups(T(:, {'case_name', 'REQ_M', 'region_name'}));
Tc.N = splitapply(@numel, T.abs_error_pct, G);
Tc.rho_abs_error_vs_grad_q = splitapply(@spearman_pair, T.abs_error_pct, T.grad_q_mag, G);
Tc.rho_abs_error_vs_distance = splitapply(@spearman_pair, T.abs_error_pct, T.distance_to_interface_m, G);
Tc.rho_abs_error_vs_grad_sws = splitapply(@spearman_pair, T.abs_error_pct, T.grad_sws_pred_mag, G);

end

function rho = spearman_pair(x, y)

valid = isfinite(x) & isfinite(y);
if nnz(valid) < 4
    rho = NaN;
else
    rho = corr(x(valid), y(valid), 'Type', 'Spearman');
end

end

function y = cov_pct(x)

x = x(isfinite(x));
if isempty(x) || abs(mean(x)) < eps
    y = NaN;
else
    y = 100 * std(x) / mean(x);
end

end

function y = prctile_finite(x, p)

x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = prctile(x, p);
end

end

function [A, x_cm, z_cm] = map_from_table(T, value_var)

value_var = char(value_var);
nz = max(T.map_iz);
nx = max(T.map_ix);
A = nan(nz, nx);
idx = sub2ind([nz nx], T.map_iz, T.map_ix);
A(idx) = T.(value_var);
x_cm = unique(T.x_center_m, 'stable') * 100;
z_cm = unique(T.z_center_m, 'stable') * 100;

end

function plot_method_metric(T, metric_var, fig_dir, file_name, title_text)

method_order = ["AdaptiveQ_HybridLocalGlobal", "TheoryQDiscrete", ...
    "FixedQ_0p50", "FixedQ_0p60", "FixedQ_0p70", "FixedQ_0p80", "FixedQ_0p90"];
T = T(ismember(T.method_name, method_order), :);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 32 17]);
tiledlayout(2, ceil(numel(unique(T.case_name, 'stable')) / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
cases = unique(T.case_name, 'stable');
for ci = 1:numel(cases)
    ax = nexttile;
    Ti = T(T.case_name == cases(ci), :);
    M_values = unique(Ti.REQ_M, 'stable');
    Y = nan(numel(M_values), numel(method_order));
    for mi = 1:numel(M_values)
        for mj = 1:numel(method_order)
            idx = Ti.REQ_M == M_values(mi) & Ti.method_name == method_order(mj);
            if any(idx)
                Y(mi, mj) = Ti.(metric_var)(find(idx, 1, 'first'));
            end
        end
    end
    bar(ax, categorical(string(M_values)), Y);
    title(ax, cases(ci), 'Interpreter', 'none', 'FontWeight', 'normal');
    ylabel(ax, metric_label(metric_var));
    xlabel(ax, 'REQ M');
    grid(ax, 'on');
    set(ax, 'FontSize', 9);
end
lgd = legend(method_short_names(method_order), 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Interpreter', 'none', 'NumColumns', 4);
lgd.FontSize = 8;
sgtitle(title_text, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, file_name), 'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_region_metric(T, metric_var, fig_dir, file_name, title_text)

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 28 16]);
tiledlayout(2, ceil(numel(unique(T.case_name, 'stable')) / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
cases = unique(T.case_name, 'stable');
region_order = ["background_far", "inclusion_core", "interface_band"];
for ci = 1:numel(cases)
    ax = nexttile;
    Ti = T(T.case_name == cases(ci), :);
    M_values = unique(Ti.REQ_M, 'stable');
    Y = nan(numel(M_values), numel(region_order));
    for mi = 1:numel(M_values)
        for ri = 1:numel(region_order)
            idx = Ti.REQ_M == M_values(mi) & Ti.region_name == region_order(ri);
            if any(idx)
                Y(mi, ri) = Ti.(metric_var)(find(idx, 1, 'first'));
            end
        end
    end
    bar(ax, categorical(string(M_values)), Y);
    title(ax, cases(ci), 'Interpreter', 'none', 'FontWeight', 'normal');
    ylabel(ax, metric_label(metric_var));
    xlabel(ax, 'REQ M');
    grid(ax, 'on');
    set(ax, 'FontSize', 9);
end
lgd = legend(strrep(cellstr(region_order), '_', ' '), 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 3);
lgd.FontSize = 9;
sgtitle(title_text, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, file_name), 'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_adaptive_vs_theory_region(T, metric_var, fig_dir, file_name, title_text)

T = T(ismember(T.method_name, ["AdaptiveQ_HybridLocalGlobal", "TheoryQDiscrete"]), :);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 32 18]);
tiledlayout(2, ceil(numel(unique(T.case_name, 'stable')) / 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
cases = unique(T.case_name, 'stable');
region_order = ["background_far", "inclusion_core", "interface_band"];
method_order = ["AdaptiveQ_HybridLocalGlobal", "TheoryQDiscrete"];
for ci = 1:numel(cases)
    ax = nexttile;
    Ti = T(T.case_name == cases(ci), :);
    labels = strings(0, 1);
    Y = [];
    for mi = unique(Ti.REQ_M, 'stable')'
        for ri = 1:numel(region_order)
            labels(end+1, 1) = "M=" + string(mi) + " " + strrep(region_order(ri), "_", " "); %#ok<AGROW>
            row = nan(1, numel(method_order));
            for mj = 1:numel(method_order)
                idx = Ti.REQ_M == mi & Ti.region_name == region_order(ri) & ...
                    Ti.method_name == method_order(mj);
                if any(idx)
                    row(mj) = Ti.(metric_var)(find(idx, 1, 'first'));
                end
            end
            Y = [Y; row]; %#ok<AGROW>
        end
    end
    bar(ax, categorical(labels), Y);
    title(ax, cases(ci), 'Interpreter', 'none', 'FontWeight', 'normal');
    ylabel(ax, metric_label(metric_var));
    grid(ax, 'on');
    xtickangle(ax, 35);
    set(ax, 'FontSize', 8);
end
lgd = legend({'AdaptiveQ', 'TheoryQDiscrete'}, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 2);
lgd.FontSize = 9;
sgtitle(title_text, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, file_name), 'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_method_map_panels(T, KW, fig_dir)

method_order = ["AdaptiveQ_HybridLocalGlobal", "TheoryQDiscrete", ...
    "FixedQ_0p50", "FixedQ_0p60", "FixedQ_0p70", "FixedQ_0p80", "FixedQ_0p90"];
for method_i = method_order
    Ti = T(T.method_name == method_i, :);
    if isempty(Ti)
        continue;
    end
    safe_name = matlab.lang.makeValidName(char(method_i));
    plot_one_method_map(Ti, KW, fig_dir, method_i, safe_name, ...
        'cs_pred', 'SWS (m/s)', 'turbo', ...
        "level12_level07_sws_maps_" + string(safe_name) + ".png");
    plot_one_method_map(Ti, KW, fig_dir, method_i, safe_name, ...
        'cs_error_pct', 'signed SWS error (%)', 'parula', ...
        "level12_level07_error_maps_" + string(safe_name) + ".png");
    plot_one_method_map(Ti, KW, fig_dir, method_i, safe_name, ...
        'q_used', 'q used', 'turbo', ...
        "level12_level07_q_maps_" + string(safe_name) + ".png");
end

end

function plot_one_method_map(T, KW, fig_dir, method_name, ~, value_var, value_label, cmap_name, file_name)

value_var = string(value_var);
cmap_name = string(cmap_name);
cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 7.0*numel(cases) 5.6*numel(M_values)]);
tl = tiledlayout(numel(M_values), numel(cases), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Tc = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile(tl);
        if isempty(Tc)
            axis(ax, 'off');
            continue;
        end
        [A, x_cm, z_cm] = map_from_table(Tc, value_var);
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontSize', 8);
        if cmap_name == "turbo"
            colormap(ax, turbo);
        else
            colormap(ax, parula);
        end
        colorbar(ax);
        draw_interface(ax, KW);
        if value_var == "cs_pred"
            stat_text = sprintf('MAPE %.2f%%', mean(Tc.cs_abs_error_pct, 'omitnan'));
        elseif value_var == "cs_error_pct"
            stat_text = sprintf('bias %.2f%%', mean(Tc.cs_error_pct, 'omitnan'));
        else
            stat_text = sprintf('median q %.3f', median(Tc.q_used, 'omitnan'));
        end
        title(ax, sprintf('%s | M=%g | %s', cases(ci), M_values(mi), stat_text), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8);
        xlabel(ax, 'x (cm)');
        ylabel(ax, 'z (cm)');
    end
end

title(tl, sprintf('%s: %s', char(method_short_names(method_name)), value_label), ...
    'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 13);
exportgraphics(fig, fullfile(fig_dir, char(file_name)), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function out = method_short_names(methods)

methods = string(methods);
out = methods;
out(methods == "AdaptiveQ_HybridLocalGlobal") = "AdaptiveQ";
out(methods == "TheoryQDiscrete") = "TheoryQ";
out(methods == "FixedQ_0p50") = "Fixed q=0.50";
out(methods == "FixedQ_0p60") = "Fixed q=0.60";
out(methods == "FixedQ_0p70") = "Fixed q=0.70";
out(methods == "FixedQ_0p80") = "Fixed q=0.80";
out(methods == "FixedQ_0p90") = "Fixed q=0.90";

end

function out = metric_label(metric_var)

switch string(metric_var)
    case "MAPE_pct"
        out = 'MAPE (%)';
    case "bias_pct"
        out = 'bias (%)';
    case "CoV_pct"
        out = 'CoV (%)';
    otherwise
        out = strrep(char(metric_var), '_', '\_');
end

end

function plot_error_vs_distance(T, fig_dir)

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 11]);
scatter(T.distance_to_interface_m * 1000, T.abs_error_pct, 8, double(T.REQ_M), 'filled', 'MarkerFaceAlpha', 0.25);
grid on; colorbar;
xlabel('distance to interface (mm)');
ylabel('|SWS error| (%)');
title('AdaptiveQ error vs distance to interface', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, 'level12_level07_error_vs_distance_to_interface.png'), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_error_vs_q_gradient(T, fig_dir)

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 11]);
scatter(T.grad_q_mag, T.abs_error_pct, 8, double(T.REQ_M), 'filled', 'MarkerFaceAlpha', 0.25);
grid on; colorbar;
xlabel('|grad q|');
ylabel('|SWS error| (%)');
title('AdaptiveQ error vs q-map gradient', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, 'level12_level07_error_vs_q_gradient.png'), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function plot_q_gradient_maps(T, KW, fig_dir)

cases = unique(T.case_name, 'stable');
M_values = unique(T.REQ_M, 'stable');
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8*numel(cases) 5.8*numel(M_values)]);
tiledlayout(numel(M_values), numel(cases), 'TileSpacing', 'compact', 'Padding', 'compact');
for mi = 1:numel(M_values)
    for ci = 1:numel(cases)
        Ti = T(T.case_name == cases(ci) & T.REQ_M == M_values(mi), :);
        ax = nexttile;
        [A, x_cm, z_cm] = scattered_map(Ti, 'grad_q_mag');
        imagesc(ax, x_cm, z_cm, A);
        axis(ax, 'image'); set(ax, 'YDir', 'normal');
        colorbar(ax); colormap(ax, turbo);
        draw_interface(ax, KW);
        title(ax, sprintf('%s | M=%g', cases(ci), M_values(mi)), ...
            'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
    end
end
sgtitle('|grad q| maps', 'FontWeight', 'normal');
exportgraphics(fig, fullfile(fig_dir, 'level12_level07_q_gradient_maps.png'), ...
    'Resolution', 250, 'BackgroundColor', 'white');
close(fig);

end

function [A, x_cm, z_cm] = scattered_map(T, value_var)

value_var = char(value_var);
x_vals = unique(T.x_center_m, 'stable');
z_vals = unique(T.z_center_m, 'stable');
[~, ix] = ismember(T.x_center_m, x_vals);
[~, iz] = ismember(T.z_center_m, z_vals);
A = nan(numel(z_vals), numel(x_vals));
A(sub2ind(size(A), iz, ix)) = T.(value_var);
x_cm = x_vals * 100;
z_cm = z_vals * 100;

end

function draw_interface(ax, KW)

hold(ax, 'on');
th = linspace(0, 2*pi, 300);
plot(ax, KW.inclusion_center_m(1)*100 + KW.inclusion_radius_m*100*cos(th), ...
    KW.inclusion_center_m(2)*100 + KW.inclusion_radius_m*100*sin(th), ...
    'w--', 'LineWidth', 1.1);
hold(ax, 'off');

end

function print_console_summary(Tg, Tr, Tc)

[~, i_best] = min(Tg.MAPE_pct);
fprintf('\nBest global method by MAPE:\n');
disp(Tg(i_best, {'method_name', 'case_name', 'REQ_M', 'MAPE_pct', 'bias_pct'}));

fprintf('\nBest method by wavefield type:\n');
cases = unique(Tg.case_name, 'stable');
for ci = 1:numel(cases)
    Ti = Tg(Tg.case_name == cases(ci), :);
    [~, ii] = min(Ti.MAPE_pct);
    disp(Ti(ii, {'case_name', 'method_name', 'REQ_M', 'MAPE_pct'}));
end

fprintf('\nBest method by region:\n');
regions = unique(Tr.region_name, 'stable');
for ri = 1:numel(regions)
    Ti = Tr(Tr.region_name == regions(ri), :);
    [~, ii] = min(Ti.MAPE_pct);
    disp(Ti(ii, {'region_name', 'case_name', 'method_name', 'REQ_M', 'MAPE_pct', 'bias_pct'}));
end

Ta = Tg(Tg.method_name == "AdaptiveQ_HybridLocalGlobal", :);
Tt = Tg(Tg.method_name == "TheoryQDiscrete", :);
fprintf('\nAdaptiveQ mean MAPE = %.3f%%; TheoryQDiscrete mean MAPE = %.3f%%.\n', ...
    mean(Ta.MAPE_pct, 'omitnan'), mean(Tt.MAPE_pct, 'omitnan'));

Tf = Tg(startsWith(Tg.method_name, "FixedQ"), :);
best_fixed = min(Tf.MAPE_pct);
fprintf('Best fixed-q global MAPE = %.3f%%; TheoryQDiscrete mean MAPE = %.3f%%.\n', ...
    best_fixed, mean(Tt.MAPE_pct, 'omitnan'));

Tincl = Tr(Tr.method_name == "AdaptiveQ_HybridLocalGlobal" & ...
    Tr.region_name == "inclusion_core", :);
fprintf('AdaptiveQ inclusion_core mean bias = %.3f%%.\n', ...
    mean(Tincl.bias_pct, 'omitnan'));

for reg = ["background_far", "inclusion_core", "interface_band"]
    Ti = Tr(Tr.method_name == "AdaptiveQ_HybridLocalGlobal" & Tr.region_name == reg, :);
    [Gm, Tm] = findgroups(Ti(:, {'REQ_M'}));
    Tm.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), Ti.MAPE_pct, Gm);
    [~, ii_best] = min(Tm.MAPE_pct);
    [~, ii_worst] = max(Tm.MAPE_pct);
    fprintf('%s best mean M=%g (MAPE %.2f%%); worst mean M=%g (MAPE %.2f%%).\n', ...
        reg, Tm.REQ_M(ii_best), Tm.MAPE_pct(ii_best), ...
        Tm.REQ_M(ii_worst), Tm.MAPE_pct(ii_worst));
end

fprintf('Mean Spearman rho |error| vs |grad q| = %.3f.\n', ...
    mean(Tc.rho_abs_error_vs_grad_q, 'omitnan'));

end

function T = concat_tables(A, B)

if isempty(A), T = B; return; end
if isempty(B), T = A; return; end
vars_all = unique([string(A.Properties.VariableNames), string(B.Properties.VariableNames)], 'stable');
A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);
T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];

end

function T = add_missing_columns(T, vars_all)

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    if ismember(vars_all(i), vars), continue; end
    name_i = char(vars_all(i));
    if any(endsWith(vars_all(i), ["name", "label", "type", "role", "family", "region"]))
        T.(name_i) = strings(height(T), 1);
    else
        T.(name_i) = nan(height(T), 1);
    end
end

end

function make_dir_if_needed(path_i)

if ~exist(path_i, 'dir')
    mkdir(path_i);
end

end
