%% analyze_test_35_extra_maps_and_summary.m
% Supplemental maps and summaries for Test 35.
%
% This script does not train or recompute REQ. It reads Test 35 patch-level
% CSVs and exports map panels similar in spirit to Test 34. Missing sampled
% pixels are shown with a light background instead of being painted as the
% lowest colormap value, which avoids confusing sparse-sampling holes with
% physical/Theory artifacts.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST35_EXTRA_MAPS = representative | all
%   ADAPTIVE_REQ_TEST35_EXTRA_KWAVE = true | false

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.08);

CFG = struct();
CFG.MapMode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST35_EXTRA_MAPS'))));
if CFG.MapMode == "", CFG.MapMode = "representative"; end
assert(ismember(CFG.MapMode, ["representative","all"]), ...
    'ADAPTIVE_REQ_TEST35_EXTRA_MAPS must be representative or all.');
CFG.PlotKWave = env_true('ADAPTIVE_REQ_TEST35_EXTRA_KWAVE', true);
CFG.Models = ["theory_discrete", "q_spectrum_only", ...
    "q_spectrum_plus_theory_composition", ...
    "delta_q_theory_composition", "delta_logk_theory_composition"];
CFG.MainModel = "delta_logk_theory_composition";

OUT.root_dir = fullfile(root_dir,'outputs','test_35_spectral_composition_to_q_model');
OUT.table_dir = fullfile(OUT.root_dir,'tables');
OUT.figure_dir = fullfile(OUT.root_dir,'figures','extra_maps');
OUT.synthetic_dir = fullfile(OUT.figure_dir,'synthetic_by_condition');
OUT.kwave_dir = fullfile(OUT.figure_dir,'kwave_by_case');
OUT.summary_dir = fullfile(OUT.figure_dir,'summary');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end

SYN_FILE = fullfile(OUT.table_dir,'test35_synthetic_patch_level_predictions.csv');
KW_FILE = fullfile(OUT.table_dir,'test35_kwave_patch_level_predictions.csv');
assert(exist(SYN_FILE,'file') == 2, 'Missing %s', SYN_FILE);

fprintf('\nTest 35 supplemental maps\n');
fprintf('Map mode: %s | kWave: %d\n', CFG.MapMode, CFG.PlotKWave);
fprintf('Reading synthetic patch table...\n');
T = readtable(SYN_FILE, 'TextType','string');
T = T(ismember(T.model_name, CFG.Models), :);
fprintf('Synthetic rows after model filter: %d\n', height(T));

plot_synthetic_maps(T, CFG, OUT);
plot_summary_heatmaps(OUT);
explain_sparse_maps(T, OUT);

if CFG.PlotKWave && exist(KW_FILE,'file') == 2
    fprintf('Reading kWave patch table...\n');
    K = readtable(KW_FILE, 'TextType','string');
    K = K(ismember(K.model_name, CFG.Models), :);
    fprintf('kWave rows after model filter: %d\n', height(K));
    plot_kwave_maps(K, CFG, OUT);
end

fprintf('\nExtra maps saved under:\n%s\n', OUT.figure_dir);
fprintf('Test 35 supplemental maps complete.\n');

%% Map plotting

function plot_synthetic_maps(T, CFG, OUT)
keys = unique(T.condition_key, 'stable');
if CFG.MapMode == "representative"
    keys = select_representative_synthetic_keys(T, keys);
end
fprintf('Saving %d synthetic condition map panels.\n', numel(keys));
for i = 1:numel(keys)
    X = T(T.condition_key == keys(i), :);
    if isempty(X), continue; end
    row = X(1,:);
    out_dir = fullfile(OUT.synthetic_dir, sanitize(row.geometry), ...
        sanitize(row.field_regime), "M" + string(row.M), ...
        "f" + string(row.f0));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    out_file = fullfile(out_dir, "test35_extra__" + sanitize(keys(i)) + ".png");
    plot_one_synthetic_condition(X, CFG, out_file);
    if mod(i, 20) == 0 || i == numel(keys)
        fprintf('  synthetic maps %d/%d\n', i, numel(keys));
    end
end
end

function keys = select_representative_synthetic_keys(T, all_keys)
want = [
    "homogeneous_cs2", 300, "directional_2D", 2
    "homogeneous_cs2", 600, "diffuse_3D", 4
    "homogeneous_cs3", 300, "directional_2D", 2
    "homogeneous_cs3", 600, "diffuse_3D", 4
    "bilayer_2_3", 300, "directional_2D", 2
    "bilayer_2_3", 500, "diffuse_2D", 3
    "bilayer_2_3", 600, "diffuse_3D", 4
    "circular_inclusion_2_3", 300, "directional_2D", 2
    "circular_inclusion_2_3", 500, "diffuse_2D", 3
    "circular_inclusion_2_3", 600, "diffuse_3D", 4
    ];
keys = strings(0,1);
base = T(T.model_name == T.model_name(1), :);
for i = 1:size(want,1)
    idx = base.geometry == string(want(i,1)) & ...
        base.f0 == double(want(i,2)) & ...
        base.field_regime == string(want(i,3)) & ...
        base.M == double(want(i,4));
    k = unique(base.condition_key(idx), 'stable');
    if ~isempty(k), keys(end+1,1) = k(1); end %#ok<AGROW>
end
if isempty(keys)
    keys = all_keys(1:min(12,numel(all_keys)));
end
keys = unique(keys, 'stable');
end

function plot_one_synthetic_condition(X, CFG, out_file)
main = X(X.model_name == CFG.MainModel, :);
if isempty(main), main = X(X.model_name == X.model_name(1), :); end
[true_map, nz, nx] = grid_from_rows(main, main.true_SWS);
purity_map = grid_from_rows(main, main.predicted_patch_purity, nz, nx);
mixed_map = grid_from_rows(main, main.p_mixed, nz, nx);
sample_map = isfinite(true_map);

fig = figure('Color','w','Units','centimeters','Position',[1 1 30 18]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), true_map, 'True SWS');
for model = CFG.Models
    M = X(X.model_name == model, :);
    if isempty(M), continue; end
    sws = grid_from_rows(M, M.sws_pred, nz, nx);
    plot_map(nexttile(tl), sws, pretty_model(model));
end
err = grid_from_rows(main, main.sws_abs_error_pct, nz, nx);
plot_map(nexttile(tl), err, "Abs error " + pretty_model(CFG.MainModel) + " (%)");
plot_map(nexttile(tl), purity_map, 'Predicted purity');
plot_map(nexttile(tl), mixed_map, 'Predicted p(mixed)');
plot_map(nexttile(tl), double(sample_map), 'Sampled-pixel mask');
title(tl, main.condition_key(1), 'Interpreter','none');
export_fig(fig, out_file);
end

function plot_kwave_maps(K, CFG, OUT)
keys = unique(K.condition_key, 'stable');
if CFG.MapMode == "representative"
    keys = keys(1:min(12,numel(keys)));
end
fprintf('Saving %d kWave map panels.\n', numel(keys));
for i = 1:numel(keys)
    X = K(K.condition_key == keys(i), :);
    if isempty(X), continue; end
    row = X(1,:);
    out_dir = fullfile(OUT.kwave_dir, sanitize(row.case_name), "M" + string(row.M));
    if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end
    out_file = fullfile(out_dir, "test35_kwave_extra__" + sanitize(keys(i)) + ".png");
    plot_one_kwave_condition(X, CFG, out_file);
end
end

function plot_one_kwave_condition(X, CFG, out_file)
main = X(X.model_name == CFG.MainModel, :);
if isempty(main), main = X(X.model_name == X.model_name(1), :); end
[true_map, nz, nx] = grid_from_rows(main, main.true_SWS);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 16]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
plot_map(nexttile(tl), true_map, 'True SWS');
for model = CFG.Models
    M = X(X.model_name == model, :);
    if isempty(M), continue; end
    sws = grid_from_rows(M, M.sws_pred, nz, nx);
    plot_map(nexttile(tl), sws, pretty_model(model));
end
plot_map(nexttile(tl), grid_from_rows(main, main.sws_abs_error_pct, nz, nx), ...
    "Abs error " + pretty_model(CFG.MainModel) + " (%)");
if ismember('predicted_patch_purity', main.Properties.VariableNames)
    plot_map(nexttile(tl), grid_from_rows(main, main.predicted_patch_purity, nz, nx), ...
        'Predicted purity');
end
if ismember('p_mixed', main.Properties.VariableNames)
    plot_map(nexttile(tl), grid_from_rows(main, main.p_mixed, nz, nx), ...
        'Predicted p(mixed)');
end
plot_map(nexttile(tl), double(isfinite(true_map)), 'Sampled-pixel mask');
title(tl, main.condition_key(1), 'Interpreter','none');
export_fig(fig, out_file);
end

%% Summary figures

function plot_summary_heatmaps(OUT)
summary_file = fullfile(OUT.table_dir,'test35_q_model_summary_synthetic_heldout.csv');
geom_file = fullfile(OUT.table_dir,'test35_q_model_summary_by_geometry.csv');
regime_file = fullfile(OUT.table_dir,'test35_q_model_summary_by_regime.csv');
M_file = fullfile(OUT.table_dir,'test35_q_model_summary_by_M.csv');
purity_file = fullfile(OUT.table_dir,'test35_q_model_summary_by_purity_bin.csv');
if exist(summary_file,'file') ~= 2, return; end
S = readtable(summary_file, 'TextType','string');
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 10]);
bar(categorical(S.model_name), S.MAPE_pct); xtickangle(25); grid on;
ylabel('Held-out synthetic MAPE (%)');
title('Test 35 held-out synthetic model ranking', 'FontWeight','normal');
export_fig(fig, fullfile(OUT.summary_dir,'test35_extra_model_ranking.png'));

plot_heatmap_csv(geom_file, 'geometry', OUT, 'test35_extra_mape_by_geometry.png');
plot_heatmap_csv(regime_file, 'field_regime', OUT, 'test35_extra_mape_by_regime.png');
plot_heatmap_csv(M_file, 'M', OUT, 'test35_extra_mape_by_M.png');
plot_heatmap_csv(purity_file, 'purity_bin', OUT, 'test35_extra_mape_by_purity.png');
end

function plot_heatmap_csv(file, col, OUT, name)
if exist(file,'file') ~= 2, return; end
T = readtable(file, 'TextType','string');
models = unique(T.model_name, 'stable');
groups = unique(string(T.(col)), 'stable');
C = nan(numel(models), numel(groups));
for mi = 1:numel(models)
    for gi = 1:numel(groups)
        idx = T.model_name == models(mi) & string(T.(col)) == groups(gi);
        if any(idx), C(mi,gi) = mean(T.MAPE_pct(idx), 'omitnan'); end
    end
end
fig = figure('Color','w','Units','centimeters','Position',[2 2 22 11]);
heatmap(groups, models, C, 'ColorbarVisible','on');
title("MAPE (%) by " + string(col));
export_fig(fig, fullfile(OUT.summary_dir, name));
end

function explain_sparse_maps(T, OUT)
base = T(T.model_name == T.model_name(1), :);
[G, keys] = findgroups(base.condition_key);
counts = splitapply(@numel, base.true_SWS, G);
txt = fullfile(OUT.summary_dir,'test35_extra_map_notes.txt');
fid = fopen(txt,'w'); assert(fid > 0);
fprintf(fid, ['Test 35 maps are drawn from sampled patch rows, not from every ', ...
    'possible REQ window. Unsampled grid positions are shown as light/blank ', ...
    'pixels in the extra maps. In the original Test 35 representative maps, ', ...
    'those NaNs could appear as dark low-colormap dots, especially in True SWS ', ...
    'or Theory panels. Those dots are sampling/plotting holes, not physical ', ...
    '2 m/s inclusions inside homogeneous maps and not isolated Theory failures.\n\n']);
fprintf(fid, 'Rows per condition in the plotted table range from %d to %d.\n', ...
    min(counts), max(counts));
fprintf(fid, 'Number of synthetic conditions in table: %d.\n', numel(keys));
fclose(fid);
end

%% Helpers

function [Z,nz,nx] = grid_from_rows(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
idx = isfinite(T.map_iz) & isfinite(T.map_ix);
Z(sub2ind([nz,nx], T.map_iz(idx), T.map_ix(idx))) = values(idx);
end

function plot_map(ax, Z, ttl)
imagesc(ax, Z, 'AlphaData', isfinite(Z));
axis(ax,'image'); ax.Color = [0.94 0.94 0.94];
colorbar(ax); title(ax, ttl, 'Interpreter','none', 'FontWeight','normal');
set(ax, 'XTick', [], 'YTick', []);
end

function s = pretty_model(s)
s = string(s);
s = replace(s, "q_spectrum_plus_theory_composition", "q + theory + composition");
s = replace(s, "delta_logk_theory_composition", "delta log-k + theory + composition");
s = replace(s, "delta_q_theory_composition", "delta q + theory + composition");
s = replace(s, "q_spectrum_only", "q spectrum only");
s = replace(s, "theory_discrete", "TheoryQDiscrete");
s = replace(s, "_", " ");
end

function tf = env_true(name, default)
raw = strtrim(lower(string(getenv(name))));
if raw == "", tf = logical(default); return; end
tf = ismember(raw, ["1","true","yes","on"]);
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
