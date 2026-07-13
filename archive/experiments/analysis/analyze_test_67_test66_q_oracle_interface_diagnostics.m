%% analyze_test_67_test66_q_oracle_interface_diagnostics.m
% Test 67: q-oracle/interface diagnostics for Test 66.
%
% This script is intentionally diagnostic only. It does not retrain any q,
% composition, confidence, or correction model. It reads existing Test 66
% patch tables and REQ caches, then compares the frozen q prediction against
% the local q-oracle implied by the REQ mapping.
%
% Main question:
%   Is the low-SWS band around interfaces caused by q-model domain shift, by
%   the finite local REQ mapping itself, or by mixed/interface patches and
%   multi-source field complexity?
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST67_TEST66_SOURCE = full_a | quick | /path/to/test66/run
%   ADAPTIVE_REQ_TEST67_INCLUDE_READOUT_MEDIUM = true | false
%   ADAPTIVE_REQ_TEST67_MAX_MAP_CONDITIONS = integer
%   ADAPTIVE_REQ_TEST67_MAX_SCATTER_POINTS = integer

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir);

fprintf('\nTest 67: Test66 q-oracle interface diagnostics\n');
fprintf('Source Test66 run: %s\n', CFG.Test66RunDir);
fprintf('No model is trained and Test66 is not modified.\n');

T = load_selected_test66_rows(CFG);
assert(~isempty(T), 'No matching Test66 rows were found for Test67 filters.');
fprintf('Selected %d patch/model rows from Test66.\n', height(T));

T = attach_oracle_sws_from_req_cache(T, CFG);
T = add_diagnostic_variables(T);

writetable(T, fullfile(OUT.table_dir, 'test67_selected_patch_level_results.csv'));

SUM = struct();
SUM.full = summarize_metrics(T, ["geometry","f0","field_regime","source_seed","M", ...
    "realism_level","roi_region","distance_bin_test67","purity_bin_test67"]);
SUM.field = summarize_metrics(T, ["geometry","f0","M","field_regime"]);
SUM.purity = summarize_metrics(T, ["geometry","f0","field_regime","purity_bin_test67"]);
SUM.distance = summarize_metrics(T, ["geometry","f0","field_regime","distance_over_window_bin_test67"]);
SUM.roi = summarize_metrics(T, ["geometry","f0","field_regime","roi_region"]);

writetable(SUM.full, fullfile(OUT.table_dir, 'summary_q_pred_vs_q_oracle.csv'));
writetable(SUM.field, fullfile(OUT.table_dir, 'summary_by_field_regime.csv'));
writetable(SUM.purity, fullfile(OUT.table_dir, 'summary_by_patch_purity.csv'));
writetable(SUM.distance, fullfile(OUT.table_dir, 'summary_by_distance_over_window.csv'));
writetable(SUM.roi, fullfile(OUT.table_dir, 'summary_by_roi.csv'));

fprintf('Writing diagnostic figures...\n');
make_representative_maps(T, OUT, CFG);
make_frequency_grids(T, OUT, CFG);
make_directional_vs_diffuse_figures(T, OUT, CFG);
make_error_vs_distance_figures(T, OUT, CFG);
make_error_vs_purity_figures(T, OUT, CFG);
make_error_vs_q_gradient_figures(T, OUT, CFG);
make_q_maps_and_histograms(T, OUT, CFG);
make_local_spectra(T, OUT, CFG);
write_readme(T, SUM, OUT, CFG);

save(fullfile(OUT.data_dir, 'test67_q_oracle_interface_diagnostics.mat'), ...
    'CFG','SUM','-v7.3');

fprintf('\nTest 67 complete.\nTables: %s\nFigures: %s\nREADME: %s\n', ...
    OUT.table_dir, OUT.figure_dir, fullfile(OUT.root_dir, 'README_results.md'));

%% Configuration

function CFG = default_config(root_dir)
CFG = struct();
CFG.Geometries = ["inclusion_2_3","inclusion_2_4","bilayer_2_3","bilayer_inclusion_2_3_4"];
CFG.Realism = "clean";
if env_true('ADAPTIVE_REQ_TEST67_INCLUDE_READOUT_MEDIUM', false)
    CFG.Realism = ["clean","readout_medium"];
end
CFG.FieldRegimes = ["single_source_lateral","diffuse_like_8src"];
CFG.Frequencies = [200 300 400 500 600];
CFG.M = 2;
CFG.SourceSeeds = [1 2];
CFG.Models = ["q_spectrum_plus_composition","q_spectrum_only"];
CFG.MainModel = "q_spectrum_plus_composition";
CFG.SecondaryModel = "q_spectrum_only";
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST67_MAX_MAP_CONDITIONS', 16);
CFG.MaxScatterPoints = env_number('ADAPTIVE_REQ_TEST67_MAX_SCATTER_POINTS', 120000);

src = env_string('ADAPTIVE_REQ_TEST67_TEST66_SOURCE', "full_a");
base = fullfile(root_dir, 'outputs', 'test_66_eikonal_realistic_transfer_validation');
if isfolder(src)
    CFG.Test66RunDir = char(src);
else
    CFG.Test66RunDir = fullfile(base, char(src));
end
if ~isfolder(CFG.Test66RunDir) && strcmpi(char(src), 'full')
    CFG.Test66RunDir = fullfile(base, 'full_a');
end
CFG.PatchCsv = fullfile(CFG.Test66RunDir, 'tables', 'test66_patch_level_results.csv');
CFG.ReqCacheDir = fullfile(CFG.Test66RunDir, 'data', 'req_cache');
assert(exist(CFG.PatchCsv,'file') == 2, 'Missing Test66 patch table: %s', CFG.PatchCsv);
assert(exist(CFG.ReqCacheDir,'dir') == 7, 'Missing Test66 REQ cache folder: %s', CFG.ReqCacheDir);
end

function OUT = make_output_dirs(root_dir)
OUT = struct();
OUT.root_dir = fullfile(root_dir, 'outputs', 'test_67_q_oracle_interface_diagnostics');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_q_oracle_diagnostics');
OUT.freq_grid_dir = fullfile(OUT.figure_dir, 'frequency_grids');
OUT.field_compare_dir = fullfile(OUT.figure_dir, 'directional_vs_diffuse_clean');
OUT.dist_dir = fullfile(OUT.figure_dir, 'error_vs_interface_distance');
OUT.purity_dir = fullfile(OUT.figure_dir, 'error_vs_patch_purity');
OUT.qgrad_dir = fullfile(OUT.figure_dir, 'error_vs_q_gradient');
OUT.q_dir = fullfile(OUT.figure_dir, 'q_maps_and_histograms');
OUT.spectra_dir = fullfile(OUT.figure_dir, 'local_spectra');

dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if ~exist(dirs(i), 'dir'), mkdir(dirs(i)); end
end
end

%% Loading and diagnostics

function T = load_selected_test66_rows(CFG)
keep = ["map_iz","map_ix","cx","cz","x_center_m","z_center_m","condition_key", ...
    "geometry","realism_level","field_regime","source_seed","noise_seed","f0","M","REQ_M", ...
    "true_SWS","k_true","patch_purity","q_oracle","distance_to_interface_m", ...
    "distance_to_interface_mm","distance_to_interface_over_window_radius", ...
    "purity_bin","distance_bin","distance_over_window_bin","roi_region", ...
    "q_theory_prior","sws_theory","model_name","q_pred","k_pred","SWS_pred", ...
    "sws_signed_error_pct","sws_abs_error_pct","high_error20"];

ds = tabularTextDatastore(CFG.PatchCsv, 'Delimiter', ',');
vars = string(ds.VariableNames);
ds.SelectedVariableNames = cellstr(intersect(keep, vars, 'stable'));
ds.ReadSize = 200000;

parts = {};
while hasdata(ds)
    C = read(ds);
    C = normalize_loaded_types(C);
    m = ismember(C.geometry, CFG.Geometries) & ...
        ismember(C.realism_level, CFG.Realism) & ...
        ismember(C.field_regime, CFG.FieldRegimes) & ...
        ismember(C.f0, CFG.Frequencies) & ...
        ismember(round(C.REQ_M), CFG.M) & ...
        ismember(C.source_seed, CFG.SourceSeeds) & ...
        ismember(C.model_name, CFG.Models);
    if any(m)
        parts{end+1,1} = C(m,:); %#ok<AGROW>
    end
end

if isempty(parts)
    T = table();
else
    T = vertcat(parts{:});
end
T = normalize_loaded_types(T);
end

function T = normalize_loaded_types(T)
if isempty(T), return; end
names = string(T.Properties.VariableNames);
string_vars = intersect(names, ["condition_key","geometry","realism_level","field_regime", ...
    "purity_bin","distance_bin","distance_over_window_bin","roi_region","model_name"], 'stable');
for v = string_vars
    T.(v) = string(T.(v));
end
numeric_vars = intersect(names, ["map_iz","map_ix","cx","cz","x_center_m","z_center_m", ...
    "source_seed","noise_seed","f0","M","REQ_M","true_SWS","k_true","patch_purity", ...
    "q_oracle","distance_to_interface_m","distance_to_interface_mm", ...
    "distance_to_interface_over_window_radius","q_theory_prior","sws_theory","q_pred", ...
    "k_pred","SWS_pred","sws_signed_error_pct","sws_abs_error_pct"], 'stable');
for v = numeric_vars
    if ~isnumeric(T.(v)), T.(v) = str2double(string(T.(v))); end
end
end

function T = attach_oracle_sws_from_req_cache(T, CFG)
T.SWS_oracle = nan(height(T),1);
T.k_oracle = nan(height(T),1);
T.q_oracle_from_cache = nan(height(T),1);
T.req_cache_found = false(height(T),1);

conds = unique(T.condition_key, 'stable');
fprintf('Reconstructing oracle SWS from %d Test66 REQ caches...\n', numel(conds));
for ci = 1:numel(conds)
    cond = conds(ci);
    cache_file = fullfile(CFG.ReqCacheDir, "req__" + sanitize(cond) + ".mat");
    idx = find(T.condition_key == cond);
    if exist(cache_file, 'file') ~= 2
        warning('Missing REQ cache for %s. Falling back to true SWS as oracle.', cond);
        T.SWS_oracle(idx) = T.true_SWS(idx);
        T.k_oracle(idx) = T.k_true(idx);
        T.q_oracle_from_cache(idx) = T.q_oracle(idx);
        continue;
    end

    S = load(cache_file, 'F');
    F = S.F;
    sws_oracle = q_to_sws(F.req_mapping, F.q_oracle, F.f0);
    k_oracle = 2*pi*F.f0 ./ sws_oracle;

    fkey = double(F.map_iz) * 1e6 + double(F.map_ix);
    rkey = double(T.map_iz(idx)) * 1e6 + double(T.map_ix(idx));
    [tf, loc] = ismember(rkey, fkey);
    if any(tf)
        ii = idx(tf);
        T.SWS_oracle(ii) = sws_oracle(loc(tf));
        T.k_oracle(ii) = k_oracle(loc(tf));
        T.q_oracle_from_cache(ii) = F.q_oracle(loc(tf));
        T.req_cache_found(ii) = true;
    end
end

missing = mean(~T.req_cache_found);
if missing > 0
    warning('%.1f%% of rows did not find a matching cache row.', 100*missing);
end
end

function T = add_diagnostic_variables(T)
T.q_error_model_oracle = T.q_pred - T.q_oracle_from_cache;
T.signed_error_pred = 100*(T.SWS_pred - T.true_SWS) ./ T.true_SWS;
T.abs_error_pred = abs(T.signed_error_pred);
T.signed_error_oracle = 100*(T.SWS_oracle - T.true_SWS) ./ T.true_SWS;
T.abs_error_oracle = abs(T.signed_error_oracle);
T.high20_pred = T.abs_error_pred > 20;
T.high20_oracle = T.abs_error_oracle > 20;
T.purity_bin_test67 = discretize_label(T.patch_purity, ...
    [0 0.50 0.75 0.90 0.95 0.99 1.000001], ...
    ["<0.50","0.50-0.75","0.75-0.90","0.90-0.95","0.95-0.99","0.99-1.00"]);
T.distance_over_window_bin_test67 = discretize_label(T.distance_to_interface_over_window_radius, ...
    [0 0.25 0.5 1 2 inf], ["0-0.25","0.25-0.5","0.5-1","1-2",">2"]);
T.distance_bin_test67 = discretize_label(1e3*T.distance_to_interface_m, ...
    [0 1 2 4 8 inf], ["0-1 mm","1-2 mm","2-4 mm","4-8 mm",">8 mm"]);
T.roi_region(T.roi_region == "" | ismissing(T.roi_region)) = "all_valid";
end

function lab = discretize_label(x, edges, names)
idx = discretize(x, edges);
lab = strings(numel(x),1);
for i = 1:numel(names)
    lab(idx == i) = names(i);
end
lab(lab == "") = "unbinned";
end

function SUM = summarize_metrics(T, group_vars)
if isempty(T)
    SUM = table();
    return;
end
Tg = T(T.model_name == "q_spectrum_plus_composition" | T.model_name == "q_spectrum_only", :);
group_vars = ["model_name", group_vars(:)'];
G = findgroups(Tg(:, cellstr(group_vars)));
[keys{1:numel(group_vars)}] = splitapply_multi(@first_value, Tg, group_vars, G); %#ok<CCAT>
SUM = table();
for i = 1:numel(group_vars)
    SUM.(group_vars(i)) = keys{i};
end
SUM.N = splitapply(@numel, Tg.abs_error_pred, G);
SUM.MAPE_pred = splitapply(@nanmean_local, Tg.abs_error_pred, G);
SUM.MAPE_oracle = splitapply(@nanmean_local, Tg.abs_error_oracle, G);
SUM.median_APE_pred = splitapply(@nanmedian_local, Tg.abs_error_pred, G);
SUM.median_APE_oracle = splitapply(@nanmedian_local, Tg.abs_error_oracle, G);
SUM.bias_pred = splitapply(@nanmean_local, Tg.signed_error_pred, G);
SUM.bias_oracle = splitapply(@nanmean_local, Tg.signed_error_oracle, G);
SUM.high20_pred = 100 * splitapply(@nanmean_local, double(Tg.high20_pred), G);
SUM.high20_oracle = 100 * splitapply(@nanmean_local, double(Tg.high20_oracle), G);
SUM.mean_q_pred = splitapply(@nanmean_local, Tg.q_pred, G);
SUM.mean_q_oracle = splitapply(@nanmean_local, Tg.q_oracle_from_cache, G);
SUM.mean_q_error = splitapply(@nanmean_local, Tg.q_error_model_oracle, G);
SUM.mean_distance_over_window_radius = splitapply(@nanmean_local, Tg.distance_to_interface_over_window_radius, G);
SUM.mean_patch_purity = splitapply(@nanmean_local, Tg.patch_purity, G);
end

function varargout = splitapply_multi(fun, T, vars, G)
varargout = cell(1,numel(vars));
for i = 1:numel(vars)
    varargout{i} = splitapply(fun, T.(vars(i)), G);
end
end

function y = first_value(x)
if iscell(x), x = string(x); end
y = x(1);
end

function y = nanmean_local(x)
y = mean(x, 'omitnan');
end

function y = nanmedian_local(x)
y = median(x, 'omitnan');
end

%% Figures

function make_representative_maps(T, OUT, CFG)
Tm = T(T.model_name == CFG.MainModel & T.realism_level == "clean", :);
conds = select_representative_conditions(Tm, CFG.MaxMapConditions);
for i = 1:numel(conds)
    C = Tm(Tm.condition_key == conds(i), :);
    if isempty(C), continue; end
    fig = figure('Visible','off','Color','w','Position',[80 80 1800 1100]);
    tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
    panels = {
        'true_SWS', 'True SWS', 'm/s'
        'SWS_pred', 'Predicted SWS', 'm/s'
        'SWS_oracle', 'Oracle SWS from q_{oracle}', 'm/s'
        'signed_error_pred', 'Signed error, q model', '%'
        'signed_error_oracle', 'Signed error, q oracle', '%'
        'abs_error_pred', 'Absolute error, q model', '%'
        'abs_error_oracle', 'Absolute error, q oracle', '%'
        'q_pred', 'q predicted', 'REQ quantile q'
        'q_oracle_from_cache', 'q oracle', 'REQ quantile q'
        'q_error_model_oracle', 'q predicted - q oracle', 'Delta q'
        'patch_purity', 'Patch purity', 'fraction'
        'distance_to_interface_over_window_radius', 'Distance / window radius', 'ratio'
        };
    for p = 1:size(panels,1)
        nexttile;
        plot_map_panel(C, panels{p,1}, panels{p,2}, panels{p,3});
    end
    sgtitle(strrep(C.condition_key(1), '_', '\_'), 'Interpreter','tex');
    out = fullfile(OUT.map_dir, "test67_map_q_oracle__" + sanitize(C.condition_key(1)) + ".png");
    exportgraphics(fig, out, 'Resolution', 180);
    close(fig);
end
end

function conds = select_representative_conditions(T, max_n)
keys = unique(T.condition_key, 'stable');
want = strings(0,1);
for geo = unique(T.geometry, 'stable')'
    for reg = ["single_source_lateral","diffuse_like_8src"]
        for f = [200 500]
            m = T.geometry == geo & T.field_regime == reg & T.f0 == f & T.source_seed == min(T.source_seed);
            k = unique(T.condition_key(m), 'stable');
            if ~isempty(k), want(end+1,1) = k(1); end %#ok<AGROW>
        end
    end
end
conds = unique(want, 'stable');
if numel(conds) < min(max_n, numel(keys))
    extra = setdiff(keys, conds, 'stable');
    conds = [conds; extra(1:min(numel(extra), max_n-numel(conds)))];
end
conds = conds(1:min(numel(conds), max_n));
end

function make_frequency_grids(T, OUT, CFG)
for geo = ["inclusion_2_3","inclusion_2_4"]
    for reg = CFG.FieldRegimes
        Tm = T(T.model_name == CFG.MainModel & T.geometry == geo & ...
            T.field_regime == reg & T.realism_level == "clean" & T.source_seed == 1, :);
        if isempty(Tm), continue; end
        fig = figure('Visible','off','Color','w','Position',[80 80 2200 1300]);
        tiledlayout(numel(CFG.Frequencies), 8, 'TileSpacing','compact','Padding','compact');
        cols = {
            'true_SWS','True SWS','m/s'
            'SWS_pred','Pred SWS','m/s'
            'SWS_oracle','Oracle SWS','m/s'
            'signed_error_pred','Signed err pred','%'
            'signed_error_oracle','Signed err oracle','%'
            'q_pred','q pred','q'
            'q_oracle_from_cache','q oracle','q'
            'q_error_model_oracle','q pred - q oracle','Delta q'
            };
        for f = CFG.Frequencies
            C = Tm(Tm.f0 == f, :);
            if isempty(C)
                for c = 1:size(cols,1), nexttile; axis off; end
                continue;
            end
            for c = 1:size(cols,1)
                nexttile;
                plot_map_panel(C, cols{c,1}, sprintf('%g Hz | %s', f, cols{c,2}), cols{c,3});
            end
        end
        sgtitle(sprintf('%s, %s, clean, source 1, M=2', pretty_name(geo), pretty_name(reg)));
        out = fullfile(OUT.freq_grid_dir, sprintf('test67_frequency_grid__%s__%s.png', geo, reg));
        exportgraphics(fig, out, 'Resolution', 170);
        close(fig);
    end
end
end

function make_directional_vs_diffuse_figures(T, OUT, CFG)
for geo = CFG.Geometries
    for f = [200 500]
        Tm = T(T.model_name == CFG.MainModel & T.geometry == geo & T.f0 == f & ...
            T.realism_level == "clean" & T.source_seed == 1, :);
        if isempty(Tm), continue; end
        fig = figure('Visible','off','Color','w','Position',[100 100 1550 1450]);
        tiledlayout(5,2,'TileSpacing','compact','Padding','compact');
        rows = {
            'SWS_pred','Predicted SWS','m/s'
            'SWS_oracle','Oracle SWS','m/s'
            'signed_error_pred','Signed error pred','%'
            'q_pred','q predicted','q'
            'q_error_model_oracle','q predicted - q oracle','Delta q'
            };
        for r = 1:size(rows,1)
            for reg = ["single_source_lateral","diffuse_like_8src"]
                nexttile;
                C = Tm(Tm.field_regime == reg, :);
                if isempty(C), axis off; title(sprintf('%s missing', reg)); continue; end
                plot_map_panel(C, rows{r,1}, sprintf('%s | %s', pretty_name(reg), rows{r,2}), rows{r,3});
            end
        end
        sgtitle(sprintf('%s, %g Hz, clean, source 1, M=2', pretty_name(geo), f));
        out = fullfile(OUT.field_compare_dir, sprintf('test67_single_vs_diffuse__%s__f%g.png', geo, f));
        exportgraphics(fig, out, 'Resolution', 180);
        close(fig);
    end
end
end

function make_error_vs_distance_figures(T, OUT, CFG)
Tm = T(T.model_name == CFG.MainModel, :);
Tm = sample_rows(Tm, CFG.MaxScatterPoints);
fig = figure('Visible','off','Color','w','Position',[100 100 1300 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; scatter_metric(Tm.distance_to_interface_mm, Tm.abs_error_pred, Tm.field_regime);
xlabel('Distance to interface (mm)'); ylabel('Absolute SWS error, q model (%)'); yline(20,'k--');
nexttile; scatter_metric(Tm.distance_to_interface_mm, Tm.abs_error_oracle, Tm.field_regime);
xlabel('Distance to interface (mm)'); ylabel('Absolute SWS error, q oracle (%)'); yline(20,'k--');
nexttile; scatter_metric(Tm.distance_to_interface_over_window_radius, Tm.abs_error_pred, Tm.field_regime);
xlabel('Distance / window radius'); ylabel('Absolute SWS error, q model (%)'); yline(20,'k--');
nexttile; scatter_metric(Tm.distance_to_interface_over_window_radius, Tm.abs_error_oracle, Tm.field_regime);
xlabel('Distance / window radius'); ylabel('Absolute SWS error, q oracle (%)'); yline(20,'k--');
sgtitle('Error versus interface distance, q model versus q oracle');
exportgraphics(fig, fullfile(OUT.dist_dir, 'test67_error_vs_interface_distance.png'), 'Resolution', 180);
close(fig);

Sd = summarize_metrics(T, ["field_regime","distance_over_window_bin_test67"]);
Sd = Sd(Sd.model_name == CFG.MainModel, :);
fig = figure('Visible','off','Color','w','Position',[100 100 1050 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; grouped_line(Sd, "distance_over_window_bin_test67", "MAPE_pred", "field_regime");
ylabel('MAPE q model (%)'); xlabel('Distance / window radius bin'); title('q model');
nexttile; grouped_line(Sd, "distance_over_window_bin_test67", "MAPE_oracle", "field_regime");
ylabel('MAPE q oracle (%)'); xlabel('Distance / window radius bin'); title('q oracle');
sgtitle('Binned error versus interface distance');
exportgraphics(fig, fullfile(OUT.dist_dir, 'test67_binned_error_vs_distance_over_window.png'), 'Resolution', 180);
close(fig);
end

function make_error_vs_purity_figures(T, OUT, CFG)
Sp = summarize_metrics(T, ["field_regime","purity_bin_test67"]);
Sp = Sp(Sp.model_name == CFG.MainModel, :);
fig = figure('Visible','off','Color','w','Position',[100 100 1200 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; grouped_line(Sp, "purity_bin_test67", "MAPE_pred", "field_regime");
ylabel('MAPE q model (%)'); xlabel('Patch purity bin'); title('MAPE, q model');
nexttile; grouped_line(Sp, "purity_bin_test67", "MAPE_oracle", "field_regime");
ylabel('MAPE q oracle (%)'); xlabel('Patch purity bin'); title('MAPE, q oracle');
nexttile; grouped_line(Sp, "purity_bin_test67", "high20_pred", "field_regime");
ylabel('High-error >20%, q model (%)'); xlabel('Patch purity bin'); title('High-error, q model');
nexttile; grouped_line(Sp, "purity_bin_test67", "high20_oracle", "field_regime");
ylabel('High-error >20%, q oracle (%)'); xlabel('Patch purity bin'); title('High-error, q oracle');
sgtitle('Error versus patch purity');
exportgraphics(fig, fullfile(OUT.purity_dir, 'test67_error_vs_patch_purity.png'), 'Resolution', 180);
close(fig);
end

function make_error_vs_q_gradient_figures(T, OUT, CFG)
[Gp, Go] = collect_q_gradient_rows(T, CFG);
if isempty(Gp), return; end
Gp = sample_rows(Gp, CFG.MaxScatterPoints);
fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; scatter_metric(Gp.grad_q_pred, Gp.abs_error_pred, Gp.field_regime);
xlabel('|grad q predicted|'); ylabel('Absolute SWS error, q model (%)'); yline(20,'k--');
nexttile; scatter_metric(Gp.grad_q_oracle, Gp.abs_error_oracle, Gp.field_regime);
xlabel('|grad q oracle|'); ylabel('Absolute SWS error, q oracle (%)'); yline(20,'k--');
sgtitle('Error versus q-map gradient');
exportgraphics(fig, fullfile(OUT.qgrad_dir, 'test67_error_vs_q_gradient.png'), 'Resolution', 180);
close(fig);

fig = figure('Visible','off','Color','w','Position',[100 100 900 450]);
histogram(Go.grad_q_pred, 80, 'DisplayStyle','stairs','LineWidth',1.5); hold on;
histogram(Go.grad_q_oracle, 80, 'DisplayStyle','stairs','LineWidth',1.5);
grid on; xlabel('|grad q|'); ylabel('Count'); legend({'q predicted','q oracle'},'Location','best');
title('q-gradient distribution across selected clean interface cases');
exportgraphics(fig, fullfile(OUT.qgrad_dir, 'test67_q_gradient_histogram.png'), 'Resolution', 180);
close(fig);
end

function make_q_maps_and_histograms(T, OUT, CFG)
Tm = T(T.model_name == CFG.MainModel & T.realism_level == "clean" & ...
    ismember(T.geometry, ["inclusion_2_3","inclusion_2_4"]) & ...
    T.field_regime == "diffuse_like_8src" & T.source_seed == 1, :);
for geo = ["inclusion_2_3","inclusion_2_4"]
    for f = [200 500]
        C = Tm(Tm.geometry == geo & Tm.f0 == f, :);
        if isempty(C), continue; end
        fig = figure('Visible','off','Color','w','Position',[100 100 1300 420]);
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        nexttile; plot_map_panel(C, 'q_pred', 'q predicted', 'REQ quantile q');
        nexttile; plot_map_panel(C, 'q_oracle_from_cache', 'q oracle', 'REQ quantile q');
        nexttile; plot_map_panel(C, 'q_error_model_oracle', 'q predicted - q oracle', 'Delta q');
        sgtitle(sprintf('%s diffuse-like, %g Hz', pretty_name(geo), f));
        exportgraphics(fig, fullfile(OUT.q_dir, sprintf('test67_q_maps__%s__f%g.png', geo, f)), 'Resolution', 180);
        close(fig);
    end
end

fig = figure('Visible','off','Color','w','Position',[100 100 1150 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
for reg = CFG.FieldRegimes
    histogram(T.q_error_model_oracle(T.model_name == CFG.MainModel & T.field_regime == reg), ...
        80, 'DisplayStyle','stairs','LineWidth',1.4); hold on;
end
grid on; xlabel('q predicted - q oracle'); ylabel('Count'); title('q error histogram');
legend(cellstr(pretty_name(CFG.FieldRegimes)), 'Location','best');
nexttile;
S = sample_rows(T(T.model_name == CFG.MainModel,:), CFG.MaxScatterPoints);
scatter(S.q_oracle_from_cache, S.q_pred, 5, S.abs_error_pred, 'filled', 'MarkerFaceAlpha',0.25);
grid on; axis square; xlabel('q oracle'); ylabel('q predicted'); cb = colorbar; ylabel(cb,'Absolute SWS error (%)');
title('q predicted versus q oracle');
exportgraphics(fig, fullfile(OUT.q_dir, 'test67_q_histogram_and_scatter.png'), 'Resolution', 180);
close(fig);
end

function make_local_spectra(T, OUT, CFG)
Tm = T(T.model_name == CFG.MainModel & T.realism_level == "clean" & T.source_seed == 1 & ...
    ismember(T.geometry, ["inclusion_2_4","bilayer_2_3"]) & T.f0 == 500, :);
conds = unique(Tm.condition_key, 'stable');
for ci = 1:min(numel(conds), 4)
    C = Tm(Tm.condition_key == conds(ci), :);
    cache_file = fullfile(CFG.ReqCacheDir, "req__" + sanitize(conds(ci)) + ".mat");
    if exist(cache_file, 'file') ~= 2, continue; end
    S = load(cache_file, 'F');
    F = S.F;
    points = select_spectral_points(C);
    if isempty(points), continue; end
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 700]);
    tiledlayout(2, numel(points), 'TileSpacing','compact','Padding','compact');
    for pi = 1:numel(points)
        row = points(pi);
        fidx = find(F.map_iz == C.map_iz(row) & F.map_ix == C.map_ix(row), 1);
        if isempty(fidx), continue; end
        mapping = F.req_mapping{fidx};
        k = double(mapping.k_cent(:));
        E = double(mapping.Ecum(:));
        Srad = [0; max(diff(E), 0)];
        kvals = [C.k_true(row), C.k_pred(row), C.k_oracle(row)];
        labs = {'k true','k predicted','k oracle'};
        nexttile(pi);
        plot(k, Srad, 'LineWidth',1.2); grid on; hold on;
        add_k_lines(kvals, labs);
        xlabel('Wavenumber k (rad/m)'); ylabel('dE/dk proxy');
        title(sprintf('%s | %s', pretty_name(C.roi_region(row)), pretty_name(C.field_regime(row))));
        nexttile(pi + numel(points));
        plot(k, E, 'LineWidth',1.2); grid on; hold on;
        add_k_lines(kvals, labs);
        yline(C.q_pred(row), 'm--', 'q pred');
        yline(C.q_oracle_from_cache(row), 'g--', 'q oracle');
        xlabel('Wavenumber k (rad/m)'); ylabel('Cumulative energy q');
        title(sprintf('purity %.2f, dist/window %.2f', C.patch_purity(row), C.distance_to_interface_over_window_radius(row)));
    end
    sgtitle(strrep(conds(ci), '_', '\_'), 'Interpreter','tex');
    out = fullfile(OUT.spectra_dir, "test67_local_spectra__" + sanitize(conds(ci)) + ".png");
    exportgraphics(fig, out, 'Resolution', 180);
    close(fig);
end
end

function rows = select_spectral_points(C)
rows = [];
labels = ["hard_core","hard_inclusion_core","interface_0_1mm","background_far","soft_core"];
for lab = labels
    m = C.roi_region == lab;
    if ~any(m)
        if contains(lab, "interface")
            m = C.distance_to_interface_over_window_radius <= 0.5;
        elseif contains(lab, "background") || contains(lab, "soft")
            m = C.patch_purity >= 0.99 & C.true_SWS <= median(C.true_SWS,'omitnan');
        elseif contains(lab, "hard")
            m = C.patch_purity >= 0.99 & C.true_SWS >= median(C.true_SWS,'omitnan');
        end
    end
    idx = find(m & isfinite(C.q_pred) & isfinite(C.q_oracle_from_cache));
    if isempty(idx), continue; end
    [~, j] = min(abs(C.distance_to_interface_over_window_radius(idx) - median(C.distance_to_interface_over_window_radius(idx),'omitnan')));
    rows(end+1) = idx(j); %#ok<AGROW>
end
rows = unique(rows, 'stable');
rows = rows(1:min(numel(rows),3));
end

%% README

function write_readme(T, SUM, OUT, CFG)
main = T(T.model_name == CFG.MainModel, :);
field = SUM.field(SUM.field.model_name == CFG.MainModel, :);
dist = SUM.distance(SUM.distance.model_name == CFG.MainModel, :);
purity = SUM.purity(SUM.purity.model_name == CFG.MainModel, :);

near = main.distance_to_interface_over_window_radius <= 0.5;
far = main.distance_to_interface_over_window_radius > 1;
lowf = ismember(main.f0, [200 300]);
highf = ismember(main.f0, [500 600]);
single = main.field_regime == "single_source_lateral";
diffuse = main.field_regime == "diffuse_like_8src";

oracle_near = mean(main.abs_error_oracle(near), 'omitnan');
pred_near = mean(main.abs_error_pred(near), 'omitnan');
oracle_far = mean(main.abs_error_oracle(far), 'omitnan');
pred_far = mean(main.abs_error_pred(far), 'omitnan');
qerr_near = mean(main.q_error_model_oracle(near), 'omitnan');
qerr_far = mean(main.q_error_model_oracle(far), 'omitnan');
lowf_mape = mean(main.abs_error_pred(lowf), 'omitnan');
highf_mape = mean(main.abs_error_pred(highf), 'omitnan');
single_mape = mean(main.abs_error_pred(single), 'omitnan');
diffuse_mape = mean(main.abs_error_pred(diffuse), 'omitnan');
single_oracle = mean(main.abs_error_oracle(single), 'omitnan');
diffuse_oracle = mean(main.abs_error_oracle(diffuse), 'omitnan');

if oracle_near > 10
    cause = "Both the q-model and the q-oracle show substantial interface error, so the dominant component is likely physical/interface mixing or an operational REQ limit.";
elseif pred_near > 10 && oracle_near <= 5
    cause = "The q-model shows the interface band but the q-oracle largely removes it, so the dominant component is q-model domain shift around Eikonal/diffuse interface patches.";
else
    cause = "The interface band is moderate; the evidence points to a mixed contribution from q-model error and interface/field complexity.";
end
if single_mape < 0.7*diffuse_mape && single_oracle <= 0.7*diffuse_oracle
    field_interp = "The diffuse-like field is worse even for the oracle, consistent with multi-source/interference making local REQ less clean.";
elseif single_mape < 0.7*diffuse_mape
    field_interp = "The diffuse-like field is worse mostly for the q-model, consistent with field-regime domain shift.";
else
    field_interp = "Single-source and diffuse-like fields are comparable at the aggregate level.";
end

lines = [
    "# Test 67: q-oracle interface diagnostics"
    ""
    "This analysis reads existing Test 66 outputs. It does not retrain models, rerun simulations, or modify Test 66."
    ""
    "## Main quantitative checks"
    sprintf("- Near-interface MAPE with q model: %.2f%%", pred_near)
    sprintf("- Near-interface MAPE with q oracle: %.2f%%", oracle_near)
    sprintf("- Far-from-interface MAPE with q model: %.2f%%", pred_far)
    sprintf("- Far-from-interface MAPE with q oracle: %.2f%%", oracle_far)
    sprintf("- Mean q error near interface: %.4f", qerr_near)
    sprintf("- Mean q error far from interface: %.4f", qerr_far)
    sprintf("- MAPE at 200/300 Hz: %.2f%%", lowf_mape)
    sprintf("- MAPE at 500/600 Hz: %.2f%%", highf_mape)
    sprintf("- MAPE single-source lateral: %.2f%%", single_mape)
    sprintf("- MAPE diffuse-like 8 source: %.2f%%", diffuse_mape)
    ""
    "## Explicit answers"
    sprintf("1. Does the band appear in oracle SWS? Near-interface oracle MAPE is %.2f%%. %s", oracle_near, yes_no(oracle_near > 10))
    sprintf("2. Does q_pred deviate from q_oracle at the interface? Mean near-interface q error is %.4f versus %.4f far from interface.", qerr_near, qerr_far)
    sprintf("3. Does error drop beyond one window radius? q-model MAPE changes from %.2f%% near to %.2f%% far.", pred_near, pred_far)
    sprintf("4. Is the band worse at 200/300 Hz? Low-frequency MAPE is %.2f%% versus %.2f%% at 500/600 Hz.", lowf_mape, highf_mape)
    sprintf("5. Is the band worse in diffuse-like fields? Diffuse-like MAPE is %.2f%% versus %.2f%% single-source.", diffuse_mape, single_mape)
    sprintf("6. Main interpretation: %s %s", cause, field_interp)
    "7. Non-estimable candidates: patches with distance/window < 0.5, patch purity < 0.95, or high predicted q-oracle disagreement should be flagged first."
    ""
    "## Generated tables"
    "- `tables/test67_selected_patch_level_results.csv`"
    "- `tables/summary_q_pred_vs_q_oracle.csv`"
    "- `tables/summary_by_field_regime.csv`"
    "- `tables/summary_by_patch_purity.csv`"
    "- `tables/summary_by_distance_over_window.csv`"
    "- `tables/summary_by_roi.csv`"
    ""
    "## Generated figures"
    "- `figures/maps_q_oracle_diagnostics/`"
    "- `figures/frequency_grids/`"
    "- `figures/directional_vs_diffuse_clean/`"
    "- `figures/error_vs_interface_distance/`"
    "- `figures/error_vs_patch_purity/`"
    "- `figures/error_vs_q_gradient/`"
    "- `figures/q_maps_and_histograms/`"
    "- `figures/local_spectra/`"
    ""
    "## Notes"
    "- `q_oracle` is diagnostic. It uses true SWS to identify the q value that would select the true k from the local REQ mapping."
    "- A low oracle error means the local REQ mapping can represent the truth and the learned q model is the limiting factor."
    "- A high oracle error means the local REQ mapping itself, the window, or the field complexity is limiting the estimate."
    ];

fid = fopen(fullfile(OUT.root_dir, 'README_results.md'), 'w');
fprintf(fid, '%s\n', lines);
fclose(fid);

% Compact console preview.
disp('Top field-regime summary for main model:');
disp(field(:, intersect(["geometry","f0","field_regime","MAPE_pred","MAPE_oracle","bias_pred","bias_oracle","mean_q_error"], string(field.Properties.VariableNames), 'stable')));
disp('Purity summary for main model:');
disp(purity(:, intersect(["field_regime","purity_bin_test67","MAPE_pred","MAPE_oracle","high20_pred","high20_oracle"], string(purity.Properties.VariableNames), 'stable')));
disp('Distance/window summary for main model:');
disp(dist(:, intersect(["field_regime","distance_over_window_bin_test67","MAPE_pred","MAPE_oracle","high20_pred","high20_oracle"], string(dist.Properties.VariableNames), 'stable')));
end

function s = yes_no(tf)
if tf, s = "Yes, this suggests an oracle/REQ-interface component."; else, s = "No or only weakly, which points more toward q-model error."; end
end

%% Plot helpers

function plot_map_panel(T, value_col, title_text, cbar_label)
[A, x_mm, z_mm] = table_to_map(T, value_col);
imagesc(x_mm, z_mm, A);
axis image; set(gca,'YDir','normal'); grid off;
title(title_text, 'Interpreter','none');
xlabel('x (mm)'); ylabel('z (mm)');
cb = colorbar; ylabel(cb, cbar_label, 'Interpreter','none');
end

function [A, x_mm, z_mm] = table_to_map(T, value_col)
nr = max(T.map_iz); nc = max(T.map_ix);
A = nan(nr, nc);
ind = sub2ind([nr nc], T.map_iz, T.map_ix);
A(ind) = T.(value_col);
x_mm = nan(1,nc); z_mm = nan(nr,1);
[~, ixu] = unique(T.map_ix, 'stable'); x_mm(T.map_ix(ixu)) = 1e3*T.x_center_m(ixu);
[~, izu] = unique(T.map_iz, 'stable'); z_mm(T.map_iz(izu)) = 1e3*T.z_center_m(izu);
if any(isnan(x_mm)), x_mm = 1:nc; end
if any(isnan(z_mm)), z_mm = (1:nr)'; end
end

function scatter_metric(x, y, g)
cats = unique(g, 'stable');
hold on;
for c = cats'
    m = g == c & isfinite(x) & isfinite(y);
    scatter(x(m), y(m), 8, 'filled', 'MarkerFaceAlpha',0.22, 'DisplayName', pretty_name(c));
end
grid on; legend('Location','best');
end

function grouped_line(T, xvar, yvar, gvar)
cats = unique(T.(gvar), 'stable');
hold on;
for c = cats'
    S = T(T.(gvar) == c, :);
    [x, ord] = sort_categorical_labels(S.(xvar));
    plot(categorical(x), S.(yvar)(ord), '-o', 'LineWidth',1.3, 'DisplayName', pretty_name(c));
end
grid on; legend('Location','best'); xtickangle(30);
end

function [x, ord] = sort_categorical_labels(x)
x = string(x);
order = ["<0.50","0.50-0.75","0.75-0.90","0.90-0.95","0.95-0.99","0.99-1.00", ...
    "0-0.25","0.25-0.5","0.5-1","1-2",">2", ...
    "single_source_lateral","diffuse_like_8src"];
score = nan(numel(x),1);
for i = 1:numel(x)
    j = find(order == x(i), 1);
    if isempty(j), score(i) = i; else, score(i) = j; end
end
[~, ord] = sort(score);
x = x(ord);
end

function [G, Tall] = collect_q_gradient_rows(T, CFG)
Tm = T(T.model_name == CFG.MainModel & T.realism_level == "clean", :);
conds = unique(Tm.condition_key, 'stable');
parts = cell(numel(conds),1);
for ci = 1:numel(conds)
    C = Tm(Tm.condition_key == conds(ci), :);
    if isempty(C), continue; end
    [Qp,~,~] = table_to_map(C, 'q_pred');
    [Qo,~,~] = table_to_map(C, 'q_oracle_from_cache');
    Gp = gradient_magnitude(Qp);
    Go = gradient_magnitude(Qo);
    kp = double(C.map_iz) * 1e6 + double(C.map_ix);
    [rr,cc] = ndgrid(1:size(Qp,1),1:size(Qp,2));
    km = double(rr(:))*1e6 + double(cc(:));
    [tf, loc] = ismember(kp, km);
    C.grad_q_pred = nan(height(C),1);
    C.grad_q_oracle = nan(height(C),1);
    C.grad_q_pred(tf) = Gp(loc(tf));
    C.grad_q_oracle(tf) = Go(loc(tf));
    parts{ci} = C;
end
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    G = table(); Tall = table();
else
    G = vertcat(parts{:});
    Tall = G;
end
end

function G = gradient_magnitude(A)
[gx, gz] = gradient(A);
G = sqrt(gx.^2 + gz.^2);
end

function T = sample_rows(T, max_n)
if height(T) <= max_n, return; end
rng(6701);
idx = randperm(height(T), max_n);
T = T(idx,:);
end

function add_k_lines(kvals, labels)
cols = lines(numel(kvals));
for i = 1:numel(kvals)
    if isfinite(kvals(i))
        xline(kvals(i), '--', labels{i}, 'Color', cols(i,:), 'LineWidth',1.1);
    end
end
end

%% Utility

function y = q_to_sws(mappings, q, f0)
y = nan(numel(q),1);
for i = 1:numel(q)
    if isscalar(f0), fi = f0; else, fi = f0(i); end
    y(i) = adaptive_req.quantile.quantile_to_cs(mappings{i}, q(i), fi);
end
end

function s = sanitize(x)
s = regexprep(string(x), '[^A-Za-z0-9_=-]+', '_');
s = char(s);
end

function s = pretty_name(x)
x = string(x);
s = strrep(x, "_", " ");
s = regexprep(s, '\bq\b', 'q');
end

function tf = env_true(name, default_value)
v = string(getenv(name));
if strlength(v) == 0
    tf = default_value;
else
    tf = any(strcmpi(v, ["1","true","yes","on"]));
end
end

function v = env_number(name, default_value)
s = string(getenv(name));
if strlength(s) == 0
    v = default_value;
else
    v = str2double(s);
    if ~isfinite(v), v = default_value; end
end
end

function s = env_string(name, default_value)
s = string(getenv(name));
if strlength(s) == 0, s = string(default_value); end
end
