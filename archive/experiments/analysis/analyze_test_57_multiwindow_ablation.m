%% analyze_test_57_multiwindow_ablation.m
% Test 57: multi-window ablation for frozen q-spectrum models.
%
% This script compares fixed-M maps and simple operational multi-window
% selectors/blends. It expects a Test55-style patch-level CSV generated with
% more than one M, for example M=[1.5 2 2.5 3]. If only one M exists, it writes
% fixed-M summaries and prints the exact rerun command needed.

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
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir,'test57_configuration.json'));

fprintf('\nTest 57: multi-window ablation\n');
fprintf('Source: %s\n', CFG.SourceCsv);
fprintf('Base model: %s | no model training.\n', CFG.BaseModel);

T0 = readtable(CFG.SourceCsv);
T0 = T0(string(T0.model_name) == CFG.BaseModel, :);
assert(~isempty(T0), 'No rows found for model %s.', CFG.BaseModel);
T0 = add_operational_features(T0, CFG);
if CFG.ValidateOnly
    T0 = sample_rows(T0, min(height(T0), CFG.ValidateRows), CFG.RandomSeed);
end

Mvals = unique(T0.M)';
fprintf('Available M values in source: %s\n', mat2str(Mvals));

T_fixed = fixed_m_long(T0);
T_multi = table();
if numel(Mvals) > 1
    T_multi = multiwindow_strategies(T0, CFG);
else
    warning('Test57:SingleM', ['Only one M value found. Multi-window strategies ', ...
        'cannot be evaluated until Test55 is rerun with ADAPTIVE_REQ_TEST55_M_LIST.']);
end
T_all = [T_fixed; T_multi];

T_overall = summarize_predictions(T_all, "strategy_name");
T_by_M = summarize_predictions(T_all, ["strategy_name","source_M"]);
T_by_case = summarize_predictions(T_all, ["strategy_name","case_id"]);
T_by_freq = summarize_predictions(T_all, ["strategy_name","f0"]);
T_by_regime = summarize_predictions(T_all, ["strategy_name","field_regime_ood"]);
T_by_roi = summarize_predictions(T_all(T_all.roi_region ~= "other",:), ...
    ["strategy_name","roi_region"]);
T_by_case_roi = summarize_predictions(T_all(T_all.roi_region ~= "other",:), ...
    ["strategy_name","case_id","roi_region"]);

writetable(remove_cell_columns(T_all), fullfile(OUT.table_dir,'test57_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test57_strategy_summary_overall.csv'));
writetable(T_by_M, fullfile(OUT.table_dir,'test57_strategy_summary_by_M.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test57_strategy_summary_by_case.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir,'test57_strategy_summary_by_frequency.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test57_strategy_summary_by_regime.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir,'test57_strategy_summary_by_roi.csv'));
writetable(T_by_case_roi, fullfile(OUT.table_dir,'test57_strategy_summary_by_case_roi.csv'));

safe_plot(@() plot_test57_summary(T_overall, T_by_M, T_by_case_roi, OUT), 'summary plots');
save(fullfile(OUT.data_dir,'test57_multiwindow_results.mat'), ...
    'CFG','T_overall','T_by_M','T_by_case','T_by_freq','T_by_regime', ...
    'T_by_roi','T_by_case_roi','-v7.3');

print_summary(T_overall, T_by_M, OUT, numel(Mvals));
fprintf('\nTables: %s\nFigures: %s\nTest 57 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST57_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST57_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST57_VALIDATE_ONLY', false);
CFG.RandomSeed = 57001;
CFG.SourceCsv = env_string('ADAPTIVE_REQ_TEST57_SOURCE_CSV', ...
    fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep', ...
    'tables','test55_patch_level_results.csv'));
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST57_BASE_MODEL', ...
    "q_spectrum_plus_composition");
CFG.CsGuess = env_number('ADAPTIVE_REQ_TEST57_CS_GUESS', 3.0);
CFG.TargetMeff = env_number('ADAPTIVE_REQ_TEST57_TARGET_MEFF', 2.0);
CFG.ValidateRows = env_number('ADAPTIVE_REQ_TEST57_VALIDATE_ROWS', 30000);
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_57_multiwindow_ablation','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_57_multiwindow_ablation');
end
if CFG.ValidateOnly
    OUT.root_dir = fullfile(OUT.root_dir,'validate');
end
OUT.table_dir = fullfile(OUT.root_dir,'tables');
OUT.figure_dir = fullfile(OUT.root_dir,'figures');
OUT.data_dir = fullfile(OUT.root_dir,'data');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

%% Strategies

function T = add_operational_features(T, CFG)
T.sws_pred = clamp(T.sws_pred, 0.2, 20);
T.k_pred_firstpass = 2*pi*T.f0 ./ T.sws_pred;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.lambda_pred_m = T.sws_pred ./ T.f0;
T.M_eff_pred = T.window_length_m ./ T.lambda_pred_m;
base_condition = regexprep(string(T.condition_key), "__M[0-9.]+$", "");
if all(ismember(["x_center_m","z_center_m"], string(T.Properties.VariableNames)))
    x_um = round(1e6*T.x_center_m);
    z_um = round(1e6*T.z_center_m);
    T.align_key = base_condition + "__xum" + string(x_um) + "__zum" + string(z_um);
else
    T.align_key = base_condition + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
end
end

function T = fixed_m_long(T0)
T = make_rows(T0, "fixed_M" + string(T0.M), T0.sws_pred, T0.M, false(height(T0),1));
end

function T = multiwindow_strategies(T0, CFG)
T0 = sortrows(T0, {'align_key','M'});
[G, ~] = findgroups(T0.align_key);
if max(G) < 1
    T = table();
    return;
end

counts = splitapply(@numel, T0.sws_pred, G);
[~, first_idx] = unique(G, 'first');
keep_group = counts >= 2;
base = T0(first_idx(keep_group),:);
if isempty(base)
    T = table();
    return;
end

mean_sws = splitapply(@(x) mean(x,'omitnan'), T0.sws_pred, G);
median_sws = splitapply(@(x) median(x,'omitnan'), T0.sws_pred, G);
[lowest_mix_sws, lowest_mix_M] = splitapply(@pick_lowest_score, ...
    T0.sws_pred, T0.M, T0.p_mixed, G);
[highest_purity_sws, highest_purity_M] = splitapply(@pick_highest_score, ...
    T0.sws_pred, T0.M, T0.predicted_patch_purity, G);
[meff_sws, meff_M] = splitapply(@(s,m,meff) pick_meff_target(s,m,meff,CFG.TargetMeff), ...
    T0.sws_pred, T0.M, T0.M_eff_pred, G);

mean_sws = mean_sws(keep_group);
median_sws = median_sws(keep_group);
lowest_mix_sws = lowest_mix_sws(keep_group);
lowest_mix_M = lowest_mix_M(keep_group);
highest_purity_sws = highest_purity_sws(keep_group);
highest_purity_M = highest_purity_M(keep_group);
meff_sws = meff_sws(keep_group);
meff_M = meff_M(keep_group);

T = [ ...
    make_rows(base, "multiwindow_mean_sws", mean_sws, NaN(height(base),1), true); ...
    make_rows(base, "multiwindow_median_sws", median_sws, NaN(height(base),1), true); ...
    make_rows(base, "multiwindow_lowest_predicted_mixedness", lowest_mix_sws, lowest_mix_M, true); ...
    make_rows(base, "multiwindow_highest_predicted_purity", highest_purity_sws, highest_purity_M, true); ...
    make_rows(base, "multiwindow_meff_near_target", meff_sws, meff_M, true) ...
    ];
end

function [sws, M] = pick_lowest_score(sws_vec, M_vec, score_vec)
[~, idx] = min(score_vec);
sws = sws_vec(idx);
M = M_vec(idx);
end

function [sws, M] = pick_highest_score(sws_vec, M_vec, score_vec)
[~, idx] = max(score_vec);
sws = sws_vec(idx);
M = M_vec(idx);
end

function [sws, M] = pick_meff_target(sws_vec, M_vec, meff_vec, target)
[~, idx] = min(abs(meff_vec - target));
sws = sws_vec(idx);
M = M_vec(idx);
end

function R = make_rows(T, strategy, sws, source_M, was_combined)
keep = ["dataset","condition_key","case_id","case_family","seen_status", ...
    "field_regime","field_regime_ood","f0","M","dx","dz", ...
    "true_SWS","k_true","patch_purity","purity_bin","distance_to_boundary_mm", ...
    "distance_bin","predicted_patch_purity","p_mixed","p_strong_mixed", ...
    "material_region","roi_region","q_pred","sws_pred","M_eff_pred"];
keep = intersect(keep, string(T.Properties.VariableNames), 'stable');
R = T(:, cellstr(keep));
strategy = string(strategy(:));
if isscalar(strategy), strategy = repmat(strategy, height(R), 1); end
R.strategy_name = strategy;
if isscalar(sws), sws = repmat(sws, height(R), 1); end
if isscalar(source_M), source_M = repmat(source_M, height(R), 1); end
if isscalar(was_combined), was_combined = repmat(was_combined, height(R), 1); end
R.sws_multiwindow = sws;
R.source_M = source_M;
R.was_combined = was_combined;
R.sws_signed_error_pct = 100*(R.sws_multiwindow - R.true_SWS)./R.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
end

%% Summaries

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T), S = table(); return; end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        100*mean(X.high_error10,'omitnan'), 100*mean(X.high_error20,'omitnan'), ...
        100*mean(X.sws_signed_error_pct < 0,'omitnan'), mean(X.sws_multiwindow,'omitnan'), ...
        mean(X.true_SWS,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','high_error10_pct','high_error20_pct', ...
        'underestimate_pct','mean_pred_sws','mean_true_sws'});
end
S = [groups vertcat(rows{:})];
end

function plot_test57_summary(T_overall, T_by_M, T_roi, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 16]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
X = sortrows(T_overall,'MAPE_pct','ascend');
barh_labels(ax, pretty_strategy(X.strategy_name), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Multi-window strategy ranking','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
F = T_by_M(startsWith(string(T_by_M.strategy_name),"fixed_M"),:);
F = sortrows(F,'source_M');
bar(ax, F.source_M, F.MAPE_pct); xlabel(ax,'REQ M'); ylabel(ax,'MAPE (%)');
title(ax,'Fixed-M ablation','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
R = T_roi(ismember(T_roi.roi_region, ["homogeneous_center_8mm","soft_core_gt8mm","hard_core_gt4mm"]),:);
R = sortrows(R,'MAPE_pct','descend');
barh_labels(ax, pretty_strategy(R.strategy_name)+" | "+string(R.case_id)+" | "+pretty_region(R.roi_region), R.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Core ROI errors','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test57_multiwindow_summary.png'));
end

function print_summary(T_overall, T_by_M, OUT, nM)
fprintf('\n================ Test 57 summary ================\n');
disp(T_overall(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct'}));
if nM <= 1
    fprintf(['Only one M was available. To run the actual multi-window ablation:\n', ...
        '  setenv(''ADAPTIVE_REQ_TEST55_M_LIST'',''1.5 2 2.5 3'') and rerun Test55, then rerun Test57.\n']);
else
    X = T_overall(~startsWith(string(T_overall.strategy_name),"fixed_M"),:);
    if ~isempty(X)
        [~,i] = min(X.MAPE_pct);
        fprintf('Best multi-window strategy: %s, MAPE %.2f%%.\n', X.strategy_name(i), X.MAPE_pct(i));
    end
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

%% Helpers

function barh_labels(ax, labels, values)
y = 1:numel(values);
barh(ax, y, values);
set(ax, 'YTick', y, 'YTickLabel', cellstr(string(labels)));
end

function T = sample_rows(T, n, seed)
rng(seed);
idx = randperm(height(T), n);
T = T(idx,:);
end

function y = clamp(x, lo, hi), y = min(max(x,lo),hi); end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
end

function s = pretty_strategy(x)
s = string(x);
s = strrep(s,"fixed_M","Fixed M=");
s(s=="multiwindow_mean_sws") = "Mean SWS across windows";
s(s=="multiwindow_median_sws") = "Median SWS across windows";
s(s=="multiwindow_lowest_predicted_mixedness") = "Lowest predicted mixedness";
s(s=="multiwindow_highest_predicted_purity") = "Highest predicted purity";
s(s=="multiwindow_meff_near_target") = "M-eff closest to target";
s = strrep(s,"_"," ");
end

function s = pretty_region(x)
s = string(x);
s = strrep(s,"_"," ");
s = strrep(s,"0p5","0.5");
s = strrep(s,"gt",">");
end

function export_fig(fig, path)
folder = fileparts(path);
if exist(folder,'dir') ~= 7, mkdir(folder); end
exportgraphics(fig, path, 'Resolution', 220);
close(fig);
end

function safe_plot(fn, label)
try
    fn();
catch ME
    warning('Test57:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w'); fprintf(fid,'%s\n',txt); fclose(fid);
end

function tf = env_true(name, default)
v = lower(strtrim(string(getenv(name))));
if v == "", tf = default; else, tf = ismember(v, ["1","true","yes","on"]); end
end

function x = env_number(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x), x = default; end
end

function s = env_string(name, default)
s = string(getenv(name));
if strlength(strtrim(s)) == 0, s = string(default); end
end
