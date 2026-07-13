%% analyze_test_41_mixedness_corrected_latest_kwave.m
% Test 41: mixedness-style correction on latest clean-q models for k-Wave.
%
% Goal:
%   Apply a Test33-like residual log-k corrector to the latest frozen q models
%   already evaluated on k-Wave in Test 40, and compare against the old
%   Test33/Test34 correction stack.
%
% Notes:
%   - No base q model is retrained.
%   - The residual corrector is a post-map model trained only from operational
%     variables: model prediction, frozen confidence, predicted mixedness /
%     predicted purity, M, field regime, and operational theory/local
%     references already available in Test 40.
%   - True SWS is used only as the residual correction label and for evaluation.
%   - True side, ROI, and distance-to-interface are evaluation variables only.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST41_MODE               = quick | full
%   ADAPTIVE_REQ_TEST41_VALIDATE_ONLY      = true | false
%   ADAPTIVE_REQ_TEST41_SAVE_ALL_MAPS      = true | false
%   ADAPTIVE_REQ_TEST41_MAX_MAP_CONDITIONS = integer or Inf

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.10);

%% Configuration

CFG = struct();
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST41_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST41_MODE must be quick or full.');
CFG.RunMode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST41_VALIDATE_ONLY', false);
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST41_SAVE_ALL_MAPS', true);
CFG.MaxMapConditions = env_number('ADAPTIVE_REQ_TEST41_MAX_MAP_CONDITIONS', ternary(CFG.QuickMode, 2, Inf));
CFG.RandomSeed = 41001;
CFG.NumLearningCycles = env_number('ADAPTIVE_REQ_TEST41_NUM_TREES', ternary(CFG.QuickMode, 80, 160));
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST41_MIN_LEAF_SIZE', 80);
CFG.KFoldByCondition = env_number('ADAPTIVE_REQ_TEST41_KFOLD', 4);
CFG.PhysicalRange = [0.5 10];
CFG.OutputRoot = fullfile(root_dir, 'outputs', 'test_41_mixedness_corrected_latest_kwave');
CFG.ReferenceStrategies = ["T34_local_baseline","T34_hybrid_baseline", ...
    "T34_theory_baseline","T34_mixedness_logk_corrected", ...
    "T34_mixedness_q_candidate_selector"];
if CFG.QuickMode
    CFG.BaseNewStrategies = ["T38_q_spectrum_plus_composition", ...
        "T38_q_spectrum_only"];
else
    CFG.BaseNewStrategies = ["T38_q_spectrum_plus_composition", ...
        "T38_q_spectrum_plus_theory_composition", ...
        "T38_q_spectrum_only", ...
        "T35_q_spectrum_plus_theory_composition"];
end
CFG.MainPlotStrategies = ["T34_local_baseline","T34_theory_baseline", ...
    "T34_mixedness_logk_corrected", ...
    "T38_q_spectrum_plus_composition", ...
    "T41_T38_q_spectrum_plus_composition_logk_all", ...
    "T41_T38_q_spectrum_plus_composition_logk_mixed_gated", ...
    "T38_q_spectrum_only", ...
    "T41_T38_q_spectrum_only_logk_mixed_gated"];

OUT = make_output_dirs(CFG.OutputRoot);
write_config_json(CFG, fullfile(OUT.root_dir, 'test41_configuration.json'));
SRC = locate_sources(root_dir);

fprintf('\nTest 41: mixedness-corrected latest q models on k-Wave\n');
fprintf('Mode: %s | output: %s\n', CFG.RunMode, OUT.root_dir);
fprintf('Base models: %s\n', strjoin(CFG.BaseNewStrategies, ', '));

assert(exist(SRC.test40_patch_csv,'file') == 2, ...
    'Missing Test 40 patch table. Run Test 40 first: %s', SRC.test40_patch_csv);

T40 = readtable(SRC.test40_patch_csv);
T40 = standardize_types(T40);
T40 = attach_operational_references(T40);
fprintf('Loaded Test40 k-Wave rows: %d.\n', height(T40));

if CFG.ValidateOnly
    validate_inputs(T40, CFG);
    fprintf('Test 41 validation-only checks passed.\n');
    return;
end

rng(CFG.RandomSeed);

%% Build reference and corrected tables

keep_ref = ismember(T40.strategy_name, [CFG.ReferenceStrategies CFG.BaseNewStrategies]);
T_ref = T40(keep_ref,:);
T_ref.correction_variant = repmat("none", height(T_ref), 1);
T_ref.base_strategy_name = T_ref.strategy_name;
T_ref.mixedness_gate = nan(height(T_ref),1);
T_ref.delta_logk_pred = nan(height(T_ref),1);
T_ref.fold_id = nan(height(T_ref),1);

parts = {T_ref};
model_summaries = cell(numel(CFG.BaseNewStrategies),1);
for i = 1:numel(CFG.BaseNewStrategies)
    base = CFG.BaseNewStrategies(i);
    B = T40(T40.strategy_name == base,:);
    assert(~isempty(B), 'Base strategy missing in Test40 table: %s', base);
    fprintf('[%d/%d] Training grouped-CV mixedness log-k corrector for %s (%d rows)...\n', ...
        i, numel(CFG.BaseNewStrategies), base, height(B));
    [C, S_model] = train_apply_logk_corrector(B, base, CFG);
    parts{end+1,1} = C; %#ok<AGROW>
    model_summaries{i} = S_model;
end

T_all = vertcat_compatible(parts{:});
T_model_cv = vertcat(model_summaries{:});

%% Summaries

T_overall = summarize_long(T_all, ["model_family","strategy_name"]);
T_by_side = summarize_long(T_all, ["strategy_name","material_side"]);
T_by_roi = summarize_long(T_all, ["strategy_name","roi_name"]);
T_by_M = summarize_long(T_all, ["strategy_name","M"]);
T_by_regime = summarize_long(T_all, ["strategy_name","field_regime"]);
T_by_region = summarize_long(T_all, ["strategy_name","analysis_region"]);
T_by_distance = summarize_long(T_all, ["strategy_name","distance_bin"]);
T_by_side_regime = summarize_long(T_all, ["strategy_name","material_side","field_regime"]);
T_by_M_regime = summarize_long(T_all, ["strategy_name","M","field_regime"]);

T33_ref = load_test33_reference(SRC);

writetable(T_all, fullfile(OUT.table_dir, 'test41_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir, 'test41_summary_overall.csv'));
writetable(T_by_side, fullfile(OUT.table_dir, 'test41_summary_by_soft_hard.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir, 'test41_summary_by_roi.csv'));
writetable(T_by_M, fullfile(OUT.table_dir, 'test41_summary_by_M.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir, 'test41_summary_by_regime.csv'));
writetable(T_by_region, fullfile(OUT.table_dir, 'test41_summary_by_analysis_region.csv'));
writetable(T_by_distance, fullfile(OUT.table_dir, 'test41_summary_by_distance.csv'));
writetable(T_by_side_regime, fullfile(OUT.table_dir, 'test41_summary_by_soft_hard_regime.csv'));
writetable(T_by_M_regime, fullfile(OUT.table_dir, 'test41_summary_by_M_regime.csv'));
writetable(T_model_cv, fullfile(OUT.table_dir, 'test41_corrector_cv_summary.csv'));
if ~isempty(T33_ref)
    writetable(T33_ref, fullfile(OUT.table_dir, 'test41_reference_test33_summary.csv'));
end

save(fullfile(OUT.data_dir, 'test41_mixedness_corrected_latest_kwave.mat'), ...
    'T_all','T_overall','T_by_side','T_by_roi','T_by_M','T_by_regime', ...
    'T_by_region','T_by_distance','T_model_cv','T33_ref','CFG','SRC','-v7.3');

%% Figures

plot_overall_ranking(T_overall, CFG, OUT);
plot_test33_style_soft_hard(T_by_side, T_by_side_regime, T_by_roi, T_by_M, CFG, OUT);
plot_mape_by_M_and_regime(T_by_M_regime, CFG, OUT);
plot_distance_curves(T_by_distance, CFG, OUT);
plot_correction_gain(T_all, CFG, OUT);
plot_test33_reference(T_overall, T33_ref, CFG, OUT);
if CFG.SaveAllMaps
    plot_condition_maps(T_all, CFG, OUT);
end

print_interpretation(T_overall, T_by_side, T_by_roi, T_by_region, T_model_cv, CFG);

fprintf('\nTables: %s\nFigures: %s\nTest 41 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Source and setup helpers

function OUT = make_output_dirs(root_dir)
OUT.root_dir = root_dir;
OUT.table_dir = fullfile(root_dir, 'tables');
OUT.figure_dir = fullfile(root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(root_dir, 'data');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end
end

function SRC = locate_sources(root_dir)
SRC = struct();
SRC.test40_patch_csv = fullfile(root_dir, 'outputs', ...
    'test_40_kwave_latest_model_comparison', 'tables', ...
    'test40_kwave_patch_level_model_comparison.csv');
SRC.test33_overall = fullfile(root_dir, 'outputs', ...
    'test_33_mixedness_aware_q_correction', 'tables', ...
    'test33_strategy_summary_overall_test.csv');
end

function validate_inputs(T, CFG)
needed = ["condition_key","case_name","field_regime","M","f0","map_iz","map_ix", ...
    "true_SWS","sws_pred","sws_signed_error_pct","strategy_name", ...
    "confidence","predicted_patch_purity","p_mixed","material_side", ...
    "roi_name","analysis_region"];
assert(all(ismember(needed, string(T.Properties.VariableNames))), ...
    'Test40 table is missing required variables.');
assert(all(ismember(CFG.BaseNewStrategies, unique(T.strategy_name))), ...
    'At least one requested Test41 base model is missing in Test40.');
assert(~any(ismember(["material_side","roi_name","analysis_region", ...
    "distance_to_interface_mm","true_SWS"], operational_feature_names())), ...
    'Oracle/evaluation variable leaked into operational features.');
end

function T = standardize_types(T)
string_vars = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "roi_name","material_side","region_label","distance_bin","purity_bin", ...
    "model_family","bundle_id","strategy_name","analysis_region"];
for v = string_vars
    if ismember(v, string(T.Properties.VariableNames))
        T.(v) = string(T.(v));
    end
end
if ~ismember("analysis_region", string(T.Properties.VariableNames))
    T.analysis_region = repmat("unknown", height(T), 1);
end
if ~ismember("p_mixed", string(T.Properties.VariableNames)) && ...
        ismember("predicted_mixed_probability", string(T.Properties.VariableNames))
    T.p_mixed = T.predicted_mixed_probability;
end
end

function T = attach_operational_references(T)
T.patch_key = T.condition_key + "__iz" + string(T.map_iz) + "__ix" + string(T.map_ix);
T.sws_theory_ref = lookup_strategy_sws(T, "T34_theory_baseline");
T.sws_local_ref = lookup_strategy_sws(T, "T34_local_baseline");
T.sws_old_mixedness_ref = lookup_strategy_sws(T, "T34_mixedness_logk_corrected");
T.log_sws_theory_ref = log(max(T.sws_theory_ref, 0.1));
T.log_sws_local_ref = log(max(T.sws_local_ref, 0.1));
end

function y = lookup_strategy_sws(T, strategy)
R = T(T.strategy_name == strategy, {'patch_key','sws_pred'});
[ok, loc] = ismember(T.patch_key, R.patch_key);
y = nan(height(T),1);
y(ok) = R.sws_pred(loc(ok));
end

%% Corrector

function [C, S_model] = train_apply_logk_corrector(B, base_strategy, CFG)
[X, feature_names] = build_operational_features(B);
assert_no_oracle_features(feature_names);
valid = all(isfinite(X),2) & isfinite(B.sws_pred) & B.sws_pred > 0 & ...
    isfinite(B.true_SWS) & B.true_SWS > 0;

base_k = 2*pi*B.f0 ./ B.sws_pred;
true_k = 2*pi*B.f0 ./ B.true_SWS;
target_delta = log(true_k) - log(base_k);

[G, condition_groups] = findgroups(B.condition_key);
fold_by_group = mod((1:numel(condition_groups))' - 1, CFG.KFoldByCondition) + 1;
fold_id = fold_by_group(G);
delta_pred = nan(height(B),1);

tree = templateTree('MinLeafSize', CFG.MinLeafSize, 'MaxNumSplits', 384);
fold_rows = table();
for f = 1:CFG.KFoldByCondition
    train = valid & fold_id ~= f;
    test = valid & fold_id == f;
    if nnz(train) < 100 || nnz(test) < 1, continue; end
    model = fitrensemble(X(train,:), target_delta(train), ...
        'Method','LSBoost','Learners',tree, ...
        'NumLearningCycles', CFG.NumLearningCycles, 'LearnRate', 0.05);
    delta_pred(test) = predict(model, X(test,:));
    row = table(string(base_strategy), f, nnz(train), nnz(test), ...
        mean(abs(target_delta(test)), 'omitnan'), mean(abs(delta_pred(test)-target_delta(test)), 'omitnan'), ...
        'VariableNames', {'base_strategy','fold_id','N_train','N_test', ...
        'mean_abs_delta_logk_label','mean_abs_delta_logk_residual'});
    fold_rows = concat_tables(fold_rows, row);
end

% Fallback model for any rows not covered by a fold.
miss = valid & ~isfinite(delta_pred);
if any(miss)
    model = fitrensemble(X(valid,:), target_delta(valid), ...
        'Method','LSBoost','Learners',tree, ...
        'NumLearningCycles', CFG.NumLearningCycles, 'LearnRate', 0.05);
    delta_pred(miss) = predict(model, X(miss,:));
end
delta_pred(~isfinite(delta_pred)) = 0;

gate = mixedness_gate(B);
base_logk = log(base_k);
sws_all = clip_sws(2*pi*B.f0 ./ exp(base_logk + delta_pred), CFG);
sws_gated = clip_sws(2*pi*B.f0 ./ exp(base_logk + gate .* delta_pred), CFG);

C1 = corrected_rows(B, base_strategy, "logk_all", sws_all, delta_pred, ones(height(B),1));
C2 = corrected_rows(B, base_strategy, "logk_mixed_gated", sws_gated, delta_pred, gate);
C = vertcat_compatible(C1, C2);
S_model = fold_rows;
end

function [X, names] = build_operational_features(B)
names = operational_feature_names();
q = B.q_pred;
q(~isfinite(q)) = 0.5;
purity = B.predicted_patch_purity;
purity(~isfinite(purity)) = 0.95;
p_mixed = B.p_mixed;
p_mixed(~isfinite(p_mixed)) = 1 - purity(~isfinite(p_mixed));
confidence = B.confidence;
confidence(~isfinite(confidence)) = 0.8;
theory = B.sws_theory_ref;
local = B.sws_local_ref;
theory(~isfinite(theory)) = B.sws_pred(~isfinite(theory));
local(~isfinite(local)) = B.sws_pred(~isfinite(local));
reg = regime_code(B.field_regime);
X = [ ...
    normalize_like(B.M, 3, 1), ...
    normalize_like(B.f0, 500, 100), ...
    reg, ...
    confidence, ...
    1-confidence, ...
    purity, ...
    1-purity, ...
    p_mixed, ...
    log(max(B.sws_pred, 0.1)), ...
    q, ...
    log(max(theory, 0.1)), ...
    log(max(local, 0.1)), ...
    log(max(B.sws_pred, 0.1)) - log(max(theory, 0.1)), ...
    log(max(B.sws_pred, 0.1)) - log(max(local, 0.1)), ...
    abs(B.sws_pred - theory)./max(theory,0.25), ...
    abs(B.sws_pred - local)./max(local,0.25), ...
    abs(local - theory)./max(theory,0.25) ...
    ];
X = double(X);
X(~isfinite(X)) = 0;
end

function names = operational_feature_names()
names = ["M_norm","f0_norm","field_regime_code","confidence", ...
    "risk","predicted_patch_purity","predicted_impurity","p_mixed", ...
    "log_sws_base","q_pred","log_sws_theory_ref","log_sws_local_ref", ...
    "log_base_minus_theory", ...
    "log_base_minus_local","rel_base_theory_disagreement", ...
    "rel_base_local_disagreement","rel_local_theory_disagreement"];
end

function gate = mixedness_gate(B)
purity = B.predicted_patch_purity;
purity(~isfinite(purity)) = 0.95;
p = B.p_mixed;
p(~isfinite(p)) = 1 - purity(~isfinite(p));
risk = 1 - B.confidence;
risk(~isfinite(risk)) = 0.2;
gate = max([p, 1-purity, 0.5*risk], [], 2);
gate = min(max(gate, 0), 1);
end

function C = corrected_rows(B, base_strategy, variant, sws, delta_logk, gate)
C = B;
C.base_strategy_name = repmat(base_strategy, height(B), 1);
C.correction_variant = repmat(variant, height(B), 1);
C.model_family = repmat("T41_mixedness_corrected_latest", height(B), 1);
C.bundle_id = repmat("test41_grouped_cv", height(B), 1);
C.strategy_name = repmat("T41_" + base_strategy + "_" + variant, height(B), 1);
C.sws_pred = sws;
C.sws_signed_error_pct = 100*(sws - C.true_SWS) ./ C.true_SWS;
C.sws_abs_error_pct = abs(C.sws_signed_error_pct);
C.high_error10 = C.sws_abs_error_pct > 10;
C.high_error20 = C.sws_abs_error_pct > 20;
C.mixedness_gate = gate;
C.delta_logk_pred = delta_logk;
C.fold_id = nan(height(B),1);
end

function assert_no_oracle_features(names)
bad = ["true","oracle","material","side","roi","distance", ...
    "error","signed","abs_error","high_error","cs_true","sws_true"];
for b = bad
    hit = names(contains(lower(names), b));
    assert(isempty(hit), 'Oracle/evaluation feature leaked into corrector: %s', strjoin(hit, ', '));
end
end

%% Summaries

function S = summarize_long(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
n = splitapply(@numel, T.sws_abs_error_pct, G);
mape = splitapply(@(x) mean(x,'omitnan'), T.sws_abs_error_pct, G);
medae = splitapply(@(x) median(x,'omitnan'), T.sws_abs_error_pct, G);
signed = splitapply(@(x) mean(x,'omitnan'), T.sws_signed_error_pct, G);
medsigned = splitapply(@(x) median(x,'omitnan'), T.sws_signed_error_pct, G);
under = splitapply(@(x) 100*mean(x < 0,'omitnan'), T.sws_signed_error_pct, G);
he10 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error10, G);
he20 = splitapply(@(x) 100*mean(x,'omitnan'), T.high_error20, G);
mean_sws = splitapply(@(x) mean(x,'omitnan'), T.sws_pred, G);
std_sws = splitapply(@std_omitnan, T.sws_pred, G);
S = [groups table(n,mape,medae,signed,medsigned,under,he10,he20,mean_sws,std_sws, ...
    'VariableNames', {'N','MAPE_pct','median_abs_error_pct','mean_signed_error_pct', ...
    'median_signed_error_pct','underestimate_pct','high_error10_pct','high_error20_pct', ...
    'mean_sws_pred','std_sws_pred'})];
end

function T33 = load_test33_reference(SRC)
if exist(SRC.test33_overall,'file') ~= 2
    T33 = table();
    return;
end
T33 = readtable(SRC.test33_overall);
for v = ["strategy_name","geometry"]
    if ismember(v, string(T33.Properties.VariableNames))
        T33.(v) = string(T33.(v));
    end
end
end

%% Plots

function plot_overall_ranking(T, CFG, OUT)
keep = select_plot_strategies(T.strategy_name, CFG);
T = T(ismember(T.strategy_name, keep),:);
T = sortrows(T, 'MAPE_pct', 'ascend');
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
barh(categorical(T.strategy_name), T.MAPE_pct);
xlabel('MAPE (%)'); grid on;
title('Test 41 k-Wave MAPE ranking','Interpreter','none','FontWeight','normal');
set(gca,'TickLabelInterpreter','none');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_overall_model_ranking.png'));
end

function plot_test33_style_soft_hard(T_side, T_side_regime, T_roi, T_M, CFG, OUT)
keep = select_plot_strategies(T_side.strategy_name, CFG);
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 20]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); grouped_bar(ax, T_side, "material_side", keep, "MAPE_pct");
ylabel(ax,'MAPE (%)'); title(ax,'MAPE by side'); legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');

ax = nexttile(tl); grouped_bar(ax, T_side, "material_side", keep, "mean_signed_error_pct");
yline(ax,0,'k-'); ylabel(ax,'mean signed error (%)'); title(ax,'Signed error by side'); legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');

roi_keep = ["background_roi","inclusion_roi"];
R = T_roi(ismember(T_roi.roi_name, roi_keep),:);
ax = nexttile(tl); grouped_bar(ax, R, "roi_name", keep, "MAPE_pct");
ylabel(ax,'MAPE (%)'); title(ax,'Core ROI MAPE'); legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');

H = T_side_regime(T_side_regime.material_side == "hard",:);
ax = nexttile(tl); hold(ax,'on');
regimes = unique(H.field_regime, 'stable');
for s = keep(:)'
    y = nan(numel(regimes),1);
    for i = 1:numel(regimes)
        idx = H.strategy_name == s & H.field_regime == regimes(i);
        if any(idx), y(i) = mean(H.MAPE_pct(idx), 'omitnan'); end
    end
    plot(ax, 1:numel(regimes), y, '-o', 'LineWidth', 1.0, 'DisplayName', s);
end
set(ax,'XTick',1:numel(regimes),'XTickLabel',regimes,'TickLabelInterpreter','none');
xtickangle(ax,25); ylabel(ax,'MAPE (%)'); title(ax,'Hard side vs field regime');
legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_soft_hard_error_diagnostics.png'));

fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13]);
ax = axes(fig); hold(ax,'on');
Ms = unique(T_M.M, 'stable');
for s = keep(:)'
    y = nan(numel(Ms),1);
    for i = 1:numel(Ms)
        idx = T_M.strategy_name == s & T_M.M == Ms(i);
        if any(idx), y(i) = mean(T_M.MAPE_pct(idx), 'omitnan'); end
    end
    plot(ax, Ms, y, '-o', 'LineWidth', 1.0, 'DisplayName', s);
end
xlabel(ax,'REQ M'); ylabel(ax,'MAPE (%)'); title(ax,'Overall MAPE vs M');
legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_mape_vs_M.png'));
end

function plot_mape_by_M_and_regime(T, CFG, OUT)
keep = select_plot_strategies(T.strategy_name, CFG);
T = T(ismember(T.strategy_name, keep),:);
regimes = unique(T.field_regime, 'stable');
fig = figure('Color','w','Units','centimeters','Position',[1 1 32 18]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for r = 1:numel(regimes)
    ax = nexttile(tl); hold(ax,'on');
    X = T(T.field_regime == regimes(r),:);
    Ms = unique(X.M, 'stable');
    for s = keep(:)'
        y = nan(numel(Ms),1);
        for i = 1:numel(Ms)
            idx = X.strategy_name == s & X.M == Ms(i);
            if any(idx), y(i) = mean(X.MAPE_pct(idx), 'omitnan'); end
        end
        plot(ax, Ms, y, '-o', 'DisplayName', s);
    end
    title(ax, regimes(r), 'Interpreter','none'); xlabel(ax,'M'); ylabel(ax,'MAPE (%)'); grid(ax,'on');
end
legend(tl.Children(end), 'Location','bestoutside','Interpreter','none');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_mape_by_M_and_regime.png'));
end

function plot_distance_curves(T_dist, CFG, OUT)
keep = select_plot_strategies(T_dist.strategy_name, CFG);
order = ["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13]);
ax = axes(fig); hold(ax,'on');
for s = keep(:)'
    y = nan(size(order));
    for i = 1:numel(order)
        idx = T_dist.strategy_name == s & string(T_dist.distance_bin) == order(i);
        if any(idx), y(i) = mean(T_dist.MAPE_pct(idx), 'omitnan'); end
    end
    plot(ax, 1:numel(order), y, '-o', 'LineWidth', 1.0, 'DisplayName', s);
end
set(ax,'XTick',1:numel(order),'XTickLabel',order,'TickLabelInterpreter','none');
xtickangle(ax,25); xlabel(ax,'Distance to interface'); ylabel(ax,'MAPE (%)');
title(ax,'MAPE vs distance to interface'); legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_mape_vs_distance.png'));
end

function plot_correction_gain(T, CFG, OUT)
rows = table();
for base = CFG.BaseNewStrategies(:)'
    B = T(T.strategy_name == base,:);
    for variant = ["logk_all","logk_mixed_gated"]
        s = "T41_" + base + "_" + variant;
        C = T(T.strategy_name == s,:);
        if isempty(B) || isempty(C), continue; end
        [ok,loc] = ismember(C.patch_key, B.patch_key);
        if ~all(ok), continue; end
        gain = B.sws_abs_error_pct(loc(ok)) - C.sws_abs_error_pct(ok);
        row = table(repmat(s, sum(ok), 1), C.material_side(ok), C.analysis_region(ok), gain, ...
            'VariableNames', {'strategy_name','material_side','analysis_region','mape_gain_points'});
        rows = concat_tables(rows, row);
    end
end
if isempty(rows), return; end
[G,groups] = findgroups(rows(:, {'strategy_name','analysis_region'}));
S = [groups table(splitapply(@numel, rows.mape_gain_points, G), ...
    splitapply(@(x) mean(x,'omitnan'), rows.mape_gain_points, G), ...
    'VariableNames', {'N','mean_mape_gain_points'})];
writetable(S, fullfile(OUT.table_dir, 'test41_correction_gain_by_region.csv'));
keep = select_plot_strategies(S.strategy_name, CFG);
S = S(ismember(S.strategy_name, keep),:);
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 14]);
grouped_bar(axes(fig), S, "analysis_region", keep, "mean_mape_gain_points");
yline(0,'k-'); ylabel('MAPE gain points (baseline - corrected)'); grid on;
title('Correction gain by region','FontWeight','normal'); legend('Location','bestoutside','Interpreter','none');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_correction_gain_by_region.png'));
end

function plot_test33_reference(T_overall, T33_ref, CFG, OUT)
if isempty(T33_ref), return; end
common = ["local_baseline","test30_theory_region_levels","mixedness_logk_corrected", ...
    "mixedness_q_candidate_selector"];
T34_names = "T34_" + common;
K = T_overall(ismember(T_overall.strategy_name, T34_names),:);
if isempty(K), return; end
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 13]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile; hold(ax,'on');
for i = 1:numel(common)
    idx = T33_ref.strategy_name == common(i);
    if any(idx)
        bar(ax, i, mean(T33_ref.MAPE(idx), 'omitnan'), 0.65);
    end
end
set(ax,'XTick',1:numel(common),'XTickLabel',common,'TickLabelInterpreter','none');
xtickangle(ax,35); ylabel(ax,'MAPE (%)'); title(ax,'Test33 synthetic held-out reference'); grid(ax,'on');
ax = nexttile; hold(ax,'on');
for i = 1:numel(common)
    idx = K.strategy_name == "T34_" + common(i);
    if any(idx)
        bar(ax, i, mean(K.MAPE_pct(idx), 'omitnan'), 0.65);
    end
end
set(ax,'XTick',1:numel(common),'XTickLabel',"T34_"+common,'TickLabelInterpreter','none');
xtickangle(ax,35); ylabel(ax,'MAPE (%)'); title(ax,'Same old stack on k-Wave/Test40'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir, 'test41_test33_vs_kwave_reference.png'));
end

function plot_condition_maps(T, CFG, OUT)
keep = select_plot_strategies(T.strategy_name, CFG);
conds = unique(T(:, {'case_name','field_regime','M'}), 'rows', 'stable');
if isfinite(CFG.MaxMapConditions) && height(conds) > CFG.MaxMapConditions
    conds = conds(1:CFG.MaxMapConditions,:);
end
fprintf('Saving %d Test41 map panels under %s.\n', height(conds), OUT.map_dir);
for ci = 1:height(conds)
    idxc = T.case_name == conds.case_name(ci) & ...
        T.field_regime == conds.field_regime(ci) & T.M == conds.M(ci);
    Xc = T(idxc,:);
    if isempty(Xc), continue; end
    base = Xc(Xc.strategy_name == Xc.strategy_name(1),:);
    [true_map,nz,nx] = rows_to_grid(base, base.true_SWS);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 36 24]);
    tl = tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    plot_map(nexttile(tl), true_map, 'True SWS');
    plot_map(nexttile(tl), rows_to_grid(base, base.distance_to_interface_mm, nz, nx), 'Distance mm');
    plot_map(nexttile(tl), rows_to_grid(base, base.predicted_patch_purity, nz, nx), 'Predicted purity');
    plot_map(nexttile(tl), rows_to_grid(base, base.p_mixed, nz, nx), 'P(mixed)');
    for si = 1:numel(keep)
        s = keep(si);
        X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_pred, nz, nx), s + ' SWS');
    end
    err_keep = keep(1:min(4,numel(keep)));
    for si = 1:numel(err_keep)
        s = err_keep(si);
        X = Xc(Xc.strategy_name == s,:);
        if isempty(X), continue; end
        plot_map(nexttile(tl), rows_to_grid(X, X.sws_abs_error_pct, nz, nx), s + ' err %');
    end
    title(tl, sprintf('%s | %s | M%d', conds.case_name(ci), conds.field_regime(ci), conds.M(ci)), ...
        'Interpreter','none');
    outdir = fullfile(OUT.map_dir, sanitize(conds.case_name(ci)), sanitize(conds.field_regime(ci)), ...
        "M" + string(conds.M(ci)));
    if exist(outdir,'dir') ~= 7, mkdir(outdir); end
    export_fig(fig, fullfile(outdir, sprintf('test41_map__%s__%s__M%d.png', ...
        sanitize(conds.case_name(ci)), sanitize(conds.field_regime(ci)), conds.M(ci))));
end
end

function grouped_bar(ax, T, group_var, strategies, metric)
groups = unique(string(T.(group_var)), 'stable');
Y = nan(numel(groups), numel(strategies));
for i = 1:numel(groups)
    for j = 1:numel(strategies)
        idx = string(T.(group_var)) == groups(i) & T.strategy_name == strategies(j);
        if any(idx), Y(i,j) = mean(T.(metric)(idx), 'omitnan'); end
    end
end
bar(ax, categorical(groups), Y);
set(ax,'TickLabelInterpreter','none');
end

function keep = select_plot_strategies(strategies, CFG)
available = unique(string(strategies), 'stable');
preferred = string(CFG.MainPlotStrategies(:));
keep = preferred(ismember(preferred, available));
if isempty(keep)
    keep = available(1:min(8,numel(available)));
end
end

function print_interpretation(T_overall, T_side, T_roi, T_region, T_model_cv, CFG)
T = sortrows(T_overall, 'MAPE_pct', 'ascend');
fprintf('\n================ Test 41 summary ================\n');
disp(T(1:min(14,height(T)), {'model_family','strategy_name','N','MAPE_pct', ...
    'mean_signed_error_pct','high_error20_pct'}));
fprintf('Best overall: %s (MAPE %.2f%%, signed %.2f%%, HE20 %.2f%%).\n', ...
    T.strategy_name(1), T.MAPE_pct(1), T.mean_signed_error_pct(1), T.high_error20_pct(1));
for base = CFG.BaseNewStrategies(:)'
    b = T_overall(T_overall.strategy_name == base,:);
    c = T_overall(T_overall.strategy_name == "T41_" + base + "_logk_mixed_gated",:);
    a = T_overall(T_overall.strategy_name == "T41_" + base + "_logk_all",:);
    if ~isempty(b) && ~isempty(c)
        fprintf('%s mixed-gated delta: %.2f MAPE points; logk-all delta: %.2f points.\n', ...
            base, c.MAPE_pct-b.MAPE_pct, a.MAPE_pct-b.MAPE_pct);
    end
end
hard = sortrows(T_side(T_side.material_side=="hard",:), 'MAPE_pct', 'ascend');
soft = sortrows(T_side(T_side.material_side=="soft",:), 'MAPE_pct', 'ascend');
fprintf('Best hard side: %s (MAPE %.2f%%).\n', hard.strategy_name(1), hard.MAPE_pct(1));
fprintf('Best soft side: %s (MAPE %.2f%%).\n', soft.strategy_name(1), soft.MAPE_pct(1));
iface = sortrows(T_region(T_region.analysis_region=="interface_0_2mm",:), 'MAPE_pct', 'ascend');
fprintf('Best interface 0-2mm: %s (MAPE %.2f%%, HE20 %.2f%%).\n', ...
    iface.strategy_name(1), iface.MAPE_pct(1), iface.high_error20_pct(1));
fprintf('Corrector CV rows: %d folds logged.\n', height(T_model_cv));
fprintf('=================================================\n');
end

%% Generic utilities

function y = clip_sws(y, CFG)
y = min(max(y, CFG.PhysicalRange(1)), CFG.PhysicalRange(2));
end

function y = normalize_like(x, mu, sigma)
y = (double(x) - mu) ./ sigma;
end

function code = regime_code(regime)
regime = string(regime);
levels = ["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
code = zeros(numel(regime),1);
for i = 1:numel(levels)
    code(regime == levels(i)) = i;
end
code = normalize_like(code, 2.5, 1.2);
end

function T = vertcat_compatible(varargin)
vars = strings(1,0);
for i = 1:nargin
    vars = [vars string(varargin{i}.Properties.VariableNames)]; %#ok<AGROW>
end
vars = unique(vars, 'stable');
string_vars = ["condition_key","case_name","field_regime","geometry","geometry_type", ...
    "roi_name","material_side","region_label","distance_bin","purity_bin", ...
    "model_family","bundle_id","strategy_name","analysis_region", ...
    "base_strategy_name","correction_variant","patch_key"];
logical_vars = ["high_error10","high_error20"];
parts = cell(nargin,1);
for i = 1:nargin
    A = varargin{i};
    for v = vars
        if ismember(v, string(A.Properties.VariableNames)), continue; end
        if ismember(v, string_vars)
            A.(v) = repmat("unknown", height(A), 1);
        elseif ismember(v, logical_vars)
            A.(v) = false(height(A), 1);
        else
            A.(v) = nan(height(A), 1);
        end
    end
    parts{i} = A(:, cellstr(vars));
end
T = vertcat(parts{:});
end

function T = concat_tables(T, row)
if isempty(T)
    T = row;
else
    T = vertcat_compatible(T, row);
end
end

function tf = env_true(name, default_value)
v = lower(strtrim(string(getenv(name))));
if v == "", tf = default_value; return; end
tf = ismember(v, ["1","true","yes","y","on"]);
end

function x = env_number(name, default_value)
v = strtrim(string(getenv(name)));
if v == "", x = default_value; return; end
x = str2double(v);
if ~isfinite(x), x = default_value; end
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end

function s = std_omitnan(x)
s = std(x, 'omitnan');
end

function [Z,nz,nx] = rows_to_grid(T, values, nz, nx)
if nargin < 3
    nz = max(T.map_iz); nx = max(T.map_ix);
end
Z = nan(nz,nx);
Z(sub2ind([nz,nx], T.map_iz, T.map_ix)) = values;
end

function plot_map(ax, Z, ttl)
imagesc(ax, Z); axis(ax,'image'); axis(ax,'off'); colorbar(ax);
title(ax, ttl, 'Interpreter','none','FontWeight','normal','FontSize',7);
end

function export_fig(fig, file)
drawnow;
exportgraphics(fig, file, 'Resolution', 220);
close(fig);
end

function s = sanitize(s)
s = regexprep(string(s), '[^A-Za-z0-9_=-]+', '_');
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file,'w');
assert(fid>0, 'Could not write config file: %s', file);
fprintf(fid,'%s\n',txt);
fclose(fid);
end
