%% analyze_test_56_meff_operational_correction.m
% Test 56: operational M-effective residual correction.
%
% This script trains a lightweight second-pass correction on top of the frozen
% q_spectrum_plus_composition output. It uses only operational first-pass
% quantities as predictors. True SWS/k are used only as residual labels and
% evaluation targets.

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
write_config_json(CFG, fullfile(OUT.root_dir,'test56_configuration.json'));

fprintf('\nTest 56: M-effective operational correction\n');
fprintf('Mode: %s | source: %s\n', CFG.Mode, CFG.SourceCsv);
fprintf('Base model: %s | no oracle variables are used as predictors.\n', CFG.BaseModel);

T0 = readtable(CFG.SourceCsv);
T0 = T0(string(T0.model_name) == CFG.BaseModel, :);
assert(~isempty(T0), 'No rows found for base model %s.', CFG.BaseModel);
T0 = add_operational_meff_features(T0, CFG);
assert_no_forbidden_predictors(CFG.FeatureVars);

if CFG.ValidateOnly
    T0 = sample_rows(T0, min(height(T0), CFG.ValidateRows), CFG.RandomSeed);
end

[train_mask, test_mask] = condition_split(T0, CFG);
MODEL = train_meff_residual(T0, train_mask, CFG);
T_pred = apply_strategies(T0(test_mask,:), MODEL, CFG);

T_overall = summarize_predictions(T_pred, "strategy_name");
T_by_case = summarize_predictions(T_pred, ["strategy_name","case_id"]);
T_by_freq = summarize_predictions(T_pred, ["strategy_name","f0"]);
T_by_regime = summarize_predictions(T_pred, ["strategy_name","field_regime_ood"]);
T_by_roi = summarize_predictions(T_pred(T_pred.roi_region ~= "other",:), ...
    ["strategy_name","roi_region"]);
T_by_case_roi = summarize_predictions(T_pred(T_pred.roi_region ~= "other",:), ...
    ["strategy_name","case_id","roi_region"]);
T_hard = summarize_predictions(T_pred(contains(string(T_pred.material_region),"hard"),:), ...
    ["strategy_name","case_id","roi_region"]);

writetable(remove_cell_columns(T_pred), fullfile(OUT.table_dir,'test56_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test56_strategy_summary_overall.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test56_strategy_summary_by_case.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir,'test56_strategy_summary_by_frequency.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test56_strategy_summary_by_regime.csv'));
writetable(T_by_roi, fullfile(OUT.table_dir,'test56_strategy_summary_by_roi.csv'));
writetable(T_by_case_roi, fullfile(OUT.table_dir,'test56_strategy_summary_by_case_roi.csv'));
writetable(T_hard, fullfile(OUT.table_dir,'test56_hard_region_summary.csv'));

safe_plot(@() plot_test56_summary(T_overall, T_by_case_roi, T_by_freq, OUT), 'summary plots');
save(fullfile(OUT.data_dir,'test56_meff_correction_model.mat'), ...
    'CFG','MODEL','T_overall','T_by_case','T_by_freq','T_by_regime', ...
    'T_by_roi','T_by_case_roi','-v7.3');

print_summary(T_overall, T_by_case_roi, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 56 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST56_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST56_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST56_VALIDATE_ONLY', false);
CFG.RandomSeed = 56001;
CFG.SourceCsv = env_string('ADAPTIVE_REQ_TEST56_SOURCE_CSV', ...
    fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep', ...
    'tables','test55_patch_level_results.csv'));
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST56_BASE_MODEL', ...
    "q_spectrum_plus_composition");
CFG.TrainFraction = env_number('ADAPTIVE_REQ_TEST56_TRAIN_FRACTION', 0.70);
CFG.CsGuess = env_number('ADAPTIVE_REQ_TEST56_CS_GUESS', 3.0);
CFG.TreeLearners = env_number('ADAPTIVE_REQ_TEST56_TREE_LEARNERS', 180);
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST56_MIN_LEAF_SIZE', 10);
CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST56_MAX_TRAIN_ROWS', ...
    ternary(mode=="quick", 60000, 250000));
CFG.ValidateRows = env_number('ADAPTIVE_REQ_TEST56_VALIDATE_ROWS', 20000);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST56_USE_PARALLEL', false);
CFG.HighPurityThreshold = env_number('ADAPTIVE_REQ_TEST56_HIGH_PURITY_THRESHOLD', 0.90);
CFG.ConservativeDeltaLimit = env_number('ADAPTIVE_REQ_TEST56_CONSERVATIVE_DELTA_LIMIT', 0.12);
CFG.FeatureVars = ["q_pred";"sws_pred";"log_sws_pred";"k_pred_firstpass"; ...
    "logk_pred_firstpass";"M_eff_pred";"lambda_pred_m";"window_length_m"; ...
    "predicted_patch_purity";"p_mixed";"p_strong_mixed";"M";"f0"; ...
    "dx";"dz";"field_code";"case_family_code"];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_56_meff_operational_correction','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_56_meff_operational_correction');
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

%% Model

function T = add_operational_meff_features(T, CFG)
T.sws_pred = clamp(T.sws_pred, 0.2, 20);
T.q_pred = clamp(T.q_pred, 0, 1);
T.k_pred_firstpass = 2*pi*T.f0 ./ T.sws_pred;
T.logk_pred_firstpass = log(T.k_pred_firstpass);
T.log_sws_pred = log(T.sws_pred);
T.lambda_pred_m = T.sws_pred ./ T.f0;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.M_eff_pred = T.window_length_m ./ T.lambda_pred_m;
T.field_code = categorical_code(string(T.field_regime_ood));
T.case_family_code = categorical_code(string(T.case_family));
end

function MODEL = train_meff_residual(T, train_mask, CFG)
target = log(T.k_true) - log(T.k_pred_firstpass);
valid = train_mask & all(isfinite(table2array(T(:, cellstr(CFG.FeatureVars)))),2) & isfinite(target);
idx = find(valid);
if CFG.MaxTrainRows > 0 && numel(idx) > CFG.MaxTrainRows
    rng(CFG.RandomSeed);
    idx = idx(randperm(numel(idx), CFG.MaxTrainRows));
end
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
MODEL = struct();
MODEL.features = CFG.FeatureVars(:);
MODEL.target_definition = "delta_logk_true = log(k_true) - log(k_firstpass)";
MODEL.residual = fitrensemble(T(idx, cellstr(CFG.FeatureVars)), target(idx), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners, ...
    'Learners',template,'Options',opts);
MODEL.n_train = numel(idx);
fprintf('Trained M-eff residual model on %d rows.\n', MODEL.n_train);
end

function T_out = apply_strategies(T, MODEL, CFG)
delta = predict(MODEL.residual, T(:, cellstr(MODEL.features)));
delta = clamp(delta, -0.45, 0.45);

base = make_strategy(T, "base_q_spectrum_plus_composition", T.sws_pred, zeros(height(T),1), false(height(T),1));
allcorr = corrected_strategy(T, "meff_delta_logk_all", delta, true(height(T),1));
highpur = corrected_strategy(T, "meff_delta_logk_high_predicted_purity", delta, ...
    T.predicted_patch_purity >= CFG.HighPurityThreshold);
conservative = corrected_strategy(T, "meff_delta_logk_conservative", delta, ...
    T.predicted_patch_purity >= CFG.HighPurityThreshold & abs(delta) <= CFG.ConservativeDeltaLimit);

oracle_sws = apply_delta_logk(T, delta);
base_abs = abs(100*(T.sws_pred - T.true_SWS)./T.true_SWS);
oracle_abs = abs(100*(oracle_sws - T.true_SWS)./T.true_SWS);
oracle_mask = oracle_abs < base_abs;
oracle = corrected_strategy(T, "oracle_apply_if_improves", delta, oracle_mask);
oracle.diagnostic_only = true(height(oracle),1);

base.diagnostic_only = false(height(base),1);
allcorr.diagnostic_only = false(height(allcorr),1);
highpur.diagnostic_only = false(height(highpur),1);
conservative.diagnostic_only = false(height(conservative),1);
T_out = [base; allcorr; highpur; conservative; oracle];
end

function R = corrected_strategy(T, name, delta, mask)
sws = T.sws_pred;
sws(mask) = apply_delta_logk(T(mask,:), delta(mask));
R = make_strategy(T, name, sws, delta, mask);
end

function sws = apply_delta_logk(T, delta)
k_corr = T.k_pred_firstpass .* exp(delta);
sws = 2*pi*T.f0 ./ k_corr;
sws = clamp(sws, 0.2, 20);
end

function R = make_strategy(T, name, sws, delta, mask)
keep = ["dataset","condition_key","case_id","case_family","seen_status", ...
    "field_regime","field_regime_ood","f0","M","dx","dz", ...
    "true_SWS","k_true","patch_purity","purity_bin","distance_to_boundary_mm", ...
    "distance_bin","predicted_patch_purity","p_mixed","p_strong_mixed", ...
    "material_region","roi_region","q_pred","sws_pred","M_eff_pred"];
keep = intersect(keep, string(T.Properties.VariableNames), 'stable');
R = T(:, cellstr(keep));
R.strategy_name = repmat(string(name), height(T), 1);
R.sws_corrected = sws;
R.delta_logk_pred = delta;
R.was_corrected = mask;
R.correction_pct = 100*(sws - T.sws_pred)./T.sws_pred;
R.sws_signed_error_pct = 100*(sws - T.true_SWS)./T.true_SWS;
R.sws_abs_error_pct = abs(R.sws_signed_error_pct);
R.high_error10 = R.sws_abs_error_pct > 10;
R.high_error20 = R.sws_abs_error_pct > 20;
end

%% Summaries and plots

function S = summarize_predictions(T, group_vars)
if isstring(group_vars), group_vars = cellstr(group_vars); end
if isempty(T), S = table(); return; end
[G, groups] = findgroups(T(:, group_vars));
rows = cell(max(G),1);
for gi = 1:max(G)
    X = T(G==gi,:);
    rows{gi} = table(height(X), mean(X.sws_abs_error_pct,'omitnan'), ...
        median(X.sws_abs_error_pct,'omitnan'), mean(X.sws_signed_error_pct,'omitnan'), ...
        median(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), 100*mean(X.sws_signed_error_pct < 0,'omitnan'), ...
        100*mean(X.was_corrected,'omitnan'), mean(abs(X.correction_pct),'omitnan'), ...
        mean(X.sws_corrected,'omitnan'), mean(X.true_SWS,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','median_abs_error_pct', ...
        'mean_signed_error_pct','median_signed_error_pct','high_error10_pct', ...
        'high_error20_pct','underestimate_pct','corrected_pct', ...
        'mean_abs_correction_pct','mean_pred_sws','mean_true_sws'});
end
S = [groups vertcat(rows{:})];
end

function plot_test56_summary(T_overall, T_roi, T_freq, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
X = sortrows(T_overall,'MAPE_pct','ascend');
barh_labels(ax, pretty_strategy(X.strategy_name), X.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'M-effective correction ranking','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for s = unique(T_freq.strategy_name,'stable')'
    if s == "oracle_apply_if_improves", continue; end
    F = sortrows(T_freq(T_freq.strategy_name==s,:), 'f0');
    plot(ax,F.f0,F.MAPE_pct,'-o','DisplayName',pretty_strategy(s));
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'MAPE vs frequency','FontWeight','normal');
legend(ax,'Location','bestoutside','Interpreter','none'); grid(ax,'on');
ax = nexttile(tl);
R = T_roi(ismember(T_roi.roi_region, ["homogeneous_center_8mm","soft_core_gt8mm","hard_core_gt4mm"]),:);
R = sortrows(R,'MAPE_pct','descend');
barh_labels(ax, pretty_strategy(R.strategy_name)+" | "+string(R.case_id)+" | "+pretty_region(R.roi_region), R.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Core ROI errors','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl);
R = T_roi(contains(string(T_roi.roi_region),"interface"),:);
R = aggregate_for_plot(R, ["strategy_name","roi_region"]);
R = sortrows(R,'MAPE_pct','descend');
barh_labels(ax, pretty_strategy(R.strategy_name)+" | "+pretty_region(R.roi_region), R.MAPE_pct);
xlabel(ax,'MAPE (%)'); title(ax,'Interface errors','FontWeight','normal'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test56_meff_correction_summary.png'));
end

function print_summary(T_overall, T_roi, OUT)
fprintf('\n================ Test 56 summary ================\n');
disp(T_overall(:, {'strategy_name','N','MAPE_pct','mean_signed_error_pct','high_error20_pct','corrected_pct'}));
B = T_overall(T_overall.strategy_name=="base_q_spectrum_plus_composition",:);
X = T_overall(~ismember(T_overall.strategy_name, ["base_q_spectrum_plus_composition","oracle_apply_if_improves"]),:);
if ~isempty(B) && ~isempty(X)
    [~,i] = min(X.MAPE_pct);
    fprintf('Best operational strategy: %s, MAPE %.2f%% (baseline %.2f%%).\n', ...
        X.strategy_name(i), X.MAPE_pct(i), B.MAPE_pct);
end
H = T_roi(contains(string(T_roi.roi_region),"hard_core"),:);
if ~isempty(H)
    H = sortrows(H,'MAPE_pct','ascend');
    fprintf('Best hard-core row: %s | %s, MAPE %.2f%%, signed %.2f%%.\n', ...
        string(H.strategy_name(1)), string(H.case_id(1)), H.MAPE_pct(1), H.mean_signed_error_pct(1));
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

%% Helpers

function barh_labels(ax, labels, values)
y = 1:numel(values);
barh(ax, y, values);
set(ax, 'YTick', y, 'YTickLabel', cellstr(string(labels)));
end

function [train_mask, test_mask] = condition_split(T, CFG)
keys = unique(string(T.condition_key), 'stable');
rng(CFG.RandomSeed);
keys = keys(randperm(numel(keys)));
ntrain = max(1, round(CFG.TrainFraction*numel(keys)));
train_keys = keys(1:ntrain);
train_mask = ismember(string(T.condition_key), train_keys);
test_mask = ~train_mask;
end

function T = sample_rows(T, n, seed)
rng(seed);
idx = randperm(height(T), n);
T = T(idx,:);
end

function code = categorical_code(x)
[~,~,ic] = unique(string(x), 'stable');
code = double(ic);
end

function assert_no_forbidden_predictors(features)
low = lower(string(features));
bad_exact = ["patch_purity","true_sws","k_true","q_oracle","material_region", ...
    "roi_region","distance_to_boundary_mm","distance_bin"];
bad_patterns = ["oracle","true","error","high_error"];
hit = low(ismember(low, bad_exact));
assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
for p = bad_patterns
    hit = low(contains(low,p));
    assert(isempty(hit), 'Forbidden operational predictor: %s', strjoin(hit, ', '));
end
end

function A = aggregate_for_plot(T, group_vars)
if isempty(T), A = table(); return; end
if isstring(group_vars), group_vars = cellstr(group_vars); end
[G, groups] = findgroups(T(:, group_vars));
A = groups;
A.N = splitapply(@(x) sum(x,'omitnan'), T.N, G);
A.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), T.MAPE_pct, T.N, G);
end

function y = weighted_mean(x,w)
ok = isfinite(x) & isfinite(w) & w > 0;
if ~any(ok), y = NaN; else, y = sum(x(ok).*w(ok))/sum(w(ok)); end
end

function y = clamp(x, lo, hi), y = min(max(x,lo),hi); end
function y = ternary(tf, a, b), if tf, y = a; else, y = b; end, end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
end

function s = pretty_strategy(x)
s = string(x);
s(s=="base_q_spectrum_plus_composition") = "Base q spectrum + composition";
s(s=="meff_delta_logk_all") = "M-eff residual applied everywhere";
s(s=="meff_delta_logk_high_predicted_purity") = "M-eff residual in high predicted purity";
s(s=="meff_delta_logk_conservative") = "Conservative M-eff residual";
s(s=="oracle_apply_if_improves") = "Oracle apply-if-improves";
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
    warning('Test56:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
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
