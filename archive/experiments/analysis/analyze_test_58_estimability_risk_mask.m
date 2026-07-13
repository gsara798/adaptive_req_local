%% analyze_test_58_estimability_risk_mask.m
% Test 58: estimability / confidence risk mask.
%
% Trains an operational high-error risk mask on clean Test55 predictions and
% optionally evaluates transfer to Test54 realistic readout outputs. No SWS
% correction is applied here; the output is reliability/estimability.

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
write_config_json(CFG, fullfile(OUT.root_dir,'test58_configuration.json'));

fprintf('\nTest 58: estimability / high-error risk mask\n');
fprintf('Clean source: %s\n', CFG.CleanCsv);
fprintf('Realistic source: %s\n', CFG.RealisticCsv);
fprintf('No correction. Truth/error labels are used only for training/evaluation.\n');

Tclean = load_source(CFG.CleanCsv, CFG.BaseModel, "clean_synthetic");
Tclean = add_operational_risk_features(Tclean, CFG);
assert_no_forbidden_predictors(CFG.FeatureVars);
if CFG.ValidateOnly
    Tclean = sample_rows(Tclean, min(height(Tclean), CFG.ValidateRows), CFG.RandomSeed);
end

[train_mask, test_mask] = condition_split(Tclean, CFG);
MODEL = train_risk_models(Tclean, train_mask, CFG);
T_eval_clean = apply_risk_model(Tclean(test_mask,:), MODEL, "clean_synthetic", CFG);

parts = {T_eval_clean};
if exist(CFG.RealisticCsv,'file') == 2
    try
        Treal = load_source(CFG.RealisticCsv, CFG.BaseModel, "realistic_readout");
        Treal = add_operational_risk_features(Treal, CFG);
        parts{end+1,1} = apply_risk_model(Treal, MODEL, "realistic_readout", CFG); %#ok<AGROW>
    catch ME
        warning('Test58:RealisticLoadFailed', 'Could not evaluate realistic source: %s', ME.message);
    end
else
    warning('Test58:MissingRealistic', 'Realistic Test54 CSV not found; clean-only evaluation will be written.');
end
Tall = vertcat(parts{:});

T_overall = summarize_predictions(Tall, ["source_domain"]);
T_by_source_region = summarize_predictions(Tall, ["source_domain","roi_region"]);
T_by_case = summarize_predictions(Tall, ["source_domain","case_id"]);
T_by_freq = summarize_predictions(Tall, ["source_domain","f0"]);
T_by_regime = summarize_predictions(Tall, ["source_domain","field_regime_ood"]);
T_by_risk_bin = summarize_predictions(Tall, ["source_domain","risk_bin"]);
T_thresh = threshold_summary(Tall, CFG);
T_auc = auc_summary(Tall);
T_cal = calibration_summary(Tall, CFG);

writetable(remove_cell_columns(Tall), fullfile(OUT.table_dir,'test58_patch_level_results.csv'));
writetable(T_overall, fullfile(OUT.table_dir,'test58_estimability_summary_overall.csv'));
writetable(T_by_source_region, fullfile(OUT.table_dir,'test58_estimability_summary_by_region.csv'));
writetable(T_by_case, fullfile(OUT.table_dir,'test58_estimability_summary_by_case.csv'));
writetable(T_by_freq, fullfile(OUT.table_dir,'test58_estimability_summary_by_frequency.csv'));
writetable(T_by_regime, fullfile(OUT.table_dir,'test58_estimability_summary_by_regime.csv'));
writetable(T_by_risk_bin, fullfile(OUT.table_dir,'test58_estimability_summary_by_risk_bin.csv'));
writetable(T_thresh, fullfile(OUT.table_dir,'test58_threshold_summary.csv'));
writetable(T_auc, fullfile(OUT.table_dir,'test58_auc_summary.csv'));
writetable(T_cal, fullfile(OUT.table_dir,'test58_reliability_calibration.csv'));

safe_plot(@() plot_test58_summary(T_thresh, T_by_risk_bin, T_cal, T_auc, OUT), 'summary plots');
save(fullfile(OUT.data_dir,'test58_estimability_risk_model.mat'), ...
    'CFG','MODEL','T_overall','T_by_source_region','T_by_case','T_by_freq', ...
    'T_by_regime','T_by_risk_bin','T_thresh','T_auc','T_cal','-v7.3');

print_summary(T_auc, T_thresh, T_by_risk_bin, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 58 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST58_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["quick","full"]), 'ADAPTIVE_REQ_TEST58_MODE must be quick or full.');
CFG = struct();
CFG.Mode = mode;
CFG.QuickMode = mode == "quick";
CFG.ValidateOnly = env_true('ADAPTIVE_REQ_TEST58_VALIDATE_ONLY', false);
CFG.RandomSeed = 58001;
CFG.CleanCsv = env_string('ADAPTIVE_REQ_TEST58_CLEAN_CSV', ...
    fullfile(root_dir,'outputs','test_55_controlled_field_frequency_roi_sweep', ...
    'tables','test55_patch_level_results.csv'));
CFG.RealisticCsv = env_string('ADAPTIVE_REQ_TEST58_REALISTIC_CSV', ...
    fullfile(root_dir,'outputs','test_54_realistic_readout_snr_validation', ...
    'tables','test54_patch_level_results.csv'));
CFG.BaseModel = env_string('ADAPTIVE_REQ_TEST58_BASE_MODEL', ...
    "q_spectrum_plus_composition");
CFG.TrainFraction = env_number('ADAPTIVE_REQ_TEST58_TRAIN_FRACTION', 0.70);
CFG.CsGuess = env_number('ADAPTIVE_REQ_TEST58_CS_GUESS', 3.0);
CFG.TreeLearners = env_number('ADAPTIVE_REQ_TEST58_TREE_LEARNERS', 180);
CFG.MinLeafSize = env_number('ADAPTIVE_REQ_TEST58_MIN_LEAF_SIZE', 10);
CFG.MaxTrainRows = env_number('ADAPTIVE_REQ_TEST58_MAX_TRAIN_ROWS', ...
    ternary(mode=="quick", 80000, 250000));
CFG.ValidateRows = env_number('ADAPTIVE_REQ_TEST58_VALIDATE_ROWS', 25000);
CFG.UseParallel = env_true('ADAPTIVE_REQ_TEST58_USE_PARALLEL', false);
CFG.RiskThresholds = [0.2 0.5 0.8];
CFG.CalibrationEdges = 0:0.1:1;
CFG.FeatureVars = ["q_pred";"sws_pred";"log_sws_pred";"M_eff_pred"; ...
    "k_pred_firstpass";"predicted_patch_purity";"p_mixed";"p_strong_mixed"; ...
    "M";"f0";"dx";"dz";"field_code";"case_family_code"];
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.QuickMode
    OUT.root_dir = fullfile(root_dir,'outputs','test_58_estimability_risk_mask','quick');
else
    OUT.root_dir = fullfile(root_dir,'outputs','test_58_estimability_risk_mask');
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

%% Data and model

function T = load_source(file, base_model, source_domain)
T = readtable(file);
if ismember("model_name", string(T.Properties.VariableNames))
    T = T(string(T.model_name) == string(base_model), :);
elseif ismember("strategy_name", string(T.Properties.VariableNames))
    T = T(string(T.strategy_name) == string(base_model), :);
end
assert(~isempty(T), 'No rows found in %s for %s.', file, base_model);
T.source_domain = repmat(string(source_domain), height(T), 1);
if ~ismember("roi_region", string(T.Properties.VariableNames))
    T.roi_region = repmat("unknown", height(T), 1);
end
if ~ismember("field_regime_ood", string(T.Properties.VariableNames)) && ...
        ismember("field_regime", string(T.Properties.VariableNames))
    T.field_regime_ood = string(T.field_regime);
end
if ~ismember("case_family", string(T.Properties.VariableNames))
    T.case_family = string(T.case_id);
end
end

function T = add_operational_risk_features(T, CFG)
T.sws_pred = clamp(T.sws_pred, 0.2, 20);
T.q_pred = clamp(T.q_pred, 0, 1);
T.k_pred_firstpass = 2*pi*T.f0 ./ T.sws_pred;
T.log_sws_pred = log(T.sws_pred);
T.lambda_pred_m = T.sws_pred ./ T.f0;
T.window_length_m = T.M .* CFG.CsGuess ./ T.f0;
T.M_eff_pred = T.window_length_m ./ T.lambda_pred_m;
T.field_code = categorical_code(string(T.field_regime_ood));
T.case_family_code = categorical_code(string(T.case_family));
if ~ismember("high_error10", string(T.Properties.VariableNames))
    T.high_error10 = abs(100*(T.sws_pred - T.true_SWS)./T.true_SWS) > 10;
end
if ~ismember("high_error20", string(T.Properties.VariableNames))
    T.high_error20 = abs(100*(T.sws_pred - T.true_SWS)./T.true_SWS) > 20;
end
end

function MODEL = train_risk_models(T, train_mask, CFG)
Xok = all(isfinite(table2array(T(:, cellstr(CFG.FeatureVars)))),2);
idx = find(train_mask & Xok & isfinite(T.sws_abs_error_pct));
if CFG.MaxTrainRows > 0 && numel(idx) > CFG.MaxTrainRows
    rng(CFG.RandomSeed);
    idx = idx(randperm(numel(idx), CFG.MaxTrainRows));
end
template = templateTree('MinLeafSize', CFG.MinLeafSize);
opts = statset('UseParallel', CFG.UseParallel);
MODEL = struct();
MODEL.features = CFG.FeatureVars(:);
MODEL.high10 = fitcensemble(T(idx, cellstr(CFG.FeatureVars)), logical(T.high_error10(idx)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template,'Options',opts);
MODEL.high20 = fitcensemble(T(idx, cellstr(CFG.FeatureVars)), logical(T.high_error20(idx)), ...
    'Method','Bag','NumLearningCycles',CFG.TreeLearners,'Learners',template,'Options',opts);
MODEL.n_train = numel(idx);
fprintf('Trained estimability risk models on %d rows.\n', MODEL.n_train);
end

function T = apply_risk_model(T, MODEL, source_domain, CFG)
X = T(:, cellstr(MODEL.features));
[~,s10] = predict(MODEL.high10, X);
[~,s20] = predict(MODEL.high20, X);
T.source_domain = repmat(string(source_domain), height(T), 1);
T.risk_high_error10 = positive_score(MODEL.high10, s10);
T.risk_high_error20 = positive_score(MODEL.high20, s20);
T.reliability = 1 - T.risk_high_error20;
T.risk_bin = discretize_risk(T.risk_high_error20, CFG.CalibrationEdges);
T.estimable_risk_lt_0p2 = T.risk_high_error20 < 0.2;
T.estimable_risk_lt_0p5 = T.risk_high_error20 < 0.5;
T.estimable_risk_lt_0p8 = T.risk_high_error20 < 0.8;
end

function s = positive_score(model, score)
classes = string(model.ClassNames);
idx = find(classes == "true" | classes == "1", 1);
if isempty(idx), idx = size(score,2); end
s = score(:,idx);
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
        mean(X.sws_signed_error_pct,'omitnan'), 100*mean(X.high_error10,'omitnan'), ...
        100*mean(X.high_error20,'omitnan'), mean(X.risk_high_error20,'omitnan'), ...
        mean(X.reliability,'omitnan'), ...
        'VariableNames', {'N','MAPE_pct','mean_signed_error_pct', ...
        'high_error10_pct','high_error20_pct','mean_risk_high20','mean_reliability'});
end
S = [groups vertcat(rows{:})];
end

function S = threshold_summary(T, CFG)
parts = {};
for src = unique(T.source_domain,'stable')'
    Xsrc = T(T.source_domain==src,:);
    for th = CFG.RiskThresholds
        low = Xsrc.risk_high_error20 >= th;
        high = ~low;
        parts{end+1,1} = table(src, th, height(Xsrc), 100*mean(low,'omitnan'), ...
            mean(Xsrc.sws_abs_error_pct(high),'omitnan'), ...
            mean(Xsrc.sws_abs_error_pct(low),'omitnan'), ...
            100*mean(Xsrc.high_error20(high),'omitnan'), ...
            100*mean(Xsrc.high_error20(low),'omitnan'), ...
            100*mean(low & Xsrc.high_error20,'omitnan') / max(mean(Xsrc.high_error20,'omitnan'), eps), ...
            100*mean(Xsrc.high_error20(low),'omitnan'), ...
            'VariableNames', {'source_domain','risk_threshold','N','low_confidence_pct', ...
            'MAPE_high_confidence','MAPE_low_confidence', ...
            'high20_rate_high_confidence','high20_rate_low_confidence', ...
            'recall_high20_pct','precision_high20_pct'}); %#ok<AGROW>
    end
end
S = vertcat(parts{:});
end

function S = auc_summary(T)
parts = {};
for src = unique(T.source_domain,'stable')'
    X = T(T.source_domain==src,:);
    parts{end+1,1} = table(src, height(X), auc_binary(X.high_error20, X.risk_high_error20, "roc"), ...
        auc_binary(X.high_error20, X.risk_high_error20, "pr"), ...
        auc_binary(X.high_error10, X.risk_high_error10, "roc"), ...
        auc_binary(X.high_error10, X.risk_high_error10, "pr"), ...
        'VariableNames', {'source_domain','N','ROC_AUC_high20','PR_AUC_high20', ...
        'ROC_AUC_high10','PR_AUC_high10'}); %#ok<AGROW>
end
S = vertcat(parts{:});
end

function S = calibration_summary(T, CFG)
parts = {};
for src = unique(T.source_domain,'stable')'
    X = T(T.source_domain==src,:);
    for i = 1:numel(CFG.CalibrationEdges)-1
        lo = CFG.CalibrationEdges(i); hi = CFG.CalibrationEdges(i+1);
        idx = X.risk_high_error20 >= lo & X.risk_high_error20 < hi;
        if i == numel(CFG.CalibrationEdges)-1
            idx = X.risk_high_error20 >= lo & X.risk_high_error20 <= hi;
        end
        parts{end+1,1} = table(src, lo, hi, sum(idx), mean(X.risk_high_error20(idx),'omitnan'), ...
            100*mean(X.high_error20(idx),'omitnan'), mean(X.sws_abs_error_pct(idx),'omitnan'), ...
            'VariableNames', {'source_domain','risk_bin_low','risk_bin_high','N', ...
            'mean_predicted_risk','observed_high20_rate_pct','MAPE_pct'}); %#ok<AGROW>
    end
end
S = vertcat(parts{:});
end

function plot_test58_summary(T_thresh, T_risk_bin, T_cal, T_auc, OUT)
fig = figure('Color','w','Units','centimeters','Position',[1 1 34 18]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
bar(ax, categorical(string(T_auc.source_domain)), T_auc.PR_AUC_high20);
ylabel(ax,'PR AUC'); title(ax,'High-error risk PR AUC','FontWeight','normal'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for src = unique(T_thresh.source_domain,'stable')'
    X = sortrows(T_thresh(T_thresh.source_domain==src,:), 'risk_threshold');
    plot(ax, X.risk_threshold, X.recall_high20_pct, '-o', 'DisplayName', src);
end
xlabel(ax,'Risk threshold'); ylabel('Recall high-error >20% (%)');
title(ax,'Recall of high-error pixels','FontWeight','normal'); legend(ax,'Location','best'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for src = unique(T_risk_bin.source_domain,'stable')'
    X = T_risk_bin(T_risk_bin.source_domain==src,:);
    plot(ax, categorical(string(X.risk_bin)), X.MAPE_pct, '-o', 'DisplayName', src);
end
xlabel(ax,'Predicted risk bin'); ylabel(ax,'MAPE (%)');
title(ax,'Error vs predicted risk','FontWeight','normal'); legend(ax,'Location','best'); grid(ax,'on');
ax = nexttile(tl); hold(ax,'on');
for src = unique(T_cal.source_domain,'stable')'
    X = T_cal(T_cal.source_domain==src & T_cal.N>0,:);
    plot(ax, X.mean_predicted_risk, X.observed_high20_rate_pct/100, '-o', 'DisplayName', src);
end
plot(ax,[0 1],[0 1],'k--','DisplayName','ideal');
xlabel(ax,'Predicted high-error probability'); ylabel(ax,'Observed high-error fraction');
title(ax,'Reliability calibration','FontWeight','normal'); legend(ax,'Location','best'); grid(ax,'on');
export_fig(fig, fullfile(OUT.figure_dir,'test58_estimability_summary.png'));
end

function print_summary(T_auc, T_thresh, T_risk_bin, OUT)
fprintf('\n================ Test 58 summary ================\n');
disp(T_auc);
X = T_thresh(T_thresh.risk_threshold==0.5,:);
if ~isempty(X)
    disp(X(:, {'source_domain','risk_threshold','low_confidence_pct', ...
        'MAPE_high_confidence','MAPE_low_confidence','recall_high20_pct'}));
end
fprintf('Outputs: %s\n', OUT.root_dir);
end

%% Helpers

function [train_mask, test_mask] = condition_split(T, CFG)
keys = unique(string(T.condition_key), 'stable');
rng(CFG.RandomSeed);
keys = keys(randperm(numel(keys)));
ntrain = max(1, round(CFG.TrainFraction*numel(keys)));
train_keys = keys(1:ntrain);
train_mask = ismember(string(T.condition_key), train_keys);
test_mask = ~train_mask;
end

function y = auc_binary(ytrue, score, mode)
ytrue = logical(ytrue);
score = double(score);
ok = isfinite(score);
ytrue = ytrue(ok); score = score(ok);
if numel(unique(ytrue)) < 2, y = NaN; return; end
[score, ord] = sort(score, 'descend');
ytrue = ytrue(ord);
tp = cumsum(ytrue);
fp = cumsum(~ytrue);
P = sum(ytrue); N = sum(~ytrue);
rec = tp / max(P, eps);
fpr = fp / max(N, eps);
prec = tp ./ max(tp+fp, eps);
switch string(mode)
    case "roc"
        y = trapz([0; fpr; 1], [0; rec; 1]);
    otherwise
        y = trapz([0; rec], [1; prec]);
end
end

function labels = discretize_risk(r, edges)
labels = strings(size(r));
for i = 1:numel(edges)-1
    idx = r >= edges(i) & r < edges(i+1);
    if i == numel(edges)-1, idx = r >= edges(i) & r <= edges(i+1); end
    labels(idx) = sprintf('%.1f-%.1f', edges(i), edges(i+1));
end
labels(labels=="") = "unknown";
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

function y = clamp(x, lo, hi), y = min(max(x,lo),hi); end
function y = ternary(tf, a, b), if tf, y = a; else, y = b; end, end

function T = remove_cell_columns(T)
vars = string(T.Properties.VariableNames);
drop = false(size(vars));
for i = 1:numel(vars), drop(i) = iscell(T.(vars(i))); end
T(:, cellstr(vars(drop))) = [];
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
    warning('Test58:PlotFailed', 'Plot failed (%s): %s', label, ME.message);
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
