%% analyze_test_11_level04_grouped_generalization.m
% Test 11 Level 04: clean grouped generalization tests.
%
% This analysis does not regenerate Test 11. It loads the saved
% test_11_global_req_features table and evaluates leave-one-group-out splits:
%   1. leave-one-frequency-out  (SIM_f0)
%   2. leave-one-M-out          (REQ_M)
%   3. leave-one-wave-model-out (SIM_WaveModel)
%   4. leave-one-aperture-out   (step_idx, when available)

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

addpath(root_dir);
root_dir = setup_adaptive_req();

%% Load Test 11

[T_feat, MC, PATHS] = adaptive_req.analysis.load_mc_results( ...
    'test_11_global_req_features', ...
    'RootDir', root_dir, ...
    'Verbose', true);

analysis_dir = fullfile(PATHS.analysis_dir, ...
    'level_04_grouped_generalization');
fig_dir = fullfile(analysis_dir, 'figures');
table_dir = fullfile(analysis_dir, 'tables');
model_dir = fullfile(analysis_dir, 'models');

make_dir_if_needed(analysis_dir);
make_dir_if_needed(fig_dir);
make_dir_if_needed(table_dir);
make_dir_if_needed(model_dir);

%% Required variables and deployment features

required_vars = [
    "q_theory"
    "req_mapping"
    "SIM_f0"
    "REQ_M"
    "SIM_WaveModel"
    "condition_id"
    "step_idx"
    "realization_idx"
    "patch_idx"];

assert(all(ismember(required_vars, string(T_feat.Properties.VariableNames))), ...
    'Test 11 data is missing required variables: %s', ...
    strjoin(required_vars(~ismember(required_vars, ...
    string(T_feat.Properties.VariableNames))), ', '));

T_feat = add_local_ecum_features_from_mapping(T_feat);
T_feat.row_key = make_row_key(T_feat);

%% Predictor sets, matching Level 01 but filtered for leakage

base_local = [
    "radial_entropy"
    "radial_peak_width"
    "radial_k_peak_norm"
    "radial_centroid_norm"
    "radial_std_norm"
    "ang_entropy"
    "ang_resultant_R"
    "ang_resultant_R2"
    "ang_moment_1"
    "ang_moment_2"
    "ang_moment_4"
    "ang_peak_count_rel"
    "ang_top2_to_top1"
    "ang_peak_separation_deg"
    "ecum_width_50_rel"
    "ecum_width_80_rel"
    "ecum_asymmetry_25_75"
    "ecum_width_ratio_80_50"
    "ecum_increment_entropy"
    "ecum_increment_peak_frac"
    "ecum_increment_gini"
    "ecum_slope_max"
    "ecum_slope_peak_to_mean"
    "ecum_slope_iqr_to_median"
    "srad_proxy_centroid_k_norm"
    "srad_proxy_std_k_norm"
    "srad_proxy_skewness"
    "srad_proxy_kurtosis"
    "srad_proxy_peak_k_norm"
    "srad_proxy_peak_to_centroid"
    "srad_proxy_low_side_frac"
    "srad_proxy_high_side_frac"
    "REQ_M"
    "SIM_f0"
    "REQ_Nbins_effective"];

base_global = [
    "global_radial_entropy"
    "global_radial_peak_width"
    "global_radial_k_peak_norm"
    "global_radial_centroid_norm"
    "global_radial_std_norm"
    "global_ang_entropy"
    "global_ang_resultant_R"
    "global_ang_resultant_R2"
    "global_ang_moment_1"
    "global_ang_moment_2"
    "global_ang_moment_4"
    "global_ang_peak_count_rel"
    "global_ang_top2_to_top1"
    "global_ang_peak_separation_deg"
    "global_ecum_width_50_rel"
    "global_ecum_width_80_rel"
    "global_ecum_asymmetry_25_75"
    "global_ecum_width_ratio_80_50"
    "global_ecum_increment_entropy"
    "global_ecum_increment_peak_frac"
    "global_ecum_increment_gini"
    "global_ecum_slope_max"
    "global_ecum_slope_peak_to_mean"
    "global_ecum_slope_iqr_to_median"
    "global_srad_proxy_centroid_k_norm"
    "global_srad_proxy_std_k_norm"
    "global_srad_proxy_skewness"
    "global_srad_proxy_kurtosis"
    "global_srad_proxy_peak_k_norm"
    "global_srad_proxy_peak_to_centroid"
    "global_srad_proxy_low_side_frac"
    "global_srad_proxy_high_side_frac"
    "REQ_M"
    "SIM_f0"
    "global_REQ_Nbins_effective"];

local_predictors = filter_operational_predictors( ...
    existing_numeric_predictors(T_feat, base_local));
global_predictors = filter_operational_predictors( ...
    existing_numeric_predictors(T_feat, base_global));
hybrid_predictors = filter_operational_predictors( ...
    unique([local_predictors; global_predictors], 'stable'));

model_specs = struct([]);
model_specs(1).name = "LocalOnly";
model_specs(1).predictors = local_predictors;
model_specs(1).role = "operational";
model_specs(2).name = "GlobalOnly";
model_specs(2).predictors = global_predictors;
model_specs(2).role = "operational";
model_specs(3).name = "HybridLocalGlobal";
model_specs(3).predictors = hybrid_predictors;
model_specs(3).role = "operational";

for i = 1:numel(model_specs)
    assert_no_leakage_predictors(model_specs(i).predictors, ...
        model_specs(i).name);
end

%% Grouped generalization experiments

split_specs = struct([]);
split_specs(1).name = "leave_one_frequency";
split_specs(1).var = "SIM_f0";
split_specs(2).name = "leave_one_M";
split_specs(2).var = "REQ_M";
split_specs(3).name = "leave_one_wave_model";
split_specs(3).var = "SIM_WaveModel";
split_specs(4).name = "leave_one_aperture";
split_specs(4).var = "step_idx";

model_types = ["linear", "boosted_trees", "bagged_trees"];

T_all_pred = table();
T_q_metrics = table();
MODELS = struct([]);
model_counter = 0;

for sidx = 1:numel(split_specs)
    split_name = split_specs(sidx).name;
    heldout_var = split_specs(sidx).var;

    if ~ismember(heldout_var, string(T_feat.Properties.VariableNames))
        fprintf('Skipping %s because %s is missing.\n', ...
            split_name, heldout_var);
        continue;
    end

    heldout_values = unique(T_feat.(char(heldout_var)), 'stable');
    if numel(heldout_values) < 2
        fprintf('Skipping %s because %s has < 2 groups.\n', ...
            split_name, heldout_var);
        continue;
    end

    for hidx = 1:numel(heldout_values)
        heldout_value = heldout_values(hidx);
        test_mask = is_group_value(T_feat.(char(heldout_var)), heldout_value);
        train_mask = ~test_mask;

        assert_grouped_split(T_feat, train_mask, test_mask, ...
            heldout_var, heldout_value);

        for midx = 1:numel(model_specs)
            spec = model_specs(midx);
            fprintf('\n=== %s | %s = %s | %s ===\n', ...
                split_name, heldout_var, value_to_string(heldout_value), ...
                spec.name);

            [MODEL_i, T_pred_i, T_metrics_i] = ...
                adaptive_req.analysis.train_q_model_fixed_split( ...
                    T_feat, spec.predictors, train_mask, test_mask, ...
                    'ModelName', spec.name, ...
                    'ModelRole', spec.role, ...
                    'ModelTypes', model_types, ...
                    'NumLearningCycles', 200, ...
                    'MinLeafSize', 8, ...
                    'Verbose', true);

            T_pred_i = add_split_metadata(T_pred_i, split_name, ...
                heldout_var, heldout_value, train_mask, test_mask, spec.role);
            T_metrics_i = add_split_metadata(T_metrics_i, split_name, ...
                heldout_var, heldout_value, train_mask, test_mask, spec.role);

            model_counter = model_counter + 1;
            MODELS(model_counter).generalization_test = split_name;
            MODELS(model_counter).heldout_var = heldout_var;
            MODELS(model_counter).heldout_value = value_to_string(heldout_value);
            MODELS(model_counter).name = spec.name;
            MODELS(model_counter).role = spec.role;
            MODELS(model_counter).predictors = spec.predictors;
            MODELS(model_counter).model = MODEL_i;

            T_all_pred = concat_tables_with_missing(T_all_pred, T_pred_i);
            T_q_metrics = concat_tables_with_missing(T_q_metrics, T_metrics_i);
        end
    end
end

assert(~isempty(T_all_pred), 'No grouped predictions were generated.');
assert(~isempty(T_q_metrics), 'No q metrics were generated.');
assert(all(T_all_pred.model_role == "operational"), ...
    'Operational Level 04 predictions must be marked operational.');

%% Convert q predictions to local SWS and summarize

T_sws = add_sws_metrics(T_all_pred, T_feat);
T_sws_test = T_sws(T_sws.split == "test", :);

assert(~isempty(T_sws_test), 'No test SWS predictions were generated.');
assert(mean(isfinite(T_sws_test.q_pred)) > 0.95, ...
    'Too many non-finite q_pred values.');
assert(mean(isfinite(T_sws_test.cs_pred)) > 0.95, ...
    'Too many non-finite cs_pred values.');

T_sws_metrics = summarize_sws_metrics(T_sws_test, [
    "generalization_test"
    "heldout_var"
    "heldout_value"
    "model_name"
    "model_type"]);

T_sws_by_model = summarize_sws_metrics(T_sws_test, [
    "generalization_test"
    "model_name"
    "model_type"]);

T_sws_by_heldout = summarize_sws_metrics(T_sws_test, [
    "heldout_var"
    "heldout_value"
    "model_name"
    "model_type"]);

T_high_error = summarize_high_error_rates(T_sws_test, [
    "generalization_test"
    "heldout_var"
    "heldout_value"
    "model_name"
    "model_type"], 20);

assert(~isempty(T_sws_metrics), 'T_sws_metrics is empty.');
assert(~isempty(T_sws_by_model), 'T_sws_by_model is empty.');
assert(~isempty(T_sws_by_heldout), 'T_sws_by_heldout is empty.');
assert(~isempty(T_high_error), 'T_high_error is empty.');

%% Save tables and models

writetable(remove_cell_columns(T_sws), fullfile(table_dir, ...
    'level11_level04_grouped_predictions.csv'));
writetable(T_q_metrics, fullfile(table_dir, ...
    'level11_level04_q_metrics.csv'));
writetable(T_sws_metrics, fullfile(table_dir, ...
    'level11_level04_sws_metrics.csv'));
writetable(T_sws_by_heldout, fullfile(table_dir, ...
    'level11_level04_sws_metrics_by_heldout.csv'));
writetable(T_sws_by_model, fullfile(table_dir, ...
    'level11_level04_sws_metrics_by_model.csv'));
writetable(T_high_error, fullfile(table_dir, ...
    'level11_level04_high_error_gt20.csv'));

save(fullfile(model_dir, ...
    'level11_level04_grouped_generalization_models.mat'), ...
    'MODELS', 'T_q_metrics', 'T_sws_metrics', ...
    'T_sws_by_heldout', 'T_sws_by_model', 'T_high_error', ...
    'MC', 'PATHS', '-v7.3');

%% Figures

plot_mape_by_heldout(T_sws_metrics, "leave_one_frequency", ...
    "SIM_f0", fig_dir, 'level11_level04_mape_leave_one_frequency.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_M", ...
    "REQ_M", fig_dir, 'level11_level04_mape_leave_one_M.png');
plot_mape_by_heldout(T_sws_metrics, "leave_one_wave_model", ...
    "SIM_WaveModel", fig_dir, ...
    'level11_level04_mape_leave_one_wave_model.png');
plot_q_true_vs_pred_by_test(T_sws_test, fig_dir);
plot_sws_error_box_by_test(T_sws_test, fig_dir);
plot_mape_summary_by_generalization_test(T_sws_by_model, fig_dir);

%% Final console report

fprintf('\nTest 11 Level 04 grouped generalization complete.\n');
fprintf('Analysis folder:\n%s\n', analysis_dir);

T_console = T_sws_metrics(T_sws_metrics.model_type == "bagged_trees", ...
    {'generalization_test', 'heldout_value', 'model_name', 'MAPE_pct'});
T_console = sortrows(T_console, ...
    {'generalization_test', 'heldout_value', 'MAPE_pct'});
fprintf('\nBagged-trees MAPE by grouped holdout.\n');
disp(T_console);

T_best = best_model_by_generalization(T_sws_by_model);
fprintf('\nBest operational model by generalization test.\n');
disp(T_best);

%% Local functions

function T = add_local_ecum_features_from_mapping(T)

if ~ismember('req_mapping', T.Properties.VariableNames)
    return;
end

feat0 = adaptive_req.quantile.extract_ecum_shape_features( ...
    T.req_mapping{find_first_mapping(T.req_mapping)});
names = fieldnames(feat0);

for j = 1:numel(names)
    if ~ismember(names{j}, T.Properties.VariableNames)
        T.(names{j}) = nan(height(T), 1);
    end
end

for i = 1:height(T)
    if isempty(T.req_mapping{i})
        continue;
    end
    feat_i = adaptive_req.quantile.extract_ecum_shape_features( ...
        T.req_mapping{i});
    for j = 1:numel(names)
        T.(names{j})(i) = feat_i.(names{j});
    end
end

end

function idx = find_first_mapping(C)

idx = find(~cellfun(@isempty, C), 1, 'first');
if isempty(idx)
    error('No non-empty REQ mapping found.');
end

end

function predictors = existing_numeric_predictors(T, candidates)

vars = string(T.Properties.VariableNames);
candidates = string(candidates(:));
predictors = strings(0, 1);

for i = 1:numel(candidates)
    name_i = candidates(i);
    if ismember(name_i, vars) && isnumeric(T.(char(name_i)))
        predictors(end + 1, 1) = name_i; %#ok<AGROW>
    end
end

end

function predictors = filter_operational_predictors(predictors)

predictors = string(predictors(:));
lower_names = lower(predictors);

exact_banned = lower([
    "q_theory"
    "q_global_theory"
    "q_local_minus_global"
    "cs_true"
    "cs_pred"
    "sws_error"
    "abs_sws_error"
    "sws_error_pct"
    "abs_sws_error_pct"
    "q_true"
    "q_pred"
    "q_pred_raw"
    "residual"
    "abs_error"]);

bad = ismember(lower_names, exact_banned) | ...
    contains(lower_names, "error") | ...
    contains(lower_names, "residual") | ...
    contains(lower_names, "target");

predictors = predictors(~bad);

end

function assert_no_leakage_predictors(predictors, model_name)

predictors = string(predictors(:));
filtered = filter_operational_predictors(predictors);
assert(numel(filtered) == numel(predictors), ...
    'Operational model %s contains leakage predictors.', model_name);
assert(~isempty(predictors), ...
    'Operational model %s has no predictors after leakage filtering.', ...
    model_name);

end

function key = make_row_key(T)

parts = strings(height(T), 4);
parts(:, 1) = string(T.condition_id);
parts(:, 2) = string(T.step_idx);
parts(:, 3) = string(T.realization_idx);
parts(:, 4) = string(T.patch_idx);
key = join(parts, "|", 2);

end

function tf = is_group_value(x, value)

if isstring(x) || iscategorical(x) || iscellstr(x)
    tf = string(x) == string(value);
else
    tf = x == value;
end

tf = tf(:);

end

function assert_grouped_split(T, train_mask, test_mask, heldout_var, heldout_value)

assert(any(train_mask), 'Grouped split has empty train set.');
assert(any(test_mask), 'Grouped split has empty test set.');
assert(~any(train_mask & test_mask), 'Train/test masks overlap.');

group_values = T.(char(heldout_var));
heldout_in_train = any(is_group_value(group_values(train_mask), heldout_value));
assert(~heldout_in_train, ...
    'Heldout value appears in train split for %s = %s.', ...
    heldout_var, value_to_string(heldout_value));

heldout_in_test_only = all(is_group_value(group_values(test_mask), ...
    heldout_value));
assert(heldout_in_test_only, ...
    'Test split is not a pure heldout group for %s = %s.', ...
    heldout_var, value_to_string(heldout_value));

assert(isscalar(unique(string(group_values(test_mask)))), ...
    'Test split contains multiple group values; this is not leave-one-group-out.');

end

function T = add_split_metadata(T, generalization_test, heldout_var, ...
    heldout_value, train_mask, test_mask, model_role)

T.generalization_test = repmat(string(generalization_test), height(T), 1);
T.heldout_var = repmat(string(heldout_var), height(T), 1);
T.heldout_value = repmat(value_to_string(heldout_value), height(T), 1);
T.N_train = sum(train_mask) * ones(height(T), 1);
T.N_test = sum(test_mask) * ones(height(T), 1);
T.model_role = repmat(string(model_role), height(T), 1);

end

function s = value_to_string(x)

s = string(x);

end

function T_sws = add_sws_metrics(T_pred, T_ref)

[tf, loc] = ismember(string(T_pred.row_key), string(T_ref.row_key));
if ~all(tf)
    error('Could not match all prediction rows back to Test 11 reference rows.');
end

T_sws = T_pred;
T_sws.cs_true = T_ref.SIM_cs_bg(loc);
T_sws.cs_pred = q_to_cs(T_sws.q_pred, T_ref.req_mapping(loc), ...
    T_ref.SIM_f0(loc));
T_sws.sws_error = T_sws.cs_pred - T_sws.cs_true;
T_sws.abs_sws_error = abs(T_sws.sws_error);
T_sws.sws_error_pct = 100 * T_sws.sws_error ./ T_sws.cs_true;
T_sws.abs_sws_error_pct = abs(T_sws.sws_error_pct);

end

function cs = q_to_cs(q, mappings, f0)

q = double(q(:));
f0 = double(f0(:));
cs = nan(numel(q), 1);

for i = 1:numel(q)
    mapping_i = mappings{i};
    if isempty(mapping_i) || ~isfinite(q(i))
        continue;
    end

    k = adaptive_req.quantile.quantile_to_k(mapping_i, q(i));
    cs(i) = 2*pi*f0(i) ./ k;
end

end

function T_sum = summarize_sws_metrics(T, group_vars)

group_vars = string(group_vars);
[G, T_keys] = findgroups(T(:, cellstr(group_vars)));

T_sum = T_keys;
T_sum.N = splitapply(@numel, T.abs_sws_error_pct, G);
T_sum.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), ...
    T.sws_error_pct, G);
T_sum.MedAE_pct = splitapply(@(x) median(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.bias_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.sws_error_pct, G);
T_sum.P95_abs_error_pct = splitapply(@(x) prctile(x, 95), ...
    T.abs_sws_error_pct, G);
T_sum.Max_abs_error_pct = splitapply(@(x) max(x, [], 'omitnan'), ...
    T.abs_sws_error_pct, G);
T_sum.HighError_gt20_pct = splitapply(@(x) 100 * mean(x > 20, 'omitnan'), ...
    T.abs_sws_error_pct, G);

T_sum = sortrows(T_sum, 'MAPE_pct', 'ascend');

end

function T_sum = summarize_high_error_rates(T, group_vars, threshold_pct)

group_vars = string(group_vars);
[G, T_keys] = findgroups(T(:, cellstr(group_vars)));
is_high = T.abs_sws_error_pct > threshold_pct;

T_sum = T_keys;
T_sum.N = splitapply(@numel, is_high, G);
T_sum.N_high_error = splitapply(@sum, is_high, G);
T_sum.HighError_gt20_pct = 100 * T_sum.N_high_error ./ T_sum.N;
T_sum.MAPE_all_pct = splitapply(@(x) mean(x, 'omitnan'), ...
    T.abs_sws_error_pct, G);

T_sum = sortrows(T_sum, 'HighError_gt20_pct', 'descend');

end

function plot_mape_by_heldout(T_metrics, generalization_test, heldout_var, ...
    fig_dir, file_name)

T = T_metrics(T_metrics.generalization_test == generalization_test & ...
    T_metrics.model_type == "bagged_trees", :);
if isempty(T)
    return;
end

plot_grouped_mape(T, "heldout_value", ...
    sprintf('MAPE SWS: %s', generalization_test), heldout_var);
exportgraphics(gcf, fullfile(fig_dir, file_name), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end

function plot_grouped_mape(T, x_var, title_text, x_label)

models = unique(T.model_name, 'stable');
x_values = unique(T.(char(x_var)), 'stable');
Y = nan(numel(x_values), numel(models));

for i = 1:numel(x_values)
    for j = 1:numel(models)
        idx = T.(char(x_var)) == x_values(i) & T.model_name == models(j);
        if any(idx)
            Y(i, j) = mean(T.MAPE_pct(idx), 'omitnan');
        end
    end
end

figure('Color', 'w', 'Position', [100 100 920 520]);
bar(categorical(x_values), Y);
ylabel('MAPE SWS (%)');
xlabel(strrep(string(x_label), '_', '\_'));
title(title_text, 'Interpreter', 'none');
legend(models, 'Location', 'best', 'Interpreter', 'none');
grid on;

end

function plot_q_true_vs_pred_by_test(T, fig_dir)

T = T(T.model_type == "bagged_trees", :);
tests = unique(T.generalization_test, 'stable');
models = unique(T.model_name, 'stable');

figure('Color', 'w', 'Position', [100 100 1200 760]);
tl = tiledlayout(numel(tests), numel(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(tests)
    for j = 1:numel(models)
        ax = nexttile(tl);
        Ti = T(T.generalization_test == tests(i) & ...
            T.model_name == models(j), :);
        scatter(ax, Ti.q_true, Ti.q_pred, 8, 'filled', ...
            'MarkerFaceAlpha', 0.35);
        hold(ax, 'on');
        plot(ax, [0 1], [0 1], 'k--');
        axis(ax, 'equal');
        xlim(ax, [0 1]);
        ylim(ax, [0 1]);
        grid(ax, 'on');
        title(ax, sprintf('%s | %s', tests(i), models(j)), ...
            'Interpreter', 'none');
        xlabel(ax, 'q true');
        ylabel(ax, 'q predicted');
    end
end

exportgraphics(gcf, fullfile(fig_dir, ...
    'level11_level04_q_true_vs_pred_by_test.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end

function plot_sws_error_box_by_test(T, fig_dir)

T = T(T.model_type == "bagged_trees", :);
labels = T.generalization_test + " | " + T.model_name;

figure('Color', 'w', 'Position', [100 100 1250 540]);
boxchart(categorical(labels), T.abs_sws_error_pct);
yline(20, 'r--', '20%');
ylabel('|SWS error| (%)');
title('Grouped generalization SWS error, bagged trees');
xtickangle(25);
grid on;

exportgraphics(gcf, fullfile(fig_dir, ...
    'level11_level04_sws_error_box_by_test.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end

function plot_mape_summary_by_generalization_test(T_by_model, fig_dir)

T = T_by_model(T_by_model.model_type == "bagged_trees", :);
plot_grouped_mape(T, "generalization_test", ...
    'MAPE summary by grouped generalization test', ...
    'generalization_test');
exportgraphics(gcf, fullfile(fig_dir, ...
    'level11_level04_mape_summary_by_generalization_test.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
close(gcf);

end

function T_best = best_model_by_generalization(T_by_model)

T = T_by_model(T_by_model.model_type == "bagged_trees", :);
tests = unique(T.generalization_test, 'stable');
T_best = table();

for i = 1:numel(tests)
    Ti = T(T.generalization_test == tests(i), :);
    Ti = sortrows(Ti, 'MAPE_pct', 'ascend');
    if ~isempty(Ti)
        T_best = [T_best; Ti(1, {'generalization_test', ...
            'model_name', 'MAPE_pct', 'HighError_gt20_pct'})]; %#ok<AGROW>
    end
end

end

function T = remove_cell_columns(T)

vars = T.Properties.VariableNames;
remove = false(size(vars));
for i = 1:numel(vars)
    remove(i) = iscell(T.(vars{i}));
end
T(:, remove) = [];

end

function T = concat_tables_with_missing(A, B)

if isempty(A)
    T = B;
    return;
end
if isempty(B)
    T = A;
    return;
end

vars_all = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');

A = add_missing_columns(A, vars_all);
B = add_missing_columns(B, vars_all);

T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];

end

function T = add_missing_columns(T, vars_all)

vars = string(T.Properties.VariableNames);
for i = 1:numel(vars_all)
    name_i = char(vars_all(i));
    if ismember(vars_all(i), vars)
        continue;
    end

    string_like = any(endsWith(vars_all(i), ...
        ["name", "type", "role", "source", "label", "split", ...
        "test", "var", "value"])) || startsWith(vars_all(i), "SIM_WaveModel");

    if string_like
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
