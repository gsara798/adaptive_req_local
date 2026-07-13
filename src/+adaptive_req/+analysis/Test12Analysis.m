classdef Test12Analysis
    %TEST12ANALYSIS Shared utilities for Test 12 adaptive-REQ analysis.

    methods (Static)

        function paths = makeOutputDirs(PATHS, level_name)
            paths.analysis_dir = fullfile(PATHS.analysis_dir, level_name);
            paths.fig_dir = fullfile(paths.analysis_dir, 'figures');
            paths.table_dir = fullfile(paths.analysis_dir, 'tables');
            paths.model_dir = fullfile(paths.analysis_dir, 'models');
            dirs = string(struct2cell(paths));
            for i = 1:numel(dirs)
                if ~exist(dirs(i), 'dir')
                    mkdir(dirs(i));
                end
            end
        end

        function T = prepareFeatureTable(T)
            T = adaptive_req.analysis.Test12Analysis.addLocalEcumFeatures(T);
            if ~ismember('row_key', string(T.Properties.VariableNames))
                T.row_key = adaptive_req.analysis.Test12Analysis.makeRowKey(T);
            end
        end

        function requireVars(T, required_vars, context)
            required_vars = string(required_vars(:));
            vars = string(T.Properties.VariableNames);
            missing = required_vars(~ismember(required_vars, vars));
            assert(isempty(missing), '%s is missing required variables: %s', ...
                context, strjoin(missing, ', '));
        end

        function specs = buildModelSpecs(T)
            base_local = [
                "radial_entropy"
                "radial_peak_width"
                "radial_k_peak_norm"
                "radial_centroid_norm"
                "radial_std_norm"
                "width_75_25_rel"
                "width_90_50_rel"
                "width_90_10_rel"
                "lowk_frac_rel"
                "midband_frac_rel"
                "highk_frac_rel"
                "ang_entropy"
                "circ_var"
                "dom_dir_frac"
                "window_max_frac"
                "window_cf"
                "ang_moment_1"
                "ang_moment_2"
                "ang_moment_4"
                "ang_peak_count_rel"
                "ang_top1_window_frac"
                "ang_top2_window_frac"
                "ang_top3_window_frac"
                "ang_top2_to_top1"
                "ang_peak_separation_deg"
                "ecum_k10_norm_k50"
                "ecum_k25_norm_k50"
                "ecum_k75_norm_k50"
                "ecum_k90_norm_k50"
                "ecum_width_50_rel"
                "ecum_width_80_rel"
                "ecum_lower_tail_rel"
                "ecum_upper_tail_rel"
                "ecum_asymmetry_10_90"
                "ecum_asymmetry_25_75"
                "ecum_width_ratio_80_50"
                "ecum_lower_upper_width_ratio"
                "ecum_auc_norm"
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
                "global_width_75_25_rel"
                "global_width_90_50_rel"
                "global_width_90_10_rel"
                "global_lowk_frac_rel"
                "global_midband_frac_rel"
                "global_highk_frac_rel"
                "global_ang_entropy"
                "global_circ_var"
                "global_dom_dir_frac"
                "global_window_max_frac"
                "global_window_cf"
                "global_ang_moment_1"
                "global_ang_moment_2"
                "global_ang_moment_4"
                "global_ang_peak_count_rel"
                "global_ang_top1_window_frac"
                "global_ang_top2_window_frac"
                "global_ang_top3_window_frac"
                "global_ang_top2_to_top1"
                "global_ang_peak_separation_deg"
                "global_ecum_k10_norm_k50"
                "global_ecum_k25_norm_k50"
                "global_ecum_k75_norm_k50"
                "global_ecum_k90_norm_k50"
                "global_ecum_width_50_rel"
                "global_ecum_width_80_rel"
                "global_ecum_lower_tail_rel"
                "global_ecum_upper_tail_rel"
                "global_ecum_asymmetry_10_90"
                "global_ecum_asymmetry_25_75"
                "global_ecum_width_ratio_80_50"
                "global_ecum_lower_upper_width_ratio"
                "global_ecum_auc_norm"
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

            local = adaptive_req.analysis.Test12Analysis.existingPredictors(T, base_local);
            global_predictors = adaptive_req.analysis.Test12Analysis.existingPredictors(T, base_global);
            hybrid = unique([local; global_predictors], 'stable');

            families = struct([]);
            families(1).name = "LocalOnly";
            families(1).predictors = local;
            families(2).name = "GlobalOnly";
            families(2).predictors = global_predictors;
            families(3).name = "HybridLocalGlobal";
            families(3).predictors = hybrid;

            variants = [
                "NoCsGuess"
                "WithCsGuess"
                "WithMeffGuess"
                "WithCsGuessAndMeffGuess"
                "DiagnosticWithMeffTrue"];

            specs = struct([]);
            n = 0;
            for i = 1:numel(families)
                for j = 1:numel(variants)
                    add = strings(0, 1);
                    role = "operational";
                    if variants(j) == "WithCsGuess"
                        add = "REQ_cs_guess";
                    elseif variants(j) == "WithMeffGuess"
                        add = "M_eff_guess";
                    elseif variants(j) == "WithCsGuessAndMeffGuess"
                        add = ["REQ_cs_guess"; "M_eff_guess"];
                    elseif variants(j) == "DiagnosticWithMeffTrue"
                        add = ["REQ_cs_guess"; "M_eff_guess"; "M_eff_true_diag"];
                        role = "diagnostic_only";
                    end

                    predictors = adaptive_req.analysis.Test12Analysis.existingPredictors( ...
                        T, unique([families(i).predictors; add(:)], 'stable'));

                    if role == "operational"
                        predictors = adaptive_req.analysis.Test12Analysis.filterOperationalPredictors(predictors);
                        adaptive_req.analysis.Test12Analysis.assertNoLeakage(predictors, families(i).name);
                    end

                    n = n + 1;
                    specs(n).model_name = families(i).name;
                    specs(n).feature_set = variants(j);
                    specs(n).model_role = role;
                    specs(n).predictors = predictors;
                end
            end
        end

        function [train_mask, test_mask] = conditionSplit(T, train_fraction, seed)
            rng(seed);
            condition_values = unique(T.condition_id, 'stable');
            idx = randperm(numel(condition_values));
            n_train = max(1, min(numel(condition_values) - 1, ...
                round(train_fraction * numel(condition_values))));
            train_conditions = condition_values(idx(1:n_train));
            train_mask = ismember(T.condition_id, train_conditions);
            test_mask = ~train_mask;
            assert(any(train_mask) && any(test_mask), ...
                'Condition split produced an empty train or test set.');
        end

        function tf = isGroupValue(x, value)
            if isstring(x) || iscategorical(x) || iscellstr(x)
                tf = string(x) == string(value);
            else
                tf = x == value;
            end
            tf = tf(:);
        end

        function assertGroupedSplit(T, train_mask, test_mask, heldout_var, heldout_value)
            assert(any(train_mask), 'Grouped split has empty train set.');
            assert(any(test_mask), 'Grouped split has empty test set.');
            assert(~any(train_mask & test_mask), 'Train/test masks overlap.');

            x = T.(char(heldout_var));
            assert(~any(adaptive_req.analysis.Test12Analysis.isGroupValue( ...
                x(train_mask), heldout_value)), ...
                'Heldout value appears in train split for %s = %s.', ...
                heldout_var, adaptive_req.analysis.Test12Analysis.valueToString(heldout_value));
            assert(all(adaptive_req.analysis.Test12Analysis.isGroupValue( ...
                x(test_mask), heldout_value)), ...
                'Test split is not a pure heldout group.');
        end

        function T = addModelMetadata(T, spec, varargin)
            p = inputParser;
            addParameter(p, 'GeneralizationTest', "", @(x) ischar(x) || isstring(x));
            addParameter(p, 'HeldoutVar', "", @(x) ischar(x) || isstring(x));
            addParameter(p, 'HeldoutValue', "", @(x) ischar(x) || isstring(x));
            addParameter(p, 'NTrain', NaN, @isnumeric);
            addParameter(p, 'NTest', NaN, @isnumeric);
            parse(p, varargin{:});

            T.model_name = repmat(string(spec.model_name), height(T), 1);
            T.feature_set = repmat(string(spec.feature_set), height(T), 1);
            T.model_role = repmat(string(spec.model_role), height(T), 1);

            if strlength(string(p.Results.GeneralizationTest)) > 0
                T.generalization_test = repmat(string(p.Results.GeneralizationTest), height(T), 1);
                T.heldout_var = repmat(string(p.Results.HeldoutVar), height(T), 1);
                T.heldout_value = repmat(string(p.Results.HeldoutValue), height(T), 1);
                T.N_train = double(p.Results.NTrain) * ones(height(T), 1);
                T.N_test = double(p.Results.NTest) * ones(height(T), 1);
            end
        end

        function T_sws = addSwsMetrics(T_pred, T_ref)
            [tf, loc] = ismember(string(T_pred.row_key), string(T_ref.row_key));
            if ~all(tf)
                error('Could not match all prediction rows back to reference rows.');
            end

            T_sws = T_pred;
            T_sws.SIM_f0 = T_ref.SIM_f0(loc);
            T_sws.SIM_cs_bg = T_ref.SIM_cs_bg(loc);
            if ismember('SIM_WaveModel', string(T_ref.Properties.VariableNames))
                T_sws.SIM_WaveModel = string(T_ref.SIM_WaveModel(loc));
            end
            T_sws.REQ_M = T_ref.REQ_M(loc);
            if ismember('REQ_cs_guess', string(T_ref.Properties.VariableNames))
                T_sws.REQ_cs_guess = T_ref.REQ_cs_guess(loc);
            end
            if ismember('M_eff_guess', string(T_ref.Properties.VariableNames))
                T_sws.M_eff_guess = T_ref.M_eff_guess(loc);
            end
            if ismember('M_eff_true_diag', string(T_ref.Properties.VariableNames))
                T_sws.M_eff_true_diag = T_ref.M_eff_true_diag(loc);
            end
            if ismember('step_idx', string(T_ref.Properties.VariableNames))
                T_sws.step_idx = T_ref.step_idx(loc);
            end
            if ismember('Omega_sr', string(T_ref.Properties.VariableNames))
                T_sws.Omega_sr = T_ref.Omega_sr(loc);
            end
            if ismember('q_local_minus_global', string(T_ref.Properties.VariableNames))
                T_sws.q_local_minus_global = T_ref.q_local_minus_global(loc);
                T_sws.abs_q_local_minus_global = abs(T_sws.q_local_minus_global);
            end

            T_sws.cs_true = T_ref.SIM_cs_bg(loc);
            T_sws.cs_pred = adaptive_req.analysis.Test12Analysis.qToCs( ...
                T_sws.q_pred, T_ref.req_mapping(loc), T_ref.SIM_f0(loc));
            T_sws.sws_error = T_sws.cs_pred - T_sws.cs_true;
            T_sws.abs_sws_error = abs(T_sws.sws_error);
            T_sws.sws_error_pct = 100 * T_sws.sws_error ./ T_sws.cs_true;
            T_sws.abs_sws_error_pct = abs(T_sws.sws_error_pct);
        end

        function cs = qToCs(q, mappings, f0)
            q = double(q(:));
            f0 = double(f0(:));
            cs = nan(numel(q), 1);
            for i = 1:numel(q)
                if isempty(mappings{i}) || ~isfinite(q(i))
                    continue;
                end
                k = adaptive_req.quantile.quantile_to_k(mappings{i}, q(i));
                cs(i) = 2*pi*f0(i) ./ k;
            end
        end

        function T_sum = summarizeSws(T, group_vars)
            group_vars = string(group_vars(:));
            T = T(:, unique([group_vars; [
                "sws_error_pct"
                "abs_sws_error_pct"]], 'stable'));
            [G, T_sum] = findgroups(T(:, cellstr(group_vars)));
            T_sum.N = splitapply(@numel, T.abs_sws_error_pct, G);
            T_sum.MAPE_pct = splitapply(@(x) mean(x, 'omitnan'), T.abs_sws_error_pct, G);
            T_sum.RMSE_pct = splitapply(@(x) sqrt(mean(x.^2, 'omitnan')), T.sws_error_pct, G);
            T_sum.MedAE_pct = splitapply(@(x) median(x, 'omitnan'), T.abs_sws_error_pct, G);
            T_sum.bias_pct = splitapply(@(x) mean(x, 'omitnan'), T.sws_error_pct, G);
            T_sum.P95_abs_error_pct = splitapply(@(x) prctile(x, 95), T.abs_sws_error_pct, G);
            T_sum.Max_abs_error_pct = splitapply(@(x) max(x, [], 'omitnan'), T.abs_sws_error_pct, G);
            T_sum.HighError_gt10_pct = splitapply(@(x) 100 * mean(x > 10, 'omitnan'), T.abs_sws_error_pct, G);
            T_sum.HighError_gt20_pct = splitapply(@(x) 100 * mean(x > 20, 'omitnan'), T.abs_sws_error_pct, G);
            T_sum = sortrows(T_sum, 'MAPE_pct', 'ascend');
        end

        function T_delta = deltaVsNoCsGuess(T_metrics, group_vars)
            group_vars = string(group_vars(:));
            keys = unique(T_metrics(:, cellstr(group_vars)), 'rows');
            T_delta = table();
            for i = 1:height(keys)
                idx_key = true(height(T_metrics), 1);
                for j = 1:numel(group_vars)
                    v = group_vars(j);
                    idx_key = idx_key & string(T_metrics.(char(v))) == string(keys.(char(v))(i));
                end
                Ti = T_metrics(idx_key, :);
                base = Ti(Ti.feature_set == "NoCsGuess", :);
                if isempty(base)
                    continue;
                end
                variants = ["WithCsGuess"; "WithMeffGuess"; "WithCsGuessAndMeffGuess"];
                for k = 1:numel(variants)
                    Tv = Ti(Ti.feature_set == variants(k), :);
                    if isempty(Tv)
                        continue;
                    end
                    row = Tv(1, :);
                    row.baseline_feature_set = "NoCsGuess";
                    row.Delta_MAPE_vs_NoCsGuess = Tv.MAPE_pct(1) - base.MAPE_pct(1);
                    row.Delta_HighError_gt20_vs_NoCsGuess = Tv.HighError_gt20_pct(1) - base.HighError_gt20_pct(1);
                    row.Delta_RMSE_vs_NoCsGuess = Tv.RMSE_pct(1) - base.RMSE_pct(1);
                    T_delta = adaptive_req.analysis.Test12Analysis.concatTables(T_delta, row);
                end
            end
        end

        function T = addBins(T)
            if ismember('M_eff_guess', string(T.Properties.VariableNames))
                T.M_eff_guess_bin = discretize(T.M_eff_guess, ...
                    [0 1.5 2 2.5 3 3.5 4 4.5 5 6 Inf], ...
                    'categorical', {'[0,1.5)', '[1.5,2)', '[2,2.5)', ...
                    '[2.5,3)', '[3,3.5)', '[3.5,4)', '[4,4.5)', ...
                    '[4.5,5)', '[5,6)', '[6,Inf)'});
            end
            if ismember('M_eff_true_diag', string(T.Properties.VariableNames))
                T.M_eff_true_diag_bin = discretize(T.M_eff_true_diag, ...
                    [0 1.5 2 2.5 3 3.5 4 4.5 5 6 Inf], ...
                    'categorical', {'[0,1.5)', '[1.5,2)', '[2,2.5)', ...
                    '[2.5,3)', '[3,3.5)', '[3.5,4)', '[4,4.5)', ...
                    '[4.5,5)', '[5,6)', '[6,Inf)'});
            end
        end

        function T = removeCellColumns(T)
            vars = T.Properties.VariableNames;
            remove = false(size(vars));
            for i = 1:numel(vars)
                remove(i) = iscell(T.(vars{i}));
            end
            T(:, remove) = [];
        end

        function T = concatTables(A, B)
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
            A = adaptive_req.analysis.Test12Analysis.addMissingColumns(A, vars_all);
            B = adaptive_req.analysis.Test12Analysis.addMissingColumns(B, vars_all);
            T = [A(:, cellstr(vars_all)); B(:, cellstr(vars_all))];
        end

        function plotMetricByTwoGroups(T, x_var, series_var, metric_var, title_text, file_path)
            if isempty(T)
                return;
            end
            x = string(T.(char(x_var)));
            series = string(T.(char(series_var)));
            x_values = unique(x, 'stable');
            series_values = unique(series, 'stable');
            Y = nan(numel(x_values), numel(series_values));
            for i = 1:numel(x_values)
                for j = 1:numel(series_values)
                    idx = x == x_values(i) & series == series_values(j);
                    if any(idx)
                        Y(i, j) = mean(T.(char(metric_var))(idx), 'omitnan');
                    end
                end
            end
            figure('Color', 'w', 'Position', [100 100 1050 560]);
            bar(categorical(x_values), Y);
            ylabel(strrep(string(metric_var), '_', '\_'));
            xlabel(strrep(string(x_var), '_', '\_'));
            title(title_text, 'Interpreter', 'none');
            legend(series_values, 'Location', 'best', 'Interpreter', 'none');
            grid on;
            exportgraphics(gcf, file_path, 'Resolution', 300, 'BackgroundColor', 'white');
            close(gcf);
        end

        function plotQScatter(T, title_text, file_path)
            T = T(T.model_type == "bagged_trees", :);
            if isempty(T)
                return;
            end
            models = unique(T.model_name, 'stable');
            features = unique(T.feature_set, 'stable');
            figure('Color', 'w', 'Position', [100 100 1350 850]);
            tl = tiledlayout(numel(models), numel(features), 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:numel(models)
                for j = 1:numel(features)
                    ax = nexttile(tl);
                    idx = T.model_name == models(i) & T.feature_set == features(j);
                    scatter(ax, T.q_true(idx), T.q_pred(idx), 6, 'filled', 'MarkerFaceAlpha', 0.25);
                    hold(ax, 'on');
                    plot(ax, [0 1], [0 1], 'k--');
                    axis(ax, 'equal');
                    xlim(ax, [0 1]);
                    ylim(ax, [0 1]);
                    grid(ax, 'on');
                    title(ax, sprintf('%s | %s', models(i), features(j)), 'Interpreter', 'none');
                    xlabel(ax, 'q true');
                    ylabel(ax, 'q predicted');
                end
            end
            title(tl, title_text, 'Interpreter', 'none');
            exportgraphics(gcf, file_path, 'Resolution', 300, 'BackgroundColor', 'white');
            close(gcf);
        end

        function plotErrorBox(T, title_text, file_path)
            T = T(T.model_type == "bagged_trees", :);
            if isempty(T)
                return;
            end
            labels = T.model_name + " | " + T.feature_set;
            figure('Color', 'w', 'Position', [100 100 1300 560]);
            boxchart(categorical(labels), T.abs_sws_error_pct);
            yline(20, 'r--', '20%');
            ylabel('|SWS error| (%)');
            title(title_text, 'Interpreter', 'none');
            xtickangle(30);
            grid on;
            exportgraphics(gcf, file_path, 'Resolution', 300, 'BackgroundColor', 'white');
            close(gcf);
        end

        function s = valueToString(x)
            s = string(x);
        end

        function key = makeRowKey(T)
            required = ["condition_id", "step_idx", "realization_idx", "patch_idx"];
            adaptive_req.analysis.Test12Analysis.requireVars(T, required, 'row_key');
            parts = strings(height(T), 4);
            parts(:, 1) = string(T.condition_id);
            parts(:, 2) = string(T.step_idx);
            parts(:, 3) = string(T.realization_idx);
            parts(:, 4) = string(T.patch_idx);
            key = join(parts, "|", 2);
        end

        function predictors = existingPredictors(T, candidates)
            vars = string(T.Properties.VariableNames);
            candidates = string(candidates(:));
            predictors = strings(0, 1);
            for i = 1:numel(candidates)
                name_i = candidates(i);
                if ismember(name_i, vars)
                    x = T.(char(name_i));
                    if isnumeric(x) || islogical(x) || isstring(x) || iscategorical(x) || iscellstr(x)
                        predictors(end + 1, 1) = name_i; %#ok<AGROW>
                    end
                end
            end
            predictors = unique(predictors, 'stable');
        end

        function predictors = filterOperationalPredictors(predictors)
            predictors = string(predictors(:));
            lower_names = lower(predictors);
            exact_banned = lower([
                "q_theory"
                "q_true"
                "q_global_theory"
                "q_local_minus_global"
                "abs_q_local_minus_global"
                "m_eff_true_diag"
                "cs_true"
                "cs_pred"
                "sws_error"
                "abs_sws_error"
                "sws_error_pct"
                "abs_sws_error_pct"
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

        function assertNoLeakage(predictors, model_name)
            filtered = adaptive_req.analysis.Test12Analysis.filterOperationalPredictors(predictors);
            assert(numel(filtered) == numel(predictors), ...
                'Operational model %s contains leakage predictors.', model_name);
            assert(~isempty(predictors), ...
                'Operational model %s has no predictors.', model_name);
        end
    end

    methods (Static, Access = private)
        function T = addLocalEcumFeatures(T)
            if ~ismember('req_mapping', string(T.Properties.VariableNames))
                return;
            end
            idx = find(~cellfun(@isempty, T.req_mapping), 1, 'first');
            if isempty(idx)
                return;
            end
            feat0 = adaptive_req.quantile.extract_ecum_shape_features(T.req_mapping{idx});
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
                feat_i = adaptive_req.quantile.extract_ecum_shape_features(T.req_mapping{i});
                for j = 1:numel(names)
                    T.(names{j})(i) = feat_i.(names{j});
                end
            end
        end

        function T = addMissingColumns(T, vars_all)
            vars = string(T.Properties.VariableNames);
            for i = 1:numel(vars_all)
                name_i = char(vars_all(i));
                if ismember(vars_all(i), vars)
                    continue;
                end
                string_like = any(endsWith(vars_all(i), ...
                    ["name", "type", "role", "source", "label", "split", ...
                    "test", "var", "value", "set", "guess"])) || ...
                    startsWith(vars_all(i), "SIM_WaveModel") || ...
                    contains(vars_all(i), "guess");
                if string_like
                    T.(name_i) = strings(height(T), 1);
                else
                    T.(name_i) = nan(height(T), 1);
                end
            end
        end
    end
end
