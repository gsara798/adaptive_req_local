function FIGS = plot_q_model_performance(T_pred, T_metrics, varargin)
%PLOT_Q_MODEL_PERFORMANCE Plot true versus predicted q for trained models.
%
% Usage
%   FIGS = adaptive_req.analysis.plot_q_model_performance(T_pred, T_metrics);

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.plot_q_model_performance';

addRequired(p, 'T_pred', @istable);
addRequired(p, 'T_metrics', @istable);

addParameter(p, 'Split', 'test', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveFigure', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'OutputDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'BaseName', 'model1_features_only', @(x) ischar(x) || isstring(x));
addParameter(p, 'SavePNG', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveFIG', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Resolution', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CloseAfterSave', false, @(x) islogical(x) || isnumeric(x));

parse(p, T_pred, T_metrics, varargin{:});

split_name = string(p.Results.Split);
visible = logical(p.Results.Visible);

save_figure = logical(p.Results.SaveFigure);
output_dir = p.Results.OutputDir;
base_name = char(p.Results.BaseName);

model_types = unique(T_pred.model_type, 'stable');
n_models = numel(model_types);

%% Figure 1: predicted versus true q

if visible
    fig1 = figure('Color', 'w');
else
    fig1 = figure('Color', 'w', 'Visible', 'off');
end

tl = tiledlayout(fig1, 1, n_models, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

plot_label = "q prediction";

if ismember('model_name', T_pred.Properties.VariableNames)
    model_names = unique(T_pred.model_name, 'stable');
    if numel(model_names) == 1
        plot_label = model_names(1);
    end
end

title(tl, sprintf('%s, split = %s', plot_label, split_name), ...
    'Interpreter', 'none');

for i = 1:n_models

    model_i = model_types(i);

    ax = nexttile(tl);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    mask = T_pred.model_type == model_i & T_pred.split == split_name;

    q_true = T_pred.q_true(mask);
    q_pred = T_pred.q_pred(mask);

    scatter(ax, q_true, q_pred, 24, 'filled', ...
        'MarkerFaceAlpha', 0.45, ...
        'MarkerEdgeAlpha', 0.45);

    plot(ax, [0 1], [0 1], 'k-', 'LineWidth', 1.2);

    xlabel(ax, 'true q');
    ylabel(ax, 'predicted q');

    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    axis(ax, 'square');

    metric_mask = T_metrics.model_type == model_i & T_metrics.split == split_name;

    if any(metric_mask)

        M = T_metrics(metric_mask, :);

        txt = sprintf('RMSE = %.4f\nMAE = %.4f\nR^2 = %.4f\n\\rho_s = %.4f', ...
            M.RMSE(1), ...
            M.MAE(1), ...
            M.R2(1), ...
            M.spearman_rho(1));

        text(ax, 0.05, 0.95, txt, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', ...
            'Margin', 3);
    end

    title(ax, char(model_i), 'Interpreter', 'none');

end

%% Figure 2: residuals

if visible
    fig2 = figure('Color', 'w');
else
    fig2 = figure('Color', 'w', 'Visible', 'off');
end

tl2 = tiledlayout(fig2, 1, n_models, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tl2, sprintf('%s residuals, split = %s', plot_label, split_name), ...
    'Interpreter', 'none');

for i = 1:n_models

    model_i = model_types(i);

    ax = nexttile(tl2);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    mask = T_pred.model_type == model_i & T_pred.split == split_name;

    q_true = T_pred.q_true(mask);
    residual = T_pred.residual(mask);

    scatter(ax, q_true, residual, 24, 'filled', ...
        'MarkerFaceAlpha', 0.45, ...
        'MarkerEdgeAlpha', 0.45);

    yline(ax, 0, 'k-', 'LineWidth', 1.2);

    xlabel(ax, 'true q');
    ylabel(ax, 'prediction residual');

    xlim(ax, [0 1]);

    title(ax, char(model_i), 'Interpreter', 'none');

end

FIGS = struct();
FIGS.predicted_vs_true = fig1;
FIGS.residuals = fig2;

if save_figure

    adaptive_req.analysis.save_analysis_figure( ...
        fig1, ...
        output_dir, ...
        [base_name '_predicted_vs_true_' char(split_name)], ...
        'SavePNG', p.Results.SavePNG, ...
        'SavePDF', p.Results.SavePDF, ...
        'SaveFIG', p.Results.SaveFIG, ...
        'Resolution', p.Results.Resolution, ...
        'CloseAfterSave', p.Results.CloseAfterSave);

    adaptive_req.analysis.save_analysis_figure( ...
        fig2, ...
        output_dir, ...
        [base_name '_residuals_' char(split_name)], ...
        'SavePNG', p.Results.SavePNG, ...
        'SavePDF', p.Results.SavePDF, ...
        'SaveFIG', p.Results.SaveFIG, ...
        'Resolution', p.Results.Resolution, ...
        'CloseAfterSave', p.Results.CloseAfterSave);

end

end
