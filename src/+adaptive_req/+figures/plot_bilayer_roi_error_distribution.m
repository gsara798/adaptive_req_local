function fig = plot_bilayer_roi_error_distribution(T_pred, varargin)
%PLOT_BILAYER_ROI_ERROR_DISTRIBUTION Boxplots of ROI patch errors.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.plot_bilayer_roi_error_distribution';
addRequired(p, 'T_pred', @istable);
addParameter(p, 'Title', 'ROI patch error distributions', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureSize', [28 10], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FontSize', 12, @(x) isnumeric(x) && isscalar(x));
parse(p, T_pred, varargin{:});

T = T_pred(T_pred.roi_name ~= "outside_roi" & ...
    isfinite(T_pred.cs_abs_error_pct), :);
labels = readable_model_name(T.model_name) + " | " + ...
    readable_roi_name(T.roi_name);
labels = categorical(labels);
labels = reordercats(labels, categories(labels));

fig = figure('Color', 'w');
ax = axes(fig);
boxchart(ax, labels, T.cs_abs_error_pct, ...
    'BoxFaceColor', [0.65 0.82 0.95], ...
    'MarkerColor', [0 0.4470 0.7410]);
ylabel(ax, '|c_s error| (%)', 'Interpreter', 'tex');
title(ax, string(p.Results.Title), 'Interpreter', 'tex');
grid(ax, 'on');
xtickangle(ax, 25);

adaptive_req.templates.apply_paper_style(fig, ax, ...
    'Times New Roman', p.Results.FontSize);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = p.Results.FigureSize;

end

function out = readable_model_name(models)

out = string(models);
out = replace(out, "TheoryDiffuse3D", "Theory");
out = replace(out, "GlobalQSingleModel", "Global q");
out = replace(out, "LocalOnly", "Local");
out = replace(out, "HybridLocalGlobal", "Hybrid");

end

function out = readable_roi_name(rois)

out = string(rois);
out = replace(out, "soft_roi", "soft");
out = replace(out, "hard_roi", "hard");
out = replace(out, "soft", "soft");
out = replace(out, "hard", "hard");

end
