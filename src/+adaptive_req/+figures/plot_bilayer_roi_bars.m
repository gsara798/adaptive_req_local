function fig = plot_bilayer_roi_bars(T_roi, varargin)
%PLOT_BILAYER_ROI_BARS Plot ROI accuracy and precision summaries.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.plot_bilayer_roi_bars';
addRequired(p, 'T_roi', @istable);
addParameter(p, 'Title', 'ROI summary', @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureSize', [28 11], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FontSize', 12, @(x) isnumeric(x) && isscalar(x));
parse(p, T_roi, varargin{:});

T = T_roi;
T = T(T.roi_name ~= "outside_roi" & T.roi_name ~= "all", :);
[~, order] = sortrows(table(T.frequency_hz, string(T.model_name), ...
    string(T.roi_name)), [1 2 3]);
T = T(order, :);

label_text = readable_label(T.model_name, T.roi_name, T.frequency_hz);
labels = categorical(label_text);
labels = reordercats(labels, cellstr(unique(label_text, 'stable')));

fig = figure('Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
bar(ax, labels, T.MAPE_pct, 'FaceColor', [0 0.4470 0.7410]);
ylabel(ax, 'MAPE (%)');
title(ax, 'Accuracy', 'Interpreter', 'none');
grid(ax, 'on');
xtickangle(ax, 25);

ax = nexttile(tl);
bar(ax, labels, T.CoV_pct, 'FaceColor', [0 0.4470 0.7410]);
ylabel(ax, 'CoV (%)');
title(ax, 'Precision', 'Interpreter', 'none');
grid(ax, 'on');
xtickangle(ax, 25);

title(tl, string(p.Results.Title), 'Interpreter', 'tex', ...
    'FontSize', p.Results.FontSize + 2);

adaptive_req.templates.apply_paper_style(fig, [], ...
    'Times New Roman', p.Results.FontSize);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = p.Results.FigureSize;

end

function labels = readable_label(models, rois, frequencies)

models = string(models);
rois = string(rois);
labels = readable_model_name(models) + " | " + readable_roi_name(rois) + ...
    " | " + string(frequencies) + " Hz";

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
