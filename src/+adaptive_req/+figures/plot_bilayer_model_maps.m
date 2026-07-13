function fig = plot_bilayer_model_maps(T_pred, roi_specs, varargin)
%PLOT_BILAYER_MODEL_MAPS Plot model maps from patch-center predictions.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.plot_bilayer_model_maps';
addRequired(p, 'T_pred', @istable);
addRequired(p, 'roi_specs');
addParameter(p, 'ValueVar', 'cs_pred', @(x) ischar(x) || isstring(x));
addParameter(p, 'Models', string.empty(0, 1), @(x) isstring(x) || iscellstr(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ColorLimits', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ColorbarLabel', 'c_s (m/s)', @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureSize', [34 11], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FontSize', 13, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ShowRoiLabels', true, @(x) islogical(x) || isnumeric(x));
parse(p, T_pred, roi_specs, varargin{:});

models = string(p.Results.Models);
if isempty(models)
    models = unique(string(T_pred.model_name), 'stable');
end

fig = figure('Color', 'w');
tl = tiledlayout(fig, 1, numel(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(models)
    ax = nexttile(tl);
    Ti = T_pred(string(T_pred.model_name) == models(i), :);
    [map_i, xc, zc] = adaptive_req.figures.table_to_map_grid( ...
        Ti, p.Results.ValueVar);
    imagesc(ax, xc * 100, zc * 100, map_i);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    xlim(ax, [min(xc), max(xc)] * 100);
    ylim(ax, [min(zc), max(zc)] * 100);
    colormap(ax, parula);
    if ~isempty(p.Results.ColorLimits)
        clim(ax, p.Results.ColorLimits);
    end
    cb = colorbar(ax);
    cb.Label.String = char(p.Results.ColorbarLabel);
    cb.Label.Interpreter = 'tex';
    title(ax, readable_model_name(models(i)), ...
        'Interpreter', 'none', 'FontSize', p.Results.FontSize + 1);
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');
    adaptive_req.figures.draw_roi_boxes(ax, roi_specs, ...
        'FontSize', max(8, p.Results.FontSize - 3), ...
        'ShowLabels', p.Results.ShowRoiLabels);
end

if strlength(string(p.Results.Title)) > 0
    title(tl, string(p.Results.Title), 'Interpreter', 'tex', ...
        'FontSize', p.Results.FontSize + 2);
end

adaptive_req.templates.apply_paper_style(fig, [], ...
    'Times New Roman', p.Results.FontSize);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = p.Results.FigureSize;

end

function name = readable_model_name(name)

name = string(name);
name = replace(name, "TheoryDiffuse3D", "Theory diffuse 3D");
name = replace(name, "GlobalQSingleModel", "Global q model");
name = replace(name, "LocalOnly", "Local only");
name = replace(name, "HybridLocalGlobal", "Hybrid");

end
