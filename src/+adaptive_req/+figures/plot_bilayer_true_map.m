function fig = plot_bilayer_true_map(sim, roi_specs, varargin)
%PLOT_BILAYER_TRUE_MAP Plot true bilayer map with valid centers and ROIs.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.plot_bilayer_true_map';
addRequired(p, 'sim', @isstruct);
addRequired(p, 'roi_specs');
addParameter(p, 'XCentersM', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ZCentersM', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Title', 'Bilayer true c_s map', @(x) ischar(x) || isstring(x));
addParameter(p, 'ColorLimits', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'FigureSize', [16 13], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FontSize', 13, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ShowPatchCenters', false, @(x) islogical(x) || isnumeric(x));
parse(p, sim, roi_specs, varargin{:});

fig = figure('Color', 'w');
ax = axes(fig);

imagesc(ax, sim.x * 100, sim.z * 100, sim.cs_map);
axis(ax, 'image');
set(ax, 'YDir', 'reverse');
colormap(ax, parula);
if ~isempty(p.Results.ColorLimits)
    clim(ax, p.Results.ColorLimits);
end
cb = colorbar(ax);
cb.Label.String = 'c_s (m/s)';
cb.Label.Interpreter = 'tex';
title(ax, string(p.Results.Title), 'Interpreter', 'tex', ...
    'FontSize', p.Results.FontSize + 2);
xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');

x_centers_m = p.Results.XCentersM;
z_centers_m = p.Results.ZCentersM;
if ~isempty(x_centers_m) && ~isempty(z_centers_m)
    xlim(ax, [min(x_centers_m), max(x_centers_m)] * 100);
    ylim(ax, [min(z_centers_m), max(z_centers_m)] * 100);
end

if logical(p.Results.ShowPatchCenters) && ...
        ~isempty(x_centers_m) && ~isempty(z_centers_m)
    hold(ax, 'on');
    [Xc, Zc] = meshgrid(unique(x_centers_m), unique(z_centers_m));
    scatter(ax, Xc(:) * 100, Zc(:) * 100, 5, 'w', 'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'MarkerEdgeAlpha', 0.45, ...
        'MarkerFaceAlpha', 0.75);
end

adaptive_req.figures.draw_roi_boxes(ax, roi_specs, ...
    'FontSize', max(8, p.Results.FontSize - 3), ...
    'ShowLabels', true);

adaptive_req.templates.apply_paper_style(fig, ax, ...
    'Times New Roman', p.Results.FontSize);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = p.Results.FigureSize;

end
