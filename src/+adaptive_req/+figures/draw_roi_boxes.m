function draw_roi_boxes(ax, roi_specs, varargin)
%DRAW_ROI_BOXES Draw ROI rectangles using meter-valued roi specs.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.draw_roi_boxes';
addRequired(p, 'ax');
addRequired(p, 'roi_specs');
addParameter(p, 'Color', 'w');
addParameter(p, 'LineWidth', 1.4);
addParameter(p, 'FontSize', 9);
addParameter(p, 'ShowLabels', true, @(x) islogical(x) || isnumeric(x));
parse(p, ax, roi_specs, varargin{:});

hold(ax, 'on');

for i = 1:numel(roi_specs)
    x_cm = roi_specs(i).xlim_m * 100;
    z_cm = roi_specs(i).zlim_m * 100;

    rectangle(ax, 'Position', [x_cm(1), z_cm(1), ...
        diff(x_cm), diff(z_cm)], ...
        'EdgeColor', p.Results.Color, ...
        'LineWidth', p.Results.LineWidth, ...
        'LineStyle', '-');

    if logical(p.Results.ShowLabels)
        text(ax, mean(x_cm), z_cm(2) + 0.06, string(roi_specs(i).label), ...
            'Color', p.Results.Color, ...
            'FontWeight', 'bold', ...
            'FontSize', p.Results.FontSize, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none');
    end
end

end
