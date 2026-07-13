function apply_style(fig, ax, varargin)
%APPLY_STYLE Apply consistent figure style for adaptive_req plots.

p = inputParser;

addRequired(p, 'fig');
addRequired(p, 'ax');

addParameter(p, 'FontSize', 14, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'LineWidth', 1.2, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FontName', 'Avenir', @(x) ischar(x) || isstring(x));

parse(p, fig, ax, varargin{:});

fs = p.Results.FontSize;
lw = p.Results.LineWidth;
font_name = char(p.Results.FontName);

try
    adaptive_req.figures.apply_paper_style(fig, ax);
catch
    set(ax, ...
        'FontSize', fs, ...
        'FontName', font_name, ...
        'LineWidth', lw, ...
        'Box', 'on');
end

if isgraphics(fig, 'figure')
    set(fig, 'Color', 'w');
end

end