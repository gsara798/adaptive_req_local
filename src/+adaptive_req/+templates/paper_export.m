function saved = paper_export(fig, fig_dir, out_name, varargin)
%PAPER_EXPORT Export PNG/PDF with stable paper-like sizing.

p = inputParser;
p.FunctionName = 'adaptive_req.templates.paper_export';

addRequired(p, 'fig');
addRequired(p, 'fig_dir', @(x) ischar(x) || isstring(x));
addRequired(p, 'out_name', @(x) ischar(x) || isstring(x));
addParameter(p, 'DPI', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'VectorPDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'WidthCm', 18, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'HeightCm', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, fig, fig_dir, out_name, varargin{:});

fig_dir = char(p.Results.fig_dir);
out_name = string(p.Results.out_name);

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

set(fig, 'Units', 'centimeters');
fig.Position(3:4) = [p.Results.WidthCm, p.Results.HeightCm];
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [p.Results.WidthCm, p.Results.HeightCm]);

png_path = fullfile(fig_dir, out_name + ".png");
pdf_path = fullfile(fig_dir, out_name + ".pdf");

exportgraphics(fig, png_path, ...
    'Resolution', p.Results.DPI, ...
    'BackgroundColor', 'white');

if logical(p.Results.VectorPDF)
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
else
    exportgraphics(fig, pdf_path, 'Resolution', p.Results.DPI);
end

saved = [png_path; pdf_path];

end
