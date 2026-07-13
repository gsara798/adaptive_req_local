function saved_files = save_analysis_figure(fig, output_dir, base_name, varargin)
%SAVE_ANALYSIS_FIGURE Save an analysis figure as PNG, PDF, and optionally FIG.
%
% Usage
%   saved_files = adaptive_req.analysis.save_analysis_figure( ...
%       fig, output_dir, base_name, ...
%       'SavePNG', true, ...
%       'SavePDF', true, ...
%       'SaveFIG', false);
%
% Inputs
%   fig:
%       Figure handle.
%
%   output_dir:
%       Folder where the figure will be saved.
%
%   base_name:
%       File name without extension.
%
% Name-value options
%   SavePNG:
%       Save PNG file.
%
%   SavePDF:
%       Save PDF file.
%
%   SaveFIG:
%       Save MATLAB FIG file.
%
%   Resolution:
%       PNG resolution in dpi.
%
%   CloseAfterSave:
%       Close figure after saving.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.save_analysis_figure';

addRequired(p, 'fig', @(x) ishghandle(x, 'figure'));
addRequired(p, 'output_dir', @(x) ischar(x) || isstring(x));
addRequired(p, 'base_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'SavePNG', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveFIG', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Resolution', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CloseAfterSave', false, @(x) islogical(x) || isnumeric(x));

parse(p, fig, output_dir, base_name, varargin{:});

output_dir = char(p.Results.output_dir);
base_name = char(p.Results.base_name);

save_png = logical(p.Results.SavePNG);
save_pdf = logical(p.Results.SavePDF);
save_fig = logical(p.Results.SaveFIG);
resolution = p.Results.Resolution;
close_after_save = logical(p.Results.CloseAfterSave);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

base_name = sanitize_file_name(base_name);

saved_files = strings(0, 1);

if save_png
    png_file = fullfile(output_dir, [base_name '.png']);
    exportgraphics(fig, png_file, 'Resolution', resolution);
    saved_files(end + 1, 1) = string(png_file);
end

if save_pdf
    pdf_file = fullfile(output_dir, [base_name '.pdf']);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
    saved_files(end + 1, 1) = string(pdf_file);
end

if save_fig
    fig_file = fullfile(output_dir, [base_name '.fig']);
    savefig(fig, fig_file);
    saved_files(end + 1, 1) = string(fig_file);
end

if close_after_save
    close(fig);
end

end

function name = sanitize_file_name(name)

name = string(name);
name = regexprep(name, '[^A-Za-z0-9_\-]+', '_');
name = char(name);

end
