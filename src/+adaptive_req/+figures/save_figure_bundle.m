function saved_files = save_figure_bundle(figs, output_dir, base_name, varargin)
%SAVE_FIGURE_BUNDLE Save one figure or a struct of figures as PNG, PDF, and FIG.
%
% Usage:
%   adaptive_req.figures.save_figure_bundle(fig, output_dir, 'my_figure');
%
%   adaptive_req.figures.save_figure_bundle(figs, output_dir, 'step_01', ...
%       'SavePNG', true, ...
%       'SavePDF', true, ...
%       'SaveFIG', true);

p = inputParser;

addRequired(p, 'figs');
addRequired(p, 'output_dir', @(x) ischar(x) || isstring(x));
addRequired(p, 'base_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'SavePNG', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveFIG', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Resolution', 300, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'CloseAfterSave', false, @(x) islogical(x) || isnumeric(x));

parse(p, figs, output_dir, base_name, varargin{:});

output_dir = char(output_dir);
base_name = sanitize_filename(char(base_name));

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

saved_files = strings(0, 1);

if isstruct(figs)
    names = fieldnames(figs);

    for i = 1:numel(names)
        name_i = names{i};
        fig_i = figs.(name_i);

        if isgraphics(fig_i, 'figure')
            base_i = sprintf('%s_%s', base_name, name_i);

            files_i = save_one_figure(fig_i, output_dir, base_i, p.Results);
            saved_files = [saved_files; files_i(:)]; %#ok<AGROW>

            if logical(p.Results.CloseAfterSave)
                close(fig_i);
            end
        end
    end

elseif isgraphics(figs, 'figure')

    files_i = save_one_figure(figs, output_dir, base_name, p.Results);
    saved_files = [saved_files; files_i(:)];

    if logical(p.Results.CloseAfterSave)
        close(figs);
    end

else
    error('Input figs must be a figure handle or a struct containing figure handles.');
end

end

function saved_files = save_one_figure(fig, output_dir, base_name, opt)

saved_files = strings(0, 1);

if logical(opt.SavePNG)
    png_file = fullfile(output_dir, [base_name '.png']);
    exportgraphics(fig, png_file, 'Resolution', opt.Resolution);
    saved_files(end+1, 1) = string(png_file);
end

if logical(opt.SavePDF)
    pdf_file = fullfile(output_dir, [base_name '.pdf']);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
    saved_files(end+1, 1) = string(pdf_file);
end

if logical(opt.SaveFIG)
    fig_file = fullfile(output_dir, [base_name '.fig']);
    savefig(fig, fig_file);
    saved_files(end+1, 1) = string(fig_file);
end

end

function safe_name = sanitize_filename(name)

safe_name = char(string(name));
safe_name = regexprep(safe_name, '[^a-zA-Z0-9_\-.]', '_');
safe_name = regexprep(safe_name, '_+', '_');

end