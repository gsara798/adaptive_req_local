function fig = plot_q_vs_omega_grouped(T_step, varargin)
%PLOT_Q_VS_OMEGA_GROUPED Plot q versus aperture solid angle.
%
% Recommended input
%   T_step:
%       Condition-step summary table produced by
%       adaptive_req.analysis.make_step_summary_table.
%
% Default plot
%   x = omega_mean
%   y = q_mean
%   error bars = q_std
%   color/group = REQ_M
%   facets = SIM_WaveModel
%
% Usage
%   fig = adaptive_req.analysis.plot_q_vs_omega_grouped(T_step);
%
%   fig = adaptive_req.analysis.plot_q_vs_omega_grouped( ...
%       T_step, ...
%       'ColorBy', 'REQ_M', ...
%       'FacetBy', 'SIM_WaveModel');

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.plot_q_vs_omega_grouped';

addRequired(p, 'T_step', @istable);

addParameter(p, 'XVar', 'omega_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'YVar', 'q_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'YErrVar', 'q_std', @(x) ischar(x) || isstring(x));

addParameter(p, 'ColorBy', 'REQ_M', @(x) ischar(x) || isstring(x));
addParameter(p, 'FacetBy', 'SIM_WaveModel', @(x) ischar(x) || isstring(x));

addParameter(p, 'UseErrorBars', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Title', 'q versus aperture', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveFigure', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'OutputDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'BaseName', 'q_vs_omega_grouped', @(x) ischar(x) || isstring(x));
addParameter(p, 'SavePNG', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveFIG', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Resolution', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CloseAfterSave', false, @(x) islogical(x) || isnumeric(x));

parse(p, T_step, varargin{:});

x_var = char(p.Results.XVar);
y_var = char(p.Results.YVar);
yerr_var = char(p.Results.YErrVar);

color_by = char(p.Results.ColorBy);
facet_by = char(p.Results.FacetBy);

use_error_bars = logical(p.Results.UseErrorBars);
visible = logical(p.Results.Visible);

validate_table_variable(T_step, x_var);
validate_table_variable(T_step, y_var);

has_yerr = ismember(yerr_var, T_step.Properties.VariableNames);

if use_error_bars && ~has_yerr
    warning('YErrVar not found. Plotting without error bars: %s', yerr_var);
    use_error_bars = false;
end

has_color = ~isempty(color_by) && ismember(color_by, T_step.Properties.VariableNames);
has_facet = ~isempty(facet_by) && ismember(facet_by, T_step.Properties.VariableNames);

if has_facet
    facet_values = unique_values(T_step.(facet_by));
else
    facet_values = "__all__";
end

n_facets = numel(facet_values);

n_cols = ceil(sqrt(n_facets));
n_rows = ceil(n_facets / n_cols);

if visible
    fig = figure('Color', 'w');
else
    fig = figure('Color', 'w', 'Visible', 'off');
end

tl = tiledlayout(fig, n_rows, n_cols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tl, char(p.Results.Title), 'Interpreter', 'none');

for fidx = 1:n_facets

    ax = nexttile(tl);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    if has_facet
        facet_mask = equals_value(T_step.(facet_by), facet_values(fidx));
        facet_label = sprintf('%s = %s', facet_by, value_to_label(facet_values(fidx)));
    else
        facet_mask = true(height(T_step), 1);
        facet_label = '';
    end

    if has_color
        color_values = unique_values(T_step.(color_by));
    else
        color_values = "__all__";
    end

    legend_entries = strings(0, 1);

    for cidx = 1:numel(color_values)

        if has_color
            color_mask = equals_value(T_step.(color_by), color_values(cidx));
            group_label = sprintf('%s = %s', color_by, value_to_label(color_values(cidx)));
        else
            color_mask = true(height(T_step), 1);
            group_label = "all";
        end

        mask = facet_mask & color_mask;

        if ~any(mask)
            continue;
        end

        x = T_step.(x_var)(mask);
        y = T_step.(y_var)(mask);

        [x, sort_idx] = sort(x);
        y = y(sort_idx);

        if use_error_bars
            yerr = T_step.(yerr_var)(mask);
            yerr = yerr(sort_idx);

            errorbar(ax, x, y, yerr, ...
                '-o', ...
                'LineWidth', 1.8, ...
                'MarkerSize', 6, ...
                'CapSize', 8);
        else
            plot(ax, x, y, '-o', ...
                'LineWidth', 1.8, ...
                'MarkerSize', 6);
        end

        legend_entries(end + 1, 1) = string(group_label);
    end

    xlabel(ax, format_axis_label(x_var), 'Interpreter', 'none');
    ylabel(ax, format_axis_label(y_var), 'Interpreter', 'none');

    if ~isempty(facet_label)
        title(ax, facet_label, 'Interpreter', 'none');
    end

    ylim(ax, [0 1]);

    if ~isempty(legend_entries)
        legend(ax, legend_entries, ...
            'Location', 'best', ...
            'Interpreter', 'none');
    end
end

if logical(p.Results.SaveFigure)
    adaptive_req.analysis.save_analysis_figure( ...
        fig, ...
        p.Results.OutputDir, ...
        p.Results.BaseName, ...
        'SavePNG', p.Results.SavePNG, ...
        'SavePDF', p.Results.SavePDF, ...
        'SaveFIG', p.Results.SaveFIG, ...
        'Resolution', p.Results.Resolution, ...
        'CloseAfterSave', p.Results.CloseAfterSave);
end

end

function validate_table_variable(T, var_name)

if ~ismember(var_name, T.Properties.VariableNames)
    error('Variable not found in table: %s', var_name);
end

end

function values = unique_values(x)

if iscell(x)
    values = unique(string(x), 'stable');
elseif isstring(x)
    values = unique(x, 'stable');
elseif iscategorical(x)
    values = unique(x, 'stable');
else
    values = unique(x, 'stable');
end

end

function mask = equals_value(x, value)

if iscell(x)
    mask = string(x) == string(value);
elseif isstring(x)
    mask = x == string(value);
elseif iscategorical(x)
    mask = x == value;
else
    mask = x == value;
end

end

function label = value_to_label(value)

if isstring(value) || ischar(value)
    label = char(string(value));
elseif iscategorical(value)
    label = char(string(value));
elseif isnumeric(value) || islogical(value)
    label = sprintf('%.6g', value);
else
    label = char(string(value));
end

end

function label = format_axis_label(var_name)

label = strrep(char(var_name), '_', ' ');

end
