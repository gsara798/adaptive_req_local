function fig = plot_q_vs_feature_grouped(T_step, feature_name, varargin)
%PLOT_Q_VS_FEATURE_GROUPED Plot q versus a feature.
%
% Recommended input
%   T_step:
%       Condition-step summary table.
%
%   feature_name:
%       Either the base feature name or the full summary column name.
%
% Examples
%   plot_q_vs_feature_grouped(T_step, 'angular_entropy')
%
%   plot_q_vs_feature_grouped(T_step, 'angular_entropy_mean')
%
% Default plot
%   x = <feature>_mean
%   y = q_mean
%   color/group = REQ_M
%   facets = SIM_WaveModel

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.plot_q_vs_feature_grouped';

addRequired(p, 'T_step', @istable);
addRequired(p, 'feature_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'YVar', 'q_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'ColorBy', 'REQ_M', @(x) ischar(x) || isstring(x));
addParameter(p, 'FacetBy', 'SIM_WaveModel', @(x) ischar(x) || isstring(x));

addParameter(p, 'ShowFit', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveFigure', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'OutputDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'BaseName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SavePNG', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SavePDF', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveFIG', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Resolution', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CloseAfterSave', false, @(x) islogical(x) || isnumeric(x));

parse(p, T_step, feature_name, varargin{:});

feature_name = char(feature_name);

x_var = resolve_feature_mean_variable(T_step, feature_name);
y_var = char(p.Results.YVar);

color_by = char(p.Results.ColorBy);
facet_by = char(p.Results.FacetBy);

show_fit = logical(p.Results.ShowFit);
visible = logical(p.Results.Visible);

validate_table_variable(T_step, x_var);
validate_table_variable(T_step, y_var);

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

plot_title = char(p.Results.Title);

if isempty(plot_title)
    plot_title = sprintf('q versus %s', strrep(x_var, '_', ' '));
end

title(tl, plot_title, 'Interpreter', 'none');

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

        scatter(ax, x, y, 45, 'filled', ...
            'MarkerFaceAlpha', 0.75, ...
            'MarkerEdgeAlpha', 0.75);

        legend_entries(end + 1, 1) = string(group_label);

        if show_fit
            add_linear_fit(ax, x, y);
            legend_entries(end + 1, 1) = string(group_label) + " fit";
        end
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

base_name = char(p.Results.BaseName);

if isempty(base_name)
    base_name = ['q_vs_feature_' x_var];
end

if logical(p.Results.SaveFigure)
    adaptive_req.analysis.save_analysis_figure( ...
        fig, ...
        p.Results.OutputDir, ...
        base_name, ...
        'SavePNG', p.Results.SavePNG, ...
        'SavePDF', p.Results.SavePDF, ...
        'SaveFIG', p.Results.SaveFIG, ...
        'Resolution', p.Results.Resolution, ...
        'CloseAfterSave', p.Results.CloseAfterSave);
end

end

function add_linear_fit(ax, x, y)

mask = isfinite(x) & isfinite(y);

x = x(mask);
y = y(mask);

if numel(x) < 3 || numel(unique(x)) < 2
    return;
end

p = polyfit(x, y, 1);

xfit = linspace(min(x), max(x), 100);
yfit = polyval(p, xfit);

plot(ax, xfit, yfit, '--', 'LineWidth', 1.5);

end

function x_var = resolve_feature_mean_variable(T, feature_name)

vars = string(T.Properties.VariableNames);
feature_name = string(feature_name);

candidates = [
    feature_name
    feature_name + "_mean"
];

idx = find(ismember(vars, candidates), 1);

if isempty(idx)
    error('Could not resolve feature mean variable for: %s', feature_name);
end

x_var = char(vars(idx));

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
