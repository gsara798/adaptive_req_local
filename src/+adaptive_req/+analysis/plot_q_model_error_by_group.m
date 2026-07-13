function fig = plot_q_model_error_by_group(T_err, group_var, varargin)
%PLOT_Q_MODEL_ERROR_BY_GROUP Plot model error metric by group.
%
% Usage
%   fig = adaptive_req.analysis.plot_q_model_error_by_group( ...
%       T_err, ...
%       'REQ_M', ...
%       'Metric', 'RMSE');

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.plot_q_model_error_by_group';

addRequired(p, 'T_err', @istable);
addRequired(p, 'group_var', @(x) ischar(x) || isstring(x));

addParameter(p, 'Metric', 'RMSE', @(x) ischar(x) || isstring(x));
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

parse(p, T_err, group_var, varargin{:});

group_var = char(group_var);
metric = char(p.Results.Metric);

if ~ismember(group_var, T_err.Properties.VariableNames)
    error('Group variable not found: %s', group_var);
end

if ~ismember(metric, T_err.Properties.VariableNames)
    error('Metric variable not found: %s', metric);
end

if ~ismember('model_type', T_err.Properties.VariableNames)
    error('T_err must contain model_type.');
end

if logical(p.Results.Visible)
    fig = figure('Color', 'w');
else
    fig = figure('Color', 'w', 'Visible', 'off');
end

ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

model_types = unique(T_err.model_type, 'stable');

x_raw = T_err.(group_var);

is_numeric_group = isnumeric(x_raw) || islogical(x_raw);

legend_entries = strings(numel(model_types), 1);

if is_numeric_group

    for i = 1:numel(model_types)

        model_i = model_types(i);
        mask = T_err.model_type == model_i;

        x = T_err.(group_var)(mask);
        y = T_err.(metric)(mask);

        [x, idx] = sort(x);
        y = y(idx);

        plot(ax, x, y, '-o', ...
            'LineWidth', 2.0, ...
            'MarkerSize', 7);

        legend_entries(i) = model_i;

    end

    xlabel(ax, strrep(group_var, '_', ' '), 'Interpreter', 'none');

else

    group_values = unique(string(T_err.(group_var)), 'stable');
    x = 1:numel(group_values);

    for i = 1:numel(model_types)

        model_i = model_types(i);

        y = NaN(size(x));

        for j = 1:numel(group_values)

            mask = T_err.model_type == model_i & ...
                   string(T_err.(group_var)) == group_values(j);

            if any(mask)
                y(j) = T_err.(metric)(find(mask, 1));
            end
        end

        plot(ax, x, y, '-o', ...
            'LineWidth', 2.0, ...
            'MarkerSize', 7);

        legend_entries(i) = model_i;

    end

    xticks(ax, x);
    xticklabels(ax, group_values);
    xtickangle(ax, 30);

    xlabel(ax, strrep(group_var, '_', ' '), 'Interpreter', 'none');

end

ylabel(ax, metric, 'Interpreter', 'none');

plot_title = char(p.Results.Title);

if isempty(plot_title)
    plot_title = sprintf('%s by %s', metric, group_var);
end

title(ax, plot_title, 'Interpreter', 'none');

legend(ax, legend_entries, ...
    'Location', 'best', ...
    'Interpreter', 'none');

base_name = char(p.Results.BaseName);

if isempty(base_name)
    base_name = sprintf('model_error_%s_by_%s', metric, group_var);
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
