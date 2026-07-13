function fig = plot_q_vs_feature_controlled(T, feature_name, varargin)
%PLOT_Q_VS_FEATURE_CONTROLLED Plot q versus one feature for controlled groups.
%
% This function is intended for controlled visual analysis.
% The input table T should already be filtered to the subset of interest.
%
% Example:
%   fig = adaptive_req.analysis.plot_q_vs_feature_controlled( ...
%       T_sub, ...
%       'ang_entropy', ...
%       'GroupVars', ["SIM_WaveModel", "REQ_M"]);
%
% Each panel corresponds to one group.
% Inside each panel:
%   x-axis  = selected feature
%   y-axis  = q_mean
%   color   = omega_mean
%   line    = points connected in order of omega_mean
%   text    = Spearman rho within that panel

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.plot_q_vs_feature_controlled';

addRequired(p, 'T', @istable);
addRequired(p, 'feature_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'QVar', 'q_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'OmegaVar', 'omega_mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'GroupVars', ["SIM_WaveModel", "REQ_M"], ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

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

parse(p, T, feature_name, varargin{:});

feature_name = string(feature_name);
q_var = char(p.Results.QVar);
omega_var = char(p.Results.OmegaVar);
group_vars = string(p.Results.GroupVars);

x_var = resolve_feature_mean_variable(T, feature_name);

validate_table_variable(T, x_var);
validate_table_variable(T, q_var);
validate_table_variable(T, omega_var);

group_vars = group_vars(ismember(group_vars, string(T.Properties.VariableNames)));

if isempty(group_vars)
    error('No valid grouping variables were found.');
end

[G, T_groups] = findgroups(T(:, cellstr(group_vars)));
n_groups = height(T_groups);

n_cols = ceil(sqrt(n_groups));
n_rows = ceil(n_groups / n_cols);

if logical(p.Results.Visible)
    fig = figure('Color', 'w');
else
    fig = figure('Color', 'w', 'Visible', 'off');
end

tl = tiledlayout(fig, n_rows, n_cols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

plot_title = string(p.Results.Title);

if strlength(plot_title) == 0
    plot_title = sprintf('q versus %s', strrep(x_var, '_', ' '));
end

title(tl, plot_title, 'Interpreter', 'none');

omega_all = T.(omega_var);
omega_all = omega_all(isfinite(omega_all));

if isempty(omega_all)
    clim_global = [];
else
    clim_global = [min(omega_all) max(omega_all)];
end

for gidx = 1:n_groups

    ax = nexttile(tl);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');

    mask = G == gidx;

    x = T.(x_var)(mask);
    y = T.(q_var)(mask);
    om = T.(omega_var)(mask);

    valid = isfinite(x) & isfinite(y) & isfinite(om);

    x = x(valid);
    y = y(valid);
    om = om(valid);

    if isempty(x)
        title(ax, make_group_title(T_groups(gidx, :)), 'Interpreter', 'none');
        hold(ax, 'off');
        continue;
    end

    [om, idx_sort] = sort(om);
    x = x(idx_sort);
    y = y(idx_sort);

    scatter(ax, x, y, 55, om, 'filled', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5);

    plot(ax, x, y, '-', 'LineWidth', 1.2, 'Color', [0.25 0.25 0.25]);

    rho = compute_spearman(x, y);

    xlabel(ax, strrep(x_var, '_', ' '), 'Interpreter', 'none');
    ylabel(ax, strrep(q_var, '_', ' '), 'Interpreter', 'none');

    title(ax, make_group_title(T_groups(gidx, :)), 'Interpreter', 'none');

    ylim(ax, [0 1]);

    if ~isempty(clim_global)
        clim(ax, clim_global);
    end

    text(ax, 0.04, 0.93, sprintf('\\rho_s = %.3f', rho), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 2);

    cb = colorbar(ax);
    cb.Label.String = strrep(omega_var, '_', ' ');
end

base_name = char(p.Results.BaseName);

if isempty(base_name)
    base_name = ['q_vs_feature_controlled_' char(feature_name)];
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

function x_var = resolve_feature_mean_variable(T, feature_name)

vars = string(T.Properties.VariableNames);

candidates = [
    feature_name
    feature_name + "_mean"
];

idx = find(ismember(vars, candidates), 1);

if isempty(idx)
    error('Could not resolve feature variable for: %s', feature_name);
end

x_var = char(vars(idx));

end

function validate_table_variable(T, var_name)

if ~ismember(var_name, T.Properties.VariableNames)
    error('Variable not found in table: %s', var_name);
end

end

function rho = compute_spearman(x, y)

rho = NaN;

if numel(x) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    return;
end

try
    rho = corr(x(:), y(:), 'Type', 'Spearman', 'Rows', 'complete');
catch
    rho = NaN;
end

end

function txt = make_group_title(Tg)

vars = string(Tg.Properties.VariableNames);

parts = strings(1, numel(vars));

for i = 1:numel(vars)
    value_i = Tg.(vars(i));
    value_i = value_i(1);

    if isnumeric(value_i)
        value_txt = sprintf('%.4g', value_i);
    elseif isstring(value_i) || ischar(value_i)
        value_txt = char(string(value_i));
    elseif iscategorical(value_i)
        value_txt = char(string(value_i));
    else
        value_txt = char(string(value_i));
    end

    parts(i) = sprintf('%s = %s', vars(i), value_txt);
end

txt = strjoin(parts, ', ');

end