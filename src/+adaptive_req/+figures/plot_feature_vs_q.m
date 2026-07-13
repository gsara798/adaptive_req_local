function [fig, ax, fit_info] = plot_feature_vs_q(T_raw, feature_name, varargin)
%PLOT_FEATURE_VS_Q Plot q_theory versus one feature.
%
% Usage:
%   adaptive_req.figures.plot_feature_vs_q(T_raw, 'ang_entropy');

p = inputParser;

addRequired(p, 'T_raw', @istable);
addRequired(p, 'feature_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'ColorBy', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ShowLinearFit', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'FigurePosition', [200 200 430 330], @(x) isnumeric(x) && numel(x) == 4);
addParameter(p, 'MarkerSize', 36, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));

parse(p, T_raw, feature_name, varargin{:});

feature_name = char(feature_name);
color_by = char(p.Results.ColorBy);
show_fit = logical(p.Results.ShowLinearFit);

assert(ismember(feature_name, T_raw.Properties.VariableNames), ...
    'Feature not found in table.');

x = T_raw.(feature_name);
y = T_raw.q_theory;

valid = isfinite(x) & isfinite(y);

fig = figure('Color', 'w', 'Position', p.Results.FigurePosition);
ax = axes('Parent', fig);
hold(ax, 'on');

if ~isempty(color_by) && ismember(color_by, T_raw.Properties.VariableNames)
    c = T_raw.(color_by);

    valid = valid & isfinite(c);

    scatter(ax, x(valid), y(valid), ...
        p.Results.MarkerSize, c(valid), ...
        'filled', ...
        'MarkerFaceAlpha', 0.80, ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.3);

    cb = colorbar(ax);
    cb.Label.String = pretty_feature_name(color_by);
    cb.Label.Interpreter = 'tex';
    colormap(ax, parula);
else
    plot(ax, x(valid), y(valid), 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [0.00 0.45 0.74], ...
        'LineWidth', 1.0);
end

fit_info = struct();
fit_info.feature_name = feature_name;
fit_info.r = NaN;
fit_info.rho = NaN;
fit_info.rmse_linear = NaN;

if show_fit && nnz(valid) >= 3 && range(x(valid)) > 0 && range(y(valid)) > 0

    xv = x(valid);
    yv = y(valid);

    pfit = polyfit(xv, yv, 1);

    xgrid = linspace(min(xv), max(xv), 200);
    ygrid = polyval(pfit, xgrid);
    yhat = polyval(pfit, xv);

    plot(ax, xgrid, ygrid, '-', ...
        'LineWidth', 2.0, ...
        'Color', [0.85 0.33 0.10]);

    C = corrcoef(xv, yv);
    r = C(1, 2);
    rho = corr(xv, yv, 'Type', 'Spearman');
    rmse_linear = sqrt(mean((yv - yhat).^2, 'omitnan'));

    fit_info.r = r;
    fit_info.rho = rho;
    fit_info.rmse_linear = rmse_linear;
    fit_info.pfit = pfit;

    if strlength(string(p.Results.Title)) == 0
        title(ax, sprintf('%s | r = %.3f | \\rho = %.3f', ...
            strrep(feature_name, '_', '\_'), r, rho), ...
            'Interpreter', 'tex');
    end
end

xlabel(ax, pretty_feature_name(feature_name), 'Interpreter', 'tex');
ylabel(ax, 'q_{theory}');

if strlength(string(p.Results.Title)) > 0
    title(ax, p.Results.Title, 'Interpreter', 'tex');
end

grid(ax, 'on');
box(ax, 'on');

adaptive_req.figures.apply_style(fig, ax);

end

function label = pretty_feature_name(name)

switch char(name)
    case 'ang_entropy'
        label = 'Angular entropy';
    case 'radial_entropy'
        label = 'Radial entropy';
    case 'q_theory'
        label = 'q_{theory}';
    case 'Omega_sr'
        label = '\Omega (sr)';
    otherwise
        label = strrep(char(name), '_', '\_');
end

end