function [fig, ax] = plot_feature_space(T_raw, varargin)
%PLOT_FEATURE_SPACE Plot two features colored by q_theory.
%
% Default:
%   x = ang_entropy
%   y = radial_entropy
%   color = q_theory

p = inputParser;

addRequired(p, 'T_raw', @istable);

addParameter(p, 'XFeature', 'ang_entropy', @(x) ischar(x) || isstring(x));
addParameter(p, 'YFeature', 'radial_entropy', @(x) ischar(x) || isstring(x));
addParameter(p, 'ColorVariable', 'q_theory', @(x) ischar(x) || isstring(x));
addParameter(p, 'FigurePosition', [200 200 430 360], @(x) isnumeric(x) && numel(x) == 4);
addParameter(p, 'MarkerSize', 42, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));

parse(p, T_raw, varargin{:});

xvar = char(p.Results.XFeature);
yvar = char(p.Results.YFeature);
cvar = char(p.Results.ColorVariable);

assert(ismember(xvar, T_raw.Properties.VariableNames), 'XFeature not found in table.');
assert(ismember(yvar, T_raw.Properties.VariableNames), 'YFeature not found in table.');
assert(ismember(cvar, T_raw.Properties.VariableNames), 'ColorVariable not found in table.');

valid = isfinite(T_raw.(xvar)) & isfinite(T_raw.(yvar)) & isfinite(T_raw.(cvar));

fig = figure('Color', 'w', 'Position', p.Results.FigurePosition);
ax = axes('Parent', fig);
hold(ax, 'on');

scatter(ax, T_raw.(xvar)(valid), T_raw.(yvar)(valid), ...
    p.Results.MarkerSize, T_raw.(cvar)(valid), ...
    'filled', ...
    'MarkerFaceAlpha', 0.90, ...
    'MarkerEdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 0.35);

xlabel(ax, pretty_feature_name(xvar), 'Interpreter', 'tex');
ylabel(ax, pretty_feature_name(yvar), 'Interpreter', 'tex');

if strlength(string(p.Results.Title)) > 0
    title(ax, p.Results.Title, 'Interpreter', 'tex');
end

grid(ax, 'on');
box(ax, 'on');

cb = colorbar(ax);
cb.Label.String = pretty_feature_name(cvar);
cb.Label.Interpreter = 'tex';

colormap(ax, parula);

adaptive_req.figures.apply_style(fig, ax);
adaptive_req.figures.format_colorbar_ticks(cb, 2);

end

function label = pretty_feature_name(name)

switch char(name)
    case 'ang_entropy'
        label = 'Angular entropy';
    case 'radial_entropy'
        label = 'Radial entropy';
    case 'q_theory'
        label = 'Reference q^*';
    case 'Omega_sr'
        label = '\Omega (sr)';
    otherwise
        label = strrep(char(name), '_', '\_');
end

end