function [fig, fit_table] = plot_features_vs_q_grid(T_raw, varargin)
%PLOT_FEATURES_VS_Q_GRID Plot q_theory versus multiple features.

p = inputParser;

addRequired(p, 'T_raw', @istable);

addParameter(p, 'FeatureNames', {}, @iscell);
addParameter(p, 'NCols', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'ColorBy', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ShowLinearFit', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'FigurePosition', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));

parse(p, T_raw, varargin{:});

feature_names = p.Results.FeatureNames;

if isempty(feature_names)
    feature_names = default_feature_names(T_raw);
end

nF = numel(feature_names);
ncols = p.Results.NCols;
nrows = ceil(nF / ncols);

if isempty(p.Results.FigurePosition)
    fig_pos = [80 60 430*ncols 310*nrows];
else
    fig_pos = p.Results.FigurePosition;
end

fig = figure('Color', 'w', 'Position', fig_pos);
tl = tiledlayout(fig, nrows, ncols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

rows = struct([]);

for i = 1:nF

    fn = feature_names{i};
    ax = nexttile(tl, i);
    hold(ax, 'on');

    if ~ismember(fn, T_raw.Properties.VariableNames)
        title(ax, sprintf('%s not found', fn), 'Interpreter', 'none');
        axis(ax, 'off');
        continue;
    end

    x = T_raw.(fn);
    y = T_raw.q_theory;

    valid = isfinite(x) & isfinite(y);

    color_by = char(p.Results.ColorBy);

    if ~isempty(color_by) && ismember(color_by, T_raw.Properties.VariableNames)
        c = T_raw.(color_by);
        valid = valid & isfinite(c);

        scatter(ax, x(valid), y(valid), 22, c(valid), ...
            'filled', ...
            'MarkerFaceAlpha', 0.75, ...
            'MarkerEdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 0.25);
        colormap(ax, parula);
    else
        plot(ax, x(valid), y(valid), 'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', [0.00 0.45 0.74], ...
            'LineWidth', 0.9);
    end

    r = NaN;
    rho = NaN;
    rmse_linear = NaN;

    if logical(p.Results.ShowLinearFit) && nnz(valid) >= 3 && ...
            range(x(valid)) > 0 && range(y(valid)) > 0

        xv = x(valid);
        yv = y(valid);

        pfit = polyfit(xv, yv, 1);
        xgrid = linspace(min(xv), max(xv), 200);
        ygrid = polyval(pfit, xgrid);
        yhat = polyval(pfit, xv);

        plot(ax, xgrid, ygrid, '-', ...
            'LineWidth', 1.8, ...
            'Color', [0.85 0.33 0.10]);

        C = corrcoef(xv, yv);
        r = C(1, 2);
        rho = corr(xv, yv, 'Type', 'Spearman');
        rmse_linear = sqrt(mean((yv - yhat).^2, 'omitnan'));
    end

    xlabel(ax, pretty_feature_name(fn), 'Interpreter', 'tex');
    ylabel(ax, 'q_{theory}');

    title(ax, sprintf('%s | r=%.2f | \\rho=%.2f', ...
        strrep(fn, '_', '\_'), r, rho), ...
        'Interpreter', 'tex');

    grid(ax, 'on');
    box(ax, 'on');

    adaptive_req.figures.apply_style(fig, ax, 'FontSize', 12);

    rows(i).feature = string(fn);
    rows(i).r = r;
    rows(i).rho = rho;
    rows(i).rmse_linear = rmse_linear;
    rows(i).n = nnz(valid);
end

fit_table = struct2table(rows);

if ~isempty(char(p.Results.ColorBy))
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = pretty_feature_name(char(p.Results.ColorBy));
    cb.Label.Interpreter = 'tex';
end

end

function feature_names = default_feature_names(T_raw)

candidate = { ...
    'radial_entropy', ...
    'width_75_25_rel', ...
    'width_90_50_rel', ...
    'width_90_10_rel', ...
    'lowk_frac_rel', ...
    'midband_frac_rel', ...
    'highk_frac_rel', ...
    'ang_entropy', ...
    'circ_var', ...
    'dom_dir_frac', ...
    'window_max_frac', ...
    'window_cf'};

feature_names = {};

for i = 1:numel(candidate)
    if ismember(candidate{i}, T_raw.Properties.VariableNames)
        feature_names{end+1} = candidate{i}; %#ok<AGROW>
    end
end

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