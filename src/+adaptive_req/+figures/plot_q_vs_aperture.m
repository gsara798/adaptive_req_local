function [fig, ax, T_summary] = plot_q_vs_aperture(T_raw, sweep, varargin)
%PLOT_Q_VS_APERTURE Plot q_theory versus aperture solid angle.
%
% Usage:
%   adaptive_req.figures.plot_q_vs_aperture(T_raw, sweep);
%
% Options:
%   'UseOmega', true
%   'ShowRaw', true
%   'ShowErrorbar', true

p = inputParser;

addRequired(p, 'T_raw', @istable);
addRequired(p, 'sweep', @isstruct);

addParameter(p, 'UseOmega', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowRaw', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowErrorbar', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'FigurePosition', [200 200 420 330], @(x) isnumeric(x) && numel(x) == 4);
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));

parse(p, T_raw, sweep, varargin{:});

use_omega = logical(p.Results.UseOmega);
show_raw = logical(p.Results.ShowRaw);
show_errorbar = logical(p.Results.ShowErrorbar);

T_summary = summarize_by_step(T_raw);

if use_omega
    x_raw = T_raw.Omega_sr;
    x_sum = T_summary.Omega_sr;
    xlab = 'Wavefield solid angle \Omega (sr)';
else
    x_raw = T_raw.aperture_value;
    x_sum = T_summary.aperture_value;

    if isfield(sweep, 'schedule') && isfield(sweep.schedule, 'name')
        xlab = sprintf('%s (%s)', sweep.schedule.name, sweep.schedule.unit);
    else
        xlab = 'Aperture';
    end
end

fig = figure('Color', 'w', 'Position', p.Results.FigurePosition);
ax = axes('Parent', fig);
hold(ax, 'on');

if show_raw
    plot(ax, x_raw, T_raw.q_theory, 'o', ...
        'MarkerSize', 4, ...
        'MarkerEdgeColor', [0.65 0.65 0.65], ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 0.8, ...
        'DisplayName', 'Local patches');
end

if show_errorbar
    errorbar(ax, x_sum, T_summary.q_theory_mean, T_summary.q_theory_std, ...
        'o-', ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', 'w', ...
        'CapSize', 5, ...
        'DisplayName', 'Mean \pm std');
else
    plot(ax, x_sum, T_summary.q_theory_mean, 'o-', ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', 'Mean');
end

xlabel(ax, xlab);
ylabel(ax, 'q_{theory}');

if strlength(string(p.Results.Title)) > 0
    title(ax, p.Results.Title, 'Interpreter', 'tex');
end

grid(ax, 'on');
box(ax, 'on');

if use_omega
    xmin = min(x_sum, [], 'omitnan');
    xmax = max(x_sum, [], 'omitnan');

    if isfinite(xmin) && isfinite(xmax)
        xlim(ax, [xmin xmax]);
    end
end

adaptive_req.figures.apply_style(fig, ax);
legend(ax, 'Location', 'best', 'Box', 'on');

end

function T_summary = summarize_by_step(T_raw)

step_ids = unique(T_raw.step_idx, 'stable');

rows = struct([]);

for i = 1:numel(step_ids)

    mask = T_raw.step_idx == step_ids(i);

    rows(i).step_idx = step_ids(i);
    rows(i).aperture_value = T_raw.aperture_value(find(mask, 1, 'first'));
    rows(i).Omega_sr = T_raw.Omega_sr(find(mask, 1, 'first'));
    rows(i).n_samples = sum(mask);

    qvals = T_raw.q_theory(mask);
    rows(i).q_theory_mean = mean(qvals, 'omitnan');
    rows(i).q_theory_std = std(qvals, 0, 'omitnan');
end

T_summary = struct2table(rows);

end