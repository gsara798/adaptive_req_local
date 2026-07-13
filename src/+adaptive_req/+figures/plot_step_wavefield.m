function [fig, ax] = plot_step_wavefield(sweep, step_idx, varargin)
%PLOT_STEP_WAVEFIELD Plot one stored simulated wavefield.
%
% Usage:
%   adaptive_req.figures.plot_step_wavefield(sweep, 3);
%
% Requires:
%   StoreWavefields = true when running the sweep.

p = inputParser;

addRequired(p, 'sweep', @isstruct);
addRequired(p, 'step_idx', @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'Realization', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Component', 'real', @(x) ischar(x) || isstring(x));
addParameter(p, 'Normalize', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowPatches', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'NumPatches', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'FigurePosition', [200 200 420 330], @(x) isnumeric(x) && numel(x) == 4);

parse(p, sweep, step_idx, varargin{:});

r = round(p.Results.Realization);

if ~isfield(sweep, 'sims') || isempty(sweep.sims)
    error('sweep.sims is empty. Run the sweep with StoreWavefields = true.');
end

sim = sweep.sims{step_idx, r};

cfg_i = sim.cfg;

if isfield(sweep, 'feat_cfg_base')
    feat_cfg = sweep.feat_cfg_base;
else
    error('sweep.feat_cfg_base was not found.');
end

[~, feat_cfg_i] = adaptive_req.config.default_req_config(cfg_i, feat_cfg);

if isempty(p.Results.NumPatches)
    if isfield(sweep, 'num_patches')
        num_patches = sweep.num_patches;
    else
        num_patches = 1;
    end
else
    num_patches = p.Results.NumPatches;
end

patch_pack = adaptive_req.simulate.build_patch_windows( ...
    cfg_i, feat_cfg_i, ...
    'NumPatches', num_patches);

Uplot = select_component(sim.Uxz, char(p.Results.Component));

if logical(p.Results.Normalize)
    Uplot = Uplot ./ max(abs(Uplot(:)) + eps);
end

fig = figure('Color', 'w', 'Position', p.Results.FigurePosition);
ax = axes('Parent', fig);

imagesc(ax, sim.x * 100, sim.z * 100, Uplot);
set(ax, 'YDir', 'normal');
axis(ax, 'image');

colormap(ax, jet);

if logical(p.Results.Normalize)
    clim(ax, [-1 1]);
end

cb = colorbar(ax);
cb.Label.String = 'Normalized Re\{u(x,z)\}';
cb.Label.Interpreter = 'tex';
cb.Ticks = [-1 0 1];

xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');

if isfield(sweep, 'schedule')
    ttl = sprintf('Step %d: \\Omega = %.2f sr', ...
        step_idx, sweep.schedule.Omega_sr(step_idx));
else
    ttl = sprintf('Step %d', step_idx);
end

title(ax, ttl, 'Interpreter', 'tex');

hold(ax, 'on');

if logical(p.Results.ShowPatches)
    for pidx = 1:patch_pack.n_patches

        xi = patch_pack.x_idx_list{pidx};
        zi = patch_pack.z_idx_list{pidx};

        x0 = sim.x(min(xi)) * 100;
        z0 = sim.z(min(zi)) * 100;

        w = (sim.x(max(xi)) - sim.x(min(xi))) * 100;
        h = (sim.z(max(zi)) - sim.z(min(zi))) * 100;

        rectangle(ax, 'Position', [x0 z0 w h], ...
            'EdgeColor', 'w', ...
            'LineWidth', 1.1);

        xc = sim.x(patch_pack.cx_list(pidx)) * 100;
        zc = sim.z(patch_pack.cz_list(pidx)) * 100;

        text(ax, xc, zc, sprintf('%d', pidx), ...
            'Color', 'w', ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
end

hold(ax, 'off');

adaptive_req.figures.apply_style(fig, ax);
adaptive_req.figures.format_axis_ticks(ax, 1);
adaptive_req.figures.format_colorbar_ticks(cb, 1);

end

function Uplot = select_component(U, component)

switch lower(component)
    case 'real'
        Uplot = real(U);
    case 'imag'
        Uplot = imag(U);
    case 'abs'
        Uplot = abs(U);
    case 'phase'
        Uplot = angle(U);
    otherwise
        error('Unknown component: %s', component);
end

end