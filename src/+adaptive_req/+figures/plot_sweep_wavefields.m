function fig = plot_sweep_wavefields(sweep, varargin)
%PLOT_SWEEP_WAVEFIELDS Plot one stored wavefield per aperture step.

p = inputParser;

addRequired(p, 'sweep', @isstruct);

addParameter(p, 'Realization', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Component', 'real', @(x) ischar(x) || isstring(x));
addParameter(p, 'Normalize', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'StepIndices', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'FigurePosition', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 4));

parse(p, sweep, varargin{:});

if ~isfield(sweep, 'sims') || isempty(sweep.sims)
    error('sweep.sims is empty. Run the sweep with StoreWavefields = true.');
end

r = round(p.Results.Realization);

if isempty(p.Results.StepIndices)
    step_indices = sweep.step_indices;
else
    step_indices = p.Results.StepIndices(:).';
end

n_steps_plot = numel(step_indices);

if isempty(p.Results.FigurePosition)
    fig_pos = [100 100 280*n_steps_plot 330];
else
    fig_pos = p.Results.FigurePosition;
end

fig = figure('Color', 'w', 'Position', fig_pos);

tl = tiledlayout(fig, 1, n_steps_plot, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for ii = 1:n_steps_plot

    step_idx = step_indices(ii);
    sim = sweep.sims{step_idx, r};

    ax = nexttile(tl);
    Uplot = select_component(sim.Uxz, char(p.Results.Component));

    if logical(p.Results.Normalize)
        Uplot = Uplot ./ max(abs(Uplot(:)) + eps);
    end

    imagesc(ax, sim.x * 100, sim.z * 100, Uplot);
    set(ax, 'YDir', 'normal');
    axis(ax, 'image');

    colormap(ax, jet);

    if logical(p.Results.Normalize)
        clim(ax, [-1 1]);
    end

    cb = colorbar(ax);
    cb.Ticks = [-1 0 1];

    xlabel(ax, 'x (cm)');
    ylabel(ax, 'z (cm)');

    if isfield(sweep, 'schedule')
        title(ax, sprintf('\\Omega = %.2f sr', sweep.schedule.Omega_sr(step_idx)), ...
            'Interpreter', 'tex');
    else
        title(ax, sprintf('Step %d', step_idx));
    end

    adaptive_req.figures.apply_style(fig, ax, 'FontSize', 12);
    adaptive_req.figures.format_axis_ticks(ax, 1);
    adaptive_req.figures.format_colorbar_ticks(cb, 1);
end

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