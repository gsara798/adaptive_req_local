function figs = plot_step_diagnostics(sim, patch_pack, rep, varargin)
%PLOT_STEP_DIAGNOSTICS Plot wavefield, 2D spectrum, and radial REQ diagnostics.
%
% Usage:
%   figs = adaptive_req.figures.plot_step_diagnostics(sim, patch_pack, rep);
%
% This function is intended to be called from:
%
%   adaptive_req.studies.run_aperture_sweep
%
% Inputs:
%   sim:
%       Output from adaptive_req.simulate.run_single_simulation.
%
%   patch_pack:
%       Output from adaptive_req.simulate.build_patch_windows.
%
%   rep:
%       Representative REQ diagnostic structure, typically req_curve from
%       adaptive_req.quantile.compute_quantile_from_patch.

p = inputParser;

addRequired(p, 'sim', @isstruct);
addRequired(p, 'patch_pack', @isstruct);
addRequired(p, 'rep', @isstruct);

addParameter(p, 'StepIndex', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Realization', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'PatchIndex', 1, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'OmegaSr', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'StepValue', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'StepName', 'Aperture', @(x) ischar(x) || isstring(x));
addParameter(p, 'StepUnit', '', @(x) ischar(x) || isstring(x));

addParameter(p, 'QMean', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'QValue', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'KReal', NaN, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'ShowPatchCenters', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowSelectedPatchBox', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowPatchLabels', false, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'FieldComponent', 'real', @(x) ischar(x) || isstring(x));
addParameter(p, 'NormalizeField', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'KShow', 1.35, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'TickDecimals', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Visible', true, @(x) islogical(x) || isnumeric(x));

parse(p, sim, patch_pack, rep, varargin{:});

opt = p.Results;

patch_idx = round(opt.PatchIndex);
nTickDecimals = round(opt.TickDecimals);

if patch_idx < 1 || patch_idx > patch_pack.n_patches
    error('PatchIndex must be between 1 and %d.', patch_pack.n_patches);
end

visible_state = logical_to_visible(opt.Visible);

x = sim.x;
z = sim.z;
Uxz = sim.Uxz;

if isfield(sim, 'cfg') && isfield(sim.cfg, 'f0') && isfield(sim.cfg, 'cs_bg')
    k0_default = 2*pi*sim.cfg.f0/sim.cfg.cs_bg;
else
    k0_default = NaN;
end

k_real = choose_first_finite( ...
    opt.KReal, ...
    getfield_with_default(rep, 'k_real', NaN), ...
    k0_default);

q_val = choose_first_finite( ...
    opt.QValue, ...
    getfield_with_default(rep, 'q', NaN));

figs = struct();

%% =========================================================================
% Figure 1: wavefield
% =========================================================================

figs.field = figure( ...
    'Name', 'REQ Diagnostics Field', ...
    'Color', 'w', ...
    'Visible', visible_state, ...
    'Position', [200 200 390 300]);

ax1 = axes('Parent', figs.field);

field_plot = select_field_component(Uxz, char(opt.FieldComponent));

if logical(opt.NormalizeField)
    field_plot = field_plot ./ max(abs(field_plot(:)) + eps);
end

imagesc(ax1, x*1e2, z*1e2, field_plot);
set(ax1, 'YDir', 'normal');
axis(ax1, 'image');
colormap(ax1, jet);

if logical(opt.NormalizeField) && ~strcmpi(opt.FieldComponent, 'abs')
    caxis(ax1, [-1 1]);
end

cb1 = colorbar(ax1);
cb1.Label.String = make_field_colorbar_label(char(opt.FieldComponent), logical(opt.NormalizeField));
cb1.Label.Interpreter = 'tex';

if logical(opt.NormalizeField) && ~strcmpi(opt.FieldComponent, 'abs')
    cb1.Ticks = [-1 0 1];
end

hold(ax1, 'on');

if logical(opt.ShowPatchCenters)
    for pidx = 1:patch_pack.n_patches
        plot(ax1, ...
            x(patch_pack.cx_list(pidx))*1e2, ...
            z(patch_pack.cz_list(pidx))*1e2, ...
            'wo', ...
            'MarkerSize', 5, ...
            'LineWidth', 1.0);

        if logical(opt.ShowPatchLabels)
            text(ax1, ...
                x(patch_pack.cx_list(pidx))*1e2, ...
                z(patch_pack.cz_list(pidx))*1e2, ...
                sprintf('%d', pidx), ...
                'Color', 'w', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end
end

if logical(opt.ShowSelectedPatchBox)
    xi = patch_pack.x_idx_list{patch_idx};
    zi = patch_pack.z_idx_list{patch_idx};

    x0 = x(min(xi))*1e2;
    z0 = z(min(zi))*1e2;
    w = (x(max(xi)) - x(min(xi)))*1e2;
    h = (z(max(zi)) - z(min(zi)))*1e2;

    rectangle(ax1, ...
        'Position', [x0 z0 w h], ...
        'EdgeColor', 'w', ...
        'LineWidth', 1.3);
end

hold(ax1, 'off');

xlabel(ax1, 'x (cm)');
ylabel(ax1, 'z (cm)');

title(ax1, make_step_title(opt, patch_idx), 'Interpreter', 'tex');

adaptive_req.figures.apply_style(figs.field, ax1);
adaptive_req.figures.format_axis_ticks(ax1, nTickDecimals);

if logical(opt.NormalizeField) && ~strcmpi(opt.FieldComponent, 'abs')
    cb1.Ticks = [-1 0 1];
end

adaptive_req.figures.format_colorbar_ticks(cb1, nTickDecimals);

%% =========================================================================
% Figure 2: normalized 2D spectrum
% =========================================================================

figs.spectrum = figure( ...
    'Name', 'REQ Diagnostics Spectrum', ...
    'Color', 'w', ...
    'Visible', visible_state, ...
    'Position', [200 200 390 300]);

ax2 = axes('Parent', figs.spectrum);

has_spectrum = all(isfield(rep, {'kx', 'kz'})) && ...
               (isfield(rep, 'Ssm_norm') || isfield(rep, 'Ssm')) && ...
               isfinite(k_real) && k_real > 0;

if has_spectrum

    if isfield(rep, 'Ssm_norm')
        Ssm_norm = rep.Ssm_norm;
    else
        Ssm_norm = rep.Ssm ./ max(rep.Ssm(:) + eps);
    end

    kx_n = rep.kx ./ k_real;
    kz_n = rep.kz ./ k_real;

    imagesc(ax2, kx_n, kz_n, Ssm_norm);
    set(ax2, 'YDir', 'normal');
    axis(ax2, 'image');
    colormap(ax2, turbo);
    caxis(ax2, [0 1]);

    hold(ax2, 'on');

    th = linspace(0, 2*pi, 500);
    plot(ax2, cos(th), sin(th), '--w', 'LineWidth', 1.5);
    plot(ax2, 0, 0, 'w+', 'MarkerSize', 7, 'LineWidth', 1.0);

    ang = pi/6;
    quiver(ax2, 0, 0, 0.92*cos(ang), 0.92*sin(ang), 0, ...
        'Color', 'w', ...
        'LineWidth', 1.5, ...
        'MaxHeadSize', 0.35);

    text(ax2, 1.05*cos(ang), 1.05*sin(ang), 'k_0', ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom');

    kshow = opt.KShow;
    xlim(ax2, [-kshow kshow]);
    ylim(ax2, [-kshow kshow]);

    hold(ax2, 'off');

    xlabel(ax2, 'k_x / k_0');
    ylabel(ax2, 'k_z / k_0');

    cb2 = colorbar(ax2);
    cb2.Label.String = 'Normalized power';
    cb2.Ticks = [0 0.2 0.4 0.6 0.8 1];

    adaptive_req.figures.apply_style(figs.spectrum, ax2);
    adaptive_req.figures.format_axis_ticks(ax2, nTickDecimals);
    adaptive_req.figures.format_colorbar_ticks(cb2, nTickDecimals);

    cb2.Ticks = [0 0.2 0.4 0.6 0.8 1];

else
    axis(ax2, 'off');
    text(ax2, 0.5, 0.5, ...
        {'2D spectrum diagnostics not available.', ...
         'req_curve should contain kx, kz, and Ssm or Ssm_norm.'}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 12);
end

%% =========================================================================
% Figure 3: radial spectrum and cumulative energy
% =========================================================================

figs.radial = figure( ...
    'Name', 'REQ Diagnostics Radial', ...
    'Color', 'w', ...
    'Visible', visible_state, ...
    'Position', [200 200 360 300]);

ax3 = axes('Parent', figs.radial);

k_cent = getfield_with_default(rep, 'k_cent', []);
Ecum = getfield_with_default(rep, 'Ecum', []);
Srad_norm = get_radial_power_norm(rep);

has_radial = ~isempty(k_cent) && ~isempty(Ecum) && ~isempty(Srad_norm) && ...
             isfinite(k_real) && k_real > 0;

if has_radial

    k_cent_n = k_cent(:) ./ k_real;
    Srad_norm = Srad_norm(:);
    Ecum = Ecum(:);

    if ~isfinite(q_val)
        q_val = interp1(k_cent_n, Ecum, 1, 'linear', 'extrap');
    end

    hold(ax3, 'on');

    yyaxis(ax3, 'left');

    h1 = plot(ax3, k_cent_n, Srad_norm, '-', ...
        'LineWidth', 1.8, ...
        'Color', [0.15 0.25 0.95]);

    ylabel(ax3, 'Normalized $\bar{S}_{\mathrm{rad}}$', ...
        'Interpreter', 'latex');

    ylim(ax3, [0 1.05]);

    yyaxis(ax3, 'right');

    h2 = plot(ax3, k_cent_n, Ecum, '-', ...
        'LineWidth', 1.8, ...
        'Color', [0 0 0]);

    h3 = plot(ax3, 1, q_val, 'o', ...
        'MarkerSize', 6, ...
        'LineWidth', 1.2, ...
        'MarkerFaceColor', [1 1 1], ...
        'Color', [1 0.25 0.25]);

    ylabel(ax3, '$E_{\mathrm{cum}}(k)$', ...
        'Interpreter', 'latex');

    ylim(ax3, [0 1]);

    idx_last = find(Srad_norm > 0.01, 1, 'last');

    if isempty(idx_last)
        k_last_n = max(k_cent_n);
    else
        k_last_n = k_cent_n(idx_last);
    end

    xmax_plot = min(max(k_cent_n), 1.15*max(1, k_last_n));
    xlim(ax3, [0 xmax_plot]);

    xlabel(ax3, 'k / k_0');
    grid(ax3, 'on');

    if isfinite(opt.QMean)
        title(ax3, sprintf('q = %.3f, mean q = %.3f', q_val, opt.QMean), ...
            'Interpreter', 'tex');
    else
        title(ax3, sprintf('q = %.3f', q_val), ...
            'Interpreter', 'tex');
    end

    hold(ax3, 'off');

    adaptive_req.figures.apply_style(figs.radial, ax3);

    legend(ax3, [h1 h2 h3], ...
        {'$\bar{S}_{rad}(k)$', '$E(k)$', 'selected $q$'}, ...
        'Location', 'southeast', ...
        'Box', 'on', ...
        'Interpreter', 'latex', ...
        'FontSize', 12);

    yyaxis(ax3, 'left');
    ax3.YColor = [0.15 0.25 0.95];

    yyaxis(ax3, 'right');
    ax3.YColor = [0 0 0];

    adaptive_req.figures.format_axis_ticks(ax3, nTickDecimals);

else
    axis(ax3, 'off');
    text(ax3, 0.5, 0.5, ...
        {'Radial diagnostics not available.', ...
         'req_curve should contain k_cent, Srad, and Ecum.'}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 12);
end

end

%% =========================================================================
% Local helper functions
% =========================================================================

function Uplot = select_field_component(U, component)

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
        error('Unknown field component: %s', component);
end

end

function label = make_field_colorbar_label(component, normalize_field)

switch lower(component)
    case 'real'
        base = 'Re\{u(x,z)\}';
    case 'imag'
        base = 'Im\{u(x,z)\}';
    case 'abs'
        base = '|u(x,z)|';
    case 'phase'
        base = 'phase(u)';
    otherwise
        base = 'u(x,z)';
end

if normalize_field
    label = ['Normalized ' base];
else
    label = base;
end

end

function title_txt = make_step_title(opt, patch_idx)

parts = strings(0);

if isfinite(opt.StepIndex)
    parts(end+1) = sprintf('Step %d', opt.StepIndex);
end

if isfinite(opt.Realization)
    parts(end+1) = sprintf('realization %d', opt.Realization);
end

parts(end+1) = sprintf('patch %d', patch_idx);

if isfinite(opt.OmegaSr)
    parts(end+1) = sprintf('\\Omega = %.2f sr', opt.OmegaSr);
elseif isfinite(opt.StepValue)
    parts(end+1) = sprintf('%s = %.2f %s', ...
        char(opt.StepName), opt.StepValue, char(opt.StepUnit));
end

title_txt = strjoin(parts, ', ');

end

function Srad_norm = get_radial_power_norm(rep)

if isfield(rep, 'Srad_norm')
    Srad_norm = rep.Srad_norm;
elseif isfield(rep, 'Srad')
    Srad_norm = rep.Srad ./ max(rep.Srad(:) + eps);
else
    Srad_norm = [];
end

end

function val = getfield_with_default(S, fieldName, defaultVal)

if isstruct(S) && isfield(S, fieldName)
    val = S.(fieldName);
else
    val = defaultVal;
end

end

function val = choose_first_finite(varargin)

val = NaN;

for i = 1:nargin
    candidate = varargin{i};
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        val = candidate;
        return;
    end
end

end

function visible_state = logical_to_visible(tf)

if logical(tf)
    visible_state = 'on';
else
    visible_state = 'off';
end

end