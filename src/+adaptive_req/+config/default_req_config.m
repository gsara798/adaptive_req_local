function [req_cfg, feat_cfg] = default_req_config(cfg, feat_cfg, varargin)
%DEFAULT_REQ_CONFIG Build REQ configuration derived from cfg and feat_cfg.
%
% Usage:
%   req_cfg = adaptive_req.config.default_req_config(cfg, feat_cfg);
%
%   [req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
%       cfg, feat_cfg, ...
%       'Nbins', 'auto', ...
%       'Nbins_auto_oversample', 1, ...
%       'Nbins_min', 16, ...
%       'smooth_sigma', 1.0);
%
% This function defines and derives all REQ-related parameters:
%   local window size
%   circular Hann window
%   zero-padding
%   radial binning
%   smoothing
%   true wavenumber
%   optional donut filter settings

validate_required_fields(cfg, {'f0', 'cs_bg', 'dx', 'dz'});
validate_required_fields(feat_cfg, {'M', 'cs_guess', 'gamma_win', 'pad_factor'});

p = inputParser;
p.FunctionName = 'adaptive_req.config.default_req_config';

% -------------------------------------------------------------------------
% Optional overrides for parameters that usually come from feat_cfg
% -------------------------------------------------------------------------

addParameter(p, 'M', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x > 0));

addParameter(p, 'cs_guess', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x > 0));

addParameter(p, 'gamma_win', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x > 0));

addParameter(p, 'pad_factor', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x >= 0));

% -------------------------------------------------------------------------
% REQ spectral parameters
% -------------------------------------------------------------------------

addParameter(p, 'smooth_sigma', 1.0, @(x) ...
    isnumeric(x) && isscalar(x) && x >= 0);

addParameter(p, 'Nbins', 'auto', @(x) ...
    is_valid_nbins(x));

addParameter(p, 'Nbins_auto_oversample', 1, @(x) ...
    isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'Nbins_min', 16, @(x) ...
    isnumeric(x) && isscalar(x) && x >= 1);

% Legacy alias. Prefer Nbins_auto_oversample in new scripts.
addParameter(p, 'oversample_factor', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x > 0));

% -------------------------------------------------------------------------
% Window construction options
% -------------------------------------------------------------------------

addParameter(p, 'cs_guess_mode', 'fixed', @(x) ...
    ischar(x) || isstring(x));

addParameter(p, 'window_mode', 'current', @(x) ...
    ischar(x) || isstring(x));

addParameter(p, 'fixed_win_size', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && x >= 3));

% -------------------------------------------------------------------------
% Optional donut filter
% -------------------------------------------------------------------------

addParameter(p, 'use_donut', false, @(x) ...
    islogical(x) || isnumeric(x));

addParameter(p, 'donut_cs_min', 1.0, @(x) ...
    isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'donut_cs_max', 5.0, @(x) ...
    isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'donut_taper_rel', 0.06, @(x) ...
    isnumeric(x) && isscalar(x) && x >= 0);

addParameter(p, 'apply_donut_to_final_map', false, @(x) ...
    islogical(x) || isnumeric(x));

parse(p, varargin{:});

opt = p.Results;

opt.cs_guess_mode = lower(char(opt.cs_guess_mode));
opt.window_mode = lower(char(opt.window_mode));

if ~isempty(opt.oversample_factor)
    opt.Nbins_auto_oversample = opt.oversample_factor;
end

if opt.donut_cs_min >= opt.donut_cs_max
    error('donut_cs_min must be smaller than donut_cs_max.');
end

% -------------------------------------------------------------------------
% Apply optional overrides to feat_cfg
% -------------------------------------------------------------------------

if ~isempty(opt.M)
    feat_cfg.M = opt.M;
end

if ~isempty(opt.cs_guess)
    feat_cfg.cs_guess = opt.cs_guess;
end

if ~isempty(opt.gamma_win)
    feat_cfg.gamma_win = opt.gamma_win;
end

if ~isempty(opt.pad_factor)
    feat_cfg.pad_factor = opt.pad_factor;
end

% -------------------------------------------------------------------------
% Select cs_guess used for the local window
% -------------------------------------------------------------------------

switch opt.cs_guess_mode

    case 'fixed'
        cs_guess_used = feat_cfg.cs_guess;

    case 'match_true'
        cs_guess_used = cfg.cs_bg;

    otherwise
        error('Unknown cs_guess_mode: %s', opt.cs_guess_mode);
end

% -------------------------------------------------------------------------
% Compute local window size
% -------------------------------------------------------------------------

switch opt.window_mode

    case 'current'
        lambda_guess = cs_guess_used / cfg.f0;
        win_size = round(feat_cfg.M * lambda_guess / cfg.dx);

    case 'fixed_pixels'
        if isempty(opt.fixed_win_size)
            error('fixed_win_size must be provided when window_mode is fixed_pixels.');
        end

        win_size = round(opt.fixed_win_size);
        lambda_guess = win_size * cfg.dx / feat_cfg.M;

    otherwise
        error('Unknown window_mode: %s', opt.window_mode);
end

if win_size < 3
    win_size = 3;
end

if mod(win_size, 2) == 0
    win_size = win_size + 1;
end

half_win = floor(win_size / 2);

lambda_true = cfg.cs_bg / cfg.f0;
M_eff = (win_size * cfg.dx) / lambda_true;

% -------------------------------------------------------------------------
% Update feature configuration with derived local-window quantities
% -------------------------------------------------------------------------

feat_cfg.win_size = win_size;
feat_cfg.half_win = half_win;
feat_cfg.cs_guess_used = cs_guess_used;
feat_cfg.lambda_guess_used = lambda_guess;
feat_cfg.lambda_true = lambda_true;
feat_cfg.M_eff = M_eff;

% -------------------------------------------------------------------------
% Build REQ configuration
% -------------------------------------------------------------------------

req_cfg = struct();

req_cfg.M = feat_cfg.M;
req_cfg.gamma_win = feat_cfg.gamma_win;
req_cfg.pad_factor = feat_cfg.pad_factor;

req_cfg.win_size = win_size;
req_cfg.half_win = half_win;

req_cfg.cs_guess_used = cs_guess_used;
req_cfg.lambda_guess_used = lambda_guess;
req_cfg.lambda_true = lambda_true;
req_cfg.M_eff = M_eff;

req_cfg.k0_true = 2*pi*cfg.f0 / cfg.cs_bg;
req_cfg.cs_true = cfg.cs_bg;
req_cfg.f0 = cfg.f0;

req_cfg.smooth_sigma = opt.smooth_sigma;

req_cfg.Nbins = normalize_nbins(opt.Nbins);
req_cfg.Nbins_auto_oversample = opt.Nbins_auto_oversample;
req_cfg.Nbins_min = round(opt.Nbins_min);

req_cfg.use_donut = logical(opt.use_donut);
req_cfg.donut_cs_min = opt.donut_cs_min;
req_cfg.donut_cs_max = opt.donut_cs_max;
req_cfg.donut_taper_rel = opt.donut_taper_rel;
req_cfg.apply_donut_to_final_map = logical(opt.apply_donut_to_final_map);

req_cfg.W2 = circular_hann_window( ...
    req_cfg.win_size, ...
    req_cfg.win_size, ...
    req_cfg.gamma_win, ...
    cfg.dx, ...
    cfg.dz);

req_cfg.PAD = max(0, round([ ...
    req_cfg.pad_factor * req_cfg.win_size, ...
    req_cfg.pad_factor * req_cfg.win_size]));

end

% =========================================================================
% Local helper functions
% =========================================================================

function W = circular_hann_window(nx, nz, gamma, dx, dz)

x = ((1:nx) - (nx + 1)/2) * dx;
z = ((1:nz) - (nz + 1)/2) * dz;

[X, Z] = meshgrid(x, z);

r = sqrt(X.^2 + Z.^2);
R = gamma * min(max(abs(x)), max(abs(z)));

W = zeros(nz, nx);

if R <= 0
    return;
end

mask = r <= R;
W(mask) = 0.5 * (1 + cos(pi * r(mask) / R));

end

function validate_required_fields(s, names)

for i = 1:numel(names)
    if ~isfield(s, names{i})
        error('Missing required field: %s', names{i});
    end
end

end

function tf = is_valid_nbins(x)

tf = is_auto_nbins(x) || ...
     (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);

end

function tf = is_auto_nbins(x)

tf = (ischar(x) && strcmpi(x, 'auto')) || ...
     (isstring(x) && isscalar(x) && strcmpi(x, "auto"));

end

function xout = normalize_nbins(x)

if is_auto_nbins(x)
    xout = 'auto';
else
    xout = round(double(x));
end

end