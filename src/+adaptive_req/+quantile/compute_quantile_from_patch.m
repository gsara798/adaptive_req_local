function [q_val, req_curve] = compute_quantile_from_patch(V_raw, cfg, req)
%COMPUTE_QUANTILE_FROM_PATCH Compute REQ reference quantile from one patch.
%
% Output fields include radial and 2D spectral diagnostics:
%   req_curve.kx
%   req_curve.kz
%   req_curve.Sxz
%   req_curve.Ssm
%   req_curve.Ssm_norm
%   req_curve.Srad
%   req_curve.Srad_norm
%   req_curve.Ecum
%   req_curve.k_cent
%   req_curve.q
%   req_curve.k_real
%   req_curve.Nbins_requested
%   req_curve.Nbins_effective

validate_required_fields(cfg, {'dx', 'dz', 'f0', 'cs_bg'});
validate_required_fields(req, { ...
    'W2', ...
    'PAD', ...
    'k0_true', ...
    'Nbins', ...
    'Nbins_auto_oversample', ...
    'Nbins_min', ...
    'smooth_sigma'});

if getfield_with_default(req, 'use_donut', false)
    validate_required_fields(req, { ...
        'donut_cs_min', ...
        'donut_cs_max', ...
        'donut_taper_rel'});
end

if ~isequal(size(V_raw), size(req.W2))
    error(['V_raw and req.W2 must have the same size. ', ...
           'V_raw is %s and W2 is %s.'], ...
           mat2str(size(V_raw)), mat2str(size(req.W2)));
end

% -------------------------------------------------------------------------
% Window, demean, and normalize patch
% -------------------------------------------------------------------------

V_win = V_raw .* req.W2;

finite_mask = isfinite(V_win);

if ~any(finite_mask(:))
    q_val = NaN;
    req_curve = empty_req_curve(req, cfg);
    return;
end

mu = mean(V_win(finite_mask));
sigma = std(V_win(finite_mask), 0);

V_win(~finite_mask) = mu;
V_win = V_win - mu;
V_win = V_win ./ (sigma + eps);

% -------------------------------------------------------------------------
% Zero pad and spectral geometry
% -------------------------------------------------------------------------

pre = get_precomputed_geometry(req);

if isempty(pre)
    V_pad = padarray(V_win, req.PAD, 0, 'both');

    [Nz2, Nx2] = size(V_pad);

    dkx = 2*pi / (Nx2 * cfg.dx);
    dkz = 2*pi / (Nz2 * cfg.dz);

    kx = (-floor(Nx2/2):ceil(Nx2/2)-1) * dkx;
    kz = (-floor(Nz2/2):ceil(Nz2/2)-1) * dkz;

    [KX, KZ] = meshgrid(kx, kz);
    KR = sqrt(KX.^2 + KZ.^2);
    TH = mod(atan2(KZ, KX), 2*pi);

    kmax = min(max(abs(kx)), max(abs(kz)));

    Nbins_requested = req.Nbins;

    if is_auto_nbins(Nbins_requested)
        dk_DFT = min(dkx, dkz);
        oversample_factor = req.Nbins_auto_oversample;

        Nbins_effective = round(kmax / (dk_DFT / oversample_factor));
        Nbins_effective = max(req.Nbins_min, Nbins_effective);
    else
        Nbins_effective = round(double(Nbins_requested));
    end

    Nbins_effective = max(1, Nbins_effective);

    k_edges = linspace(0, kmax, Nbins_effective + 1);
    k_cent = 0.5 * (k_edges(1:end-1) + k_edges(2:end));
    radial_bin = discretize(KR(:), k_edges);
else
    V_pad = zeros(pre.Nz2, pre.Nx2, 'like', V_win);
    z_dst = (pre.pad_z + 1):(pre.pad_z + size(V_win, 1));
    x_dst = (pre.pad_x + 1):(pre.pad_x + size(V_win, 2));
    V_pad(z_dst, x_dst) = V_win;

    dkx = pre.dkx;
    dkz = pre.dkz;
    kx = pre.kx;
    kz = pre.kz;
    KR = pre.KR;
    TH = pre.TH;
    kmax = pre.kmax;
    Nbins_requested = pre.Nbins_requested;
    Nbins_effective = pre.Nbins_effective;
    k_edges = pre.k_edges;
    k_cent = pre.k_cent;
    radial_bin = pre.radial_bin;
end

% -------------------------------------------------------------------------
% 2D spectrum
% -------------------------------------------------------------------------

Uf = fftshift(fft2(V_pad));

if getfield_with_default(req, 'use_donut', false)

    Hd = build_donut_mask_local( ...
        KR, ...
        cfg.f0, ...
        req.donut_cs_min, ...
        req.donut_cs_max, ...
        req.donut_taper_rel);

    Uf = Uf .* Hd;

else

    Hd = [];

end

Sxz = abs(Uf).^2;

idx0z = floor(size(Sxz, 1)/2) + 1;
idx0x = floor(size(Sxz, 2)/2) + 1;

Sxz(idx0z, idx0x) = 0;

if req.smooth_sigma > 0
    Ssm = imgaussfilt(Sxz, req.smooth_sigma);
else
    Ssm = Sxz;
end

Ssm_norm = Ssm ./ max(Ssm(:) + eps);

% -------------------------------------------------------------------------
% Radial spectrum
% -------------------------------------------------------------------------

if isempty(pre)
    use_parfor = getfield_with_default(cfg, 'UseParfor', false);

    Srad = radial_average_power( ...
        Ssm, ...
        KR, ...
        k_edges, ...
        Nbins_effective, ...
        use_parfor);
else
    Srad = radial_average_power_precomputed( ...
        Ssm, radial_bin, Nbins_effective);
end

Srad = max(Srad, 0);
Srad(~isfinite(Srad)) = 0;

Srad_norm = Srad ./ max(Srad(:) + eps);

Ecum_raw = cumsum(Srad);

if Ecum_raw(end) > 0
    Ecum = Ecum_raw ./ Ecum_raw(end);
else
    Ecum = nan(size(Ecum_raw));
end

% -------------------------------------------------------------------------
% Evaluate q at the true wavenumber
% -------------------------------------------------------------------------

k_real = req.k0_true;

valid_interp = isfinite(k_cent) & isfinite(Ecum);

if nnz(valid_interp) >= 2

    q_val = interp1( ...
        k_cent(valid_interp), ...
        Ecum(valid_interp), ...
        k_real, ...
        'linear', ...
        'extrap');

else

    q_val = NaN;

end

q_val = max(0, min(1, q_val));

% -------------------------------------------------------------------------
% Output diagnostics
% -------------------------------------------------------------------------

req_curve = struct();

req_curve.Sxz = Sxz;
req_curve.Ssm = Ssm;
req_curve.Ssm_norm = Ssm_norm;

req_curve.Srad = Srad(:);
req_curve.Srad_norm = Srad_norm(:);
req_curve.Ecum = Ecum(:);

req_curve.q = q_val;
req_curve.q_theory = q_val;

req_curve.kx = kx(:).';
req_curve.kz = kz(:).';
req_curve.KR = KR;
req_curve.TH = TH;
req_curve.k_cent = k_cent(:);

req_curve.k_real = k_real;
req_curve.k0_true = k_real;

req_curve.f0 = cfg.f0;
req_curve.cs_true = cfg.cs_bg;

req_curve.Hd = Hd;

req_curve.Nbins_requested = Nbins_requested;
req_curve.Nbins_effective = Nbins_effective;

req_curve.dkx = dkx;
req_curve.dkz = dkz;
req_curve.kmax = kmax;

req_curve.M = getfield_with_default(req, 'M', NaN);
req_curve.win_size = getfield_with_default(req, 'win_size', NaN);
req_curve.half_win = getfield_with_default(req, 'half_win', NaN);
req_curve.pad_factor = getfield_with_default(req, 'pad_factor', NaN);
req_curve.smooth_sigma = req.smooth_sigma;

req_curve.use_donut = getfield_with_default(req, 'use_donut', false);
req_curve.donut_cs_min = getfield_with_default(req, 'donut_cs_min', NaN);
req_curve.donut_cs_max = getfield_with_default(req, 'donut_cs_max', NaN);
req_curve.donut_taper_rel = getfield_with_default(req, 'donut_taper_rel', NaN);

end

% =========================================================================
% Local helper functions
% =========================================================================

function pre = get_precomputed_geometry(req)

if isfield(req, 'precomputed') && isstruct(req.precomputed)
    pre = req.precomputed;
else
    pre = [];
end

end

function Srad = radial_average_power_precomputed(Ssm, radial_bin, Nbins)

w = Ssm(:);
valid = ~isnan(radial_bin) & isfinite(w);
Srad = accumarray(radial_bin(valid), w(valid), [Nbins, 1], @mean, 0);
Srad = Srad(:).';

end

function Srad = radial_average_power(Ssm, KR, k_edges, Nbins, use_parfor)

Srad = zeros(1, Nbins);

if use_parfor

    parfor b = 1:Nbins
        mask = (KR >= k_edges(b)) & (KR < k_edges(b+1));
        vals = Ssm(mask);

        vals = vals(isfinite(vals));

        if isempty(vals)
            Srad(b) = 0;
        else
            Srad(b) = mean(vals);
        end
    end

else

    for b = 1:Nbins
        mask = (KR >= k_edges(b)) & (KR < k_edges(b+1));
        vals = Ssm(mask);

        vals = vals(isfinite(vals));

        if isempty(vals)
            Srad(b) = 0;
        else
            Srad(b) = mean(vals);
        end
    end

end

end

function Hd = build_donut_mask_local(KR, f0, cs_min, cs_max, taper_rel)

k_low = 2*pi*f0 / cs_max;
k_high = 2*pi*f0 / cs_min;

Hd = double((KR >= k_low) & (KR <= k_high));

if nargin >= 5 && ~isempty(taper_rel) && taper_rel > 0

    bw = taper_rel * (k_high - k_low);

    k1 = k_low;
    k2 = k_low + bw;
    k3 = k_high - bw;
    k4 = k_high;

    Hd = zeros(size(KR));

    Hd(KR >= k2 & KR <= k3) = 1;

    idx = (KR >= k1 & KR < k2);
    x = (KR(idx) - k1) / max(k2 - k1, eps);
    Hd(idx) = 0.5 * (1 - cos(pi*x));

    idx = (KR > k3 & KR <= k4);
    x = (KR(idx) - k3) / max(k4 - k3, eps);
    Hd(idx) = 0.5 * (1 + cos(pi*x));

end

end

function tf = is_auto_nbins(x)

tf = (ischar(x) && strcmpi(x, 'auto')) || ...
     (isstring(x) && isscalar(x) && strcmpi(x, "auto"));

end

function validate_required_fields(s, names)

for i = 1:numel(names)
    if ~isfield(s, names{i})
        error('Missing required field: %s', names{i});
    end
end

end

function val = getfield_with_default(S, fieldName, defaultVal)

if isstruct(S) && isfield(S, fieldName)
    val = S.(fieldName);
else
    val = defaultVal;
end

end

function req_curve = empty_req_curve(req, cfg)

req_curve = struct();

req_curve.Sxz = [];
req_curve.Ssm = [];
req_curve.Ssm_norm = [];

req_curve.Srad = [];
req_curve.Srad_norm = [];
req_curve.Ecum = [];

req_curve.q = NaN;
req_curve.q_theory = NaN;

req_curve.kx = [];
req_curve.kz = [];
req_curve.k_cent = [];

req_curve.k_real = getfield_with_default(req, 'k0_true', NaN);
req_curve.k0_true = req_curve.k_real;

req_curve.f0 = getfield_with_default(cfg, 'f0', NaN);
req_curve.cs_true = getfield_with_default(cfg, 'cs_bg', NaN);

req_curve.Hd = [];

req_curve.Nbins_requested = getfield_with_default(req, 'Nbins', NaN);
req_curve.Nbins_effective = NaN;

req_curve.dkx = NaN;
req_curve.dkz = NaN;
req_curve.kmax = NaN;

end
