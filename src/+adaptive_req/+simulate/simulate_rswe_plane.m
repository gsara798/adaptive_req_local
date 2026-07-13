function [Uxz, x, z, cs_map, k_map, diag] = simulate_rswe_plane(varargin)
% SIMULATE_RSWE_PLANE
% Synthetic narrowband reverberant/shear-like field on imaging plane y = y0.
%
% Output convention:
%   Uxz(z,x), cs_map(z,x), k_map(z,x)
%
% Supports:
%   - Angular support:
%       * SourceSampling = 'ranges' or 'cone'
%   - Angular point generation:
%       * AngularSamplingMethod = 'random' or 'fibonacci'
%   - Wave models:
%       * WaveModel = 'spherical'
%       * WaveModel = 'planewave'   (homogeneous medium only)
%
% Notes:
%   - 'spherical' uses a shell of point-like emitters.
%   - 'planewave' uses true plane-wave superposition.
%   - In Is2D=true mode, strict in-plane SV-like polarization is imposed.

% =========================================================================
% Optional struct expansion
% =========================================================================
if ~isempty(varargin) && isstruct(varargin{1})
    cfg = varargin{1};
    varargin(1) = [];
    fn = fieldnames(cfg);
    for k = 1:numel(fn)
        varargin(end+1:end+2) = {fn{k}, cfg.(fn{k})};
    end
end

% =========================================================================
% Parse inputs
% =========================================================================
P = inputParser;
P.KeepUnmatched = false;

addParameter(P, 'MaskConfig', []);
addParameter(P, 'MaskType', 'homogeneous');
addParameter(P, 'cs_bg', 2.0);
addParameter(P, 'cs_inc', 3.0);
addParameter(P, 'MaskParams', struct());

addParameter(P, 'Lx', 0.05);
addParameter(P, 'Lz', 0.05);
addParameter(P, 'dx', 1e-4);
addParameter(P, 'dz', 1e-4);

addParameter(P, 'f0', 700);
addParameter(P, 'y0', 0.0);

addParameter(P, 'Nwaves', 30);
addParameter(P, 'rmin', []);
addParameter(P, 'rmax', []);

addParameter(P, 'AmpJitter', 0.10);
addParameter(P, 'DecayAlpha', 0);
addParameter(P, 'Seed', 1);
addParameter(P, 'UseParfor', true);

addParameter(P, 'PhiRange', [0, 2*pi]);
addParameter(P, 'ThetaRange', [0, pi]);
addParameter(P, 'SNR', Inf);

addParameter(P, 'Is2D', false, @(x) islogical(x) || isnumeric(x));
addParameter(P, 'AngleRange2D', [0, 2*pi]);

addParameter(P, 'SourceSampling', 'ranges');          % 'ranges', 'cone', or 'band'
addParameter(P, 'AngularSamplingMethod', 'random');   % 'random' or 'fibonacci'
addParameter(P, 'ForceInPlaneWave', false);
addParameter(P, 'ConeAxis', [-1 0 0]);
addParameter(P, 'ConeHalfAngleDeg', 180);
addParameter(P, 'BandAxis', [-1 0 0]);
addParameter(P, 'BandHalfWidthDeg', 90);


addParameter(P, 'WaveModel', 'spherical');            % 'spherical' or 'planewave'

parse(P, varargin{:});
P = P.Results;

if isempty(P.dz)
    P.dz = P.dx;
end

if P.ConeHalfAngleDeg < 0 || P.ConeHalfAngleDeg > 180
    error('ConeHalfAngleDeg must be in [0, 180].');
end

if P.BandHalfWidthDeg < 0 || P.BandHalfWidthDeg > 90
    error('BandHalfWidthDeg must be in [0, 90].');
end

waveModel = lower(string(P.WaveModel));
if ~ismember(waveModel, ["spherical", "planewave"])
    error('WaveModel must be ''spherical'' or ''planewave''.');
end

angMethod = lower(string(P.AngularSamplingMethod));
if ~ismember(angMethod, ["random", "fibonacci"])
    error('AngularSamplingMethod must be ''random'' or ''fibonacci''.');
end

omega = 2*pi*P.f0;

% =========================================================================
% Grid
% =========================================================================
Nx = round(P.Lx / P.dx) + 1;
Nz = round(P.Lz / P.dz) + 1;

x = linspace(0, P.Lx, Nx);
z = linspace(0, P.Lz, Nz);

% =========================================================================
% Shell radii defaults
% =========================================================================
if isempty(P.rmin), P.rmin = 0.9 * P.Lx; end
if isempty(P.rmax), P.rmax = 1.0 * P.Lx; end

% =========================================================================
% Build velocity / wavenumber maps
% Internal convention: (x,z) => [Nx, Nz]
% =========================================================================
if ~isempty(P.MaskConfig)
    mask_cfg = P.MaskConfig;
else
    mask_cfg = legacy_to_maskcfg(P);
end

[cs_map_xz, k_map_xz] = make_cs_mask_multi(x, z, omega, mask_cfg);

% =========================================================================
% Generate wave directions
% =========================================================================
rng(P.Seed);

[ux, uy, uz] = generate_wave_directions(P, angMethod);

if logical(P.ForceInPlaneWave)
    [ux, uy, uz] = force_in_plane_wave_direction(ux, uy, uz, P);
end

% =========================================================================
% Radial distances / source locations
% Only relevant for spherical mode
% =========================================================================
rho_u = rand(1, P.Nwaves, 'single');

if P.Is2D
    r_src = sqrt((single(P.rmax)^2 - single(P.rmin)^2) .* rho_u + single(P.rmin)^2);
else
    r_src = ((single(P.rmax)^3 - single(P.rmin)^3) .* rho_u + single(P.rmin)^3).^(1/3);
end

x_src = r_src .* ux + single(P.Lx / 2);
y_src = r_src .* uy + single(P.y0);
z_src = r_src .* uz + single(P.Lz / 2);

% =========================================================================
% Polarization projected onto z
% =========================================================================
pz = compute_pz(ux, uy, uz, P);

% =========================================================================
% Random amplitudes / phases
% =========================================================================
amp_jit = single(P.AmpJitter);
vamp = single(1 + amp_jit * randn(1, P.Nwaves, 'single'));
phi0 = 2*pi * rand(1, P.Nwaves, 'single');

gamma = (vamp .* pz .* exp(1i * phi0)) / sqrt(P.Nwaves);

% =========================================================================
% Field synthesis
% Internal convention: Uxz_xz(x,z)
% =========================================================================
Uxz_xz = complex(zeros(Nx, Nz, 'single'));
y0 = single(P.y0);

if P.UseParfor
    try
        if isempty(gcp('nocreate'))
            parpool('threads');
        end
    catch
    end
end

switch waveModel

    case "spherical"
        parfor (j = 1:Nz, bool2par(P.UseParfor))
            u_col = complex(zeros(Nx, 1, 'single'));
            xj = single(x(:));
            zj = single(z(j));
            kj = single(k_map_xz(:, j));

            for n = 1:P.Nwaves
                R = sqrt((xj - x_src(n)).^2 + ...
                         (y0 - y_src(n)).^2 + ...
                         (zj - z_src(n)).^2);

                if P.DecayAlpha ~= 0
                    Adec = 1 ./ max(R, 1e-6).^P.DecayAlpha;
                else
                    Adec = 1;
                end

                spatial = gamma(n) .* exp(1i * (kj .* R)) .* Adec;
                u_col = u_col + spatial;
            end

            Uxz_xz(:, j) = u_col;
        end

    case "planewave"
        k_pw = double(k_map_xz(1,1));
        rel_var = max(abs(double(k_map_xz(:)) - k_pw)) / max(abs(k_pw), eps);

        if rel_var > 1e-9
            error(['WaveModel = ''planewave'' currently requires a homogeneous medium. ', ...
                   'Detected nonuniform k_map. Relative variation = ', num2str(rel_var)]);
        end

        parfor (j = 1:Nz, bool2par(P.UseParfor))
            u_col = complex(zeros(Nx, 1, 'single'));
            xj = single(x(:));
            zj = single(z(j));

            for n = 1:P.Nwaves
                phase = k_pw * (ux(n) * xj + uy(n) * y0 + uz(n) * zj);
                u_col = u_col + gamma(n) .* exp(1i * phase);
            end

            Uxz_xz(:, j) = u_col;
        end
end

% =========================================================================
% Convert to output convention: (z,x)
% =========================================================================
% Uxz   = double(Uxz_xz);
% cs_map = cs_map_xz;
% k_map  = k_map_xz;
Uxz   = double(Uxz_xz.');
cs_map = cs_map_xz.';
k_map  = k_map_xz.';

% =========================================================================
% Add complex AWGN
% =========================================================================
if isfinite(P.SNR)
    sig_power = mean(abs(Uxz(:)).^2);
    noise_power = sig_power / (10^(P.SNR/10));

    dataType = class(Uxz);
    noise_real = sqrt(noise_power/2) * randn(size(Uxz), dataType);
    noise_imag = sqrt(noise_power/2) * randn(size(Uxz), dataType);

    Uxz = Uxz + complex(noise_real, noise_imag);
end

% =========================================================================
% Diagnostics
% =========================================================================
diag = struct();
diag.Params = P;

diag.x = x;
diag.z = z;
diag.dx = P.dx;
diag.dz = P.dz;
diag.f0 = P.f0;
diag.omega = omega;

diag.waveModel = char(waveModel);
diag.angularSamplingMethod = char(angMethod);

diag.waveDirs.ux = double(ux);
diag.waveDirs.uy = double(uy);
diag.waveDirs.uz = double(uz);
diag.wavePol.pz  = double(pz);
diag.inPlaneCoverage = ...
    adaptive_req.simulate.summarize_wave_direction_plane_coverage( ...
        diag.waveDirs, 'Tolerance', 1e-6);

if waveModel == "spherical"
    diag.sources.x = double(x_src);
    diag.sources.y = double(y_src);
    diag.sources.z = double(z_src);
else
    diag.sources.x = [];
    diag.sources.y = [];
    diag.sources.z = [];
end

diag.sources.w = double(gamma);
diag.outputConvention = 'Uxz(z,x), cs_map(z,x), k_map(z,x)';

end
% =========================== END MAIN FUNCTION ===========================


% =========================================================================
% Direction generation
% =========================================================================
function [ux, uy, uz] = generate_wave_directions(P, angMethod)

switch lower(P.SourceSampling)

    case 'ranges'
        if P.Is2D
            a0 = single(P.AngleRange2D(1));
            a1 = single(P.AngleRange2D(2));

            switch angMethod
                case "random"
                    u = rand(1, P.Nwaves, 'single');
                    alpha = a0 + u .* (a1 - a0);

                case "fibonacci"
                    alpha = linspace(a0, a1, P.Nwaves + 1);
                    alpha(end) = [];
                    alpha = single(alpha);
            end

            ux = cos(alpha);
            uy = zeros(1, P.Nwaves, 'single');
            uz = sin(alpha);

        else
            fullSphere = abs(P.PhiRange(1)) < 1e-12 && ...
                         abs(P.PhiRange(2) - 2*pi) < 1e-12 && ...
                         abs(P.ThetaRange(1)) < 1e-12 && ...
                         abs(P.ThetaRange(2) - pi) < 1e-12;

            if angMethod == "fibonacci" && fullSphere
                [ux, uy, uz] = fibonacci_sphere(P.Nwaves);
                ux = single(ux); uy = single(uy); uz = single(uz);
            else
                th_min = single(P.ThetaRange(1));
                th_max = single(P.ThetaRange(2));
                ph_min = single(P.PhiRange(1));
                ph_max = single(P.PhiRange(2));

                switch angMethod
                    case "random"
                        u1 = rand(1, P.Nwaves, 'single');
                        u2 = rand(1, P.Nwaves, 'single');

                        cos_th_min = cos(th_min);
                        cos_th_max = cos(th_max);

                        cos_theta = cos_th_min - u1 .* (cos_th_min - cos_th_max);
                        theta = acos(cos_theta);
                        phi = ph_min + u2 .* (ph_max - ph_min);

                    case "fibonacci"
                        t = ((0:P.Nwaves-1) + 0.5) / P.Nwaves;
                        t = single(t);

                        cos_th_min = cos(th_min);
                        cos_th_max = cos(th_max);

                        cos_theta = cos_th_min - t .* (cos_th_min - cos_th_max);
                        theta = acos(cos_theta);

                        phi_span = ph_max - ph_min;
                        phi = ph_min + mod((0:P.Nwaves-1) * (2*pi / ((1+sqrt(5))/2)), phi_span);
                        phi = single(phi);
                end

                ux = sin(theta) .* cos(phi);
                uy = sin(theta) .* sin(phi);
                uz = cos(theta);
            end
        end

    case 'cone'
        if P.Is2D
            axis_xz = single([P.ConeAxis(1), P.ConeAxis(3)]);
            na2 = norm(axis_xz);

            if na2 < 1e-12
                error('For Is2D=true, ConeAxis must have a nonzero x-z projection.');
            end

            axis_xz = axis_xz / na2;
            alpha0 = atan2(axis_xz(2), axis_xz(1));
            half_ang = deg2rad(single(P.ConeHalfAngleDeg));

            switch angMethod
                case "random"
                    u = rand(1, P.Nwaves, 'single');
                    alpha = alpha0 + (2*u - 1) .* half_ang;

                case "fibonacci"
                    t = ((0:P.Nwaves-1) + 0.5) / P.Nwaves;
                    t = single(t);
                    alpha = alpha0 + (-1 + 2*t) .* half_ang;
            end

            ux = cos(alpha);
            uy = zeros(1, P.Nwaves, 'single');
            uz = sin(alpha);

        else
            axis_hat = single(P.ConeAxis(:));
            na = norm(axis_hat);
            if na < 1e-12
                error('ConeAxis must be nonzero.');
            end
            axis_hat = axis_hat / na;

            half_ang = deg2rad(single(P.ConeHalfAngleDeg));

            switch angMethod
                case "random"
                    [ux, uy, uz] = sample_unit_vectors_in_cone(axis_hat, half_ang, P.Nwaves);

                case "fibonacci"
                    [ux, uy, uz] = fibonacci_cap(axis_hat, half_ang, P.Nwaves);
            end
        end

    case 'band'
        if P.Is2D
            axis_xz = single([P.BandAxis(1), P.BandAxis(3)]);
            na2 = norm(axis_xz);

            if na2 < 1e-12
                error('For Is2D=true, BandAxis must have a nonzero x-z projection.');
            end

            axis_xz = axis_xz / na2;
            alpha0 = atan2(axis_xz(2), axis_xz(1));
            half_w = deg2rad(single(P.BandHalfWidthDeg));

            % In 2D, a band centered on the equator of the axis becomes
            % two symmetric sectors around alpha0 +/- pi/2.
            c1 = alpha0 + pi/2;
            c2 = alpha0 - pi/2;

            switch angMethod
                case "random"
                    u = rand(1, P.Nwaves, 'single');
                    side = rand(1, P.Nwaves, 'single') > 0.5;

                    alpha = zeros(1, P.Nwaves, 'single');
                    alpha(~side) = c1 + (-1 + 2*u(~side)) .* half_w;
                    alpha(side)  = c2 + (-1 + 2*u(side))  .* half_w;

                case "fibonacci"
                    t = ((0:P.Nwaves-1) + 0.5) / P.Nwaves;
                    t = single(t);

                    alpha = zeros(1, P.Nwaves, 'single');
                    idx1 = 1:2:P.Nwaves;
                    idx2 = 2:2:P.Nwaves;

                    alpha(idx1) = c1 + (-1 + 2*t(idx1)) .* half_w;
                    alpha(idx2) = c2 + (-1 + 2*t(idx2)) .* half_w;
            end

            ux = cos(alpha);
            uy = zeros(1, P.Nwaves, 'single');
            uz = sin(alpha);

        else
            axis_hat = single(P.BandAxis(:));
            na = norm(axis_hat);
            if na < 1e-12
                error('BandAxis must be nonzero.');
            end
            axis_hat = axis_hat / na;

            half_w = deg2rad(single(P.BandHalfWidthDeg));

            switch angMethod
                case "random"
                    [ux, uy, uz] = sample_unit_vectors_in_band(axis_hat, half_w, P.Nwaves);

                case "fibonacci"
                    [ux, uy, uz] = fibonacci_band(axis_hat, half_w, P.Nwaves);
            end
        end

    otherwise
        error('Unknown SourceSampling mode: %s', P.SourceSampling);
end

end

% =========================================================================
% Polarization
% =========================================================================
function pz = compute_pz(ux, uy, uz, P)

if P.Is2D
    pz = ux;
    return
end

N = numel(ux);
pz = zeros(1, N, 'single');

for n = 1:N
    k_dir = [ux(n), uy(n), uz(n)];

    rand_vec = randn(1, 3, 'single');
    p_vec = rand_vec - dot(rand_vec, k_dir) * k_dir;
    np = norm(p_vec);

    if np < 1e-10
        rand_vec = randn(1, 3, 'single');
        p_vec = rand_vec - dot(rand_vec, k_dir) * k_dir;
        np = norm(p_vec);

        if np < 1e-10
            if abs(k_dir(1)) < 0.9
                p_vec = [1, 0, 0] - dot([1, 0, 0], k_dir) * k_dir;
            else
                p_vec = [0, 1, 0] - dot([0, 1, 0], k_dir) * k_dir;
            end
            np = norm(p_vec);
        end
    end

    p_vec = p_vec / np;
    pz(n) = p_vec(3);
end

end

function [ux, uy, uz] = sample_unit_vectors_in_band(axis_hat, half_w, N)
% Uniform sampling over a spherical band centered at 90 deg from axis_hat:
% beta in [pi/2 - half_w, pi/2 + half_w]

axis_hat = axis_hat(:);
axis_hat = axis_hat / norm(axis_hat);

if abs(dot(axis_hat, [0;0;1])) < 0.99
    tmp = [0;0;1];
else
    tmp = [0;1;0];
end

e1 = cross(tmp, axis_hat);
e1 = e1 / norm(e1);

e2 = cross(axis_hat, e1);
e2 = e2 / norm(e2);

% beta in [pi/2-half_w, pi/2+half_w]
% so mu = cos(beta) in [-sin(half_w), +sin(half_w)]
u = rand(1, N, 'single');
mu_max = sin(half_w);
mu = -mu_max + 2*mu_max*u;

sin_beta = sqrt(max(0, 1 - mu.^2));
psi = 2*pi * rand(1, N, 'single');

dirs = axis_hat * mu + ...
       e1 * (sin_beta .* cos(psi)) + ...
       e2 * (sin_beta .* sin(psi));

ux = dirs(1,:);
uy = dirs(2,:);
uz = dirs(3,:);
end

function [ux, uy, uz] = fibonacci_band(axis_hat, half_w, N)
% Quasi-uniform sampling over a spherical band centered at 90 deg from axis_hat

axis_hat = axis_hat(:);
axis_hat = axis_hat / norm(axis_hat);

if abs(dot(axis_hat, [0;0;1])) < 0.99
    tmp = [0;0;1];
else
    tmp = [0;1;0];
end

e1 = cross(tmp, axis_hat);
e1 = e1 / norm(e1);

e2 = cross(axis_hat, e1);
e2 = e2 / norm(e2);

g = pi * (3 - sqrt(5));
k = 0:N-1;

t = (k + 0.5) / N;

mu_max = sin(half_w);
mu = -mu_max + 2*mu_max*t;   % quasi-uniform in cos(beta)

sin_beta = sqrt(max(0, 1 - mu.^2));
psi = g * k;

dirs = axis_hat * mu + ...
       e1 * (sin_beta .* cos(psi)) + ...
       e2 * (sin_beta .* sin(psi));

ux = dirs(1,:);
uy = dirs(2,:);
uz = dirs(3,:);
end

% =========================================================================
% Helper utilities
% =========================================================================
function cfg = legacy_to_maskcfg(P)
cfg = struct();
cfg.cs_bg = P.cs_bg;
cfg.CombineMode = 'overlay';

if strcmpi(P.MaskType, 'homogeneous')
    cfg.Masks = {struct('Type', 'homogeneous', 'cs_inc', P.cs_bg, 'Params', struct())};
else
    cfg.Masks = {struct('Type', P.MaskType, 'cs_inc', P.cs_inc, 'Params', P.MaskParams)};
end
end

function tf = bool2par(flag)
tf = 1;
if ~flag, tf = 0; end
end

function [cs_map, k_map, alpha_stack] = make_cs_mask_multi(x, z, omega, cfg)

if ~isfield(cfg, 'cs_bg'),       cfg.cs_bg = 2.0; end
if ~isfield(cfg, 'Masks'),       cfg.Masks = {};  end
if ~isfield(cfg, 'CombineMode'), cfg.CombineMode = 'overlay'; end

Nx = numel(x);
Nz = numel(z);
[X, Z] = ndgrid(x, z);

if isempty(cfg.Masks)
    cs_map = cfg.cs_bg * ones(Nx, Nz);
    k_map  = omega ./ cs_map;
    alpha_stack = {};
    return
end

M = numel(cfg.Masks);
alpha_stack = cell(1, M);
cs_inc_list = zeros(1, M);

for m = 1:M
    ms = cfg.Masks{m};
    if ~isfield(ms, 'Params'), ms.Params = struct(); end
    alpha_stack{m} = alpha_from_mask(X, Z, ms.Type, ms.Params);
    cs_inc_list(m) = ms.cs_inc;
end

switch lower(cfg.CombineMode)
    case 'overlay'
        cs_map = cfg.cs_bg * ones(Nx, Nz);
        for m = 1:M
            a = alpha_stack{m};
            cs_map = (1 - a) .* cs_map + a .* cs_inc_list(m);
        end

    case 'blend'
        num = zeros(Nx, Nz);
        den = zeros(Nx, Nz);
        for m = 1:M
            a = alpha_stack{m};
            num = num + a .* cs_inc_list(m);
            den = den + a;
        end
        a_tot = min(den, 1);
        cs_map = (1 - a_tot) .* cfg.cs_bg + a_tot .* (num ./ max(den, eps));

    case 'max'
        cs_map = cfg.cs_bg * ones(Nx, Nz);
        for m = 1:M
            a = alpha_stack{m};
            c_local = (1 - a) .* cfg.cs_bg + a .* cs_inc_list(m);
            cs_map = max(cs_map, c_local);
        end

    case 'min'
        cs_map = cfg.cs_bg * ones(Nx, Nz);
        for m = 1:M
            a = alpha_stack{m};
            c_local = (1 - a) .* cfg.cs_bg + a .* cs_inc_list(m);
            cs_map = min(cs_map, c_local);
        end

    otherwise
        error('CombineMode "%s" not supported.', cfg.CombineMode);
end

k_map = omega ./ cs_map;
end

function a = alpha_from_mask(X, Z, type, P)

if ~isfield(P, 'SigmaEdge'), P.SigmaEdge = 0; end

switch lower(type)
    case 'homogeneous'
        a = zeros(size(X));

    case 'circle'
        if ~isfield(P, 'Center'), P.Center = [0, 0]; end
        if ~isfield(P, 'Radius'), P.Radius = 1e-2; end

        xc = P.Center(1);
        zc = P.Center(2);
        r  = sqrt((X - xc).^2 + (Z - zc).^2);
        M0 = (r <= P.Radius);
        a  = soften(M0, P.SigmaEdge, X, Z);

    case 's-curve'
        if ~isfield(P, 'Center'),       P.Center = [0, 0]; end
        if ~isfield(P, 'S_Amplitude'),  P.S_Amplitude = 6e-3; end
        if ~isfield(P, 'S_Wavelength'), P.S_Wavelength = 25e-3; end
        if ~isfield(P, 'S_Thickness'),  P.S_Thickness = 10e-3; end

        X0 = X - P.Center(1);
        Z0 = Z - P.Center(2);
        zS = P.S_Amplitude .* sin(2*pi*X0 / max(P.S_Wavelength, eps));
        M0 = abs(Z0 - zS) <= (P.S_Thickness/2);
        a = soften(M0, P.SigmaEdge, X, Z);

    case 'bilayer'
        if ~isfield(P, 'Bi_Angle'),  P.Bi_Angle = 0; end
        if ~isfield(P, 'Bi_Offset'), P.Bi_Offset = 0; end

        n = [cos(P.Bi_Angle); sin(P.Bi_Angle)];
        d = X*n(1) + Z*n(2) - P.Bi_Offset;
        M0 = (d > 0);
        a = soften(M0, P.SigmaEdge, X, Z);

    case 'custom'
        if ~isfield(P, 'CustomMask')
            error('Custom requires Params.CustomMask.');
        end
        M0 = logical(P.CustomMask);

        if ~isequal(size(M0), size(X))
            error('CustomMask must match grid size.');
        end
        a = soften(M0, P.SigmaEdge, X, Z);

    otherwise
        error('Unknown mask Type: %s', type);
end
end

function a = soften(M0, sigma_edge, X, Z)

if sigma_edge <= 0
    a = double(M0);
    return
end

dx = mean(diff(unique(X(:,1))));
dz = mean(diff(unique(Z(1,:))));
sig_px = max(0.5, sigma_edge / max(dx, dz));
fsz = 2*ceil(3*sig_px) + 1;

a = imgaussfilt(double(M0), sig_px, 'FilterSize', fsz, 'Padding', 'replicate');
a = a / max(a(:) + eps);
end

function [ux, uy, uz] = sample_unit_vectors_in_cone(axis_hat, half_ang, N)
axis_hat = axis_hat(:);
axis_hat = axis_hat / norm(axis_hat);

if abs(dot(axis_hat, [0;0;1])) < 0.99
    tmp = [0;0;1];
else
    tmp = [0;1;0];
end

e1 = cross(tmp, axis_hat);
e1 = e1 / norm(e1);

e2 = cross(axis_hat, e1);
e2 = e2 / norm(e2);

u = rand(1, N, 'single');
cos_beta = 1 - u .* (1 - cos(half_ang));
sin_beta = sqrt(max(0, 1 - cos_beta.^2));

psi = 2*pi * rand(1, N, 'single');

dirs = axis_hat * cos_beta + ...
       e1 * (sin_beta .* cos(psi)) + ...
       e2 * (sin_beta .* sin(psi));

ux = dirs(1,:);
uy = dirs(2,:);
uz = dirs(3,:);
end

function [ux, uy, uz] = fibonacci_sphere(N)
g = pi * (3 - sqrt(5));
k = 0:N-1;

z = 1 - 2*(k + 0.5)/N;
r = sqrt(max(0, 1 - z.^2));
phi = g * k;

ux = r .* cos(phi);
uy = r .* sin(phi);
uz = z;
end

function [ux, uy, uz] = fibonacci_cap(axis_hat, half_ang, N)
axis_hat = axis_hat(:);
axis_hat = axis_hat / norm(axis_hat);

if abs(dot(axis_hat, [0;0;1])) < 0.99
    tmp = [0;0;1];
else
    tmp = [0;1;0];
end

e1 = cross(tmp, axis_hat);
e1 = e1 / norm(e1);
e2 = cross(axis_hat, e1);
e2 = e2 / norm(e2);

g = pi * (3 - sqrt(5));
k = 0:N-1;

t = (k + 0.5) / N;
cos_beta = 1 - t .* (1 - cos(half_ang));
sin_beta = sqrt(max(0, 1 - cos_beta.^2));
psi = g * k;

dirs = axis_hat * cos_beta + ...
       e1 * (sin_beta .* cos(psi)) + ...
       e2 * (sin_beta .* sin(psi));

ux = dirs(1,:);
uy = dirs(2,:);
uz = dirs(3,:);
end

function [ux, uy, uz] = force_in_plane_wave_direction(ux, uy, uz, P)

switch lower(P.SourceSampling)
    case 'cone'
        axis_hat = double(P.ConeAxis(:));
        axis_hat = axis_hat / norm(axis_hat);
        target = [axis_hat(1); 0; axis_hat(3)];

        if norm(target) < eps
            return;
        end

        target = target / norm(target);
        sep = acos(max(-1, min(1, dot(axis_hat, target))));

        if sep > deg2rad(double(P.ConeHalfAngleDeg)) + 1e-12
            warning(['ForceInPlaneWave requested, but the cone does not ', ...
                'intersect the imaging plane. No direction was changed.']);
            return;
        end

    case 'ranges'
        target = [1; 0; 0];

    case 'band'
        axis_hat = double(P.BandAxis(:));
        axis_hat = axis_hat / norm(axis_hat);
        target = cross(axis_hat, [0; 1; 0]);

        if norm(target) < eps
            target = [1; 0; 0];
        else
            target = target / norm(target);
        end

    otherwise
        return;
end

[~, idx] = min(abs(double(uy)));
ux(idx) = single(target(1));
uy(idx) = single(0);
uz(idx) = single(target(3));

end
