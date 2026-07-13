function out = q_theory_REQ_discrete_shearUZ(dx, dz, f0, cs_guess, varargin)
% q_theory_REQ_discrete_shearUZ
% Discrete-theory quantile q_th = Ecum(k0) matched to REQ implementation.
%
% This version also computes:
%   - Srad_true : ideal radial spectrum on the same discrete k grid
%   - Ecum_true : cumulative ideal radial spectrum
%   - q_ideal   : ideal quantile evaluated at k0
%
% TheoryMode:
%   - 'S2D'  : Build full ideal 2D spectrum, optionally apply donut filter
%              in 2D, convolve in 2D, then radialize
%   - 'SRAD' : Build ideal radial spectrum directly, then apply a radial
%              transfer operator derived from the exact 2D windowing pipeline
%
% FieldType:
%   - 'Diffuse3D'
%   - 'Diffuse2D'
%   - 'SingleWave'
%
% Notes
% -----
% - The donut filter is applied ONLY in the S2D route.
% - q_ideal remains the unfiltered ideal quantile.
% - q_th reflects the filtered discrete/windowed theory when
%   UseDonutFilter = true.

% ---------------- Parse ----------------
p = inputParser;
p.KeepUnmatched = false;

addParameter(p, 'M', 4);
addParameter(p, 'Gamma', 1.0);
addParameter(p, 'PadFactor', 1);
addParameter(p, 'Nbins', 120);
addParameter(p, 'SmoothSigma', 1);
addParameter(p, 'CraterPower', 1.0);
addParameter(p, 'FieldType', 'Diffuse3D');
addParameter(p, 'TheoryMode', 'S2D');
addParameter(p, 'Plot', true);
addParameter(p, 'BuildTransferMatrix', true);

% ---- Donut filter options (used only in S2D) ----
addParameter(p, 'UseDonutFilter', false);
addParameter(p, 'DonutCsMin', 1.0);
addParameter(p, 'DonutCsMax', 10.0);
addParameter(p, 'DonutTaperRel', 0.0);

parse(p, varargin{:});
P = p.Results;

omega   = 2*pi*f0;
lambda0 = cs_guess / f0;
k0      = omega / cs_guess;

% ----- Window sizes and padding (same logic as REQ) -----
nx = round(P.M * lambda0 / dx);
nz = round(P.M * lambda0 / dz);
if mod(nx,2)==0, nx = nx+1; end
if mod(nz,2)==0, nz = nz+1; end

W = make_window_hann_circular_m(nx, nz, P.Gamma, dx, dz);

if P.PadFactor > 0
    PAD  = max(0, round(double([P.PadFactor * nz, P.PadFactor * nx])));
    Wpad = padarray(W, PAD, 0, 'both');
else
    Wpad = W;
end

% ----- FFT grid exactly like REQ -----
[Z2, X2] = size(Wpad);

dkx = 2*pi/(X2*dx);
dkz = 2*pi/(Z2*dz);

kx = (-floor(X2/2):ceil(X2/2)-1) * dkx;
kz = (-floor(Z2/2):ceil(Z2/2)-1) * dkz;

[KX, KZ] = meshgrid(kx, kz);
KR = sqrt(KX.^2 + KZ.^2);

kmax   = min(max(abs(kx)), max(abs(kz)));
dk_rad = min(dkx, dkz);

% ----- Radial bins exactly like REQ -----
if ischar(P.Nbins) && strcmpi(P.Nbins, 'auto')
    oversample_factor = 1;
    Nbins = round(kmax / (dk_rad / oversample_factor));
    Nbins = max(16, Nbins);
else
    Nbins = max(16, P.Nbins);
end

k_edges = linspace(0, kmax, Nbins+1);
k_cent  = 0.5 * (k_edges(1:end-1) + k_edges(2:end));

% ----- Window power spectrum -----
Wk = fftshift(fft2(Wpad));
Wpow2D = abs(Wk).^2;
Wpow2D = Wpow2D / (sum(Wpow2D(:)) + eps);

% ----- Radialized window profile (diagnostic only) -----
Wpow_rad = radial_average_mean(Wpow2D, KR, k_edges);
Wpow_rad = Wpow_rad / (sum(Wpow_rad) + eps);

% ----- Donut mask (diagnostic / S2D only) -----
Hd_donut = [];
if P.UseDonutFilter
    Hd_donut = build_donut_mask_2D(KR, f0, P.DonutCsMin, P.DonutCsMax, P.DonutTaperRel);
end

% Initialize outputs
Strue2D   = [];
Smean2D   = [];
Srad_true = [];
Srad_mean = [];
T_rad     = [];
Srad      = [];
Ecum      = [];
Ecum_true = [];
q_th      = NaN;
q_ideal   = NaN;

% ================================================================
% ALWAYS BUILD THE IDEAL RADIAL SPECTRUM ON THE SAME k GRID
% (kept unfiltered on purpose so q_ideal remains the pure ideal limit)
% ================================================================
Srad_true = build_Srad_true(k_cent, k0, dk_rad, P);
Srad_true(~isfinite(Srad_true)) = 0;
Srad_true = max(Srad_true, 0);

if sum(Srad_true) > eps
    switch upper(P.FieldType)
        case {'SINGLEWAVE','DIFFUSE2D'}
            Ecum_true = zeros(size(k_cent));
            [~, idx0] = min(abs(k_cent - k0));
            Ecum_true(1:idx0-1) = 0;
            Ecum_true(idx0) = 0.5;
            Ecum_true(idx0+1:end) = 1.0;
            q_ideal = 0.5;

        otherwise
            Ecum_true = cumsum(Srad_true);
            Ecum_true = Ecum_true / (Ecum_true(end) + eps);
            q_ideal = interp1_monotone(k_cent, Ecum_true, k0);
    end
else
    Ecum_true = nan(size(Srad_true));
    q_ideal = NaN;
end

switch upper(P.TheoryMode)

    case 'S2D'
        % ================================================================
        % ROUTE A: BUILD FULL Strue2D, OPTIONAL DONUT IN 2D, WINDOW IN 2D,
        % THEN RADIALIZE
        % ================================================================
        Strue2D = build_Strue2D(KX, KZ, KR, k0, dk_rad, P);

        % Optional donut filter in k-space:
        % filtered power spectrum = |Hd|^2 * Strue2D
        if P.UseDonutFilter
            Strue2D = Strue2D .* (Hd_donut.^2);
        end

        Strue2D(~isfinite(Strue2D)) = 0;
        Strue2D = max(Strue2D, 0);
        Strue2D = Strue2D / (sum(Strue2D(:)) + eps);

        % Exact 2D convolution with window power
        Smean2D = fftconv2_same(Strue2D, Wpow2D);

        % Mimic REQ: remove DC bin
        idx_DC_z = floor(Z2/2) + 1;
        idx_DC_x = floor(X2/2) + 1;
        Smean2D(idx_DC_z, idx_DC_x) = 0;

        if P.SmoothSigma > 0
            Smean2D = imgaussfilt(Smean2D, P.SmoothSigma);
        end

        Smean2D(~isfinite(Smean2D)) = 0;
        Smean2D = max(Smean2D, 0);

        % Radial averaging EXACTLY like REQ
        Srad = radial_average_mean(Smean2D, KR, k_edges);

        % Keep for plotting convenience
        Srad_mean = Srad;

    case 'SRAD'
        % ================================================================
        % ROUTE B: UNCHANGED
        % ================================================================
        Srad_shape = Srad_true;

        if P.SmoothSigma > 0
            Srad_shape = imgaussfilt(Srad_shape, P.SmoothSigma);
        end

        Srad_shape = Srad_shape / (sum(Srad_shape) + eps);

        if P.BuildTransferMatrix
            T_rad = build_radial_transfer_matrix(KR, k_edges, Wpow2D, P.SmoothSigma);
            Srad_mean = T_rad * Srad_shape(:);
            Srad_mean = Srad_mean(:).';
        else
            kernel1D = radial_average_sum(Wpow2D, KR, k_edges);
            kernel1D = kernel1D / (sum(kernel1D) + eps);
            Srad_mean = conv(Srad_shape, kernel1D, 'same');
        end

        Srad_mean(~isfinite(Srad_mean)) = 0;
        Srad_mean = max(Srad_mean, 0);

        if sum(Srad_mean) > eps
            Srad_mean = Srad_mean / sum(Srad_mean);
        end

        Srad = Srad_mean;

    otherwise
        error('Invalid TheoryMode. Choose ''S2D'' or ''SRAD''.');
end

% ---------------- Ecum and q_th ----------------
Srad(~isfinite(Srad)) = 0;
Srad = max(Srad, 0);

if sum(Srad) <= eps
    warning('Srad collapsed to zero. Returning NaN outputs.');
    Ecum = nan(size(Srad));
    q_th = NaN;
else
    Ecum = cumsum(Srad);
    Ecum = Ecum / (Ecum(end) + eps);
    q_th = interp1_monotone(k_cent, Ecum, k0);
end

% ---------------- Pack outputs ----------------
out = struct();
out.q_th       = q_th;
out.q_ideal    = q_ideal;
out.k0         = k0;
out.k_cent     = k_cent;
out.k_edges    = k_edges;
out.Srad       = Srad;
out.Srad_true  = Srad_true;
out.Srad_mean  = Srad_mean;
out.Ecum       = Ecum;
out.Ecum_true  = Ecum_true;
out.kx         = kx;
out.kz         = kz;
out.KR         = KR;
out.Wpow2D     = Wpow2D;
out.Wpow_rad   = Wpow_rad;
out.Strue2D    = Strue2D;
out.Smean2D    = Smean2D;
out.T_rad      = T_rad;
out.Hd_donut   = Hd_donut;
out.params     = P;

% ---------------- Optional quick plot ----------------
if P.Plot
    figure('Color', 'w', 'Position', [100 100 1500 420]);

    subplot(1,4,1);
    imagesc(kx, kz, Strue2D);
    axis image;
    title(sprintf('Ideal S_{2D} (%s)', P.FieldType), 'Interpreter', 'none');
    xlabel('k_x'); ylabel('k_z');
    xlim([-2*k0 2*k0]); ylim([-2*k0 2*k0]);
    colorbar;

    subplot(1,4,2);
    imagesc(kx, kz, Smean2D);
    axis image;
    title('Discrete/windowed S_{2D}');
    xlabel('k_x'); ylabel('k_z');
    xlim([-2*k0 2*k0]); ylim([-2*k0 2*k0]);
    colorbar;

    subplot(1,4,3);
    hold on;
    plot(k_cent, Srad_true./(max(Srad_true)+eps), 'k-', 'LineWidth', 2);
    plot(k_cent, Srad./(max(Srad)+eps), 'b-', 'LineWidth', 2);
    xline(k0, 'r--', 'k_0');
    xlabel('k');
    ylabel('Normalized \bar{S}_{rad}');
    title('Radial spectrum');
    legend('Ideal', 'Discrete', 'Location', 'best');
    grid on;
    xlim([0 2.5*k0]);

    subplot(1,4,4);
    hold on;
    plot(k_cent, Ecum_true, 'k-', 'LineWidth', 2);
    plot(k_cent, Ecum, 'b-', 'LineWidth', 2);
    xline(k0, 'r--', 'k_0');
    if isfinite(q_ideal), plot(k0, q_ideal, 'ko', 'MarkerFaceColor', 'k'); end
    if isfinite(q_th),    plot(k0, q_th,    'bo', 'MarkerFaceColor', 'b'); end
    xlabel('k');
    ylabel('E_{cum}(k)');
    title(sprintf('Ecum: q_{ideal}=%.3f, q_{th}=%.3f', q_ideal, q_th));
    legend('Ideal', 'Discrete', 'Location', 'best');
    grid on;
    xlim([0 2.5*k0]);
    ylim([0 1]);
end

end

% =====================================================================
% LOCAL HELPERS
% =====================================================================

function q = interp1_monotone(k_cent, Ecum, k0)
    valid = isfinite(k_cent) & isfinite(Ecum);
    kc_tmp = k_cent(valid);
    Ec_tmp = Ecum(valid);

    [kc_u, iu] = unique(kc_tmp, 'stable');
    Ec_u = Ec_tmp(iu);

    if numel(kc_u) < 2
        q = NaN;
        return;
    end

    q = interp1(kc_u, Ec_u, k0, 'linear', 'extrap');
end

function kq = invert_ecum_to_k(k_cent, Ecum, q)
    valid = isfinite(k_cent) & isfinite(Ecum);
    kc_tmp = k_cent(valid);
    Ec_tmp = Ecum(valid);

    [Ec_u, iu] = unique(Ec_tmp, 'stable');
    kc_u = kc_tmp(iu);

    if numel(Ec_u) < 2
        kq = NaN;
        return;
    end

    kq = interp1(Ec_u, kc_u, q, 'linear', 'extrap');
end

function Strue2D = build_Strue2D(KX, KZ, KR, k0, dk_rad, P)
    Strue2D = zeros(size(KR));

    switch upper(P.FieldType)
        case 'DIFFUSE3D'
            inside = (KR < k0);
            edge_eps = max(dk_rad^2, 1e-12);

            Strue2D(inside) = (1 - (KZ(inside).^2)/(k0^2)) ./ ...
                              sqrt(max(k0^2 - KR(inside).^2, edge_eps));

            if ~isempty(P.CraterPower) && P.CraterPower ~= 1
                Strue2D(inside) = Strue2D(inside) .^ P.CraterPower;
            end

        case 'DIFFUSE2D'
            tol = 0.5 * dk_rad;
            ring_shape = double(abs(KR - k0) <= tol);
            polarization_weight = (KX ./ (KR + eps)).^2;
            Strue2D = polarization_weight .* ring_shape;

        case 'SINGLEWAVE'
            kx = KX(1,:);
            kz = KZ(:,1);
            [~, idx_x_pos] = min(abs(kx - k0));
            [~, idx_x_neg] = min(abs(kx + k0));
            [~, idx_z_0]   = min(abs(kz - 0));
            Strue2D(idx_z_0, idx_x_pos) = 1;
            Strue2D(idx_z_0, idx_x_neg) = 1;

        otherwise
            error('Invalid FieldType.');
    end
end

function Srad_true = build_Srad_true(k_cent, k0, dk_rad, P)
    Srad_true = zeros(size(k_cent));

    switch upper(P.FieldType)
        case 'DIFFUSE3D'
            inside = (k_cent < k0);
            edge_eps = max(dk_rad^2, 1e-12);

            Srad_true(inside) = ...
                (2*k0^2 - k_cent(inside).^2) ./ ...
                sqrt(max(k0^2 - k_cent(inside).^2, edge_eps));

            if ~isempty(P.CraterPower) && P.CraterPower ~= 1
                Srad_true(inside) = Srad_true(inside) .^ P.CraterPower;
            end

        case 'DIFFUSE2D'
            Srad_true = zeros(size(k_cent));
            [~, idx0] = min(abs(k_cent - k0));
            Srad_true(idx0) = 1;

        case 'SINGLEWAVE'
            Srad_true = zeros(size(k_cent));
            [~, idx0] = min(abs(k_cent - k0));
            Srad_true(idx0) = 1;

        otherwise
            error('Invalid FieldType.');
    end
end

function T = build_radial_transfer_matrix(KR, k_edges, Wpow2D, smoothSigma)
    Nbins = numel(k_edges) - 1;
    T = zeros(Nbins, Nbins);

    ring_masks = cell(1, Nbins);
    for b = 1:Nbins
        ring_masks{b} = (KR >= k_edges(b) & KR < k_edges(b+1));
    end

    for j = 1:Nbins
        Bj = zeros(size(KR));
        Bj(ring_masks{j}) = 1;

        prof_in = zeros(1, Nbins);
        for b = 1:Nbins
            vals = Bj(ring_masks{b});
            if isempty(vals)
                prof_in(b) = 0;
            else
                prof_in(b) = mean(vals, 'omitnan');
            end
        end

        gain_in = prof_in(j);
        if gain_in <= eps
            gain_in = 1;
        end

        Bj_conv = fftconv2_same(Bj, Wpow2D);

        if smoothSigma > 0
            Bj_conv = imgaussfilt(Bj_conv, smoothSigma);
        end

        Bj_conv(~isfinite(Bj_conv)) = 0;
        Bj_conv = max(Bj_conv, 0);

        prof_out = zeros(1, Nbins);
        for b = 1:Nbins
            vals = Bj_conv(ring_masks{b});
            if isempty(vals)
                prof_out(b) = 0;
            else
                prof_out(b) = mean(vals, 'omitnan');
            end
        end

        T(:, j) = prof_out(:) / gain_in;
    end
end

function prof = radial_average_mean(S, KR, k_edges)
    Nbins = numel(k_edges) - 1;
    prof = zeros(1, Nbins);

    for b = 1:Nbins
        mask = (KR >= k_edges(b) & KR < k_edges(b+1));
        vals = S(mask);
        if isempty(vals)
            prof(b) = 0;
        else
            prof(b) = mean(vals, 'omitnan');
        end
    end

    prof(~isfinite(prof)) = 0;
end

function prof = radial_average_sum(S, KR, k_edges)
    Nbins = numel(k_edges) - 1;
    prof = zeros(1, Nbins);

    for b = 1:Nbins
        mask = (KR >= k_edges(b) & KR < k_edges(b+1));
        vals = S(mask);
        if isempty(vals)
            prof(b) = 0;
        else
            prof(b) = sum(vals, 'omitnan');
        end
    end

    prof(~isfinite(prof)) = 0;
end

function W = make_window_hann_circular_m(nx, nz, gamma, dx, dz)
    if nargin < 3 || isempty(gamma), gamma = 1.0; end

    x = ((1:nx) - (nx+1)/2) * dx;
    z = ((1:nz) - (nz+1)/2) * dz;
    [X, Z] = meshgrid(x, z);

    r = sqrt(X.^2 + Z.^2);

    Rx = max(abs(x));
    Rz = max(abs(z));
    R  = gamma * min(Rx, Rz);

    W = zeros(nz, nx);
    mask = (r <= R);
    W(mask) = 0.5 * (1 + cos(pi * r(mask) / R));
end

function Y = fftconv2_same(A, B)
    [Ma, Na] = size(A);
    [Mb, Nb] = size(B);

    M = Ma + Mb - 1;
    N = Na + Nb - 1;

    FA = fft2(A, M, N);
    FB = fft2(B, M, N);

    C = real(ifft2(FA .* FB));

    row0 = floor((Mb-1)/2) + 1;
    col0 = floor((Nb-1)/2) + 1;

    Y = C(row0:row0+Ma-1, col0:col0+Na-1);
end

function Hd = build_donut_mask_2D(KR, f0, csMin, csMax, taperRel)
% Build 2D donut mask in k-space.
%
% Passband:
%   kLow  = 2*pi*f0/csMax
%   kHigh = 2*pi*f0/csMin
%
% If taperRel > 0, use raised-cosine ramps at both edges.

    kLow  = 2*pi*f0 / csMax;
    kHigh = 2*pi*f0 / csMin;

    Hd = double((KR >= kLow) & (KR <= kHigh));

    if nargin >= 5 && ~isempty(taperRel) && taperRel > 0
        bw = taperRel * (kHigh - kLow);

        k1 = kLow;
        k2 = kLow + bw;
        k3 = kHigh - bw;
        k4 = kHigh;

        Hd = zeros(size(KR));

        % Full passband
        Hd(KR >= k2 & KR <= k3) = 1;

        % Rising taper
        idx = (KR >= k1 & KR < k2);
        x = (KR(idx) - k1) / max(k2 - k1, eps);
        Hd(idx) = 0.5 * (1 - cos(pi*x));

        % Falling taper
        idx = (KR > k3 & KR <= k4);
        x = (KR(idx) - k3) / max(k4 - k3, eps);
        Hd(idx) = 0.5 * (1 + cos(pi*x));
    end
end