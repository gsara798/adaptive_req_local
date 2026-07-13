function spec = compute_local_spectrum(V_raw, dx, dz, f0, cs_ref, opt)

k0 = 2*pi*f0 / cs_ref;

V0 = V_raw - mean(V_raw(:));

W = adaptive_req.spectrum.hann2_circular_shrink(size(V0,2), size(V0,1), opt.gamma_win, dx, dz);
Vw = V0 .* W;
Vw = Vw / (std(Vw(:)) + eps);

pad_xy = round(opt.pad_factor * [size(Vw,1), size(Vw,2)]);
Vpad = padarray(Vw, pad_xy, 0, 'both');

[Nz2, Nx2] = size(Vpad);

dkx = 2*pi / (Nx2 * dx);
dkz = 2*pi / (Nz2 * dz);

kx = (-floor(Nx2/2):ceil(Nx2/2)-1) * dkx;
kz = (-floor(Nz2/2):ceil(Nz2/2)-1) * dkz;

[KX, KZ] = meshgrid(kx, kz);
KR = sqrt(KX.^2 + KZ.^2);
TH = mod(atan2(KZ, KX), 2*pi);

Uf = fftshift(fft2(Vpad));
S2D = abs(Uf).^2;

icx = floor(size(S2D,2)/2) + 1;
icz = floor(size(S2D,1)/2) + 1;
S2D(icz, icx) = 0;

if opt.smooth_sigma_2d > 0
    S2D = imgaussfilt(S2D, opt.smooth_sigma_2d);
end

S2D = max(S2D, 0);
S2D = S2D / (sum(S2D(:)) + eps);

spec = struct();
spec.V_raw = V_raw;
spec.Vw    = Vw;
spec.Vpad  = Vpad;
spec.S2D   = S2D;

spec.kx    = kx;
spec.kz    = kz;
spec.k0    = k0;

spec.KX    = KX;
spec.KZ    = KZ;
spec.KR    = KR;
spec.TH    = TH;

end