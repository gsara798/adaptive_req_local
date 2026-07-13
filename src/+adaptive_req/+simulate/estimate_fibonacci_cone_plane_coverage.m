function T = estimate_fibonacci_cone_plane_coverage(varargin)
%ESTIMATE_FIBONACCI_CONE_PLANE_COVERAGE Fast analytic plane-coverage audit.

p = inputParser;
p.FunctionName = 'adaptive_req.simulate.estimate_fibonacci_cone_plane_coverage';

addParameter(p, 'Nwaves', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'ConeAxis', [1 0 0], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'FullApertureDeg', linspace(0, 360, 10), @isnumeric);
addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'ForceInPlaneWave', false, @(x) islogical(x) || isnumeric(x));

parse(p, varargin{:});

N = round(p.Results.Nwaves);
axis_hat = double(p.Results.ConeAxis(:));
axis_hat = axis_hat / norm(axis_hat);
full_apertures = double(p.Results.FullApertureDeg(:));
tol = p.Results.Tolerance;
force_in_plane = logical(p.Results.ForceInPlaneWave);

rows(numel(full_apertures)) = struct();

for i = 1:numel(full_apertures)

    half_angle_deg = full_apertures(i) / 2;
    [ux, uy, uz] = fibonacci_cap_local(axis_hat, deg2rad(half_angle_deg), N);

    if force_in_plane
        target = [axis_hat(1); 0; axis_hat(3)];

        if norm(target) > eps
            target = target / norm(target);
            [~, idx] = min(abs(uy));
            ux(idx) = target(1);
            uy(idx) = 0;
            uz(idx) = target(3);
        end
    end

    summary = adaptive_req.simulate.summarize_wave_direction_plane_coverage( ...
        struct('ux', ux, 'uy', uy, 'uz', uz), ...
        'Tolerance', tol);

    rows(i).step_idx = i;
    rows(i).full_aperture_deg = full_apertures(i);
    rows(i).cone_half_angle_deg = half_angle_deg;
    rows(i).Omega_sr = 2 * pi * (1 - cos(deg2rad(half_angle_deg)));
    rows(i).Nwaves = N;
    rows(i).tolerance = tol;
    rows(i).force_in_plane_wave = force_in_plane;
    rows(i).min_abs_uy = summary.min_abs_uy;
    rows(i).median_abs_uy = summary.median_abs_uy;
    rows(i).p01_abs_uy = summary.p01_abs_uy;
    rows(i).count_in_plane = summary.count_in_plane;
    rows(i).has_in_plane_wave = summary.has_in_plane_wave;

end

T = struct2table(rows);

end

function [ux, uy, uz] = fibonacci_cap_local(axis_hat, half_ang, N)

if abs(dot(axis_hat, [0; 0; 1])) < 0.99
    tmp = [0; 0; 1];
else
    tmp = [0; 1; 0];
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

ux = dirs(1, :).';
uy = dirs(2, :).';
uz = dirs(3, :).';

end
