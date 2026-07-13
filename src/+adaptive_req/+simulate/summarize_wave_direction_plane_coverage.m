function summary = summarize_wave_direction_plane_coverage(wave_dirs, varargin)
%SUMMARIZE_WAVE_DIRECTION_PLANE_COVERAGE Check if waves include y=0 plane.
%
% A wave is exactly in the imaging plane when its y component is zero.

p = inputParser;
p.FunctionName = 'adaptive_req.simulate.summarize_wave_direction_plane_coverage';

addRequired(p, 'wave_dirs', @isstruct);
addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);

parse(p, wave_dirs, varargin{:});

tol = p.Results.Tolerance;

if ~isfield(wave_dirs, 'uy')
    error('wave_dirs must contain field uy.');
end

uy = double(wave_dirs.uy(:));
uy = uy(isfinite(uy));

if isempty(uy)
    summary = empty_summary(tol);
    return;
end

abs_uy = abs(uy);

summary = struct();
summary.tolerance = tol;
summary.n_waves = numel(uy);
summary.min_abs_uy = min(abs_uy);
summary.median_abs_uy = median(abs_uy);
summary.p01_abs_uy = percentile(abs_uy, 1);
summary.p05_abs_uy = percentile(abs_uy, 5);
summary.count_in_plane = sum(abs_uy <= tol);
summary.has_in_plane_wave = summary.count_in_plane > 0;
summary.nearest_idx = find(abs_uy == summary.min_abs_uy, 1, 'first');

end

function summary = empty_summary(tol)

summary = struct();
summary.tolerance = tol;
summary.n_waves = 0;
summary.min_abs_uy = NaN;
summary.median_abs_uy = NaN;
summary.p01_abs_uy = NaN;
summary.p05_abs_uy = NaN;
summary.count_in_plane = 0;
summary.has_in_plane_wave = false;
summary.nearest_idx = NaN;

end

function value = percentile(x, pct)

x = sort(x(:));

if isempty(x)
    value = NaN;
    return;
end

idx = max(1, min(numel(x), ceil(pct / 100 * numel(x))));
value = x(idx);

end
