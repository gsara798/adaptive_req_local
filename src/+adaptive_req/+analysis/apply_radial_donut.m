function [mapping, info] = apply_radial_donut(req_curve, varargin)
%APPLY_RADIAL_DONUT Restrict an existing REQ radial spectrum in k-space.
%
% [MAPPING, INFO] = adaptive_req.analysis.apply_radial_donut(CURVE, ...)
% applies a physical speed range, a prior-relative k range, or their
% intersection to CURVE.Srad. The returned MAPPING can be passed directly
% to quantile_to_k or quantile_to_cs. No truth or error variable is used.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.apply_radial_donut';
addRequired(p, 'req_curve', @isstruct);
addParameter(p, 'F0', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ApplyPhysical', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PhysicalSpeedRange', [0.5 10], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
addParameter(p, 'ApplyPrior', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PriorK', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'RelativeRange', [0.7 1.3], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
addParameter(p, 'TaperRelative', 0.05, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
parse(p, req_curve, varargin{:});

required = {'k_cent', 'Srad'};
missing = required(~isfield(req_curve, required));
assert(isempty(missing), 'REQ curve missing fields: %s', ...
    strjoin(string(missing), ', '));

k = double(req_curve.k_cent(:));
s = max(double(req_curve.Srad(:)), 0);
s(~isfinite(s)) = 0;
assert(numel(k) == numel(s), 'k_cent and Srad must have equal length.');

k_low = -Inf;
k_high = Inf;
if logical(p.Results.ApplyPhysical)
    f0 = p.Results.F0;
    speeds = sort(double(p.Results.PhysicalSpeedRange(:)));
    assert(isfinite(f0) && f0 > 0, ...
        'F0 must be positive when ApplyPhysical is true.');
    k_low = max(k_low, 2*pi*f0/speeds(2));
    k_high = min(k_high, 2*pi*f0/speeds(1));
end
if logical(p.Results.ApplyPrior)
    prior_k = p.Results.PriorK;
    rel = sort(double(p.Results.RelativeRange(:)));
    assert(isfinite(prior_k) && prior_k > 0, ...
        'PriorK must be positive when ApplyPrior is true.');
    k_low = max(k_low, rel(1)*prior_k);
    k_high = min(k_high, rel(2)*prior_k);
end
assert(k_low < k_high, 'The requested donut intersection is empty.');

weight = radial_taper(k, k_low, k_high, p.Results.TaperRelative);
s_filtered = s .* weight.^2;
if sum(s_filtered) <= eps
    s_filtered = s;
    weight = ones(size(s));
    used_fallback = true;
else
    used_fallback = false;
end

ecum = cumsum(s_filtered);
if ecum(end) > 0
    ecum = ecum/ecum(end);
else
    ecum(:) = NaN;
end

mapping = req_curve;
mapping.Srad = s_filtered;
mapping.Srad_norm = s_filtered/max(s_filtered + eps);
mapping.Ecum = ecum;

info = struct('k_low', k_low, 'k_high', k_high, ...
    'weight', weight, 'retained_energy_fraction', ...
    sum(s_filtered)/(sum(s) + eps), 'used_fallback', used_fallback);
end

function w = radial_taper(k, low, high, taper_relative)
w = ones(size(k));
if isfinite(low), w(k < low) = 0; end
if isfinite(high), w(k > high) = 0; end
if taper_relative <= 0, return; end

if isfinite(low) && low > 0
    width = taper_relative*low;
    idx = k >= low & k < low + width;
    w(idx) = sin(0.5*pi*(k(idx)-low)/max(width, eps));
end
if isfinite(high) && high > 0
    width = taper_relative*high;
    idx = k <= high & k > high - width;
    w(idx) = min(w(idx), ...
        sin(0.5*pi*(high-k(idx))/max(width, eps)));
end
end
