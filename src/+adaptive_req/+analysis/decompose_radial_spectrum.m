function D = decompose_radial_spectrum(curve, varargin)
%DECOMPOSE_RADIAL_SPECTRUM Detect two separated radial spectral components.
%   D = decompose_radial_spectrum(CURVE) uses only CURVE.k_cent and
%   CURVE.Srad. Expected material speeds, true SWS, and interface labels are
%   intentionally not inputs. The detector is conservative: two components
%   require two resolved maxima, adequate energy on both sides of the
%   intervening valley, and a sufficiently deep valley.

p = inputParser;
addParameter(p, 'KRange', [0 Inf], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p, 'SmoothBins', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p, 'MinPeakHeightRatio', 0.10, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MinSeparationFraction', 0.12, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MaxValleyRatio', 0.82, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MinComponentWeight', 0.08, @(x)isnumeric(x)&&isscalar(x));
parse(p, varargin{:});
opt = p.Results;

D = empty_result();
if ~isstruct(curve) || ~all(isfield(curve, {'k_cent','Srad'}))
    return;
end
k = double(curve.k_cent(:));
s = max(double(curve.Srad(:)), 0);
valid = isfinite(k) & isfinite(s) & k >= opt.KRange(1) & k <= opt.KRange(2);
k = k(valid); s = s(valid);
if numel(k) < 5 || sum(s) <= 0
    return;
end
[k, order] = sort(k); s = s(order);
s = s / sum(s);
span = max(k) - min(k);
if span <= 0
    return;
end

smooth_bins = max(1, round(opt.SmoothBins));
ss = movmean(s, smooth_bins, 'Endpoints', 'shrink');
loc = find(ss(2:end-1) >= ss(1:end-2) & ss(2:end-1) > ss(3:end)) + 1;
if isempty(loc)
    [~, loc] = max(ss);
end
loc = loc(ss(loc) >= opt.MinPeakHeightRatio * max(ss));

D.valid = true;
D.dominant_k = k(find(ss == max(ss), 1));
D.num_candidate_peaks = numel(loc);
D.radial_width_10_90 = spectral_quantile(k, s, .90) - ...
    spectral_quantile(k, s, .10);
if numel(loc) < 2
    return;
end

best_score = -Inf;
best = [];
for a = 1:(numel(loc)-1)
    for b = (a+1):numel(loc)
        ia = loc(a); ib = loc(b);
        separation = (k(ib)-k(ia)) / max(0.5*(k(ib)+k(ia)), eps);
        if separation < opt.MinSeparationFraction
            continue;
        end
        [valley, rel] = min(ss(ia:ib));
        iv = ia + rel - 1;
        weak_height = min(ss(ia),ss(ib)) / max(max(ss),eps);
        valley_ratio = valley / max(min(ss(ia),ss(ib)),eps);
        low_weight = sum(s(1:iv));
        high_weight = sum(s((iv+1):end));
        balance = 2*min(low_weight,high_weight) / max(low_weight+high_weight,eps);
        score = weak_height * max(1-valley_ratio,0) * balance * separation;
        if score > best_score
            best_score = score;
            best = [ia ib iv separation valley_ratio low_weight high_weight ...
                weak_height balance];
        end
    end
end
if isempty(best)
    return;
end

D.low_k = k(best(1));
D.high_k = k(best(2));
D.valley_k = k(best(3));
D.separation_fraction = best(4);
D.valley_ratio = best(5);
D.low_component_weight = best(6);
D.high_component_weight = best(7);
D.weak_peak_height_ratio = best(8);
D.component_balance = best(9);
D.two_component_score = max(best_score,0);
D.two_component_detected = ...
    D.valley_ratio <= opt.MaxValleyRatio && ...
    min(D.low_component_weight,D.high_component_weight) >= opt.MinComponentWeight;
end

function q = spectral_quantile(k,s,p)
E = cumsum(s) / max(sum(s),eps);
[Eu,ia] = unique(E,'stable');
if numel(Eu) < 2
    q = NaN;
else
    q = interp1(Eu,k(ia),p,'linear','extrap');
end
end

function D = empty_result()
D = struct('valid',false,'two_component_detected',false, ...
    'dominant_k',NaN,'low_k',NaN,'high_k',NaN,'valley_k',NaN, ...
    'num_candidate_peaks',0,'separation_fraction',NaN, ...
    'valley_ratio',NaN,'low_component_weight',NaN, ...
    'high_component_weight',NaN,'weak_peak_height_ratio',NaN, ...
    'component_balance',NaN,'two_component_score',0, ...
    'radial_width_10_90',NaN);
end
