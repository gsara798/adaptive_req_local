function angular = extract_angular_features(S2D, TH, KR, radial, opt)

k_center = radial.k_q50;
iqr_k    = radial.iqr_k;

half_band = max(opt.angular_band_rel_iqr * iqr_k, opt.angular_band_min_abs + eps);
k_min = max(0, k_center - half_band);
k_max = k_center + half_band;

th_edges = linspace(0, 2*pi, opt.Nang + 1);
th_cent  = 0.5 * (th_edges(1:end-1) + th_edges(2:end));

Pang = zeros(1, opt.Nang);

mask_ring = (KR >= k_min) & (KR <= k_max);

for i = 1:opt.Nang
    mask_th = (TH >= th_edges(i)) & (TH < th_edges(i+1));
    mask = mask_ring & mask_th;
    Pang(i) = sum(S2D(mask), 'omitnan');
end

Pang = max(Pang, 0);

if sum(Pang) <= 0
    Pang = nan(1, opt.Nang);
    angular = struct();
    angular.ang_entropy = NaN;
    angular.circ_var = NaN;
    angular.dom_dir_frac = NaN;
    angular.window_max_frac = NaN;
    angular.window_cf = NaN;
    shape = adaptive_req.features.extract_angular_shape_features( ...
        Pang, th_cent, nan(size(Pang)), opt);
    angular = merge_structs(angular, shape);
    angular.theta_pref_deg = NaN;
    angular.theta_pref_rad = NaN;
    angular.k_ang_min = k_min;
    angular.k_ang_max = k_max;
    angular.k_ang_center = k_center;
    angular.window_uniform_baseline = min(1, opt.window_deg / 360);
    angular.th_cent = th_cent;
    angular.th_deg = mod(rad2deg(th_cent), 360);
    angular.Pang = Pang;
    angular.win_sum = nan(size(Pang));
    return;
end

Pang = Pang / (sum(Pang) + eps);

if opt.angular_smooth_bins > 1
    Pang = movmean(Pang, opt.angular_smooth_bins);
    Pang = Pang / (sum(Pang) + eps);
end

th_deg = mod(rad2deg(th_cent), 360);
[th_deg, ord] = sort(th_deg);
th_cent = th_cent(ord);
Pang = Pang(ord);

[~, idx_max] = max(Pang);
theta_pref_deg = th_deg(idx_max);
theta_pref_rad = th_cent(idx_max);

ang_entropy = -sum(Pang .* log(Pang + eps)) / log(numel(Pang));

R = abs(sum(Pang .* exp(1i * th_cent)));
circ_var = 1 - R;

dom_dir_frac = max(Pang);

window_bins = max(1, round(opt.window_deg / (360 / opt.Nang)));
win_sum = adaptive_req.spectrum.circular_window_sum(Pang, window_bins);

window_max_frac = max(win_sum);
window_uniform_baseline = min(1, opt.window_deg / 360);
window_cf = window_max_frac / (window_uniform_baseline + eps);
shape = adaptive_req.features.extract_angular_shape_features( ...
    Pang, th_cent, win_sum, opt);

angular = struct();

% Features for regression
angular.ang_entropy     = ang_entropy;
angular.circ_var        = circ_var;
angular.dom_dir_frac    = dom_dir_frac;
angular.window_max_frac = window_max_frac;
angular.window_cf       = window_cf;
angular = merge_structs(angular, shape);

% Diagnostics
angular.theta_pref_deg = theta_pref_deg;
angular.theta_pref_rad = theta_pref_rad;
angular.k_ang_min      = k_min;
angular.k_ang_max      = k_max;
angular.k_ang_center   = k_center;
angular.window_uniform_baseline = window_uniform_baseline;

angular.th_cent = th_cent;
angular.th_deg  = th_deg;
angular.Pang    = Pang;
angular.win_sum = win_sum;

end

function out = merge_structs(out, extra)

names = fieldnames(extra);

for i = 1:numel(names)
    out.(names{i}) = extra.(names{i});
end

end
