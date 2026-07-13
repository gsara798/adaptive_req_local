function feat = extract_ecum_shape_features(req_curve)
%EXTRACT_ECUM_SHAPE_FEATURES Shape descriptors from the REQ cumulative curve.
%
% The descriptors are operational: they use only Ecum(k) and k_cent. They do
% not use true SWS, q_true, or a guessed wave speed.

feat = empty_features();

if isempty(req_curve) || ~isstruct(req_curve) || ...
        ~isfield(req_curve, 'Ecum') || ~isfield(req_curve, 'k_cent')
    return;
end

E = double(req_curve.Ecum(:));
k = double(req_curve.k_cent(:));

valid = isfinite(E) & isfinite(k) & k > 0;
E = E(valid);
k = k(valid);

if numel(E) < 3
    return;
end

[k, ord] = sort(k, 'ascend');
E = E(ord);

E = cummax(E);

if max(E) <= min(E)
    return;
end

E = (E - min(E)) ./ (max(E) - min(E));

kmax = max(k);
if ~isfinite(kmax) || kmax <= 0
    return;
end

k10 = k_at_q(k, E, 0.10);
k25 = k_at_q(k, E, 0.25);
k50 = k_at_q(k, E, 0.50);
k75 = k_at_q(k, E, 0.75);
k90 = k_at_q(k, E, 0.90);

k50_safe = max(k50, eps);
width_50 = max(k75 - k25, eps);
width_80 = max(k90 - k10, eps);
lower_width_50 = max(k50 - k25, eps);
upper_width_50 = max(k75 - k50, eps);

feat.ecum_k10_norm_k50 = k10 / k50_safe;
feat.ecum_k25_norm_k50 = k25 / k50_safe;
feat.ecum_k75_norm_k50 = k75 / k50_safe;
feat.ecum_k90_norm_k50 = k90 / k50_safe;

feat.ecum_width_50_rel = (k75 - k25) / k50_safe;
feat.ecum_width_80_rel = (k90 - k10) / k50_safe;
feat.ecum_lower_tail_rel = (k50 - k10) / k50_safe;
feat.ecum_upper_tail_rel = (k90 - k50) / k50_safe;
feat.ecum_asymmetry_10_90 = (k90 + k10 - 2 * k50) / width_80;
feat.ecum_asymmetry_25_75 = (k75 + k25 - 2 * k50) / width_50;
feat.ecum_width_ratio_80_50 = width_80 / width_50;
feat.ecum_lower_upper_width_ratio = lower_width_50 / upper_width_50;

k_norm = k ./ kmax;
feat.ecum_auc_norm = trapz(k_norm, E) / max(range(k_norm), eps);

dE = diff(E);
dk = diff(k_norm);
valid_step = isfinite(dE) & isfinite(dk) & dk > 0;
dE = max(dE(valid_step), 0);
dk = dk(valid_step);

if isempty(dE) || sum(dE) <= 0
    return;
end

p = dE ./ sum(dE);
slope = dE ./ dk;
slope = slope(isfinite(slope));
ds_centers = 0.5 * (k_norm(1:end-1) + k_norm(2:end));
ds_centers = ds_centers(valid_step);

feat.ecum_increment_entropy = -sum(p .* log(p + eps)) / log(numel(p));
feat.ecum_increment_peak_frac = max(p);
feat.ecum_increment_gini = 1 - sum(p.^2);

if numel(ds_centers) == numel(p)
    mu_k = sum(p .* ds_centers);
    sigma_k = sqrt(sum(p .* (ds_centers - mu_k).^2));
    sigma_safe = max(sigma_k, eps);
    [~, idx_peak] = max(p);

    feat.srad_proxy_centroid_k_norm = mu_k;
    feat.srad_proxy_std_k_norm = sigma_k;
    feat.srad_proxy_skewness = sum(p .* (ds_centers - mu_k).^3) / ...
        sigma_safe^3;
    feat.srad_proxy_kurtosis = sum(p .* (ds_centers - mu_k).^4) / ...
        sigma_safe^4;
    feat.srad_proxy_peak_k_norm = ds_centers(idx_peak);
    feat.srad_proxy_peak_to_centroid = ds_centers(idx_peak) / ...
        max(mu_k, eps);
    feat.srad_proxy_low_side_frac = sum(p(ds_centers < mu_k));
    feat.srad_proxy_high_side_frac = sum(p(ds_centers >= mu_k));
end

if ~isempty(slope)
    feat.ecum_slope_max = max(slope);
    feat.ecum_slope_peak_to_mean = max(slope) / (mean(slope, 'omitnan') + eps);
    feat.ecum_slope_iqr_to_median = iqr(slope) / (median(slope, 'omitnan') + eps);
end

end

function feat = empty_features()

feat = struct();
feat.ecum_k10_norm_k50 = NaN;
feat.ecum_k25_norm_k50 = NaN;
feat.ecum_k75_norm_k50 = NaN;
feat.ecum_k90_norm_k50 = NaN;
feat.ecum_width_50_rel = NaN;
feat.ecum_width_80_rel = NaN;
feat.ecum_lower_tail_rel = NaN;
feat.ecum_upper_tail_rel = NaN;
feat.ecum_asymmetry_10_90 = NaN;
feat.ecum_asymmetry_25_75 = NaN;
feat.ecum_width_ratio_80_50 = NaN;
feat.ecum_lower_upper_width_ratio = NaN;
feat.ecum_auc_norm = NaN;
feat.ecum_increment_entropy = NaN;
feat.ecum_increment_peak_frac = NaN;
feat.ecum_increment_gini = NaN;
feat.ecum_slope_max = NaN;
feat.ecum_slope_peak_to_mean = NaN;
feat.ecum_slope_iqr_to_median = NaN;
feat.srad_proxy_centroid_k_norm = NaN;
feat.srad_proxy_std_k_norm = NaN;
feat.srad_proxy_skewness = NaN;
feat.srad_proxy_kurtosis = NaN;
feat.srad_proxy_peak_k_norm = NaN;
feat.srad_proxy_peak_to_centroid = NaN;
feat.srad_proxy_low_side_frac = NaN;
feat.srad_proxy_high_side_frac = NaN;

end

function kq = k_at_q(k, E, q)

[Euniq, ia] = unique(E, 'stable');
kuniq = k(ia);

if numel(Euniq) < 2
    kq = NaN;
    return;
end

q = min(max(q, min(Euniq)), max(Euniq));
kq = interp1(Euniq, kuniq, q, 'linear', 'extrap');

end
