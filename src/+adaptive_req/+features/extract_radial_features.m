function radial = extract_radial_features(S2D, KR, opt)

k = KR(:);
w = S2D(:);

valid = isfinite(k) & isfinite(w) & (w > 0);
k = k(valid);
w = w(valid);

w = w / (sum(w) + eps);

[k_sorted, ord] = sort(k);
w_sorted = w(ord);

% combine repeated k values before building the CDF
[k_unique, ~, ic] = unique(k_sorted);
w_unique = accumarray(ic, w_sorted);

F_sorted = cumsum(w_unique);
F_sorted = F_sorted / (F_sorted(end) + eps);

k_q10  = weighted_quantile_sorted(k_unique, F_sorted, 0.10);
k_q25  = weighted_quantile_sorted(k_unique, F_sorted, 0.25);
k_q50  = weighted_quantile_sorted(k_unique, F_sorted, 0.50);
k_q75  = weighted_quantile_sorted(k_unique, F_sorted, 0.75);
k_q90  = weighted_quantile_sorted(k_unique, F_sorted, 0.90);
k_q995 = weighted_quantile_sorted(k_unique, F_sorted, opt.k_plot_qhi);

k50_safe = max(k_q50, eps);
iqr_k = max(k_q75 - k_q25, eps);

lowk_frac_rel = sum(w(k <= opt.lowk_rel_alpha * k50_safe));

mid_half_width = max(opt.midband_rel_iqr * iqr_k, opt.midband_min_abs + eps);
midband_frac_rel = sum(w(abs(k - k_q50) <= mid_half_width));

highk_frac_rel = sum(w(k >= opt.highk_rel_alpha * k50_safe));

k_plot_hi = max(k_q995, k_q90 + 2 * iqr_k);
k_plot_hi = max(k_plot_hi, k_q50 + 3 * iqr_k);

Nrad_plot = opt.Nrad_plot;
k_edges_plot = linspace(0, k_plot_hi, Nrad_plot + 1);
k_cent_plot  = 0.5 * (k_edges_plot(1:end-1) + k_edges_plot(2:end));

Prad_plot = weighted_histogram_1d(k, w, k_edges_plot);
if opt.radial_smooth_bins > 1
    Prad_plot = movmean(Prad_plot, opt.radial_smooth_bins);
end
Prad_plot = Prad_plot / (sum(Prad_plot) + eps);

Frad_plot = zeros(size(k_cent_plot));
for i = 1:numel(k_cent_plot)
    Frad_plot(i) = sum(w(k <= k_cent_plot(i)));
end

k_edges_entropy = [k_edges_plot, inf];
P_entropy = weighted_histogram_1d(k, w, k_edges_entropy);
P_entropy = P_entropy / (sum(P_entropy) + eps);
radial_entropy = -sum(P_entropy .* log(P_entropy + eps)) / log(numel(P_entropy));

width_75_25_rel = (k_q75 - k_q25) / k50_safe;
width_90_50_rel = (k_q90 - k_q50) / k50_safe;
width_90_10_rel = (k_q90 - k_q10) / k50_safe;

radial = struct();

% Features for regression
radial.radial_entropy   = radial_entropy;
radial.width_75_25_rel  = width_75_25_rel;
radial.width_90_50_rel  = width_90_50_rel;
radial.width_90_10_rel  = width_90_10_rel;
radial.lowk_frac_rel    = lowk_frac_rel;
radial.midband_frac_rel = midband_frac_rel;
radial.highk_frac_rel   = highk_frac_rel;

% Diagnostics
radial.k_q10   = k_q10;
radial.k_q25   = k_q25;
radial.k_q50   = k_q50;
radial.k_q75   = k_q75;
radial.k_q90   = k_q90;
radial.k_q995  = k_q995;
radial.iqr_k   = iqr_k;

radial.k_plot_hi   = k_plot_hi;
radial.k_cent_plot = k_cent_plot;
radial.Prad_plot   = Prad_plot;
radial.Frad_plot   = Frad_plot;

end

function xq = weighted_quantile_sorted(x_sorted, F_sorted, q)

q = max(0, min(1, q));

valid = isfinite(x_sorted) & isfinite(F_sorted);
x_sorted = x_sorted(valid);
F_sorted = F_sorted(valid);

if isempty(x_sorted)
    xq = NaN;
    return;
end

Fend = F_sorted(end);
if ~isfinite(Fend) || Fend <= 0
    xq = NaN;
    return;
end

F_sorted = F_sorted / Fend;

idx = find(F_sorted >= q, 1, 'first');

if isempty(idx)
    xq = x_sorted(end);
    return;
end

if idx == 1
    xq = x_sorted(1);
    return;
end

x1 = x_sorted(idx-1);
x2 = x_sorted(idx);
F1 = F_sorted(idx-1);
F2 = F_sorted(idx);

if F2 <= F1
    xq = x2;
else
    t = (q - F1) / (F2 - F1);
    xq = x1 + t * (x2 - x1);
end

end

function h = weighted_histogram_1d(x, w, edges)

bin = discretize(x, edges);
valid = ~isnan(bin);

h = accumarray(bin(valid), w(valid), [numel(edges)-1, 1], @sum, 0);
h = h(:).';

end