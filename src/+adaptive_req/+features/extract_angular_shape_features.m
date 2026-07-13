function shape = extract_angular_shape_features(Pang, th_cent, win_sum, opt)
%EXTRACT_ANGULAR_SHAPE_FEATURES Rotation-invariant angular spectrum features.
%
% These features distinguish uniform, unidirectional, bidirectional, and
% multilobe spectra without using the simulated wave-field label.

Pang = Pang(:).';
th_cent = th_cent(:).';
win_sum = win_sum(:).';

shape = empty_shape_features();

valid = isfinite(Pang) & isfinite(th_cent) & isfinite(win_sum);

if ~all(valid) || isempty(Pang) || sum(Pang) <= 0
    return;
end

Pang = max(Pang, 0);
Pang = Pang / (sum(Pang) + eps);

for order = [1 2 4]
    value = abs(sum(Pang .* exp(1i * order * th_cent)));
    shape.(sprintf('ang_moment_%d', order)) = value;
end

n = numel(win_sum);
window_bins = max(1, round(opt.window_deg / (360 / n)));
suppression_bins = max(window_bins, ...
    round(opt.angular_peak_suppression_deg / (360 / n)));

[peak_values, peak_indices] = select_separated_peaks( ...
    win_sum, suppression_bins, opt.angular_peak_max_count);

if isempty(peak_values)
    return;
end

peak_values = peak_values(:).';
peak_indices = peak_indices(:).';

shape.ang_top1_window_frac = peak_values(1);
shape.ang_top2_window_frac = sum(peak_values(1:min(2, end)));
shape.ang_top3_window_frac = sum(peak_values(1:min(3, end)));

if numel(peak_values) >= 2
    shape.ang_top2_to_top1 = peak_values(2) / (peak_values(1) + eps);
    shape.ang_peak_separation_deg = circular_separation_deg( ...
        th_cent(peak_indices(1)), th_cent(peak_indices(2)));
else
    shape.ang_top2_to_top1 = 0;
    shape.ang_peak_separation_deg = 0;
end

relative_threshold = opt.angular_peak_min_rel * peak_values(1);
shape.ang_peak_count_rel = sum(peak_values >= relative_threshold);

end

function shape = empty_shape_features()

shape = struct();
shape.ang_moment_1 = NaN;
shape.ang_moment_2 = NaN;
shape.ang_moment_4 = NaN;
shape.ang_peak_count_rel = NaN;
shape.ang_top1_window_frac = NaN;
shape.ang_top2_window_frac = NaN;
shape.ang_top3_window_frac = NaN;
shape.ang_top2_to_top1 = NaN;
shape.ang_peak_separation_deg = NaN;

end

function [values, indices] = select_separated_peaks(x, suppression_bins, max_count)

x = x(:).';
n = numel(x);
available = true(1, n);

values = zeros(1, 0);
indices = zeros(1, 0);

for peak_idx = 1:min(max_count, n)

    candidates = x;
    candidates(~available) = -Inf;

    [value, idx] = max(candidates);

    if ~isfinite(value)
        break;
    end

    values(end + 1) = value; %#ok<AGROW>
    indices(end + 1) = idx; %#ok<AGROW>

    offsets = -suppression_bins:suppression_bins;
    suppress_idx = mod(idx - 1 + offsets, n) + 1;
    available(suppress_idx) = false;
end

end

function separation = circular_separation_deg(theta_a, theta_b)

delta = abs(rad2deg(atan2( ...
    sin(theta_a - theta_b), ...
    cos(theta_a - theta_b))));

separation = min(delta, 360 - delta);

end
