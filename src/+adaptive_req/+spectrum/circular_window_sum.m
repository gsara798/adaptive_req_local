function y = circular_window_sum(x, window_bins)

x = x(:).';
N = numel(x);

left_bins  = floor((window_bins - 1) / 2);
right_bins = ceil((window_bins - 1) / 2);

xwrap = [x, x, x];
mid = N + (1:N);

y = zeros(1, N);
for i = 1:N
    idx = (mid(i) - left_bins):(mid(i) + right_bins);
    y(i) = sum(xwrap(idx));
end

end