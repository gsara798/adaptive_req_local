function kq = quantile_to_k(req_curve, q)

k = double(req_curve.k_cent(:));
E = double(req_curve.Ecum(:));

valid = isfinite(k) & isfinite(E);
k = k(valid);
E = E(valid);

if numel(k) < 2 || numel(E) < 2 || ~isfinite(q)
    kq = NaN;
    return;
end

[Euniq, ia] = unique(E, 'stable');
kuniq = k(ia);

if numel(Euniq) < 2
    kq = NaN;
    return;
end

q_clamped = min(max(q, min(Euniq)), max(Euniq));

kq = interp1(Euniq, kuniq, q_clamped, 'linear', 'extrap');

end
