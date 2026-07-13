function cs = quantile_to_cs(req_curve, q, f0)

kq = adaptive_req.quantile.quantile_to_k(req_curve, q);

if ~isfinite(kq) || kq <= 0
    cs = NaN;
else
    cs = 2*pi*f0 ./ kq;
end

end