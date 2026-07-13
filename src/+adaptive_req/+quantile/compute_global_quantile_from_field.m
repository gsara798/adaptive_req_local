function [q_val, req_curve, feat_global] = compute_global_quantile_from_field( ...
    Uxz, cfg, req, feat_cfg)
%COMPUTE_GLOBAL_QUANTILE_FROM_FIELD Compute one REQ curve from the full field.
%
% This is the global counterpart to compute_quantile_from_patch. It uses the
% full simulated field as the analysis window and returns:
%   q_val      : Ecum(k0_true) from the global spectrum
%   req_curve  : lightweight/global REQ diagnostics
%   feat_global: scalar spectral features from the global 2D spectrum

req_global = req;
req_global.W2 = adaptive_req.spectrum.hann2_circular_shrink( ...
    size(Uxz, 2), size(Uxz, 1), feat_cfg.gamma_win, cfg.dx, cfg.dz);
req_global.PAD = round(feat_cfg.pad_factor * [size(Uxz, 1), size(Uxz, 2)]);
req_global.win_size = size(Uxz, 2);
req_global.half_win = floor(size(Uxz, 2) / 2);
req_global.pad_factor = feat_cfg.pad_factor;

[q_val, req_curve] = adaptive_req.quantile.compute_quantile_from_patch( ...
    Uxz, cfg, req_global);

feat_pack = adaptive_req.features.req_extract_patch_features( ...
    Uxz, cfg.dx, cfg.dz, cfg.f0, cfg.cs_bg, feat_cfg);

feat_global = feat_pack.scalar;

end
