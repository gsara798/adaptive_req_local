function out = req_extract_patch_features_from_curve(req_curve, opt)
%REQ_EXTRACT_PATCH_FEATURES_FROM_CURVE Extract local features from REQ output.
%
% This avoids recomputing the windowed FFT after REQ has already built the
% smoothed 2D spectrum for the same patch.

if isfield(req_curve, 'Ssm')
    S2D = double(req_curve.Ssm);
elseif isfield(req_curve, 'Sxz')
    S2D = double(req_curve.Sxz);
else
    error('req_curve must contain Ssm or Sxz.');
end

S2D = max(S2D, 0);
S2D(~isfinite(S2D)) = 0;
S2D = S2D ./ (sum(S2D(:)) + eps);

kx = double(req_curve.kx(:).');
kz = double(req_curve.kz(:));

if isfield(req_curve, 'KR') && isfield(req_curve, 'TH') && ...
        isequal(size(req_curve.KR), size(S2D)) && ...
        isequal(size(req_curve.TH), size(S2D))
    KR = req_curve.KR;
    TH = req_curve.TH;
    KX = [];
    KZ = [];
else
    [KX, KZ] = meshgrid(kx, kz);
    KR = sqrt(KX.^2 + KZ.^2);
    TH = mod(atan2(KZ, KX), 2*pi);
end

radial = adaptive_req.features.extract_radial_features(S2D, KR, opt);
angular = adaptive_req.features.extract_angular_features( ...
    S2D, TH, KR, radial, opt);

out = struct();
out.scalar = adaptive_req.features.combine_feature_structs(radial, angular);
out.radial = radial;
out.angular = angular;
out.diag = struct( ...
    'S2D', S2D, ...
    'kx', kx, ...
    'kz', kz, ...
    'KX', KX, ...
    'KZ', KZ, ...
    'KR', KR, ...
    'TH', TH);

end
