function out = req_extract_patch_features(V_raw, dx, dz, f0, cs_ref, opt)

spec = adaptive_req.spectrum.compute_local_spectrum(V_raw, dx, dz, f0, cs_ref, opt);

radial  = adaptive_req.features.extract_radial_features(spec.S2D, spec.KR, opt);
angular = adaptive_req.features.extract_angular_features(spec.S2D, spec.TH, spec.KR, radial, opt);

out = struct();
out.scalar  = adaptive_req.features.combine_feature_structs(radial, angular);
out.radial  = radial;
out.angular = angular;
out.diag    = spec;

end