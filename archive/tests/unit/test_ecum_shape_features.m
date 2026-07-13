function tests = test_ecum_shape_features
%TEST_ECUM_SHAPE_FEATURES Verify cumulative-curve descriptors.

tests = functiontests(localfunctions);

end

function testSharpCurveHasHigherPeakConcentration(testCase)

k = linspace(1, 100, 100).';

sharp = struct();
sharp.k_cent = k;
sharp.Ecum = 1 ./ (1 + exp(-(k - 50) / 1.5));

diffuse = struct();
diffuse.k_cent = k;
diffuse.Ecum = 1 ./ (1 + exp(-(k - 50) / 12));

feat_sharp = adaptive_req.quantile.extract_ecum_shape_features(sharp);
feat_diffuse = adaptive_req.quantile.extract_ecum_shape_features(diffuse);

verifyGreaterThan(testCase, feat_sharp.ecum_slope_peak_to_mean, ...
    feat_diffuse.ecum_slope_peak_to_mean);
verifyLessThan(testCase, feat_sharp.ecum_width_80_rel, ...
    feat_diffuse.ecum_width_80_rel);

end

function testFeaturesAreScaleInvariant(testCase)

k = linspace(1, 100, 100).';
curve_a = struct('k_cent', k, 'Ecum', (k / max(k)).^2);
curve_b = struct('k_cent', 3 * k, 'Ecum', curve_a.Ecum);

feat_a = adaptive_req.quantile.extract_ecum_shape_features(curve_a);
feat_b = adaptive_req.quantile.extract_ecum_shape_features(curve_b);

verifyEqual(testCase, feat_b.ecum_width_80_rel, ...
    feat_a.ecum_width_80_rel, 'AbsTol', 1e-12);
verifyEqual(testCase, feat_b.ecum_auc_norm, ...
    feat_a.ecum_auc_norm, 'AbsTol', 1e-12);

end
