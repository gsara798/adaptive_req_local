function tests = test_req_mapping
%TEST_REQ_MAPPING Verify that compact REQ storage preserves q-to-SWS results.

tests = functiontests(localfunctions);

end

function testMappingPreservesQuantileConversion(testCase)

curve = struct();
curve.Ecum = [0; 0.1; 0.4; 0.75; 1.0];
curve.k_cent = [100; 200; 300; 400; 500];
curve.f0 = 600;
curve.Nbins_effective = 5;

mapping = adaptive_req.quantile.make_req_mapping(curve);
q = 0.6;

k_full = adaptive_req.quantile.quantile_to_k(curve, q);
k_light = adaptive_req.quantile.quantile_to_k(mapping, q);

cs_full = adaptive_req.quantile.quantile_to_cs(curve, q, curve.f0);
cs_light = adaptive_req.quantile.quantile_to_cs(mapping, q, curve.f0);

verifyEqual(testCase, k_light, k_full, 'AbsTol', 1e-4);
verifyEqual(testCase, cs_light, cs_full, 'AbsTol', 1e-6);

verifyEqual(testCase, class(mapping.Ecum), 'single');
verifyEqual(testCase, class(mapping.k_cent), 'single');
verifyFalse(testCase, isfield(mapping, 'Sxz'));
verifyFalse(testCase, isfield(mapping, 'Ssm'));

end
