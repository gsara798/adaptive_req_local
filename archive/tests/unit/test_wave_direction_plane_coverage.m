function tests = test_wave_direction_plane_coverage
%TEST_WAVE_DIRECTION_PLANE_COVERAGE Verify in-plane direction diagnostics.

tests = functiontests(localfunctions);

end

function testDetectsInPlaneWave(testCase)

dirs = struct();
dirs.uy = [0.3, 0, -0.2];

summary = adaptive_req.simulate.summarize_wave_direction_plane_coverage( ...
    dirs, 'Tolerance', 1e-9);

verifyTrue(testCase, summary.has_in_plane_wave);
verifyEqual(testCase, summary.count_in_plane, 1);
verifyEqual(testCase, summary.min_abs_uy, 0);

end

function testRejectsNearButNotExactWave(testCase)

dirs = struct();
dirs.uy = [1e-4, -2e-4, 3e-4];

summary = adaptive_req.simulate.summarize_wave_direction_plane_coverage( ...
    dirs, 'Tolerance', 1e-6);

verifyFalse(testCase, summary.has_in_plane_wave);
verifyEqual(testCase, summary.min_abs_uy, 1e-4, 'AbsTol', 1e-12);

end
