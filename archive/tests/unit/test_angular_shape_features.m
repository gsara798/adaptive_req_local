function tests = test_angular_shape_features
%TEST_ANGULAR_SHAPE_FEATURES Verify rotation-invariant angular descriptors.

tests = functiontests(localfunctions);

end

function testUniformSpectrumHasLowMoments(testCase)

[theta, opt] = make_inputs();
Pang = ones(size(theta));
Pang = Pang / sum(Pang);
win_sum = adaptive_req.spectrum.circular_window_sum(Pang, 4);

shape = adaptive_req.features.extract_angular_shape_features( ...
    Pang, theta, win_sum, opt);

verifyLessThan(testCase, shape.ang_moment_1, 1e-12);
verifyLessThan(testCase, shape.ang_moment_2, 1e-12);
verifyLessThan(testCase, shape.ang_moment_4, 1e-12);

end

function testOpposingLobesAreCapturedBySecondMoment(testCase)

[theta, opt] = make_inputs();
Pang = circular_lobe(theta, 0, 0.12) + circular_lobe(theta, pi, 0.12);
Pang = Pang / sum(Pang);
win_sum = adaptive_req.spectrum.circular_window_sum(Pang, 4);

shape = adaptive_req.features.extract_angular_shape_features( ...
    Pang, theta, win_sum, opt);

verifyLessThan(testCase, shape.ang_moment_1, 0.05);
verifyGreaterThan(testCase, shape.ang_moment_2, 0.90);
verifyEqual(testCase, shape.ang_peak_separation_deg, 180, 'AbsTol', 5);

end

function testFeaturesAreRotationInvariant(testCase)

[theta, opt] = make_inputs();
Pang = circular_lobe(theta, 0.4, 0.18) + ...
    0.6 * circular_lobe(theta, 2.2, 0.22);
Pang = Pang / sum(Pang);

shift_bins = 11;
Pang_shifted = circshift(Pang, shift_bins);

shape_a = adaptive_req.features.extract_angular_shape_features( ...
    Pang, theta, adaptive_req.spectrum.circular_window_sum(Pang, 4), opt);
shape_b = adaptive_req.features.extract_angular_shape_features( ...
    Pang_shifted, theta, ...
    adaptive_req.spectrum.circular_window_sum(Pang_shifted, 4), opt);

names = fieldnames(shape_a);

for i = 1:numel(names)
    verifyEqual(testCase, shape_b.(names{i}), shape_a.(names{i}), ...
        'AbsTol', 1e-10);
end

end

function [theta, opt] = make_inputs()

n = 72;
theta = (0:n-1) * 2 * pi / n;
opt = struct();
opt.window_deg = 20;
opt.angular_peak_min_rel = 0.35;
opt.angular_peak_suppression_deg = 30;
opt.angular_peak_max_count = 8;

end

function y = circular_lobe(theta, center, sigma)

delta = atan2(sin(theta - center), cos(theta - center));
y = exp(-0.5 * (delta / sigma).^2);

end
