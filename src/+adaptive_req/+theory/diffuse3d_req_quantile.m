function out = diffuse3d_req_quantile(dx, dz, f0, cs_guess, varargin)
%DIFFUSE3D_REQ_QUANTILE Discrete REQ q for projected diffuse 3D theory.
%
% This is a lightweight wrapper around the project-local discrete theory
% implementation adaptive_req.theory.q_theory_REQ_discrete_shearUZ.

p = inputParser;
p.FunctionName = 'adaptive_req.theory.diffuse3d_req_quantile';

addRequired(p, 'dx', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'dz', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'f0', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'cs_guess', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'M', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Gamma', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PadFactor', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Nbins', 'auto', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'SmoothSigma', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'TheoryMode', 'S2D', @(x) ischar(x) || isstring(x));

parse(p, dx, dz, f0, cs_guess, varargin{:});
P = p.Results;

out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    dx, dz, f0, cs_guess, ...
    'M', P.M, ...
    'Gamma', P.Gamma, ...
    'PadFactor', P.PadFactor, ...
    'Nbins', P.Nbins, ...
    'SmoothSigma', P.SmoothSigma, ...
    'TheoryMode', P.TheoryMode, ...
    'FieldType', 'Diffuse3D', ...
    'Plot', false);
out.source = "adaptive_req.theory.q_theory_REQ_discrete_shearUZ";

end
