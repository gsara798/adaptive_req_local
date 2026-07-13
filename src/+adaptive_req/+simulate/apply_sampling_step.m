function cfg_i = apply_sampling_step(cfg, sampling_mode, step_value, varargin)
%APPLY_SAMPLING_STEP Apply one aperture value to a simulation config.
%
% Usage:
%   cfg_i = adaptive_req.simulate.apply_sampling_step(cfg, 'cone', 180);
%
% For cone mode:
%   step_value is full aperture in degrees.
%   cfg_i.ConeHalfAngleDeg = step_value/2.
%
% For band mode:
%   step_value is band half-width in degrees.
%
% For ranges mode:
%   step_value is angular range in degrees.

p = inputParser;

addRequired(p, 'cfg', @isstruct);
addRequired(p, 'sampling_mode', @(x) ischar(x) || isstring(x));
addRequired(p, 'step_value', @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'RangeCenterDeg', 0, @(x) isnumeric(x) && isscalar(x));

parse(p, cfg, sampling_mode, step_value, varargin{:});

cfg_i = cfg;
sampling_mode = lower(char(p.Results.sampling_mode));

switch sampling_mode

    case 'cone'
        cfg_i.SourceSampling = 'cone';
        cfg_i.ConeHalfAngleDeg = step_value / 2;

        if ~isfield(cfg_i, 'ConeAxis') || isempty(cfg_i.ConeAxis)
            cfg_i.ConeAxis = [1 0 0];
        end

        cfg_i.ConeHalfAngleDeg = max(0, min(180, cfg_i.ConeHalfAngleDeg));

    case 'band'
        cfg_i.SourceSampling = 'band';
        cfg_i.BandHalfWidthDeg = step_value;

        if ~isfield(cfg_i, 'BandAxis') || isempty(cfg_i.BandAxis)
            cfg_i.BandAxis = [0 1 0];
        end

        cfg_i.BandHalfWidthDeg = max(0, min(90, cfg_i.BandHalfWidthDeg));

    case 'ranges'
        cfg_i.SourceSampling = 'ranges';

        range_center = deg2rad(p.Results.RangeCenterDeg);
        half_range = deg2rad(step_value / 2);

        if isfield(cfg_i, 'Is2D') && logical(cfg_i.Is2D)
            cfg_i.AngleRange2D = range_center + [-half_range, half_range];
        else
            cfg_i.PhiRange = [0, 2*pi];
            cfg_i.ThetaRange = [0, min(pi, half_range)];
        end

    otherwise
        error('Unknown sampling_mode: %s', sampling_mode);
end

end