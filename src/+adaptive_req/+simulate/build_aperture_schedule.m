function schedule = build_aperture_schedule(sampling_mode, num_steps, varargin)
%BUILD_APERTURE_SCHEDULE Build angular aperture schedule for simulation sweeps.
%
% Usage:
%   schedule = adaptive_req.simulate.build_aperture_schedule('cone', 10);
%   schedule = adaptive_req.simulate.build_aperture_schedule('band', 10);
%
% For cone mode:
%   schedule.values is the full cone aperture in degrees, from 0 to 360.
%   The corresponding cone half-angle is values/2.
%
% For band mode:
%   schedule.values is the band half-width in degrees, from 0 to 90.
%
% Output:
%   schedule.values
%   schedule.Omega_sr
%   schedule.name
%   schedule.unit
%   schedule.mode

p = inputParser;

addRequired(p, 'sampling_mode', @(x) ischar(x) || isstring(x));
addRequired(p, 'num_steps', @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'ConeApertureDeg', [0 360], ...
    @(x) isnumeric(x) && numel(x) == 2);

addParameter(p, 'BandHalfWidthDeg', [0 90], ...
    @(x) isnumeric(x) && numel(x) == 2);

addParameter(p, 'RangeApertureDeg', [0 360], ...
    @(x) isnumeric(x) && numel(x) == 2);

parse(p, sampling_mode, num_steps, varargin{:});

sampling_mode = lower(char(p.Results.sampling_mode));
num_steps = p.Results.num_steps;

schedule = struct();
schedule.mode = sampling_mode;
schedule.num_steps = num_steps;

switch sampling_mode

    case 'cone'
        vals = linspace(p.Results.ConeApertureDeg(1), ...
                        p.Results.ConeApertureDeg(2), ...
                        num_steps);

        half_angle_rad = deg2rad(vals / 2);

        schedule.values = vals;
        schedule.Omega_sr = 2*pi*(1 - cos(half_angle_rad));
        schedule.name = 'Cone aperture';
        schedule.unit = 'deg';
        schedule.description = 'values are full cone aperture; ConeHalfAngleDeg = values/2';

    case 'band'
        vals = linspace(p.Results.BandHalfWidthDeg(1), ...
                        p.Results.BandHalfWidthDeg(2), ...
                        num_steps);

        delta_rad = deg2rad(vals);

        schedule.values = vals;
        schedule.Omega_sr = 4*pi*sin(delta_rad);
        schedule.name = 'Band half-width';
        schedule.unit = 'deg';
        schedule.description = 'values are band half-width around the equator of BandAxis';

    case 'ranges'
        vals = linspace(p.Results.RangeApertureDeg(1), ...
                        p.Results.RangeApertureDeg(2), ...
                        num_steps);

        half_angle_rad = deg2rad(vals / 2);

        schedule.values = vals;
        schedule.Omega_sr = 2*pi*(1 - cos(half_angle_rad));
        schedule.name = 'Angular range';
        schedule.unit = 'deg';
        schedule.description = 'legacy range mode; values are angular aperture in degrees';

    otherwise
        error('Unknown sampling_mode: %s', sampling_mode);
end

end