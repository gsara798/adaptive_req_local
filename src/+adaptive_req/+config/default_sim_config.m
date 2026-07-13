function cfg = default_sim_config(varargin)
%DEFAULT_SIM_CONFIG Default simulation configuration for adaptive REQ studies.
%
% Usage:
%   cfg = adaptive_req.config.default_sim_config();
%   cfg = adaptive_req.config.default_sim_config('f0', 300, 'cs_bg', 2.5);

p = inputParser;

addParameter(p, 'cs_bg', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MaskType', 'homogeneous', @(x) ischar(x) || isstring(x));

addParameter(p, 'Lx', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Lz', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'dx', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'dz', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'f0', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Nwaves', 500, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'AmpJitter', 0.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'SNR', inf, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'SourceSampling', 'cone', @(x) ischar(x) || isstring(x));
addParameter(p, 'AngularSamplingMethod', 'random', @(x) ischar(x) || isstring(x));
addParameter(p, 'ForceInPlaneWave', false, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'ConeAxis', [1 0 0], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'ConeHalfAngleDeg', 180, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'BandAxis', [0 1 0], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'BandHalfWidthDeg', 90, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'WaveModel', 'planewave', @(x) ischar(x) || isstring(x));
addParameter(p, 'UseParfor', true, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'Seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

parse(p, varargin{:});

cfg = p.Results;

cfg.MaskType = char(cfg.MaskType);
cfg.SourceSampling = char(cfg.SourceSampling);
cfg.AngularSamplingMethod = char(cfg.AngularSamplingMethod);
cfg.WaveModel = char(cfg.WaveModel);
cfg.UseParfor = logical(cfg.UseParfor);
cfg.ForceInPlaneWave = logical(cfg.ForceInPlaneWave);

cfg.Nx = round(cfg.Lx / cfg.dx) + 1;
cfg.Nz = round(cfg.Lz / cfg.dz) + 1;

cfg.x = ((1:cfg.Nx) - 1) * cfg.dx;
cfg.z = ((1:cfg.Nz) - 1) * cfg.dz;

cfg.k0_true = 2*pi*cfg.f0 / cfg.cs_bg;
cfg.lambda_true = cfg.cs_bg / cfg.f0;

end
