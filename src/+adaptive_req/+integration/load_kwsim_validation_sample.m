function data = load_kwsim_validation_sample(sample_file)
%LOAD_KWSIM_VALIDATION_SAMPLE Load and validate a KWSIM REQ sample.
%
% The input file must contain:
%
%   req_validation_sample
%
% The public wavefield contract is always:
%
%   wavefield_complex_zx : [Nz,Nx]
%   axes.x_m
%   axes.z_m
%   truth.cs_m_s_zx
%   frequency_hz
%
% No transpose is applied in this adapter.

arguments
    sample_file {mustBeTextScalar}
end

sample_file = string(sample_file);

if ~isfile(sample_file)
    error( ...
        "adaptive_req:KwsimSampleNotFound", ...
        "KWSIM validation sample was not found: %s", ...
        sample_file);
end

loaded = load( ...
    sample_file, ...
    "req_validation_sample");

if ~isfield(loaded, "req_validation_sample")
    error( ...
        "adaptive_req:InvalidKwsimSample", ...
        ["The MAT file must contain a variable named " ...
         "'req_validation_sample'."]);
end

sample = loaded.req_validation_sample;

required_fields = [
    "sample_schema_version"
    "source_dimension"
    "wavefield_complex_zx"
    "axes"
    "truth"
    "frequency_hz"
    "orientation"
];

for field_name = required_fields.'
    if ~isfield(sample, field_name)
        error( ...
            "adaptive_req:InvalidKwsimSample", ...
            "KWSIM sample is missing field '%s'.", ...
            field_name);
    end
end

if string(sample.sample_schema_version) ~= "1.0"
    error( ...
        "adaptive_req:UnsupportedKwsimSampleSchema", ...
        "Unsupported KWSIM sample schema: %s", ...
        string(sample.sample_schema_version));
end

if string(sample.orientation) ~= "[Nz,Nx]"
    error( ...
        "adaptive_req:InvalidKwsimOrientation", ...
        "Expected KWSIM orientation [Nz,Nx], received %s.", ...
        string(sample.orientation));
end

wavefield_zx = sample.wavefield_complex_zx;

if ~isnumeric(wavefield_zx) || ...
        ~ismatrix(wavefield_zx) || ...
        isempty(wavefield_zx)
    error( ...
        "adaptive_req:InvalidKwsimWavefield", ...
        "wavefield_complex_zx must be a nonempty numeric matrix.");
end

if ~all(isfinite(wavefield_zx), "all")
    error( ...
        "adaptive_req:InvalidKwsimWavefield", ...
        "wavefield_complex_zx contains non-finite values.");
end

if ~isfield(sample.axes, "x_m") || ...
        ~isfield(sample.axes, "z_m")
    error( ...
        "adaptive_req:InvalidKwsimAxes", ...
        "KWSIM sample must contain x_m and z_m axes.");
end

x_m = double(sample.axes.x_m(:));
z_m = double(sample.axes.z_m(:));

expected_size = [
    numel(z_m), ...
    numel(x_m)
];

if ~isequal(size(wavefield_zx), expected_size)
    error( ...
        "adaptive_req:InvalidKwsimWavefieldSize", ...
        ["Wavefield size %s is inconsistent with " ...
         "the x/z axes %s."], ...
        mat2str(size(wavefield_zx)), ...
        mat2str(expected_size));
end

if ~isfield(sample.truth, "cs_m_s_zx")
    error( ...
        "adaptive_req:InvalidKwsimTruth", ...
        "KWSIM sample is missing truth.cs_m_s_zx.");
end

cs_truth_zx = double(sample.truth.cs_m_s_zx);

if ~isequal(size(cs_truth_zx), expected_size)
    error( ...
        "adaptive_req:InvalidKwsimTruth", ...
        "Truth SWS map size is inconsistent with the wavefield.");
end

if ~all(isfinite(cs_truth_zx), "all") || ...
        any(cs_truth_zx <= 0, "all")
    error( ...
        "adaptive_req:InvalidKwsimTruth", ...
        "Truth SWS map must contain finite positive values.");
end

frequency_hz = double(sample.frequency_hz);

if ~isscalar(frequency_hz) || ...
        ~isfinite(frequency_hz) || ...
        frequency_hz <= 0
    error( ...
        "adaptive_req:InvalidKwsimFrequency", ...
        "KWSIM frequency must be a finite positive scalar.");
end

dx_m = resolve_spacing( ...
    sample, ...
    "dx_m", ...
    x_m);

dz_m = resolve_spacing( ...
    sample, ...
    "dz_m", ...
    z_m);

validate_uniform_axis( ...
    x_m, ...
    dx_m, ...
    "x");

validate_uniform_axis( ...
    z_m, ...
    dz_m, ...
    "z");

data = struct();

data.sample_file = sample_file;
data.sample_schema_version = ...
    string(sample.sample_schema_version);

data.source_dimension = ...
    double(sample.source_dimension);

data.wavefield_zx = ...
    wavefield_zx;

data.axes = struct();
data.axes.x_m = x_m;
data.axes.z_m = z_m;

data.dx_m = dx_m;
data.dz_m = dz_m;
data.frequency_hz = frequency_hz;

data.truth = struct();
data.truth.cs_m_s_zx = cs_truth_zx;

if isfield(sample.truth, "material_id_zx")
    data.truth.material_id_zx = ...
        sample.truth.material_id_zx;
else
    data.truth.material_id_zx = [];
end

if isfield(sample, "quantity")
    data.quantity = string(sample.quantity);
else
    data.quantity = "";
end

if isfield(sample, "units")
    data.units = string(sample.units);
else
    data.units = "";
end

if isfield(sample, "phasor_convention")
    data.phasor_convention = ...
        string(sample.phasor_convention);
else
    data.phasor_convention = "";
end

if isfield(sample, "extraction")
    data.extraction = sample.extraction;
else
    data.extraction = struct();
end

if isfield(sample, "simulation_valid")
    data.simulation_valid = ...
        sample.simulation_valid;
else
    data.simulation_valid = [];
end

if isfield(sample, "req_readiness")
    data.req_readiness = ...
        sample.req_readiness;
else
    data.req_readiness = struct();
end

data.raw_sample = sample;

end


function spacing_m = ...
    resolve_spacing(sample, field_name, axis_m)

spacing_m = NaN;

if isfield(sample, "spacing") && ...
        isfield(sample.spacing, field_name)
    spacing_m = double( ...
        sample.spacing.(field_name));
end

if ~isscalar(spacing_m) || ...
        ~isfinite(spacing_m) || ...
        spacing_m <= 0

    if numel(axis_m) < 2
        error( ...
            "adaptive_req:InvalidKwsimAxes", ...
            "At least two spatial samples are required.");
    end

    spacing_m = median(diff(axis_m));
end

if ~isfinite(spacing_m) || ...
        spacing_m <= 0
    error( ...
        "adaptive_req:InvalidKwsimSpacing", ...
        "Spatial sampling must be positive.");
end

end


function validate_uniform_axis( ...
    axis_m, spacing_m, axis_name)

if numel(axis_m) < 2
    error( ...
        "adaptive_req:InvalidKwsimAxes", ...
        "Axis %s must contain at least two samples.", ...
        axis_name);
end

differences = diff(axis_m);

tolerance = max( ...
    1e-12, ...
    1e-6 * abs(spacing_m));

if max(abs(differences - spacing_m)) > tolerance
    error( ...
        "adaptive_req:NonuniformKwsimAxis", ...
        "Axis %s is not uniformly sampled.", ...
        axis_name);
end

end
