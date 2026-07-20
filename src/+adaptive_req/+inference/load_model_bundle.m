function [bundle, metadata] = load_model_bundle(bundle_file, options)
%LOAD_MODEL_BUNDLE Load and cache one frozen adaptive-REQ model bundle.
%
% The MAT file must contain MODEL_BUNDLE. The bundle is cached for the
% current MATLAB session because trained ensemble files may be very large.

arguments
    bundle_file {mustBeTextScalar}
    options.ForceReload (1,1) logical = false
end

persistent cached_path
persistent cached_bytes
persistent cached_datenum
persistent cached_bundle

bundle_file = absolute_path(string(bundle_file));

if ~isfile(bundle_file)
    error( ...
        "adaptive_req:ModelBundleNotFound", ...
        "Model bundle was not found: %s", ...
        bundle_file);
end

file_info = dir(bundle_file);

cache_matches = ...
    ~options.ForceReload && ...
    ~isempty(cached_bundle) && ...
    string(cached_path) == bundle_file && ...
    cached_bytes == file_info.bytes && ...
    cached_datenum == file_info.datenum;

if ~cache_matches
    fprintf("Loading frozen model bundle:\n%s\n", ...
        bundle_file);

    timer = tic;

    loaded = load( ...
        bundle_file, ...
        "MODEL_BUNDLE");

    if ~isfield(loaded, "MODEL_BUNDLE")
        error( ...
            "adaptive_req:InvalidModelBundle", ...
            "MAT file does not contain MODEL_BUNDLE.");
    end

    validate_bundle(loaded.MODEL_BUNDLE);

    cached_path = bundle_file;
    cached_bytes = file_info.bytes;
    cached_datenum = file_info.datenum;
    cached_bundle = loaded.MODEL_BUNDLE;

    fprintf("Model bundle loaded in %.2f s.\n", ...
        toc(timer));
end

bundle = cached_bundle;

metadata = struct();
metadata.path = bundle_file;
metadata.bytes_on_disk = file_info.bytes;
metadata.loaded_from_cache = cache_matches;
metadata.model_names = ...
    string(bundle.MODELS.model_names);

end


function validate_bundle(bundle)

required_bundle_fields = [
    "MODELS"
    "BASE_FEATURES"
];

for field_name = required_bundle_fields.'
    if ~isfield(bundle, field_name)
        error( ...
            "adaptive_req:InvalidModelBundle", ...
            "MODEL_BUNDLE is missing '%s'.", ...
            field_name);
    end
end

required_model_fields = [
    "composition"
    "q"
    "model_names"
];

for field_name = required_model_fields.'
    if ~isfield(bundle.MODELS, field_name)
        error( ...
            "adaptive_req:InvalidModelBundle", ...
            "MODEL_BUNDLE.MODELS is missing '%s'.", ...
            field_name);
    end
end

required_names = [
    "q_spectrum_only"
    "q_spectrum_plus_composition"
];

available_names = ...
    string(bundle.MODELS.model_names);

if ~all(ismember(required_names, available_names))
    error( ...
        "adaptive_req:InvalidModelBundle", ...
        "Bundle does not contain the required frozen q models.");
end

end


function path_value = absolute_path(path_value)

[status, attributes] = fileattrib(path_value);

if status
    path_value = string(attributes.Name);
end

end
