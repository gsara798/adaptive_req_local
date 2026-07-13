function mapping = make_req_mapping(req_curve, varargin)
%MAKE_REQ_MAPPING Keep only the 1D REQ mapping needed for q-to-k conversion.
%
% The full req_curve contains large 2D diagnostic arrays. Converting a
% predicted quantile to wavenumber and SWS only requires Ecum and k_cent.

p = inputParser;
p.FunctionName = 'adaptive_req.quantile.make_req_mapping';

addRequired(p, 'req_curve', @isstruct);
addParameter(p, 'StorageType', 'single', @(x) ischar(x) || isstring(x));

parse(p, req_curve, varargin{:});

storage_type = lower(string(p.Results.StorageType));

if ~ismember(storage_type, ["single", "double"])
    error('StorageType must be ''single'' or ''double''.');
end

if ~isfield(req_curve, 'Ecum') || ~isfield(req_curve, 'k_cent')
    error('req_curve must contain Ecum and k_cent.');
end

Ecum = req_curve.Ecum(:);
k_cent = req_curve.k_cent(:);

if storage_type == "single"
    Ecum = single(Ecum);
    k_cent = single(k_cent);
else
    Ecum = double(Ecum);
    k_cent = double(k_cent);
end

mapping = struct();
mapping.Ecum = Ecum;
mapping.k_cent = k_cent;
mapping.f0 = getfield_with_default(req_curve, 'f0', NaN);
mapping.Nbins_effective = getfield_with_default(req_curve, ...
    'Nbins_effective', numel(k_cent));
mapping.storage_type = char(storage_type);

end

function val = getfield_with_default(S, field_name, default_val)

if isfield(S, field_name)
    val = S.(field_name);
else
    val = default_val;
end

end
