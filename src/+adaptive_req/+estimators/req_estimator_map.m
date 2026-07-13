function OUT = req_estimator_map(Uxz, cfg, feat_cfg, varargin)
%REQ_ESTIMATOR_MAP Estimate local REQ q and SWS maps from a 2D wavefield.
%
% This project-native estimator mirrors the role of the reference REQ map
% estimators, but uses the adaptive_req quantile implementation directly.
% It is intended for heterogeneous tests where we want spatial maps:
%
%   OUT = adaptive_req.estimators.req_estimator_map(Uxz, cfg, feat_cfg)
%
% Required inputs
%   Uxz      : nz-by-nx wavefield
%   cfg      : simulation/physical config with dx, dz, f0
%   feat_cfg : REQ/window config, usually default_feature_config output
%
% Optional name-value inputs
%   StepX/StepZ       : x/z stride in pixels. Default half window.
%   EdgeMode          : currently 'valid'. Full window must stay in field.
%   QuantileMode      : 'local_req', 'fixed', 'provided', or
%                       'theory_discrete'.
%   TheoryFieldType   : 'Diffuse3D', 'Diffuse2D', or 'SingleWave'.
%   FixedQuantile     : scalar q used when QuantileMode='fixed'.
%   ProvidedQuantiles : nZ-by-nX q map used when QuantileMode='provided'.
%   ReqOptions        : options passed to default_req_config.
%   ReturnFeatures    : store a per-window scalar-feature cell grid. Default true.
%   ReturnFeatureTable: return rows with scalar features. Default true.
%   StoreReqCurves    : store full per-window REQ curves. Default true.
%   ReuseReqSpectrumForFeatures : avoid a second FFT for features. Default true.
%   UseWindowParfor   : parallelize across REQ windows. Default false.
%   Verbose           : print progress. Default false.

p = inputParser;
p.FunctionName = 'adaptive_req.estimators.req_estimator_map';

addRequired(p, 'Uxz', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'cfg', @isstruct);
addRequired(p, 'feat_cfg', @isstruct);
addParameter(p, 'StridePixels', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
addParameter(p, 'StepX', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
addParameter(p, 'StepZ', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
addParameter(p, 'EdgeMode', 'valid', @(x) ischar(x) || isstring(x));
addParameter(p, 'QuantileMode', 'local_req', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'TheoryFieldType', 'Diffuse3D', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'TheoryMode', 'S2D', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'FixedQuantile', NaN, ...
    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ProvidedQuantiles', [], ...
    @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ReqOptions', {}, @iscell);
addParameter(p, 'ReturnFeatures', true, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ReturnFeatureTable', true, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'StoreReqCurves', true, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ReuseReqSpectrumForFeatures', true, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'UseWindowParfor', false, ...
    @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', false, ...
    @(x) islogical(x) || isnumeric(x));

parse(p, Uxz, cfg, feat_cfg, varargin{:});

stride_pixels = p.Results.StridePixels;
step_x = p.Results.StepX;
step_z = p.Results.StepZ;
edge_mode = lower(string(p.Results.EdgeMode));
quantile_mode = lower(string(p.Results.QuantileMode));
theory_field_type = string(p.Results.TheoryFieldType);
theory_mode = string(p.Results.TheoryMode);
fixed_quantile = p.Results.FixedQuantile;
provided_quantiles = p.Results.ProvidedQuantiles;
return_features = logical(p.Results.ReturnFeatures);
return_feature_table = logical(p.Results.ReturnFeatureTable);
compute_features = return_features || return_feature_table;
store_req_curves = logical(p.Results.StoreReqCurves);
reuse_req_spectrum_for_features = logical(p.Results.ReuseReqSpectrumForFeatures);
use_window_parfor = logical(p.Results.UseWindowParfor);
verbose = logical(p.Results.Verbose);

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, p.Results.ReqOptions{:});
req_cfg.precomputed = build_req_precomputed_geometry(cfg, req_cfg);

win_size = feat_cfg.win_size;
half_win = floor(win_size / 2);

if edge_mode ~= "valid"
    error('req_estimator_map currently supports EdgeMode=''valid'' only.');
end

if isempty(step_x) || isempty(step_z)
    if isempty(stride_pixels)
        stride_pixels = max(1, half_win);
    else
        stride_pixels = round(stride_pixels);
    end

    if isempty(step_x)
        step_x = stride_pixels;
    end

    if isempty(step_z)
        step_z = stride_pixels;
    end
end

step_x = round(step_x);
step_z = round(step_z);

[nz, nx] = size(Uxz);
x_centers = (1 + half_win):step_x:(nx - half_win);
z_centers = (1 + half_win):step_z:(nz - half_win);

if isempty(x_centers) || isempty(z_centers)
    error(['Field is too small for the requested REQ window. ', ...
           'Field size is %d-by-%d, win_size is %d.'], ...
           nz, nx, win_size);
end

nz_map = numel(z_centers);
nx_map = numel(x_centers);

q_map = nan(nz_map, nx_map);
k_map = nan(nz_map, nx_map);
cs_map = nan(nz_map, nx_map);
valid_map = false(nz_map, nx_map);
req_mappings = cell(nz_map, nx_map);
if store_req_curves
    req_curves = cell(nz_map, nx_map);
else
    req_curves = {};
end

if return_features
    features = cell(nz_map, nx_map);
else
    features = {};
end

if return_feature_table
    rows(nz_map * nx_map) = struct();
else
    rows = struct([]);
end

row_idx = 0;

if quantile_mode == "provided"
    if ~isequal(size(provided_quantiles), [nz_map, nx_map])
        error(['ProvidedQuantiles must have size nZ-by-nX = %d-by-%d. ', ...
               'Received %d-by-%d.'], ...
            nz_map, nx_map, size(provided_quantiles, 1), ...
            size(provided_quantiles, 2));
    end
elseif quantile_mode == "fixed"
    if ~isfinite(fixed_quantile)
        error('FixedQuantile must be finite when QuantileMode=''fixed''.');
    end
elseif quantile_mode == "theory_discrete"
    theory_out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
        cfg.dx, cfg.dz, cfg.f0, feat_cfg.cs_guess_used, ...
        'M', feat_cfg.M, ...
        'Gamma', feat_cfg.gamma_win, ...
        'PadFactor', feat_cfg.pad_factor, ...
        'Nbins', req_cfg.Nbins, ...
        'SmoothSigma', req_cfg.smooth_sigma, ...
        'TheoryMode', theory_mode, ...
        'FieldType', theory_field_type, ...
        'Plot', false);
    fixed_quantile = theory_out.q_th;
elseif quantile_mode ~= "local_req"
    error('Unknown QuantileMode: %s', quantile_mode);
end

if use_window_parfor
    nwin = nz_map * nx_map;
    q_vec = nan(nwin, 1);
    k_vec = nan(nwin, 1);
    cs_vec = nan(nwin, 1);
    valid_vec = false(nwin, 1);
    mapping_vec = cell(nwin, 1);
    curve_vec = cell(nwin, 1);
    feature_vec = cell(nwin, 1);
    row_vec = cell(nwin, 1);

    parfor wi = 1:nwin
        [iz, ix] = ind2sub([nz_map, nx_map], wi);
        provided_q_i = NaN;
        if quantile_mode == "provided"
            provided_q_i = provided_quantiles(iz, ix);
        end

        [q_vec(wi), k_vec(wi), cs_vec(wi), valid_vec(wi), ...
            mapping_vec{wi}, curve_vec{wi}, feature_vec{wi}, row_vec{wi}] = ...
            process_req_window( ...
            Uxz, cfg, req_cfg, feat_cfg, iz, ix, wi, ...
            x_centers(ix), z_centers(iz), half_win, quantile_mode, ...
            fixed_quantile, provided_q_i, compute_features, return_features, ...
            return_feature_table, store_req_curves, ...
            reuse_req_spectrum_for_features);
    end

    for wi = 1:nwin
        [iz, ix] = ind2sub([nz_map, nx_map], wi);
        q_map(iz, ix) = q_vec(wi);
        k_map(iz, ix) = k_vec(wi);
        cs_map(iz, ix) = cs_vec(wi);
        valid_map(iz, ix) = valid_vec(wi);
        req_mappings{iz, ix} = mapping_vec{wi};
        if store_req_curves
            req_curves{iz, ix} = curve_vec{wi};
        end
        if return_features
            features{iz, ix} = feature_vec{wi};
        end
    end

    if return_feature_table
        rows = [row_vec{:}];
        row_idx = numel(rows);
    end
else
    for iz = 1:nz_map
        for ix = 1:nx_map
            wi = sub2ind([nz_map, nx_map], iz, ix);
            provided_q_i = NaN;
            if quantile_mode == "provided"
                provided_q_i = provided_quantiles(iz, ix);
            end

            [q_i, k_i, cs_i, valid_i, mapping_i, curve_i, feature_i, row_i] = ...
                process_req_window( ...
                Uxz, cfg, req_cfg, feat_cfg, iz, ix, wi, ...
                x_centers(ix), z_centers(iz), half_win, quantile_mode, ...
                fixed_quantile, provided_q_i, compute_features, ...
                return_features, return_feature_table, store_req_curves, ...
                reuse_req_spectrum_for_features);

            q_map(iz, ix) = q_i;
            k_map(iz, ix) = k_i;
            cs_map(iz, ix) = cs_i;
            valid_map(iz, ix) = valid_i;
            req_mappings{iz, ix} = mapping_i;
            if store_req_curves
                req_curves{iz, ix} = curve_i;
            end
            if return_features
                features{iz, ix} = feature_i;
            end
            if return_feature_table
                row_idx = row_idx + 1;
                if row_idx == 1
                    rows = repmat(row_i, nz_map * nx_map, 1);
                else
                    rows(row_idx) = row_i;
                end
            end
        end

        if verbose
            fprintf('REQ map row %d / %d completed.\n', iz, nz_map);
        end
    end
end

x_m = (x_centers - 1) * cfg.dx;
z_m = (z_centers - 1) * cfg.dz;

OUT = struct();
OUT.q_map = q_map;
OUT.k_map = k_map;
OUT.cs_map = cs_map;
OUT.valid_map = valid_map;
OUT.x_idx = x_centers(:).';
OUT.z_idx = z_centers(:);
OUT.x_m = x_m(:).';
OUT.z_m = z_m(:);
OUT.req_mappings = req_mappings;
OUT.req_curves = req_curves;
OUT.features = features;
if return_feature_table
    OUT.feature_table = struct2table(rows(1:row_idx));
else
    OUT.feature_table = table();
end
OUT.win_size = win_size;
OUT.half_win = half_win;
OUT.stride_pixels = [];
OUT.step_x = step_x;
OUT.step_z = step_z;
OUT.edge_mode = edge_mode;
OUT.quantile_mode = quantile_mode;
if quantile_mode == "theory_discrete"
    OUT.theory = theory_out;
else
    OUT.theory = struct();
end
OUT.cfg = cfg;
OUT.feat_cfg = feat_cfg;
OUT.req_cfg = req_cfg;

end

function [q_i, k_i, cs_i, valid_i, mapping_i, curve_out, feature_out, row] = ...
    process_req_window( ...
    Uxz, cfg, req_cfg, feat_cfg, iz, ix, wi, cx, cz, half_win, ...
    quantile_mode, fixed_quantile, provided_q_i, compute_features, ...
    store_feature_grid, ...
    return_feature_table, store_req_curves, reuse_req_spectrum_for_features)

x_idx = (cx - half_win):(cx + half_win);
z_idx = (cz - half_win):(cz + half_win);
patch = Uxz(z_idx, x_idx);

[q_local_i, req_curve_i] = adaptive_req.quantile.compute_quantile_from_patch( ...
    patch, cfg, req_cfg);

switch quantile_mode
    case "local_req"
        q_i = q_local_i;
    case "fixed"
        q_i = fixed_quantile;
    case "theory_discrete"
        q_i = fixed_quantile;
    case "provided"
        q_i = provided_q_i;
    otherwise
        q_i = NaN;
end

k_i = adaptive_req.quantile.quantile_to_k(req_curve_i, q_i);
cs_i = adaptive_req.quantile.quantile_to_cs(req_curve_i, q_i, cfg.f0);
valid_i = isfinite(q_i) && isfinite(cs_i);
mapping_i = adaptive_req.quantile.make_req_mapping(req_curve_i);

if store_req_curves
    curve_out = req_curve_i;
else
    curve_out = [];
end

feature_out = [];
feat_scalar_i = struct();
shape_i = adaptive_req.quantile.extract_ecum_shape_features(req_curve_i);

if compute_features
    if reuse_req_spectrum_for_features
        feat_i = adaptive_req.features.req_extract_patch_features_from_curve( ...
            req_curve_i, feat_cfg);
    else
        feat_i = adaptive_req.features.req_extract_patch_features( ...
            patch, cfg.dx, cfg.dz, cfg.f0, feat_cfg.cs_guess_used, feat_cfg);
    end
    feat_scalar_i = feat_i.scalar;
    if store_feature_grid
        feature_out = feat_scalar_i;
    end
end

if return_feature_table
    row = struct();
    row.map_iz = iz;
    row.map_ix = ix;
    row.patch_idx = wi;
    row.cx = cx;
    row.cz = cz;
    row.x_center_m = (cx - 1) * cfg.dx;
    row.z_center_m = (cz - 1) * cfg.dz;
    row.q_local_req = q_local_i;
    row.q_pred = q_i;
    row.q_theory_discrete = NaN;
    if quantile_mode == "theory_discrete"
        row.q_theory_discrete = fixed_quantile;
    end
    row.k_pred = k_i;
    row.cs_pred = cs_i;
    row.REQ_M = feat_cfg.M;
    row.M = feat_cfg.M;
    row.SIM_f0 = cfg.f0;
    row.SIM_cs_bg = cfg.cs_bg;
    row.SIM_WaveModel = get_string_field(cfg, 'WaveModel', "");
    row.REQ_Nbins_effective = get_numeric_field( ...
        req_curve_i, 'Nbins_effective', NaN);
    row.req_mapping = {mapping_i};
    row = assign_scalar_fields_to_struct(row, feat_scalar_i, "");
    row = assign_scalar_fields_to_struct(row, shape_i, "");
else
    row = struct();
end

end

function rows = assign_scalar_fields_to_row(rows, row_idx, values, prefix)

if ~isstruct(values)
    return;
end

names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if isnumeric(value_i) && isscalar(value_i)
        rows(row_idx).(char(string(prefix) + string(names{i}))) = ...
            double(value_i);
    elseif islogical(value_i) && isscalar(value_i)
        rows(row_idx).(char(string(prefix) + string(names{i}))) = ...
            double(value_i);
    end
end

end

function row = assign_scalar_fields_to_struct(row, values, prefix)

if ~isstruct(values)
    return;
end

names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if isnumeric(value_i) && isscalar(value_i)
        row.(char(string(prefix) + string(names{i}))) = double(value_i);
    elseif islogical(value_i) && isscalar(value_i)
        row.(char(string(prefix) + string(names{i}))) = double(value_i);
    end
end

end

function pre = build_req_precomputed_geometry(cfg, req)

win_size = size(req.W2);
pad_z = req.PAD(1);
pad_x = req.PAD(2);
Nz2 = win_size(1) + 2 * pad_z;
Nx2 = win_size(2) + 2 * pad_x;

dkx = 2*pi / (Nx2 * cfg.dx);
dkz = 2*pi / (Nz2 * cfg.dz);

kx = (-floor(Nx2/2):ceil(Nx2/2)-1) * dkx;
kz = (-floor(Nz2/2):ceil(Nz2/2)-1) * dkz;

[KX, KZ] = meshgrid(kx, kz);
KR = sqrt(KX.^2 + KZ.^2);
TH = mod(atan2(KZ, KX), 2*pi);

kmax = min(max(abs(kx)), max(abs(kz)));
Nbins_requested = req.Nbins;

if ischar(Nbins_requested) || isstring(Nbins_requested)
    if lower(string(Nbins_requested)) == "auto"
        dk_DFT = min(dkx, dkz);
        oversample_factor = req.Nbins_auto_oversample;
        Nbins_effective = round(kmax / (dk_DFT / oversample_factor));
        Nbins_effective = max(req.Nbins_min, Nbins_effective);
    else
        error('Unknown Nbins value: %s', string(Nbins_requested));
    end
else
    Nbins_effective = round(double(Nbins_requested));
end

Nbins_effective = max(1, Nbins_effective);
k_edges = linspace(0, kmax, Nbins_effective + 1);
k_cent = 0.5 * (k_edges(1:end-1) + k_edges(2:end));
radial_bin = discretize(KR(:), k_edges);

pre = struct();
pre.Nz2 = Nz2;
pre.Nx2 = Nx2;
pre.pad_z = pad_z;
pre.pad_x = pad_x;
pre.dkx = dkx;
pre.dkz = dkz;
pre.kx = kx;
pre.kz = kz;
pre.KR = KR;
pre.TH = TH;
pre.kmax = kmax;
pre.Nbins_requested = Nbins_requested;
pre.Nbins_effective = Nbins_effective;
pre.k_edges = k_edges;
pre.k_cent = k_cent;
pre.radial_bin = radial_bin;

end

function val = get_numeric_field(S, name, default_val)

if isstruct(S) && isfield(S, name) && isnumeric(S.(name)) && ...
        isscalar(S.(name))
    val = double(S.(name));
else
    val = default_val;
end

end

function val = get_string_field(S, name, default_val)

if isstruct(S) && isfield(S, name)
    val = string(S.(name));
else
    val = string(default_val);
end

end
