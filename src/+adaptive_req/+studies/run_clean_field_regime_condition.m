function [T_condition, REP, META] = run_clean_field_regime_condition( ...
    CFG, regime, f0, cs_bg, dx, condition_id, varargin)
%RUN_CLEAN_FIELD_REGIME_CONDITION Generate one Test 18 physical condition.
%
% One simulated realization is reused across every requested REQ.M value.
% This avoids repeating wavefield synthesis when only the REQ window changes.

p = inputParser;
p.FunctionName = 'adaptive_req.studies.run_clean_field_regime_condition';
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
verbose = logical(p.Results.Verbose);

T_condition = table();
REP = struct([]);
META = struct();
META.condition_id = condition_id;
META.field_regime_label = string(regime.field_regime_label);
META.field_regime_variant = string(regime.field_regime_variant);
META.f0 = f0;
META.cs_bg = cs_bg;
META.dx = dx;
META.dz = dx;
META.success = false;
META.elapsed_s = NaN;

t_condition = tic;

for realization_idx = 1:CFG.EXP.num_realizations
    cfg = build_sim_cfg(CFG, regime, f0, cs_bg, dx, ...
        condition_id, realization_idx);

    if verbose
        fprintf('    realization %d/%d: simulate %s, N=%d\n', ...
            realization_idx, CFG.EXP.num_realizations, ...
            regime.field_regime_variant, regime.Nwaves);
    end
    sim = adaptive_req.simulate.run_single_simulation(cfg);
    verify_in_plane_requirement(sim, regime);

    for M = CFG.GRID.REQ_M(:)'
        [T_M, req_meta, central_patch] = extract_all_patches( ...
            sim, cfg, CFG, regime, M, condition_id, realization_idx);
        T_condition = concat_tables(T_condition, T_M);

        if is_representative_case(CFG, regime, f0, cs_bg, dx, M, realization_idx)
            REP = struct();
            REP.field_regime_label = string(regime.field_regime_label);
            REP.field_regime_variant = string(regime.field_regime_variant);
            REP.Nwaves = regime.Nwaves;
            REP.Is2D = regime.Is2D;
            REP.f0 = f0;
            REP.cs_bg = cs_bg;
            REP.dx = dx;
            REP.REQ_M = M;
            REP.x = sim.x;
            REP.z = sim.z;
            REP.Uxz = sim.Uxz;
            REP.central_patch = central_patch;
            REP.central_power_spectrum = ...
                log1p(abs(fftshift(fft2(central_patch))).^2);
            REP.req_meta = req_meta;
            REP.in_plane_fraction = in_plane_fraction(sim);
        end
    end
end

META.n_rows = height(T_condition);
META.elapsed_s = toc(t_condition);
META.success = true;

end

function cfg = build_sim_cfg(CFG, regime, f0, cs_bg, dx, condition_id, realization_idx)
cfg = adaptive_req.config.default_sim_config( ...
    'cs_bg', cs_bg, ...
    'MaskType', 'homogeneous', ...
    'Lx', CFG.SIM.Lx, ...
    'Lz', CFG.SIM.Lz, ...
    'dx', dx, ...
    'dz', dx, ...
    'f0', f0, ...
    'Nwaves', regime.Nwaves, ...
    'AmpJitter', CFG.SIM.AmpJitter, ...
    'SNR', CFG.SIM.SNR, ...
    'SourceSampling', 'ranges', ...
    'AngularSamplingMethod', char(regime.AngularSamplingMethod), ...
    'ForceInPlaneWave', regime.ForceInPlaneWave, ...
    'WaveModel', char(regime.WaveModel), ...
    'UseParfor', CFG.SIM.UseParfor, ...
    'Seed', CFG.EXP.seed_base + 10000*condition_id + realization_idx);
cfg.Is2D = logical(regime.Is2D);
cfg.AngleRange2D = [0 2*pi];
cfg.PhiRange = [0 2*pi];
cfg.ThetaRange = [0 pi];
cfg.DecayAlpha = CFG.SIM.DecayAlpha;
end

function [T, req_meta, central_patch] = extract_all_patches( ...
    sim, cfg, CFG, regime, M, condition_id, realization_idx)

feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', M, ...
    'cs_guess', CFG.GRID.REQ_cs_guess, ...
    'gamma_win', CFG.REQ.gamma_win, ...
    'pad_factor', CFG.REQ.pad_factor);
req_options = { ...
    'Nbins', CFG.REQ.Nbins, ...
    'Nbins_auto_oversample', CFG.REQ.Nbins_auto_oversample, ...
    'Nbins_min', CFG.REQ.Nbins_min, ...
    'smooth_sigma', CFG.REQ.smooth_sigma, ...
    'use_donut', CFG.REQ.use_donut, ...
    'donut_cs_min', CFG.REQ.donut_cs_min, ...
    'donut_cs_max', CFG.REQ.donut_cs_max, ...
    'donut_taper_rel', CFG.REQ.donut_taper_rel, ...
    'apply_donut_to_final_map', CFG.REQ.apply_donut_to_final_map};
[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg, req_options{:});

[q_global, global_curve, global_features] = ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz, cfg, req_cfg, feat_cfg);
global_shape = adaptive_req.quantile.extract_ecum_shape_features(global_curve);
global_mapping = adaptive_req.quantile.make_req_mapping(global_curve);

patch_pack = adaptive_req.simulate.build_patch_windows( ...
    cfg, feat_cfg, 'NumPatches', CFG.EXP.num_patches);
rows = cell(patch_pack.n_patches, 1);
central_patch = [];
central_distance = Inf;

for patch_idx = 1:patch_pack.n_patches
    xi = patch_pack.x_idx_list{patch_idx};
    zi = patch_pack.z_idx_list{patch_idx};
    patch = sim.Uxz(zi, xi);

    feat = adaptive_req.features.req_extract_patch_features( ...
        patch, cfg.dx, cfg.dz, cfg.f0, cfg.cs_bg, feat_cfg);
    [q_theory, req_curve] = adaptive_req.quantile.compute_quantile_from_patch( ...
        patch, cfg, req_cfg);
    local_shape = adaptive_req.quantile.extract_ecum_shape_features(req_curve);

    row = struct();
    row.condition_id = condition_id;
    row.field_regime_label = string(regime.field_regime_label);
    row.field_regime_variant = string(regime.field_regime_variant);
    row.realization_idx = realization_idx;
    row.patch_idx = patch_idx;
    row.patch_label = string(patch_pack.patch_labels(patch_idx));
    row.SIM_f0 = cfg.f0;
    row.SIM_cs_bg = cfg.cs_bg;
    row.SIM_Nwaves = cfg.Nwaves;
    row.SIM_Is2D = logical(regime.Is2D);
    row.SIM_ForceInPlaneWave = logical(regime.ForceInPlaneWave);
    row.SIM_WaveModel = string(regime.WaveModel);
    row.SIM_AngularSamplingMethod = string(regime.AngularSamplingMethod);
    row.SIM_SourceSampling = string(cfg.SourceSampling);
    row.SIM_SNR = cfg.SNR;
    row.SIM_dx = cfg.dx;
    row.SIM_dz = cfg.dz;
    row.Seed = cfg.Seed;
    row.REQ_M = feat_cfg.M;
    row.REQ_cs_guess = feat_cfg.cs_guess_used;
    row.REQ_win_size = feat_cfg.win_size;
    row.REQ_pad_factor = feat_cfg.pad_factor;
    row.REQ_gamma_win = feat_cfg.gamma_win;
    row.REQ_Nbins_effective = numeric_field(req_curve, 'Nbins_effective');
    row.cx = patch_pack.cx_list(patch_idx);
    row.cz = patch_pack.cz_list(patch_idx);
    row.x_center_m = sim.x(row.cx);
    row.z_center_m = sim.z(row.cz);
    row.patch_nx = numel(xi);
    row.patch_nz = numel(zi);
    row.lambda_true = cfg.cs_bg/cfg.f0;
    row.lambda_guess = feat_cfg.cs_guess_used/cfg.f0;
    row.window_length_x_m = feat_cfg.win_size*cfg.dx;
    row.window_length_z_m = feat_cfg.win_size*cfg.dz;
    row.M_eff_true_diag = row.window_length_x_m/row.lambda_true;
    row.M_eff_guess = row.window_length_x_m/row.lambda_guess;
    row.pixels_per_wavelength = row.lambda_true/cfg.dx;
    row.pixels_per_window = feat_cfg.win_size;
    row.k0_over_knyquist = (2*pi*cfg.f0/cfg.cs_bg)/(pi/cfg.dx);
    row.q_theory = q_theory;
    row.q_global_theory = q_global;
    row.q_local_minus_global = q_theory-q_global;
    row.k_req_q_theory = adaptive_req.quantile.quantile_to_k(req_curve, q_theory);
    row.cs_req_q_theory = adaptive_req.quantile.quantile_to_cs( ...
        req_curve, q_theory, cfg.f0);
    row.req_mapping = {adaptive_req.quantile.make_req_mapping(req_curve)};
    row.global_req_mapping = {global_mapping};
    row = assign_numeric_fields(row, feat.scalar, "");
    row = assign_numeric_fields(row, local_shape, "");
    row = assign_numeric_fields(row, global_features, "global_");
    row = assign_numeric_fields(row, global_shape, "global_");
    rows{patch_idx} = row;

    d = hypot(row.x_center_m-CFG.SIM.Lx/2, row.z_center_m-CFG.SIM.Lz/2);
    if d < central_distance
        central_distance = d;
        central_patch = patch;
    end
end

T = struct2table(vertcat(rows{:}));
req_meta = struct('win_size', feat_cfg.win_size, ...
    'Nbins_effective', numeric_field(global_curve, 'Nbins_effective'), ...
    'q_global', q_global);
end

function row = assign_numeric_fields(row, values, prefix)
if ~isstruct(values), return; end
names = fieldnames(values);
for i = 1:numel(names)
    value_i = values.(names{i});
    if (isnumeric(value_i) || islogical(value_i)) && isscalar(value_i)
        row.(char(string(prefix)+string(names{i}))) = double(value_i);
    end
end
end

function value = numeric_field(S, name)
if isfield(S, name) && isnumeric(S.(name)) && isscalar(S.(name))
    value = double(S.(name));
else
    value = NaN;
end
end

function tf = is_representative_case(CFG, regime, f0, cs_bg, dx, M, realization_idx)
canonical = ...
    (regime.field_regime_label == "directional_2D" && regime.Nwaves == 1) || ...
    (regime.field_regime_label == "diffuse_2D" && regime.Nwaves == 128) || ...
    (regime.field_regime_label == "partial_3D" && regime.Nwaves == 8) || ...
    (regime.field_regime_label == "diffuse_3D" && regime.Nwaves == 128);
tf = canonical && f0 == 500 && cs_bg == 3 && ...
    abs(dx-0.5e-3) < eps && M == 3 && realization_idx == 1 && ...
    ismember(M, CFG.GRID.REQ_M);
end

function verify_in_plane_requirement(sim, regime)
if regime.field_regime_label == "partial_3D"
    assert(min(abs(sim.diag.waveDirs.uy)) < 1e-7, ...
        'Partial 3D condition does not contain an in-plane source.');
end
end

function frac = in_plane_fraction(sim)
coverage = sim.diag.inPlaneCoverage;
if isfield(coverage, 'fraction_in_plane')
    frac = coverage.fraction_in_plane;
elseif isfield(coverage, 'count_in_plane') && isfield(coverage, 'n_waves')
    frac = coverage.count_in_plane/max(coverage.n_waves, 1);
else
    frac = NaN;
end
end

function T = concat_tables(A, B)
if isempty(A), T = B; return; end
if isempty(B), T = A; return; end
vars = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)], 'stable');
A = add_missing(A, vars);
B = add_missing(B, vars);
T = [A(:, cellstr(vars)); B(:, cellstr(vars))];
end

function T = add_missing(T, vars)
for i = 1:numel(vars)
    if ismember(vars(i), string(T.Properties.VariableNames)), continue; end
    T.(char(vars(i))) = nan(height(T), 1);
end
end
