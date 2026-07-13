function [T_raw, sweep] = run_aperture_sweep(cfg, feat_cfg, varargin)
%RUN_APERTURE_SWEEP Run a full aperture sweep for adaptive REQ studies.
%
% This function loops over aperture steps, realizations, and local patches.
% For each simulated wavefield, it extracts local patches and computes:
%
%   1. Local spectral features.
%   2. Reference REQ quantile q_theory.
%   3. k and cs recovered from q_theory.
%   4. REQ metadata needed for later q-to-cs analysis.
%
% Usage
%   [T_raw, sweep] = adaptive_req.studies.run_aperture_sweep(cfg, feat_cfg);
%
%   [T_raw, sweep] = adaptive_req.studies.run_aperture_sweep(RUN);
%
% Output
%   T_raw:
%       One row per aperture step, realization, and patch.
%
%   sweep:
%       Metadata, schedule, patch information, and optional wavefields.

%% Allow config-driven calling mode

if is_run_struct(cfg)

    RUN = cfg;

    if nargin >= 2
        extra_args = [{feat_cfg}, varargin];
    else
        extra_args = {};
    end

    [cfg, feat_cfg, varargin] = unpack_run_for_aperture_sweep(RUN, extra_args);

end

%% Parse inputs

p = inputParser;
p.FunctionName = 'adaptive_req.studies.run_aperture_sweep';

addRequired(p, 'cfg', @isstruct);
addRequired(p, 'feat_cfg', @isstruct);

addParameter(p, 'SamplingMode', 'cone', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'NumSteps', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'StepIndices', [], ...
    @(x) isempty(x) || isnumeric(x));

addParameter(p, 'NumRealizations', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'NumPatches', 9, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'SeedBase', 1000, ...
    @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'PatchOptions', {}, @iscell);
addParameter(p, 'ReqOptions', {}, @iscell);

addParameter(p, 'StoreWavefields', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StoreReqCurve', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StoreReqMapping', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'ComputeGlobalReq', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StoreGlobalReqMapping', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StoreReqMetadata', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StoreFeatureStruct', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'Verbose', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'PlotStepDiagnostics', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveStepDiagnostics', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'StepDiagnosticPatchIndex', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, 'StepDiagnosticDir', '', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'StepDiagnosticVisible', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'CloseStepDiagnosticsAfterSave', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveDiagnosticPNG', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveDiagnosticPDF', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'SaveDiagnosticFIG', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'DiagnosticResolution', 300, ...
    @(x) isnumeric(x) && isscalar(x));

parse(p, cfg, feat_cfg, varargin{:});

opt = p.Results;

sampling_mode = lower(char(opt.SamplingMode));
num_steps = round(opt.NumSteps);
num_realizations = round(opt.NumRealizations);
num_patches = round(opt.NumPatches);

store_wavefields = logical(opt.StoreWavefields);
store_req_curve = logical(opt.StoreReqCurve);
store_req_mapping = logical(opt.StoreReqMapping);
compute_global_req = logical(opt.ComputeGlobalReq);
store_global_req_mapping = logical(opt.StoreGlobalReqMapping);
store_req_metadata = logical(opt.StoreReqMetadata);
store_feature_struct = logical(opt.StoreFeatureStruct);
verbose = logical(opt.Verbose);

plot_step_diagnostics = logical(opt.PlotStepDiagnostics);
save_step_diagnostics = logical(opt.SaveStepDiagnostics);

step_diagnostic_patch_idx = round(opt.StepDiagnosticPatchIndex);
step_diagnostic_visible = logical(opt.StepDiagnosticVisible);
close_step_diagnostics_after_save = logical(opt.CloseStepDiagnosticsAfterSave);

step_diagnostic_dir = char(opt.StepDiagnosticDir);

if save_step_diagnostics && isempty(step_diagnostic_dir)
    step_diagnostic_dir = fullfile(pwd, 'figures', 'step_diagnostics');
end

if save_step_diagnostics && ~exist(step_diagnostic_dir, 'dir')
    mkdir(step_diagnostic_dir);
end

%% Build aperture schedule

schedule = adaptive_req.simulate.build_aperture_schedule( ...
    sampling_mode, ...
    num_steps);

if isempty(opt.StepIndices)
    step_indices = 1:num_steps;
else
    step_indices = opt.StepIndices(:).';
end

if any(step_indices < 1) || any(step_indices > num_steps)
    error('StepIndices contains values outside the valid range 1:NumSteps.');
end

%% Feature names

if isfield(feat_cfg, 'feature_names')
    feature_names = feat_cfg.feature_names;
else
    error('feat_cfg must contain feature_names. Check default_feature_config.m.');
end

%% Initialize sweep metadata

sweep = struct();

sweep.cfg_base = cfg;
sweep.feat_cfg_base = feat_cfg;
sweep.req_options = opt.ReqOptions;
sweep.patch_options = opt.PatchOptions;

sweep.schedule = schedule;
sweep.sampling_mode = sampling_mode;
sweep.num_steps = num_steps;
sweep.step_indices = step_indices;
sweep.num_realizations = num_realizations;
sweep.num_patches = num_patches;

sweep.feature_names = feature_names;

sweep.store_wavefields = store_wavefields;
sweep.store_req_curve = store_req_curve;
sweep.store_req_mapping = store_req_mapping;
sweep.compute_global_req = compute_global_req;
sweep.store_global_req_mapping = store_global_req_mapping;
sweep.store_req_metadata = store_req_metadata;
sweep.store_feature_struct = store_feature_struct;

if store_wavefields
    sweep.sims = cell(num_steps, num_realizations);
else
    sweep.sims = {};
end

sweep.plot_step_diagnostics = plot_step_diagnostics;
sweep.save_step_diagnostics = save_step_diagnostics;
sweep.step_diagnostic_patch_idx = step_diagnostic_patch_idx;
sweep.step_diagnostic_dir = step_diagnostic_dir;

%% Main sweep

n_rows_est = numel(step_indices) * num_realizations * num_patches;

rows(n_rows_est) = struct();
row_idx = 0;

if verbose
    fprintf('\nRunning aperture sweep.\n');
    fprintf('Sampling mode      : %s\n', sampling_mode);
    fprintf('Number of steps    : %d\n', numel(step_indices));
    fprintf('Realizations/step  : %d\n', num_realizations);
    fprintf('Patches/simulation : %d\n', num_patches);
end

for sidx = 1:numel(step_indices)

    step_idx = step_indices(sidx);
    step_value = schedule.values(step_idx);

    cfg_step = adaptive_req.simulate.apply_sampling_step( ...
        cfg, ...
        sampling_mode, ...
        step_value);

    if verbose
        fprintf('\nStep %d / %d\n', step_idx, num_steps);
        fprintf('  Aperture value = %.3f %s\n', step_value, schedule.unit);
        fprintf('  Omega = %.4f sr\n', schedule.Omega_sr(step_idx));
    end

    for r = 1:num_realizations

        cfg_i = cfg_step;
        cfg_i.Seed = opt.SeedBase + 10000 * step_idx + r;

        [req_cfg_i, feat_cfg_i] = adaptive_req.config.default_req_config( ...
            cfg_i, ...
            feat_cfg, ...
            opt.ReqOptions{:});

        sim = adaptive_req.simulate.run_single_simulation(cfg_i);

        global_req = struct( ...
            'q_theory', NaN, ...
            'mapping', [], ...
            'curve', [], ...
            'features', struct(), ...
            'shape_features', struct());

        if compute_global_req
            [global_req.q_theory, global_req.curve, global_req.features] = ...
                adaptive_req.quantile.compute_global_quantile_from_field( ...
                    sim.Uxz, cfg_i, req_cfg_i, feat_cfg_i);

            global_req.shape_features = ...
                adaptive_req.quantile.extract_ecum_shape_features( ...
                    global_req.curve);

            global_req.mapping = adaptive_req.quantile.make_req_mapping( ...
                global_req.curve);
        end

        patch_pack = adaptive_req.simulate.build_patch_windows( ...
            cfg_i, ...
            feat_cfg_i, ...
            'NumPatches', num_patches, ...
            opt.PatchOptions{:});

        if store_wavefields
            sweep.sims{step_idx, r} = sim;
        end

        q_local_this_realization = nan(1, patch_pack.n_patches);

        rep_diag = struct();
        rep_diag_found = false;

        for pidx = 1:patch_pack.n_patches

            xi = patch_pack.x_idx_list{pidx};
            zi = patch_pack.z_idx_list{pidx};

            V_patch = sim.Uxz(zi, xi);

            feat = adaptive_req.features.req_extract_patch_features( ...
                V_patch, ...
                cfg_i.dx, ...
                cfg_i.dz, ...
                cfg_i.f0, ...
                cfg_i.cs_bg, ...
                feat_cfg_i);

            [q_theory, req_curve] = adaptive_req.quantile.compute_quantile_from_patch( ...
                V_patch, ...
                cfg_i, ...
                req_cfg_i);

            req_shape_features = ...
                adaptive_req.quantile.extract_ecum_shape_features(req_curve);

            q_local_this_realization(pidx) = q_theory;

            if pidx == step_diagnostic_patch_idx

                rep_diag = req_curve;
                rep_diag_found = true;

                if ~isfield(rep_diag, 'q')
                    rep_diag.q = q_theory;
                end

                if ~isfield(rep_diag, 'k_real')
                    rep_diag.k_real = 2*pi*cfg_i.f0/cfg_i.cs_bg;
                end
            end

            k_req_q_theory = adaptive_req.quantile.quantile_to_k( ...
                req_curve, ...
                q_theory);

            cs_req_q_theory = adaptive_req.quantile.quantile_to_cs( ...
                req_curve, ...
                q_theory, ...
                cfg_i.f0);

            row_idx = row_idx + 1;

            %% Row identifiers

            rows(row_idx).step_idx = step_idx;
            rows(row_idx).realization = r;
            rows(row_idx).realization_idx = r;
            rows(row_idx).patch_idx = pidx;
            rows(row_idx).patch_label = string(patch_pack.patch_labels(pidx));

            %% Aperture metadata

            rows(row_idx).aperture_value = step_value;
            rows(row_idx).Omega_sr = schedule.Omega_sr(step_idx);
            rows(row_idx).omega_sr = schedule.Omega_sr(step_idx);
            rows(row_idx).omega_mean = schedule.Omega_sr(step_idx);

            rows(row_idx).aperture_name = string(schedule.name);
            rows(row_idx).aperture_unit = string(schedule.unit);

            rows(row_idx).ConeHalfAngleDeg = get_optional_numeric( ...
                cfg_i, ...
                'ConeHalfAngleDeg', ...
                NaN);

            rows(row_idx).BandHalfWidthDeg = get_optional_numeric( ...
                cfg_i, ...
                'BandHalfWidthDeg', ...
                NaN);

            %% Simulation metadata

            rows(row_idx).SIM_WaveModel = get_optional_string( ...
                cfg_i, ...
                'WaveModel', ...
                "");

            rows(row_idx).SIM_SourceSampling = get_optional_string( ...
                cfg_i, ...
                'SourceSampling', ...
                "");

            rows(row_idx).SIM_AngularSamplingMethod = get_optional_string( ...
                cfg_i, ...
                'AngularSamplingMethod', ...
                "");

            rows(row_idx).SIM_f0 = cfg_i.f0;
            rows(row_idx).SIM_cs_bg = cfg_i.cs_bg;
            rows(row_idx).SIM_Nwaves = cfg_i.Nwaves;
            rows(row_idx).SIM_SNR = cfg_i.SNR;

            rows(row_idx).SIM_AmpJitter = get_optional_numeric( ...
                cfg_i, ...
                'AmpJitter', ...
                NaN);

            rows(row_idx).SIM_dx = cfg_i.dx;
            rows(row_idx).SIM_dz = cfg_i.dz;

            rows(row_idx).SIM_Lx = get_optional_numeric( ...
                cfg_i, ...
                'Lx', ...
                NaN);

            rows(row_idx).SIM_Lz = get_optional_numeric( ...
                cfg_i, ...
                'Lz', ...
                NaN);

            rows(row_idx).Seed = cfg_i.Seed;

            % Backward-compatible aliases.
            rows(row_idx).f0 = cfg_i.f0;
            rows(row_idx).cs_true = cfg_i.cs_bg;
            rows(row_idx).k0_true = 2*pi*cfg_i.f0/cfg_i.cs_bg;

            %% Patch metadata

            rows(row_idx).cx = patch_pack.cx_list(pidx);
            rows(row_idx).cz = patch_pack.cz_list(pidx);

            rows(row_idx).x_idx_start = xi(1);
            rows(row_idx).x_idx_end = xi(end);
            rows(row_idx).z_idx_start = zi(1);
            rows(row_idx).z_idx_end = zi(end);

            rows(row_idx).patch_nx = numel(xi);
            rows(row_idx).patch_nz = numel(zi);

            %% REQ window metadata

            lambda_true_i = cfg_i.cs_bg / cfg_i.f0;
            lambda_guess_i = feat_cfg_i.cs_guess_used / cfg_i.f0;

            window_length_x_i = feat_cfg_i.win_size * cfg_i.dx;
            window_length_z_i = feat_cfg_i.win_size * cfg_i.dz;

            M_eff_true_i = window_length_x_i / lambda_true_i;
            M_eff_guess_i = window_length_x_i / lambda_guess_i;

            rows(row_idx).REQ_M = feat_cfg_i.M;
            rows(row_idx).REQ_cs_guess = feat_cfg_i.cs_guess_used;
            rows(row_idx).REQ_gamma_win = feat_cfg_i.gamma_win;
            rows(row_idx).REQ_pad_factor = feat_cfg_i.pad_factor;

            rows(row_idx).REQ_win_size = feat_cfg_i.win_size;
            rows(row_idx).REQ_half_win = feat_cfg_i.half_win;

            rows(row_idx).lambda_true = lambda_true_i;
            rows(row_idx).lambda_guess = lambda_guess_i;
            rows(row_idx).lambda_guess_used = feat_cfg_i.lambda_guess_used;

            rows(row_idx).window_length_x_m = window_length_x_i;
            rows(row_idx).window_length_z_m = window_length_z_i;

            % Diagnostic effective window size.
            % This uses cs_bg, so it is diagnostic, not an experimental input.
            rows(row_idx).M_eff_true_diag = M_eff_true_i;
            rows(row_idx).M_eff_guess = M_eff_guess_i;

            rows(row_idx).cs_guess_to_true_ratio = ...
                feat_cfg_i.cs_guess_used / cfg_i.cs_bg;

            rows(row_idx).cs_true_to_guess_ratio = ...
                cfg_i.cs_bg / feat_cfg_i.cs_guess_used;

            % Backward-compatible aliases.
            rows(row_idx).M = feat_cfg_i.M;
            rows(row_idx).cs_guess_used = feat_cfg_i.cs_guess_used;
            rows(row_idx).M_eff = M_eff_true_i;
            rows(row_idx).win_size = feat_cfg_i.win_size;

            %% q, k, and cs target quantities

            rows(row_idx).q_theory = q_theory;
            rows(row_idx).q_reference = q_theory;

            if compute_global_req
                rows(row_idx).q_global_theory = global_req.q_theory;
                rows(row_idx).q_local_minus_global = ...
                    q_theory - global_req.q_theory;
            end

            rows(row_idx).k_req_q_theory = k_req_q_theory;
            rows(row_idx).cs_req_q_theory = cs_req_q_theory;

            rows(row_idx).cs_req_q_theory_error = ...
                cs_req_q_theory - cfg_i.cs_bg;

            rows(row_idx).cs_req_q_theory_abs_error = ...
                abs(cs_req_q_theory - cfg_i.cs_bg);

            rows(row_idx).cs_req_q_theory_rel_error = ...
                (cs_req_q_theory - cfg_i.cs_bg) / cfg_i.cs_bg;

            rows(row_idx).cs_req_q_theory_abs_rel_error = ...
                abs(rows(row_idx).cs_req_q_theory_rel_error);

            rows(row_idx).cs_req_q_theory_rel_error_pct = ...
                100 * rows(row_idx).cs_req_q_theory_rel_error;

            rows(row_idx).cs_req_q_theory_abs_rel_error_pct = ...
                100 * rows(row_idx).cs_req_q_theory_abs_rel_error;

            %% REQ curve metadata

            if store_req_metadata

                rows(row_idx).REQ_Nbins_requested = ...
                    get_curve_field_string(req_curve, 'Nbins_requested');

                rows(row_idx).REQ_Nbins_effective = ...
                    get_curve_field_numeric(req_curve, 'Nbins_effective');

                rows(row_idx).REQ_Nbins_auto_oversample = ...
                    get_optional_numeric(req_cfg_i, 'Nbins_auto_oversample', NaN);

                rows(row_idx).REQ_Nbins_min = ...
                    get_optional_numeric(req_cfg_i, 'Nbins_min', NaN);

                rows(row_idx).REQ_smooth_sigma = ...
                    get_optional_numeric(req_cfg_i, 'smooth_sigma', NaN);

                rows(row_idx).REQ_use_donut = ...
                    get_optional_logical(req_cfg_i, 'use_donut', false);

                rows(row_idx).REQ_donut_cs_min = ...
                    get_optional_numeric(req_cfg_i, 'donut_cs_min', NaN);

                rows(row_idx).REQ_donut_cs_max = ...
                    get_optional_numeric(req_cfg_i, 'donut_cs_max', NaN);

                rows(row_idx).REQ_donut_taper_rel = ...
                    get_optional_numeric(req_cfg_i, 'donut_taper_rel', NaN);

                rows(row_idx).req_dkx = ...
                    get_curve_field_numeric(req_curve, 'dkx');

                rows(row_idx).req_dkz = ...
                    get_curve_field_numeric(req_curve, 'dkz');

                rows(row_idx).req_kmax = ...
                    get_curve_field_numeric(req_curve, 'kmax');

                rows(row_idx).req_k_real = ...
                    get_curve_field_numeric(req_curve, 'k_real');

                rows(row_idx).req_k0_true = ...
                    get_curve_field_numeric(req_curve, 'k0_true');

                rows(row_idx).req_win_size = ...
                    get_curve_field_numeric(req_curve, 'win_size');

                rows(row_idx).req_half_win = ...
                    get_curve_field_numeric(req_curve, 'half_win');

                rows(row_idx).req_pad_factor = ...
                    get_curve_field_numeric(req_curve, 'pad_factor');

                rows(row_idx).req_smooth_sigma = ...
                    get_curve_field_numeric(req_curve, 'smooth_sigma');

                if compute_global_req
                    rows(row_idx).global_REQ_Nbins_effective = ...
                        get_curve_field_numeric( ...
                            global_req.curve, ...
                            'Nbins_effective');

                    rows(row_idx).global_req_dkx = ...
                        get_curve_field_numeric(global_req.curve, 'dkx');

                    rows(row_idx).global_req_dkz = ...
                        get_curve_field_numeric(global_req.curve, 'dkz');

                    rows(row_idx).global_req_kmax = ...
                        get_curve_field_numeric(global_req.curve, 'kmax');
                end

            end

            %% Local scalar features

            for fn_idx = 1:numel(feature_names)

                fn = get_feature_name(feature_names, fn_idx);

                if isfield(feat.scalar, fn)
                    rows(row_idx).(fn) = feat.scalar.(fn);
                elseif isfield(req_shape_features, fn)
                    rows(row_idx).(fn) = req_shape_features.(fn);
                else
                    rows(row_idx).(fn) = NaN;
                end
            end

            if compute_global_req
                rows = assign_prefixed_numeric_fields_to_row( ...
                    rows, ...
                    row_idx, ...
                    global_req.features, ...
                    "global_");

                rows = assign_prefixed_numeric_fields_to_row( ...
                    rows, ...
                    row_idx, ...
                    global_req.shape_features, ...
                    "global_");
            end

            %% Optional heavy outputs

            if store_req_curve
                rows(row_idx).req_curve = {req_curve};
            end

            if store_req_mapping
                rows(row_idx).req_mapping = { ...
                    adaptive_req.quantile.make_req_mapping(req_curve)};
            end

            if compute_global_req && store_global_req_mapping
                rows(row_idx).global_req_mapping = {global_req.mapping};
            end

            if store_feature_struct
                rows(row_idx).feat = {feat};
            end

        end

        %% Step diagnostics, once per realization

        if (plot_step_diagnostics || save_step_diagnostics) && rep_diag_found

            q_mean_i = mean(q_local_this_realization, 'omitnan');

            figs_diag = adaptive_req.figures.plot_step_diagnostics( ...
                sim, ...
                patch_pack, ...
                rep_diag, ...
                'StepIndex', step_idx, ...
                'Realization', r, ...
                'PatchIndex', step_diagnostic_patch_idx, ...
                'OmegaSr', schedule.Omega_sr(step_idx), ...
                'StepValue', step_value, ...
                'StepName', schedule.name, ...
                'StepUnit', schedule.unit, ...
                'QMean', q_mean_i, ...
                'KReal', 2*pi*cfg_i.f0/cfg_i.cs_bg, ...
                'ShowPatchCenters', true, ...
                'ShowSelectedPatchBox', true, ...
                'ShowPatchLabels', false, ...
                'Visible', step_diagnostic_visible);

            if save_step_diagnostics

                base_name = sprintf( ...
                    'step_%03d_realization_%03d_patch_%03d', ...
                    step_idx, ...
                    r, ...
                    step_diagnostic_patch_idx);

                adaptive_req.figures.save_figure_bundle( ...
                    figs_diag, ...
                    step_diagnostic_dir, ...
                    base_name, ...
                    'SavePNG', opt.SaveDiagnosticPNG, ...
                    'SavePDF', opt.SaveDiagnosticPDF, ...
                    'SaveFIG', opt.SaveDiagnosticFIG, ...
                    'Resolution', opt.DiagnosticResolution, ...
                    'CloseAfterSave', close_step_diagnostics_after_save);
            end

        elseif (plot_step_diagnostics || save_step_diagnostics) && ~rep_diag_found

            warning( ...
                'Step diagnostics requested, but diagnostic patch %d was not found.', ...
                step_diagnostic_patch_idx);

        end

        if verbose
            fprintf('  Realization %d / %d completed.\n', r, num_realizations);
        end
    end
end

%% Finalize table

rows = rows(1:row_idx);

T_raw = struct2table(rows);

sweep.n_rows = height(T_raw);

if verbose
    fprintf('\nAperture sweep completed.\n');
    fprintf('Generated %d rows.\n', height(T_raw));
end

end

%% ========================================================================
% Local helper functions
% ========================================================================

function tf = is_run_struct(S)

tf = isstruct(S) && ...
     isfield(S, 'cfg') && ...
     isfield(S, 'feat_cfg') && ...
     isfield(S, 'EXP') && ...
     isfield(S, 'PLOT') && ...
     isfield(S, 'SAVE');

end

function [cfg, feat_cfg, args] = unpack_run_for_aperture_sweep(RUN, extra_args)

cfg = RUN.cfg;
feat_cfg = RUN.feat_cfg;

args_default = { ...
    'SamplingMode', RUN.EXP.sampling_mode, ...
    'NumSteps', RUN.EXP.num_steps, ...
    'StepIndices', RUN.EXP.step_indices, ...
    'NumRealizations', RUN.EXP.num_realizations, ...
    'NumPatches', RUN.EXP.num_patches, ...
    'SeedBase', RUN.EXP.seed_base, ...
    'ReqOptions', RUN.req_options, ...
    'StoreWavefields', RUN.PLOT.store_wavefields, ...
    'StoreReqCurve', false, ...
    'StoreReqMapping', true, ...
    'ComputeGlobalReq', false, ...
    'StoreGlobalReqMapping', true, ...
    'StoreReqMetadata', true, ...
    'StoreFeatureStruct', false, ...
    'PlotStepDiagnostics', RUN.PLOT.show_step_diagnostics || RUN.PLOT.save_step_diagnostics, ...
    'SaveStepDiagnostics', RUN.PLOT.save_step_diagnostics, ...
    'StepDiagnosticPatchIndex', RUN.EXP.selected_patch, ...
    'StepDiagnosticDir', RUN.SAVE.step_diag_dir, ...
    'StepDiagnosticVisible', RUN.PLOT.step_diagnostic_visible, ...
    'CloseStepDiagnosticsAfterSave', RUN.PLOT.close_step_diagnostics_after_save, ...
    'SaveDiagnosticPNG', RUN.SAVE.save_png, ...
    'SaveDiagnosticPDF', RUN.SAVE.save_pdf, ...
    'SaveDiagnosticFIG', RUN.SAVE.save_fig, ...
    'DiagnosticResolution', RUN.SAVE.png_resolution, ...
    'Verbose', true};

if isfield(RUN, 'PatchOptions')
    args_default = merge_name_value_args( ...
        args_default, ...
        {'PatchOptions', RUN.PatchOptions});
end

if isfield(RUN, 'OUTPUT')

    if isfield(RUN.OUTPUT, 'store_req_curve')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'StoreReqCurve', RUN.OUTPUT.store_req_curve});
    end

    if isfield(RUN.OUTPUT, 'store_req_mapping')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'StoreReqMapping', RUN.OUTPUT.store_req_mapping});
    end

    if isfield(RUN.OUTPUT, 'compute_global_req')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'ComputeGlobalReq', RUN.OUTPUT.compute_global_req});
    end

    if isfield(RUN.OUTPUT, 'store_global_req_mapping')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'StoreGlobalReqMapping', RUN.OUTPUT.store_global_req_mapping});
    end

    if isfield(RUN.OUTPUT, 'store_req_metadata')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'StoreReqMetadata', RUN.OUTPUT.store_req_metadata});
    end

    if isfield(RUN.OUTPUT, 'store_feature_struct')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'StoreFeatureStruct', RUN.OUTPUT.store_feature_struct});
    end

    if isfield(RUN.OUTPUT, 'verbose')
        args_default = merge_name_value_args( ...
            args_default, ...
            {'Verbose', RUN.OUTPUT.verbose});
    end
end

args = merge_name_value_args(args_default, extra_args);

end

function args_out = merge_name_value_args(args_base, args_override)

if isempty(args_override)
    args_out = args_base;
    return;
end

if mod(numel(args_override), 2) ~= 0
    error('Optional arguments must be provided as name-value pairs.');
end

args_out = args_base;

for i = 1:2:numel(args_override)

    name_i = args_override{i};
    value_i = args_override{i + 1};

    if ~(ischar(name_i) || isstring(name_i))
        error('Argument names must be character vectors or strings.');
    end

    names_existing = args_out(1:2:end);
    match_idx = find(strcmpi(string(names_existing), string(name_i)), 1);

    if isempty(match_idx)
        args_out = [args_out, {char(name_i), value_i}]; %#ok<AGROW>
    else
        value_position = 2 * match_idx;
        args_out{value_position} = value_i;
    end
end

end

function rows = assign_prefixed_numeric_fields_to_row( ...
    rows, row_idx, values, prefix)

if ~isstruct(values)
    return;
end

names = fieldnames(values);

for i = 1:numel(names)

    value_i = values.(names{i});

    if isnumeric(value_i) && isscalar(value_i)
        rows(row_idx).(char(prefix + string(names{i}))) = double(value_i);
    elseif islogical(value_i) && isscalar(value_i)
        rows(row_idx).(char(prefix + string(names{i}))) = double(value_i);
    end
end

end

function val = get_optional_numeric(S, field_name, default_val)

if isstruct(S) && isfield(S, field_name)

    x = S.(field_name);

    if isnumeric(x) && isscalar(x)
        val = double(x);
    elseif islogical(x) && isscalar(x)
        val = double(x);
    else
        val = default_val;
    end
else
    val = default_val;
end

end

function val = get_optional_logical(S, field_name, default_val)

if isstruct(S) && isfield(S, field_name)
    x = S.(field_name);

    if islogical(x) && isscalar(x)
        val = x;
    elseif isnumeric(x) && isscalar(x)
        val = logical(x);
    else
        val = logical(default_val);
    end
else
    val = logical(default_val);
end

end

function txt = get_optional_string(S, field_name, default_val)

if isstruct(S) && isfield(S, field_name)

    x = S.(field_name);

    if isstring(x)
        txt = x(1);
    elseif ischar(x)
        txt = string(x);
    elseif iscategorical(x)
        txt = string(x);
    elseif isnumeric(x) && isscalar(x)
        txt = string(sprintf('%.12g', x));
    else
        txt = string(x);
    end
else
    txt = string(default_val);
end

end

function val = get_curve_field_numeric(curve, field_name)

if isstruct(curve) && isfield(curve, field_name)

    x = curve.(field_name);

    if isnumeric(x) && isscalar(x)
        val = double(x);
    elseif islogical(x) && isscalar(x)
        val = double(x);
    else
        val = NaN;
    end
else
    val = NaN;
end

end

function txt = get_curve_field_string(curve, field_name)

if isstruct(curve) && isfield(curve, field_name)

    x = curve.(field_name);

    if isstring(x)
        txt = x(1);
    elseif ischar(x)
        txt = string(x);
    elseif isnumeric(x) && isscalar(x)
        txt = string(sprintf('%.12g', x));
    elseif islogical(x) && isscalar(x)
        txt = string(logical(x));
    else
        txt = string(missing);
    end
else
    txt = string(missing);
end

end

function fn = get_feature_name(feature_names, idx)

if iscell(feature_names)
    fn = char(string(feature_names{idx}));
elseif isstring(feature_names)
    fn = char(feature_names(idx));
elseif ischar(feature_names)
    fn = char(string(feature_names));
else
    fn = char(string(feature_names(idx)));
end

end
