%% test_50_controlled_kwave_phase1.m
% Phase 1 controlled k-Wave simulations for adaptive REQ.
%
% Purpose:
%   Create small, interpretable, reproducible k-Wave fields under our own
%   configuration before launching larger validations. The simulation code is
%   intentionally kept in tests/, while reusable k-Wave helpers live in
%   src/+adaptive_req/+kwave/.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST50_MODE = validate | quick | full
%   ADAPTIVE_REQ_KWAVE_PATH  = /Users/sara/Documents/k-wave-toolbox-version-1.4.1
%   ADAPTIVE_REQ_TEST50_RUN_REQ = true | false
%   ADAPTIVE_REQ_TEST50_TARGET_STEP_M = 0.001
%   ADAPTIVE_REQ_TEST50_USE_PARFOR = true | false
%   ADAPTIVE_REQ_TEST50_SAVE_TIME_SERIES = true | false
%   ADAPTIVE_REQ_TEST50_NX / ADAPTIVE_REQ_TEST50_NZ / ADAPTIVE_REQ_TEST50_T_END
%   ADAPTIVE_REQ_TEST50_VELOCITY_COMPONENT = axial_shear | lateral_shear | shear_magnitude
%   ADAPTIVE_REQ_TEST50_SOURCE_SIDE = left | right | top | bottom
%   ADAPTIVE_REQ_TEST50_SOURCE_POLARIZATION = axial | lateral | radial | transverse
%   ADAPTIVE_REQ_TEST50_ANALYSIS_ROI = exclude_source_buffer | full
%   ADAPTIVE_REQ_TEST50_ANALYSIS_BUFFER_M = 0.012
%   ADAPTIVE_REQ_TEST50_ANALYSIS_MARGIN_M = 0.002
%
% validate:
%   No k-Wave simulation. Checks toolbox path, material-map construction,
%   source configuration, and output folders.
%
% quick:
%   Small homogeneous and inclusion cases with one source, M=2.
%
% full:
%   Homogeneous 2/3, bilayer, inclusion; single and multi-source cases;
%   M=[2 3]. This is still a controlled phase-1 set, not a large sweep.

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));
run(fullfile(root_dir, 'setup_adaptive_req.m'));

CFG = struct();
CFG.Mode = lower(env_string('ADAPTIVE_REQ_TEST50_MODE', "validate"));
CFG.KWavePath = env_string('ADAPTIVE_REQ_KWAVE_PATH', ...
    "/Users/sara/Documents/k-wave-toolbox-version-1.4.1");
CFG.RunREQ = env_true('ADAPTIVE_REQ_TEST50_RUN_REQ', CFG.Mode ~= "validate");
CFG.TargetStepM = env_number('ADAPTIVE_REQ_TEST50_TARGET_STEP_M', 1.0e-3);
CFG.UseWindowParfor = env_true('ADAPTIVE_REQ_TEST50_USE_PARFOR', false);
CFG.SaveTimeSeries = env_true('ADAPTIVE_REQ_TEST50_SAVE_TIME_SERIES', false);
CFG.VelocityComponent = env_string('ADAPTIVE_REQ_TEST50_VELOCITY_COMPONENT', "axial_shear");
CFG.SourceSide = env_string('ADAPTIVE_REQ_TEST50_SOURCE_SIDE', "left");
CFG.SourcePolarization = env_string('ADAPTIVE_REQ_TEST50_SOURCE_POLARIZATION', "axial");
CFG.AnalysisROIMode = env_string('ADAPTIVE_REQ_TEST50_ANALYSIS_ROI', "exclude_source_buffer");
CFG.AnalysisBufferM = env_number('ADAPTIVE_REQ_TEST50_ANALYSIS_BUFFER_M', 12e-3);
CFG.AnalysisMarginM = env_number('ADAPTIVE_REQ_TEST50_ANALYSIS_MARGIN_M', 2e-3);
CFG.OverrideNx = env_number_or_nan('ADAPTIVE_REQ_TEST50_NX');
CFG.OverrideNz = env_number_or_nan('ADAPTIVE_REQ_TEST50_NZ');
CFG.OverrideTEnd = env_number_or_nan('ADAPTIVE_REQ_TEST50_T_END');
CFG.OutputRoot = fullfile(root_dir, 'outputs', 'test_50_controlled_kwave_phase1', char(CFG.Mode));
CFG.DataDir = fullfile(CFG.OutputRoot, 'data');
CFG.FigureDir = fullfile(CFG.OutputRoot, 'figures');
if ~exist(CFG.DataDir, 'dir'), mkdir(CFG.DataDir); end
if ~exist(CFG.FigureDir, 'dir'), mkdir(CFG.FigureDir); end

fprintf('\nTest 50: controlled k-Wave phase 1\n');
fprintf('Mode: %s | Run REQ: %d | target REQ step %.3f mm\n', ...
    CFG.Mode, CFG.RunREQ, 1e3 * CFG.TargetStepM);
fprintf('k-Wave path request: %s\n', CFG.KWavePath);

kwave_root = adaptive_req.kwave.locate_kwave_toolbox(CFG.KWavePath);
fprintf('k-Wave root: %s\n', kwave_root);
fprintf('pstdElastic2D: %s\n', which('pstdElastic2D'));

cases = build_cases(CFG.Mode);
cases = apply_case_overrides(cases, CFG);
M_list = build_M_list(CFG.Mode);

write_text_config(CFG, kwave_root);

if CFG.Mode == "validate"
    C0 = cases(1);
    sim_cfg = adaptive_req.kwave.default_controlled_config( ...
        'Geometry', C0.geometry, 'SourceMode', C0.source_mode, ...
        'Seed', C0.seed, 'KWavePath', CFG.KWavePath);
    MAT = adaptive_req.kwave.make_material_map_2d(sim_cfg);
    assert(isequal(size(MAT.cs_xz), [sim_cfg.Nx, sim_cfg.Nz]));
    assert(all(isfinite(MAT.cs_xz(:))));
    fprintf('Validation-only checks passed. No k-Wave time stepping was run.\n');
    fprintf('Outputs: %s\n', CFG.OutputRoot);
    return;
end

Tall = table();
n_total = numel(cases) * numel(M_list);
counter = 0;

for ci = 1:numel(cases)
    C = cases(ci);
    sim_cfg = adaptive_req.kwave.default_controlled_config( ...
        'Geometry', C.geometry, 'SourceMode', C.source_mode, ...
        'Seed', C.seed, 'KWavePath', CFG.KWavePath, ...
        'Nx', C.Nx, 'Nz', C.Nz, 'dx', C.dx, 'dz', C.dz, ...
        'f0', C.f0, 't_end', C.t_end, ...
        'CompressionMode', C.compression_mode, ...
        'alpha_shear', C.alpha_shear, ...
        'VelocityComponent', CFG.VelocityComponent, ...
        'SourceSide', CFG.SourceSide, ...
        'SourcePolarization', CFG.SourcePolarization, ...
        'AnalysisROIMode', CFG.AnalysisROIMode, ...
        'AnalysisBufferM', CFG.AnalysisBufferM, ...
        'AnalysisMarginM', CFG.AnalysisMarginM);

    key = sprintf('%s__%s__%s_%s__%s__f%g__dx%gum', ...
        sanitize(C.geometry), sanitize(C.source_mode), ...
        sanitize(CFG.SourceSide), sanitize(CFG.SourcePolarization), ...
        sanitize(CFG.AnalysisROIMode), ...
        C.f0, round(1e6*C.dx));
    fprintf('\nRunning k-Wave field: %s\n', key);
    t_sim = tic;
    S = adaptive_req.kwave.run_controlled_2d(sim_cfg);
    fprintf('  k-Wave completed in %.1f s | analysis field %s | full field %s | Nt=%d\n', ...
        toc(t_sim), mat2str(size(S.Uxz)), mat2str(size(S.full_Uxz)), numel(S.kgrid_t_array));

    save_field(CFG, key, S);
    plot_field_summary(S, CFG, key);

    for mi = 1:numel(M_list)
        counter = counter + 1;
        M = M_list(mi);
        fprintf('[%d/%d] REQ diagnostics for %s | M=%g\n', ...
            counter, n_total, key, M);

        if CFG.RunREQ
            T = extract_req_summary(S, sim_cfg, C, key, M, CFG);
            Tall = vertcat_compatible(Tall, T);
        else
            fprintf('  skipped REQ extraction because ADAPTIVE_REQ_TEST50_RUN_REQ=false.\n');
        end
    end
end

if ~isempty(Tall)
    writetable(Tall, fullfile(CFG.DataDir, 'test50_req_patch_summary.csv'));
end

fprintf('\nTest 50 complete.\nData: %s\nFigures: %s\n', CFG.DataDir, CFG.FigureDir);

%% Local helpers

function cases = build_cases(mode)
base = struct('geometry', "", 'source_mode', "", 'seed', 1001, ...
    'Nx', 96, 'Nz', 96, 'dx', 0.5e-3, 'dz', 0.5e-3, ...
    'f0', 500, 't_end', 0.040, 'compression_mode', "matched_shear", ...
    'alpha_shear', 20);

if mode == "quick"
    cases = repmat(base, 2, 1);
    cases(1).geometry = "homogeneous_cs2";
    cases(1).source_mode = "single_sine";
    cases(1).seed = 1101;
    cases(2).geometry = "inclusion_2_3";
    cases(2).source_mode = "single_sine";
    cases(2).seed = 1102;
else
    specs = [
        "homogeneous_cs2", "single_sine", 1201
        "homogeneous_cs3", "single_sine", 1202
        "bilayer_2_3", "single_sine", 1203
        "inclusion_2_3", "single_sine", 1204
        "inclusion_2_3", "sources8_sine", 1205
        "inclusion_2_3", "sources128_sine", 1206
        "inclusion_2_3", "single_square", 1207
        ];
    cases = repmat(base, size(specs, 1), 1);
    for i = 1:size(specs, 1)
        cases(i).geometry = specs(i, 1);
        cases(i).source_mode = specs(i, 2);
        cases(i).seed = str2double(specs(i, 3));
    end
end
end

function M_list = build_M_list(mode)
if mode == "quick"
    M_list = 2;
else
    M_list = [2 3];
end
end

function cases = apply_case_overrides(cases, CFG)
for i = 1:numel(cases)
    if isfinite(CFG.OverrideNx)
        cases(i).Nx = round(CFG.OverrideNx);
    end
    if isfinite(CFG.OverrideNz)
        cases(i).Nz = round(CFG.OverrideNz);
    end
    if isfinite(CFG.OverrideTEnd)
        cases(i).t_end = CFG.OverrideTEnd;
    end
end
end

function T = extract_req_summary(S, sim_cfg, C, key, M, CFG)
cfg_req = struct();
cfg_req.dx = sim_cfg.dx;
cfg_req.dz = sim_cfg.dz;
cfg_req.f0 = sim_cfg.f0;
cfg_req.cs_bg = sim_cfg.cs_soft;
cfg_req.WaveModel = char(C.source_mode);

feat = adaptive_req.config.default_feature_config( ...
    'M', M, 'cs_guess', 3.0, 'gamma_win', 1, 'pad_factor', 1);
step_x = max(1, round(CFG.TargetStepM / sim_cfg.dx));
step_z = max(1, round(CFG.TargetStepM / sim_cfg.dz));
req_options = {'Nbins', 'auto', 'Nbins_auto_oversample', 1, ...
    'Nbins_min', 16, 'smooth_sigma', 1};

t_req = tic;
O = adaptive_req.estimators.req_estimator_map(S.Uxz, cfg_req, feat, ...
    'StepX', step_x, 'StepZ', step_z, 'EdgeMode', 'valid', ...
    'QuantileMode', 'local_req', 'ReqOptions', req_options, ...
    'ReturnFeatures', false, 'ReturnFeatureTable', true, ...
    'StoreReqCurves', false, 'UseWindowParfor', CFG.UseWindowParfor, ...
    'Verbose', false);
fprintf('  REQ extracted %d windows in %.1f s | StepX=%d StepZ=%d\n', ...
    height(O.feature_table), toc(t_req), step_x, step_z);

T = O.feature_table;
T.condition_key = repmat(string(key), height(T), 1);
T.geometry = repmat(string(C.geometry), height(T), 1);
T.source_mode = repmat(string(C.source_mode), height(T), 1);
T.f0 = sim_cfg.f0 * ones(height(T), 1);
T.dx = sim_cfg.dx * ones(height(T), 1);
T.dz = sim_cfg.dz * ones(height(T), 1);
T.M = M * ones(height(T), 1);
T.REQ_StepX = step_x * ones(height(T), 1);
T.REQ_StepZ = step_z * ones(height(T), 1);
T.TargetStepM = CFG.TargetStepM * ones(height(T), 1);
T.true_SWS = sample_true_sws(S.cs_map, T.x_center_m, T.z_center_m, sim_cfg);
T.abs_error_local_req_pct = 100 * abs(T.cs_pred - T.true_SWS) ./ T.true_SWS;
T.signed_error_local_req_pct = 100 * (T.cs_pred - T.true_SWS) ./ T.true_SWS;
end

function true_sws = sample_true_sws(cs_map_zx, x_m, z_m, cfg)
ix = min(max(round(x_m ./ cfg.dx) + 1, 1), size(cs_map_zx, 2));
iz = min(max(round(z_m ./ cfg.dz) + 1, 1), size(cs_map_zx, 1));
idx = sub2ind(size(cs_map_zx), iz, ix);
true_sws = cs_map_zx(idx);
end

function save_field(CFG, key, S)
file = fullfile(CFG.DataDir, sprintf('%s_field.mat', key));
if CFG.SaveTimeSeries
    save(file, 'S', '-v7.3');
else
    Slim = rmfield_if_present(S, {'selected_shear_time_xzt'});
    save(file, 'Slim', '-v7.3');
end
end

function S = rmfield_if_present(S, names)
for i = 1:numel(names)
    if isfield(S, names{i})
        S = rmfield(S, names{i});
    end
end
end

function plot_field_summary(S, CFG, key)
fig = figure('Color', 'w', 'Position', [100 100 1280 720], 'Visible', 'off');
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, S.cs_map); axis image; set(gca,'YDir','normal');
title('True SWS'); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='m/s';

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, real(S.Uxz)); axis image; set(gca,'YDir','normal');
title(sprintf('Real harmonic field: %s', S.velocity_component)); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='a.u.';

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, abs(S.Uxz)); axis image; set(gca,'YDir','normal');
title('Harmonic amplitude'); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='a.u.';

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, angle(S.Uxz)); axis image; set(gca,'YDir','normal');
title('Harmonic phase'); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='rad';

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, S.source_mask_zx); axis image; set(gca,'YDir','normal');
title('Source mask'); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='source';

nexttile; imagesc(S.x_m*1e3, S.z_m*1e3, S.material_id); axis image; set(gca,'YDir','normal');
title('Material label'); xlabel('x (mm)'); ylabel('z (mm)'); cb=colorbar; cb.Label.String='label';

sgtitle(sprintf('Test 50 controlled k-Wave: %s', strrep(key, '_', '\_')));
exportgraphics(fig, fullfile(CFG.FigureDir, sprintf('%s_field_summary.png', key)), 'Resolution', 180);
close(fig);
end

function write_text_config(CFG, kwave_root)
fid = fopen(fullfile(CFG.OutputRoot, 'test50_configuration.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Test 50 controlled k-Wave phase 1\n');
fprintf(fid, 'Mode: %s\n', CFG.Mode);
fprintf(fid, 'k-Wave root: %s\n', kwave_root);
fprintf(fid, 'RunREQ: %d\n', CFG.RunREQ);
fprintf(fid, 'TargetStepM: %.9g\n', CFG.TargetStepM);
fprintf(fid, 'UseWindowParfor: %d\n', CFG.UseWindowParfor);
fprintf(fid, 'SaveTimeSeries: %d\n', CFG.SaveTimeSeries);
fprintf(fid, 'VelocityComponent: %s\n', CFG.VelocityComponent);
fprintf(fid, 'SourceSide: %s\n', CFG.SourceSide);
fprintf(fid, 'SourcePolarization: %s\n', CFG.SourcePolarization);
fprintf(fid, 'AnalysisROIMode: %s\n', CFG.AnalysisROIMode);
fprintf(fid, 'AnalysisBufferM: %.9g\n', CFG.AnalysisBufferM);
fprintf(fid, 'AnalysisMarginM: %.9g\n', CFG.AnalysisMarginM);
end

function out = env_string(name, default_val)
val = string(getenv(name));
if strlength(val) == 0
    out = string(default_val);
else
    out = val;
end
end

function out = env_number(name, default_val)
val = string(getenv(name));
if strlength(val) == 0
    out = default_val;
else
    out = str2double(val);
    if ~isfinite(out)
        error('Environment variable %s must be numeric.', name);
    end
end
end

function out = env_number_or_nan(name)
val = string(getenv(name));
if strlength(val) == 0
    out = NaN;
else
    out = str2double(val);
    if ~isfinite(out)
        error('Environment variable %s must be numeric.', name);
    end
end
end

function out = env_true(name, default_val)
val = lower(string(getenv(name)));
if strlength(val) == 0
    out = logical(default_val);
else
    out = any(val == ["1", "true", "yes", "y", "on"]);
end
end

function s = sanitize(x)
s = regexprep(lower(char(string(x))), '[^a-z0-9]+', '_');
s = regexprep(s, '^_|_$', '');
end

function T = vertcat_compatible(A, B)
if isempty(A)
    T = B;
elseif isempty(B)
    T = A;
else
    all_vars = unique([string(A.Properties.VariableNames), string(B.Properties.VariableNames)], 'stable');
    A = add_missing_vars(A, all_vars);
    B = add_missing_vars(B, all_vars);
    T = [A(:, all_vars); B(:, all_vars)];
end
end

function T = add_missing_vars(T, vars)
for v = vars(:).'
    if ~ismember(v, string(T.Properties.VariableNames))
        T.(v) = missing_column(height(T));
    end
end
end

function x = missing_column(n)
x = nan(n, 1);
end
