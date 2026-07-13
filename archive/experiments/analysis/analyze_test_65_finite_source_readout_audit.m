clear; close all; clc;

% Test 65: finite-source and readout geometry audit before adaptive REQ.
%
% This analysis lives in adaptive_req_local because it is a pre-REQ validation
% gate. It uses wave_sim_project as the simulation backend, but writes all
% diagnostic outputs under outputs/test_65_finite_source_readout_audit/.
%
% It does not run REQ and does not train any model.

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = struct();
CFG.Mode = string(getenv_default('ADAPTIVE_REQ_TEST65_MODE', 'quick')); % validate|quick|full
CFG.WaveSimPath = string(getenv_default('ADAPTIVE_REQ_TEST65_WAVE_SIM_PATH', ...
    '/Users/sara/Documents/wave_sim_project'));
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST65_SAVE_ALL_MAPS', true);

if exist(fullfile(CFG.WaveSimPath, 'src'), 'dir') ~= 7
    error('wave_sim_project src folder not found: %s', fullfile(CFG.WaveSimPath, 'src'));
end
addpath(fullfile(CFG.WaveSimPath, 'src'));

CFG.Lx = 0.050;
CFG.Lz = 0.050;
CFG.dx = 0.0002;
CFG.dz = 0.0002;
CFG.sourceX = -0.012;
CFG.sourceY = 0;
CFG.sourceZ = 0.025;
CFG.sourceNumPoints = 11;
CFG.sourceMotionAxis = [0 0 1];
CFG.measurementAxis = [0 0 1];

switch CFG.Mode
    case "validate"
        CFG.caseList = "homogeneous_cs2";
        CFG.f0List = 300;
        CFG.sourceAperturesMM = [0 8];
        CFG.realismList = ["clean", "readout_depth_plus_acoustic_shadow"];
        CFG.maxDetailedMaps = inf;
    case "quick"
        CFG.caseList = ["homogeneous_cs2", "inclusion_2_3", "bilayer_2_3"];
        CFG.f0List = 500;
        CFG.sourceAperturesMM = [0 4 8];
        CFG.realismList = ["clean", "readout_depth_plus_acoustic_shadow"];
        CFG.maxDetailedMaps = inf;
    otherwise
        CFG.caseList = ["homogeneous_cs2", "inclusion_2_3", "bilayer_2_3"];
        CFG.f0List = [300 500];
        CFG.sourceAperturesMM = [0 4 8];
        CFG.realismList = ["clean", "readout_depth_plus_acoustic_shadow"];
        CFG.maxDetailedMaps = inf;
end

OUT = struct();
OUT.root = fullfile(root_dir, 'outputs', ...
    'test_65_finite_source_readout_audit', char(CFG.Mode));
OUT.figures = fullfile(OUT.root, 'figures');
OUT.maps = fullfile(OUT.figures, 'maps_by_condition');
OUT.tables = fullfile(OUT.root, 'tables');
ensure_dir(OUT.figures); ensure_dir(OUT.maps); ensure_dir(OUT.tables);

fprintf('\nTest 65: finite-source and readout geometry audit before adaptive REQ\n');
fprintf('Mode: %s | wave_sim_project: %s\n', CFG.Mode, CFG.WaveSimPath);
fprintf('No REQ extraction and no model training.\n');

rows = {};
rowIdx = 0;
mapCount = 0;
nTotal = numel(CFG.caseList) * numel(CFG.f0List) * ...
    numel(CFG.sourceAperturesMM) * numel(CFG.realismList);
iCond = 0;

for ic = 1:numel(CFG.caseList)
    for if0 = 1:numel(CFG.f0List)
        for ia = 1:numel(CFG.sourceAperturesMM)
            for ir = 1:numel(CFG.realismList)
                iCond = iCond + 1;
                C = struct();
                C.case_id = CFG.caseList(ic);
                C.f0 = CFG.f0List(if0);
                C.aperture_mm = CFG.sourceAperturesMM(ia);
                C.realism = CFG.realismList(ir);
                C.key = sprintf('%s__f%d__src%gmm__%s', C.case_id, C.f0, ...
                    C.aperture_mm, C.realism);

                fprintf('[%d/%d] %s\n', iCond, nTotal, C.key);
                cfg = make_config(C, CFG);
                out = simcore.simulateEikonalVectorReadout2p5D(cfg);

                S = summarize_condition(out, cfg, C);
                rowIdx = rowIdx + 1;
                rows(rowIdx,:) = struct_to_row(S); %#ok<AGROW>

                if CFG.SaveAllMaps && mapCount < CFG.maxDetailedMaps
                    plot_condition_maps(out, cfg, C, S, OUT);
                    mapCount = mapCount + 1;
                end
            end
        end
    end
end

T = cell2table(rows, 'VariableNames', summary_names());
writetable(T, fullfile(OUT.tables, 'test65_finite_source_summary.csv'));
plot_summary(T, OUT);
save(fullfile(OUT.root, 'test65_config.mat'), 'CFG');

fprintf('\nTest 65 complete.\n');
fprintf('Conditions audited: %d\n', height(T));
fprintf('Tables: %s\n', OUT.tables);
fprintf('Figures: %s\n', OUT.figures);

%% Local functions

function cfg = make_config(C, CFG)
cfg = struct();
cfg.Lx = CFG.Lx;
cfg.Lz = CFG.Lz;
cfg.dx = CFG.dx;
cfg.dz = CFG.dz;
cfg.f0 = C.f0;
cfg.cs_bg = 2.0;
cfg.y0 = 0;

cfg.Measurement.Axis = CFG.measurementAxis;
cfg.Acquisition.Nt = 36;
cfg.Acquisition.PRF = 5000;

cfg.Eikonal.MaxIterations = 120;
cfg.Eikonal.Tolerance = 1e-8;
cfg.Eikonal.NormalizeBoundaryTime = true;
cfg.Eikonal.DecayPower = 0.5;
cfg.Eikonal.StorePerSourceFields = false;
cfg.Eikonal.UseAttenuation = false;
cfg.Eikonal.AttenuationMode = 'straight_ray_integral';

cfg.Medium = make_medium(C);

src = simcore.propagation.makeFiniteSource2p5D( ...
    [CFG.sourceX CFG.sourceY CFG.sourceZ], C.aperture_mm*1e-3, ...
    source_point_count(C.aperture_mm, CFG), ...
    'Axis', [0 0 1], ...
    'MotionAxis', CFG.sourceMotionAxis, ...
    'Amplitude', 1, ...
    'Phase', 0, ...
    'Normalize', true);

cfg.Sources.PositionXYZ = src.PositionXYZ;
cfg.Sources.Amplitude = src.Amplitude;
cfg.Sources.Phase = src.Phase;
cfg.Sources.MotionAxis = src.MotionAxis;
cfg.SourceDiagnostic = src;

cfg.Ultrasound.FrequencyMHz = 5;
cfg.Ultrasound.AcousticAttenuation.Enabled = false;
cfg.Ultrasound.ReadoutNoise.Enabled = false;

if C.realism == "readout_depth_plus_acoustic_shadow"
    % Shear attenuation is physical attenuation of the shear wave and is
    % resolved at cfg.f0 by makeSimpleMedium2D/resolveMaterialAtFrequency.
    cfg.Eikonal.UseAttenuation = true;

    cfg.Ultrasound.AcousticAttenuation.Enabled = true;
    cfg.Ultrasound.AcousticAttenuation.Mode = 'plane_wave_compounding';
    cfg.Ultrasound.AcousticAttenuation.FrequencyMHz = cfg.Ultrasound.FrequencyMHz;
    cfg.Ultrasound.AcousticAttenuation.SteeringAnglesDeg = [-8 0 8];
    cfg.Ultrasound.AcousticAttenuation.CombineMode = 'mean_amplitude';
    cfg.Ultrasound.AcousticAttenuation.TwoWay = true;
    cfg.Ultrasound.AcousticAttenuation.CouplingToReadoutNoise = true;
    cfg.Ultrasound.AcousticAttenuation.NoiseCouplingMode = 'snr_from_echo';
    cfg.Ultrasound.AcousticAttenuation.IncludeBackscatter = true;
    cfg.Ultrasound.AcousticAttenuation.ReferenceSNRdB = 45;
    cfg.Ultrasound.AcousticAttenuation.EchoReferenceAmplitude = 1;
    cfg.Ultrasound.AcousticAttenuation.MaxNoiseGain = 6;
    cfg.Ultrasound.AcousticAttenuation.MinNoiseGain = 0.25;

    cfg.Ultrasound.ReadoutNoise.Enabled = true;
    cfg.Ultrasound.ReadoutNoise.Model = 'boukraa_time_domain';
    cfg.Ultrasound.ReadoutNoise.Level = 'low';
    % In this realism level, depth degradation enters through the acoustic
    % echo/SNR map. Boukraa then converts lower SNR into tracking noise.
    cfg.Ultrasound.ReadoutNoise.DepthProfile = 'none';
    cfg.Ultrasound.ReadoutNoise.DepthExponent = 3;
    cfg.Ultrasound.ReadoutNoise.AmplitudeNormalization = 'max_abs';
    cfg.Ultrasound.ReadoutNoise.CompensateTemporalExtraction = false;
    cfg.Ultrasound.ReadoutNoise.FrequencyNoiseGain = 0.5;
    cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.Enabled = true;
    cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaX_mm = 0.45;
    cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaZ_mm = 0.45;
    cfg.Ultrasound.ReadoutNoise.Seed = deterministic_seed(C);
end
end

function n = source_point_count(apertureMM, CFG)
if apertureMM <= 0
    n = 1;
else
    n = CFG.sourceNumPoints;
end
end

function seed = deterministic_seed(C)
txt = sprintf('%s|%g|%g|%s', C.case_id, C.f0, C.aperture_mm, C.realism);
vals = double(char(txt));
seed = 1000 + mod(sum(vals .* (1:numel(vals))), 900000);
end

function M = make_medium(C)
soft = material(2.0, 4, 0.50, 1.00, 0);
hard = material(3.0, 30, 1.50, 1.15, 6);
M = struct();
M.Material = soft;
M.cs_bg = 2.0;

switch C.case_id
    case "homogeneous_cs2"
        M.Masks = {};
    case "inclusion_2_3"
        M.Masks = {struct('Type','circle', ...
            'Params', struct('Center', [0.026 0.025], 'Radius', 0.007, 'SigmaEdge', 0), ...
            'Material', hard)};
    case "bilayer_2_3"
        M.Masks = {struct('Type','bilayer', ...
            'Params', struct('Bi_Angle', 0, 'Bi_Offset', 0.026, 'SigmaEdge', 0), ...
            'Material', hard)};
end
end

function S = material(c, alpha_np_m, us_alpha, us_power, backscatter_db)
S = struct('Type','isotropic', 'c', c, 'Density', 1000, ...
    'Attenuation', struct('Model','power_law', 'alpha0', alpha_np_m, ...
        'Units','Np/m', 'fRef',400, 'exponent',1.0), ...
    'USAttenuation', struct('alpha0', us_alpha, 'power', us_power, ...
        'Units','dB/cm/MHz'), ...
    'Backscatter', struct('GainDB', backscatter_db));
end

function S = summarize_condition(out, cfg, C)
U = out.U;
Uc = out.U_clean;
noise = abs(U - Uc);
amp = abs(U);
phaseErr = angle(U .* conj(Uc)) ./ (2*pi);
roi = center_roi(out.x, out.z);

NW = noise_weights(out);
sourceMap = source_map_2d(out.x, out.z, cfg.SourceDiagnostic);

S = struct();
S.case_id = string(C.case_id);
S.f0 = C.f0;
S.aperture_mm = C.aperture_mm;
S.realism = string(C.realism);
S.num_source_points = size(cfg.Sources.PositionXYZ,1);
S.source_weight_sum = sum(cfg.Sources.Amplitude(:));
S.source_weight_peak = max(cfg.Sources.Amplitude(:));
S.mean_amp_roi = mean(amp(roi), 'omitnan');
S.std_amp_roi = std(amp(roi), 0, 'omitnan');
S.mean_noise_roi = mean(noise(roi), 'omitnan');
S.mean_phase_error_cycles_roi = mean(abs(phaseErr(roi)), 'omitnan');
S.mean_depth_weight_roi = mean(NW.depth(roi), 'omitnan');
S.mean_echo_snr_noise_weight_roi = mean(NW.acoustic(roi), 'omitnan');
S.mean_total_noise_weight_roi = mean(NW.total(roi), 'omitnan');
S.max_extra_shadow_loss_db = max(NW.extra_shadow_loss(:), [], 'omitnan');
S.max_source_map = max(sourceMap(:));
S.mean_abs_measurement_weight = source_measurement_weight(out);
S.output_key = string(C.key);
end

function names = summary_names()
names = {'case_id','f0','aperture_mm','realism','num_source_points', ...
    'source_weight_sum','source_weight_peak','mean_amp_roi','std_amp_roi', ...
    'mean_noise_roi','mean_phase_error_cycles_roi','mean_depth_weight_roi', ...
    'mean_echo_snr_noise_weight_roi','mean_total_noise_weight_roi', ...
    'max_extra_shadow_loss_db','max_source_map','mean_abs_measurement_weight', ...
    'output_key'};
end

function row = struct_to_row(S)
names = summary_names();
row = cell(1, numel(names));
for i = 1:numel(names)
    row{i} = S.(names{i});
end
end

function roi = center_roi(x, z)
[X,Z] = meshgrid(x, z);
roi = abs(X - mean(x)) <= 0.006 & abs(Z - mean(z)) <= 0.006;
end

function w = source_measurement_weight(out)
w = NaN;
try
    info = out.diag.SourceInfo;
    vals = nan(numel(info),1);
    for i = 1:numel(info)
        vals(i) = info{i}.MeanAbsMeasurementWeight;
    end
    w = mean(vals, 'omitnan');
catch
end
end

function NW = noise_weights(out)
nanMap = NaN(size(out.cs_map));
NW = struct('depth', nanMap, 'acoustic', nanMap, 'total', nanMap, ...
    'alpha_us', nanMap, 'total_loss', nanMap, 'background_loss', nanMap, ...
    'extra_shadow_loss', nanMap, 'echo_amplitude', nanMap, ...
    'relative_echo_amplitude', nanMap, 'tracking_snr_db', nanMap);
if ~isfield(out.diag, 'timeReadoutNoise')
    return
end
D = out.diag.timeReadoutNoise;
if isfield(D, 'depthWeight') && ~isempty(D.depthWeight), NW.depth = D.depthWeight; end
if isfield(D, 'acousticNoiseWeight') && ~isempty(D.acousticNoiseWeight), NW.acoustic = D.acousticNoiseWeight; end
if isfield(D, 'totalNoiseWeight') && ~isempty(D.totalNoiseWeight), NW.total = D.totalNoiseWeight; end
if isfield(D, 'acoustic') && isstruct(D.acoustic) && isfield(D.acoustic, 'maps')
    A = D.acoustic;
    if isfield(A.maps, 'alpha_dB_cm_zx'), NW.alpha_us = A.maps.alpha_dB_cm_zx; end
    if isfield(A, 'totalLossDB'), NW.total_loss = A.totalLossDB; end
    if isfield(A, 'backgroundLossDB'), NW.background_loss = A.backgroundLossDB; end
    if isfield(A, 'extraShadowLossDB'), NW.extra_shadow_loss = A.extraShadowLossDB; end
    if isfield(A, 'echoAmplitude'), NW.echo_amplitude = A.echoAmplitude; end
    if isfield(A, 'relativeEchoAmplitude'), NW.relative_echo_amplitude = A.relativeEchoAmplitude; end
    if isfield(A, 'trackingSNRdB'), NW.tracking_snr_db = A.trackingSNRdB; end
end
end

function plot_condition_maps(out, cfg, C, S, OUT)
NW = noise_weights(out);
sourceMap = source_map_2d(out.x, out.z, cfg.SourceDiagnostic);
phaseErr = angle(out.U .* conj(out.U_clean)) ./ (2*pi);
roi = center_roi(out.x, out.z);

fig = figure('Color','w', 'Position', [40 40 2050 1250], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 4, 5, 'TileSpacing','compact', 'Padding','compact');
plot_map(out.x, out.z, out.cs_map, 'True SWS', 'm/s', roi, cfg.SourceDiagnostic);
plot_map(out.x, out.z, sourceMap, 'Finite source weights', 'relative', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, abs(out.U_clean), 'Clean amplitude |U|', 'a.u.', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, angle(out.U_clean), 'Clean phase', 'rad', [], cfg.SourceDiagnostic);

plot_map(out.x, out.z, abs(out.U), 'Readout amplitude |U|', 'a.u.', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, angle(out.U), 'Readout phase', 'rad', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, abs(out.U - out.U_clean), 'Readout/noise difference', 'a.u.', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, phaseErr, 'Phase error vs clean', 'cycles', [], cfg.SourceDiagnostic);

plot_map(out.x, out.z, out.alpha_map, 'Shear attenuation alpha(f0)', 'Np/m', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.alpha_us, 'US attenuation alpha(fUS)', 'dB/cm', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.echo_amplitude, 'US echo amplitude', 'relative', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.tracking_snr_db, 'Tracking SNR from echo', 'dB', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.extra_shadow_loss, 'US extra shadow loss', 'dB', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.total_loss, 'US total loss', 'dB', [], cfg.SourceDiagnostic);

plot_map(out.x, out.z, NW.depth, 'Boukraa depth weight', 'relative', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.acoustic, 'Echo/SNR noise weight', 'relative', [], cfg.SourceDiagnostic);
plot_map(out.x, out.z, NW.total, 'Total Boukraa noise weight', 'relative', [], cfg.SourceDiagnostic);
nexttile; axis off;
text(0, 0.95, sprintf(['Case: %s\\nf0: %g Hz\\nAperture: %.1f mm (%d pts)\\nRealism: %s\\n', ...
    'ROI mean |phase error|: %.4f cycles\\nROI mean noise: %.4g\\nMax shadow: %.2f dB'], ...
    C.case_id, C.f0, C.aperture_mm, S.num_source_points, C.realism, ...
    S.mean_phase_error_cycles_roi, S.mean_noise_roi, S.max_extra_shadow_loss_db), ...
    'Interpreter','none', 'FontSize', 11);

sgtitle(sprintf('Test 65: %s, f=%g Hz, source %.1f mm, %s', ...
    C.case_id, C.f0, C.aperture_mm, C.realism), 'Interpreter','none');
file = fullfile(OUT.maps, sprintf('test65__%s.png', sanitize(C.key)));
exportgraphics(fig, file, 'Resolution', 170);
close(fig);

plot_profiles(out, NW, C, OUT);
end

function plot_profiles(out, NW, C, OUT)
ix = round(numel(out.x)/2);
iz = round(numel(out.z)/2);
fig = figure('Color','w', 'Position', [80 80 1700 980], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');

nexttile; hold on; grid on;
plot(1e3*out.x, abs(out.U_clean(iz,:)), 'LineWidth', 1.6);
plot(1e3*out.x, abs(out.U(iz,:)), 'LineWidth', 1.6);
xlabel('x (mm)'); ylabel('|U|'); title('Lateral amplitude profile');
legend({'clean','readout'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.z, abs(out.U_clean(:,ix)), 'LineWidth', 1.6);
plot(1e3*out.z, abs(out.U(:,ix)), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('|U|'); title('Axial/depth amplitude profile');
legend({'clean','readout'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.x, unwrap(angle(out.U_clean(iz,:))), 'LineWidth', 1.6);
plot(1e3*out.x, unwrap(angle(out.U(iz,:))), 'LineWidth', 1.6);
xlabel('x (mm)'); ylabel('phase (rad)'); title('Lateral phase profile');
legend({'clean','readout'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.z, unwrap(angle(out.U_clean(:,ix))), 'LineWidth', 1.6);
plot(1e3*out.z, unwrap(angle(out.U(:,ix))), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('phase (rad)'); title('Axial/depth phase profile');
legend({'clean','readout'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.z, NW.extra_shadow_loss(:,ix), 'LineWidth', 1.6);
plot(1e3*out.z, NW.total_loss(:,ix), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('loss (dB)'); title('US acoustic loss at center x');
legend({'extra shadow','total'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.z, NW.echo_amplitude(:,ix), 'LineWidth', 1.6);
plot(1e3*out.z, NW.relative_echo_amplitude(:,ix), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('relative amplitude'); title('US echo amplitude at center x');
legend({'absolute echo','relative to background'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*out.z, NW.tracking_snr_db(:,ix), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('SNR (dB)'); title('Tracking SNR at center x');

nexttile; hold on; grid on;
plot(1e3*out.z, NW.depth(:,ix), 'LineWidth', 1.6);
plot(1e3*out.z, NW.acoustic(:,ix), 'LineWidth', 1.6);
plot(1e3*out.z, NW.total(:,ix), 'LineWidth', 1.6);
xlabel('z (mm)'); ylabel('weight'); title('Boukraa noise weights at center x');
legend({'depth','echo/SNR','total'}, 'Location','best');

sgtitle(sprintf('Test 65 profiles: %s, f=%g Hz, source %.1f mm, %s', ...
    C.case_id, C.f0, C.aperture_mm, C.realism), 'Interpreter','none');
file = fullfile(OUT.maps, sprintf('test65_profiles__%s.png', sanitize(C.key)));
exportgraphics(fig, file, 'Resolution', 170);
close(fig);
end

function plot_map(x, z, A, ttl, cblabel, roi, src)
nexttile;
imagesc(1e3*x, 1e3*z, A);
axis image; set(gca,'YDir','normal');
title(ttl, 'Interpreter','none'); xlabel('x (mm)'); ylabel('z (mm)');
cb = colorbar; ylabel(cb, cblabel, 'Interpreter','none');
hold on;
if ~isempty(roi)
    [rr, cc] = find(roi);
    rectangle('Position', [1e3*x(min(cc)) 1e3*z(min(rr)) ...
        1e3*(x(max(cc))-x(min(cc))) 1e3*(z(max(rr))-z(min(rr)))], ...
        'EdgeColor', [1 1 1], 'LineStyle','--', 'LineWidth', 1.8);
end
if nargin >= 7 && ~isempty(src)
    pos = src.PositionXYZ;
    amp = src.Amplitude(:) ./ max(src.Amplitude(:));
    scatter(1e3*pos(:,1), 1e3*pos(:,3), 24 + 120*amp, ...
        'w', 'filled', 'MarkerEdgeColor','k');
end
end

function M = source_map_2d(x, z, src)
[X,Z] = meshgrid(x, z);
M = zeros(numel(z), numel(x));
sigma = 0.00045;
for i = 1:size(src.PositionXYZ,1)
    dx = X - src.PositionXYZ(i,1);
    dz = Z - src.PositionXYZ(i,3);
    M = M + src.Amplitude(i) * exp(-0.5*(dx.^2 + dz.^2) / sigma^2);
end
if max(M(:)) > 0
    M = M ./ max(M(:));
end
end

function plot_summary(T, OUT)
fig = figure('Color','w', 'Position', [90 90 1400 720], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile;
gscatter(T.aperture_mm, T.mean_phase_error_cycles_roi, T.realism);
grid on; xlabel('Source aperture (mm)'); ylabel('ROI mean |phase error| (cycles)');
title('Phase error vs source size');

nexttile;
gscatter(T.aperture_mm, T.mean_noise_roi, T.realism);
grid on; xlabel('Source aperture (mm)'); ylabel('ROI mean |U-U clean|');
title('Readout/noise difference vs source size');

nexttile;
gscatter(T.aperture_mm, T.mean_amp_roi, T.case_id);
grid on; xlabel('Source aperture (mm)'); ylabel('ROI mean amplitude');
title('Amplitude vs source size');

nexttile;
gscatter(T.aperture_mm, T.max_extra_shadow_loss_db, T.case_id);
grid on; xlabel('Source aperture (mm)'); ylabel('Max extra shadow (dB)');
title('Acoustic shadow by case');

sgtitle('Test 65 finite-source summary');
exportgraphics(fig, fullfile(OUT.figures, 'test65_summary_metrics.png'), 'Resolution', 180);
close(fig);
end

function ensure_dir(p)
if exist(p, 'dir') ~= 7
    mkdir(p);
end
end

function v = getenv_default(name, defaultVal)
v = getenv(name);
if isempty(v), v = defaultVal; end
end

function tf = env_true(name, defaultVal)
v = getenv(name);
if isempty(v), tf = defaultVal; return; end
tf = any(strcmpi(v, {'1','true','yes','on'}));
end

function s = sanitize(s)
s = char(string(s));
s = regexprep(s, '[^a-zA-Z0-9_\-]+', '_');
end
