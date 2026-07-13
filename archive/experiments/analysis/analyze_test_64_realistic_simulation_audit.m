%% analyze_test_64_realistic_simulation_audit.m
% Test 64: audit wave_sim_project realistic simulations before REQ validation.
%
% This script does not train a q model and does not run REQ. It checks that
% the realistic simulation source is interpretable and reproducible before it
% becomes a validation set:
%   - homogeneous phase/k sanity checks;
%   - inclusion and bilayer mask/ROI checks;
%   - source/polarization/axial-readout diagnostics;
%   - clean/shear-attenuation/readout-depth/acoustic-shadow layer comparisons;
%   - compact outputs only, with time-series saved for diagnostic cases only.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST64_MODE          = validate | quick | full
%   ADAPTIVE_REQ_TEST64_WAVE_SIM_PATH = /Users/sara/Documents/wave_sim_project
%   ADAPTIVE_REQ_TEST64_SAVE_TIMESERIES_DIAG = true | false
%   ADAPTIVE_REQ_TEST64_SAVE_ALL_MAPS = true | false

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test64_configuration.json'));

fprintf('\nTest 64: realistic simulation audit before REQ validation\n');
fprintf('Mode: %s | wave_sim_project: %s\n', CFG.Mode, CFG.WaveSimPath);
fprintf('Primary simulator: simcore.simulateEikonalVectorReadout2p5D\n');
fprintf('Output used for REQ later: axial complex phasor out.U.\n');
fprintf('No q model training and no REQ extraction in this audit.\n');

assert(exist(fullfile(CFG.WaveSimPath, 'src'), 'dir') == 7, ...
    'Missing wave_sim_project src folder: %s', CFG.WaveSimPath);
addpath(fullfile(CFG.WaveSimPath, 'src'));

CONDITIONS = build_conditions(CFG);
fprintf('Audit conditions: %d\n', numel(CONDITIONS));

summary_parts = cell(numel(CONDITIONS), 1);
roi_parts = cell(numel(CONDITIONS), 1);
mask_parts = cell(numel(CONDITIONS), 1);
capability_parts = {};

for ci = 1:numel(CONDITIONS)
    C = CONDITIONS(ci);
    fprintf('[%d/%d] %s\n', ci, numel(CONDITIONS), C.condition_key);
    t0 = tic;
    [SIM, cfg_sim, capability] = run_condition(C, CFG, OUT);
    summary_parts{ci} = summarize_condition(SIM, cfg_sim, C, capability, CFG);
    [roi_parts{ci}, ROI] = audit_rois(SIM, C, CFG);
    mask_parts{ci} = audit_mask(SIM, C, CFG);
    capability_parts{end+1,1} = capability; %#ok<AGROW>
    if CFG.SaveAllMaps
        plot_audit_maps(SIM, C, ROI, capability, OUT, CFG);
    end
    save_compact_condition(SIM, cfg_sim, C, ROI, capability, OUT, CFG);
    fprintf('  done in %.1f s | simulator=%s | field=%s\n', ...
        toc(t0), capability.simulator_used, mat2str(size(SIM.U)));
end

T_summary = vertcat(summary_parts{:});
T_roi = vertcat(roi_parts{:});
T_mask = vertcat(mask_parts{:});
T_capability = struct2table(vertcat(capability_parts{:}));

writetable(T_summary, fullfile(OUT.table_dir, 'test64_simulation_audit_summary.csv'));
writetable(T_roi, fullfile(OUT.table_dir, 'test64_roi_checks.csv'));
writetable(T_mask, fullfile(OUT.table_dir, 'test64_mask_checks.csv'));
writetable(T_capability, fullfile(OUT.table_dir, 'test64_capability_gaps.csv'));

plot_summary_figures(T_summary, T_roi, T_mask, T_capability, OUT);

save(fullfile(OUT.data_dir, 'test64_compact_results.mat'), ...
    'CFG', 'T_summary', 'T_roi', 'T_mask', 'T_capability', '-v7.3');

print_console_summary(T_summary, T_roi, T_mask, T_capability, OUT);
fprintf('\nTables: %s\nFigures: %s\nTest 64 complete.\n', OUT.table_dir, OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST64_MODE'))));
if mode == "", mode = "quick"; end
assert(ismember(mode, ["validate","quick","full"]), ...
    'ADAPTIVE_REQ_TEST64_MODE must be validate, quick, or full.');

CFG = struct();
CFG.Mode = mode;
CFG.ValidateOnly = mode == "validate" || env_true('ADAPTIVE_REQ_TEST64_VALIDATE_ONLY', false);
CFG.QuickMode = mode == "quick";
CFG.FullMode = mode == "full";
CFG.WaveSimPath = char(env_string('ADAPTIVE_REQ_TEST64_WAVE_SIM_PATH', ...
    "/Users/sara/Documents/wave_sim_project"));
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST64_SAVE_ALL_MAPS', true);
CFG.SaveTimeSeriesDiag = env_true('ADAPTIVE_REQ_TEST64_SAVE_TIMESERIES_DIAG', false);
CFG.RandomSeed = 64001;

CFG.dx = 0.2e-3;
CFG.dz = 0.2e-3;
CFG.Lx = 0.05;
CFG.Lz = 0.05;
CFG.USFrequencyMHz = 7;
CFG.SourceXYZ = [-0.020 0.000 0.025];
CFG.MotionAxis = [0 0 1];
CFG.MeasurementAxis = [0 0 1];
CFG.AcquisitionPRF = 4000;
CFG.AcquisitionCycles = 4;

CFG.InclusionCenter = [0.025 0.025];
CFG.InclusionRadius = 0.008;
CFG.BilayerOffset = 0.025;
CFG.BilayerAngleRad = 0;
CFG.RoiSizeMm = 8;
CFG.InclusionHardRoiSizeMm = 5;
CFG.ThinRoiTolerance = 1e-12;
CFG.RoiPurityTolerance = 0.999;

if CFG.ValidateOnly
    CFG.dx = 0.5e-3;
    CFG.dz = 0.5e-3;
    CFG.Lx = 0.03;
    CFG.Lz = 0.03;
    CFG.SourceXYZ = [-0.015 0.000 0.015];
    CFG.InclusionCenter = [0.015 0.015];
    CFG.InclusionRadius = 0.005;
    CFG.BilayerOffset = 0.015;
    CFG.RoiSizeMm = 4;
    CFG.InclusionHardRoiSizeMm = 3;
end

CFG.OutputRoot = fullfile(root_dir, 'outputs', 'test_64_realistic_simulation_audit');
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.FullMode
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_64_realistic_simulation_audit');
else
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_64_realistic_simulation_audit', char(CFG.Mode));
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
OUT.condition_dir = fullfile(OUT.data_dir, 'conditions');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if exist(dirs(i), 'dir') ~= 7, mkdir(dirs(i)); end
end
end

function C = build_conditions(CFG)
cases = [
    case_spec("homogeneous_cs2", "homogeneous", 2, 2)
    case_spec("inclusion_2_3", "inclusion", 2, 3)
    case_spec("bilayer_2_3", "bilayer", 2, 3)
    ];
levels = ["clean","shear_attenuation","readout_depth_noise", ...
    "readout_depth_plus_acoustic_shadow"];
freqs = [300 500];
if CFG.FullMode
    freqs = [200 300 400 500 600];
end
if CFG.ValidateOnly
    cases = cases(1:2);
    levels = "clean";
    freqs = 300;
end

C = repmat(empty_condition(), 0, 1);
for ci = 1:numel(cases)
    for fi = 1:numel(freqs)
        for li = 1:numel(levels)
            X = empty_condition();
            X.case_id = cases(ci).case_id;
            X.case_family = cases(ci).case_family;
            X.cs_soft = cases(ci).cs_soft;
            X.cs_hard = cases(ci).cs_hard;
            X.f0 = freqs(fi);
            X.realism_level = levels(li);
            X.condition_key = string(sprintf('%s__f%g__%s', X.case_id, X.f0, X.realism_level));
            C(end+1,1) = X; %#ok<AGROW>
        end
    end
end
end

function C = case_spec(id, family, cs_soft, cs_hard)
C = struct('case_id', string(id), 'case_family', string(family), ...
    'cs_soft', cs_soft, 'cs_hard', cs_hard);
end

function X = empty_condition()
X = struct('case_id',"", 'case_family',"", 'cs_soft',NaN, 'cs_hard',NaN, ...
    'f0',NaN, 'realism_level',"", 'condition_key',"");
end

%% Simulation

function [SIM, cfg_sim, capability] = run_condition(C, CFG, OUT)
capability = base_capability(C);
switch C.case_family
    case {"homogeneous","inclusion","bilayer"}
        cfg_sim = build_vector_eikonal_cfg(C, CFG);
        out = simcore.simulateEikonalVectorReadout2p5D(cfg_sim);
        capability.simulator_used = "simulateEikonalVectorReadout2p5D";
        capability.vector_2p5d_supported = true;
    otherwise
        error('Unknown case family: %s', C.case_family);
end
SIM = standardize_output(out, C, CFG, capability);
if CFG.SaveTimeSeriesDiag
    file = fullfile(OUT.condition_dir, "timeseries_diag__" + sanitize(C.condition_key) + ".mat");
    if isfield(out, 'Ut_clean') || isfield(out, 'Ut_noisy')
        Ut_clean = getfield_or(out, 'Ut_clean', []);
        Ut_noisy = getfield_or(out, 'Ut_noisy', []);
        save(file, 'Ut_clean', 'Ut_noisy', '-v7.3');
    end
end
end

function cap = base_capability(C)
cap = struct();
cap.condition_key = C.condition_key;
cap.case_id = C.case_id;
cap.case_family = C.case_family;
cap.realism_level = C.realism_level;
cap.f0 = C.f0;
cap.primary_requested = "simulateEikonalVectorReadout2p5D";
cap.simulator_used = "";
cap.vector_2p5d_supported = false;
cap.axial_complex_output_available = true;
cap.time_series_saved = false;
cap.capability_gap = "";
end

function cfg = build_vector_eikonal_cfg(C, CFG)
cfg = struct();
cfg.Lx = CFG.Lx; cfg.Lz = CFG.Lz; cfg.dx = CFG.dx; cfg.dz = CFG.dz;
cfg.y0 = 0;
cfg.f0 = C.f0;
cfg.cs_bg = C.cs_soft;
cfg.Sources.PositionXYZ = CFG.SourceXYZ;
cfg.Sources.MotionAxis = CFG.MotionAxis;
cfg.Sources.Amplitude = 1;
cfg.Sources.Phase = 0;
cfg.Measurement.Axis = CFG.MeasurementAxis;
cfg.Medium = medium_config(C, CFG, "vector");
cfg.Eikonal.MaxIterations = 180;
cfg.Eikonal.Tolerance = 1e-8;
cfg.Eikonal.Verbose = false;
cfg.Eikonal.DecayPower = 0.5;
cfg.Eikonal.AttenuationMode = 'eikonal_ray_integral';
cfg.Eikonal.AttenuationSamples = 96;
cfg.Eikonal.RayPath.StepLength = 2*CFG.dx;
cfg.Eikonal.RayPath.MaxSteps = 1200;
cfg.Eikonal.Transport.Enabled = false;
cfg.Eikonal.StorePerSourceFields = true;
cfg.Shear.InterfaceRT.Enabled = false;
cfg.Shear.Scattering.Enabled = false;
cfg.Acquisition.PRF = max(CFG.AcquisitionPRF, 8*C.f0);
cfg.Acquisition.Nt = max(16, ceil(CFG.AcquisitionCycles * cfg.Acquisition.PRF / C.f0));
cfg.Acquisition.StoreTimeSeries = false;
cfg.Ultrasound.FrequencyMHz = CFG.USFrequencyMHz;
cfg = apply_realism_level(cfg, C, "vector");
end

function M = medium_config(C, CFG, ~)
M = struct();
M.cs_bg = C.cs_soft;
M.Material = material(C.cs_soft, 3, 0.50, 1.00, 0);
switch C.case_family
    case "homogeneous"
        M.Masks = {};
    case "inclusion"
        M.Masks = {struct('Type','circle', 'Params', struct( ...
            'Center', CFG.InclusionCenter, 'Radius', CFG.InclusionRadius, 'SigmaEdge', 0), ...
            'Material', material(C.cs_hard, 30, 1.50, 1.15, 6))};
    case "bilayer"
        M.Masks = {struct('Type','bilayer', 'Params', struct( ...
            'Bi_Angle', CFG.BilayerAngleRad, 'Bi_Offset', CFG.BilayerOffset, ...
            'SigmaEdge', 0), ...
            'Material', material(C.cs_hard, 30, 1.50, 1.15, 6))};
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

function cfg = apply_realism_level(cfg, C, simulator_family)
level = string(C.realism_level);
cfg.Eikonal.UseAttenuation = false;
cfg.Shear.Attenuation.Enabled = false;
cfg.Shear.Attenuation.Mode = 'ray_integral';
cfg.Shear.Attenuation.Model = 'none';
cfg.Shear.Attenuation.RaySamples = 48;
cfg.Ultrasound.AcousticAttenuation.Enabled = false;
cfg.Ultrasound.ReadoutNoise.Enabled = false;
switch string(level)
    case "clean"
        return;
    case {"shear_attenuation","attenuated_clean"}
        cfg = enable_shear_attenuation(cfg, simulator_family);
    case {"readout_depth_noise","readout_noise"}
        cfg = enable_shear_attenuation(cfg, simulator_family);
        cfg = enable_readout_noise(cfg, C, 'medium', 1.0, false);
    case {"readout_depth_plus_acoustic_shadow","shadow_noise"}
        cfg = enable_shear_attenuation(cfg, simulator_family);
        cfg = enable_readout_noise(cfg, C, 'high', 1.5, true);
        cfg = enable_acoustic_shadow(cfg);
    otherwise
        error('Unknown realism level: %s', level);
end
end

function cfg = enable_shear_attenuation(cfg, simulator_family)
if simulator_family == "vector"
    cfg.Eikonal.UseAttenuation = true;
else
    cfg.Shear.Attenuation.Enabled = true;
end
end

function cfg = enable_readout_noise(cfg, C, level, freq_gain, acoustic_shadow_enabled)
cfg.Ultrasound.ReadoutNoise.Enabled = true;
cfg.Ultrasound.ReadoutNoise.Model = 'boukraa_time_domain';
cfg.Ultrasound.ReadoutNoise.Level = level;
if acoustic_shadow_enabled
    % In the acoustic-shadow layer, depth/SNR degradation is computed from
    % the ultrasound echo-amplitude map. Boukraa converts that local SNR into
    % tracking noise, so we do not add a second phenomenological depth weight.
    cfg.Ultrasound.ReadoutNoise.DepthProfile = 'none';
else
    cfg.Ultrasound.ReadoutNoise.DepthProfile = 'exponential_0_to_1';
end
cfg.Ultrasound.ReadoutNoise.DepthExponent = 3;
cfg.Ultrasound.ReadoutNoise.AmplitudeNormalization = 'max_abs';
cfg.Ultrasound.ReadoutNoise.CompensateTemporalExtraction = false;
cfg.Ultrasound.ReadoutNoise.FrequencyNoiseGain = 0.5 * freq_gain;
cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.Enabled = true;
cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaX_mm = 0.45;
cfg.Ultrasound.ReadoutNoise.SpatialCorrelation.SigmaZ_mm = 0.45;
cfg.Ultrasound.ReadoutNoise.Seed = deterministic_seed(C, acoustic_shadow_enabled);
end

function seed = deterministic_seed(C, acoustic_shadow_enabled)
txt = sprintf('%s|%g|%s|shadow%d', ...
    C.case_id, C.f0, C.realism_level, acoustic_shadow_enabled);
vals = double(char(txt));
seed = 1000 + mod(sum(vals .* (1:numel(vals))), 900000);
end

function cfg = enable_acoustic_shadow(cfg)
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
end

function SIM = standardize_output(out, C, CFG, capability)
SIM = struct();
SIM.U = out.U;
SIM.U_clean = get_clean_field(out);
SIM.x = out.x(:).';
SIM.z = out.z(:);
SIM.cs_map = out.cs_map;
SIM.alpha_map = getfield_or(out, 'alpha_map', zeros(size(out.cs_map)));
SIM.amp_map = abs(SIM.U);
SIM.phase_map = angle(SIM.U);
SIM.phase_unwrapped = unwrap(unwrap(SIM.phase_map, [], 2), [], 1);
SIM.signed_distance_mm = signed_distance_map(SIM.x, SIM.z, C, CFG);
SIM.abs_distance_mm = abs(SIM.signed_distance_mm);
SIM.material_side = material_side_map(SIM.cs_map, C);
SIM.diag = getfield_or(out, 'diag', struct());
SIM.capability = capability;
end

function Uc = get_clean_field(out)
if isfield(out, 'U_clean') && ~isempty(out.U_clean)
    Uc = out.U_clean;
elseif isfield(out, 'U_clean_frequency') && ~isempty(out.U_clean_frequency)
    Uc = out.U_clean_frequency;
elseif isfield(out, 'U_recovered_clean') && ~isempty(out.U_recovered_clean)
    Uc = out.U_recovered_clean;
else
    Uc = out.U;
end
end

%% Audits

function T = summarize_condition(SIM, cfg, C, capability, CFG)
U = SIM.U; Uc = SIM.U_clean;
noise = U - Uc;
nsr = abs(noise) ./ (abs(Uc) + eps);
phase_err_cycles = angle(U .* conj(Uc)) ./ (2*pi);
[k_med, sws_phase_med] = phase_k_sanity(SIM, C, CFG);
[mw_mean, mw_max] = measurement_weight_summary(SIM);
T = table( ...
    string(C.condition_key), string(C.case_id), string(C.case_family), string(C.realism_level), C.f0, ...
    string(capability.simulator_used), capability.vector_2p5d_supported, ...
    size(U,1), size(U,2), cfg.dx, cfg.dz, ...
    mean(abs(U(:)), 'omitnan'), median(abs(U(:)), 'omitnan'), ...
    20*log10(norm(Uc(:)) / max(norm(noise(:)), eps)), ...
    median(20*log10(nsr(:) + eps), 'omitnan'), ...
    prctile(abs(phase_err_cycles(:)), 90), ...
    k_med, sws_phase_med, 100*(sws_phase_med - C.cs_soft)/C.cs_soft, ...
    mw_mean, mw_max, ...
    'VariableNames', {'condition_key','case_id','case_family','realism_level','f0', ...
    'simulator_used','vector_2p5d_supported','Nz','Nx','dx_m','dz_m', ...
    'mean_abs_U','median_abs_U','frequency_snr_db','median_NSR_db', ...
    'p90_phase_error_cycles','median_phase_k_rad_m','median_phase_sws_m_s', ...
    'homogeneous_phase_sws_error_pct','mean_abs_measurement_weight','max_abs_measurement_weight'});
end

function [k_med, sws_med] = phase_k_sanity(SIM, C, CFG)
k_med = NaN; sws_med = NaN;
if C.case_family ~= "homogeneous"
    return;
end
phi = SIM.phase_unwrapped;
dx = mean(diff(SIM.x)); dz = mean(diff(SIM.z));
[dphidz, dphidx] = gradient(phi, dz, dx);
k = hypot(dphidx, dphidz);
ROI = center_square_roi(SIM.x, SIM.z, CFG.RoiSizeMm*1e-3);
vals = k(ROI.mask);
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals), return; end
k_med = median(vals, 'omitnan');
sws_med = 2*pi*C.f0 / k_med;
end

function [mw_mean, mw_max] = measurement_weight_summary(SIM)
mw_mean = NaN; mw_max = NaN;
if isfield(SIM.diag, 'SourceInfo') && ~isempty(SIM.diag.SourceInfo)
    S = SIM.diag.SourceInfo;
    vals_mean = nan(numel(S),1);
    vals_max = nan(numel(S),1);
    for i = 1:numel(S)
        vals_mean(i) = getfield_or(S{i}, 'MeanAbsMeasurementWeight', NaN);
        vals_max(i) = getfield_or(S{i}, 'MaxAbsMeasurementWeight', NaN);
    end
    mw_mean = mean(vals_mean, 'omitnan');
    mw_max = max(vals_max, [], 'omitnan');
end
end

function [T, ROI] = audit_rois(SIM, C, CFG)
ROI = build_rois(SIM, C, CFG);
rows = cell(numel(ROI), 1);
for i = 1:numel(ROI)
    mask = ROI(i).mask;
    cs = SIM.cs_map(mask);
    side = SIM.material_side(mask);
    expected = ROI(i).expected_side;
    purity = mean(side == expected, 'omitnan');
    rows{i} = table(string(C.condition_key), string(C.case_id), string(C.case_family), string(C.realism_level), C.f0, ...
        string(ROI(i).name), string(expected), nnz(mask), ...
        ROI(i).width_mm, ROI(i).height_mm, ROI(i).center_x_mm, ROI(i).center_z_mm, ...
        mean(cs, 'omitnan'), std(cs, 'omitnan'), purity, purity >= CFG.RoiPurityTolerance, ...
        'VariableNames', {'condition_key','case_id','case_family','realism_level','f0', ...
        'roi_name','expected_side','N_pixels','width_mm','height_mm','center_x_mm','center_z_mm', ...
        'mean_true_sws','std_true_sws','material_purity','passes_purity_check'});
end
T = vertcat(rows{:});
end

function T = audit_mask(SIM, C, CFG)
cs = SIM.cs_map;
soft_mask = SIM.material_side == "soft";
hard_mask = SIM.material_side == "hard";
interface_mask = abs(SIM.signed_distance_mm) <= max(2*CFG.dx*1e3, 0.5);
speeds = unique(round(cs(:), 6));
expected_hard = C.case_family ~= "homogeneous";
T = table(string(C.condition_key), string(C.case_id), string(C.case_family), string(C.realism_level), C.f0, ...
    numel(speeds), min(cs(:)), max(cs(:)), ...
    mean(soft_mask(:)), mean(hard_mask(:)), mean(interface_mask(:)), ...
    expected_hard, any(hard_mask(:)) == expected_hard, ...
    mean(abs(cs(soft_mask) - C.cs_soft) < 1e-9, 'omitnan'), ...
    mean(abs(cs(hard_mask) - C.cs_hard) < 1e-9, 'omitnan'), ...
    'VariableNames', {'condition_key','case_id','case_family','realism_level','f0', ...
    'N_unique_speeds','min_sws','max_sws','soft_fraction','hard_fraction', ...
    'interface_band_fraction','expects_hard_region','passes_hard_presence_check', ...
    'soft_speed_match_fraction','hard_speed_match_fraction'});
end

function ROI = build_rois(SIM, C, CFG)
switch C.case_family
    case "homogeneous"
        R = center_square_roi(SIM.x, SIM.z, CFG.RoiSizeMm*1e-3);
        R.name = "center ROI"; R.expected_side = "soft";
        ROI = R;
    case "inclusion"
        hard_size = min(CFG.InclusionHardRoiSizeMm*1e-3, 1.4*CFG.InclusionRadius);
        R1 = rectangular_roi(SIM.x, SIM.z, CFG.InclusionCenter, [hard_size hard_size], "hard core ROI", "hard");
        bg_center = [max(min(SIM.x)+0.25*(max(SIM.x)-min(SIM.x)), 0.006), CFG.InclusionCenter(2)];
        bg_size = CFG.RoiSizeMm*1e-3;
        if hypot(bg_center(1)-CFG.InclusionCenter(1), bg_center(2)-CFG.InclusionCenter(2)) < CFG.InclusionRadius + 0.5*bg_size
            bg_center = [max(SIM.x)-0.25*(max(SIM.x)-min(SIM.x)), CFG.InclusionCenter(2)];
        end
        R2 = rectangular_roi(SIM.x, SIM.z, bg_center, [bg_size bg_size], "background soft ROI", "soft");
        ROI = [R1 R2];
    case "bilayer"
        size_m = CFG.RoiSizeMm*1e-3;
        x_soft = CFG.BilayerOffset - 0.5*size_m - 0.004;
        x_hard = CFG.BilayerOffset + 0.5*size_m + 0.004;
        z_mid = 0.5*(min(SIM.z)+max(SIM.z));
        R1 = rectangular_roi(SIM.x, SIM.z, [x_soft z_mid], [size_m size_m], "soft ROI", "soft");
        R2 = rectangular_roi(SIM.x, SIM.z, [x_hard z_mid], [size_m size_m], "hard ROI", "hard");
        ROI = [R1 R2];
    otherwise
        ROI = center_square_roi(SIM.x, SIM.z, CFG.RoiSizeMm*1e-3);
end
end

function R = center_square_roi(x, z, side_m)
R = rectangular_roi(x, z, [mean(x) mean(z)], [side_m side_m], "center ROI", "soft");
end

function R = rectangular_roi(x, z, center, size_m, name, side)
[X,Z] = meshgrid(x, z);
cx = center(1); cz = center(2);
wx = size_m(1); wz = size_m(2);
mask = abs(X - cx) <= wx/2 & abs(Z - cz) <= wz/2;
R = struct('name', string(name), 'expected_side', string(side), 'mask', mask, ...
    'center_x_mm', 1e3*cx, 'center_z_mm', 1e3*cz, ...
    'width_mm', 1e3*wx, 'height_mm', 1e3*wz, ...
    'x0_mm', 1e3*(cx-wx/2), 'z0_mm', 1e3*(cz-wz/2));
end

function S = signed_distance_map(x, z, C, CFG)
[X,Z] = meshgrid(x, z);
switch C.case_family
    case "homogeneous"
        S = NaN(size(X));
    case "inclusion"
        r = hypot(X - CFG.InclusionCenter(1), Z - CFG.InclusionCenter(2));
        S = 1e3*(CFG.InclusionRadius - r); % <0 background soft, >0 inclusion hard
    case "bilayer"
        n = [cos(CFG.BilayerAngleRad), sin(CFG.BilayerAngleRad)];
        S = 1e3*(X*n(1) + Z*n(2) - CFG.BilayerOffset); % <0 soft, >0 hard
    otherwise
        S = NaN(size(X));
end
end

function side = material_side_map(cs_map, C)
side = strings(size(cs_map));
if C.case_family == "homogeneous"
    side(:) = "soft";
    return;
end
mid = 0.5*(C.cs_soft + C.cs_hard);
side(cs_map < mid) = "soft";
side(cs_map >= mid) = "hard";
end

%% Plots and saving

function plot_audit_maps(SIM, C, ROI, capability, OUT, CFG)
fig = figure('Color','w','Position',[40 40 2050 1250], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 5, 5, 'TileSpacing','compact', 'Padding','compact');
plot_map(SIM.x, SIM.z, SIM.cs_map, 'True SWS', 'SWS (m/s)', ROI);
plot_map(SIM.x, SIM.z, SIM.amp_map, 'Axial amplitude |U|', 'amplitude', []);
plot_map(SIM.x, SIM.z, SIM.phase_map, 'Axial phase', 'phase (rad)', []);
plot_map(SIM.x, SIM.z, SIM.signed_distance_mm, 'Signed distance', 'distance (mm)', ROI);

plot_map(SIM.x, SIM.z, abs(SIM.U - SIM.U_clean), 'Readout/noise difference', '|U-U clean|', []);
phase_err = angle(SIM.U .* conj(SIM.U_clean)) ./ (2*pi);
plot_map(SIM.x, SIM.z, phase_err, 'Phase error vs clean', 'cycles', []);
plot_map(SIM.x, SIM.z, SIM.alpha_map, 'Shear attenuation map', 'Np/m', []);
plot_map(SIM.x, SIM.z, double(SIM.material_side == "hard"), 'Material mask', 'hard=1', ROI);

NW = noise_weight_maps(SIM);
plot_map(SIM.x, SIM.z, NW.alpha_us, 'US attenuation alpha(fUS)', 'dB/cm', []);
plot_map(SIM.x, SIM.z, NW.echo_amplitude, 'US echo amplitude', 'relative', []);
plot_map(SIM.x, SIM.z, NW.tracking_snr_db, 'Tracking SNR from echo', 'dB', []);
plot_map(SIM.x, SIM.z, NW.total_loss, 'US total acoustic loss', 'dB', []);
plot_map(SIM.x, SIM.z, NW.background_loss, 'US background loss', 'dB', []);
plot_map(SIM.x, SIM.z, NW.extra_shadow_loss, 'US extra shadow loss', 'dB', []);

plot_map(SIM.x, SIM.z, NW.depth, 'Boukraa depth weight', 'relative weight', []);
plot_map(SIM.x, SIM.z, NW.acoustic, 'Echo/SNR noise weight', 'relative weight', []);
plot_map(SIM.x, SIM.z, NW.total, 'Total Boukraa noise weight', 'relative weight', []);
plot_map(SIM.x, SIM.z, SIM.phase_unwrapped, 'Unwrapped phase', 'phase (rad)', []);
axis off; text(0, 0.9, sprintf('Simulator: %s\\nVector 2.5D supported: %d\\nGap: %s', ...
    capability.simulator_used, capability.vector_2p5d_supported, capability.capability_gap), ...
    'Interpreter','none', 'FontSize',10);

sgtitle(sprintf('Test 64 audit: %s, f=%g Hz, %s', ...
    C.case_id, C.f0, C.realism_level), 'Interpreter','none');
case_dir = fullfile(OUT.map_dir, char(C.case_id));
if exist(case_dir, 'dir') ~= 7, mkdir(case_dir); end
file = fullfile(case_dir, "test64_audit__" + sanitize(C.condition_key) + ".png");
exportgraphics(fig, file, 'Resolution', 180);
close(fig);

if CFG.SaveAllMaps
    plot_roi_overlay(SIM, C, ROI, OUT);
end
end

function NW = noise_weight_maps(SIM)
nan_map = NaN(size(SIM.cs_map));
NW = struct('depth', nan_map, 'acoustic', nan_map, 'total', nan_map, ...
    'alpha_us', nan_map, 'total_loss', nan_map, 'background_loss', nan_map, ...
    'extra_shadow_loss', nan_map, 'echo_amplitude', nan_map, ...
    'relative_echo_amplitude', nan_map, 'tracking_snr_db', nan_map);

if isfield(SIM.diag, 'timeReadoutNoise')
    D = SIM.diag.timeReadoutNoise;
elseif isfield(SIM.diag, 'readoutNoise')
    D = SIM.diag.readoutNoise;
else
    return;
end

if isfield(D, 'depthWeight') && ~isempty(D.depthWeight)
    NW.depth = D.depthWeight;
end
if isfield(D, 'acousticNoiseWeight') && ~isempty(D.acousticNoiseWeight)
    NW.acoustic = D.acousticNoiseWeight;
end
if isfield(D, 'totalNoiseWeight') && ~isempty(D.totalNoiseWeight)
    NW.total = D.totalNoiseWeight;
end
if isfield(D, 'acoustic') && isstruct(D.acoustic) && isfield(D.acoustic, 'maps')
    A = D.acoustic;
    if isfield(A.maps, 'alpha_dB_cm_zx') && ~isempty(A.maps.alpha_dB_cm_zx)
        NW.alpha_us = A.maps.alpha_dB_cm_zx;
    elseif isfield(A.maps, 'alpha_dB_cm_MHz_zx') && ~isempty(A.maps.alpha_dB_cm_MHz_zx)
        NW.alpha_us = A.maps.alpha_dB_cm_MHz_zx;
    end
    if isfield(A, 'totalLossDB') && ~isempty(A.totalLossDB)
        NW.total_loss = A.totalLossDB;
    elseif isfield(A, 'lossDB') && ~isempty(A.lossDB)
        NW.total_loss = A.lossDB;
    end
    if isfield(A, 'backgroundLossDB') && ~isempty(A.backgroundLossDB)
        NW.background_loss = A.backgroundLossDB;
    end
    if isfield(A, 'extraShadowLossDB') && ~isempty(A.extraShadowLossDB)
        NW.extra_shadow_loss = A.extraShadowLossDB;
    elseif isfield(A, 'relativeLossDB') && ~isempty(A.relativeLossDB)
        NW.extra_shadow_loss = A.relativeLossDB;
    end
    if isfield(A, 'echoAmplitude') && ~isempty(A.echoAmplitude)
        NW.echo_amplitude = A.echoAmplitude;
    end
    if isfield(A, 'relativeEchoAmplitude') && ~isempty(A.relativeEchoAmplitude)
        NW.relative_echo_amplitude = A.relativeEchoAmplitude;
    end
    if isfield(A, 'trackingSNRdB') && ~isempty(A.trackingSNRdB)
        NW.tracking_snr_db = A.trackingSNRdB;
    end
end
end

function plot_map(x, z, A, ttl, cblabel, ROI)
nexttile;
imagesc(1e3*x, 1e3*z, A);
axis image; set(gca,'YDir','normal');
title(ttl, 'Interpreter','none');
xlabel('x (mm)'); ylabel('z (mm)');
cb = colorbar; ylabel(cb, cblabel, 'Interpreter','none');
if ~isempty(ROI)
    hold on;
    draw_roi_boxes(ROI);
end
end

function draw_roi_boxes(ROI)
for i = 1:numel(ROI)
    color = [0 0.35 1];
    if ROI(i).expected_side == "hard", color = [1 0 0]; end
    rectangle('Position', [ROI(i).x0_mm ROI(i).z0_mm ROI(i).width_mm ROI(i).height_mm], ...
        'EdgeColor', color, 'LineWidth', 2, 'LineStyle','--');
    text(ROI(i).center_x_mm, ROI(i).center_z_mm, ROI(i).name, ...
        'Color','w', 'FontWeight','bold', 'HorizontalAlignment','center', ...
        'BackgroundColor',[0 0 0 0.4], 'Interpreter','none');
end
end

function plot_roi_overlay(SIM, C, ROI, OUT)
fig = figure('Color','w','Position',[80 80 900 700]);
imagesc(1e3*SIM.x, 1e3*SIM.z, SIM.cs_map);
axis image; set(gca,'YDir','normal'); hold on;
draw_roi_boxes(ROI);
xlabel('x (mm)'); ylabel('z (mm)');
title(sprintf('%s ROI/mask QA, f=%g Hz, %s', C.case_id, C.f0, C.realism_level), ...
    'Interpreter','none');
cb = colorbar; ylabel(cb, 'True SWS (m/s)');
case_dir = fullfile(OUT.figure_dir, 'roi_overlays', char(C.case_id));
if exist(case_dir, 'dir') ~= 7, mkdir(case_dir); end
file = fullfile(case_dir, "test64_roi_overlay__" + sanitize(C.condition_key) + ".png");
exportgraphics(fig, file, 'Resolution', 180);
close(fig);
end

function plot_summary_figures(T_summary, T_roi, T_mask, T_cap, OUT)
fig = figure('Color','w','Position',[100 100 1300 650]);
tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
nexttile;
gscatter(categorical(T_summary.realism_level), T_summary.p90_phase_error_cycles, categorical(T_summary.case_family));
grid on; ylabel('p90 phase error (cycles)'); title('Readout phase error by realism level');
nexttile;
roi_names = unique(string(T_roi.roi_name), 'stable');
roi_mean = nan(numel(roi_names), 1);
for i = 1:numel(roi_names)
    roi_mean(i) = mean(T_roi.material_purity(string(T_roi.roi_name) == roi_names(i)), 'omitnan');
end
bar(categorical(roi_names), roi_mean); grid on;
ylabel('Mean ROI material purity'); title('ROI purity checks'); ylim([0 1.05]);
nexttile;
case_names = unique(string(T_mask.case_id), 'stable');
hard_mean = nan(numel(case_names), 1);
for i = 1:numel(case_names)
    hard_mean(i) = mean(T_mask.hard_fraction(string(T_mask.case_id) == case_names(i)), 'omitnan');
end
bar(categorical(case_names), hard_mean); grid on;
ylabel('Mean hard material fraction'); title('Mask hard fraction by case');
nexttile;
case_names_cap = unique(string(T_cap.case_id), 'stable');
support_mean = nan(numel(case_names_cap), 1);
for i = 1:numel(case_names_cap)
    support_mean(i) = mean(double(T_cap.vector_2p5d_supported(string(T_cap.case_id) == case_names_cap(i))), 'omitnan');
end
bar(categorical(case_names_cap), support_mean); grid on;
ylabel('Fraction using vector 2.5D'); title('Capability support');
file = fullfile(OUT.figure_dir, 'test64_audit_summary.png');
exportgraphics(fig, file, 'Resolution', 180);
close(fig);
end

function save_compact_condition(SIM, cfg_sim, C, ROI, capability, OUT, CFG)
S = struct();
S.U = single(SIM.U);
S.U_clean = single(SIM.U_clean);
S.x = SIM.x;
S.z = SIM.z;
S.cs_map = single(SIM.cs_map);
S.alpha_map = single(SIM.alpha_map);
S.signed_distance_mm = single(SIM.signed_distance_mm);
S.material_side = SIM.material_side;
S.condition = C;
S.capability = capability;
S.roi = rmfield(ROI, 'mask');
S.cfg_sim = compact_config(cfg_sim);
if CFG.SaveTimeSeriesDiag
    S.note = 'Time series saved separately only when requested.';
else
    S.note = 'Time series intentionally omitted; this is a compact audit output.';
end
file = fullfile(OUT.condition_dir, "compact__" + sanitize(C.condition_key) + ".mat");
save(file, 'S', '-v7.3');
end

function cfg2 = compact_config(cfg)
cfg2 = cfg;
if isfield(cfg2, 'Acquisition')
    cfg2.Acquisition.StoreTimeSeries = false;
end
end

%% Console

function print_console_summary(T_summary, T_roi, T_mask, T_cap, OUT)
fprintf('\nInterpretive simulation-audit summary:\n');
bad_roi = T_roi(T_roi.passes_purity_check == false, :);
bad_mask = T_mask(T_mask.passes_hard_presence_check == false, :);
fprintf('  Conditions audited: %d.\n', height(T_summary));
fprintf('  ROI purity failures: %d.\n', height(bad_roi));
fprintf('  Mask hard-region presence failures: %d.\n', height(bad_mask));
fprintf('  Vector 2.5D unsupported cases: %d.\n', sum(~T_cap.vector_2p5d_supported));
H = T_summary(T_summary.case_family == "homogeneous", :);
if ~isempty(H)
    fprintf('  Homogeneous median phase-SWS error: %.2f%%.\n', ...
        median(abs(H.homogeneous_phase_sws_error_pct), 'omitnan'));
end
if any(~T_cap.vector_2p5d_supported)
    fprintf('  Main gap: bilayer needs primary vector-2.5D medium support before realistic bilayer validation.\n');
    fprintf('  Recommended next step: use passing inclusion/homogeneous cases for REQ readout stress-test, and add bilayer support or keep bilayer in 2D fallback until fixed.\n');
else
    fprintf('  Primary vector 2.5D path supports all audited geometries.\n');
    fprintf('  Recommended next step: use homogeneous, inclusion, and bilayer cases for REQ readout stress-test.\n');
end
fprintf('  Outputs under: %s\n', OUT.root_dir);
end

%% Small utilities

function val = getfield_or(S, field, default_val)
if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
    val = S.(field);
else
    val = default_val;
end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file, 'w');
assert(fid > 0, 'Could not write %s', file);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function tf = env_true(name, default_val)
v = lower(strtrim(string(getenv(name))));
if v == ""
    tf = default_val;
else
    tf = ismember(v, ["1","true","yes","y","on"]);
end
end

function v = env_string(name, default_val)
s = string(getenv(name));
if strlength(strtrim(s)) == 0
    v = string(default_val);
else
    v = s;
end
end

function s = sanitize(x)
s = regexprep(char(string(x)), '[^A-Za-z0-9_\-]+', '_');
end
