%% analyze_test_66_shear_shadow_source_integral_comparison.m
% Test 66: compare finite-source and attenuation-integral choices.
%
% Motivation:
% Test 65 can show a double-looking shear-amplitude shadow behind an
% inclusion. This diagnostic isolates whether that structure comes from:
%   1) point vs finite source aperture;
%   2) straight-ray vs eikonal-ray shear attenuation integration.
%
% This test does not run REQ and does not train/apply q models. It only
% generates clean shear fields and compares each inclusion field against its
% matched homogeneous reference.

clear; close all; clc;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = default_config(root_dir);
OUT = make_output_dirs(root_dir, CFG);
write_config_json(CFG, fullfile(OUT.root_dir, 'test66_configuration.json'));

fprintf('\nTest 66: shear shadow source/integral comparison\n');
fprintf('Mode: %s | wave_sim_project: %s\n', CFG.Mode, CFG.WaveSimPath);
fprintf('No REQ extraction. No q/confidence model training or inference.\n');

rows = {};
rowIdx = 0;
nTotal = numel(CFG.FrequenciesHz) * numel(CFG.SourceAperturesMM) * numel(CFG.AttenuationModes);
iCond = 0;

allResults = {};

for if0 = 1:numel(CFG.FrequenciesHz)
    for ia = 1:numel(CFG.SourceAperturesMM)
        for im = 1:numel(CFG.AttenuationModes)
            iCond = iCond + 1;
            C = struct();
            C.f0 = CFG.FrequenciesHz(if0);
            C.aperture_mm = CFG.SourceAperturesMM(ia);
            C.attenuation_mode = CFG.AttenuationModes(im);
            C.key = sprintf('f%d__src%gmm__%s', C.f0, C.aperture_mm, C.attenuation_mode);

            fprintf('[%d/%d] %s\n', iCond, nTotal, C.key);
            ticCond = tic;

            cfgIncl = make_config(C, CFG, "inclusion");
            cfgHom = make_config(C, CFG, "homogeneous");

            outIncl = simcore.simulateEikonalVector2p5D(cfgIncl);
            outHom = simcore.simulateEikonalVector2p5D(cfgHom);

            R = analyze_pair(outIncl, outHom, cfgIncl, C, CFG);
            R.elapsed_s = toc(ticCond);
            allResults{end+1} = R; %#ok<SAGROW>

            S = summary_row(R, C);
            rowIdx = rowIdx + 1;
            rows(rowIdx,:) = struct2row(S, summary_names()); %#ok<SAGROW>

            if CFG.SaveAllMaps
                plot_condition_pair(R, C, OUT);
            end
        end
    end
end

T = cell2table(rows, 'VariableNames', summary_names());
writetable(T, fullfile(OUT.table_dir, 'test66_shadow_comparison_summary.csv'));

plot_summary_grid([allResults{:}], OUT, CFG);

save(fullfile(OUT.data_dir, 'test66_compact_results.mat'), 'CFG', 'T', 'allResults', '-v7.3');

fprintf('\nTest 66 complete.\n');
fprintf('Conditions: %d\n', height(T));
fprintf('Tables: %s\n', OUT.table_dir);
fprintf('Figures: %s\n', OUT.figure_dir);

%% Configuration

function CFG = default_config(root_dir)
mode = lower(string(getenv_default('ADAPTIVE_REQ_TEST66_MODE', 'quick')));
CFG = struct();
CFG.Mode = mode;
CFG.WaveSimPath = char(env_string('ADAPTIVE_REQ_TEST66_WAVE_SIM_PATH', ...
    "/Users/sara/Documents/wave_sim_project"));
CFG.SaveAllMaps = env_true('ADAPTIVE_REQ_TEST66_SAVE_ALL_MAPS', true);

if exist(fullfile(CFG.WaveSimPath, 'src'), 'dir') ~= 7
    error('wave_sim_project src folder not found: %s', fullfile(CFG.WaveSimPath, 'src'));
end
addpath(fullfile(CFG.WaveSimPath, 'src'));

CFG.Lx = 0.050;
CFG.Lz = 0.050;
CFG.dx = 0.2e-3;
CFG.dz = 0.2e-3;
CFG.SourceXYZ = [-0.012 0 0.025];
CFG.SourceNumPoints = 11;
CFG.SourceAxis = [0 0 1];
CFG.MotionAxis = [0 0 1];
CFG.MeasurementAxis = [0 0 1];
CFG.DecayPower = 0.5;
CFG.InclusionCenter = [0.026 0.025];
CFG.InclusionRadius = 0.007;
CFG.FrequenciesHz = 500;
CFG.SourceAperturesMM = [0 4];
CFG.AttenuationModes = ["straight_ray_integral", "eikonal_ray_integral"];
CFG.ShadowThresholdDB = 3;

if mode == "validate"
    CFG.dx = 0.5e-3;
    CFG.dz = 0.5e-3;
elseif mode == "full"
    CFG.FrequenciesHz = [300 500];
    CFG.SourceAperturesMM = [0 4 8];
end

CFG.OutputRoot = fullfile(root_dir, 'outputs', 'test_66_shear_shadow_source_integral_comparison');
end

function OUT = make_output_dirs(root_dir, CFG)
if CFG.Mode == "full"
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_66_shear_shadow_source_integral_comparison');
else
    OUT.root_dir = fullfile(root_dir, 'outputs', 'test_66_shear_shadow_source_integral_comparison', char(CFG.Mode));
end
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.map_dir = fullfile(OUT.figure_dir, 'maps_by_condition');
OUT.data_dir = fullfile(OUT.root_dir, 'data');
dirs = string(struct2cell(OUT));
for i = 1:numel(dirs)
    if exist(dirs(i), 'dir') ~= 7
        mkdir(dirs(i));
    end
end
end

function cfg = make_config(C, CFG, caseKind)
cfg = struct();
cfg.Lx = CFG.Lx;
cfg.Lz = CFG.Lz;
cfg.dx = CFG.dx;
cfg.dz = CFG.dz;
cfg.f0 = C.f0;
cfg.cs_bg = 2.0;
cfg.y0 = 0;
cfg.Measurement.Axis = CFG.MeasurementAxis;

cfg.Eikonal.MaxIterations = 180;
cfg.Eikonal.Tolerance = 1e-8;
cfg.Eikonal.Verbose = false;
cfg.Eikonal.NormalizeBoundaryTime = true;
cfg.Eikonal.DecayPower = CFG.DecayPower;
cfg.Eikonal.StorePerSourceFields = true;
cfg.Eikonal.UseAttenuation = true;
cfg.Eikonal.AttenuationMode = char(C.attenuation_mode);
cfg.Eikonal.AttenuationSamples = 96;
cfg.Eikonal.RayPath.StepLength = 2*CFG.dx;
cfg.Eikonal.RayPath.MaxSteps = 1200;
cfg.Eikonal.Transport.Enabled = false;

cfg.Medium = make_medium(caseKind, CFG);

src = simcore.propagation.makeFiniteSource2p5D( ...
    CFG.SourceXYZ, C.aperture_mm*1e-3, source_point_count(C.aperture_mm, CFG), ...
    'Axis', CFG.SourceAxis, ...
    'MotionAxis', CFG.MotionAxis, ...
    'Amplitude', 1, ...
    'Phase', 0, ...
    'Normalize', true);

cfg.Sources.PositionXYZ = src.PositionXYZ;
cfg.Sources.Amplitude = src.Amplitude;
cfg.Sources.Phase = src.Phase;
cfg.Sources.MotionAxis = src.MotionAxis;
cfg.SourceDiagnostic = src;
end

function n = source_point_count(apertureMM, CFG)
if apertureMM <= 0
    n = 1;
else
    n = CFG.SourceNumPoints;
end
end

function M = make_medium(caseKind, CFG)
soft = material(2.0, 4, 0.50, 1.00, 0);
hard = material(3.0, 30, 1.50, 1.15, 6);
M = struct();
M.Material = soft;
M.cs_bg = 2.0;
if caseKind == "inclusion"
    M.Masks = {struct('Type','circle', ...
        'Params', struct('Center', CFG.InclusionCenter, ...
        'Radius', CFG.InclusionRadius, 'SigmaEdge', 0), ...
        'Material', hard)};
else
    M.Masks = {};
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

%% Analysis

function R = analyze_pair(outIncl, outHom, cfg, C, CFG)
ampIncl = abs(outIncl.U);
ampHom = abs(outHom.U);

% Ratio against a matched homogeneous reference removes source spreading and
% leaves the inclusion-induced shadow/lensing pattern.
ampRatio = ampIncl ./ max(ampHom, eps);
shadowLossDB = -20*log10(max(ampRatio, eps));

x = outIncl.x(:).';
z = outIncl.z(:);
[X, Z] = meshgrid(x, z);
cx = CFG.InclusionCenter(1);
cz = CFG.InclusionCenter(2);
r = CFG.InclusionRadius;
inclusionMask = (X - cx).^2 + (Z - cz).^2 <= r^2;
behindMask = X > (cx + r);
farX = cx + 2*r;
[~, ixFar] = min(abs(x - farX));
[~, izCenter] = min(abs(z - cz));

profileVerticalDB = shadowLossDB(:, ixFar);
profileLateralDB = shadowLossDB(izCenter, :);
profileAmpIncl = ampIncl(izCenter, :);
profileAmpHom = ampHom(izCenter, :);

thr = CFG.ShadowThresholdDB;
pixelAreaMM2 = (mean(diff(x))*1e3) * (mean(diff(z))*1e3);
shadowAreaBehindMM2 = nnz(behindMask & shadowLossDB > thr) * pixelAreaMM2;
verticalWidthMM = nnz(profileVerticalDB > thr) * mean(diff(z))*1e3;

R = struct();
R.key = string(C.key);
R.f0 = C.f0;
R.aperture_mm = C.aperture_mm;
R.attenuation_mode = string(C.attenuation_mode);
R.x = x;
R.z = z;
R.cs_map = outIncl.cs_map;
R.alpha_map = outIncl.alpha_map;
R.U_inclusion = outIncl.U;
R.U_homogeneous = outHom.U;
R.amp_inclusion = ampIncl;
R.amp_homogeneous = ampHom;
R.amp_ratio = ampRatio;
R.shadow_loss_db = shadowLossDB;
R.phase_inclusion = angle(outIncl.U);
R.phase_homogeneous = angle(outHom.U);
R.source_map = source_map_2d(x, z, cfg.SourceDiagnostic);
R.inclusion_mask = inclusionMask;
R.behind_mask = behindMask;
R.ix_far = ixFar;
R.iz_center = izCenter;
R.profile_vertical_db = profileVerticalDB;
R.profile_lateral_db = profileLateralDB;
R.profile_amp_inclusion = profileAmpIncl;
R.profile_amp_homogeneous = profileAmpHom;
R.max_shadow_loss_db = max(shadowLossDB(behindMask), [], 'omitnan');
R.mean_shadow_loss_behind_db = mean(shadowLossDB(behindMask), 'omitnan');
R.shadow_area_gt3db_behind_mm2 = shadowAreaBehindMM2;
R.vertical_shadow_width_gt3db_mm = verticalWidthMM;
R.mean_amp_inclusion_roi = mean(ampIncl(inclusionMask), 'omitnan');
R.mean_amp_behind = mean(ampIncl(behindMask), 'omitnan');
end

function S = summary_row(R, C)
S = struct();
S.f0 = C.f0;
S.aperture_mm = C.aperture_mm;
S.attenuation_mode = string(C.attenuation_mode);
S.max_shadow_loss_db = R.max_shadow_loss_db;
S.mean_shadow_loss_behind_db = R.mean_shadow_loss_behind_db;
S.shadow_area_gt3db_behind_mm2 = R.shadow_area_gt3db_behind_mm2;
S.vertical_shadow_width_gt3db_mm = R.vertical_shadow_width_gt3db_mm;
S.mean_amp_inclusion_roi = R.mean_amp_inclusion_roi;
S.mean_amp_behind = R.mean_amp_behind;
S.elapsed_s = R.elapsed_s;
S.output_key = string(C.key);
end

function names = summary_names()
names = {'f0','aperture_mm','attenuation_mode','max_shadow_loss_db', ...
    'mean_shadow_loss_behind_db','shadow_area_gt3db_behind_mm2', ...
    'vertical_shadow_width_gt3db_mm','mean_amp_inclusion_roi', ...
    'mean_amp_behind','elapsed_s','output_key'};
end

%% Plots

function plot_condition_pair(R, C, OUT)
fig = figure('Color','w', 'Position', [40 40 1800 1120], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 3, 4, 'TileSpacing','compact', 'Padding','compact');

plot_map(R.x, R.z, R.cs_map, 'True SWS', 'm/s', R.inclusion_mask);
plot_map(R.x, R.z, R.source_map, 'Source weights', 'relative', []);
plot_map(R.x, R.z, R.amp_inclusion, 'Inclusion amplitude |U|', 'a.u.', []);
plot_map(R.x, R.z, R.amp_homogeneous, 'Matched homogeneous |U|', 'a.u.', []);

plot_map(R.x, R.z, R.amp_ratio, 'Amplitude ratio incl/hom', 'ratio', []);
plot_map(R.x, R.z, R.shadow_loss_db, 'Inclusion shadow loss', 'dB', []);
plot_map(R.x, R.z, R.phase_inclusion, 'Inclusion phase', 'rad', []);
plot_map(R.x, R.z, R.alpha_map, 'Shear attenuation alpha(f0)', 'Np/m', []);

nexttile; hold on; grid on;
plot(1e3*R.x, R.profile_amp_homogeneous, 'LineWidth', 1.5);
plot(1e3*R.x, R.profile_amp_inclusion, 'LineWidth', 1.5);
xlabel('x (mm)'); ylabel('|U|'); title('Center lateral amplitude');
legend({'homogeneous','inclusion'}, 'Location','best');

nexttile; hold on; grid on;
plot(1e3*R.x, R.profile_lateral_db, 'LineWidth', 1.5);
xline(1e3*mean(R.x(R.inclusion_mask(R.iz_center,:))), '--k', 'inclusion');
xlabel('x (mm)'); ylabel('loss (dB)'); title('Center lateral shadow loss');

nexttile; hold on; grid on;
plot(1e3*R.z, R.profile_vertical_db, 'LineWidth', 1.5);
xlabel('z (mm)'); ylabel('loss (dB)'); title('Vertical shadow profile behind inclusion');

nexttile; axis off;
txt = sprintf(['f0: %g Hz\nsource aperture: %.1f mm\nattenuation mode: %s\n', ...
    'max behind loss: %.2f dB\nmean behind loss: %.2f dB\n', ...
    'area >3 dB: %.1f mm^2\nvertical width >3 dB: %.1f mm'], ...
    C.f0, C.aperture_mm, C.attenuation_mode, R.max_shadow_loss_db, ...
    R.mean_shadow_loss_behind_db, R.shadow_area_gt3db_behind_mm2, ...
    R.vertical_shadow_width_gt3db_mm);
text(0.02, 0.95, txt, 'VerticalAlignment','top', 'Interpreter','none');

sgtitle(sprintf('Test 66: %s', R.key), 'Interpreter','none');
file = fullfile(OUT.map_dir, sprintf('test66__%s.png', sanitize(R.key)));
exportgraphics(fig, file, 'Resolution', 170);
close(fig);
end

function plot_summary_grid(R, OUT, CFG)
if isempty(R), return; end
fig = figure('Color','w', 'Position', [40 40 1900 1180], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig, 4, numel(R), 'TileSpacing','compact', 'Padding','compact');

for i = 1:numel(R)
    plot_map(R(i).x, R(i).z, R(i).amp_inclusion, ...
        sprintf('Amp |U|\n%.1f mm, %s', R(i).aperture_mm, pretty_mode(R(i).attenuation_mode)), ...
        'a.u.', []);
end
for i = 1:numel(R)
    plot_map(R(i).x, R(i).z, R(i).shadow_loss_db, ...
        sprintf('Shadow loss\n%.1f mm, %s', R(i).aperture_mm, pretty_mode(R(i).attenuation_mode)), ...
        'dB', []);
end
for i = 1:numel(R)
    nexttile; hold on; grid on;
    plot(1e3*R(i).x, R(i).profile_lateral_db, 'LineWidth', 1.3);
    xlabel('x (mm)'); ylabel('loss (dB)');
    title(sprintf('Lateral profile\n%.1f mm, %s', R(i).aperture_mm, pretty_mode(R(i).attenuation_mode)), ...
        'Interpreter','none');
end
for i = 1:numel(R)
    nexttile; hold on; grid on;
    plot(1e3*R(i).z, R(i).profile_vertical_db, 'LineWidth', 1.3);
    xlabel('z (mm)'); ylabel('loss (dB)');
    title(sprintf('Vertical behind inclusion\n%.1f mm, %s', R(i).aperture_mm, pretty_mode(R(i).attenuation_mode)), ...
        'Interpreter','none');
end

sgtitle(sprintf('Test 66 shadow comparison, f=%g Hz', R(1).f0));
file = fullfile(OUT.figure_dir, sprintf('test66_summary_shadow_comparison_f%d.png', R(1).f0));
exportgraphics(fig, file, 'Resolution', 170);
close(fig);

T = struct2table(R);
fig2 = figure('Color','w', 'Position', [80 80 1300 680], ...
    'Toolbar','none', 'MenuBar','none');
tiledlayout(fig2, 1, 3, 'TileSpacing','compact', 'Padding','compact');
labels = strings(numel(R),1);
for i = 1:numel(R)
    labels(i) = sprintf('%.1f mm\n%s', R(i).aperture_mm, pretty_mode(R(i).attenuation_mode));
end
nexttile; bar(categorical(labels), T.max_shadow_loss_db); grid on;
ylabel('max loss behind inclusion (dB)'); title('Maximum shadow loss');
nexttile; bar(categorical(labels), T.shadow_area_gt3db_behind_mm2); grid on;
ylabel('area (mm^2)'); title('Behind-inclusion area >3 dB');
nexttile; bar(categorical(labels), T.vertical_shadow_width_gt3db_mm); grid on;
ylabel('width (mm)'); title('Vertical width >3 dB');
sgtitle('Test 66 shadow metrics');
exportgraphics(fig2, fullfile(OUT.figure_dir, 'test66_shadow_metric_bars.png'), 'Resolution', 170);
close(fig2);
end

function plot_map(x, z, A, ttl, cblabel, overlayMask)
nexttile;
imagesc(1e3*x, 1e3*z, A);
axis image; set(gca,'YDir','normal');
title(ttl, 'Interpreter','none');
xlabel('x (mm)'); ylabel('z (mm)');
cb = colorbar; ylabel(cb, cblabel, 'Interpreter','none');
if ~isempty(overlayMask)
    hold on;
    contour(1e3*x, 1e3*z, double(overlayMask), [0.5 0.5], 'w--', 'LineWidth', 1.2);
end
end

function M = source_map_2d(x, z, src)
[X, Z] = meshgrid(x, z);
M = zeros(size(X));
sigma = 0.6e-3;
for i = 1:size(src.PositionXYZ,1)
    sx = src.PositionXYZ(i,1);
    sz = src.PositionXYZ(i,3);
    M = M + src.Amplitude(i) .* exp(-((X-sx).^2 + (Z-sz).^2)/(2*sigma^2));
end
if max(M(:)) > 0
    M = M ./ max(M(:));
end
end

function p = pretty_mode(s)
s = string(s);
switch s
    case "straight_ray_integral"
        p = "straight ray";
    case "eikonal_ray_integral"
        p = "eikonal ray";
    otherwise
        p = s;
end
end

%% Utilities

function row = struct2row(S, names)
row = cell(1, numel(names));
for i = 1:numel(names)
    row{i} = S.(names{i});
end
end

function write_config_json(CFG, file)
txt = jsonencode(CFG, PrettyPrint=true);
fid = fopen(file, 'w');
if fid < 0, error('Could not write %s', file); end
fwrite(fid, txt, 'char');
fclose(fid);
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_=-]+', '_');
end

function v = getenv_default(name, defaultVal)
v = getenv(name);
if isempty(v)
    v = defaultVal;
end
end

function v = env_string(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    v = string(defaultVal);
else
    v = string(raw);
end
end

function tf = env_true(name, defaultVal)
raw = getenv(name);
if isempty(raw)
    tf = logical(defaultVal);
else
    tf = any(strcmpi(raw, {'1','true','yes','on'}));
end
end
