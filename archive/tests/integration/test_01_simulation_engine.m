%% test_01_simulation_engine.m
% Test one synthetic simulation from the adaptive_req package.
%
% This script verifies that:
%   1. The project setup works.
%   2. The default simulation, feature, and REQ configurations are valid.
%   3. One aperture condition can be selected from an aperture schedule.
%   4. A synthetic wavefield can be generated.
%   5. Local analysis patches can be placed on the wavefield.
%   6. The patch size is consistent with M * lambda_guess.
%
% This test runs only one simulation. It does not perform a full aperture
% sweep. The full sweep is tested in test_02_aperture_sweep.m.

clear; clc; close all;
format compact;

%% ========================================================================
% 1. Locate project root and set up package path
% ========================================================================
% This script is assumed to live in:
%
%   adaptive_req/tests/integration/
%
% Therefore:
%   this_dir -> adaptive_req/tests/integration
%   root_dir -> adaptive_req
%
% The setup file adds:
%
%   adaptive_req/src/
%
% to the MATLAB path.

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

run(fullfile(root_dir, 'setup_adaptive_req.m'));

%% ========================================================================
% 2. Test settings
% ========================================================================

EXP = struct();

EXP.name = 'test_01_simulation_engine';

% Aperture schedule settings.
% num_steps defines how many aperture conditions are available.
% step_idx selects only one of those conditions for this test.
EXP.sampling_mode = 'cone';
EXP.num_steps = 5;
EXP.step_idx = 2;

% Local patch settings.
% NumPatches controls how many local analysis windows are drawn.
% M controls the size of each window through M * lambda_guess.
EXP.num_patches = 4;
EXP.M = 2;

% Reproducibility.
EXP.seed_base = 1000;

%% ========================================================================
% 3. Create simulation, feature, and REQ configurations
% ========================================================================

cfg = adaptive_req.config.default_sim_config();

feat_cfg = adaptive_req.config.default_feature_config( ...
    'M', EXP.M);

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config( ...
    cfg, feat_cfg); %#ok<ASGLU>

%% ========================================================================
% 4. Build aperture schedule and select one step
% ========================================================================
% For sampling_mode = 'cone', the schedule defines a sequence of cone
% apertures and their corresponding solid angles.
%
% In this test, we only select one aperture step. The full aperture sweep is
% handled later by adaptive_req.studies.run_aperture_sweep.

schedule = adaptive_req.simulate.build_aperture_schedule( ...
    EXP.sampling_mode, EXP.num_steps);

if EXP.step_idx < 1 || EXP.step_idx > EXP.num_steps
    error('EXP.step_idx must be between 1 and EXP.num_steps.');
end

cfg_i = adaptive_req.simulate.apply_sampling_step( ...
    cfg, EXP.sampling_mode, schedule.values(EXP.step_idx));

cfg_i.Seed = EXP.seed_base + EXP.step_idx;

%% ========================================================================
% 5. Run one synthetic simulation
% ========================================================================
% sim contains:
%   sim.Uxz     complex wavefield with convention Uxz(z,x)
%   sim.x       x coordinates in meters
%   sim.z       z coordinates in meters
%   sim.cs_map  shear wave speed map with convention cs_map(z,x)
%   sim.k_map   wavenumber map with convention k_map(z,x)
%   sim.diag    diagnostic information from the simulation engine

sim = adaptive_req.simulate.run_single_simulation(cfg_i);

%% ========================================================================
% 6. Build local patch windows
% ========================================================================
% The patches define the local regions where spectra, radial energy curves,
% features, and q_theory will later be computed.
%
% Important:
%   NumPatches changes the number and placement of local windows.
%   M changes the physical size of each local window.

patch_pack = adaptive_req.simulate.build_patch_windows( ...
    cfg_i, feat_cfg, ...
    'NumPatches', EXP.num_patches);

%% ========================================================================
% 7. Print simulation diagnostics
% ========================================================================

fprintf('\nSimulation completed.\n');
fprintf('Uxz size = %d z-pixels by %d x-pixels\n', ...
    size(sim.Uxz, 1), size(sim.Uxz, 2));

fprintf('\nPhysical parameters.\n');
fprintf('f0 = %.1f Hz\n', cfg_i.f0);
fprintf('cs_bg = %.3f m/s\n', cfg_i.cs_bg);
fprintf('k0 = %.3f rad/m\n', 2*pi*cfg_i.f0/cfg_i.cs_bg);

fprintf('\nAperture condition.\n');
fprintf('sampling_mode = %s\n', EXP.sampling_mode);
fprintf('step_idx = %d / %d\n', EXP.step_idx, EXP.num_steps);
fprintf('aperture value = %.3f %s\n', ...
    schedule.values(EXP.step_idx), schedule.unit);
fprintf('Omega = %.3f sr\n', schedule.Omega_sr(EXP.step_idx));

if isfield(cfg_i, 'ConeHalfAngleDeg')
    fprintf('ConeHalfAngleDeg = %.2f deg\n', cfg_i.ConeHalfAngleDeg);
end

fprintf('\nPatch diagnostics.\n');
fprintf('Number of patches = %d\n', patch_pack.n_patches);
fprintf('M = %.2f\n', feat_cfg.M);
fprintf('cs_guess_used = %.3f m/s\n', feat_cfg.cs_guess_used);
fprintf('lambda_guess = %.3f cm\n', feat_cfg.lambda_guess_used * 100);
fprintf('M * lambda_guess = %.3f cm\n', ...
    feat_cfg.M * feat_cfg.lambda_guess_used * 100);
fprintf('win_size = %d pixels\n', patch_pack.win_size);
fprintf('patch size in x = %.3f cm\n', patch_pack.win_size * cfg_i.dx * 100);
fprintf('patch size in z = %.3f cm\n', patch_pack.win_size * cfg_i.dz * 100);

%% ========================================================================
% 8. Visualize simulated wavefield and patch locations
% ========================================================================
% The plotted wavefield is the real part of the complex narrowband field.
% The white rectangles show local analysis patches.

fig = figure('Color', 'w', 'Position', [100 100 620 500]);
ax = axes('Parent', fig);

imagesc(ax, sim.x * 100, sim.z * 100, real(sim.Uxz));
set(ax, 'YDir', 'normal');

axis(ax, 'image');
colormap(ax, parula);

cb = colorbar(ax);
cb.Label.String = 'Re\{u(x,z)\}';
cb.Label.Interpreter = 'tex';

xlabel(ax, 'x (cm)');
ylabel(ax, 'z (cm)');

title(ax, sprintf('Single simulated wavefield, step %d', EXP.step_idx), ...
    'Interpreter', 'tex');

hold(ax, 'on');

for pidx = 1:patch_pack.n_patches

    xi = patch_pack.x_idx_list{pidx};
    zi = patch_pack.z_idx_list{pidx};

    x0 = sim.x(min(xi)) * 100;
    z0 = sim.z(min(zi)) * 100;

    w = (sim.x(max(xi)) - sim.x(min(xi))) * 100;
    h = (sim.z(max(zi)) - sim.z(min(zi))) * 100;

    rectangle(ax, ...
        'Position', [x0 z0 w h], ...
        'EdgeColor', 'w', ...
        'LineWidth', 1.2);

    xc = sim.x(patch_pack.cx_list(pidx)) * 100;
    zc = sim.z(patch_pack.cz_list(pidx)) * 100;

    text(ax, xc, zc, sprintf('%d', pidx), ...
        'Color', 'w', ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

hold(ax, 'off');

grid(ax, 'on');
box(ax, 'on');

try
    adaptive_req.figures.apply_style(fig, ax);
catch
    set(ax, 'FontSize', 14, 'LineWidth', 1.2);
end

drawnow;

fprintf('\nTest 01 completed successfully.\n');