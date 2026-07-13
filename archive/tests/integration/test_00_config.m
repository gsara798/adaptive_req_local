%% test_00_config.m
% Basic configuration test for the adaptive_req project.
%
% This script verifies that:
%   1. The project setup file correctly adds src/ to the MATLAB path.
%   2. The default simulation configuration can be created.
%   3. The default feature configuration can be created.
%   4. The derived REQ configuration is consistent with cfg and feat_cfg.
%   5. The circular Hann window used by REQ can be visualized.

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
%   this_dir  -> adaptive_req/tests/integration
%   root_dir  -> adaptive_req
%
% The setup file adds:
%
%   adaptive_req/src/
%
% to the MATLAB path, which makes the package +adaptive_req available.

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
root_dir = fileparts(fileparts(this_dir));

run(fullfile(root_dir, 'setup_adaptive_req.m'));

%% ========================================================================
% 2. Create default simulation configuration
% ========================================================================
% cfg contains the physical and numerical simulation settings, such as:
%   f0, cs_bg, Lx, Lz, dx, dz, Nwaves, SNR, and angular sampling settings.

cfg = adaptive_req.config.default_sim_config();

%% ========================================================================
% 3. Create default feature configuration
% ========================================================================
% feat_cfg contains settings used for local spectral feature extraction.
%
% Important fields include:
%   M          : local window size in wavelengths
%   cs_guess   : guessed shear wave speed used to define lambda_guess
%   gamma_win  : circular Hann shrink factor
%   pad_factor : zero-padding factor for spectral analysis

feat_cfg = adaptive_req.config.default_feature_config();

%% ========================================================================
% 4. Create derived REQ configuration
% ========================================================================
% default_req_config combines cfg and feat_cfg to compute derived REQ
% quantities, such as:
%
%   win_size
%   half_win
%   W2
%   PAD
%   k0_true
%   lambda_guess_used
%   Nbins settings
%
% The second output updates feat_cfg with derived quantities such as
% win_size, half_win, lambda_guess_used, and M_eff.

[req_cfg, feat_cfg] = adaptive_req.config.default_req_config(cfg, feat_cfg);

%% ========================================================================
% 5. Display configurations
% ========================================================================

disp(' ');
disp('================ cfg ================');
disp(cfg);

disp(' ');
disp('============= feat_cfg ==============');
disp(feat_cfg);

disp(' ');
disp('============= req_cfg ===============');
disp(req_cfg);

%% ========================================================================
% 6. Print key derived quantities
% ========================================================================

fprintf('\nConfiguration summary.\n');
fprintf('f0 = %.1f Hz\n', cfg.f0);
fprintf('cs_bg = %.3f m/s\n', cfg.cs_bg);
fprintf('k0_true = %.3f rad/m\n', req_cfg.k0_true);

fprintf('\nWindow summary.\n');
fprintf('M = %.2f\n', feat_cfg.M);
fprintf('cs_guess_used = %.3f m/s\n', feat_cfg.cs_guess_used);
fprintf('lambda_guess = %.3f cm\n', feat_cfg.lambda_guess_used * 100);
fprintf('win_size = %d pixels\n', feat_cfg.win_size);
fprintf('patch size = %.3f cm\n', feat_cfg.win_size * cfg.dx * 100);
fprintf('M_eff = %.3f\n', feat_cfg.M_eff);

fprintf('\nREQ spectral settings.\n');

if isfield(req_cfg, 'Nbins')
    fprintf('Nbins requested = %s\n', string(req_cfg.Nbins));
end

if isfield(req_cfg, 'Nbins_auto_oversample')
    fprintf('Nbins auto oversample = %.2f\n', req_cfg.Nbins_auto_oversample);
end

if isfield(req_cfg, 'Nbins_min')
    fprintf('Nbins min = %d\n', req_cfg.Nbins_min);
end

if isfield(req_cfg, 'smooth_sigma')
    fprintf('smooth_sigma = %.3f\n', req_cfg.smooth_sigma);
end

%% ========================================================================
% 7. Visualize REQ circular Hann window
% ========================================================================
% W2 is the spatial window applied to each local patch before computing the
% local 2D Fourier spectrum.

figure('Color', 'w', 'Position', [200 200 430 360]);

imagesc(req_cfg.W2);
axis image;
colormap parula;
colorbar;

xlabel('x pixel');
ylabel('z pixel');
title('REQ circular Hann window');

set(gca, 'YDir', 'normal');
grid on;
box on;

fprintf('\nTest 00 completed successfully.\n');