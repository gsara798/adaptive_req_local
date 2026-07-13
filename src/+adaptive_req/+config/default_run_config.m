function CFG = default_run_config()
%DEFAULT_RUN_CONFIG Default high-level configuration for adaptive_req runs.
%
% This function defines the user-facing configuration.
% It does not call run_aperture_sweep and does not generate figures.

CFG = struct();

%% Experiment settings

CFG.EXP = struct();

CFG.EXP.name = 'default_adaptive_req_run';
CFG.EXP.sampling_mode = 'cone';

CFG.EXP.num_steps = 10;
CFG.EXP.step_indices = [];

CFG.EXP.num_realizations = 1;
CFG.EXP.num_patches = 3;

CFG.EXP.seed_base = 1000;

CFG.EXP.selected_step = 1;
CFG.EXP.selected_patch = 1;

%% Simulation overrides

CFG.SIM = struct();

CFG.SIM.f0 = 600;
CFG.SIM.cs_bg = 3.0;

CFG.SIM.Nwaves = 2000;
CFG.SIM.SNR = Inf;
CFG.SIM.AmpJitter = 0;

CFG.SIM.Lx = 0.06;
CFG.SIM.Lz = 0.06;
CFG.SIM.dx = 2.5e-4;
CFG.SIM.dz = 2.5e-4;

CFG.SIM.WaveModel = 'planewave';
CFG.SIM.SourceSampling = 'cone';
CFG.SIM.AngularSamplingMethod = 'fibonacci';

%% REQ settings

CFG.REQ = struct();

CFG.REQ.M = 3;
CFG.REQ.cs_guess = 3.0;
CFG.REQ.gamma_win = 1.0;
CFG.REQ.pad_factor = 2.0;

CFG.REQ.Nbins = 'auto';
CFG.REQ.Nbins_auto_oversample = 1;
CFG.REQ.Nbins_min = 16;

CFG.REQ.smooth_sigma = 0.001;

CFG.REQ.use_donut = false;
CFG.REQ.donut_cs_min = 1.0;
CFG.REQ.donut_cs_max = 5.0;
CFG.REQ.donut_taper_rel = 0.06;
CFG.REQ.apply_donut_to_final_map = false;

%% Output settings

CFG.OUTPUT = struct();

CFG.OUTPUT.store_req_curve = false;
CFG.OUTPUT.store_req_mapping = true;
CFG.OUTPUT.compute_global_req = false;
CFG.OUTPUT.store_global_req_mapping = true;
CFG.OUTPUT.store_req_metadata = true;
CFG.OUTPUT.store_feature_struct = false;
CFG.OUTPUT.verbose = true;

CFG.OUTPUT.save_condition_table = true;
CFG.OUTPUT.save_condition_mat = false;
CFG.OUTPUT.save_condition_summary_figures = false;

%% Plot settings

CFG.PLOT = struct();

CFG.PLOT.show_step_diagnostics = false;
CFG.PLOT.save_step_diagnostics = false;
CFG.PLOT.step_diagnostic_visible = false;
CFG.PLOT.close_step_diagnostics_after_save = true;

CFG.PLOT.show_q_vs_aperture = true;
CFG.PLOT.show_feature_space = true;
CFG.PLOT.show_feature_vs_q = true;
CFG.PLOT.show_feature_grid = true;

CFG.PLOT.show_selected_step_wavefield = false;
CFG.PLOT.show_all_step_wavefields = false;

CFG.PLOT.save_summary_figures = false;

%% Save settings

CFG.SAVE = struct();

CFG.SAVE.save_png = true;
CFG.SAVE.save_pdf = true;
CFG.SAVE.save_fig = false;

CFG.SAVE.png_resolution = 300;
CFG.SAVE.close_after_save = false;

end
