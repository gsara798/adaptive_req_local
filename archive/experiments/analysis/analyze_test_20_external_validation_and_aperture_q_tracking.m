%% analyze_test_20_external_validation_and_aperture_q_tracking.m
% Test 20, Level 01: external validation and homogeneous aperture q tracking.
%
% Part A applies the Test 19 and historical Test 12 deployments to the
% external homogeneous/inclusion/bilayer cases introduced in Test 17.
% Part B follows local q through a ten-step homogeneous cone-aperture sweep.
% No model is trained in this script.

clear; clc; close all;
format compact;

set(groot,'defaultAxesFontSize',11);
set(groot,'defaultTextFontSize',11);
set(groot,'defaultLegendFontSize',9);

%% Project setup and fixed Test 20 design

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
% setup_style uses large presentation defaults; Test 20 has dense panels.
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');
set(groot,'defaultAxesFontSize',8);
set(groot,'defaultTextFontSize',8);
set(groot,'defaultLegendFontSize',7);

CFG = adaptive_req.config.load_profile_config( ...
    'test_17_synthetic_inclusion_kWave_like','RootDir',root_dir);
CFG.REQ.M_list = [2 3 4];
CFG.REQ.cs_guess = 3;
CFG.REQ.StepX = 1;
CFG.REQ.StepZ = 1;
CFG.REQ.EdgeMode = 'valid';
CFG.FIG.Resolution = 220;

AP = struct();
AP.cs_values = [2 3 4];
AP.f0 = 500;
AP.dx = 0.5e-3;
AP.dz = 0.5e-3;
AP.Lx = 0.05;
AP.Lz = 0.05;
AP.Nwaves = 2000;
AP.num_steps = 10;
AP.REQ_M = [2 3 4];
AP.cs_guess = 3;
AP.seed_base = 20001;
% [x,z] sample indices: center, left, right, top, bottom.
AP.patch_centers = [51 51; 26 51; 76 51; 51 26; 51 76];
AP.patch_names = ["center";"left";"right";"top";"bottom"];

RUN_MODE=lower(string(getenv('ADAPTIVE_REQ_TEST20_MODE')));
if RUN_MODE=="", RUN_MODE="all"; end
assert(ismember(RUN_MODE,["all","part_a_only","part_b_only"]), ...
    'ADAPTIVE_REQ_TEST20_MODE must be all, part_a_only, or part_b_only.');

OUT = make_output_dirs(root_dir);
GEOMETRIES = build_geometry_specs(CFG);
MODEL_SPECS = load_all_models(root_dir);
validate_model_predictors(MODEL_SPECS);

fprintf('\nTest 20, Level 01\n');
fprintf('Models loaded: %d learned + TheoryQDiscrete.\n',numel(MODEL_SPECS));
disp(struct2table(rmfield(MODEL_SPECS,'model')));

%% Part A: external heterogeneous validation

if RUN_MODE=="part_b_only"
    external_metrics_file=fullfile(OUT.table_dir, ...
        'level20_external_sws_metrics.csv');
    assert(exist(external_metrics_file,'file')==2, ...
        'Part A metrics not found: %s',external_metrics_file);
    T_external_sws=readtable(external_metrics_file, ...
        'TextType','string','VariableNamingRule','preserve');
    fprintf('\nSkipping completed Part A and loading metrics from:\n%s\n', ...
        external_metrics_file);
else
external_checkpoint = fullfile(OUT.data_dir,'level20_external_checkpoint.mat');
T_external = table();
EXTERNAL_CASES = struct([]);
completed_external = strings(0,1);
if exist(external_checkpoint,'file') == 2
    S = load(external_checkpoint,'T_external','EXTERNAL_CASES', ...
        'completed_external');
    T_external = S.T_external;
    EXTERNAL_CASES = S.EXTERNAL_CASES;
    completed_external = S.completed_external;
    fprintf('\nResuming Part A: %d of %d conditions complete.\n', ...
        numel(completed_external),numel(GEOMETRIES)*numel(CFG.CASES)*numel(CFG.REQ.M_list));
end

condition_id = numel(completed_external);
for gi = 1:numel(GEOMETRIES)
    geometry = GEOMETRIES(gi);
    for wi = 1:numel(CFG.CASES)
        wave_case = CFG.CASES(wi);
        regime = canonical_regime_label(wave_case);
        for mi = 1:numel(CFG.REQ.M_list)
            M = CFG.REQ.M_list(mi);
            key = geometry.geometry_id + "__" + regime + "__M" + string(M);
            if any(completed_external == key)
                continue;
            end

            condition_id = condition_id + 1;
            fprintf('\nPart A %d: %s | %s | M=%g\n',condition_id, ...
                geometry.geometry_id,regime,M);
            [T_feat,sim,req_out,cfg_sim] = load_or_generate_test17_condition( ...
                root_dir,CFG,geometry,wave_case,M,gi,wi,condition_id);
            T_feat = prepare_external_features(T_feat,cfg_sim,geometry, ...
                wave_case,regime,M,condition_id);

            T_condition = predict_all_models( ...
                MODEL_SPECS,T_feat,cfg_sim,regime,true);
            T_condition = add_prediction_errors(T_condition);
            T_condition = add_q_gradients(T_condition,cfg_sim.dx,cfg_sim.dz);
            T_external = concat_tables(T_external,T_condition);

            k = numel(EXTERNAL_CASES)+1;
            EXTERNAL_CASES(k).key = key;
            EXTERNAL_CASES(k).geometry = geometry;
            EXTERNAL_CASES(k).regime = regime;
            EXTERNAL_CASES(k).REQ_M = M;
            EXTERNAL_CASES(k).sim = sim;
            EXTERNAL_CASES(k).req_out = req_out;
            EXTERNAL_CASES(k).T_pred = T_condition;

            completed_external(end+1,1) = key; %#ok<SAGROW>
            save(external_checkpoint,'CFG','GEOMETRIES','MODEL_SPECS', ...
                'T_external','EXTERNAL_CASES','completed_external','-v7.3');
        end
    end
end

assert(~isempty(T_external),'Part A produced no external predictions.');
external_groups = ["model_name","feature_set","model_role","model_type", ...
    "geometry_id","field_regime_label","REQ_M"];
T_external_sws = summarize_sws(T_external,external_groups);
T_external_q = summarize_q(T_external,external_groups);
T_external_distance = summarize_binned(T_external( ...
    T_external.geometry_type=="bilayer" | T_external.geometry_type=="inclusion",:), ...
    'distance_to_interface_m','distance_bin_mm');
T_external_qgrad = summarize_binned(T_external( ...
    T_external.model_role~="diagnostic_only",:),'q_gradient_mag','q_gradient_bin');

if RUN_MODE~="part_a_only"
    writetable(remove_cell_columns(T_external),fullfile(OUT.table_dir, ...
        'level20_external_predictions.csv'));
    writetable(T_external_sws,fullfile(OUT.table_dir, ...
        'level20_external_sws_metrics.csv'));
    writetable(T_external_q,fullfile(OUT.table_dir, ...
        'level20_external_q_metrics.csv'));
    writetable(T_external_distance,fullfile(OUT.table_dir, ...
        'level20_external_error_vs_interface_distance.csv'));
    writetable(T_external_qgrad,fullfile(OUT.table_dir, ...
        'level20_external_error_vs_q_gradient.csv'));

    save(fullfile(OUT.data_dir,'level20_external_validation.mat'), ...
        'CFG','GEOMETRIES','MODEL_SPECS','T_external','T_external_sws', ...
        'T_external_q','T_external_distance','T_external_qgrad', ...
        'EXTERNAL_CASES','-v7.3');
end

plot_external_metric(T_external_sws,'MAPE_pct','SWS MAPE (%)', ...
    'level20_external_mape_by_model_case_regime_M.png',OUT);
plot_external_metric(T_external_sws,'HighError_gt20_pct', ...
    'High-error >20% (%)', ...
    'level20_external_high_error_by_model_case_regime_M.png',OUT);
plot_external_maps(EXTERNAL_CASES,OUT);
plot_binned_error(T_external_distance,'distance_bin_mm', ...
    'Distance to interface bin (mm)','External error vs interface distance', ...
    'level20_external_error_vs_interface_distance.png',OUT);
plot_binned_error(T_external_qgrad,'q_gradient_bin', ...
    '|grad q| bin','External error vs |grad q|', ...
    'level20_external_error_vs_q_gradient.png',OUT);

if RUN_MODE=="part_a_only"
    fprintf('\nTest 20 Part A figures regenerated from checkpoints.\n');
    fprintf('Figure folder:\n%s\n',OUT.fig_dir);
    return;
end
end

%% Part B: homogeneous aperture sweep q tracking

aperture_checkpoint = fullfile(OUT.data_dir,'level20_aperture_checkpoint.mat');
T_aperture = table();
APERTURE_SWEEPS = struct([]);
completed_aperture = strings(0,1);
if exist(aperture_checkpoint,'file') == 2
    S = load(aperture_checkpoint,'T_aperture','APERTURE_SWEEPS', ...
        'completed_aperture');
    T_aperture = S.T_aperture;
    APERTURE_SWEEPS = S.APERTURE_SWEEPS;
    completed_aperture = S.completed_aperture;
    fprintf('\nResuming Part B: %d of %d sweeps complete.\n', ...
        numel(completed_aperture),numel(AP.cs_values)*numel(AP.REQ_M));
end

for ci = 1:numel(AP.cs_values)
    cs = AP.cs_values(ci);
    for mi = 1:numel(AP.REQ_M)
        M = AP.REQ_M(mi);
        key = "cs"+string(cs)+"__M"+string(M);
        if any(completed_aperture==key)
            continue;
        end

        fprintf('\nPart B: homogeneous cs=%g m/s | M=%g\n',cs,M);
        [T_feat,sweep] = run_homogeneous_aperture_sweep(AP,cs,M,ci,mi);
        T_feat = prepare_aperture_features(T_feat,AP,cs,M);
        % Unknown is deliberate: aperture must not leak through the user guess.
        T_condition = predict_all_models( ...
            MODEL_SPECS,T_feat,aperture_cfg(AP,cs,ci,mi),"unknown",false);
        T_condition = add_prediction_errors(T_condition);
        T_aperture = concat_tables(T_aperture,T_condition);

        k = numel(APERTURE_SWEEPS)+1;
        APERTURE_SWEEPS(k).key = key;
        APERTURE_SWEEPS(k).cs_true = cs;
        APERTURE_SWEEPS(k).REQ_M = M;
        APERTURE_SWEEPS(k).sweep = sweep;
        APERTURE_SWEEPS(k).T_pred = T_condition;
        completed_aperture(end+1,1) = key; %#ok<SAGROW>
        save(aperture_checkpoint,'AP','MODEL_SPECS','T_aperture', ...
            'APERTURE_SWEEPS','completed_aperture','-v7.3');
    end
end

assert(~isempty(T_aperture),'Part B produced no aperture predictions.');
aperture_groups = ["model_name","feature_set","model_role","model_type", ...
    "SIM_cs_bg","REQ_M"];
T_aperture_q = summarize_q(T_aperture,aperture_groups);
T_aperture_sws = summarize_sws(T_aperture,aperture_groups);

writetable(remove_cell_columns(T_aperture),fullfile(OUT.table_dir, ...
    'level20_aperture_sweep_predictions.csv'));
writetable(T_aperture_q,fullfile(OUT.table_dir, ...
    'level20_aperture_sweep_q_metrics.csv'));
writetable(T_aperture_sws,fullfile(OUT.table_dir, ...
    'level20_aperture_sweep_sws_metrics.csv'));
save(fullfile(OUT.data_dir,'level20_aperture_q_tracking.mat'), ...
    'AP','MODEL_SPECS','T_aperture','T_aperture_q', ...
    'T_aperture_sws','APERTURE_SWEEPS','-v7.3');

plot_q_vs_aperture(T_aperture,OUT);
plot_aperture_error(T_aperture,'q_error','q error', ...
    'q_error_vs_aperture.png',OUT);
plot_aperture_error(T_aperture,'sws_error_pct','SWS signed error (%)', ...
    'sws_error_vs_aperture.png',OUT);
plot_aperture_q_scatter(T_aperture,OUT);
plot_aperture_patch_map(AP,OUT);

%% Final console summary

fprintf('\nTest 20, Level 01 complete.\nAnalysis folder:\n%s\n',OUT.analysis_dir);
print_external_summary(T_external_sws);
print_old_new_summary(T_external_sws);
print_aperture_trend_summary(T_aperture);
print_user_guess_summary(T_external_sws,T_aperture_sws);

%% Local functions

function OUT = make_output_dirs(root_dir)
OUT.analysis_dir = fullfile(root_dir,'outputs', ...
    'test_20_external_validation_and_aperture_q_tracking','analysis');
OUT.fig_dir = fullfile(OUT.analysis_dir,'figures');
OUT.map_dir = fullfile(OUT.fig_dir,'external_maps');
OUT.table_dir = fullfile(OUT.analysis_dir,'tables');
OUT.data_dir = fullfile(OUT.analysis_dir,'data');
dirs = string(struct2cell(OUT));
for i=1:numel(dirs)
    if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end
end
end

function specs = build_geometry_specs(CFG)
specs = [ ...
    make_geometry("homogeneous_cs2","Homogeneous c_s=2 m/s","homogeneous",2,2,CFG)
    make_geometry("homogeneous_cs3","Homogeneous c_s=3 m/s","homogeneous",3,3,CFG)
    make_geometry("bilayer_2_3","Bilayer 2/3 m/s","bilayer",2,3,CFG)
    make_geometry("inclusion_2_3","Circular inclusion 2/3 m/s","inclusion",2,3,CFG)];
end

function g = make_geometry(id,name,type,low,high,CFG)
g = struct('geometry_id',string(id),'geometry_name',string(name), ...
    'geometry_type',string(type),'cs_low',low,'cs_high',high, ...
    'center',CFG.INCLUSION.Center,'radius',CFG.INCLUSION.Radius, ...
    'interface_half_width',CFG.INCLUSION.InterfaceHalfWidth, ...
    'far_margin',CFG.INCLUSION.BackgroundFarMargin);
if string(type)=="bilayer", g.center=[0.025 0.025]; g.radius=NaN; end
end

function SPECS = load_all_models(root_dir)
new_dir = fullfile(root_dir,'outputs','model_registry','test19_clean_field_regime');
old_dir = fullfile(root_dir,'outputs','model_registry','test12_hybrid_baseline');
new_req = [ ...
    request("LocalOnly_T18","CleanFieldRegime_noUser","LocalOnly_T18","operational")
    request("GlobalOnly_T18","CleanFieldRegime_noUser","GlobalOnly_T18","operational")
    request("HybridLocalGlobal_T18_noUserRegime","CleanFieldRegime_noUser", ...
        "HybridLocalGlobal_T18_noUserRegime","operational")
    request("HybridLocalGlobal_T18_withUserRegimeGuess", ...
        "CleanFieldRegime_withUserGuess", ...
        "HybridLocalGlobal_T18_withUserRegimeGuess","user_informed")];
old_req = [ ...
    request("LocalOnly","NoCsGuess","LocalOnly_old","external_old_model")
    request("GlobalOnly","NoCsGuess","GlobalOnly_old","external_old_model")
    request("HybridLocalGlobal","WithCsGuess","HybridLocalGlobal_old", ...
        "external_old_model")];
n_models = numel(new_req)+numel(old_req);
SPECS = repmat(load_one(new_dir,new_req(1)),n_models,1);
for i=2:numel(new_req)
    SPECS(i)=load_one(new_dir,new_req(i));
end
offset=numel(new_req);
for i=1:numel(old_req)
    SPECS(offset+i)=load_one(old_dir,old_req(i));
end
end

function r=request(stored,feature,report,role)
r=struct('stored',string(stored),'feature',string(feature), ...
    'report',string(report),'role',string(role));
end

function s=load_one(folder,r)
[M,I,F]=adaptive_req.analysis.load_q_model_deployment(folder, ...
    'ModelName',r.stored,'FeatureSet',r.feature,'ModelType','bagged_trees');
s=struct('model_name',r.report,'feature_set',string(I.feature_set), ...
    'model_type',"bagged_trees",'model_role',r.role, ...
    'model_file',string(F),'model',M);
end

function validate_model_predictors(specs)
forbidden=lower(["q_true";"q_theory";"cs_true";"M_eff_true_diag"; ...
    "aperture_weight";"true_aperture_weight";"solid_angle_weight"]);
for i=1:numel(specs)
    p=lower(string({specs(i).model.encoder.entries.name}));
    bad=intersect(p,forbidden);
    assert(isempty(bad),'Forbidden predictors in %s: %s', ...
        specs(i).model_name,strjoin(bad,', '));
    assert(~any(contains(p,"error")|contains(p,"residual")), ...
        'Error/residual predictor found in %s.',specs(i).model_name);
end
end

function label=canonical_regime_label(wave_case)
switch string(wave_case.wave_label)
    case "directional_2d", label="directional_2D";
    case "diffuse_2d", label="diffuse_2D";
    case "partial_diffuse_3d", label="partial_3D";
    case "diffuse_3d", label="diffuse_3D";
    otherwise, error('Unknown Test 17 wave case: %s',wave_case.wave_label);
end
end

function [T,sim,req_out,cfg] = load_or_generate_test17_condition( ...
    root_dir,CFG,g,wave_case,M,gi,wi,condition_id)
file=fullfile(root_dir,'outputs','test_17_model_comparison_heterogeneous_cases', ...
    'analysis','data','conditions',sprintf('level17_compare_%s_%s_M%g.mat', ...
    g.geometry_id,sanitize_filename(wave_case.wave_label),M));
if exist(file,'file')==2
    S=load(file,'T_feat','sim','req_out','cfg_sim');
    T=S.T_feat; sim=S.sim; req_out=S.req_out; cfg=S.cfg_sim;
    fprintf('Reusing Test 17 condition: %s\n',file);
    return;
end

cfg=build_external_sim_cfg(CFG,g,wave_case,gi,wi);
sim=adaptive_req.simulate.run_single_simulation(cfg);
[feat_cfg,req_options]=build_req_settings(CFG,M);
[T,req_out]=extract_external_feature_table( ...
    sim,cfg,feat_cfg,req_options,g,wave_case,CFG,condition_id);
end

function cfg=build_external_sim_cfg(CFG,g,w,gi,wi)
cfg=struct('Lx',CFG.SIM.Lx,'Lz',CFG.SIM.Lz,'dx',CFG.SIM.dx, ...
    'dz',CFG.SIM.dz,'f0',500,'cs_bg',g.cs_low,'cs_inc',g.cs_high, ...
    'Nwaves',w.Nwaves,'Is2D',logical(w.Is2D), ...
    'WaveModel',char(w.WaveModel),'AngularSamplingMethod', ...
    char(w.AngularSamplingMethod),'ForceInPlaneWave',logical(w.ForceInPlaneWave), ...
    'SNR',Inf,'AmpJitter',0,'DecayAlpha',0, ...
    'Seed',CFG.SIM.Seed+1000*gi+101*wi,'UseParfor',CFG.SIM.UseParfor, ...
    'PhiRange',[0 2*pi],'ThetaRange',[0 pi], ...
    'AngleRange2D',[0 2*pi],'SourceSampling','ranges');
switch g.geometry_type
    case "homogeneous", masks={};
    case "inclusion"
        masks={struct('Type','circle','cs_inc',g.cs_high,'Params', ...
            struct('Center',g.center,'Radius',g.radius,'SigmaEdge',1e-6))};
    case "bilayer"
        masks={struct('Type','bilayer','cs_inc',g.cs_high,'Params', ...
            struct('Bi_Angle',0,'Bi_Offset',g.center(1),'SigmaEdge',1e-6))};
end
cfg.MaskConfig=struct('cs_bg',g.cs_low,'CombineMode','overlay','Masks',{masks});
end

function [feat_cfg,req_options]=build_req_settings(CFG,M)
feat_cfg=adaptive_req.config.default_feature_config('M',M, ...
    'cs_guess',CFG.REQ.cs_guess,'gamma_win',CFG.REQ.Gamma, ...
    'pad_factor',CFG.REQ.PadFactor);
req_options={'Nbins',CFG.REQ.Nbins,'Nbins_auto_oversample', ...
    CFG.REQ.Nbins_auto_oversample,'Nbins_min',CFG.REQ.Nbins_min, ...
    'smooth_sigma',CFG.REQ.SmoothSigma};
end

function [T,OUT]=extract_external_feature_table( ...
    sim,cfg,feat_cfg,req_options,g,w,CFG,condition_id)
[req_cfg,feat_cfg]=adaptive_req.config.default_req_config(cfg,feat_cfg,req_options{:});
cfg_req=cfg; cfg_req.UseParfor=false;
[qg,curve,fg]=adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz,cfg_req,req_cfg,feat_cfg);
sg=adaptive_req.quantile.extract_ecum_shape_features(curve);
OUT=adaptive_req.estimators.req_estimator_map(sim.Uxz,cfg_req,feat_cfg, ...
    'StepX',CFG.REQ.StepX,'StepZ',CFG.REQ.StepZ,'EdgeMode',CFG.REQ.EdgeMode, ...
    'QuantileMode','local_req','ReqOptions',req_options, ...
    'ReturnFeatures',true,'ReturnFeatureTable',true, ...
    'ReuseReqSpectrumForFeatures',true,'UseWindowParfor',true, ...
    'StoreReqCurves',false,'Verbose',false);
T=OUT.feature_table;
n=height(T);
T.global_req_mapping=repmat({adaptive_req.quantile.make_req_mapping(curve)},n,1);
T.q_global_req=qg*ones(n,1);
T.global_REQ_Nbins_effective=curve.Nbins_effective*ones(n,1);
T=assign_prefixed_fields(T,fg,"global_");
T=assign_prefixed_fields(T,sg,"global_");
T.condition_id=condition_id*ones(n,1);
T.geometry_id=repmat(g.geometry_id,n,1);
T.geometry_name=repmat(g.geometry_name,n,1);
T.geometry_type=repmat(g.geometry_type,n,1);
T.wave_label=repmat(string(w.wave_label),n,1);
end

function T=assign_prefixed_fields(T,S,prefix)
names=fieldnames(S);
for i=1:numel(names)
    v=S.(names{i});
    if isnumeric(v)&&isscalar(v)
        T.(char(prefix+string(names{i})))=repmat(double(v),height(T),1);
    end
end
end

function T=prepare_external_features(T,cfg,g,w,regime,M,condition_id)
n=height(T);
T.condition_id=condition_id*ones(n,1);
T.condition_label=repmat(g.geometry_id+"__"+regime+"__M"+string(M),n,1);
T.geometry_id=repmat(g.geometry_id,n,1);
T.geometry_name=repmat(g.geometry_name,n,1);
T.geometry_type=repmat(g.geometry_type,n,1);
T.field_regime_label=repmat(regime,n,1);
T.field_regime_variant=repmat(regime+"_N"+string(w.Nwaves),n,1);
T.user_field_guess=map_user_guess(T.field_regime_label);
T.SIM_f0=cfg.f0*ones(n,1); T.SIM_dx=cfg.dx*ones(n,1);
T.SIM_dz=cfg.dz*ones(n,1); T.SIM_cs_bg=g.cs_low*ones(n,1);
T.SIM_cs_inc=g.cs_high*ones(n,1); T.SIM_Nwaves=cfg.Nwaves*ones(n,1);
T.SIM_Is2D=repmat(logical(cfg.Is2D),n,1);
T.SIM_ForceInPlaneWave=repmat(logical(cfg.ForceInPlaneWave),n,1);
T.REQ_M=M*ones(n,1); T.REQ_cs_guess=3*ones(n,1);
T.M_eff_guess=M*ones(n,1); T.step_idx=M*ones(n,1);
T.realization_idx=ones(n,1);
if ~ismember('patch_idx',T.Properties.VariableNames), T.patch_idx=(1:n)'; end
if ~ismember('map_ix',T.Properties.VariableNames), T.map_ix=T.patch_idx; end
if ~ismember('map_iz',T.Properties.VariableNames), T.map_iz=ones(n,1); end
[cs,region,material,distance]=geometry_at_points(T.x_center_m,T.z_center_m,g);
T.cs_true=cs; T.region_name=region; T.material_region=material;
T.distance_to_interface_m=distance;
T.q_true=q_at_known_k(T.req_mapping,T.SIM_f0,T.cs_true);
if ~ismember('global_REQ_Nbins_effective',T.Properties.VariableNames)
    T.global_REQ_Nbins_effective=cellfun(@mapping_nbins,T.global_req_mapping);
end
T=adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T);
T=adaptive_req.analysis.Test12Analysis.addBins(T);
end

function [cs,region,material,distance]=geometry_at_points(x,z,g)
n=numel(x); cs=g.cs_low*ones(n,1); region=strings(n,1); material=strings(n,1);
switch g.geometry_type
    case "homogeneous"
        distance=nan(n,1); region(:)="homogeneous"; material(:)="homogeneous";
    case "inclusion"
        d=hypot(x-g.center(1),z-g.center(2))-g.radius; distance=abs(d);
        inside=d<0; cs(inside)=g.cs_high; material(~inside)="background";
        material(inside)="inclusion"; region(:)="transition_excluded";
        region(d<=-g.far_margin)="inclusion_core";
        region(abs(d)<=g.interface_half_width)="interface_band";
        region(d>=g.far_margin)="background_far";
    case "bilayer"
        d=x-g.center(1); distance=abs(d); high=d>=0; cs(high)=g.cs_high;
        material(~high)="layer_1"; material(high)="layer_2";
        region(:)="transition_excluded";
        region(d<=-g.far_margin)="layer_1_far";
        region(d>=g.far_margin)="layer_2_far";
        region(abs(d)<=g.interface_half_width)="interface_band";
end
end

function q=q_at_known_k(mappings,f0,cs)
q=nan(numel(cs),1);
for i=1:numel(cs)
    m=mappings{i}; k=double(m.k_cent(:)); E=double(m.Ecum(:));
    ok=isfinite(k)&isfinite(E);
    if nnz(ok)>=2
        [ku,ia]=unique(k(ok),'stable'); Eu=E(ok); Eu=Eu(ia);
        q(i)=interp1(ku,Eu,2*pi*f0(i)/cs(i),'linear','extrap');
        q(i)=min(max(q(i),0),1);
    end
end
end

function T=predict_all_models(specs,T_feat,cfg,regime,theory_all_rows)
T=table();
for i=1:numel(specs)
    P=adaptive_req.analysis.predict_q_model_from_table(specs(i).model,T_feat, ...
        'ModelType',specs(i).model_type,'ModelName',specs(i).model_name);
    Ti=prediction_base(T_feat);
    Ti.model_name=repmat(specs(i).model_name,height(Ti),1);
    Ti.feature_set=repmat(specs(i).feature_set,height(Ti),1);
    Ti.model_role=repmat(specs(i).model_role,height(Ti),1);
    Ti.model_type=repmat(specs(i).model_type,height(Ti),1);
    Ti.q_pred_raw=P.q_pred_raw; Ti.q_pred=min(max(P.q_pred,0.001),0.999);
    Ti.cs_pred=q_to_cs(Ti.q_pred,T_feat.req_mapping,cfg.f0);
    T=concat_tables(T,Ti);
end
q=theory_q_for_rows(T_feat,regime,theory_all_rows);
Ti=prediction_base(T_feat);
Ti.model_name=repmat("TheoryQDiscrete",height(Ti),1);
Ti.feature_set=repmat("field_matched_discrete_theory",height(Ti),1);
Ti.model_role=repmat("diagnostic_only",height(Ti),1);
Ti.model_type=repmat("theory_no_ml",height(Ti),1);
Ti.q_pred_raw=q; Ti.q_pred=q; Ti.cs_pred=q_to_cs(q,T_feat.req_mapping,cfg.f0);
T=concat_tables(T,Ti);
end

function T=prediction_base(F)
vars=["condition_id","condition_label","geometry_id","geometry_name", ...
    "geometry_type","field_regime_label","field_regime_variant", ...
    "SIM_f0","SIM_cs_bg","SIM_cs_inc","SIM_dx","SIM_dz","SIM_Nwaves", ...
    "REQ_M","REQ_cs_guess","step_idx","realization_idx","patch_idx", ...
    "patch_label","map_iz","map_ix","cx","cz","x_center_m","z_center_m", ...
    "aperture_value","normalized_aperture","Omega_sr","region_name", ...
    "material_region","distance_to_interface_m","q_true","cs_true"];
vars=vars(ismember(vars,string(F.Properties.VariableNames)));
T=F(:,cellstr(vars));
end

function q=theory_q_for_rows(T,regime,all_rows)
q=nan(height(T),1); cache=containers.Map('KeyType','char','ValueType','double');
for i=1:height(T)
    if ~all_rows && T.step_idx(i)~=1 && T.step_idx(i)~=10, continue; end
    label=string(regime);
    if label=="unknown"
        if T.step_idx(i)==1, label="directional_2D"; else, label="diffuse_3D"; end
    end
    key=sprintf('%s_dx%.12g_f%g_M%g',label,T.SIM_dx(i),T.SIM_f0(i),T.REQ_M(i));
    if isKey(cache,key), q(i)=cache(key); continue; end
    if label=="partial_3D"
        qi=.5*(theory_one(T,i,"Diffuse2D")+theory_one(T,i,"Diffuse3D"));
    elseif label=="directional_2D"
        qi=theory_one(T,i,"SingleWave");
    elseif label=="diffuse_2D"
        qi=theory_one(T,i,"Diffuse2D");
    else
        qi=theory_one(T,i,"Diffuse3D");
    end
    cache(key)=qi; q(i)=qi;
end
end

function q=theory_one(T,i,type)
out=adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    T.SIM_dx(i),T.SIM_dz(i),T.SIM_f0(i),T.REQ_cs_guess(i), ...
    'M',T.REQ_M(i),'Gamma',1,'PadFactor',1,'Nbins','auto', ...
    'SmoothSigma',1,'TheoryMode','S2D','FieldType',type,'Plot',false);
q=out.q_th;
end

function cs=q_to_cs(q,mappings,f0)
cs=nan(numel(q),1);
for i=1:numel(q)
    if isfinite(q(i))&&~isempty(mappings{i})
        cs(i)=adaptive_req.quantile.quantile_to_cs(mappings{i},q(i),f0);
    end
end
end

function T=add_prediction_errors(T)
T.q_error=T.q_pred-T.q_true; T.abs_q_error=abs(T.q_error);
T.sws_error=T.cs_pred-T.cs_true;
T.sws_error_pct=100*T.sws_error./T.cs_true;
T.abs_sws_error_pct=abs(T.sws_error_pct);
T.is_high_error_20=T.abs_sws_error_pct>20;
end

function T=add_q_gradients(T,dx,dz)
T.q_gradient_mag=nan(height(T),1);
G=unique(T(:,{'condition_id','model_name'}),'rows');
for i=1:height(G)
    idx=T.condition_id==G.condition_id(i)&T.model_name==G.model_name(i);
    Ti=T(idx,:);
    if ~all(ismember({'map_iz','map_ix'},Ti.Properties.VariableNames)), continue; end
    Q=map_from_table(Ti,'q_pred');
    [dqz,dqx]=gradient(Q,dz,dx);
    lin=sub2ind(size(Q),Ti.map_iz,Ti.map_ix);
    grad=hypot(dqx,dqz); T.q_gradient_mag(idx)=grad(lin);
end
end

function Q=map_from_table(T,var)
Q=nan(max(T.map_iz),max(T.map_ix));
Q(sub2ind(size(Q),T.map_iz,T.map_ix))=T.(var);
end

function [T,sweep]=run_homogeneous_aperture_sweep(AP,cs,M,ci,mi)
cfg=aperture_cfg(AP,cs,ci,mi);
feat=adaptive_req.config.default_feature_config('M',M,'cs_guess',AP.cs_guess, ...
    'gamma_win',1,'pad_factor',1);
req={'Nbins','auto','Nbins_auto_oversample',1,'Nbins_min',16,'smooth_sigma',1};
[T,sweep]=adaptive_req.studies.run_aperture_sweep(cfg,feat, ...
    'SamplingMode','cone','NumSteps',AP.num_steps,'NumRealizations',1, ...
    'NumPatches',5,'SeedBase',AP.seed_base+1000*ci+100*mi, ...
    'PatchOptions',{'Pattern','manual','ManualCenters',AP.patch_centers}, ...
    'ReqOptions',req,'StoreWavefields',true,'StoreReqCurve',false, ...
    'StoreReqMapping',true,'ComputeGlobalReq',true, ...
    'StoreGlobalReqMapping',true,'StoreReqMetadata',true, ...
    'StoreFeatureStruct',false,'Verbose',true);
end

function cfg=aperture_cfg(AP,cs,ci,mi)
cfg=adaptive_req.config.default_sim_config('Lx',AP.Lx,'Lz',AP.Lz, ...
    'dx',AP.dx,'dz',AP.dz,'f0',AP.f0,'cs_bg',cs,'Nwaves',AP.Nwaves, ...
    'WaveModel','planewave','SourceSampling','cone', ...
    'AngularSamplingMethod','fibonacci','ForceInPlaneWave',true, ...
    'SNR',Inf,'AmpJitter',0, ...
    'Seed',AP.seed_base+1000*ci+100*mi,'UseParfor',true);
cfg.Is2D=false;
cfg.DecayAlpha=0;
end

function T=prepare_aperture_features(T,AP,cs,M)
n=height(T); T.condition_id=(100000+100*cs+M)*ones(n,1);
T.condition_label=repmat("aperture_cs"+string(cs)+"_M"+string(M),n,1);
T.geometry_id=repmat("homogeneous_cs"+string(cs),n,1);
T.geometry_name=repmat("Homogeneous c_s="+string(cs)+" m/s",n,1);
T.geometry_type=repmat("homogeneous_aperture_sweep",n,1);
T.field_regime_label=repmat("aperture_continuous",n,1);
T.field_regime_variant=repmat("cone_aperture_N2000",n,1);
T.user_field_guess=categorical(repmat("unknown",n,1), ...
    ["directional_like","partially_diffuse","diffuse_like","unknown"]);
T.SIM_cs_bg=cs*ones(n,1); T.SIM_cs_inc=cs*ones(n,1);
T.REQ_M=M*ones(n,1); T.REQ_cs_guess=AP.cs_guess*ones(n,1);
T.M_eff_guess=M*ones(n,1); T.cs_true=cs*ones(n,1);
T.normalized_aperture=T.aperture_value/360;
T.patch_label=AP.patch_names(T.patch_idx);
T.x_center_m=(T.cx-1)*AP.dx; T.z_center_m=(T.cz-1)*AP.dz;
T.q_true=T.q_theory;
if ~ismember('global_REQ_Nbins_effective',T.Properties.VariableNames)
    T.global_REQ_Nbins_effective=cellfun(@mapping_nbins,T.global_req_mapping);
end
T=adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T);
T=adaptive_req.analysis.Test12Analysis.addBins(T);
end

function labels=map_user_guess(regime)
r=string(regime); labels=strings(size(r));
labels(r=="directional_2D")="directional_like";
labels(r=="partial_3D")="partially_diffuse";
labels(r=="diffuse_2D"|r=="diffuse_3D")="diffuse_like";
labels(labels=="")="unknown";
labels=categorical(labels,["directional_like","partially_diffuse", ...
    "diffuse_like","unknown"]);
end

function n=mapping_nbins(m)
if isfield(m,'Nbins_effective'), n=double(m.Nbins_effective); else, n=numel(m.Ecum); end
end

function Tsum=summarize_sws(T,groups)
[G,Tsum]=findgroups(T(:,cellstr(groups)));
Tsum.N=splitapply(@finite_count,T.abs_sws_error_pct,G);
Tsum.MAPE_pct=splitapply(@(x)mean(x,'omitnan'),T.abs_sws_error_pct,G);
Tsum.RMSE_pct=splitapply(@(x)sqrt(mean(x.^2,'omitnan')),T.sws_error_pct,G);
Tsum.MedAE_pct=splitapply(@(x)median(x,'omitnan'),T.abs_sws_error_pct,G);
Tsum.bias_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_error_pct,G);
Tsum.HighError_gt20_pct=splitapply(@high_error_rate, ...
    T.abs_sws_error_pct,G);
Tsum=sortrows(Tsum,'MAPE_pct');
end

function Tsum=summarize_q(T,groups)
[G,Tsum]=findgroups(T(:,cellstr(groups)));
Tsum.N=splitapply(@finite_pair_count,T.q_true,T.q_pred,G);
Tsum.MAE_q=splitapply(@(x)mean(x,'omitnan'),T.abs_q_error,G);
Tsum.RMSE_q=splitapply(@(x)sqrt(mean(x.^2,'omitnan')),T.q_error,G);
Tsum.bias_q=splitapply(@(x)mean(x,'omitnan'),T.q_error,G);
Tsum=sortrows(Tsum,'MAE_q');
end

function Tsum=summarize_binned(T,var,bin_name)
x=T.(var); finite_x=x(isfinite(x));
if isempty(finite_x), Tsum=table(); return; end
edges=linspace(min(finite_x),max(finite_x),9);
if edges(1)==edges(end), edges=[edges(1)-eps(edges(1)) edges(1)+eps(edges(1))]; end
T.(bin_name)=discretize(x,edges,'IncludedEdge','right');
groups=["model_name","geometry_id","field_regime_label","REQ_M",string(bin_name)];
[G,Tsum]=findgroups(T(:,cellstr(groups)));
Tsum.N=splitapply(@numel,T.abs_sws_error_pct,G);
Tsum.MAPE_pct=splitapply(@(v)mean(v,'omitnan'),T.abs_sws_error_pct,G);
Tsum.bias_pct=splitapply(@(v)mean(v,'omitnan'),T.sws_error_pct,G);
Tsum.bin_center=splitapply(@(v)mean(v,'omitnan'),x,G);
if string(bin_name)=="distance_bin_mm", Tsum.bin_center=Tsum.bin_center*1e3; end
end

function plot_external_metric(T,var,ylab,file,OUT)
models=ordered_model_names(T.model_name); cases=unique(T.geometry_id,'stable');
finite_all=T.(var)(isfinite(T.(var)));
global_top=max(finite_all); if global_top<=0, global_top=1; end
fig=figure('Color','w','Units','centimeters','Position',[2 2 30 17]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for ci=1:numel(cases)
    ax=nexttile(tl); Tc=T(T.geometry_id==cases(ci),:);
    Tc.x_group=Tc.field_regime_label+" | M="+string(Tc.REQ_M);
    xgroups=ordered_external_groups(Tc);
    C=nan(numel(models),numel(xgroups));
    for mi=1:numel(models)
        for ri=1:numel(xgroups)
            idx=Tc.model_name==models(mi)&Tc.x_group==xgroups(ri);
            C(mi,ri)=mean(Tc.(var)(idx),'omitnan');
        end
    end
    imagesc(ax,C); colormap(ax,parula(256)); cb=colorbar(ax);
    cb.FontSize=6; cb.Label.String=ylab; cb.Label.FontSize=7;
    xticks(ax,1:numel(xgroups)); xticklabels(ax,short_external_groups(Tc,xgroups));
    yticks(ax,1:numel(models)); yticklabels(ax,short_model_names(models));
    xtickangle(ax,45); ax.FontSize=6.5; ax.TickLabelInterpreter='none';
    title(ax,short_geometry_name(cases(ci)),'FontSize',9,'FontWeight','normal');
    xlabel(ax,'Field regime and M','FontSize',7);
    clim(ax,[0 global_top]);
end
title(tl,ylab+" by external case, regime, and M", ...
    'FontSize',10,'FontWeight','normal');
export_fig(fig,fullfile(OUT.fig_dir,file),OUT);
end

function plot_external_maps(CASES,OUT)
for ci=1:numel(CASES)
    T=CASES(ci).T_pred; models=unique(T.model_name,'stable');
    fig=figure('Color','w','Units','centimeters','Position',[2 2 21 27]);
    tl=tiledlayout(numel(models),3,'TileSpacing','tight','Padding','compact');
    q_limits=[0 1];
    true_cs=T.cs_true(isfinite(T.cs_true));
    sws_limits=[max(0,min(true_cs)-0.5),max(true_cs)+0.5];
    err_abs=abs(T.sws_error_pct(isfinite(T.sws_error_pct)));
    err_limit=max(10,prctile(err_abs,98));
    for mi=1:numel(models)
        Ti=T(T.model_name==models(mi),:);
        vals={Ti.q_pred,Ti.cs_pred,Ti.sws_error_pct};
        limits={q_limits,sws_limits,[-err_limit err_limit]};
        labels={'q','SWS (m/s)','signed error (%)'};
        for j=1:3
            ax=nexttile(tl); A=map_values(Ti,vals{j}); imagesc(ax,A);
            axis(ax,'image','off'); set(ax,'YDir','normal'); clim(ax,limits{j});
            if j==3, colormap(ax,turbo(256)); else, colormap(ax,parula(256)); end
            if mi==1
                title(ax,labels{j},'FontSize',8,'FontWeight','normal');
                cb=colorbar(ax); cb.FontSize=5;
            end
            if j==1
                text(ax,-0.04,0.5,short_model_name(models(mi)), ...
                    'Units','normalized','HorizontalAlignment','right', ...
                    'VerticalAlignment','middle','FontSize',6.5, ...
                    'Interpreter','none','Clipping','off');
            end
        end
    end
    title(tl,strrep(CASES(ci).key,'_',' '), ...
        'Interpreter','none','FontSize',10,'FontWeight','normal');
    export_fig(fig,fullfile(OUT.map_dir,"level20_maps__"+ ...
        sanitize_filename(CASES(ci).key)+".png"),OUT);
end
end

function A=map_values(T,v)
A=nan(max(T.map_iz),max(T.map_ix)); A(sub2ind(size(A),T.map_iz,T.map_ix))=v;
end

function plot_binned_error(T,xvar,xlab,ttl,file,OUT)
if isempty(T), return; end
fig=figure('Color','w','Units','centimeters','Position',[2 2 19 11]);
ax=axes(fig); hold(ax,'on'); models=unique(T.model_name,'stable');
for i=1:numel(models)
    Ti=T(T.model_name==models(i),:);
    S=groupsummary(Ti,xvar,'mean',{'bin_center','MAPE_pct'});
    [x,ord]=sort(S.mean_bin_center);
    plot(ax,x,S.mean_MAPE_pct(ord),'-o','LineWidth',1, ...
        'MarkerSize',4,'DisplayName',short_model_name(models(i)));
end
xlabel(ax,xlab,'FontSize',8); ylabel(ax,'SWS MAPE (%)','FontSize',8);
title(ax,ttl,'FontSize',9,'FontWeight','normal'); ax.FontSize=7.5;
legend(ax,'Location','eastoutside','Interpreter','none','FontSize',6.5); grid(ax,'on');
export_fig(fig,fullfile(OUT.fig_dir,file),OUT);
end

function plot_q_vs_aperture(T,OUT)
models=unique(T.model_name,'stable');
fig=figure('Color','w','Units','centimeters','Position',[2 2 30 22]);
tl=tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for cs=[2 3 4]
    for M=[2 3 4]
        ax=nexttile(tl); hold(ax,'on'); S=T(T.SIM_cs_bg==cs&T.REQ_M==M,:);
        truth=S(S.model_name==models(1),:);
        truth=groupsummary(truth,'normalized_aperture','mean','q_true');
        plot(ax,truth.normalized_aperture,truth.mean_q_true,'k-','LineWidth',2, ...
            'DisplayName','q true');
        for i=1:numel(models)
            Si=S(S.model_name==models(i),:);
            y=groupsummary(Si,'normalized_aperture','mean','q_pred');
            plot(ax,y.normalized_aperture,y.mean_q_pred,'-o','DisplayName',models(i));
        end
        title(ax,sprintf('c_s=%g, M=%g',cs,M)); xlabel(ax,'normalized aperture');
        ylabel(ax,'q'); grid(ax,'on');
    end
end
legend(nexttile(tl,1),'Location','eastoutside','Interpreter','none');
title(tl,'q versus aperture (mean over five fixed patches)','FontWeight','normal');
export_fig(fig,fullfile(OUT.fig_dir,'q_vs_aperture_by_patch.png'),OUT);
end

function plot_aperture_error(T,var,ylab,file,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 25 15]);
ax=axes(fig); hold(ax,'on'); models=unique(T.model_name,'stable');
for i=1:numel(models)
    Ti=T(T.model_name==models(i),:);
    S=groupsummary(Ti,'normalized_aperture','mean',var);
    plot(ax,S.normalized_aperture,S.("mean_"+var),'-o','DisplayName',models(i));
end
xlabel(ax,'normalized aperture'); ylabel(ax,ylab); yline(ax,0,'k:'); grid(ax,'on');
legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.fig_dir,file),OUT);
end

function plot_aperture_q_scatter(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 24 18]);
models=unique(T.model_name,'stable'); tl=tiledlayout(3,3,'TileSpacing','compact');
for i=1:numel(models)
    ax=nexttile(tl); Ti=T(T.model_name==models(i)&isfinite(T.q_pred),:);
    scatter(ax,Ti.q_true,Ti.q_pred,8,Ti.normalized_aperture,'filled'); hold(ax,'on');
    plot(ax,[0 1],[0 1],'k--'); axis(ax,'square'); xlim(ax,[0 1]); ylim(ax,[0 1]);
    title(ax,models(i),'Interpreter','none'); xlabel(ax,'q true'); ylabel(ax,'q predicted');
end
title(tl,'Aperture sweep q tracking','FontWeight','normal');
export_fig(fig,fullfile(OUT.fig_dir,'q_true_vs_q_pred_aperture_sweep.png'),OUT);
end

function plot_aperture_patch_map(AP,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 14 12]); ax=axes(fig);
rectangle(ax,'Position',[0 0 AP.Lx*1e3 AP.Lz*1e3]); hold(ax,'on'); axis(ax,'equal');
x=(AP.patch_centers(:,1)-1)*AP.dx*1e3; z=(AP.patch_centers(:,2)-1)*AP.dz*1e3;
scatter(ax,x,z,80,'filled'); text(ax,x+1,z,AP.patch_names,'Interpreter','none');
set(ax,'YDir','reverse'); xlim(ax,[0 AP.Lx*1e3]); ylim(ax,[0 AP.Lz*1e3]);
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)'); title(ax,'Five fixed aperture-sweep patches');
export_fig(fig,fullfile(OUT.fig_dir,'aperture_sweep_five_patches.png'),OUT);
end

function print_external_summary(T)
fprintf('\nBest model by external case:\n');
cases=unique(T.geometry_id,'stable');
for i=1:numel(cases)
    Ti=T(T.geometry_id==cases(i)&T.model_role~="diagnostic_only",:);
    S=groupsummary(Ti,'model_name','mean','MAPE_pct'); [v,j]=min(S.mean_MAPE_pct);
    fprintf('  %-22s %-48s %.3f%% MAPE\n',cases(i),S.model_name(j),v);
end
fprintf('\nBest model by field regime:\n');
reg=unique(T.field_regime_label,'stable');
for i=1:numel(reg)
    Ti=T(T.field_regime_label==reg(i)&T.model_role~="diagnostic_only",:);
    S=groupsummary(Ti,'model_name','mean','MAPE_pct'); [v,j]=min(S.mean_MAPE_pct);
    fprintf('  %-16s %-48s %.3f%% MAPE\n',reg(i),S.model_name(j),v);
end
end

function print_old_new_summary(T)
is_new=contains(T.model_name,"_T18"); is_old=endsWith(T.model_name,"_old");
fprintf('\nT18 versus old deployments (mean grouped MAPE):\n');
fprintf('  T18 models: %.3f%%\n',mean(T.MAPE_pct(is_new),'omitnan'));
fprintf('  Old models: %.3f%%\n',mean(T.MAPE_pct(is_old),'omitnan'));
end

function print_aperture_trend_summary(T)
fprintf('\nAperture q-trend correlations (mean trajectory):\n');
models=unique(T.model_name,'stable');
for i=1:numel(models)
    Ti=T(T.model_name==models(i)&isfinite(T.q_pred),:);
    S=groupsummary(Ti,'step_idx','mean',{'q_true','q_pred'});
    r=corr(S.mean_q_true,S.mean_q_pred,'Rows','complete');
    fprintf('  %-48s r=%.3f\n',models(i),r);
end
end

function print_user_guess_summary(Te,Ta)
no="HybridLocalGlobal_T18_noUserRegime";
yes="HybridLocalGlobal_T18_withUserRegimeGuess";
fprintf('\nUser regime guess contribution:\n');
fprintf('  External delta MAPE (with - without): %.3f percentage points\n', ...
    mean(Te.MAPE_pct(Te.model_name==yes),'omitnan')- ...
    mean(Te.MAPE_pct(Te.model_name==no),'omitnan'));
fprintf('  Aperture delta MAPE with unknown guess: %.3f percentage points\n', ...
    mean(Ta.MAPE_pct(Ta.model_name==yes),'omitnan')- ...
    mean(Ta.MAPE_pct(Ta.model_name==no),'omitnan'));
end

function labels=short_model_names(names)
labels=strings(size(names));
for i=1:numel(names), labels(i)=short_model_name(names(i)); end
end

function names=ordered_model_names(available)
preferred=["LocalOnly_T18";"GlobalOnly_T18"; ...
    "HybridLocalGlobal_T18_noUserRegime"; ...
    "HybridLocalGlobal_T18_withUserRegimeGuess"; ...
    "LocalOnly_old";"GlobalOnly_old";"HybridLocalGlobal_old"; ...
    "TheoryQDiscrete"];
available=unique(string(available),'stable');
names=preferred(ismember(preferred,available));
names=[names;available(~ismember(available,names))];
end

function label=short_model_name(name)
switch string(name)
    case "LocalOnly_T18", label="Local T18";
    case "GlobalOnly_T18", label="Global T18";
    case "HybridLocalGlobal_T18_noUserRegime", label="Hybrid T18";
    case "HybridLocalGlobal_T18_withUserRegimeGuess", label="Hybrid T18 + guess";
    case "LocalOnly_old", label="Local old";
    case "GlobalOnly_old", label="Global old";
    case "HybridLocalGlobal_old", label="Hybrid old";
    case "TheoryQDiscrete", label="Discrete theory";
    otherwise, label=string(name);
end
end

function labels=short_external_groups(T,xgroups)
labels=strings(size(xgroups));
for i=1:numel(xgroups)
    row=find(T.x_group==xgroups(i),1);
    switch string(T.field_regime_label(row))
        case "directional_2D", regime="Dir 2D";
        case "diffuse_2D", regime="Diff 2D";
        case "partial_3D", regime="Partial 3D";
        case "diffuse_3D", regime="Diff 3D";
        otherwise, regime=string(T.field_regime_label(row));
    end
    labels(i)=regime+" M"+string(T.REQ_M(row));
end
end


function groups=ordered_external_groups(T)
regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
M_values=sort(unique(T.REQ_M));
groups=strings(0,1);
for ri=1:numel(regimes)
    for mi=1:numel(M_values)
        candidate=regimes(ri)+" | M="+string(M_values(mi));
        if any(T.x_group==candidate), groups(end+1,1)=candidate; end %#ok<AGROW>
    end
end
extra=unique(T.x_group(~ismember(T.x_group,groups)),'stable');
groups=[groups;extra(:)];
end

function label=short_geometry_name(name)
switch string(name)
    case "homogeneous_cs2", label="Homogeneous 2 m/s";
    case "homogeneous_cs3", label="Homogeneous 3 m/s";
    case "bilayer_2_3", label="Bilayer 2/3 m/s";
    case "inclusion_2_3", label="Circular inclusion 2/3 m/s";
    otherwise, label=strrep(string(name),'_',' ');
end
end

function n=finite_count(x)
n=sum(isfinite(x));
end

function n=finite_pair_count(a,b)
n=sum(isfinite(a)&isfinite(b));
end

function pct=high_error_rate(x)
x=x(isfinite(x));
if isempty(x), pct=NaN; else, pct=100*mean(x>20); end
end

function export_fig(fig,file,~)
exportgraphics(fig,file,'Resolution',220); close(fig);
end

function T=remove_cell_columns(T)
is_cell=varfun(@iscell,T,'OutputFormat','uniform'); T(:,is_cell)=[];
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
if isempty(B), T=A; return; end
vars=unique([string(A.Properties.VariableNames),string(B.Properties.VariableNames)],'stable');
A=add_missing(A,vars); B=add_missing(B,vars);
T=[A(:,cellstr(vars));B(:,cellstr(vars))];
end

function T=add_missing(T,vars)
for i=1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)), continue; end
    T.(char(vars(i)))=nan(height(T),1);
end
end

function name=sanitize_filename(name)
name=regexprep(char(string(name)),'[^A-Za-z0-9_-]+','_');
end
