%% analyze_test_22_confidence_external_validation.m
% Test 22: external validation of the frozen Test 21 confidence detectors.
%
% This script DOES NOT train or recalibrate any detector or q model. It loads
% the seven frozen Test 21 detectors and applies them to clean external
% simulations spanning geometry, frequency, field regime, REQ M, and spatial
% resolution. Oracle quantities are added only after detector inference.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST22_SIZE      = pilot | medium | full (default pilot)
%   ADAPTIVE_REQ_TEST22_FAST_MODE = true | false (pilot-only debug subset)
%   ADAPTIVE_REQ_TEST22_MODE      = full | validate_only

clear; clc; close all;
format compact;

%% Fixed external-validation design

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontName','Helvetica', ...
    'defaultTextFontName','Helvetica','defaultAxesFontSize',8, ...
    'defaultTextFontSize',8,'defaultLegendFontSize',7);

TEST_SIZE=lower(string(getenv('ADAPTIVE_REQ_TEST22_SIZE')));
if TEST_SIZE=="", TEST_SIZE="pilot"; end
switch TEST_SIZE
    case "pilot"
        DESIGN.f0=500; DESIGN.M=3;
        DESIGN.dx=[0.2 0.5]*1e-3; TARGET_STEP_M=1.0e-3;
    case "medium"
        DESIGN.f0=[400 500 600]; DESIGN.M=[2 3 4];
        DESIGN.dx=[0.2 0.5]*1e-3; TARGET_STEP_M=0.75e-3;
    case "full"
        DESIGN.f0=[300 400 500 600]; DESIGN.M=[2 3 4];
        DESIGN.dx=[0.2 0.3 0.5]*1e-3; TARGET_STEP_M=0.5e-3;
    otherwise
        error('ADAPTIVE_REQ_TEST22_SIZE must be pilot, medium, or full.');
end
FAST_MODE=any(strcmpi(strtrim(getenv('ADAPTIVE_REQ_TEST22_FAST_MODE')), ...
    {'true','1','yes','on'}));
if FAST_MODE && TEST_SIZE~="pilot"
    warning('Test22:FastModeIgnored', ...
        'ADAPTIVE_REQ_TEST22_FAST_MODE is only active for pilot.');
    FAST_MODE=false;
end
DESIGN.cs_guess=3;
DESIGN.Lx=0.05; DESIGN.Lz=0.05;
DESIGN.seed=22001;
DESIGN.risk_thresholds=[0.2 0.5 0.8];
DESIGN.confidence_edges=0:0.1:1;
DESIGN.distance_edges_mm=[0 0.25 0.5 1 2 3 5 8 12 Inf];
DESIGN.primary_model="HybridLocalGlobal_T18_noUserRegime";
DESIGN.optional_model="HybridLocalGlobal_T18_withUserRegimeGuess";
DESIGN.test_size=TEST_SIZE;
DESIGN.fast_mode=FAST_MODE;
RUN_MODE=lower(string(getenv('ADAPTIVE_REQ_TEST22_MODE')));
if RUN_MODE=="", RUN_MODE="full"; end
assert(ismember(RUN_MODE,["full","validate_only"]), ...
    'ADAPTIVE_REQ_TEST22_MODE must be full or validate_only.');

CFG17=adaptive_req.config.load_profile_config( ...
    'test_17_synthetic_inclusion_kWave_like','RootDir',root_dir);
CFG17.REQ.cs_guess=DESIGN.cs_guess;
CFG17.REQ.EdgeMode='valid';
GEOMETRIES=build_geometry_specs(CFG17);
REGIMES=CFG17.CASES;
OUT=make_output_dirs(root_dir,TEST_SIZE,FAST_MODE);

fprintf('\nTest 22: frozen confidence-detector external validation\n');
fprintf('No detector or q model is trained in this script.\n');
fprintf('Size: %s | fast mode: %d | target REQ step: %.2f mm.\n', ...
    TEST_SIZE,FAST_MODE,1e3*TARGET_STEP_M);
fprintf('Design: %d simulation fields and %d REQ conditions.\n', ...
    numel(GEOMETRIES)*numel(DESIGN.f0)*numel(REGIMES)*numel(DESIGN.dx), ...
    numel(GEOMETRIES)*numel(DESIGN.f0)*numel(REGIMES)* ...
    numel(DESIGN.dx)*numel(DESIGN.M));

%% Load frozen models and detectors, continuing when optional artifacts lack

if FAST_MODE
    requested_models=DESIGN.primary_model;
else
    requested_models=strings(0,1);
end
MODEL_SPECS=load_available_q_models(root_dir,requested_models);
fprintf('Available q models: %s\n',strjoin(string({MODEL_SPECS.model_name}),', '));
if ~any([MODEL_SPECS.model_name]==DESIGN.primary_model)
    warning('Test22:PrimaryModelMissing', ...
        ['Primary model %s is unavailable. Test 22 cannot form the trained ', ...
        'confidence features; no simulations will be run.'],DESIGN.primary_model);
    return;
end
DETECTORS=load_frozen_detectors(root_dir);
if FAST_MODE
    keep=ismember(string({DETECTORS.detector_name}), ...
        ["ML_bagged_trees","ML_logistic_regression","Rule_q_gradient"]);
    DETECTORS=DETECTORS(keep);
end
fprintf('Available frozen detectors: %s\n', ...
    strjoin(string({DETECTORS.detector_name}),', '));
if isempty(DETECTORS)
    warning('Test22:NoDetectors','No frozen Test 21 detectors are available.');
    return;
end
validate_q_model_predictors(MODEL_SPECS);
validate_detector_predictors(DETECTORS);
fprintf('Loaded %d q models and %d frozen confidence detectors.\n', ...
    numel(MODEL_SPECS),numel(DETECTORS));
if FAST_MODE
    fprintf(['Fast pilot: unavailable disagreement inputs are NaN and use ', ...
        'the frozen detectors'' training-median imputation. Debug only.\n']);
end
if RUN_MODE=="validate_only"
    fprintf(['Test 22 validation-only check passed: artifacts loaded, ', ...
        'predictor policies passed, and no simulation was generated.\n']);
    return;
end

%% Simulate, predict q/SWS, and apply detectors

T_timing=table();
num_conditions=numel(GEOMETRIES)*numel(DESIGN.f0)*numel(REGIMES)* ...
    numel(DESIGN.dx)*numel(DESIGN.M);
EXTERNAL_PARTS=cell(num_conditions,1);
condition_dir=fullfile(OUT.data_dir,'conditions');
if exist(condition_dir,'dir')~=7, mkdir(condition_dir); end
field_cache_dir=fullfile(root_dir,'outputs', ...
    'test_22_confidence_external_validation','analysis','data','field_cache');
if exist(field_cache_dir,'dir')~=7, mkdir(field_cache_dir); end

condition_id=0;
for gi=1:numel(GEOMETRIES)
    g=GEOMETRIES(gi);
    for fi=1:numel(DESIGN.f0)
        f0=DESIGN.f0(fi);
        for ri=1:numel(REGIMES)
            w=REGIMES(ri); regime=canonical_regime_label(w);
            for di=1:numel(DESIGN.dx)
                dx=DESIGN.dx(di);
                CFG17.REQ.StepX=max(1,round(TARGET_STEP_M/dx));
                CFG17.REQ.StepZ=max(1,round(TARGET_STEP_M/dx));
                sim_key=sprintf('%s__f%g__%s__dx%gum',g.geometry_id,f0, ...
                    lower(regime),round(1e6*dx));
                sim_file=fullfile(field_cache_dir,"field__"+sanitize(sim_key)+".mat");
                sim_tic=tic;
                if exist(sim_file,'file')==2
                    S=load(sim_file,'sim','cfg_sim'); sim=S.sim; cfg_sim=S.cfg_sim;
                    sim_source="cache";
                else
                    cfg_sim=build_sim_cfg(DESIGN,g,w,f0,dx,gi,fi,ri,di);
                    fprintf('\nGenerating %s\n',sim_key);
                    sim=adaptive_req.simulate.run_single_simulation(cfg_sim);
                    save(sim_file,'sim','cfg_sim','-v7.3');
                    sim_source="generated";
                end
                sim_seconds=toc(sim_tic);
                fprintf('Field %s in %.2f s; size %d x %d.\n', ...
                    sim_source,sim_seconds,size(sim.Uxz,1),size(sim.Uxz,2));

                for mi=1:numel(DESIGN.M)
                    M=DESIGN.M(mi); condition_id=condition_id+1;
                    key=sim_key+"__M"+string(M)+"__step"+ ...
                        string(CFG17.REQ.StepX)+"px";
                    result_file=fullfile(condition_dir,"result__"+sanitize(key)+".mat");
                    if exist(result_file,'file')==2
                        S=load(result_file,'T_condition');
                        EXTERNAL_PARTS{condition_id}=S.T_condition;
                        fprintf('Reused %s\n',key);
                        continue;
                    end

                    condition_tic=tic;
                    fprintf(['Condition %d: %s | dx=%.3f mm | M=%g | ', ...
                        'StepX=%d | StepZ=%d | field=%dx%d\n'],condition_id,key, ...
                        1e3*dx,M,CFG17.REQ.StepX,CFG17.REQ.StepZ, ...
                        size(sim.Uxz,1),size(sim.Uxz,2));
                    stage_tic=tic;
                    [T_feat,~]=extract_feature_table(sim,cfg_sim,CFG17,M);
                    extract_seconds=toc(stage_tic);
                    T_feat=prepare_features(T_feat,cfg_sim,g,w,regime,M,condition_id);
                    n_windows=height(T_feat);
                    fprintf('  REQ/features: %d windows in %.2f s.\n', ...
                        n_windows,extract_seconds);
                    stage_tic=tic;
                    T_models=predict_available_models(MODEL_SPECS,T_feat,cfg_sim, ...
                        regime,~FAST_MODE);
                    q_seconds=toc(stage_tic);
                    stage_tic=tic;
                    D=build_operational_detector_table(T_feat,T_models,DESIGN);
                    detector_feature_seconds=toc(stage_tic);
                    stage_tic=tic;
                    R=apply_frozen_detectors(DETECTORS,D);
                    detector_seconds=toc(stage_tic);
                    T_condition=make_external_predictions(D,R,T_models,DESIGN);
                    stage_tic=tic;
                    save(result_file,'T_condition','-v7.3');
                    save_seconds=toc(stage_tic);
                    total_seconds=toc(condition_tic);
                    fprintf(['  q models %.2f s | detector features %.2f s | ', ...
                        'detectors %.2f s | save %.2f s | total %.2f s.\n'], ...
                        q_seconds,detector_feature_seconds,detector_seconds, ...
                        save_seconds,total_seconds);
                    timing_row=table(condition_id,g.geometry_id,f0,regime,M,dx, ...
                        CFG17.REQ.StepX,CFG17.REQ.StepZ,size(sim.Uxz,1), ...
                        size(sim.Uxz,2),n_windows,sim_seconds,extract_seconds, ...
                        q_seconds,detector_feature_seconds,detector_seconds, ...
                        save_seconds,total_seconds, ...
                        'VariableNames',{'condition_id','geometry','f0', ...
                        'field_regime','M','dx','StepX','StepZ','field_nz', ...
                        'field_nx','n_req_windows','simulation_seconds', ...
                        'extract_seconds','q_prediction_seconds', ...
                        'detector_feature_seconds','detector_seconds', ...
                        'save_seconds','total_condition_seconds'});
                    T_timing=concat_tables(T_timing,timing_row);
                    EXTERNAL_PARTS{condition_id}=T_condition;
                    clear T_feat T_models D R T_condition;
                end
                clear sim;
            end
        end
    end
end

assert(all(~cellfun(@isempty,EXTERNAL_PARTS)), ...
    'One or more Test 22 conditions did not produce predictions.');
T_external=vertcat(EXTERNAL_PARTS{:});
clear EXTERNAL_PARTS;
assert(~isempty(T_external),'Test 22 produced no external predictions.');
writetable(T_external,fullfile(OUT.table_dir,'level22_external_predictions.csv'));
if ~isempty(T_timing)
    writetable(T_timing,fullfile(OUT.table_dir,'level22_condition_timing.csv'));
end
save(fullfile(OUT.data_dir,'level22_external_predictions.mat'), ...
    'T_external','DESIGN','GEOMETRIES','-v7.3');

%% Metrics and diagnostic summaries

T_metrics=evaluate_grouped(T_external,strings(0,1),DESIGN.risk_thresholds);
T_by_geometry=evaluate_grouped(T_external,"geometry",DESIGN.risk_thresholds);
T_by_frequency=evaluate_grouped(T_external,"f0",DESIGN.risk_thresholds);
T_by_regime=evaluate_grouped(T_external,"field_regime",DESIGN.risk_thresholds);
T_by_M=evaluate_grouped(T_external,"M",DESIGN.risk_thresholds);
T_by_dx=evaluate_grouped(T_external,"dx",DESIGN.risk_thresholds);
T_by_geometry_regime_M=evaluate_grouped(T_external, ...
    ["geometry" "field_regime" "M"],DESIGN.risk_thresholds);
T_threshold=threshold_summary(T_external,DESIGN.risk_thresholds);
T_confidence=confidence_bin_summary(T_external,DESIGN.confidence_edges);
T_distance=distance_summary(T_external,DESIGN.distance_edges_mm);

writetable(T_metrics,fullfile(OUT.table_dir,'level22_detector_metrics_overall.csv'));
writetable(T_by_geometry,fullfile(OUT.table_dir,'level22_detector_metrics_by_geometry.csv'));
writetable(T_by_frequency,fullfile(OUT.table_dir,'level22_detector_metrics_by_frequency.csv'));
writetable(T_by_regime,fullfile(OUT.table_dir,'level22_detector_metrics_by_regime.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'level22_detector_metrics_by_M.csv'));
writetable(T_by_dx,fullfile(OUT.table_dir,'level22_detector_metrics_by_dx.csv'));
writetable(T_by_geometry_regime_M,fullfile(OUT.table_dir, ...
    'level22_detector_metrics_by_geometry_regime_M.csv'));
writetable(T_distance,fullfile(OUT.table_dir,'level22_distance_to_interface_summary.csv'));
writetable(T_confidence,fullfile(OUT.table_dir,'level22_confidence_bin_summary.csv'));
writetable(T_threshold,fullfile(OUT.table_dir,'level22_threshold_summary.csv'));
save(fullfile(OUT.data_dir,'level22_analysis_results.mat'),'T_metrics', ...
    'T_by_geometry','T_by_frequency','T_by_regime','T_by_M','T_by_dx', ...
    'T_by_geometry_regime_M','T_distance','T_confidence','T_threshold','-v7.3');

%% Figures

plot_representative_maps(T_external,OUT);
plot_confidence_curves(T_confidence,OUT,'mean_abs_error_pct', ...
    'Mean absolute SWS error (%)','level22_error_vs_confidence.png');
plot_confidence_curves(T_confidence,OUT,'high_error_rate_gt20', ...
    'Observed high-error >20 rate','level22_high_error_rate_vs_confidence.png');
plot_binary_curves(T_external,OUT,"roc");
plot_binary_curves(T_external,OUT,"pr");
plot_auc_summary(T_metrics,T_threshold,OUT);
plot_heatmaps(T_external,OUT);
plot_distance_analysis(T_distance,OUT);
plot_calibration(T_confidence,OUT);

%% Automatic scientific summary

print_summary(T_external,T_metrics,T_by_geometry,T_by_frequency, ...
    T_by_regime,T_distance,T_threshold,DESIGN);
fprintf('\nTest 22 complete. Analysis folder:\n%s\n',OUT.analysis_dir);

%% Local functions

function OUT=make_output_dirs(root_dir,test_size,fast_mode)
OUT.analysis_dir=fullfile(root_dir,'outputs', ...
    'test_22_confidence_external_validation','analysis');
if test_size~="full"
    suffix=test_size;
    if fast_mode, suffix=suffix+"_fast"; end
    OUT.analysis_dir=fullfile(OUT.analysis_dir,char(suffix));
end
OUT.table_dir=fullfile(OUT.analysis_dir,'tables');
OUT.figure_dir=fullfile(OUT.analysis_dir,'figures');
OUT.data_dir=fullfile(OUT.analysis_dir,'data');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs), if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end, end
end

function specs=build_geometry_specs(CFG)
specs=[make_geometry("homogeneous_cs2","homogeneous",2,2,CFG); ...
    make_geometry("homogeneous_cs3","homogeneous",3,3,CFG); ...
    make_geometry("bilayer_2_3","bilayer",2,3,CFG); ...
    make_geometry("circular_inclusion_2_3","inclusion",2,3,CFG)];
end

function g=make_geometry(id,type,low,high,CFG)
g=struct('geometry_id',string(id),'geometry_type',string(type), ...
    'cs_low',low,'cs_high',high,'center',CFG.INCLUSION.Center, ...
    'radius',CFG.INCLUSION.Radius,'interface_half_width', ...
    CFG.INCLUSION.InterfaceHalfWidth,'far_margin', ...
    CFG.INCLUSION.BackgroundFarMargin);
if string(type)=="bilayer", g.center=[0.025 0.025]; g.radius=NaN; end
end

function specs=load_available_q_models(root_dir,requested_names)
new_dir=fullfile(root_dir,'outputs','model_registry','test19_clean_field_regime');
old_dir=fullfile(root_dir,'outputs','model_registry','test12_hybrid_baseline');
requests=[ ...
    request(new_dir,"LocalOnly_T18","CleanFieldRegime_noUser","LocalOnly_T18","operational"); ...
    request(new_dir,"GlobalOnly_T18","CleanFieldRegime_noUser","GlobalOnly_T18","operational"); ...
    request(new_dir,"HybridLocalGlobal_T18_noUserRegime","CleanFieldRegime_noUser", ...
        "HybridLocalGlobal_T18_noUserRegime","operational"); ...
    request(new_dir,"HybridLocalGlobal_T18_withUserRegimeGuess", ...
        "CleanFieldRegime_withUserGuess","HybridLocalGlobal_T18_withUserRegimeGuess","user_informed"); ...
    request(old_dir,"LocalOnly","NoCsGuess","LocalOnly_old","external_old_model"); ...
    request(old_dir,"GlobalOnly","NoCsGuess","GlobalOnly_old","external_old_model"); ...
    request(old_dir,"HybridLocalGlobal","WithCsGuess","HybridLocalGlobal_old","external_old_model")];
specs=struct('model_name',{},'feature_set',{},'model_type',{}, ...
    'model_role',{},'model_file',{},'model',{});
for i=1:numel(requests)
    if ~isempty(requested_names) && ...
            ~ismember(requests(i).report,string(requested_names))
        continue;
    end
    try
        [M,I,F]=adaptive_req.analysis.load_q_model_deployment(requests(i).folder, ...
            'ModelName',requests(i).stored,'FeatureSet',requests(i).feature, ...
            'ModelType','bagged_trees');
        specs(end+1)=struct('model_name',requests(i).report, ...
            'feature_set',string(I.feature_set),'model_type',"bagged_trees", ...
            'model_role',requests(i).role,'model_file',string(F),'model',M); %#ok<AGROW>
    catch ME
        warning('Test22:ModelUnavailable','Skipping q model %s: %s', ...
            requests(i).report,ME.message);
    end
end
end

function r=request(folder,stored,feature,report,role)
r=struct('folder',string(folder),'stored',string(stored), ...
    'feature',string(feature),'report',string(report),'role',string(role));
end

function D=load_frozen_detectors(root_dir)
file=fullfile(root_dir,'outputs','test_21_interface_confidence_detector', ...
    'analysis','models','level21_detector_training_checkpoint.mat');
D=struct([]);
if exist(file,'file')~=2
    warning('Test22:DetectorFileMissing','Frozen detector checkpoint missing: %s',file);
    return;
end
S=load(file,'DETECTORS');
if ~isfield(S,'DETECTORS')
    warning('Test22:DetectorVariableMissing','DETECTORS not found in %s.',file); return;
end
D=S.DETECTORS;
expected=["Rule_q_gradient","Rule_q_hybrid_minus_local", ...
    "Rule_q_hybrid_minus_theory","Rule_max_local_theory_disagreement", ...
    "ML_logistic_regression","ML_bagged_trees","ML_boosted_trees"];
found=string({D.detector_name});
for i=1:numel(expected)
    if ~ismember(expected(i),found)
        warning('Test22:DetectorUnavailable','Frozen detector unavailable: %s',expected(i));
    end
end
end

function validate_q_model_predictors(specs)
forbidden=lower(["q_true","q_theory","cs_true","distance_to_interface", ...
    "sws_error","abs_sws_error","high_error","m_eff_true_diag"]);
for i=1:numel(specs)
    p=lower(string({specs(i).model.encoder.entries.name}));
    assert(isempty(intersect(p,forbidden)) && ...
        ~any(contains(p,"error")|contains(p,"target")), ...
        'Oracle predictor in q model %s.',specs(i).model_name);
end
end

function validate_detector_predictors(detectors)
banned=lower(["q_true","q_theory","cs_true","true_sws", ...
    "distance_to_interface","sws_error","abs_sws_error", ...
    "high_error_gt20","high_error_gt10","mixed_window_fraction"]);
for i=1:numel(detectors)
    p=lower(string(detectors(i).predictors));
    assert(isempty(intersect(p,banned)) && ~any(contains(p,"error")), ...
        'Oracle predictor in frozen detector %s.',detectors(i).detector_name);
end
end

function label=canonical_regime_label(w)
switch string(w.wave_label)
    case "directional_2d", label="directional_2D";
    case "diffuse_2d", label="diffuse_2D";
    case "partial_diffuse_3d", label="partial_3D";
    case "diffuse_3d", label="diffuse_3D";
    otherwise, error('Unknown regime %s.',w.wave_label);
end
end

function cfg=build_sim_cfg(D,g,w,f0,dx,gi,~,ri,~)
cfg=struct('Lx',D.Lx,'Lz',D.Lz,'dx',dx,'dz',dx,'f0',f0, ...
    'cs_bg',g.cs_low,'cs_inc',g.cs_high,'Nwaves',w.Nwaves, ...
    'Is2D',logical(w.Is2D),'WaveModel',char(w.WaveModel), ...
    'AngularSamplingMethod',char(w.AngularSamplingMethod), ...
    'ForceInPlaneWave',logical(w.ForceInPlaneWave),'SNR',Inf, ...
    'AmpJitter',0,'DecayAlpha',0,'UseParfor',true, ...
    'Seed',D.seed+100000*gi+10000*ri+10*round(f0)+round(1e6*dx), ...
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

function [T,O]=extract_feature_table(sim,cfg,CFG,M)
feat=adaptive_req.config.default_feature_config('M',M, ...
    'cs_guess',CFG.REQ.cs_guess,'gamma_win',CFG.REQ.Gamma, ...
    'pad_factor',CFG.REQ.PadFactor);
req={'Nbins',CFG.REQ.Nbins,'Nbins_auto_oversample', ...
    CFG.REQ.Nbins_auto_oversample,'Nbins_min',CFG.REQ.Nbins_min, ...
    'smooth_sigma',CFG.REQ.SmoothSigma};
[req_cfg,feat]=adaptive_req.config.default_req_config(cfg,feat,req{:});
cfg_req=cfg; cfg_req.UseParfor=false;
[qg,curve,fg]=adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz,cfg_req,req_cfg,feat);
sg=adaptive_req.quantile.extract_ecum_shape_features(curve);
O=adaptive_req.estimators.req_estimator_map(sim.Uxz,cfg_req,feat, ...
    'StepX',CFG.REQ.StepX,'StepZ',CFG.REQ.StepZ, ...
    'EdgeMode',CFG.REQ.EdgeMode,'QuantileMode','local_req', ...
    'ReqOptions',req,'ReturnFeatures',false,'ReturnFeatureTable',true, ...
    'ReuseReqSpectrumForFeatures',true,'UseWindowParfor',true, ...
    'StoreReqCurves',false,'Verbose',false);
T=O.feature_table; n=height(T);
T.global_req_mapping=repmat({adaptive_req.quantile.make_req_mapping(curve)},n,1);
T.q_global_req=qg*ones(n,1);
T.global_REQ_Nbins_effective=curve.Nbins_effective*ones(n,1);
T=assign_prefixed(T,fg,"global_"); T=assign_prefixed(T,sg,"global_");
end

function T=assign_prefixed(T,S,prefix)
names=fieldnames(S);
for i=1:numel(names)
    v=S.(names{i});
    if isnumeric(v)&&isscalar(v), T.(char(prefix+string(names{i})))=repmat(double(v),height(T),1); end
end
end

function T=prepare_features(T,cfg,g,w,regime,M,id)
n=height(T); T.condition_id=id*ones(n,1);
T.geometry_id=repmat(g.geometry_id,n,1);
T.geometry_type=repmat(g.geometry_type,n,1);
T.field_regime_label=repmat(regime,n,1);
T.field_regime_variant=repmat(regime+"_N"+string(w.Nwaves),n,1);
T.user_field_guess=map_user_guess(T.field_regime_label);
T.SIM_f0=cfg.f0*ones(n,1); T.SIM_dx=cfg.dx*ones(n,1); T.SIM_dz=cfg.dz*ones(n,1);
T.SIM_cs_bg=g.cs_low*ones(n,1); T.SIM_cs_inc=g.cs_high*ones(n,1);
T.SIM_Nwaves=cfg.Nwaves*ones(n,1); T.SIM_Is2D=repmat(logical(cfg.Is2D),n,1);
T.SIM_ForceInPlaneWave=repmat(logical(cfg.ForceInPlaneWave),n,1);
T.REQ_M=M*ones(n,1); T.REQ_cs_guess=3*ones(n,1); T.M_eff_guess=M*ones(n,1);
T.step_idx=M*ones(n,1); T.realization_idx=ones(n,1);
if ~ismember('patch_idx',T.Properties.VariableNames), T.patch_idx=(1:n)'; end
if ~ismember('map_ix',T.Properties.VariableNames), T.map_ix=T.patch_idx; end
if ~ismember('map_iz',T.Properties.VariableNames), T.map_iz=ones(n,1); end
[cs,region,distance]=geometry_at_points(T.x_center_m,T.z_center_m,g);
T.cs_true=cs; T.region_name=region; T.distance_to_interface_m=distance;
T.q_true=q_at_known_k(T.req_mapping,T.SIM_f0,T.cs_true);
T=adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T);
T=adaptive_req.analysis.Test12Analysis.addBins(T);
end

function labels=map_user_guess(regime)
labels=strings(size(regime)); labels(regime=="directional_2D")="directional_like";
labels(regime=="partial_3D")="partially_diffuse";
labels(regime=="diffuse_2D"|regime=="diffuse_3D")="diffuse_like";
labels(labels=="")="unknown";
labels=categorical(labels,["directional_like","partially_diffuse","diffuse_like","unknown"]);
end

function [cs,region,distance]=geometry_at_points(x,z,g)
n=numel(x); cs=g.cs_low*ones(n,1); region=strings(n,1);
switch g.geometry_type
    case "homogeneous", distance=nan(n,1); region(:)="homogeneous_no_interface";
    case "inclusion"
        d=hypot(x-g.center(1),z-g.center(2))-g.radius; distance=abs(d);
        cs(d<0)=g.cs_high; region(:)="transition";
        region(d<=-g.far_margin)="inclusion_core";
        region(abs(d)<=g.interface_half_width)="interface_band";
        region(d>=g.far_margin)="background_far";
    case "bilayer"
        d=x-g.center(1); distance=abs(d); cs(d>=0)=g.cs_high; region(:)="transition";
        region(d<=-g.far_margin)="layer_1_far";
        region(d>=g.far_margin)="layer_2_far";
        region(abs(d)<=g.interface_half_width)="interface_band";
end
end

function q=q_at_known_k(mappings,f0,cs)
q=nan(numel(cs),1);
for i=1:numel(cs)
    m=mappings{i}; k=double(m.k_cent(:)); E=double(m.Ecum(:)); ok=isfinite(k)&isfinite(E);
    if nnz(ok)>=2
        [ku,ia]=unique(k(ok),'stable'); Eu=E(ok); Eu=Eu(ia);
        q(i)=min(max(interp1(ku,Eu,2*pi*f0(i)/cs(i),'linear','extrap'),0),1);
    end
end
end

function P=predict_available_models(specs,F,cfg,regime,include_theory)
P=table();
for i=1:numel(specs)
    try
        Q=adaptive_req.analysis.predict_q_model_from_table(specs(i).model,F, ...
            'ModelType',specs(i).model_type,'ModelName',specs(i).model_name);
        T=prediction_base(F); T.model_name=repmat(specs(i).model_name,height(F),1);
        T.q_pred=min(max(Q.q_pred,0.001),0.999);
        T.sws_pred=q_to_cs(T.q_pred,F.req_mapping,cfg.f0);
        P=concat_tables(P,T);
    catch ME
        warning('Test22:PredictionFailed','Skipping %s for this condition: %s', ...
            specs(i).model_name,ME.message);
    end
end
if include_theory
    try
        q=theory_q(F,regime); T=prediction_base(F);
        T.model_name=repmat("TheoryQDiscrete",height(F),1); T.q_pred=q;
        T.sws_pred=q_to_cs(q,F.req_mapping,cfg.f0); P=concat_tables(P,T);
    catch ME
        warning('Test22:TheoryUnavailable','TheoryQDiscrete skipped: %s',ME.message);
    end
end
end

function T=prediction_base(F)
T=F(:,{'condition_id','map_iz','map_ix','q_true','cs_true'});
end

function q=theory_q(T,regime)
switch string(regime)
    case "directional_2D", type="SingleWave";
    case "diffuse_2D", type="Diffuse2D";
    case "diffuse_3D", type="Diffuse3D";
    otherwise, type="Partial3D";
end
if type=="Partial3D"
    q=.5*(theory_one(T,"Diffuse2D")+theory_one(T,"Diffuse3D"));
else
    q=theory_one(T,type);
end
q=repmat(q,height(T),1);
end

function q=theory_one(T,type)
o=adaptive_req.theory.q_theory_REQ_discrete_shearUZ(T.SIM_dx(1),T.SIM_dz(1), ...
    T.SIM_f0(1),T.REQ_cs_guess(1),'M',T.REQ_M(1),'Gamma',1, ...
    'PadFactor',1,'Nbins','auto','SmoothSigma',1,'TheoryMode','S2D', ...
    'FieldType',type,'Plot',false); q=o.q_th;
end

function cs=q_to_cs(q,mappings,f0)
cs=nan(size(q));
for i=1:numel(q)
    if isfinite(q(i)), cs(i)=adaptive_req.quantile.quantile_to_cs(mappings{i},q(i),f0); end
end
end

function D=build_operational_detector_table(F,P,DESIGN)
H=sortrows(P(P.model_name==DESIGN.primary_model,:),{'map_iz','map_ix'});
assert(height(H)==height(F),'Primary prediction rows are incomplete.');
F=sortrows(F,{'map_iz','map_ix'}); D=F;
D.cs_guess=D.REQ_cs_guess;
D.dx=D.SIM_dx;
D.dz=D.SIM_dz;
D.q_pred_hybrid_T18=H.q_pred; D.sws_pred_hybrid_T18=H.sws_pred;
map_names=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","HybridLocalGlobal_old","TheoryQDiscrete"];
qvars=["q_pred_local_T18","q_pred_global_T18","q_pred_hybrid_T18_with_guess", ...
    "q_pred_hybrid_old","q_theory_discrete"];
svars=["sws_pred_local_T18","sws_pred_global_T18","sws_pred_hybrid_T18_with_guess", ...
    "sws_pred_hybrid_old","sws_theory_discrete"];
for i=1:numel(map_names)
    Ti=sortrows(P(P.model_name==map_names(i),:),{'map_iz','map_ix'});
    if height(Ti)==height(D), D.(qvars(i))=Ti.q_pred; D.(svars(i))=Ti.sws_pred;
    else, D.(qvars(i))=nan(height(D),1); D.(svars(i))=nan(height(D),1); end
end
D.abs_q_hybrid_minus_local=abs(D.q_pred_hybrid_T18-D.q_pred_local_T18);
D.abs_q_hybrid_minus_global=abs(D.q_pred_hybrid_T18-D.q_pred_global_T18);
D.abs_q_hybrid_minus_theory=abs(D.q_pred_hybrid_T18-D.q_theory_discrete);
D.abs_q_hybrid_minus_hybrid_guess=abs(D.q_pred_hybrid_T18-D.q_pred_hybrid_T18_with_guess);
D.abs_q_hybrid_T18_minus_old=abs(D.q_pred_hybrid_T18-D.q_pred_hybrid_old);
D.abs_sws_hybrid_minus_local=abs(D.sws_pred_hybrid_T18-D.sws_pred_local_T18);
D.abs_sws_hybrid_minus_global=abs(D.sws_pred_hybrid_T18-D.sws_pred_global_T18);
D.abs_sws_hybrid_minus_theory=abs(D.sws_pred_hybrid_T18-D.sws_theory_discrete);
D.abs_sws_hybrid_T18_minus_old=abs(D.sws_pred_hybrid_T18-D.sws_pred_hybrid_old);
D.max_local_theory_disagreement=max(D.abs_q_hybrid_minus_local,D.abs_q_hybrid_minus_theory);
D=add_map_features(D,3);
end

function D=add_map_features(D,nhood)
Q=map_from_rows(D,D.q_pred_hybrid_T18); C=map_from_rows(D,D.sws_pred_hybrid_T18);
[dqz,dqx]=gradient(Q,D.SIM_dz(1),D.SIM_dx(1));
[dcz,dcx]=gradient(C,D.SIM_dz(1),D.SIM_dx(1));
[qs,qr]=local_stats(Q,nhood); [cs,cr]=local_stats(C,nhood);
lin=sub2ind(size(Q),D.map_iz,D.map_ix); qg=hypot(dqx,dqz); cg=hypot(dcx,dcz);
D.grad_q_hybrid=qg(lin); D.grad_sws_hybrid=cg(lin);
D.local_std_q_hybrid=qs(lin); D.local_std_sws_hybrid=cs(lin);
D.local_range_q_hybrid=qr(lin); D.local_range_sws_hybrid=cr(lin);
end

function [S,R]=local_stats(A,nhood)
K=ones(nhood); valid=isfinite(A); A0=A; A0(~valid)=0;
n=conv2(double(valid),K,'same'); mu=conv2(A0,K,'same')./max(n,1);
v=conv2(A0.^2,K,'same')./max(n,1)-mu.^2; S=sqrt(max(v,0));
mx=movmax(movmax(A,nhood,1,'omitnan'),nhood,2,'omitnan');
mn=movmin(movmin(A,nhood,1,'omitnan'),nhood,2,'omitnan'); R=mx-mn;
end

function R=apply_frozen_detectors(detectors,D)
R=table();
for i=1:numel(detectors)
    det=detectors(i);
    try
        if det.detector_family=="ml"
            names=string(det.predictors); require_variables(D,names,det.detector_name);
            X=double(D{:,cellstr(names)}); med=det.impute_median;
            for j=1:size(X,2), bad=~isfinite(X(:,j)); X(bad,j)=med(j); end
            X=(X-det.predictor_mu)./det.predictor_sigma;
        else
            require_variables(D,det.rule_feature,det.detector_name);
            X=double(D.(det.rule_feature)); X(~isfinite(X))=det.impute_median;
        end
        [~,score]=predict(det.model,X); p=positive_score(det.model,score);
        T=table(repmat(string(det.detector_name),height(D),1),p,1-p, ...
            'VariableNames',{'detector','risk','confidence'});
        R=concat_tables(R,T);
    catch ME
        warning('Test22:DetectorPredictionFailed','Skipping detector %s: %s', ...
            det.detector_name,ME.message);
    end
end
end

function require_variables(T,names,label)
missing=setdiff(string(names),string(T.Properties.VariableNames));
assert(isempty(missing),'%s missing predictors: %s',label,strjoin(missing,', '));
end

function p=positive_score(model,scores)
classes=string(model.ClassNames); j=find(classes=="1"|classes=="true",1);
if isempty(j), j=size(scores,2); end
p=min(max(double(scores(:,j)),0),1);
end

function T=make_external_predictions(D,R,P,DESIGN)
T=table(); detectors=unique(R.detector,'stable');
eval_models=[DESIGN.primary_model DESIGN.optional_model];
for mi=1:numel(eval_models)
    E=sortrows(P(P.model_name==eval_models(mi),:),{'map_iz','map_ix'});
    if height(E)~=height(D)
        continue;
    end
    err=100*(E.sws_pred-D.cs_true)./D.cs_true; ae=abs(err);
    for di=1:numel(detectors)
        Ri=R(R.detector==detectors(di),:);
        B=table(); B.condition_id=D.condition_id; B.geometry=D.geometry_id;
        B.geometry_type=D.geometry_type; B.f0=D.SIM_f0;
        B.field_regime=D.field_regime_label; B.M=D.REQ_M; B.dx=D.SIM_dx; B.dz=D.SIM_dz;
        B.map_iz=D.map_iz; B.map_ix=D.map_ix; B.x=D.x_center_m; B.z=D.z_center_m;
        B.sws_model=repmat(eval_models(mi),height(D),1);
        B.detector=repmat(detectors(di),height(D),1);
        B.q_pred=E.q_pred; B.SWS_pred=E.sws_pred; B.true_SWS=D.cs_true;
        B.SWS_error_pct=err; B.abs_SWS_error_pct=ae;
        B.high_error_gt20=ae>20; B.high_error_gt10=ae>10;
        B.risk=Ri.risk; B.confidence=Ri.confidence;
        B.low_confidence_risk_gt020=B.risk>0.2;
        B.low_confidence_risk_gt050=B.risk>0.5;
        B.low_confidence_risk_gt080=B.risk>0.8;
        B.q_gradient=D.grad_q_hybrid;
        B.distance_to_interface_mm=1e3*D.distance_to_interface_m;
        T=concat_tables(T,B);
    end
end
end

function T=evaluate_grouped(P,groups,thresholds)
base=["sws_model" "detector"]; groups=string(groups(:))'; allgroups=[base groups];
[G,K]=findgroups(P(:,cellstr(allgroups))); T=K;
T.N=splitapply(@numel,P.risk,G);
T.ROC_AUC_gt20=splitapply(@(y,p)auc_value(y,p,"roc"),P.high_error_gt20,P.risk,G);
T.PR_AUC_gt20=splitapply(@(y,p)auc_value(y,p,"pr"),P.high_error_gt20,P.risk,G);
T.ROC_AUC_gt10=splitapply(@(y,p)auc_value(y,p,"roc"),P.high_error_gt10,P.risk,G);
T.PR_AUC_gt10=splitapply(@(y,p)auc_value(y,p,"pr"),P.high_error_gt10,P.risk,G);
for i=1:numel(thresholds)
    tag=sprintf('%03d',round(100*thresholds(i)));
    T.("recall_risk_gt"+tag)=splitapply(@(y,p)recall_at(y,p,thresholds(i)),P.high_error_gt20,P.risk,G);
    T.("precision_risk_gt"+tag)=splitapply(@(y,p)precision_at(y,p,thresholds(i)),P.high_error_gt20,P.risk,G);
    T.("fraction_low_confidence_risk_gt"+tag)=splitapply(@(p)mean(p>thresholds(i),'omitnan'),P.risk,G);
end
T.mean_error_high_confidence=splitapply(@(e,p)mean(e(p<=0.2),'omitnan'),P.abs_SWS_error_pct,P.risk,G);
T.mean_error_low_confidence=splitapply(@(e,p)mean(e(p>0.2),'omitnan'),P.abs_SWS_error_pct,P.risk,G);
T.high_error_rate_high_confidence=splitapply(@(y,p)mean(y(p<=0.2),'omitnan'),P.high_error_gt20,P.risk,G);
T.high_error_rate_low_confidence=splitapply(@(y,p)mean(y(p>0.2),'omitnan'),P.high_error_gt20,P.risk,G);
end

function A=auc_value(y,p,type)
valid=isfinite(p); y=logical(y(valid)); p=p(valid);
if isempty(y)||~any(y)||~any(~y), A=NaN; return; end
[fpr,tpr,r,prec]=binary_curves(y,p);
if type=="roc", A=trapz(fpr,tpr); else, A=trapz(r,prec); end
end

function [fpr,tpr,r,prec]=binary_curves(y,p)
[~,ord]=sort(p,'descend'); y=double(y(ord)); tp=cumsum(y); fp=cumsum(1-y);
tpr=[0;tp/max(sum(y),1);1]; fpr=[0;fp/max(sum(1-y),1);1];
r=[0;tp/max(sum(y),1)]; prec=[1;tp./max(tp+fp,1)];
end

function v=recall_at(y,p,t), pred=p>t; v=safe_div(sum(pred&y),sum(y)); end
function v=precision_at(y,p,t), pred=p>t; v=safe_div(sum(pred&y),sum(pred)); end
function v=safe_div(a,b), if b==0, v=NaN; else, v=a/b; end, end

function T=threshold_summary(P,thresholds)
T=table(); names=unique(P(:,{'sws_model','detector'}),'rows');
for i=1:height(names)
    X=P(P.sws_model==names.sws_model(i)&P.detector==names.detector(i),:);
    for j=1:numel(thresholds)
        t=thresholds(j); low=X.risk>t; R=names(i,:); R.risk_threshold=t; R.N=height(X);
        R.recall=recall_at(X.high_error_gt20,X.risk,t); R.precision=precision_at(X.high_error_gt20,X.risk,t);
        R.fraction_low_confidence=mean(low); R.mean_error_high_confidence=mean(X.abs_SWS_error_pct(~low),'omitnan');
        R.mean_error_low_confidence=mean(X.abs_SWS_error_pct(low),'omitnan');
        R.high_error_rate_high_confidence=mean(X.high_error_gt20(~low),'omitnan');
        R.high_error_rate_low_confidence=mean(X.high_error_gt20(low),'omitnan');
        T=concat_tables(T,R);
    end
end
end

function T=confidence_bin_summary(P,edges)
P.confidence_bin=discretize(P.confidence,edges,'IncludedEdge','right');
P.confidence_bin(P.confidence==0)=1;
groups={'sws_model','detector','confidence_bin'}; [G,T]=findgroups(P(:,groups));
T.N=splitapply(@numel,P.confidence,G);
T.mean_confidence=splitapply(@(x)mean(x,'omitnan'),P.confidence,G);
T.mean_predicted_risk=splitapply(@(x)mean(x,'omitnan'),P.risk,G);
T.mean_abs_error_pct=splitapply(@(x)mean(x,'omitnan'),P.abs_SWS_error_pct,G);
T.high_error_rate_gt20=splitapply(@(x)mean(x,'omitnan'),P.high_error_gt20,G);
T.high_error_rate_gt10=splitapply(@(x)mean(x,'omitnan'),P.high_error_gt10,G);
end

function T=distance_summary(P,edges)
P=P(ismember(P.geometry_type,["bilayer","inclusion"]),:);
P.distance_bin=discretize(P.distance_to_interface_mm,edges);
groups={'sws_model','detector','geometry','distance_bin'}; [G,T]=findgroups(P(:,groups));
T.N=splitapply(@numel,P.risk,G);
T.mean_distance_mm=splitapply(@(x)mean(x,'omitnan'),P.distance_to_interface_mm,G);
T.mean_abs_error_pct=splitapply(@(x)mean(x,'omitnan'),P.abs_SWS_error_pct,G);
T.mean_confidence=splitapply(@(x)mean(x,'omitnan'),P.confidence,G);
T.high_error_rate_gt20=splitapply(@(x)mean(x,'omitnan'),P.high_error_gt20,G);
end

function plot_representative_maps(P,OUT)
cases=["homogeneous_cs2","homogeneous_cs3","bilayer_2_3","bilayer_2_3", ...
    "circular_inclusion_2_3","circular_inclusion_2_3"];
regimes=["directional_2D","diffuse_3D","directional_2D","diffuse_3D", ...
    "directional_2D","diffuse_3D"];
primary=P.sws_model=="HybridLocalGlobal_T18_noUserRegime";
[~,best]=max_group_auc(P(primary,:));
fig=figure('Color','w','Units','centimeters','Position',[1 1 30 25]);
tl=tiledlayout(6,6,'TileSpacing','tight','Padding','compact');
for i=1:6
    X=P(primary&P.detector==best&P.geometry==cases(i)&P.field_regime==regimes(i)& ...
        P.f0==500&P.M==3&abs(P.dx-0.2e-3)<eps,:);
    if isempty(X), continue; end
    vals={X.true_SWS,X.SWS_pred,X.abs_SWS_error_pct,X.confidence, ...
        X.distance_to_interface_mm,X.q_gradient};
    labels={'true SWS','predicted SWS','absolute error (%)','confidence','distance (mm)','|grad q|'};
    for j=1:6
        ax=nexttile(tl); imagesc(ax,map_from_rows(X,vals{j})); axis(ax,'image','off'); colorbar(ax,'FontSize',5);
        if i==1, title(ax,labels{j},'FontWeight','normal','FontSize',7); end
    end
end
title(tl,'Representative external confidence maps','FontWeight','normal');
export_fig(fig,fullfile(OUT.figure_dir,'level22_representative_confidence_maps.png'));
end

function [A,best]=max_group_auc(P)
names=unique(P.detector,'stable'); A=nan(size(names));
for i=1:numel(names), X=P(P.detector==names(i),:); A(i)=auc_value(X.high_error_gt20,X.risk,"pr"); end
[~,j]=max(A); best=names(j);
end

function plot_confidence_curves(T,OUT,var,ylab,file)
fig=figure('Color','w','Units','centimeters','Position',[2 2 20 13]); ax=axes(fig); hold(ax,'on');
names=unique(T.detector,'stable'); model="HybridLocalGlobal_T18_noUserRegime";
for i=1:numel(names)
    X=T(T.detector==names(i)&T.sws_model==model,:); plot(ax,X.mean_confidence,X.(var),'-o','MarkerSize',3,'DisplayName',short_detector(names(i)));
end
xlabel(ax,'Confidence'); ylabel(ax,ylab); grid(ax,'on'); legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,file));
end

function plot_binary_curves(P,OUT,type)
P=P(P.sws_model=="HybridLocalGlobal_T18_noUserRegime",:);
fig=figure('Color','w','Units','centimeters','Position',[2 2 19 13]); ax=axes(fig); hold(ax,'on'); names=unique(P.detector,'stable');
for i=1:numel(names)
    X=P(P.detector==names(i),:); [fpr,tpr,r,prec]=binary_curves(X.high_error_gt20,X.risk);
    if type=="roc", x=fpr; y=tpr; a=trapz(fpr,tpr); else, x=r; y=prec; a=trapz(r,prec); end
    plot(ax,x,y,'LineWidth',1.2,'DisplayName',short_detector(names(i))+sprintf(' (%.3f)',a));
end
if type=="roc", plot(ax,[0 1],[0 1],'k:'); xlabel(ax,'False-positive rate'); ylabel(ax,'Recall'); file='level22_external_roc_curves.png';
else, xlabel(ax,'Recall'); ylabel(ax,'Precision'); file='level22_external_pr_curves.png'; end
grid(ax,'on'); xlim(ax,[0 1]); ylim(ax,[0 1]); legend(ax,'Location','eastoutside','Interpreter','none'); export_fig(fig,fullfile(OUT.figure_dir,file));
end

function plot_auc_summary(M,T,OUT)
M=M(M.sws_model=="HybridLocalGlobal_T18_noUserRegime",:); T=T(T.sws_model=="HybridLocalGlobal_T18_noUserRegime"&T.risk_threshold==0.2,:);
[~,loc]=ismember(M.detector,T.detector); C=[M.PR_AUC_gt20 M.ROC_AUC_gt20 T.recall(loc) T.fraction_low_confidence(loc)];
fig=figure('Color','w','Units','centimeters','Position',[2 2 24 13]); bar(C); xticks(1:height(M)); xticklabels(short_detector(M.detector)); xtickangle(30); ylim([0 1]); grid on;
legend({'PR AUC','ROC AUC','recall risk>0.2','fraction low confidence'},'Location','eastoutside');
export_fig(fig,fullfile(OUT.figure_dir,'level22_auc_summary.png'));
end

function plot_heatmaps(P,OUT)
P=P(P.sws_model=="HybridLocalGlobal_T18_noUserRegime",:); geoms=unique(P.geometry,'stable');
for gi=1:numel(geoms)
    X=P(P.geometry==geoms(gi),:); specs={["field_regime" "M"],["f0" "M"],["dx" "M"]}; tags=["regime_M","f0_M","dx_M"];
    for si=1:3
        row=string(X.detector); col=string(X.(specs{si}(1)))+" | M="+string(X.M); rn=unique(row,'stable'); cn=unique(col,'stable'); C=nan(numel(rn),numel(cn));
        for r=1:numel(rn), for c=1:numel(cn), idx=row==rn(r)&col==cn(c); C(r,c)=auc_value(X.high_error_gt20(idx),X.risk(idx),"pr"); end, end
        fig=figure('Color','w','Units','centimeters','Position',[2 2 24 12]); imagesc(C,[0 1]); colorbar; yticks(1:numel(rn)); yticklabels(short_detector(rn)); xticks(1:numel(cn)); xticklabels(cn); xtickangle(45); title(geoms(gi)+' PR AUC');
        export_fig(fig,fullfile(OUT.figure_dir,"level22_heatmap_"+geoms(gi)+"_"+tags(si)+".png"));
    end
end
end

function plot_distance_analysis(T,OUT)
model="HybridLocalGlobal_T18_noUserRegime"; names=unique(T.detector,'stable'); geoms=unique(T.geometry,'stable');
for gi=1:numel(geoms)
    fig=figure('Color','w','Units','centimeters','Position',[2 2 25 9]); tl=tiledlayout(1,3,'TileSpacing','compact'); vars=["mean_abs_error_pct","mean_confidence","high_error_rate_gt20"];
    for vi=1:3, ax=nexttile(tl); hold(ax,'on'); for di=1:numel(names), X=T(T.sws_model==model&T.geometry==geoms(gi)&T.detector==names(di),:); plot(ax,X.mean_distance_mm,X.(vars(vi)),'-o','MarkerSize',3,'DisplayName',short_detector(names(di))); end, xlabel(ax,'Distance to interface (mm)'); ylabel(ax,strrep(vars(vi),'_',' ')); grid(ax,'on'); end
    legend(nexttile(tl,3),'Location','eastoutside','Interpreter','none'); title(tl,geoms(gi),'Interpreter','none'); export_fig(fig,fullfile(OUT.figure_dir,"level22_distance_analysis_"+geoms(gi)+".png"));
end
end

function plot_calibration(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 20 13]); ax=axes(fig); hold(ax,'on'); names=unique(T.detector,'stable');
for i=1:numel(names), X=T(T.detector==names(i)&T.sws_model=="HybridLocalGlobal_T18_noUserRegime",:); plot(ax,X.mean_predicted_risk,X.high_error_rate_gt20,'-o','MarkerSize',3,'DisplayName',short_detector(names(i))); end
plot(ax,[0 1],[0 1],'k:'); xlabel(ax,'Predicted risk'); ylabel(ax,'Observed high-error rate'); grid(ax,'on'); legend(ax,'Location','eastoutside','Interpreter','none'); export_fig(fig,fullfile(OUT.figure_dir,'level22_risk_calibration.png'));
end

function print_summary(P,M,Mg,Mf,Mr,D,T,DESIGN)
M=M(M.sws_model==DESIGN.primary_model,:); [~,i]=max(M.PR_AUC_gt20); best=M(i,:);
T=T(T.sws_model==DESIGN.primary_model&T.risk_threshold==0.2,:); [~,i]=max(T.recall); bestrec=T(i,:);
fprintf('\n================ Test 22 summary ================\n');
fprintf('Best global PR AUC: %s (%.4f).\n',short_detector(best.detector),best.PR_AUC_gt20);
fprintf('Best recall at risk > 0.2: %s (%.4f).\n',short_detector(bestrec.detector),bestrec.recall);
for f=[400 600]
    if any(Mf.f0==f)
        report_domain(Mf,Mf.f0==f,best.detector,"f0="+string(f));
    else
        fprintf('f0=%g: not included in %s profile.\n',f,DESIGN.test_size);
    end
end
for r=["diffuse_2D","partial_3D","diffuse_3D"], report_domain(Mr,Mr.field_regime==r,best.detector,r); end
for g=["bilayer_2_3","circular_inclusion_2_3"], report_domain(Mg,Mg.geometry==g,best.detector,g); end
X=P(P.sws_model==DESIGN.primary_model&P.detector==best.detector,:);
fail=summarize_failures(X); fprintf('Worst case: %s.\n',fail);
Db=D(D.sws_model==DESIGN.primary_model&D.detector==best.detector,:);
near=mean(Db.mean_confidence(Db.mean_distance_mm<=1),'omitnan'); far=mean(Db.mean_confidence(Db.mean_distance_mm>=3),'omitnan');
fprintf('Interface concentration: mean confidence near=%.3f, far=%.3f.\n',near,far);
bycase=evaluate_grouped(X,["geometry" "f0" "field_regime" "M" "dx"],DESIGN.risk_thresholds);
isint=ismember(bycase.geometry,["bilayer_2_3","circular_inclusion_2_3"]);
ready=mean(bycase.PR_AUC_gt20(isint)>0.80,'omitnan')>0.5 && ...
    best.mean_error_high_confidence<best.mean_error_low_confidence && near<far;
homlow=mean(X.risk(contains(X.geometry,"homogeneous"))>0.2,'omitnan');
if ready && homlow<0.5
    recommendation="ready for regularization";
elseif mean(bycase.PR_AUC_gt20(isint)>0.80,'omitnan')<=0.5
    recommendation="needs retraining with broader interface data";
else
    recommendation="works only for some regimes";
end
fprintf('Automatic recommendation: %s.\n',recommendation);
if DESIGN.test_size~="full"
    fprintf('Recommendation scope: preliminary %s profile, not full validation.\n', ...
        DESIGN.test_size);
end
fprintf('Homogeneous fraction marked low confidence at risk>0.2: %.3f.\n',homlow);
fprintf('No oracle or diagnostic variable was used for detector inference.\n');
fprintf('=================================================\n');
end

function report_domain(T,idx,det,label)
X=T(idx&T.detector==det,:); fprintf('%s: PR AUC %.3f, recall@0.2 %.3f.\n',label,mean(X.PR_AUC_gt20,'omitnan'),mean(X.recall_risk_gt020,'omitnan'));
end

function s=summarize_failures(X)
T=evaluate_grouped(X,["geometry" "f0" "field_regime" "M" "dx"],[0.2 0.5 0.8]); T=sortrows(T,'PR_AUC_gt20','ascend');
if isempty(T), s="not available"; else, s=sprintf('%s, f0=%g, %s, M=%g, dx=%.1f mm (PR AUC %.3f)',T.geometry(1),T.f0(1),T.field_regime(1),T.M(1),1e3*T.dx(1),T.PR_AUC_gt20(1)); end
end

function A=map_from_rows(T,v)
A=nan(max(T.map_iz),max(T.map_ix)); A(sub2ind(size(A),T.map_iz,T.map_ix))=v;
end

function label=short_detector(x)
label=replace(string(x),["Rule_q_gradient","Rule_q_hybrid_minus_local", ...
    "Rule_q_hybrid_minus_theory","Rule_max_local_theory_disagreement", ...
    "ML_logistic_regression","ML_bagged_trees","ML_boosted_trees"], ...
    ["q-gradient rule","hybrid-local rule","hybrid-theory rule", ...
    "max-disagreement rule","logistic ML","bagged trees ML","boosted trees ML"]);
end

function export_fig(fig,file), exportgraphics(fig,file,'Resolution',220); close(fig); end
function s=sanitize(x), s=regexprep(char(x),'[^A-Za-z0-9_-]','_'); end

function T=concat_tables(A,B)
if isempty(A), T=B; return; elseif isempty(B), T=A; return; end
vars=unique([string(A.Properties.VariableNames),string(B.Properties.VariableNames)],'stable');
A=add_missing(A,vars); B=add_missing(B,vars); T=[A(:,cellstr(vars));B(:,cellstr(vars))];
end

function T=add_missing(T,vars)
for i=1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)), continue; end
    T.(vars(i))=nan(height(T),1);
end
end
