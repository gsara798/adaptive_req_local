%% analyze_test_26_confidence_gated_corrections.m
% Test 26: confidence-gated spectral and spatial corrections.
%
% This script never trains or modifies a q model or confidence detector. It
% consumes frozen Test 22 condition checkpoints, reconstructs each local REQ
% radial spectrum once, and compares operational corrections. The frozen
% q_pred is held fixed when a donut changes the q-to-k inversion. True SWS,
% material purity, and interface distance are attached only after inference.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST26_MODE          = quick | full (preferred)
%   ADAPTIVE_REQ_TEST26_QUICK_MODE    = true | false (legacy; default false)
%   ADAPTIVE_REQ_TEST26_VALIDATE_ONLY = true | false (default false)
%   ADAPTIVE_REQ_TEST26_INCLUDE_GLOBAL = true | false (default false)
%   ADAPTIVE_REQ_TEST26_INCLUDE_BOOSTED = true | false (default false)
%
% QUICK_MODE uses f0=500 Hz, M=[2 3], directional_2D and diffuse_3D.

clear; clc; close all;
format compact;

%% Configuration

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();

CORR=struct();
run_mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST26_MODE'))));
if run_mode=="quick"
    CORR.QuickMode=true;
elseif run_mode=="full"
    CORR.QuickMode=false;
elseif run_mode==""
    CORR.QuickMode=env_true('ADAPTIVE_REQ_TEST26_QUICK_MODE',false);
else
    error('ADAPTIVE_REQ_TEST26_MODE must be quick or full.');
end
CORR.RunMode=ternary(CORR.QuickMode,"quick","full");
CORR.ValidateOnly=env_true('ADAPTIVE_REQ_TEST26_VALIDATE_ONLY',false);
CORR.SourceRoot=fullfile(root_dir,'outputs', ...
    'test_22_confidence_external_validation','analysis');
CORR.Models=["HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","LocalOnly_T18"];
CORR.IncludeGlobalBaseline=env_true( ...
    'ADAPTIVE_REQ_TEST26_INCLUDE_GLOBAL',false);
if CORR.IncludeGlobalBaseline, CORR.Models(end+1)="GlobalOnly_T18"; end
CORR.Detectors="ML_bagged_trees";
CORR.IncludeBoostedDetector=env_true( ...
    'ADAPTIVE_REQ_TEST26_INCLUDE_BOOSTED',false);
if CORR.IncludeBoostedDetector, CORR.Detectors(end+1)="ML_boosted_trees"; end
CORR.SaveAllConditionMaps=env_true( ...
    'ADAPTIVE_REQ_TEST26_SAVE_ALL_MAPS',true);
CORR.PrimaryModel="HybridLocalGlobal_T18_noUserRegime";
CORR.PrimaryDetector="ML_bagged_trees";
CORR.PhysicalDonut=struct('c_min',0.5,'c_max',10,'taper_relative',0.05);
CORR.PriorDonutRanges=[0.7 1.3; 0.5 1.5];
CORR.PriorMedianWindow=[3 3];
CORR.ConfidenceThresholds=[0.5 0.7 0.8];
CORR.PeakMinRelativeHeight=0.15;
CORR.PurityEdges=[0 0.70 0.90 0.99 1+eps];
CORR.DistanceEdgesMm=[0 1 2 4 8 Inf];
CORR.RandomSeed=26001;
CORR.Dx=0.2e-3;
CORR.CompactPatchStrategies=["baseline","physical_donut_all", ...
    "prior_donut_all_a070_130_physical", ...
    "confidence_gated_prior_donut_c070_a070_130_physical", ...
    "confidence_gated_interpolation_c080", ...
    "peak_candidate_low_confidence_c080"];
CORR.RepresentativeStrategy= ...
    "confidence_gated_interpolation_c080";
if CORR.QuickMode
    CORR.Frequencies=500;
    CORR.M=[2 3];
    CORR.Regimes=["directional_2D","diffuse_3D"];
    CORR.TargetStepM=0.5e-3;
else
    CORR.Frequencies=[300 400 500 600];
    CORR.M=[2 3 4];
    CORR.Regimes=["directional_2D","diffuse_2D", ...
        "partial_3D","diffuse_3D"];
    CORR.TargetStepM=0.5e-3;
end

OUT=make_output_dirs(root_dir,CORR.QuickMode);
write_config_json(CORR,fullfile(OUT.root_dir,'test26_configuration.json'));

fprintf('\nTest 26: frozen confidence-gated corrections\n');
fprintf('No training. No oracle variables are used for correction.\n');
fprintf('Run mode: %s | Quick mode: %d | source: %s\n', ...
    CORR.RunMode,CORR.QuickMode,CORR.SourceRoot);
fprintf('Target correction-map step: %.2f mm.\n',1e3*CORR.TargetStepM);
fprintf('Spatial resolutions: %s mm.\n', ...
    strjoin(compose('%.1f',1e3*CORR.Dx),', '));

%% Discover frozen Test 22 conditions

condition_dir=fullfile(CORR.SourceRoot,'data','conditions');
field_dir=fullfile(CORR.SourceRoot,'data','field_cache');
assert(exist(condition_dir,'dir')==7, ...
    'Test 22 condition cache not found: %s',condition_dir);
assert(exist(field_dir,'dir')==7, ...
    'Test 22 field cache not found: %s',field_dir);
FILES=discover_conditions(condition_dir,CORR);
assert(~isempty(FILES),'No Test 22 conditions match the Test 26 design.');
fprintf('Matched %d cached Test 22 REQ conditions.\n',numel(FILES));

[available_models,available_detectors]=inspect_availability(FILES(1).path);
cached_models=intersect(CORR.Models,available_models,'stable');
missing_model_names=setdiff(CORR.Models,cached_models,'stable');
MISSING_MODEL_SPECS=load_missing_q_models(root_dir,missing_model_names);
reconstructed_models=string({MISSING_MODEL_SPECS.model_name});
models=[cached_models(:); reconstructed_models(:)];
detectors=intersect(CORR.Detectors,available_detectors,'stable');
report_missing(CORR.Models,models,'frozen q model');
report_missing(CORR.Detectors,detectors,'confidence detector');
assert(ismember(CORR.PrimaryModel,models), ...
    'Primary frozen q model is unavailable in Test 22 cache.');
assert(ismember(CORR.PrimaryDetector,detectors), ...
    'Primary frozen confidence detector is unavailable.');
fprintf('Models evaluated: %s\n',strjoin(models,', '));
fprintf('Detectors evaluated: %s\n',strjoin(detectors,', '));
if ~isempty(MISSING_MODEL_SPECS)
    fprintf(['Missing Test 22 predictions will be reconstructed from frozen ', ...
        'models once per condition (%s).\n'],strjoin(reconstructed_models,', '));
end
if CORR.ValidateOnly
    validate_helper();
    fprintf('Validation-only check passed. No condition was evaluated.\n');
    return;
end

%% Apply corrections condition by condition

PATCH_PARTS=cell(numel(FILES),1);
STAT_PARTS=cell(numel(FILES),1);
for ci=1:numel(FILES)
    condition_tic=tic;
    F=FILES(ci);
    checkpoint=fullfile(OUT.condition_dir, ...
        "corrected__"+sanitize(F.key)+".mat");
    if exist(checkpoint,'file')==2
        checkpoint_info=whos('-file',checkpoint);
        variables=string({checkpoint_info.name});
        if ismember("RESULT",variables)
            S=load(checkpoint,'RESULT'); RESULT=S.RESULT;
            required=["T_compact","T_stats","prediction_cube", ...
                "strategy_specs","models","detectors","step_info"];
            if all(isfield(RESULT,cellstr(required)))
                RESULT.T_compact=refresh_compact_result(RESULT,CORR);
                PATCH_PARTS{ci}=RESULT.T_compact;
                STAT_PARTS{ci}=RESULT.T_stats;
                save(checkpoint,'RESULT','-v7.3');
                fprintf('[%d/%d] Reused %s (%d compact patches).\n', ...
                    ci,numel(FILES),F.key,height(RESULT.T_compact));
                continue;
            end
        end
        warning('Test26:LegacyCheckpoint', ...
            ['Ignoring incompatible checkpoint %s. It was produced by an ', ...
            'older Test 26 layout and will be replaced after this condition.'], ...
            checkpoint);
    end

    S=load(F.path,'T_condition');
    T=S.T_condition; clear S;
    T=T(ismember(T.detector,detectors),:);
    assert(~isempty(T),'No requested frozen predictions in %s.',F.path);
    reference_full=sortrows(T(T.sws_model==CORR.PrimaryModel & ...
        T.detector==CORR.PrimaryDetector,:),{'map_iz','map_ix'});
    assert(height(reference_full)==height(unique(reference_full(:,{'map_iz','map_ix'}),'rows')), ...
        'Duplicate primary patch rows in %s.',F.path);
    [keep_grid,step_info]=spatial_subsample(reference_full,CORR.TargetStepM);
    reference=reference_full(keep_grid,:);
    T=T(ismember([T.map_iz T.map_ix], ...
        [reference.map_iz reference.map_ix],'rows'),:);

    field_file=find_field_file(field_dir,reference(1,:));
    SF=load(field_file,'sim','cfg_sim');
    if ~isempty(MISSING_MODEL_SPECS)
        prediction_file=fullfile(OUT.spectral_dir, ...
            "missing_models__"+sanitize(F.key)+".mat");
        if exist(prediction_file,'file')==2
            PM=load(prediction_file,'T_missing_models');
            T_missing_models=PM.T_missing_models;
            expected_rows=height(reference)*numel(MISSING_MODEL_SPECS)*numel(detectors);
            if height(T_missing_models)~=expected_rows
                warning('Test26:StaleMissingModelCache', ...
                    'Rebuilding stale missing-model cache for %s.',F.key);
                T_missing_models=reconstruct_missing_model_rows(reference,T, ...
                    SF.sim,SF.cfg_sim,MISSING_MODEL_SPECS,detectors, ...
                    step_info.source_step_x_px,step_info.source_step_z_px);
                save(prediction_file,'T_missing_models','-v7.3');
            end
        else
            fprintf('  Reconstructing missing frozen q-model features for %s.\n',F.key);
            T_missing_models=reconstruct_missing_model_rows(reference,T, ...
                SF.sim,SF.cfg_sim,MISSING_MODEL_SPECS,detectors, ...
                step_info.source_step_x_px,step_info.source_step_z_px);
            save(prediction_file,'T_missing_models','-v7.3');
        end
        T=concat_tables(T,T_missing_models);
    end
    T=T(ismember(T.sws_model,models)&ismember(T.detector,detectors),:);
    spectral_file=fullfile(OUT.spectral_dir, ...
        "spectra__"+sanitize(F.key)+".mat");
    if exist(spectral_file,'file')==2
        Q=load(spectral_file,'SPEC'); SPEC=Q.SPEC;
        if ~spectral_cache_matches(SPEC,reference,SF.cfg_sim)
            warning('Test26:StaleSpectralCache', ...
                'Rebuilding stale spectral cache for %s.',F.key);
            SPEC=extract_patch_spectra(reference,SF.sim,SF.cfg_sim);
            save(spectral_file,'SPEC','-v7.3');
        end
    else
        SPEC=extract_patch_spectra(reference,SF.sim,SF.cfg_sim);
        save(spectral_file,'SPEC','-v7.3');
    end

    prediction_cube=[]; STRATEGY_SPECS=table(); T_stats=table();
    for mi=1:numel(models)
        for di=1:numel(detectors)
            X=sortrows(T(T.sws_model==models(mi)&T.detector==detectors(di),:), ...
                {'map_iz','map_ix'});
            if height(X)~=height(reference)
                warning('Test26:IncompleteRows', ...
                    'Skipping %s/%s in %s: incomplete rows.', ...
                    models(mi),detectors(di),F.key);
                continue;
            end
            rows=evaluate_strategy_family(X,SPEC,SF.cfg_sim,CORR,F.key);
            [cube_slice,specs]=pack_strategy_rows(rows,height(reference));
            if isempty(prediction_cube)
                prediction_cube=nan(height(reference),height(specs), ...
                    numel(models),numel(detectors),'single');
                STRATEGY_SPECS=specs;
            else
                assert(isequal(STRATEGY_SPECS.strategy_name,specs.strategy_name), ...
                    'Strategy ordering changed within condition.');
            end
            prediction_cube(:,:,mi,di)=single(cube_slice); %#ok<SAGROW>
            stats=condition_statistics(rows,reference(1,:));
            T_stats=concat_tables(T_stats,stats);
        end
    end
    assert(all(isfinite(prediction_cube(:))&prediction_cube(:)>0), ...
        'Correction produced invalid SWS values.');
    T_meta=make_patch_metadata(reference,SPEC,F.key,step_info);
    T_compact=make_compact_patch_table(T_meta,prediction_cube, ...
        STRATEGY_SPECS,models,detectors,CORR);
    RESULT=struct('prediction_cube',prediction_cube, ...
        'strategy_specs',STRATEGY_SPECS,'models',models, ...
        'detectors',detectors,'T_stats',T_stats,'T_compact',T_compact, ...
        'step_info',step_info);
    save(checkpoint,'RESULT','-v7.3');
    PATCH_PARTS{ci}=T_compact; STAT_PARTS{ci}=T_stats;
    fprintf(['[%d/%d] %s: %d/%d patches at %.2f x %.2f mm -> ', ...
        '%d strategies in %.1f s.\n'],ci,numel(FILES),F.key, ...
        height(reference),height(reference_full),1e3*step_info.actual_step_x_m, ...
        1e3*step_info.actual_step_z_m,height(STRATEGY_SPECS),toc(condition_tic));
    clear T RESULT T_meta T_compact T_stats prediction_cube SPEC SF;
end

T_patch=vertcat(PATCH_PARTS{:}); clear PATCH_PARTS;
T_condition_stats=vertcat(STAT_PARTS{:}); clear STAT_PARTS;
writetable(T_patch,fullfile(OUT.table_dir,'test26_patch_level_results.csv'));
save(fullfile(OUT.data_dir,'test26_compact_results.mat'), ...
    'T_patch','T_condition_stats','CORR','-v7.3');

%% Metrics and strategy selection

T_overall=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector"],"overall");
T_geometry=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector","geometry"],"overall");
T_regime=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector","field_regime"],"overall");
T_frequency_M=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector","f0","M"],"overall");
T_purity=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector","geometry","purity_bin"],"purity");
T_distance=aggregate_condition_stats(T_condition_stats, ...
    ["strategy_name","sws_model","detector","geometry","distance_bin"],"distance");
T_best=best_strategy_candidates_from_summary(T_overall,T_regime,T_frequency_M);

writetable(T_overall,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_overall.csv'));
writetable(T_geometry,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_by_geometry.csv'));
writetable(T_regime,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_by_regime.csv'));
writetable(T_frequency_M,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_by_frequency_M.csv'));
writetable(T_purity,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_by_purity_bin.csv'));
writetable(T_distance,fullfile(OUT.table_dir, ...
    'test26_strategy_summary_by_distance_bin.csv'));
writetable(T_best,fullfile(OUT.table_dir, ...
    'test26_best_strategy_candidates.csv'));

%% Figures

T_plot=compact_to_plot_table(T_patch,CORR,200000);
plot_summary_bars(T_plot,CORR,OUT);
plot_error_by_distance(T_plot,CORR,OUT);
plot_error_by_purity(T_plot,CORR,OUT);
plot_representative_maps(T_plot,CORR,OUT);
plot_baseline_vs_corrected(T_plot,CORR,OUT);
plot_strategy_heatmap(T_plot,CORR,OUT);
if CORR.SaveAllConditionMaps
    plot_all_condition_maps(T_patch,CORR,OUT);
end

%% Interpretation

print_interpretation_summary(T_overall,CORR);
fprintf('\nTables: %s\nFigures: %s\n',OUT.table_dir,OUT.figure_dir);
fprintf('Test 26 complete.\n');

%% Local functions

function OUT=make_output_dirs(root_dir,quick)
OUT.root_dir=fullfile(root_dir,'outputs', ...
    'test_26_confidence_gated_corrections');
if quick, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables');
OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data');
OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
OUT.spectral_dir=fullfile(OUT.data_dir,'spectral_checkpoints');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs)
    if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end
end
end

function tf=env_true(name,default_value)
value=strtrim(getenv(name));
if isempty(value), tf=default_value; return; end
tf=any(strcmpi(value,{'true','1','yes','on'}));
end

function write_config_json(CORR,file)
copy=CORR;
fields=fieldnames(copy);
for i=1:numel(fields)
    if isstring(copy.(fields{i})), copy.(fields{i})=cellstr(copy.(fields{i})); end
end
fid=fopen(file,'w'); assert(fid>=0,'Cannot write %s.',file);
cleaner=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(copy,'PrettyPrint',true),'char');
end

function FILES=discover_conditions(folder,CORR)
d=dir(fullfile(folder,'result__*.mat'));
FILES=struct('path',{},'key',{});
for i=1:numel(d)
    name=string(d(i).name); lower_name=lower(name);
    keep_f=any(arrayfun(@(f)contains(lower_name,"__f"+string(f)+"__"), ...
        CORR.Frequencies));
    keep_m=any(arrayfun(@(m)contains(lower_name,"__m"+string(m)+"__"),CORR.M));
    keep_dx=any(arrayfun(@(dx)contains(lower_name, ...
        "__dx"+string(round(1e6*dx))+"um__"),CORR.Dx));
    keep_r=false;
    for r=CORR.Regimes
        keep_r=keep_r||contains(lower_name,"__"+lower(r)+"__");
    end
    if keep_f&&keep_m&&keep_dx&&keep_r
        FILES(end+1)=struct('path',fullfile(d(i).folder,d(i).name), ...
            'key',erase(name,["result__" ".mat"])); %#ok<AGROW>
    end
end
end

function [models,detectors]=inspect_availability(file)
S=load(file,'T_condition');
models=unique(string(S.T_condition.sws_model),'stable');
detectors=unique(string(S.T_condition.detector),'stable');
end

function report_missing(requested,available,label)
missing=setdiff(requested,available,'stable');
for x=missing
    warning('Test26:FrozenInputUnavailable', ...
        '%s unavailable and will be skipped: %s.',label,x);
end
end

function specs=load_missing_q_models(root_dir,names)
specs=struct('model_name',{},'model_type',{},'model',{});
if isempty(names), return; end
folder=fullfile(root_dir,'outputs','model_registry','test19_clean_field_regime');
requests=[ ...
    struct('report',"LocalOnly_T18",'stored',"LocalOnly_T18", ...
        'feature',"CleanFieldRegime_noUser"); ...
    struct('report',"GlobalOnly_T18",'stored',"GlobalOnly_T18", ...
        'feature',"CleanFieldRegime_noUser")];
for i=1:numel(requests)
    if ~ismember(requests(i).report,names), continue; end
    try
        model=adaptive_req.analysis.load_q_model_deployment(folder, ...
            'ModelName',requests(i).stored,'FeatureSet',requests(i).feature, ...
            'ModelType','bagged_trees');
        predictors=lower(string({model.encoder.entries.name}));
        banned=["q_true","q_theory","cs_true","true_sws", ...
            "distance_to_interface","sws_error","high_error"];
        assert(isempty(intersect(predictors,banned))&& ...
            ~any(contains(predictors,"error")|contains(predictors,"target")), ...
            'Oracle predictor in frozen model %s.',requests(i).report);
        specs(end+1)=struct('model_name',requests(i).report, ...
            'model_type',"bagged_trees",'model',model); %#ok<AGROW>
    catch ME
        warning('Test26:ModelUnavailable','Skipping %s: %s', ...
            requests(i).report,ME.message);
    end
end
end

function T_missing=reconstruct_missing_model_rows(reference,T_cached,sim,cfg, ...
    specs,detectors,step_x,step_z)
M=reference.M(1);
F=extract_frozen_feature_table(sim,cfg,M,step_x,step_z, ...
    reference.field_regime(1),reference.condition_id(1));
F=sortrows(F,{'map_iz','map_ix'});
F=F(ismember([F.map_iz F.map_ix], ...
    [reference.map_iz reference.map_ix],'rows'),:);
assert(height(F)==height(reference)&& ...
    isequal(F.map_iz,reference.map_iz)&&isequal(F.map_ix,reference.map_ix), ...
    'Reconstructed feature grid does not match frozen Test 22 rows.');
T_missing=table();
for si=1:numel(specs)
    Q=adaptive_req.analysis.predict_q_model_from_table(specs(si).model,F, ...
        'ModelType',specs(si).model_type,'ModelName',specs(si).model_name);
    q=min(max(Q.q_pred,0.001),0.999);
    sws=q_to_cs(q,F.req_mapping,cfg.f0);
    for di=1:numel(detectors)
        template=sortrows(T_cached(T_cached.sws_model== ...
            "HybridLocalGlobal_T18_noUserRegime"& ...
            T_cached.detector==detectors(di),:),{'map_iz','map_ix'});
        assert(height(template)==height(F), ...
            'Detector template rows are incomplete for %s.',detectors(di));
        template.sws_model(:)=specs(si).model_name;
        template.q_pred=q; template.SWS_pred=sws;
        template.SWS_error_pct=100*(sws-template.true_SWS)./template.true_SWS;
        template.abs_SWS_error_pct=abs(template.SWS_error_pct);
        template.high_error_gt10=template.abs_SWS_error_pct>10;
        template.high_error_gt20=template.abs_SWS_error_pct>20;
        T_missing=concat_tables(T_missing,template);
    end
end
end

function [keep,info]=spatial_subsample(T,target_step_m)
x=unique(T.x,'sorted'); z=unique(T.z,'sorted');
assert(numel(x)>=2&&numel(z)>=2,'Cannot infer Test 22 map stride.');
source_x=median(diff(x)); source_z=median(diff(z));
factor_x=max(1,ceil(target_step_m/source_x));
factor_z=max(1,ceil(target_step_m/source_z));
keep=mod(T.map_ix-min(T.map_ix),factor_x)==0& ...
    mod(T.map_iz-min(T.map_iz),factor_z)==0;
info=struct('source_step_x_m',source_x,'source_step_z_m',source_z, ...
    'source_step_x_px',max(1,round(source_x/T.dx(1))), ...
    'source_step_z_px',max(1,round(source_z/T.dz(1))), ...
    'subsample_factor_x',factor_x,'subsample_factor_z',factor_z, ...
    'actual_step_x_m',factor_x*source_x, ...
    'actual_step_z_m',factor_z*source_z);
end

function T=extract_frozen_feature_table(sim,cfg,M,step_x,step_z,regime,id)
feat=adaptive_req.config.default_feature_config('M',M,'cs_guess',3, ...
    'gamma_win',1,'pad_factor',1);
req_args={'Nbins','auto','Nbins_auto_oversample',1,'Nbins_min',16, ...
    'smooth_sigma',1};
[req_cfg,feat]=adaptive_req.config.default_req_config(cfg,feat,req_args{:});
cfg_req=cfg; cfg_req.UseParfor=false;
[q_global,curve,global_features]= ...
    adaptive_req.quantile.compute_global_quantile_from_field( ...
    sim.Uxz,cfg_req,req_cfg,feat);
global_shape=adaptive_req.quantile.extract_ecum_shape_features(curve);
O=adaptive_req.estimators.req_estimator_map(sim.Uxz,cfg_req,feat, ...
    'StepX',step_x,'StepZ',step_z,'EdgeMode','valid', ...
    'QuantileMode','local_req','ReqOptions',req_args, ...
    'ReturnFeatures',false,'ReturnFeatureTable',true, ...
    'ReuseReqSpectrumForFeatures',true,'UseWindowParfor',true, ...
    'StoreReqCurves',false,'Verbose',false);
T=O.feature_table; n=height(T);
T.global_req_mapping=repmat( ...
    {adaptive_req.quantile.make_req_mapping(curve)},n,1);
T.q_global_req=q_global*ones(n,1);
T.global_REQ_Nbins_effective=curve.Nbins_effective*ones(n,1);
T=assign_prefixed(T,global_features,"global_");
T=assign_prefixed(T,global_shape,"global_");
T.condition_id=id*ones(n,1);
T.field_regime_label=repmat(string(regime),n,1);
T.user_field_guess=map_user_guess(T.field_regime_label);
T.SIM_f0=cfg.f0*ones(n,1); T.SIM_dx=cfg.dx*ones(n,1);
T.SIM_dz=cfg.dz*ones(n,1); T.SIM_cs_bg=cfg.cs_bg*ones(n,1);
T.SIM_cs_inc=cfg.cs_inc*ones(n,1); T.SIM_Nwaves=cfg.Nwaves*ones(n,1);
T.SIM_Is2D=repmat(logical(cfg.Is2D),n,1);
T.SIM_ForceInPlaneWave=repmat(logical(cfg.ForceInPlaneWave),n,1);
T.REQ_M=M*ones(n,1); T.REQ_cs_guess=3*ones(n,1);
T.M_eff_guess=M*ones(n,1); T.step_idx=M*ones(n,1);
T.realization_idx=ones(n,1);
T=adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T);
T=adaptive_req.analysis.Test12Analysis.addBins(T);
end

function T=assign_prefixed(T,S,prefix)
names=fieldnames(S);
for i=1:numel(names)
    value=S.(names{i});
    if isnumeric(value)&&isscalar(value)
        T.(char(prefix+string(names{i})))=repmat(double(value),height(T),1);
    end
end
end

function labels=map_user_guess(regime)
labels=strings(size(regime));
labels(regime=="directional_2D")="directional_like";
labels(regime=="partial_3D")="partially_diffuse";
labels(regime=="diffuse_2D"|regime=="diffuse_3D")="diffuse_like";
labels(labels=="")="unknown";
labels=categorical(labels,["directional_like","partially_diffuse", ...
    "diffuse_like","unknown"]);
end

function sws=q_to_cs(q,mappings,f0)
sws=nan(size(q));
for i=1:numel(q)
    if isfinite(q(i))
        sws(i)=adaptive_req.quantile.quantile_to_cs(mappings{i},q(i),f0);
    end
end
end

function file=find_field_file(folder,row)
key=sprintf('%s__f%g__%s__dx%gum',row.geometry,row.f0, ...
    lower(string(row.field_regime)),round(1e6*row.dx));
file=fullfile(folder,"field__"+sanitize(key)+".mat");
assert(exist(file,'file')==2,'Field cache missing: %s',file);
end

function SPEC=extract_patch_spectra(T,sim,cfg)
M=T.M(1);
feat=adaptive_req.config.default_feature_config('M',M,'cs_guess',3, ...
    'gamma_win',1,'pad_factor',1);
[req,feat]=adaptive_req.config.default_req_config(cfg,feat, ...
    'Nbins','auto','Nbins_auto_oversample',1,'Nbins_min',16, ...
    'smooth_sigma',1);
half=floor(feat.win_size/2); n=height(T);
curves=cell(n,1); purity=nan(n,1); soft_fraction=nan(n,1);
hard_fraction=nan(n,1); cx=round(T.x/cfg.dx)+1; cz=round(T.z/cfg.dz)+1;
cfg_patch=cfg; cfg_patch.UseParfor=false;
Uxz=sim.Uxz; cs_map=sim.cs_map;
parfor i=1:n
    zi=(cz(i)-half):(cz(i)+half); xi=(cx(i)-half):(cx(i)+half);
    patch=Uxz(zi,xi); %#ok<PFBNS>
    [~,curve]=adaptive_req.quantile.compute_quantile_from_patch( ...
        patch,cfg_patch,req);
    curves{i}=struct('k_cent',curve.k_cent,'Srad',curve.Srad, ...
        'Ecum',curve.Ecum);
    material=cs_map(zi,xi); center=cs_map(cz(i),cx(i)); %#ok<PFBNS>
    soft_fraction(i)=mean(material(:)<2.5);
    hard_fraction(i)=mean(material(:)>=2.5);
    purity(i)=mean(abs(material(:)-center)<0.25);
end
SPEC=struct('curves',{curves},'patch_purity',purity, ...
    'soft_fraction',soft_fraction,'hard_fraction',hard_fraction, ...
    'cx',cx,'cz',cz,'win_size',feat.win_size);
end

function tf=spectral_cache_matches(SPEC,T,cfg)
required={'curves','patch_purity','soft_fraction','hard_fraction','cx','cz'};
tf=all(isfield(SPEC,required));
if ~tf, return; end
n=height(T);
tf=numel(SPEC.curves)==n&&numel(SPEC.patch_purity)==n&& ...
    numel(SPEC.soft_fraction)==n&&numel(SPEC.hard_fraction)==n;
if ~tf, return; end
expected_cx=round(T.x/cfg.dx)+1;
expected_cz=round(T.z/cfg.dz)+1;
tf=isequal(double(SPEC.cx(:)),double(expected_cx(:)))&& ...
    isequal(double(SPEC.cz(:)),double(expected_cz(:)));
end

function T=evaluate_strategy_family(X,SPEC,~,CORR,condition_key)
q=X.q_pred; baseline=X.SWS_pred;
[baseline_map,compact_iz,compact_ix]=map_from_rows(X,baseline);
prior_map=medfilt2(baseline_map,CORR.PriorMedianWindow,'symmetric');
lin=sub2ind(size(prior_map),compact_iz,compact_ix);
prior_sws=prior_map(lin);
prior_k=2*pi*X.f0./prior_sws;

physical=correct_with_donut(SPEC.curves,q,X.f0,prior_k,CORR, ...
    true,false,CORR.PriorDonutRanges(1,:),baseline);
prior=cell(size(CORR.PriorDonutRanges,1),2);
for ai=1:size(CORR.PriorDonutRanges,1)
    range=CORR.PriorDonutRanges(ai,:);
    prior{ai,1}=correct_with_donut(SPEC.curves,q,X.f0,prior_k,CORR, ...
        false,true,range,baseline);
    prior{ai,2}=correct_with_donut(SPEC.curves,q,X.f0,prior_k,CORR, ...
        true,true,range,baseline);
end

T=table();
T=append_strategy(T,make_rows(X,SPEC,baseline,"baseline",false, ...
    NaN,NaN,NaN,condition_key));
T=append_strategy(T,make_rows(X,SPEC,physical,"physical_donut_all",true, ...
    NaN,NaN,NaN,condition_key));
for ai=1:size(CORR.PriorDonutRanges,1)
    range=CORR.PriorDonutRanges(ai,:); tag=alpha_tag(range);
    T=append_strategy(T,make_rows(X,SPEC,prior{ai,1}, ...
        "prior_donut_all_"+tag,false,range(1),range(2),NaN,condition_key));
    T=append_strategy(T,make_rows(X,SPEC,prior{ai,2}, ...
        "prior_donut_all_"+tag+"_physical",true, ...
        range(1),range(2),NaN,condition_key));
end

for ti=1:numel(CORR.ConfidenceThresholds)
    threshold=CORR.ConfidenceThresholds(ti);
    low=X.confidence<threshold; ctag="c"+sprintf('%03d',round(100*threshold));
    for ai=1:size(CORR.PriorDonutRanges,1)
        range=CORR.PriorDonutRanges(ai,:); tag=alpha_tag(range);
        for physical_flag=0:1
            corrected=baseline;
            candidate=prior{ai,physical_flag+1};
            corrected(low)=candidate(low);
            name="confidence_gated_prior_donut_"+ctag+"_"+tag;
            if physical_flag, name=name+"_physical"; end
            T=append_strategy(T,make_rows(X,SPEC,corrected,name, ...
                logical(physical_flag),range(1),range(2),threshold,condition_key));
        end
    end

    for physical_flag=0:1
        source=baseline;
        if physical_flag, source=physical; end
        corrected=interpolate_low_confidence(X,source,low);
        name="confidence_gated_interpolation_"+ctag;
        if physical_flag, name=name+"_physical"; end
        T=append_strategy(T,make_rows(X,SPEC,corrected,name, ...
            logical(physical_flag),NaN,NaN,threshold,condition_key));

        peak=peak_candidate_correction(SPEC.curves,source,prior_k,X.f0, ...
            low,CORR,logical(physical_flag));
        name="peak_candidate_low_confidence_"+ctag;
        if physical_flag, name=name+"_physical"; end
        T=append_strategy(T,make_rows(X,SPEC,peak,name, ...
            logical(physical_flag),NaN,NaN,threshold,condition_key));
    end
end
end

function sws=correct_with_donut(curves,q,f0,prior_k,CORR, ...
    physical,apply_prior,range,fallback)
n=numel(curves); sws=nan(n,1);
for i=1:n
    try
        mapping=adaptive_req.analysis.apply_radial_donut(curves{i}, ...
            'F0',f0(i),'ApplyPhysical',physical, ...
            'PhysicalSpeedRange',[CORR.PhysicalDonut.c_min ...
            CORR.PhysicalDonut.c_max], ...
            'ApplyPrior',apply_prior,'PriorK',prior_k(i), ...
            'RelativeRange',range,'TaperRelative', ...
            CORR.PhysicalDonut.taper_relative);
        sws(i)=adaptive_req.quantile.quantile_to_cs(mapping,q(i),f0(i));
    catch
        sws(i)=fallback(i);
    end
end
bad=~isfinite(sws)|sws<=0; sws(bad)=fallback(bad);
end

function corrected=interpolate_low_confidence(X,source,low)
corrected=source; good=~low&isfinite(source);
if nnz(good)<3, return; end
try
    interpolant=scatteredInterpolant(X.x(good),X.z(good),source(good), ...
        'natural','nearest');
    values=interpolant(X.x(low),X.z(low));
    ok=isfinite(values); idx=find(low); corrected(idx(ok))=values(ok);
catch
    % Preserve source if the reliable-pixel geometry is degenerate.
end
end

function corrected=peak_candidate_correction(curves,source,prior_k,f0, ...
    low,CORR,physical)
corrected=source; idx=find(low);
for jj=1:numel(idx)
    i=idx(jj); curve=curves{i};
    if physical
        curve=adaptive_req.analysis.apply_radial_donut(curve,'F0',f0(i), ...
            'ApplyPhysical',true,'PhysicalSpeedRange', ...
            [CORR.PhysicalDonut.c_min CORR.PhysicalDonut.c_max], ...
            'TaperRelative',CORR.PhysicalDonut.taper_relative);
    end
    k=double(curve.k_cent(:)); s=double(curve.Srad(:));
    loc=find(s(2:end-1)>=s(1:end-2)&s(2:end-1)>s(3:end))+1;
    loc=loc(s(loc)>=CORR.PeakMinRelativeHeight*max(s));
    k_base=2*pi*f0(i)/source(i); candidates=[k_base;k(loc)];
    candidates=candidates(isfinite(candidates)&candidates>0);
    if isempty(candidates), continue; end
    [~,best]=min(abs(log(candidates/prior_k(i))));
    corrected(i)=2*pi*f0(i)/candidates(best);
end
end

function R=make_rows(X,SPEC,pred,name,physical,alo,ahi,cth,key)
R=table(); n=height(X);
R.condition_key=repmat(string(key),n,1);
copy_names={'condition_id','geometry','geometry_type','f0','field_regime', ...
    'M','dx','dz','map_iz','map_ix','x','z','sws_model','detector', ...
    'q_pred','true_SWS','risk','confidence','distance_to_interface_mm'};
for i=1:numel(copy_names), R.(copy_names{i})=X.(copy_names{i}); end
R.strategy_name=repmat(string(name),n,1);
R.physical_donut=repmat(logical(physical),n,1);
R.prior_alpha_low=repmat(alo,n,1); R.prior_alpha_high=repmat(ahi,n,1);
R.confidence_threshold=repmat(cth,n,1);
R.SWS_pred_corrected=pred;
R.SWS_error_signed_pct=100*(pred-R.true_SWS)./R.true_SWS;
R.SWS_error_abs_pct=abs(R.SWS_error_signed_pct);
R.high_error_gt10=R.SWS_error_abs_pct>10;
R.high_error_gt20=R.SWS_error_abs_pct>20;
R.is_underestimate=R.SWS_error_signed_pct<0;
R.is_overestimate=R.SWS_error_signed_pct>0;
R.patch_purity=SPEC.patch_purity;
R.patch_soft_fraction=SPEC.soft_fraction;
R.patch_hard_fraction=SPEC.hard_fraction;
R.patch_is_mixed=R.patch_soft_fraction>0.05&R.patch_hard_fraction>0.05;
R.material_side=repmat("soft",n,1);
R.material_side(R.true_SWS>=2.5)="hard";
R.material_side(R.geometry_type=="homogeneous")="homogeneous";
R.purity_bin=purity_labels(R.patch_purity,R.geometry_type);
R.distance_bin=distance_labels(R.distance_to_interface_mm,R.geometry_type);
R.low_confidence=false(n,1);
if isfinite(cth), R.low_confidence=R.confidence<cth; end
end

function labels=purity_labels(p,geometry_type)
labels=repmat("strongly_mixed",size(p));
labels(p>=0.70)="moderately_mixed";
labels(p>=0.90)="near_pure";
labels(p>=0.99)="pure";
labels(geometry_type=="homogeneous")="homogeneous";
end

function labels=distance_labels(d,geometry_type)
labels=repmat("gt_8mm",size(d));
labels(d<8)="4_8mm"; labels(d<4)="2_4mm";
labels(d<2)="1_2mm"; labels(d<1)="0_1mm";
labels(geometry_type=="homogeneous")="homogeneous";
end

function tag=alpha_tag(range)
tag="a"+sprintf('%03d',round(100*range(1)))+"_"+ ...
    sprintf('%03d',round(100*range(2)));
end

function T=append_strategy(T,R)
if isempty(T), T=R; else, T=[T;R]; end
end

function [cube,specs]=pack_strategy_rows(rows,n)
assert(mod(height(rows),n)==0,'Strategy rows are not patch-aligned.');
indices=1:n:height(rows);
spec_names={'strategy_name','physical_donut','prior_alpha_low', ...
    'prior_alpha_high','confidence_threshold'};
specs=rows(indices,spec_names);
cube=reshape(rows.SWS_pred_corrected,n,[]);
end

function T=condition_statistics(rows,reference)
base=["strategy_name","sws_model","detector"];
T=raw_stats(rows,base,"overall");
T=concat_tables(T,raw_stats(rows,[base "purity_bin"],"purity"));
T=concat_tables(T,raw_stats(rows,[base "distance_bin"],"distance"));
metadata={'geometry','field_regime','f0','M','dx','dz'};
for i=1:numel(metadata)
    T.(metadata{i})=repmat(reference.(metadata{i})(1),height(T),1);
end
end

function S=raw_stats(T,groups,scope)
G=findgroups(T(:,cellstr(groups))); ids=unique(G,'stable'); S=table();
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups));
    row.stat_scope=string(scope); row.N=height(X);
    row.sum_abs_error=sum(X.SWS_error_abs_pct,'omitnan');
    row.sum_signed_error=sum(X.SWS_error_signed_pct,'omitnan');
    row.condition_median_signed_error=median(X.SWS_error_signed_pct,'omitnan');
    row.underestimate_count=sum(X.is_underestimate);
    row.overestimate_count=sum(X.is_overestimate);
    row.high_error_10_count=sum(X.high_error_gt10);
    row.high_error_20_count=sum(X.high_error_gt20);
    descriptors=["physical_donut","prior_alpha_low", ...
        "prior_alpha_high","confidence_threshold"];
    for descriptor=descriptors, row.(descriptor)=X.(descriptor)(1); end
    classes=["homogeneous","pure","near_pure", ...
        "moderately_mixed","strongly_mixed"];
    for class=classes
        idx=X.purity_bin==class; token=char(class);
        row.("N_"+token)=sum(idx);
        row.("sum_abs_"+token)=sum(X.SWS_error_abs_pct(idx),'omitnan');
    end
    idx=X.geometry_type~="homogeneous"&X.patch_purity>=0.90;
    row.N_heterogeneous_nonmixed=sum(idx);
    row.sum_abs_heterogeneous_nonmixed= ...
        sum(X.SWS_error_abs_pct(idx),'omitnan');
    for side=["soft","hard"]
        idx=X.material_side==side; token=char(side);
        row.("N_"+token)=sum(idx);
        row.("sum_abs_"+token)=sum(X.SWS_error_abs_pct(idx),'omitnan');
    end
    S=concat_tables(S,row);
end
end

function T=make_patch_metadata(X,SPEC,key,step_info)
names={'condition_id','geometry','geometry_type','f0','field_regime', ...
    'M','dx','dz','map_iz','map_ix','x','z','q_pred','true_SWS', ...
    'risk','confidence','distance_to_interface_mm'};
T=X(:,names); T.condition_key=repmat(string(key),height(T),1);
T.patch_purity=SPEC.patch_purity; T.patch_soft_fraction=SPEC.soft_fraction;
T.patch_hard_fraction=SPEC.hard_fraction;
T.patch_is_mixed=T.patch_soft_fraction>0.05&T.patch_hard_fraction>0.05;
T.material_side=repmat("soft",height(T),1);
T.material_side(T.true_SWS>=2.5)="hard";
T.material_side(T.geometry_type=="homogeneous")="homogeneous";
T.purity_bin=purity_labels(T.patch_purity,T.geometry_type);
T.distance_bin=distance_labels(T.distance_to_interface_mm,T.geometry_type);
T.correction_step_x_mm=1e3*step_info.actual_step_x_m*ones(height(T),1);
T.correction_step_z_mm=1e3*step_info.actual_step_z_m*ones(height(T),1);
end

function T=make_compact_patch_table(meta,cube,specs,models,detectors,CORR)
T=meta; T.sws_model=repmat(CORR.PrimaryModel,height(T),1);
T.detector=repmat(CORR.PrimaryDetector,height(T),1);
mi=find(models==CORR.PrimaryModel,1); di=find(detectors==CORR.PrimaryDetector,1);
for name=CORR.CompactPatchStrategies
    si=find(specs.strategy_name==name,1);
    if isempty(si), continue; end
    token=matlab.lang.makeValidName("sws_"+name);
    pred=double(cube(:,si,mi,di)); T.(token)=pred;
end
end

function T=refresh_compact_result(RESULT,CORR)
old=RESULT.T_compact; names=string(old.Properties.VariableNames);
prediction_vars=startsWith(names,"sws_")&names~="sws_model";
meta=old(:,~prediction_vars);
meta(:,intersect({'sws_model','detector'},meta.Properties.VariableNames))=[];
T=make_compact_patch_table(meta,RESULT.prediction_cube, ...
    RESULT.strategy_specs,string(RESULT.models),string(RESULT.detectors),CORR);
end

function S=aggregate_condition_stats(T,groups,scope)
T=T(T.stat_scope==scope,:); G=findgroups(T(:,cellstr(groups)));
ids=unique(G,'stable'); S=table();
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups));
    row.N=sum(X.N); row.MAPE=sum(X.sum_abs_error)/max(row.N,1);
    row.mean_signed_error=sum(X.sum_signed_error)/max(row.N,1);
    row.median_signed_error=weighted_median( ...
        X.condition_median_signed_error,X.N);
    row.underestimate_pct=100*sum(X.underestimate_count)/max(row.N,1);
    row.overestimate_pct=100*sum(X.overestimate_count)/max(row.N,1);
    row.high_error_10_pct=100*sum(X.high_error_10_count)/max(row.N,1);
    row.high_error_20_pct=100*sum(X.high_error_20_count)/max(row.N,1);
    descriptors=["physical_donut","prior_alpha_low", ...
        "prior_alpha_high","confidence_threshold"];
    for descriptor=descriptors, row.(descriptor)=X.(descriptor)(1); end
    classes=["homogeneous","pure","near_pure", ...
        "moderately_mixed","strongly_mixed","heterogeneous_nonmixed", ...
        "soft","hard"];
    for class=classes
        token=char(class); count=sum(X.("N_"+token));
        row.("MAPE_"+token)=sum(X.("sum_abs_"+token))/max(count,1);
        if count==0, row.("MAPE_"+token)=NaN; end
    end
    S=concat_tables(S,row);
end
end

function value=weighted_median(x,w)
ok=isfinite(x)&isfinite(w)&w>0; x=x(ok); w=w(ok);
if isempty(x), value=NaN; return; end
[x,order]=sort(x); w=w(order); value=x(find(cumsum(w)>=0.5*sum(w),1));
end

function B=best_strategy_candidates_from_summary(overall,regime,frequency_M)
B=table();
B=concat_tables(B,best_from_one_summary(overall,strings(0,1),"overall"));
B=concat_tables(B,best_from_one_summary(regime,"field_regime","by_regime"));
B=concat_tables(B,best_from_one_summary(frequency_M,["f0","M"],"by_frequency_M"));
end

function B=best_from_one_summary(S,extra_groups,scope)
groups=["sws_model","detector" string(extra_groups)];
G=findgroups(S(:,cellstr(groups))); ids=unique(G,'stable'); B=table();
for i=1:numel(ids)
    X=S(G==ids(i),:); baseline=X(X.strategy_name=="baseline",:);
    if isempty(baseline), continue; end
    [~,im]=min(X.MAPE); [~,ih]=min(X.high_error_20_pct);
    hom_delta=X.MAPE_homogeneous-baseline.MAPE_homogeneous;
    pure_delta=X.MAPE_pure-baseline.MAPE_pure;
    mixed=mean([X.MAPE_moderately_mixed X.MAPE_strongly_mixed],2,'omitnan');
    base_mixed=mean([baseline.MAPE_moderately_mixed ...
        baseline.MAPE_strongly_mixed],2,'omitnan');
    tradeoff=(mixed-base_mixed)+2*max(hom_delta-0.75,0)+max(pure_delta,0);
    [~,it]=min(tradeoff); [~,id]=min(abs(hom_delta));
    row=X(1,cellstr(groups)); row.candidate_scope=string(scope);
    row.lowest_MAPE_strategy=X.strategy_name(im);
    row.lowest_high_error20_strategy=X.strategy_name(ih);
    row.lowest_homogeneous_degradation_strategy=X.strategy_name(id);
    row.best_tradeoff_strategy=X.strategy_name(it);
    row.best_tradeoff_score=tradeoff(it);
    B=concat_tables(B,row);
end
end

function T=compact_to_plot_table(C,CORR,max_patches)
keys=unique(C.condition_key(C.geometry_type~="homogeneous"),'stable');
must_keep=false(height(C),1);
if ~isempty(keys), must_keep=C.condition_key==keys(1); end
remaining=find(~must_keep); rng(CORR.RandomSeed);
n=max(0,min(numel(remaining),max_patches-sum(must_keep)));
if n>0, must_keep(remaining(randperm(numel(remaining),n)))=true; end
C=C(must_keep,:); T=table();
for name=CORR.CompactPatchStrategies
    token=matlab.lang.makeValidName("sws_"+name);
    if ~ismember(token,string(C.Properties.VariableNames)), continue; end
    R=C; R.strategy_name=repmat(name,height(R),1);
    R.SWS_pred_corrected=R.(token);
    R.SWS_error_signed_pct=100*(R.SWS_pred_corrected-R.true_SWS)./R.true_SWS;
    R.SWS_error_abs_pct=abs(R.SWS_error_signed_pct);
    R.high_error_gt20=R.SWS_error_abs_pct>20;
    R.low_confidence=false(height(R),1);
    if contains(name,"c080"), R.low_confidence=R.confidence<0.8;
    elseif contains(name,"c070"), R.low_confidence=R.confidence<0.7;
    elseif contains(name,"c050"), R.low_confidence=R.confidence<0.5;
    end
    variable_names=string(R.Properties.VariableNames);
    prediction_vars=startsWith(variable_names,"sws_")&variable_names~="sws_model";
    R(:,prediction_vars)=[];
    T=concat_tables(T,R);
end
end

function plot_summary_bars(T,CORR,OUT)
X=primary_rows(T,CORR); X.summary_region=X.purity_bin;
X.summary_region(ismember(X.purity_bin,["moderately_mixed","strongly_mixed"]))="mixed";
X.summary_region(ismember(X.purity_bin,["pure","near_pure"]))="pure_heterogeneous";
strategies=plot_strategy_subset(X,CORR); regions=["homogeneous","pure_heterogeneous","mixed"];
M=nan(numel(strategies),numel(regions)); H=M;
for i=1:numel(strategies)
    for j=1:numel(regions)
        idx=X.strategy_name==strategies(i)&X.summary_region==regions(j);
        M(i,j)=mean(X.SWS_error_abs_pct(idx),'omitnan');
        H(i,j)=100*mean(X.high_error_gt20(idx),'omitnan');
    end
end
fig=figure('Color','w','Position',[100 100 1250 470]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors=[0.25 0.45 0.75;0.35 0.70 0.45;0.90 0.50 0.20];
ax1=nexttile(tl); b=bar(ax1,M,'grouped'); set_bar_colors(b,colors);
title(ax1,'Mean absolute error'); ylabel(ax1,'MAPE (%)'); style_axes(ax1);
ax2=nexttile(tl); b=bar(ax2,H,'grouped'); set_bar_colors(b,colors);
title(ax2,'Large-error rate'); ylabel(ax2,'Pixels with error >20% (%)'); style_axes(ax2);
labels=short_strategy(strategies);
set([ax1 ax2],'XTick',1:numel(strategies),'XTickLabel',labels, ...
    'XTickLabelRotation',18,'TickLabelInterpreter','none');
legend(ax2,["Homogeneous","Pure / near-pure","Mixed"], ...
    'Location','northoutside','Orientation','horizontal','Box','off');
export_fig(fig,OUT,'test26_strategy_bar_summary.png');
end

function plot_error_by_distance(T,CORR,OUT)
X=primary_rows(T,CORR); X=X(X.geometry_type~="homogeneous",:);
strategies=plot_strategy_subset(X,CORR); bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
fig=figure('Color','w','Position',[100 100 1180 500]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors=strategy_colors(numel(strategies)); geometry=["bilayer_2_3","circular_inclusion_2_3"];
for gi=1:2
    ax=nexttile(tl); hold(ax,'on');
    for si=1:numel(strategies)
        y=nan(size(bins));
        for bi=1:numel(bins)
            idx=X.geometry==geometry(gi)&X.strategy_name==strategies(si)&X.distance_bin==bins(bi);
            y(bi)=mean(X.SWS_error_abs_pct(idx),'omitnan');
        end
        plot(ax,1:numel(bins),y,'-o','Color',colors(si,:), ...
            'LineWidth',1.8,'MarkerSize',5,'MarkerFaceColor','w', ...
            'DisplayName',short_strategy(strategies(si)));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',["0–1","1–2","2–4","4–8",">8"], ...
        'TickLabelInterpreter','none');
    title(ax,ternary(gi==1,'Bilayer 2/3 m/s','Circular inclusion 2/3 m/s'));
    ylabel(ax,'MAPE (%)'); xlabel(ax,'Distance to interface (mm)'); style_axes(ax);
end
lg=legend(nexttile(tl,2),'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',3,'Box','off','Interpreter','none'); lg.Layout.Tile='south';
export_fig(fig,OUT,'test26_error_vs_distance.png');
end

function plot_error_by_purity(T,CORR,OUT)
X=primary_rows(T,CORR); X=X(X.geometry_type~="homogeneous",:);
strategies=plot_strategy_subset(X,CORR); colors=strategy_colors(numel(strategies));
bins=["strongly_mixed","moderately_mixed","near_pure","pure"];
fig=figure('Color','w','Position',[100 100 1180 500]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
metrics=["SWS_error_abs_pct","high_error_gt20"];
for pi=1:2, ax=nexttile(tl); hold(ax,'on');
    for si=1:numel(strategies)
        y=nan(size(bins));
        for bi=1:numel(bins)
            idx=X.strategy_name==strategies(si)&X.purity_bin==bins(bi);
            y(bi)=mean(X.(metrics(pi))(idx),'omitnan');
        end
        if pi==2, y=100*y; end
        plot(ax,1:numel(bins),y,'-o','Color',colors(si,:), ...
            'LineWidth',1.8,'MarkerSize',5,'MarkerFaceColor','w', ...
            'DisplayName',short_strategy(strategies(si)));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel', ...
        ["Strongly mixed","Moderately mixed","Near-pure","Pure"], ...
        'XTickLabelRotation',15,'TickLabelInterpreter','none');
    title(ax,ternary(pi==1,'Mean absolute error','Large-error rate'));
    ylabel(ax,ternary(pi==1,'MAPE (%)','Pixels with error >20% (%)')); style_axes(ax);
end
lg=legend(nexttile(tl,2),'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',3,'Box','off','Interpreter','none'); lg.Layout.Tile='south';
export_fig(fig,OUT,'test26_error_vs_patch_purity.png');
end

function plot_representative_maps(T,CORR,OUT)
X=primary_rows(T,CORR); keys=unique(X.condition_key(X.geometry_type~="homogeneous"),'stable');
if isempty(keys), return; end
K=X(X.condition_key==keys(1),:); baseline=K(K.strategy_name=="baseline",:);
corrected=K(K.strategy_name==CORR.RepresentativeStrategy,:);
if isempty(corrected)
    candidates=plot_strategy_subset(K,CORR);
    corrected=K(K.strategy_name==candidates(end),:);
end
baseline=sortrows(baseline,{'map_iz','map_ix'}); corrected=sortrows(corrected,{'map_iz','map_ix'});
maps={map_from_rows(baseline,baseline.true_SWS),map_from_rows(baseline,baseline.SWS_pred_corrected), ...
    map_from_rows(corrected,corrected.SWS_pred_corrected), ...
    map_from_rows(baseline,baseline.SWS_error_abs_pct), ...
    map_from_rows(corrected,corrected.SWS_error_abs_pct), ...
    map_from_rows(baseline,baseline.confidence), ...
    map_from_rows(corrected,corrected.low_confidence), ...
    map_from_rows(corrected,corrected.SWS_pred_corrected-baseline.SWS_pred_corrected)};
titles=["True SWS","Baseline SWS","Corrected SWS","Baseline abs error (%)", ...
    "Corrected abs error (%)","Confidence","Low-confidence mask","Corrected - baseline"];
fig=figure('Color','w','Position',[40 60 1420 690]);
tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
error_max=min(60,max([prctile(maps{4}(:),98) prctile(maps{5}(:),98)]));
delta_max=max(abs(maps{8}(:)),[],'omitnan'); delta_max=max(delta_max,0.1);
for i=1:8
    ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image');
    set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',11,'FontWeight','normal');
    if i<=3, clim(ax,[1.5 3.5]); colormap(ax,turbo); colorbar(ax);
    elseif i<=5, clim(ax,[0 error_max]); colormap(ax,turbo); colorbar(ax);
    elseif i==6, clim(ax,[0 1]); colormap(ax,parula); colorbar(ax);
    elseif i==7, clim(ax,[0 1]); colormap(ax,gray); colorbar(ax,'Ticks',[0 1]);
    else, clim(ax,[-delta_max delta_max]); colormap(ax,bluewhitered_map()); colorbar(ax);
    end
end
title(tl,"Representative bilayer — "+short_strategy(CORR.RepresentativeStrategy), ...
    'Interpreter','none','FontSize',13,'FontWeight','bold');
export_fig(fig,OUT,'test26_representative_sws_maps.png');
end

function plot_all_condition_maps(T,CORR,OUT)
baseline_var=matlab.lang.makeValidName("sws_baseline");
corrected_var=matlab.lang.makeValidName("sws_"+CORR.RepresentativeStrategy);
required=[string(baseline_var) string(corrected_var)];
if ~all(ismember(required,string(T.Properties.VariableNames)))
    warning('Test26:ConditionMapsUnavailable', ...
        'Compact patch table lacks baseline or representative correction.');
    return;
end
dx_token="dx"+strjoin(string(round(1e6*CORR.Dx)),"_")+"um";
root=fullfile(OUT.figure_dir,"maps_by_condition_"+dx_token);
keys=unique(T.condition_key,'stable');
fprintf('Saving %d condition map panels under %s.\n',numel(keys),root);
for ki=1:numel(keys)
    C=sortrows(T(T.condition_key==keys(ki),:),{'map_iz','map_ix'});
    folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end
    file=fullfile(folder,"test26_maps__"+sanitize(keys(ki))+".png");
    if exist(file,'file')==2, continue; end

    baseline=C.(baseline_var); corrected=C.(corrected_var);
    signed_base=100*(baseline-C.true_SWS)./C.true_SWS;
    signed_corrected=100*(corrected-C.true_SWS)./C.true_SWS;
    low=C.confidence<0.8;
    maps={map_from_rows(C,C.true_SWS),map_from_rows(C,baseline), ...
        map_from_rows(C,corrected),map_from_rows(C,abs(signed_base)), ...
        map_from_rows(C,abs(signed_corrected)),map_from_rows(C,C.confidence), ...
        map_from_rows(C,low),map_from_rows(C,corrected-baseline)};
    titles=["True SWS","Baseline SWS","Corrected SWS", ...
        "Baseline abs error (%)","Corrected abs error (%)", ...
        "Confidence","Confidence < 0.8","Corrected - baseline"];
    x_mm=1e3*unique(C.x,'sorted'); z_mm=1e3*unique(C.z,'sorted');
    error_max=min(60,max([prctile(maps{4}(:),98) prctile(maps{5}(:),98)]));
    error_max=max(error_max,10);
    delta_max=max(abs(maps{8}(:)),[],'omitnan'); delta_max=max(delta_max,0.1);

    fig=figure('Color','w','Visible','off','Position',[40 60 1420 690]);
    tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
    for i=1:8
        ax=nexttile(tl); imagesc(ax,x_mm,z_mm,maps{i});
        axis(ax,'image'); set(ax,'YDir','normal','FontSize',8);
        title(ax,titles(i),'FontSize',10,'FontWeight','normal');
        if i<=3, clim(ax,[1.5 3.5]); colormap(ax,turbo); colorbar(ax);
        elseif i<=5, clim(ax,[0 error_max]); colormap(ax,turbo); colorbar(ax);
        elseif i==6, clim(ax,[0 1]); colormap(ax,parula); colorbar(ax);
        elseif i==7, clim(ax,[0 1]); colormap(ax,gray); colorbar(ax,'Ticks',[0 1]);
        else, clim(ax,[-delta_max delta_max]); colormap(ax,bluewhitered_map()); colorbar(ax);
        end
        if i>4, xlabel(ax,'x (mm)'); end
        if mod(i-1,4)==0, ylabel(ax,'z (mm)'); end
    end
    heading=sprintf('%s | %s | f=%g Hz | M=%g | dx=%.1f mm', ...
        pretty_geometry(C.geometry(1)),pretty_regime(C.field_regime(1)), ...
        C.f0(1),C.M(1),1e3*C.dx(1));
    title(tl,heading,'Interpreter','none','FontSize',12,'FontWeight','bold');
    hide_axes_toolbars(fig);
    exportgraphics(fig,file,'Resolution',180); close(fig);
    if mod(ki,25)==0||ki==numel(keys)
        fprintf('  Condition maps: %d/%d.\n',ki,numel(keys));
    end
end
end

function label=pretty_geometry(x)
switch string(x)
    case "bilayer_2_3", label='Bilayer 2/3 m/s';
    case "circular_inclusion_2_3", label='Circular inclusion 2/3 m/s';
    case "homogeneous_cs2", label='Homogeneous 2 m/s';
    case "homogeneous_cs3", label='Homogeneous 3 m/s';
    otherwise, label=char(x);
end
end

function label=pretty_regime(x)
switch string(x)
    case "directional_2D", label='Directional 2D';
    case "diffuse_2D", label='Diffuse 2D';
    case "partial_3D", label='Partial 3D';
    case "diffuse_3D", label='Diffuse 3D';
    otherwise, label=char(x);
end
end

function plot_baseline_vs_corrected(T,CORR,OUT)
X=primary_rows(T,CORR); chosen=X(X.strategy_name==CORR.RepresentativeStrategy,:);
base=X(X.strategy_name=="baseline",:);
keys={'condition_key','map_iz','map_ix'};
L=base(:,keys); L.baseline_error=base.SWS_error_abs_pct;
R=chosen(:,keys); R.corrected_error=chosen.SWS_error_abs_pct;
R.patch_purity=chosen.patch_purity;
J=innerjoin(L,R,'Keys',keys);
if height(J)>25000, rng(CORR.RandomSeed); J=J(randperm(height(J),25000),:); end
fig=figure('Color','w','Position',[100 100 720 620]);
scatter(J.baseline_error,J.corrected_error,11,J.patch_purity, ...
    'filled','MarkerFaceAlpha',0.22); hold on;
lim=prctile([J.baseline_error;J.corrected_error],99.5); lim=max(lim,20);
plot([0 lim],[0 lim],'k--','LineWidth',1); xlim([0 lim]); ylim([0 lim]); axis square;
xlabel('Baseline absolute error (%)'); ylabel('Corrected absolute error (%)');
title('Pixel-wise effect of gated interpolation'); style_axes(gca);
cb=colorbar; cb.Label.String='Patch purity'; clim([0 1]); colormap(turbo);
changed=abs(J.corrected_error-J.baseline_error)>1e-9;
improved=100*mean(J.corrected_error(changed)<J.baseline_error(changed));
changed_pct=100*mean(changed);
text(0.04*lim,0.91*lim,sprintf('Changed: %.1f%% | improved among changed: %.1f%%', ...
    changed_pct,improved), ...
    'FontSize',10,'BackgroundColor','w','Margin',4);
export_fig(fig,OUT,'test26_baseline_vs_corrected_scatter.png');
end

function plot_strategy_heatmap(T,CORR,OUT)
X=primary_rows(T,CORR); strategies=plot_strategy_subset(X,CORR);
geometries=["bilayer_2_3","circular_inclusion_2_3","homogeneous_cs2","homogeneous_cs3"];
geometry_titles=["Bilayer 2/3","Circular inclusion 2/3","Homogeneous 2","Homogeneous 3"];
fig=figure('Color','w','Position',[80 60 1120 680]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
all_values=[]; axes_list=gobjects(4,1);
for gi=1:4
    G=X(X.geometry==geometries(gi),:);
    groups=unique(G(:,{'field_regime','M'}),'rows','stable');
    A=nan(numel(strategies),height(groups));
    for si=1:numel(strategies)
        for j=1:height(groups)
            idx=G.strategy_name==strategies(si)& ...
                G.field_regime==groups.field_regime(j)&G.M==groups.M(j);
            A(si,j)=mean(G.SWS_error_abs_pct(idx),'omitnan');
        end
    end
    all_values=[all_values;A(:)]; %#ok<AGROW>
    ax=nexttile(tl); axes_list(gi)=ax; imagesc(ax,A);
    labels=replace(groups.field_regime,["directional_2D","diffuse_3D"],["Directional","Diffuse 3D"])+" / M"+string(groups.M);
    set(ax,'YTick',1:numel(strategies),'YTickLabel',short_strategy(strategies), ...
        'XTick',1:height(groups),'XTickLabel',labels,'XTickLabelRotation',15, ...
        'TickLabelInterpreter','none','FontSize',8);
    title(ax,geometry_titles(gi),'FontSize',11); colorbar(ax); colormap(ax,turbo);
end
limits=[0 prctile(all_values,97)];
for ax=axes_list', clim(ax,limits); end
title(tl,'MAPE by geometry, regime and M','FontSize',13,'FontWeight','bold');
export_fig(fig,OUT,'test26_strategy_condition_heatmap.png');
end

function X=primary_rows(T,CORR)
X=T(T.sws_model==CORR.PrimaryModel&T.detector==CORR.PrimaryDetector,:);
end

function strategies=plot_strategy_subset(X,CORR)
wanted=["baseline","physical_donut_all", ...
    "confidence_gated_prior_donut_c070_a070_130_physical", ...
    CORR.RepresentativeStrategy,"peak_candidate_low_confidence_c080"];
strategies=intersect(wanted,unique(X.strategy_name,'stable'),'stable');
end

function labels=short_strategy(x)
labels=strings(size(x));
for i=1:numel(x)
    switch x(i)
        case "baseline", labels(i)="Baseline";
        case "physical_donut_all", labels(i)="Physical donut";
        case "confidence_gated_prior_donut_c070_a070_130_physical", labels(i)="Gated prior donut";
        case "confidence_gated_interpolation_c080", labels(i)="Gated interpolation";
        case "peak_candidate_low_confidence_c080", labels(i)="Peak candidate";
        otherwise, labels(i)=replace(x(i),"_"," ");
    end
end
end

function colors=strategy_colors(n)
base=[0.15 0.35 0.65;0.85 0.45 0.15;0.55 0.35 0.70; ...
    0.10 0.60 0.50;0.45 0.45 0.45];
colors=base(1:n,:);
end

function set_bar_colors(bars,colors)
for i=1:numel(bars), bars(i).FaceColor=colors(i,:); end
end

function style_axes(ax)
set(ax,'FontName','Helvetica','FontSize',9,'LineWidth',0.8, ...
    'Box','off','TickDir','out');
grid(ax,'on'); ax.GridAlpha=0.16; ax.MinorGridAlpha=0.08;
end

function cmap=bluewhitered_map()
n=256; half=n/2;
blue=[linspace(0.15,1,half)' linspace(0.35,1,half)' ones(half,1)];
red=[ones(half,1) linspace(1,0.25,half)' linspace(1,0.20,half)'];
cmap=[blue;red];
end

function print_interpretation_summary(S,CORR)
X=S(S.sws_model==CORR.PrimaryModel& ...
    S.detector==CORR.PrimaryDetector,:);
baseline=X(X.strategy_name=="baseline",:);
[~,ig]=min(X.MAPE);
mixed=mean([X.MAPE_moderately_mixed X.MAPE_strongly_mixed],2,'omitnan');
[~,im]=min(mixed);
fprintf('\nInterpretive summary (%s / %s):\n', ...
    CORR.PrimaryModel,CORR.PrimaryDetector);
fprintf('  Best global MAPE: %s (%.2f%%).\n', ...
    X.strategy_name(ig),X.MAPE(ig));
fprintf('  Best mixed-patch MAPE: %s (%.2f%%).\n', ...
    X.strategy_name(im),mixed(im));
physical=X(X.strategy_name=="physical_donut_all",:);
if ~isempty(physical)
    fprintf('  Physical donut homogeneous delta: %+.2f points.\n', ...
        physical.MAPE_homogeneous-baseline.MAPE_homogeneous);
end
selected=X(X.strategy_name==CORR.RepresentativeStrategy,:);
if ~isempty(selected)
    fprintf('  Selected gated strategy pure delta: %+.2f points.\n', ...
        selected.MAPE_pure-baseline.MAPE_pure);
end
fprintf(['  Recommendation: accept a gated correction only when both ', ...
    'heterogeneous geometries improve without >0.5-1 point ', ...
    'homogeneous degradation.\n']);
end

function [A,iz,ix]=map_from_rows(T,v)
[z_values,~,iz]=unique(T.map_iz,'sorted');
[x_values,~,ix]=unique(T.map_ix,'sorted');
A=nan(numel(z_values),numel(x_values));
A(sub2ind(size(A),iz,ix))=v;
end

function export_fig(fig,OUT,name)
hide_axes_toolbars(fig);
exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',220); close(fig);
end

function hide_axes_toolbars(fig)
axes_list=findall(fig,'Type','axes');
for ax=axes_list'
    try
        ax.Toolbar.Visible='off';
    catch
    end
end
end

function validate_helper()
k=linspace(0,5000,101)'; c=struct('k_cent',k,'Srad',ones(size(k)), ...
    'Ecum',linspace(0,1,numel(k))');
[m,info]=adaptive_req.analysis.apply_radial_donut(c,'F0',500, ...
    'ApplyPhysical',true,'PhysicalSpeedRange',[0.5 10], ...
    'ApplyPrior',true,'PriorK',1500,'RelativeRange',[0.7 1.3]);
assert(abs(m.Ecum(end)-1)<1e-12&&info.k_low<info.k_high);
end

function value=ternary(condition,a,b)
if condition, value=a; else, value=b; end
end

function s=sanitize(x)
s=regexprep(char(x),'[^A-Za-z0-9_-]','_');
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
missingA=setdiff(B.Properties.VariableNames,A.Properties.VariableNames);
missingB=setdiff(A.Properties.VariableNames,B.Properties.VariableNames);
for x=missingA, A.(x{1})=missing_like(B.(x{1}),height(A)); end
for x=missingB, B.(x{1})=missing_like(A.(x{1}),height(B)); end
B=B(:,A.Properties.VariableNames); T=[A;B];
end

function x=missing_like(example,n)
if isstring(example), x=repmat("",n,1);
elseif islogical(example), x=false(n,1);
elseif isnumeric(example), x=nan(n,1);
else, x=cell(n,1);
end
end
