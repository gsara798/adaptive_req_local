%% analyze_test_29_oracle_halfwindow_graph_tv.m
% Test 29: oracle ceiling, operational half-window bank, and graph/TV map.
%
% Stage 1 is diagnostic_only: a true same-material mask estimates whether
% removing the opposite material can rescue REQ. Stage 2 evaluates a bank of
% smoothly tapered half-windows without oracle information. Stage 3 applies
% confidence-weighted bilateral graph/TV reconstruction. Full mode is blocked
% unless an operational quick strategy clearly improves mixed patches without
% degrading homogeneous or pure/near-pure regions.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST29_MODE          = quick | full
%   ADAPTIVE_REQ_TEST29_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST29_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST29_FORCE_FULL    = true | false

clear; clc; close all;
format compact;

%% Configuration

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.15, ...
    'defaultAxesLabelFontSizeMultiplier',1.05);

CFG=struct(); mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST29_MODE'))));
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]),'ADAPTIVE_REQ_TEST29_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST29_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST29_SAVE_ALL_MAPS',true);
CFG.ForceFull=env_true('ADAPTIVE_REQ_TEST29_FORCE_FULL',false);
CFG.Version=2; CFG.Dx=0.2e-3; CFG.ConfidenceThreshold=.80;
CFG.PhysicalSwsRange=[.5 10]; CFG.MaxRelativeCorrection=.50;
CFG.HalfWindowCount=8; CFG.HalfTransitionPixels=1.5; CFG.MaskFloor=.03;
CFG.HalfWidthRatioMax=.95; CFG.HalfScoreImprovementMin=.005;
CFG.HalfPriorPenalty=.15;
CFG.GraphLambdas=[.25 .75]; CFG.GraphIterations=30;
CFG.GraphEdgeSigma=.20; CFG.GraphTvEpsilon=.03;
CFG.FullGateMixedImprovement=0.50;
CFG.FullGateHomogeneousDeltaMax=0.50;
CFG.FullGatePureDeltaMax=0.25;
CFG.RandomSeed=29001;
CFG.Geometries=["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"];
if CFG.QuickMode
    CFG.Frequencies=500; CFG.M=[2 3];
    CFG.Regimes=["directional_2D","diffuse_3D"];
    CFG.RunOracle=true;
    CFG.OperationalStrategies=["hybrid_baseline","local_baseline", ...
        "confidence_switch_c080","halfwindow_local_c080", ...
        "halfwindow_switch_c080","graph_switch_l025","graph_switch_l075", ...
        "graph_halfwindow_switch_l025","graph_halfwindow_switch_l075"];
    CFG.DiagnosticStrategies=["oracle_same_material_c080","theory_discrete"];
else
    CFG.Frequencies=[300 400 500 600]; CFG.M=[2 3 4];
    CFG.Regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
    CFG.RunOracle=false;
    CFG.GraphLambdas=.75;
    CFG.OperationalStrategies=["hybrid_baseline","local_baseline", ...
        "confidence_switch_c080","halfwindow_switch_c080", ...
        "graph_halfwindow_switch_l075"];
    CFG.DiagnosticStrategies="theory_discrete";
end
CFG.AllStrategies=[CFG.OperationalStrategies CFG.DiagnosticStrategies];

TEST27=fullfile(root_dir,'outputs','test_27_adaptive_window_edge_aware');
TEST26=fullfile(root_dir,'outputs','test_26_confidence_gated_corrections');
if CFG.QuickMode
    quick27=fullfile(TEST27,'quick'); quick26=fullfile(TEST26,'quick'); use_quick=false;
    if exist(fullfile(quick27,'tables','test27_patch_level_results.csv'),'file')==2
        TEST27=quick27; use_quick=true;
    end
    if use_quick&&exist(fullfile(quick26,'data','spectral_checkpoints'),'dir')==7
        TEST26=quick26;
    end
end
SOURCE_TABLE=fullfile(TEST27,'tables','test27_patch_level_results.csv');
SPECTRAL_DIR=fullfile(TEST26,'data','spectral_checkpoints');
FIELD_DIR=fullfile(root_dir,'outputs','test_22_confidence_external_validation', ...
    'analysis','data','field_cache');
OUT=make_output_dirs(root_dir,CFG);
write_config_json(CFG,fullfile(OUT.root_dir,'test29_configuration.json'));

if ~CFG.QuickMode&&~CFG.ValidateOnly&&~CFG.ForceFull
    enforce_full_gate(root_dir,CFG);
end

fprintf('\nTest 29: oracle ceiling, half-window bank, and graph/TV\n');
fprintf('Mode: %s | oracle strategy is diagnostic_only.\n',CFG.RunMode);
assert(exist(SOURCE_TABLE,'file')==2,'Test 27 table missing: %s',SOURCE_TABLE);
T=readtable(SOURCE_TABLE,'TextType','string'); T=filter_design(T,CFG);
keys=unique(T.condition_key,'stable');
assert(~isempty(keys),'No conditions match the Test 29 design.');
fprintf('Matched %d conditions and %d patches.\n',numel(keys),height(T));

if CFG.ValidateOnly
    validate_test29(T,keys(1),SPECTRAL_DIR,FIELD_DIR,CFG);
    fprintf('Test 29 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Evaluate conditions

parts=cell(numel(keys),1);
for ci=1:numel(keys)
    timer=tic; key=keys(ci); checkpoint=fullfile(OUT.condition_dir, ...
        "test29__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        S=load(checkpoint,'RESULT');
        if isfield(S,'RESULT')&&isfield(S.RESULT,'version')&&S.RESULT.version==CFG.Version
            R=S.RESULT;
            if ~ismember('sws_theory_discrete',R.T_patch.Properties.VariableNames)
                X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
                SPEC=load_spectral_cache(SPECTRAL_DIR,key,X);
                R.T_patch.sws_theory_discrete=theory_discrete_prediction(X,SPEC);
                RESULT=R; save(checkpoint,'RESULT','-v7.3');
            end
            parts{ci}=R.T_patch;
            fprintf('[%d/%d] Reused %s (%d patches).\n',ci,numel(keys),key,height(parts{ci}));
            continue;
        end
    end
    X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
    SPEC=load_spectral_cache(SPECTRAL_DIR,key,X);
    F=load_field_cache(FIELD_DIR,X(1,:));
    [pred,diag]=evaluate_condition(X,SPEC,F,CFG,inf);
    assert_policy(pred,diag,X,CFG);
    T_patch=attach_outputs(X,pred,diag,CFG);
    RESULT=struct('version',CFG.Version,'T_patch',T_patch); save(checkpoint,'RESULT','-v7.3');
    parts{ci}=T_patch;
    fprintf('[%d/%d] %s: low %d, oracle %d, half-window %d in %.1f s.\n', ...
        ci,numel(keys),key,sum(diag.low_confidence),sum(diag.oracle_applied), ...
        sum(diag.halfwindow_accepted),toc(timer));
end
T_patch=vertcat(parts{:}); clear parts T;

%% Summaries and full gate

T_overall=summarize(T_patch,"strategy_name",CFG);
T_geometry=summarize(T_patch,["strategy_name","geometry"],CFG);
T_regime=summarize(T_patch,["strategy_name","field_regime"],CFG);
T_frequency_M=summarize(T_patch,["strategy_name","f0","M"],CFG);
T_purity=summarize(T_patch,["strategy_name","geometry","purity_bin"],CFG);
T_distance=summarize(T_patch,["strategy_name","geometry","distance_bin"],CFG);
T_acceptance=summarize_acceptance(T_patch);
T_gate=make_full_gate(T_overall,CFG);

writetable(T_patch,fullfile(OUT.table_dir,'test29_patch_level_results.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test29_strategy_summary_overall.csv'));
writetable(T_geometry,fullfile(OUT.table_dir,'test29_strategy_summary_by_geometry.csv'));
writetable(T_regime,fullfile(OUT.table_dir,'test29_strategy_summary_by_regime.csv'));
writetable(T_frequency_M,fullfile(OUT.table_dir,'test29_strategy_summary_by_frequency_M.csv'));
writetable(T_purity,fullfile(OUT.table_dir,'test29_strategy_summary_by_purity.csv'));
writetable(T_distance,fullfile(OUT.table_dir,'test29_strategy_summary_by_distance.csv'));
writetable(T_acceptance,fullfile(OUT.table_dir,'test29_halfwindow_acceptance.csv'));
writetable(T_gate,fullfile(OUT.table_dir,'test29_full_gate.csv'));
save(fullfile(OUT.data_dir,'test29_compact_results.mat'),'T_patch','T_overall','T_gate','CFG','-v7.3');

plot_summary(T_patch,CFG,OUT);
if CFG.RunOracle, plot_oracle_ceiling(T_patch,CFG,OUT); end
plot_error_bins(T_patch,CFG,OUT,"purity_bin", ...
    ["strongly_mixed","moderately_mixed","near_pure","pure"], ...
    ["Strongly mixed","Moderately mixed","Near-pure","Pure"], ...
    'Patch class','test29_error_by_purity.png');
plot_error_bins(T_patch,CFG,OUT,"distance_bin", ...
    ["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"], ...
    ["0-1","1-2","2-4","4-8",">8"], ...
    'Distance to interface (mm)','test29_error_by_distance.png');
plot_acceptance(T_patch,OUT); plot_representative_maps(T_patch,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T_patch,CFG,OUT); end
print_summary(T_overall,T_gate,T_acceptance,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 29 complete.\n',OUT.table_dir,OUT.figure_dir);

%% Setup and gate helpers

function OUT=make_output_dirs(root,CFG)
OUT.root_dir=fullfile(root,'outputs','test_29_oracle_halfwindow_graph_tv');
if CFG.QuickMode, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables'); OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data'); OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
for d=string(struct2cell(OUT))', if exist(d,'dir')~=7, mkdir(d); end, end
end

function tf=env_true(name,default)
value=strtrim(getenv(name)); if isempty(value), tf=default;
else, tf=any(strcmpi(value,{'true','1','yes','on'})); end
end

function write_config_json(CFG,file)
C=CFG; for f=string(fieldnames(C))', if isstring(C.(f)), C.(f)=cellstr(C.(f)); end, end
fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char'); clear cleanup;
end

function enforce_full_gate(root,~)
file=fullfile(root,'outputs','test_29_oracle_halfwindow_graph_tv','quick', ...
    'tables','test29_full_gate.csv');
assert(exist(file,'file')==2,["Quick gate is missing. Run Test 29 quick first: " file]);
G=readtable(file,'TextType','string'); passed=G.gate_pass&G.operational;
assert(any(passed),['Test 29 full is blocked: no quick operational strategy ', ...
    'clearly improved mixed patches without degrading homogeneous/pure regions.']);
fprintf('Full gate passed by: %s\n',strjoin(G.strategy_name(passed),', '));
end

function T=filter_design(T,CFG)
T=T(abs(T.dx-CFG.Dx)<1e-12&ismember(T.geometry,CFG.Geometries)& ...
    ismember(T.f0,CFG.Frequencies)&ismember(T.M,CFG.M)&ismember(T.field_regime,CFG.Regimes),:);
end

function SPEC=load_spectral_cache(folder,key,X)
file=fullfile(folder,"spectra__"+sanitize(key)+".mat"); assert(exist(file,'file')==2,'Missing %s',file);
S=load(file,'SPEC'); SPEC=S.SPEC; assert(numel(SPEC.curves)==height(X));
cx=round(X.x/X.dx(1))+1; cz=round(X.z/X.dz(1))+1;
assert(isequal(double(SPEC.cx(:)),double(cx(:)))&&isequal(double(SPEC.cz(:)),double(cz(:))));
end

function F=load_field_cache(folder,row)
key=sprintf('%s__f%g__%s__dx%gum',row.geometry,row.f0,lower(string(row.field_regime)),round(1e6*row.dx));
file=fullfile(folder,"field__"+sanitize(key)+".mat"); assert(exist(file,'file')==2,'Missing %s',file);
F=load(file,'sim','cfg_sim');
end

%% Spectral candidates and graph reconstruction

function [pred,D]=evaluate_condition(X,SPEC,F,CFG,max_low)
n=height(X); local=double(X.sws_local_baseline); hybrid=double(X.sws_hybrid_baseline);
switch_map=double(X.sws_confidence_switch_local_hybrid_c080); low=X.confidence<CFG.ConfidenceThreshold;
idx=find(low); if isfinite(max_low), idx=idx(1:min(numel(idx),max_low)); low(:)=false; low(idx)=true; end

oracle=local; half=local; oracle_applied=false(n,1); half_accept=false(n,1);
half_angle=nan(n,1); half_width_ratio=nan(n,1); half_candidate=nan(n,1);
feat=adaptive_req.config.default_feature_config('M',X.M(1),'cs_guess',3,'gamma_win',1,'pad_factor',1);
[req,feat]=adaptive_req.config.default_req_config(F.cfg_sim,feat,'Nbins','auto', ...
    'Nbins_auto_oversample',1,'Nbins_min',16,'smooth_sigma',1);
half_size=floor(feat.win_size/2); cfg_patch=F.cfg_sim; cfg_patch.UseParfor=false;
angles=linspace(0,2*pi,CFG.HalfWindowCount+1); angles(end)=[];

count=numel(idx); oracle_value=local(idx); oracle_used=false(count,1);
half_value=local(idx); half_used=false(count,1); angle_value=nan(count,1);
width_value=nan(count,1); candidate_value=nan(count,1);
run_oracle=CFG.RunOracle;
parfor jj=1:count
    i=idx(jj); zi=(SPEC.cz(i)-half_size):(SPEC.cz(i)+half_size);
    xi=(SPEC.cx(i)-half_size):(SPEC.cx(i)+half_size);
    if min(zi)<1||max(zi)>size(F.sim.Uxz,1)||min(xi)<1||max(xi)>size(F.sim.Uxz,2), continue; end
    patch=F.sim.Uxz(zi,xi); original=curve_quality(SPEC.curves{i});
    qlocal=q_at_sws(SPEC.curves{i},local(i),X.f0(i));

    % Diagnostic-only oracle mask.
    if run_oracle
        material=F.sim.cs_map(zi,xi); center=F.sim.cs_map(SPEC.cz(i),SPEC.cx(i));
        same=double(abs(material-center)<.25);
        if mean(same(:))<.99
            weight=smooth_mask(same,CFG.MaskFloor);
            value=masked_sws(patch,weight,cfg_patch,req,qlocal,X.f0(i));
            if candidate_ok(value,local(i),CFG)
                oracle_value(jj)=value; oracle_used(jj)=true;
            end
        end
    end

    [ZZ,XX]=ndgrid(-half_size:half_size,-half_size:half_size);
    best_objective=Inf; best_value=NaN; best_angle=NaN; best_ratio=NaN;
    for theta=angles
        signed=XX*cos(theta)+ZZ*sin(theta);
        weight=CFG.MaskFloor+(1-CFG.MaskFloor)./(1+exp(signed/CFG.HalfTransitionPixels));
        try
            [value,curve]=masked_sws(patch,weight,cfg_patch,req,qlocal,X.f0(i));
            quality=curve_quality(curve); ratio=quality.width/max(original.width,eps);
            improvement=original.mixture_score-quality.mixture_score;
            acceptable=(ratio<=CFG.HalfWidthRatioMax||improvement>=CFG.HalfScoreImprovementMin)&& ...
                candidate_ok(value,local(i),CFG);
            objective=ratio+CFG.HalfPriorPenalty*abs(log(value/local(i)))+.25*quality.mixture_score;
            if acceptable&&objective<best_objective
                best_objective=objective; best_value=value; best_angle=theta; best_ratio=ratio;
            end
        catch
        end
    end
    if isfinite(best_value)
        half_value(jj)=best_value; half_used(jj)=true; angle_value(jj)=best_angle;
        width_value(jj)=best_ratio; candidate_value(jj)=best_value;
    end
end
oracle(idx)=oracle_value; oracle_applied(idx)=oracle_used;
half(idx)=half_value; half_accept(idx)=half_used; half_angle(idx)=angle_value;
half_width_ratio(idx)=width_value; half_candidate(idx)=candidate_value;

half_switch=switch_map; half_switch(half_accept)=half(half_accept);
pred=struct('hybrid_baseline',hybrid,'local_baseline',local, ...
    'confidence_switch_c080',switch_map,'halfwindow_local_c080',half, ...
    'halfwindow_switch_c080',half_switch);
if CFG.RunOracle, pred.oracle_same_material_c080=oracle; end
pred.theory_discrete=theory_discrete_prediction(X,SPEC);
for lambda=CFG.GraphLambdas
    tag="l"+sprintf('%03d',round(100*lambda));
    pred.("graph_switch_"+tag)=graph_correct(X,switch_map,local,lambda,CFG);
    pred.("graph_halfwindow_switch_"+tag)=graph_correct(X,half_switch,local,lambda,CFG);
end
D=struct('low_confidence',low,'oracle_applied',oracle_applied, ...
    'halfwindow_accepted',half_accept,'halfwindow_angle_rad',half_angle, ...
    'halfwindow_width_ratio',half_width_ratio,'halfwindow_candidate_sws',half_candidate);
end

function corrected=graph_correct(X,observed,seed,lambda,CFG)
[Y,iz,ix]=map_from_rows(X,observed); C=map_from_rows(X,X.confidence); S=map_from_rows(X,seed);
[U,~]=adaptive_req.analysis.confidence_weighted_graph_tv(Y,C,S, ...
    'Lambda',lambda,'Iterations',CFG.GraphIterations,'EdgeSigma',CFG.GraphEdgeSigma, ...
    'TvEpsilon',CFG.GraphTvEpsilon,'HighConfidenceThreshold',CFG.ConfidenceThreshold, ...
    'PhysicalRange',CFG.PhysicalSwsRange);
corrected=U(sub2ind(size(U),iz,ix));
end

function varargout=masked_sws(patch,weight,cfg,req,q,f0)
[~,curve]=adaptive_req.quantile.compute_quantile_from_patch(patch.*sqrt(weight),cfg,req);
value=adaptive_req.quantile.quantile_to_cs(curve,q,f0); varargout={value,curve};
end

function weight=smooth_mask(mask,floor_value)
kernel=[1 2 1]'*[1 2 1]/16; weight=conv2(mask,kernel,'same'); weight=conv2(weight,kernel,'same');
weight=floor_value+(1-floor_value)*min(max(weight,0),1);
end

function Q=curve_quality(curve)
D=adaptive_req.analysis.decompose_radial_spectrum(curve,'MinPeakHeightRatio',.08, ...
    'MinSeparationFraction',.10,'MaxValleyRatio',.9,'MinComponentWeight',.05);
Q=struct('width',D.radial_width_10_90,'mixture_score',D.two_component_score);
end

function q=q_at_sws(curve,sws,f0)
k=double(curve.k_cent(:)); E=double(curve.Ecum(:)); ok=isfinite(k)&isfinite(E);
[k,ia]=unique(k(ok),'stable'); E=E(ok); E=E(ia);
if numel(k)<2, q=NaN; else, q=interp1(k,E,2*pi*f0/sws,'linear','extrap'); end
q=min(max(q,0),1);
end

function sws=theory_discrete_prediction(X,SPEC)
regime=string(X.field_regime(1));
switch regime
    case "directional_2D", type="SingleWave";
    case "diffuse_2D", type="Diffuse2D";
    case "diffuse_3D", type="Diffuse3D";
    otherwise, type="Partial3D";
end
if type=="Partial3D"
    q=.5*(theory_q_one(X,"Diffuse2D")+theory_q_one(X,"Diffuse3D"));
else
    q=theory_q_one(X,type);
end
sws=nan(height(X),1);
for i=1:height(X)
    sws(i)=adaptive_req.quantile.quantile_to_cs(SPEC.curves{i},q,X.f0(i));
end
end

function q=theory_q_one(X,type)
o=adaptive_req.theory.q_theory_REQ_discrete_shearUZ(X.dx(1),X.dz(1), ...
    X.f0(1),3,'M',X.M(1),'Gamma',1,'PadFactor',1,'Nbins','auto', ...
    'SmoothSigma',1,'TheoryMode','S2D','FieldType',type,'Plot',false);
q=o.q_th;
end

function tf=candidate_ok(value,fallback,CFG)
tf=isfinite(value)&&value>=CFG.PhysicalSwsRange(1)&&value<=CFG.PhysicalSwsRange(2)&& ...
    abs(value-fallback)/max(fallback,eps)<=CFG.MaxRelativeCorrection;
end

function assert_policy(pred,D,X,CFG)
for name=CFG.AllStrategies
    value=pred.(name); assert(numel(value)==height(X)&&all(isfinite(value))&& ...
        all(value>=CFG.PhysicalSwsRange(1)&value<=CFG.PhysicalSwsRange(2)));
end
assert(~any(D.halfwindow_accepted&~D.low_confidence));
operational=["local","hybrid","confidence","spectrum","coordinates","half_windows"];
oracle=["true_SWS","patch_purity","material_side","distance_to_interface","cs_map"];
assert(isempty(intersect(lower(operational),lower(oracle))));
end

%% Outputs and statistics

function R=attach_outputs(X,pred,D,CFG)
keep={'condition_key','geometry','geometry_type','field_regime','f0','M','dx','dz', ...
    'map_iz','map_ix','x','z','true_SWS','confidence','patch_purity','material_side', ...
    'purity_bin','distance_bin','distance_to_interface_mm'};
R=X(:,keep); for name=CFG.AllStrategies, R.("sws_"+name)=pred.(name); end
for f=string(fieldnames(D))', R.(f)=D.(f); end
end

function S=summarize(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for strategy=CFG.AllStrategies
    pred=T.("sws_"+strategy); error=100*(pred-T.true_SWS)./T.true_SWS;
    if isempty(basegroups), G=ones(height(T),1); else, G=findgroups(T(:,cellstr(basegroups))); end
    for id=unique(G,'stable')'
        idx=G==id; row=table();
        for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy; row.diagnostic_only=ismember(strategy,CFG.DiagnosticStrategies);
        row.N=sum(idx); row.MAPE=mean(abs(error(idx))); row.mean_signed_error=mean(error(idx));
        row.high_error_10_pct=100*mean(abs(error(idx))>10); row.high_error_20_pct=100*mean(abs(error(idx))>20);
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            j=idx&T.purity_bin==cls; token=char(cls); row.("N_"+token)=sum(j);
            row.("MAPE_"+token)=mean(abs(error(j)),'omitnan');
            row.("high20_"+token)=100*mean(abs(error(j))>20,'omitnan');
        end
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function S=summarize_acceptance(T)
groups=["geometry","field_regime","f0","M"]; [G,S]=findgroups(T(:,cellstr(groups)));
S.N=splitapply(@numel,T.confidence,G); S.low_confidence_pct=100*splitapply(@mean,T.low_confidence,G);
S.oracle_applied_pct=100*splitapply(@mean,T.oracle_applied,G);
S.halfwindow_accepted_pct=100*splitapply(@mean,T.halfwindow_accepted,G);
S.mean_accepted_width_ratio=splitapply(@(x)mean(x,'omitnan'),T.halfwindow_width_ratio,G);
end

function G=make_full_gate(S,CFG)
base=S(S.strategy_name=="confidence_switch_c080",:); G=S(:,{'strategy_name','diagnostic_only'});
G.operational=~G.diagnostic_only; mixed=mean([S.MAPE_moderately_mixed S.MAPE_strongly_mixed],2,'omitnan');
base_mixed=mean([base.MAPE_moderately_mixed base.MAPE_strongly_mixed],2,'omitnan');
pure=mean([S.MAPE_pure S.MAPE_near_pure],2,'omitnan'); base_pure=mean([base.MAPE_pure base.MAPE_near_pure],2,'omitnan');
G.mixed_improvement_points=base_mixed-mixed;
G.homogeneous_delta_points=S.MAPE_homogeneous-base.MAPE_homogeneous;
G.pure_near_delta_points=pure-base_pure;
G.gate_pass=G.operational&G.strategy_name~="confidence_switch_c080"& ...
    G.mixed_improvement_points>=CFG.FullGateMixedImprovement& ...
    G.homogeneous_delta_points<=CFG.FullGateHomogeneousDeltaMax& ...
    G.pure_near_delta_points<=CFG.FullGatePureDeltaMax;
end

%% Figures and reporting

function plot_summary(T,CFG,OUT)
M=nan(numel(CFG.AllStrategies),3); H=M;
for i=1:numel(CFG.AllStrategies)
    e=abs(100*(T.("sws_"+CFG.AllStrategies(i))-T.true_SWS)./T.true_SWS);
    masks={T.purity_bin=="homogeneous",ismember(T.purity_bin,["pure","near_pure"]), ...
        ismember(T.purity_bin,["moderately_mixed","strongly_mixed"])};
    for j=1:3, M(i,j)=mean(e(masks{j}),'omitnan'); H(i,j)=100*mean(e(masks{j})>20,'omitnan'); end
end
fig=figure('Color','w','Position',[50 40 1300 720]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,M); xlabel(ax,'MAPE (%)'); title(ax,'Mean absolute error'); style(ax);
ax=nexttile(tl); barh(ax,H); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style(ax);
for ax=findall(fig,'Type','axes')', set(ax,'YTick',1:numel(CFG.AllStrategies), ...
        'YTickLabel',pretty(CFG.AllStrategies),'TickLabelInterpreter','none','FontSize',8); end
legend(["Homogeneous","Pure / near-pure","Mixed"],'Location','northoutside', ...
    'Orientation','horizontal','Box','off'); export_fig(fig,OUT,'test29_strategy_summary.png');
end

function plot_oracle_ceiling(T,~,OUT)
strategies=["confidence_switch_c080","halfwindow_switch_c080", ...
    "graph_halfwindow_switch_l025","theory_discrete","oracle_same_material_c080"];
cases=["bilayer_2_3","circular_inclusion_2_3"]; values=nan(2,numel(strategies));
for c=1:2
    for s=1:numel(strategies)
        idx=T.geometry==cases(c)&ismember(T.purity_bin,["moderately_mixed","strongly_mixed"]);
        e=abs(100*(T.("sws_"+strategies(s))-T.true_SWS)./T.true_SWS);
        values(c,s)=mean(e(idx),'omitnan');
    end
end
fig=figure('Color','w','Position',[100 100 850 430]); ax=axes(fig); bar(ax,values);
set(ax,'XTick',1:2,'XTickLabel',["Bilayer","Inclusion"]); ylabel(ax,'Mixed-patch MAPE (%)');
legend(ax,pretty(strategies),'Location','northoutside','Orientation','horizontal','Box','off'); style(ax);
export_fig(fig,OUT,'test29_oracle_ceiling.png');
end

function plot_error_bins(T,CFG,OUT,var,bins,labels,xlab,file)
cases=["bilayer_2_3","circular_inclusion_2_3"];
strategies=["confidence_switch_c080","halfwindow_switch_c080", ...
    "graph_halfwindow_switch_l075"];
if CFG.QuickMode
    strategies=["confidence_switch_c080","halfwindow_switch_c080", ...
        "graph_switch_l025","graph_halfwindow_switch_l025"];
end
strategies(end+1)="theory_discrete";
if CFG.RunOracle, strategies(end+1)="oracle_same_material_c080"; end
fig=figure('Color','w','Position',[80 80 1150 470]); tl=tiledlayout(1,2,'TileSpacing','compact');
for c=1:2, ax=nexttile(tl); hold(ax,'on');
    for strategy=strategies
        e=abs(100*(T.("sws_"+strategy)-T.true_SWS)./T.true_SWS); y=nan(size(bins));
        for b=1:numel(bins), idx=T.geometry==cases(c)&T.(var)==bins(b); y(b)=mean(e(idx),'omitnan'); end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.5,'DisplayName',pretty(strategy));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',labels); xlabel(ax,xlab); ylabel(ax,'MAPE (%)'); title(ax,replace(cases(c),"_"," ")); style(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off'); export_fig(fig,OUT,file);
end

function plot_acceptance(T,OUT)
classes=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]; value=nan(numel(classes),3);
for i=1:numel(classes), idx=T.purity_bin==classes(i); value(i,:)=[100*mean(T.low_confidence(idx)), ...
        100*mean(T.oracle_applied(idx)),100*mean(T.halfwindow_accepted(idx))]; end
fig=figure('Color','w','Position',[100 100 900 430]); ax=axes(fig); bar(ax,value);
set(ax,'XTick',1:numel(classes),'XTickLabel',replace(classes,"_"," "),'XTickLabelRotation',18);
ylabel(ax,'Patches (%)'); legend(ax,["Low confidence","Oracle applied","Half-window accepted"], ...
    'Location','northoutside','Orientation','horizontal','Box','off'); style(ax);
export_fig(fig,OUT,'test29_acceptance_by_region.png');
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
plot_map(T(T.condition_key==keys(1),:),CFG,fullfile(OUT.figure_dir,'test29_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable'); root=fullfile(OUT.figure_dir,'maps_by_condition_dx200um');
for i=1:numel(keys), C=T(T.condition_key==keys(i),:); folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end; file=fullfile(folder,"test29_maps__"+sanitize(keys(i))+".png");
    if exist(file,'file')~=2, plot_map(C,CFG,file); end
end
end

function plot_map(C,CFG,file)
C=sortrows(C,{'map_iz','map_ix'});
if CFG.QuickMode, selected=C.sws_graph_halfwindow_switch_l025;
else, selected=C.sws_graph_halfwindow_switch_l075; end
e0=abs(100*(C.sws_confidence_switch_c080-C.true_SWS)./C.true_SWS); e1=abs(100*(selected-C.true_SWS)./C.true_SWS);
eTheory=abs(100*(C.sws_theory_discrete-C.true_SWS)./C.true_SWS);
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_confidence_switch_c080), ...
    map_from_rows(C,C.sws_halfwindow_switch_c080),map_from_rows(C,selected), ...
    map_from_rows(C,C.sws_theory_discrete),map_from_rows(C,e0),map_from_rows(C,e1), ...
    map_from_rows(C,eTheory), ...
    map_from_rows(C,C.confidence),map_from_rows(C,C.halfwindow_accepted), ...
    map_from_rows(C,C.halfwindow_angle_rad),map_from_rows(C,C.halfwindow_width_ratio)};
titles=["True","Switch","Half-window","Graph + half-window","Theory discrete", ...
    "Switch error","Selected error","Theory error","Confidence","Half accepted", ...
    "Half angle","Width ratio"];
if CFG.RunOracle
    maps=[maps(1:4) {map_from_rows(C,C.sws_oracle_same_material_c080)} maps(5:end)];
    titles=[titles(1:4) "Oracle ceiling" titles(5:end)];
end
fig=figure('Color','w','Visible','off','Position',[20 40 1500 820]); tl=tiledlayout(3,5,'TileSpacing','compact');
for i=1:numel(maps), ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i)); colorbar(ax); end
title(tl,C.condition_key(1),'Interpreter','none'); exportgraphics(fig,file,'Resolution',170); close(fig);
end

function print_summary(S,G,A,CFG)
oper=S(~S.diagnostic_only,:); [~,i]=min(oper.MAPE);
base=S(S.strategy_name=="confidence_switch_c080",:); mixed=@(X)mean([X.MAPE_moderately_mixed X.MAPE_strongly_mixed],2,'omitnan');
fprintf('\nInterpretation:\n  Best operational global: %s (%.2f%%).\n',oper.strategy_name(i),oper.MAPE(i));
if CFG.RunOracle
    oracle=S(S.strategy_name=="oracle_same_material_c080",:);
    fprintf('  Switch mixed MAPE %.2f%%; oracle ceiling %.2f%%.\n',mixed(base),mixed(oracle));
else
    fprintf('  Switch mixed MAPE %.2f%%. Oracle ceiling is retained in quick outputs.\n',mixed(base));
end
fprintf('  Half-window accepted %.1f%% of all patches.\n',sum(A.N.*A.halfwindow_accepted_pct)/sum(A.N));
passed=G.strategy_name(G.gate_pass); if isempty(passed), fprintf('  FULL GATE: FAIL. Do not run full yet.\n');
else, fprintf('  FULL GATE: PASS (%s).\n',strjoin(passed,', ')); end
end

%% Validation and generic helpers

function validate_test29(T,key,spectral_dir,field_dir,CFG)
Y=[2 2 2;2 NaN 3;3 3 3]; C=ones(3); C(2,2)=NaN; seed=Y;
[U,~]=adaptive_req.analysis.confidence_weighted_graph_tv(Y,C,seed,'Iterations',3);
assert(all(isfinite(U(isfinite(Y))))&&isnan(U(2,2)));
X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'}); SPEC=load_spectral_cache(spectral_dir,key,X); F=load_field_cache(field_dir,X(1,:));
[pred,D]=evaluate_condition(X,SPEC,F,CFG,3); assert_policy(pred,D,X,CFG);
assert(all(isfinite(pred.graph_halfwindow_switch_l025)));
fprintf('  Graph/TV synthetic map and three cached low-confidence patches passed.\n');
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end; B=B(:,A.Properties.VariableNames); T=[A;B];
end

function label=pretty(name)
label=replace(string(name),["hybrid_baseline","local_baseline","confidence_switch_c080", ...
    "halfwindow_local_c080","halfwindow_switch_c080","graph_switch_l025", ...
    "graph_switch_l075","graph_halfwindow_switch_l025","graph_halfwindow_switch_l075", ...
    "oracle_same_material_c080","theory_discrete"], ...
    ["Hybrid","Local","Switch","Half-window Local","Half-window Switch", ...
    "Graph Switch 0.25","Graph Switch 0.75","Graph Half-window 0.25", ...
    "Graph Half-window 0.75","Oracle same-material","Theory discrete"]);
end

function style(ax), grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9; end
function export_fig(fig,OUT,name), exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',200); close(fig); end
function value=sanitize(value), value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_'); end
