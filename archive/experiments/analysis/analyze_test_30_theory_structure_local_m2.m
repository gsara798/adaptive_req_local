%% analyze_test_30_theory_structure_local_m2.m
% Test 30: Theory structure, Local M2 levels, and boundary reconstruction.
%
% The output grid and quantitative observations always come from M=2.
% M=3/M=4 are interpolated only to estimate multiscale instability and can
% only reduce reliability. A user geometry guess selects a topology family
% (homogeneous/bilayer/inclusion/unknown), never a true mask or SWS value.
% Oracle material labels are attached after operational inference solely for
% evaluation and the explicitly diagnostic oracle_region_graph strategy.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST30_MODE          = quick | full
%   ADAPTIVE_REQ_TEST30_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST30_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST30_FORCE_FULL    = true | false

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.15, ...
    'defaultAxesLabelFontSizeMultiplier',1.05);

%% Configuration

CFG=struct(); mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST30_MODE'))));
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]),'ADAPTIVE_REQ_TEST30_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST30_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST30_SAVE_ALL_MAPS',true);
CFG.ForceFull=env_true('ADAPTIVE_REQ_TEST30_FORCE_FULL',false);
CFG.Version=2; CFG.OutputM=2; CFG.ConfidenceThreshold=.80;
CFG.UnknownMinSeparation=.25; CFG.SeedBoundaryMarginMm=2.0;
CFG.BoundaryBandMm=2.0; CFG.MinimumRegionSeeds=12;
CFG.MultiscaleSeedMax=.12; CFG.MultiscaleReliabilityScale=.12;
CFG.GraphLambda=.75; CFG.GraphIterations=35; CFG.GraphEdgeSigma=.15;
CFG.GraphTvEpsilon=.025; CFG.PhysicalRange=[.5 10];
CFG.FullGateMixedImprovement=.50; CFG.FullGateHomogeneousDelta=.50;
CFG.FullGatePureDelta=.25; CFG.RandomSeed=30001;
CFG.Geometries=["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"];
if CFG.QuickMode
    CFG.Frequencies=500; CFG.Regimes=["directional_2D","diffuse_3D"];
else
    CFG.Frequencies=[300 400 500 600];
    CFG.Regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
end
CFG.OperationalStrategies=["local_m2","switch_m2","theory_m2", ...
    "unknown_theory_calibrated","unknown_region_levels","unknown_boundary_graph", ...
    "userguess_theory_calibrated","userguess_region_levels", ...
    "userguess_boundary_graph","userguess_boundary_graph_multiscale"];
CFG.DiagnosticStrategies="oracle_region_graph";
CFG.AllStrategies=[CFG.OperationalStrategies CFG.DiagnosticStrategies];

SOURCE=fullfile(root_dir,'outputs','test_29_oracle_halfwindow_graph_tv', ...
    'data','test29_compact_results.mat');
assert(exist(SOURCE,'file')==2,'Full Test 29 compact results missing: %s',SOURCE);
OUT=make_output_dirs(root_dir,CFG); write_config_json(CFG,fullfile(OUT.root_dir,'test30_configuration.json'));
if ~CFG.QuickMode&&~CFG.ValidateOnly&&~CFG.ForceFull, enforce_full_gate(root_dir); end

fprintf('\nTest 30: Theory structure + Local M2 + boundary reconstruction\n');
fprintf('Mode: %s | output M=%d | M3/M4 are reliability-only.\n',CFG.RunMode,CFG.OutputM);
S=load(SOURCE,'T_patch'); T=S.T_patch; clear S;
T=T(ismember(T.geometry,CFG.Geometries)&ismember(T.f0,CFG.Frequencies)& ...
    ismember(T.field_regime,CFG.Regimes),:);
groups=unique(T(:,{'geometry','geometry_type','field_regime','f0','dx'}),'rows','stable');
fprintf('Matched %d physical fields (%d M-specific rows).\n',height(groups),height(T));

if CFG.ValidateOnly
    validate_test30(T,groups(1,:),CFG);
    fprintf('Test 30 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Operational maps per physical field

parts=cell(height(groups),1); ring_parts=cell(height(groups),1);
for gi=1:height(groups)
    timer=tic; key=physical_key(groups(gi,:));
    checkpoint=fullfile(OUT.condition_dir,"test30__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        C=load(checkpoint,'RESULT');
        if isfield(C,'RESULT')&&C.RESULT.version==CFG.Version
            parts{gi}=C.RESULT.T_patch; ring_parts{gi}=C.RESULT.T_ring;
            fprintf('[%d/%d] Reused %s.\n',gi,height(groups),key); continue;
        end
    end
    F=T(T.geometry==groups.geometry(gi)&T.field_regime==groups.field_regime(gi)& ...
        T.f0==groups.f0(gi)&abs(T.dx-groups.dx(gi))<1e-12,:);
    [R,T_ring]=evaluate_field(F,CFG);
    RESULT=struct('version',CFG.Version,'T_patch',R,'T_ring',T_ring);
    save(checkpoint,'RESULT','-v7.3'); parts{gi}=R; ring_parts{gi}=T_ring;
    fprintf('[%d/%d] %s: %d M2 pixels in %.2f s.\n', ...
        gi,height(groups),key,height(R),toc(timer));
end
T_patch=vertcat(parts{:}); T_ring=vertcat(ring_parts{:}); clear parts ring_parts T;

%% Summaries

T_overall=summarize_predictions(T_patch,"strategy_name",CFG);
T_geometry=summarize_predictions(T_patch,["strategy_name","geometry"],CFG);
T_regime=summarize_predictions(T_patch,["strategy_name","field_regime"],CFG);
T_frequency=summarize_predictions(T_patch,["strategy_name","f0"],CFG);
T_purity=summarize_predictions(T_patch,["strategy_name","geometry","purity_bin"],CFG);
T_gate=make_full_gate(T_overall,CFG);

writetable(T_patch,fullfile(OUT.table_dir,'test30_patch_level_results.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test30_strategy_summary_overall.csv'));
writetable(T_geometry,fullfile(OUT.table_dir,'test30_strategy_summary_by_geometry.csv'));
writetable(T_regime,fullfile(OUT.table_dir,'test30_strategy_summary_by_regime.csv'));
writetable(T_frequency,fullfile(OUT.table_dir,'test30_strategy_summary_by_frequency.csv'));
writetable(T_purity,fullfile(OUT.table_dir,'test30_strategy_summary_by_purity.csv'));
writetable(T_ring,fullfile(OUT.table_dir,'test30_inclusion_ring_core_metrics.csv'));
writetable(T_gate,fullfile(OUT.table_dir,'test30_full_gate.csv'));
save(fullfile(OUT.data_dir,'test30_compact_results.mat'),'T_patch','T_ring','T_overall','T_gate','CFG','-v7.3');

plot_strategy_summary(T_patch,CFG,OUT); plot_ring_profiles(T_patch,CFG,OUT);
plot_multiscale_instability(T_patch,OUT); plot_geometry_quality(T_patch,OUT);
plot_input_M_dependence(SOURCE,OUT); plot_representative_maps(T_patch,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T_patch,CFG,OUT); end
print_interpretation(T_overall,T_ring,T_gate,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 30 complete.\n',OUT.table_dir,OUT.figure_dir);

%% Setup helpers

function OUT=make_output_dirs(root,CFG)
OUT.root_dir=fullfile(root,'outputs','test_30_theory_structure_local_m2');
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

function enforce_full_gate(root)
file=fullfile(root,'outputs','test_30_theory_structure_local_m2','quick','tables','test30_full_gate.csv');
assert(exist(file,'file')==2,'Run Test 30 quick before full.'); G=readtable(file,'TextType','string');
passed=logical(G.gate_pass);
assert(any(passed),['Test 30 full blocked: no quick strategy improved mixed patches ', ...
    'without degrading homogeneous and pure regions.']);
fprintf('Full gate passed by: %s\n',strjoin(G.strategy_name(passed),', '));
end

function key=physical_key(G)
key=sprintf('%s__f%g__%s__dx%gum',G.geometry,G.f0,lower(G.field_regime),round(1e6*G.dx));
end

%% Field evaluation

function [R,T_ring]=evaluate_field(F,CFG)
M2=sortrows(F(F.M==2,:),{'map_iz','map_ix'}); assert(~isempty(M2),'M2 field missing.');
M3=sortrows(F(F.M==3,:),{'map_iz','map_ix'}); M4=sortrows(F(F.M==4,:),{'map_iz','map_ix'});
local2=double(M2.sws_local_baseline); switch2=double(M2.sws_confidence_switch_c080);
theory2=double(M2.sws_theory_discrete); confidence=double(M2.confidence);
local3=align_map(M3,double(M3.sws_local_baseline),M2);
local4=align_map(M4,double(M4.sws_local_baseline),M2);
instability=max(abs(local2-local3),abs(local3-local4))./max(local2,eps);
monotonic_drop=max(local2-local4,0)./max(local2,eps);

[unknown_label,unknown_info]=segment_theory_map(M2,theory2,"unknown",CFG);
guess=geometry_guess_from_case(M2.geometry_type(1));
[user_label,user_info]=segment_theory_map(M2,theory2,guess,CFG);

unknown=build_region_maps(M2,local2,theory2,confidence,instability, ...
    unknown_label,unknown_info,CFG);
user=build_region_maps(M2,local2,theory2,confidence,instability, ...
    user_label,user_info,CFG);

% Diagnostic-only ceiling: true material labels are never passed above.
oracle_label=ones(height(M2),1); oracle_label(M2.material_side=="hard")=2;
oracle_info=label_info(M2,oracle_label,"oracle",CFG);
oracle=build_region_maps(M2,local2,theory2,confidence,instability, ...
    oracle_label,oracle_info,CFG);

pred=struct('local_m2',local2,'switch_m2',switch2,'theory_m2',theory2, ...
    'unknown_theory_calibrated',unknown.theory_calibrated, ...
    'unknown_region_levels',unknown.region_levels, ...
    'unknown_boundary_graph',unknown.boundary_graph, ...
    'userguess_theory_calibrated',user.theory_calibrated, ...
    'userguess_region_levels',user.region_levels, ...
    'userguess_boundary_graph',user.boundary_graph, ...
    'userguess_boundary_graph_multiscale',user.boundary_graph_multiscale, ...
    'oracle_region_graph',oracle.boundary_graph_multiscale);
assert_operational_policy(pred,M2,CFG);
R=attach_outputs(M2,pred,CFG);
R.local_m3_aligned=local3; R.local_m4_aligned=local4;
R.multiscale_instability=instability; R.monotonic_M_drop=monotonic_drop;
R.unknown_region=unknown_label; R.userguess_region=user_label;
R.unknown_boundary_distance_mm=unknown_info.boundary_distance_mm;
R.userguess_boundary_distance_mm=user_info.boundary_distance_mm;
R.user_geometry_guess=repmat(guess,height(R),1);
R.unknown_geometry_guess=repmat("unknown",height(R),1);
R.userguess_region_iou=region_iou(user_label,oracle_label)*ones(height(R),1);
R.unknown_region_iou=region_iou(unknown_label,oracle_label)*ones(height(R),1);
T_ring=ring_metrics(R,CFG);
end

function values=align_map(source,source_values,target)
if isempty(source), values=nan(height(target),1); return; end
ok=isfinite(source_values); F=scatteredInterpolant(source.x(ok),source.z(ok), ...
    source_values(ok),'linear','nearest'); values=F(target.x,target.z);
end

function guess=geometry_guess_from_case(type)
switch string(type)
    case "homogeneous", guess="homogeneous";
    case "bilayer", guess="bilayer";
    case "inclusion", guess="inclusion";
    otherwise, guess="unknown";
end
end

%% Operational segmentation and reconstruction

function [labels,info]=segment_theory_map(meta,theory,guess,CFG)
[A,iz,ix]=map_from_rows(meta,theory); A=medfilt2(A,[3 3],'symmetric');
x=A(isfinite(A)); c=[prctile(x,25) prctile(x,75)];
for iter=1:20
    first=abs(x-c(1))<=abs(x-c(2));
    if any(first), c(1)=median(x(first)); end
    if any(~first), c(2)=median(x(~first)); end
end
c=sort(c); threshold=mean(c); separation=diff(c);
switch string(guess)
    case "homogeneous"
        L=ones(size(A));
    case "bilayer"
        col=median(A,1,'omitnan'); row=median(A,2,'omitnan');
        if range(col)>=range(row), profile=col; split=profile>threshold; L=repmat(split,size(A,1),1)+1;
        else, profile=row; split=profile>threshold; L=repmat(split,1,size(A,2))+1; end
    case "inclusion"
        raw=A>threshold; if mean(raw(:),'omitnan')>.5, raw=~raw; end
        raw=medfilt2(double(raw),[5 5],'symmetric')>.5;
        components=bwconncomp(raw,8);
        if components.NumObjects>0
            sizes=cellfun(@numel,components.PixelIdxList); [~,largest]=max(sizes);
            clean=false(size(raw)); clean(components.PixelIdxList{largest})=true;
            raw=imfill(clean,'holes');
        end
        L=ones(size(A)); L(raw)=2;
    otherwise
        if separation<CFG.UnknownMinSeparation, L=ones(size(A));
        else, raw=medfilt2(double(A>threshold),[3 3],'symmetric')>.5; L=ones(size(A)); L(raw)=2; end
end
labels=L(sub2ind(size(L),iz,ix)); info=label_info(meta,labels,guess,CFG);
info.cluster_separation=separation; info.threshold=threshold;
end

function info=label_info(meta,labels,guess,~)
[L,iz,ix]=map_from_rows(meta,labels); boundary=false(size(L));
boundary(:,1:end-1)=boundary(:,1:end-1)|(L(:,1:end-1)~=L(:,2:end));
boundary(:,2:end)=boundary(:,2:end)|(L(:,1:end-1)~=L(:,2:end));
boundary(1:end-1,:)=boundary(1:end-1,:)|(L(1:end-1,:)~=L(2:end,:));
boundary(2:end,:)=boundary(2:end,:)|(L(1:end-1,:)~=L(2:end,:));
if any(boundary(:)), distance=bwdist(boundary)*mean([median(diff(unique(meta.x))) median(diff(unique(meta.z)))]);
else, distance=inf(size(L)); end
info=struct('guess',string(guess),'boundary_distance_mm',1e3*distance(sub2ind(size(L),iz,ix)));
end

function maps=build_region_maps(meta,local,theory,confidence,instability,labels,info,CFG)
n=height(meta); region_levels=nan(n,1); theory_cal=theory;
for region=unique(labels(:))'
    member=labels==region; seed=member&confidence>=CFG.ConfidenceThreshold& ...
        info.boundary_distance_mm>=CFG.SeedBoundaryMarginMm&instability<=CFG.MultiscaleSeedMax;
    if nnz(seed)<CFG.MinimumRegionSeeds
        seed=member&confidence>=.65&info.boundary_distance_mm>=.5*CFG.SeedBoundaryMarginMm;
    end
    if nnz(seed)<CFG.MinimumRegionSeeds, seed=member&isfinite(local); end
    level=weighted_median(local(seed),max(confidence(seed),.05));
    residual=weighted_median(local(seed)-theory(seed),max(confidence(seed),.05));
    region_levels(member)=level; theory_cal(member)=theory(member)+residual;
end

effective_conf=confidence; boundary=info.boundary_distance_mm<=CFG.BoundaryBandMm;
effective_conf(boundary)=min(effective_conf(boundary),.08);
boundary_graph=graph_map(meta,local,effective_conf,region_levels,CFG);
multi_conf=effective_conf.*exp(-instability/CFG.MultiscaleReliabilityScale);
multi_conf(~isfinite(multi_conf))=effective_conf(~isfinite(multi_conf));
boundary_graph_multi=graph_map(meta,local,multi_conf,region_levels,CFG);
maps=struct('theory_calibrated',theory_cal,'region_levels',region_levels, ...
    'boundary_graph',boundary_graph,'boundary_graph_multiscale',boundary_graph_multi);
end

function corrected=graph_map(meta,observed,confidence,seed,CFG)
[Y,iz,ix]=map_from_rows(meta,observed); C=map_from_rows(meta,confidence); S=map_from_rows(meta,seed);
[U,~]=adaptive_req.analysis.confidence_weighted_graph_tv(Y,C,S, ...
    'Lambda',CFG.GraphLambda,'Iterations',CFG.GraphIterations, ...
    'EdgeSigma',CFG.GraphEdgeSigma,'TvEpsilon',CFG.GraphTvEpsilon, ...
    'HighConfidenceThreshold',CFG.ConfidenceThreshold,'PhysicalRange',CFG.PhysicalRange);
corrected=U(sub2ind(size(U),iz,ix));
end

function value=weighted_median(x,w)
ok=isfinite(x)&isfinite(w)&w>0; x=x(ok); w=w(ok);
if isempty(x), value=NaN; return; end
[x,order]=sort(x); w=w(order); value=x(find(cumsum(w)>=.5*sum(w),1));
end

function score=region_iou(labels,truth)
if numel(unique(truth))<2, score=mean(labels==mode(labels)); return; end
a=mean(labels==truth); flipped=labels; flipped(labels==1)=2; flipped(labels==2)=1;
score=max(a,mean(flipped==truth));
end

function assert_operational_policy(pred,meta,CFG)
for name=CFG.AllStrategies
    value=pred.(name); assert(numel(value)==height(meta)&&all(isfinite(value))&& ...
        all(value>=CFG.PhysicalRange(1)&value<=CFG.PhysicalRange(2)), ...
        'Invalid strategy output: %s',name);
end
operational_inputs=["local_m2","theory_m2","confidence","M_instability","geometry_guess"];
oracle=["true_SWS","patch_purity","material_side","distance_to_interface"];
assert(isempty(intersect(lower(operational_inputs),lower(oracle))));
end

%% Outputs and metrics

function R=attach_outputs(M2,pred,CFG)
keep={'condition_key','geometry','geometry_type','field_regime','f0','dx','dz', ...
    'map_iz','map_ix','x','z','true_SWS','confidence','patch_purity', ...
    'material_side','purity_bin','distance_bin','distance_to_interface_mm'};
R=M2(:,keep); for name=CFG.AllStrategies, R.("sws_"+name)=pred.(name); end
end

function S=summarize_predictions(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for strategy=CFG.AllStrategies
    pred=T.("sws_"+strategy); error=100*(pred-T.true_SWS)./T.true_SWS;
    if isempty(basegroups), G=ones(height(T),1); else, G=findgroups(T(:,cellstr(basegroups))); end
    for id=unique(G,'stable')'
        idx=G==id; row=table(); for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy; row.diagnostic_only=ismember(strategy,CFG.DiagnosticStrategies);
        row.N=sum(idx); row.MAPE=mean(abs(error(idx))); row.mean_signed_error=mean(error(idx));
        row.high_error_10_pct=100*mean(abs(error(idx))>10); row.high_error_20_pct=100*mean(abs(error(idx))>20);
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            j=idx&T.purity_bin==cls; token=char(cls); row.("N_"+token)=sum(j);
            row.("MAPE_"+token)=mean(abs(error(j)),'omitnan'); row.("high20_"+token)=100*mean(abs(error(j))>20,'omitnan');
        end
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function T=ring_metrics(R,CFG)
T=table(); strategies=CFG.AllStrategies; base=R(1,{'geometry','field_regime','f0','dx'});
for strategy=strategies
    pred=R.("sws_"+strategy); error=100*(pred-R.true_SWS)./R.true_SWS; row=base;
    row.strategy_name=strategy; row.diagnostic_only=ismember(strategy,CFG.DiagnosticStrategies);
    if R.geometry_type(1)=="inclusion"
        ring=R.material_side=="soft"&R.distance_to_interface_mm<=2;
        far=R.material_side=="soft"&R.distance_to_interface_mm>8;
        core=R.material_side=="hard"&R.distance_to_interface_mm>4;
        edge=R.distance_to_interface_mm<=2;
        row.N_ring=sum(ring); row.mean_soft_ring=mean(pred(ring),'omitnan');
        row.mean_soft_far=mean(pred(far),'omitnan'); row.soft_ring_drop=row.mean_soft_far-row.mean_soft_ring;
        row.mean_hard_core=mean(pred(core),'omitnan'); row.core_signed_error=mean(error(core),'omitnan');
        row.edge_MAPE=mean(abs(error(edge)),'omitnan'); row.inclusion_contrast=row.mean_hard_core-row.mean_soft_far;
    else
        row.N_ring=0; row.mean_soft_ring=NaN; row.mean_soft_far=NaN; row.soft_ring_drop=NaN;
        row.mean_hard_core=NaN; row.core_signed_error=NaN; row.edge_MAPE=NaN; row.inclusion_contrast=NaN;
    end
    T=concat_tables(T,row);
end
end

function G=make_full_gate(S,CFG)
base=S(S.strategy_name=="switch_m2",:); G=S(:,{'strategy_name','diagnostic_only'});
mixed=mean([S.MAPE_moderately_mixed S.MAPE_strongly_mixed],2,'omitnan');
base_mixed=mean([base.MAPE_moderately_mixed base.MAPE_strongly_mixed],2,'omitnan');
pure=mean([S.MAPE_pure S.MAPE_near_pure],2,'omitnan'); base_pure=mean([base.MAPE_pure base.MAPE_near_pure],2,'omitnan');
G.mixed_improvement_points=base_mixed-mixed;
G.homogeneous_delta_points=S.MAPE_homogeneous-base.MAPE_homogeneous;
G.pure_near_delta_points=pure-base_pure;
G.gate_pass=~G.diagnostic_only&G.strategy_name~="switch_m2"& ...
    G.mixed_improvement_points>=CFG.FullGateMixedImprovement& ...
    G.homogeneous_delta_points<=CFG.FullGateHomogeneousDelta& ...
    G.pure_near_delta_points<=CFG.FullGatePureDelta;
end

%% Figures

function plot_strategy_summary(T,CFG,OUT)
strategies=CFG.AllStrategies(~ismember(CFG.AllStrategies,["unknown_region_levels","userguess_region_levels"]));
M=nan(numel(strategies),3); H=M;
for i=1:numel(strategies)
    e=abs(100*(T.("sws_"+strategies(i))-T.true_SWS)./T.true_SWS);
    masks={T.purity_bin=="homogeneous",ismember(T.purity_bin,["pure","near_pure"]), ...
        ismember(T.purity_bin,["moderately_mixed","strongly_mixed"])};
    for j=1:3, M(i,j)=mean(e(masks{j}),'omitnan'); H(i,j)=100*mean(e(masks{j})>20,'omitnan'); end
end
fig=figure('Color','w','Position',[50 30 1350 760]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,M); xlabel(ax,'MAPE (%)'); title(ax,'Mean absolute error'); style(ax);
ax=nexttile(tl); barh(ax,H); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style(ax);
for ax=findall(fig,'Type','axes')', set(ax,'YTick',1:numel(strategies),'YTickLabel',pretty(strategies), ...
        'TickLabelInterpreter','none','FontSize',8); end
legend(["Homogeneous","Pure / near-pure","Mixed"],'Location','northoutside','Orientation','horizontal','Box','off');
export_fig(fig,OUT,'test30_strategy_summary.png');
end

function plot_ring_profiles(T,~,OUT)
T=T(T.geometry_type=="inclusion",:); bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
strategies=["switch_m2","theory_m2","unknown_boundary_graph", ...
    "userguess_boundary_graph_multiscale","oracle_region_graph"];
fig=figure('Color','w','Position',[70 70 1200 500]); tl=tiledlayout(1,2,'TileSpacing','compact');
for side=["soft","hard"]
    ax=nexttile(tl); hold(ax,'on');
    for strategy=strategies
        pred=T.("sws_"+strategy); y=nan(size(bins));
        for b=1:numel(bins), idx=T.material_side==side&T.distance_bin==bins(b); y(b)=mean(pred(idx),'omitnan'); end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.5,'DisplayName',pretty(strategy));
    end
    yline(ax,ternary(side=="soft",2,3),'k--','True'); set(ax,'XTick',1:numel(bins),'XTickLabel',["0-1","1-2","2-4","4-8",">8"]);
    xlabel(ax,'Distance to interface (mm)'); ylabel(ax,'Mean SWS (m/s)'); title(ax,"Inclusion "+side+" side"); style(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off'); export_fig(fig,OUT,'test30_inclusion_distance_profiles.png');
end

function plot_multiscale_instability(T,OUT)
fig=figure('Color','w','Position',[80 80 1100 450]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); scatter(ax,T.multiscale_instability,T.confidence,7,T.patch_purity,'filled','MarkerFaceAlpha',.15);
xlabel(ax,'M2/M3/M4 instability'); ylabel(ax,'Confidence'); colorbar(ax); style(ax);
ax=nexttile(tl); scatter(ax,T.monotonic_M_drop,100*(T.sws_local_m2-T.true_SWS)./T.true_SWS,7,T.patch_purity,'filled','MarkerFaceAlpha',.15);
xlabel(ax,'Monotonic SWS drop M2 to M4'); ylabel(ax,'Local M2 signed error (%)'); yline(ax,0,'k:'); colorbar(ax); style(ax);
export_fig(fig,OUT,'test30_multiscale_instability.png');
end

function plot_geometry_quality(T,OUT)
cases=["bilayer_2_3","circular_inclusion_2_3"]; U=nan(2,2);
for i=1:2, idx=T.geometry==cases(i); U(i,:)=[mean(T.unknown_region_iou(idx)) mean(T.userguess_region_iou(idx))]; end
fig=figure('Color','w','Position',[100 100 750 400]); ax=axes(fig); bar(ax,U); ylim(ax,[0 1]);
set(ax,'XTick',1:2,'XTickLabel',["Bilayer","Inclusion"]); ylabel(ax,'Material-region agreement');
legend(ax,["Unknown auto","User geometry guess"],'Location','northoutside','Orientation','horizontal','Box','off'); style(ax);
export_fig(fig,OUT,'test30_geometry_segmentation_quality.png');
end

function plot_input_M_dependence(source,OUT)
S=load(source,'T_patch'); T=S.T_patch; T=T(T.geometry=="circular_inclusion_2_3",:);
roi=abs(1e3*T.x-25)<=3&abs(1e3*T.z-25)<=3; T=T(roi,:); methods=["local_baseline","theory_discrete"];
fig=figure('Color','w','Position',[100 100 760 430]); ax=axes(fig); hold(ax,'on');
for method=methods
    y=nan(1,3); sd=y; for m=2:4, idx=T.M==m; p=T.("sws_"+method); y(m-1)=mean(p(idx)); sd(m-1)=std(p(idx)); end
    errorbar(ax,2:4,y,sd,'-o','LineWidth',1.5,'DisplayName',pretty(method));
end
yline(ax,3,'k--','True'); xlabel(ax,'REQ M'); ylabel(ax,'Inclusion-core SWS (m/s)'); legend(ax,'Box','off'); style(ax);
export_fig(fig,OUT,'test30_input_M_dependence_inclusion_core.png');
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
plot_map(T(T.condition_key==keys(1),:),CFG,fullfile(OUT.figure_dir,'test30_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable'); root=fullfile(OUT.figure_dir,'maps_by_condition');
for i=1:numel(keys), C=T(T.condition_key==keys(i),:); folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end; file=fullfile(folder,"test30__"+sanitize(keys(i))+".png");
    if exist(file,'file')~=2, plot_map(C,CFG,file); end
end
end

function plot_map(C,~,file)
C=sortrows(C,{'map_iz','map_ix'}); selected=C.sws_userguess_region_levels;
e0=abs(100*(C.sws_switch_m2-C.true_SWS)./C.true_SWS); e1=abs(100*(selected-C.true_SWS)./C.true_SWS);
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_theory_m2),map_from_rows(C,C.sws_local_m2), ...
    map_from_rows(C,C.sws_unknown_boundary_graph),map_from_rows(C,C.sws_userguess_region_levels), ...
    map_from_rows(C,C.sws_userguess_boundary_graph_multiscale), ...
    map_from_rows(C,C.sws_oracle_region_graph),map_from_rows(C,e0),map_from_rows(C,e1), ...
    map_from_rows(C,C.confidence),map_from_rows(C,C.multiscale_instability), ...
    map_from_rows(C,C.unknown_region),map_from_rows(C,C.userguess_region), ...
    map_from_rows(C,C.userguess_boundary_distance_mm)};
titles=["True","Theory M2","Local M2","Unknown graph","User region levels", ...
    "User-guess graph","Oracle-region graph","Switch error","Selected error","Confidence", ...
    "M instability","Unknown regions","User-guess regions","Boundary distance"];
fig=figure('Color','w','Visible','off','Position',[20 30 1500 850]); tl=tiledlayout(3,5,'TileSpacing','compact');
for i=1:numel(maps), ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i)); colorbar(ax); end
title(tl,C.condition_key(1),'Interpreter','none'); exportgraphics(fig,file,'Resolution',180); close(fig);
end

function print_interpretation(S,R,G,~)
oper=S(~S.diagnostic_only,:); [~,i]=min(oper.MAPE); base=S(S.strategy_name=="switch_m2",:);
selected=S(S.strategy_name=="userguess_region_levels",:);
fprintf('\nInterpretation:\n  Best operational global: %s (%.2f%%).\n',oper.strategy_name(i),oper.MAPE(i));
fprintf('  Selected homogeneous delta vs Switch M2: %+.2f points.\n',selected.MAPE_homogeneous-base.MAPE_homogeneous);
ring=R(R.strategy_name=="userguess_region_levels"&R.geometry=="circular_inclusion_2_3",:);
fprintf('  Selected mean soft-ring drop: %.3f m/s; hard-core SWS: %.3f m/s.\n', ...
    mean(ring.soft_ring_drop,'omitnan'),mean(ring.mean_hard_core,'omitnan'));
passed=G.strategy_name(G.gate_pass); if isempty(passed), fprintf('  FULL GATE: FAIL.\n');
else, fprintf('  FULL GATE: PASS (%s).\n',strjoin(passed,', ')); end
fprintf('  Output remains M2-based; M3/M4 only modify reliability.\n');
end

%% Validation and generic helpers

function validate_test30(T,G,CFG)
F=T(T.geometry==G.geometry&T.field_regime==G.field_regime&T.f0==G.f0,:);
[R,ring]=evaluate_field(F,CFG); assert(height(R)>0&&height(ring)==numel(CFG.AllStrategies));
assert(all(isfinite(R.sws_userguess_boundary_graph_multiscale)));
if ismember('M',R.Properties.VariableNames), assert(all(R.M==2)); end
fprintf('  M2 alignment, segmentation, regional calibration, and graph reconstruction passed.\n');
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end; B=B(:,A.Properties.VariableNames); T=[A;B];
end

function label=pretty(name)
label=replace(string(name),["local_m2","switch_m2","theory_m2", ...
    "unknown_theory_calibrated","unknown_region_levels","unknown_boundary_graph", ...
    "userguess_theory_calibrated","userguess_region_levels","userguess_boundary_graph", ...
    "userguess_boundary_graph_multiscale","oracle_region_graph", ...
    "local_baseline","theory_discrete"], ...
    ["Local M2","Switch M2","Theory M2","Unknown: calibrated Theory", ...
    "Unknown: region levels","Unknown: boundary graph","User guess: calibrated Theory", ...
    "User guess: region levels","User guess: boundary graph","User guess: graph + M reliability", ...
    "Oracle-region graph","Local","Theory discrete"]);
end

function style(ax), grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9; end
function export_fig(fig,OUT,name), exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',210); close(fig); end
function value=sanitize(value), value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_'); end
function value=ternary(condition,a,b), if condition, value=a; else, value=b; end, end
