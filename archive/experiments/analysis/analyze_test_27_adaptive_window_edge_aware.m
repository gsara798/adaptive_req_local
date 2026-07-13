%% analyze_test_27_adaptive_window_edge_aware.m
% Test 27: confidence-gated adaptive window and edge-aware correction.
%
% No model is trained here. Frozen Test 26 prediction cubes and Test 22
% detector scores are reused. Oracle quantities are attached only after all
% operational strategies are formed, except oracle_region_reference, which
% is explicitly marked diagnostic_only.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST27_MODE          = quick | full
%   ADAPTIVE_REQ_TEST27_QUICK_MODE    = true | false (legacy)
%   ADAPTIVE_REQ_TEST27_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST27_SAVE_ALL_MAPS = true | false

clear; clc; close all;
format compact;

%% Configuration

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();

CORR=struct();
mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST27_MODE'))));
if mode=="quick", CORR.QuickMode=true;
elseif mode=="full", CORR.QuickMode=false;
elseif mode=="", CORR.QuickMode=env_true('ADAPTIVE_REQ_TEST27_QUICK_MODE',true);
else, error('ADAPTIVE_REQ_TEST27_MODE must be quick or full.');
end
CORR.RunMode=ternary(CORR.QuickMode,"quick","full");
CORR.ValidateOnly=env_true('ADAPTIVE_REQ_TEST27_VALIDATE_ONLY',false);
CORR.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST27_SAVE_ALL_MAPS',true);
CORR.TargetStepM=0.2e-3;
CORR.Dx=0.2e-3;
CORR.ConfidenceThresholds=[0.5 0.7 0.8];
CORR.PrimaryDetector="ML_bagged_trees";
CORR.BaseModels=["HybridLocalGlobal_T18_noUserRegime","LocalOnly_T18"];
CORR.OptionalModel="HybridLocalGlobal_T18_withUserRegimeGuess";
CORR.SwitchThresholds=[0.7 0.8];
CORR.BlendRanges=[0.5 0.8;0.6 0.9];
CORR.PhysicalSwsRange=[0.5 10];
CORR.SmallWindowWidthRatioMax=1.10;
CORR.SmallWindowTailDeltaMax=0.05;
CORR.EdgeRadii=[2 4 6];
CORR.EdgeMinimumNeighbors=3;
CORR.PseudoRegionMinimumSeparation=0.25;
CORR.RandomSeed=27001;
CORR.SelectedStrategy="small_window_then_edge_aware_c080";
CORR.CompactStrategies=["hybrid_baseline","local_baseline", ...
    "confidence_switch_local_hybrid_c080", ...
    "confidence_blend_local_hybrid_050_080", ...
    "small_window_low_confidence_c080", ...
    "small_window_if_operationally_better_c080", ...
    "edge_aware_interpolation_low_confidence_c080", ...
    CORR.SelectedStrategy,"oracle_region_reference"];
if CORR.QuickMode
    CORR.Frequencies=500; CORR.M=[2 3];
    CORR.Regimes=["directional_2D","diffuse_3D"];
else
    CORR.Frequencies=[300 400 500 600]; CORR.M=[2 3 4];
    CORR.Regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
end
CORR.Geometries=["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"];

SOURCE=fullfile(root_dir,'outputs','test_26_confidence_gated_corrections');
if CORR.QuickMode, SOURCE=fullfile(SOURCE,'quick'); end
OUT=make_output_dirs(root_dir,CORR.QuickMode);
write_config_json(CORR,fullfile(OUT.root_dir,'test27_configuration.json'));

fprintf('\nTest 27: confidence-gated adaptive window and edge-aware correction\n');
fprintf('No training. Operational correction inputs exclude oracle variables.\n');
fprintf('Run mode: %s | source: %s\n',CORR.RunMode,SOURCE);

%% Discover frozen Test 26 conditions

source_condition_dir=fullfile(SOURCE,'data','condition_checkpoints');
source_spectral_dir=fullfile(SOURCE,'data','spectral_checkpoints');
assert(exist(source_condition_dir,'dir')==7, ...
    'Test 26 checkpoints missing: %s',source_condition_dir);
FILES=discover_conditions(source_condition_dir,CORR);
assert(~isempty(FILES),'No Test 26 checkpoints match the Test 27 design.');
fprintf('Matched %d Test 26 conditions.\n',numel(FILES));

if CORR.ValidateOnly
    validate_test27(FILES,source_spectral_dir,CORR);
    fprintf('Test 27 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Evaluate strategies condition by condition

PATCH_PARTS=cell(numel(FILES),1); STAT_PARTS=cell(numel(FILES),1);
LOCAL_PARTS=cell(numel(FILES),1); SMALL_PARTS=cell(numel(FILES),1);
EDGE_PARTS=cell(numel(FILES),1);
for ci=1:numel(FILES)
    timer=tic; F=FILES(ci);
    checkpoint=fullfile(OUT.condition_dir,"test27__"+sanitize(F.key)+".mat");
    if exist(checkpoint,'file')==2
        S=load(checkpoint,'RESULT');
        if isfield(S,'RESULT')&&result_is_compatible(S.RESULT,CORR)
            R=S.RESULT; PATCH_PARTS{ci}=R.T_compact; STAT_PARTS{ci}=R.T_stats;
            LOCAL_PARTS{ci}=R.T_local; SMALL_PARTS{ci}=R.T_small;
            EDGE_PARTS{ci}=R.T_edge;
            fprintf('[%d/%d] Reused %s (%d patches).\n', ...
                ci,numel(FILES),F.key,height(R.T_compact));
            continue;
        end
    end

    base=load_test26_condition(F.path,CORR);
    spectral=load_spectral_features(source_spectral_dir,F.key,base.meta);
    small=load_small_window(F,FILES,CORR,source_spectral_dir);
    context=build_operational_context(base,small,spectral,CORR);
    [cube,specs,diagnostics]=evaluate_strategies(base,context,CORR);
    assert_operational_outputs(cube,specs,base.meta,CORR);

    T_stats=condition_statistics(base.meta,cube,specs,base.hybrid,CORR);
    T_compact=make_compact_table(base.meta,cube,specs,diagnostics,CORR);
    T_local=local_vs_hybrid_rows(base.meta,base.local,base.hybrid);
    T_small=small_window_rows(base.meta,context,diagnostics,CORR);
    T_edge=edge_aware_rows(base.meta,diagnostics,CORR);
    RESULT=struct('T_compact',T_compact,'prediction_cube',single(cube), ...
        'strategy_specs',specs,'diagnostics',diagnostics,'T_stats',T_stats, ...
        'T_local',T_local,'T_small',T_small,'T_edge',T_edge, ...
        'config_signature',config_signature(CORR));
    save(checkpoint,'RESULT','-v7.3');
    PATCH_PARTS{ci}=T_compact; STAT_PARTS{ci}=T_stats;
    LOCAL_PARTS{ci}=T_local; SMALL_PARTS{ci}=T_small; EDGE_PARTS{ci}=T_edge;
    fprintf('[%d/%d] %s: %d patches, %d strategies in %.1f s.\n', ...
        ci,numel(FILES),F.key,height(base.meta),height(specs),toc(timer));
end

T_patch=vertcat(PATCH_PARTS{:}); T_condition_stats=vertcat(STAT_PARTS{:});
T_local_condition=vertcat(LOCAL_PARTS{:}); T_small_condition=vertcat(SMALL_PARTS{:});
T_edge_condition=vertcat(EDGE_PARTS{:});
clear PATCH_PARTS STAT_PARTS LOCAL_PARTS SMALL_PARTS EDGE_PARTS;

%% Aggregate metrics and save tables

T_overall=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only"],"overall");
T_geometry=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only","geometry"],"overall");
T_regime=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only","field_regime"],"overall");
T_frequency_M=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only","f0","M"],"overall");
T_purity=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only","geometry","purity_bin"],"purity");
T_distance=aggregate_stats(T_condition_stats, ...
    ["strategy_name","diagnostic_only","geometry","distance_bin"],"distance");
T_local=aggregate_local_rows(T_local_condition);
T_small=aggregate_simple_rows(T_small_condition, ...
    ["geometry","field_regime","f0","M","confidence_threshold"]);
T_edge=aggregate_simple_rows(T_edge_condition, ...
    ["geometry","field_regime","f0","M","confidence_threshold"]);
T_best=best_strategy_candidates(T_overall,T_regime,T_frequency_M);

writetable(T_patch,fullfile(OUT.table_dir,'test27_patch_level_results.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test27_strategy_summary_overall.csv'));
writetable(T_geometry,fullfile(OUT.table_dir,'test27_strategy_summary_by_geometry.csv'));
writetable(T_regime,fullfile(OUT.table_dir,'test27_strategy_summary_by_regime.csv'));
writetable(T_frequency_M,fullfile(OUT.table_dir,'test27_strategy_summary_by_frequency_M.csv'));
writetable(T_purity,fullfile(OUT.table_dir,'test27_strategy_summary_by_purity_bin.csv'));
writetable(T_distance,fullfile(OUT.table_dir,'test27_strategy_summary_by_distance_bin.csv'));
writetable(T_local,fullfile(OUT.table_dir,'test27_local_vs_hybrid_summary.csv'));
writetable(T_small,fullfile(OUT.table_dir,'test27_small_window_acceptance_summary.csv'));
writetable(T_edge,fullfile(OUT.table_dir,'test27_edge_aware_summary.csv'));
writetable(T_best,fullfile(OUT.table_dir,'test27_best_strategy_candidates.csv'));
save(fullfile(OUT.data_dir,'test27_compact_results.mat'), ...
    'T_patch','T_condition_stats','CORR','-v7.3');

%% Figures and interpretation

plot_strategy_summary(T_patch,CORR,OUT);
plot_error_vs_distance(T_patch,CORR,OUT);
plot_error_vs_purity(T_patch,CORR,OUT);
plot_error_scatter(T_patch,CORR,OUT);
plot_local_vs_hybrid(T_patch,CORR,OUT);
plot_small_window_diagnostics(T_patch,CORR,OUT);
plot_edge_diagnostics(T_patch,CORR,OUT);
plot_representative_maps(T_patch,CORR,OUT);
if CORR.SaveAllMaps, plot_all_maps(T_patch,CORR,OUT); end
print_interpretation(T_overall,T_local,T_small,T_edge,CORR);
fprintf('\nTables: %s\nFigures: %s\nTest 27 complete.\n', ...
    OUT.table_dir,OUT.figure_dir);

%% Setup and source helpers

function OUT=make_output_dirs(root_dir,quick)
OUT.root_dir=fullfile(root_dir,'outputs','test_27_adaptive_window_edge_aware');
if quick, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables');
OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data');
OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
for d=string(struct2cell(OUT))'
    if exist(d,'dir')~=7, mkdir(d); end
end
end

function tf=env_true(name,default_value)
value=strtrim(getenv(name));
if isempty(value), tf=default_value; else
    tf=any(strcmpi(value,{'true','1','yes','on'}));
end
end

function write_config_json(CORR,file)
C=CORR;
for f=string(fieldnames(C))'
    if isstring(C.(f)), C.(f)=cellstr(C.(f)); end
end
fid=fopen(file,'w'); assert(fid>=0,'Cannot write %s.',file);
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char');
end

function FILES=discover_conditions(folder,CORR)
d=dir(fullfile(folder,'corrected__*.mat'));
FILES=struct('path',{},'key',{},'M',{});
for i=1:numel(d)
    name=string(d(i).name); key=erase(name,["corrected__" ".mat"]); low=lower(key);
    keep_f=any(arrayfun(@(f)contains(low,"__f"+string(f)+"__"),CORR.Frequencies));
    keep_m=any(arrayfun(@(m)contains(low,"__m"+string(m)+"__"),CORR.M));
    keep_dx=contains(low,"__dx"+string(round(1e6*CORR.Dx))+"um__");
    keep_g=any(arrayfun(@(g)startsWith(low,lower(g)+"__"),CORR.Geometries));
    keep_r=false;
    for r=CORR.Regimes, keep_r=keep_r||contains(low,"__"+lower(r)+"__"); end
    if keep_f&&keep_m&&keep_dx&&keep_g&&keep_r
        token=regexp(key,'__M(\d+)__','tokens','once');
        FILES(end+1)=struct('path',fullfile(d(i).folder,d(i).name), ...
            'key',key,'M',str2double(token{1})); %#ok<AGROW>
    end
end
end

function base=load_test26_condition(file,CORR)
S=load(file,'RESULT'); assert(isfield(S,'RESULT'),'Invalid Test 26 checkpoint: %s',file);
R=S.RESULT; specs=string(R.strategy_specs.strategy_name);
si=find(specs=="baseline",1); assert(~isempty(si),'Test 26 baseline missing.');
models=string(R.models); detector=string(R.detectors);
di=find(detector==CORR.PrimaryDetector,1); assert(~isempty(di),'Bagged detector missing.');
hi=find(models==CORR.BaseModels(1),1); li=find(models==CORR.BaseModels(2),1);
assert(~isempty(hi)&&~isempty(li),'Hybrid or LocalOnly missing in Test 26 cube.');
meta=R.T_compact; names=string(meta.Properties.VariableNames);
drop=startsWith(names,"sws_")|names=="sws_model"|names=="detector";
meta(:,drop)=[];
base=struct('meta',meta,'hybrid',double(R.prediction_cube(:,si,hi,di)), ...
    'local',double(R.prediction_cube(:,si,li,di)),'source_result',R);
oi=find(models==CORR.OptionalModel,1);
if isempty(oi), base.hybrid_with_guess=nan(height(meta),1);
else, base.hybrid_with_guess=double(R.prediction_cube(:,si,oi,di)); end
end

function small=load_small_window(F,FILES,CORR,spectral_dir)
small=struct('available',false,'meta',table(),'hybrid',[],'local',[], ...
    'width',[],'tail',[],'target_M',NaN);
if F.M<=2, return; end
target=F.M-1; small_key=regexprep(F.key,'__M\d+__',"__M"+string(target)+"__");
j=find(string({FILES.key})==small_key,1);
if isempty(j), return; end
small_base=load_test26_condition(FILES(j).path,CORR);
sf=load_spectral_features(spectral_dir,small_key,small_base.meta);
small=struct('available',true,'meta',small_base.meta, ...
    'hybrid',small_base.hybrid,'local',small_base.local, ...
    'width',sf.width,'tail',sf.tail,'target_M',target);
end

function F=load_spectral_features(folder,key,meta)
file=fullfile(folder,"spectra__"+sanitize(key)+".mat");
F=struct('width',nan(height(meta),1),'tail',nan(height(meta),1));
if exist(file,'file')~=2
    warning('Test27:SpectralCacheMissing','Spectral cache unavailable: %s',file); return;
end
S=load(file,'SPEC'); SPEC=S.SPEC;
if numel(SPEC.curves)~=height(meta)
    warning('Test27:SpectralCacheMismatch','Ignoring mismatched spectral cache %s.',file); return;
end
for i=1:height(meta)
    [F.width(i),F.tail(i)]=spectral_shape(SPEC.curves{i});
end
end

function [width,tail]=spectral_shape(curve)
k=double(curve.k_cent(:)); s=max(double(curve.Srad(:)),0);
s(~isfinite(s))=0; E=cumsum(s)/(sum(s)+eps);
k10=interp_quantile(k,E,.1); k50=interp_quantile(k,E,.5); k90=interp_quantile(k,E,.9);
width=(k90-k10)/max(k50,eps); tail=sum(s(k>1.5*k50))/(sum(s)+eps);
end

function value=interp_quantile(k,E,q)
[Eu,ia]=unique(E,'stable'); ku=k(ia);
if numel(Eu)<2, value=NaN; else, value=interp1(Eu,ku,q,'linear','extrap'); end
end

%% Operational strategy construction

function C=build_operational_context(base,small,spectral,CORR)
n=height(base.meta); C=struct();
C.base_width=spectral.width; C.base_tail=spectral.tail;
C.base_disagreement=abs(base.local-base.hybrid);
C.pseudo_region=pseudo_segment(base.local,base.hybrid,base.meta,CORR);
C.small_available=false(n,1); C.small_local=nan(n,1); C.small_hybrid=nan(n,1);
C.small_width=nan(n,1); C.small_tail=nan(n,1); C.small_M=nan(n,1);
if small.available
    C.small_available(:)=true; C.small_M(:)=small.target_M;
    C.small_local=interpolate_to_grid(small.meta,small.local,base.meta);
    C.small_hybrid=interpolate_to_grid(small.meta,small.hybrid,base.meta);
    C.small_width=interpolate_to_grid(small.meta,small.width,base.meta);
    C.small_tail=interpolate_to_grid(small.meta,small.tail,base.meta);
end
valid=C.small_available&isfinite(C.small_local)& ...
    C.small_local>=CORR.PhysicalSwsRange(1)&C.small_local<=CORR.PhysicalSwsRange(2);
width_ok=C.small_width<=CORR.SmallWindowWidthRatioMax*C.base_width;
tail_ok=C.small_tail<=C.base_tail+CORR.SmallWindowTailDeltaMax;
disagreement_ok=abs(C.small_local-C.small_hybrid)<=C.base_disagreement;
spectral_missing=~isfinite(C.small_width)|~isfinite(C.base_width)| ...
    ~isfinite(C.small_tail)|~isfinite(C.base_tail);
C.small_operational_accept=valid&disagreement_ok& ...
    ((width_ok&tail_ok)|spectral_missing);
end

function values=interpolate_to_grid(source_meta,source_values,target_meta)
ok=isfinite(source_values);
if nnz(ok)<3, values=nan(height(target_meta),1); return; end
F=scatteredInterpolant(source_meta.x(ok),source_meta.z(ok),source_values(ok), ...
    'linear','nearest');
values=F(target_meta.x,target_meta.z);
end

function labels=pseudo_segment(local,hybrid,meta,CORR)
seed=.7*local+.3*hybrid; [A,iz,ix]=map_from_rows(meta,seed);
A=medfilt2(A,[3 3],'symmetric'); x=A(isfinite(A));
if isempty(x), labels=ones(size(seed)); return; end
c=[prctile(x,25) prctile(x,75)];
for iter=1:15
    group1=abs(x-c(1))<=abs(x-c(2));
    if any(group1), c(1)=median(x(group1)); end
    if any(~group1), c(2)=median(x(~group1)); end
end
if abs(diff(sort(c)))<CORR.PseudoRegionMinimumSeparation
    L=ones(size(A));
else
    [c,~]=sort(c); L=ones(size(A)); L(abs(A-c(2))<abs(A-c(1)))=2;
end
labels=L(sub2ind(size(L),iz,ix));
end

function [cube,specs,D]=evaluate_strategies(base,C,CORR)
n=height(base.meta); cube=zeros(n,0); specs=table(); D=struct();
[cube,specs]=add_strategy(cube,specs,base.hybrid,"hybrid_baseline",NaN,false,"baseline");
[cube,specs]=add_strategy(cube,specs,base.local,"local_baseline",NaN,false,"baseline");
for threshold=CORR.SwitchThresholds
    pred=base.hybrid; pred(base.meta.confidence>=threshold)=base.local(base.meta.confidence>=threshold);
    name="confidence_switch_local_hybrid_c"+sprintf('%03d',round(100*threshold));
    [cube,specs]=add_strategy(cube,specs,pred,name,threshold,false,"switch");
end
for i=1:size(CORR.BlendRanges,1)
    lo=CORR.BlendRanges(i,1); hi=CORR.BlendRanges(i,2);
    w=min(max((base.meta.confidence-lo)/(hi-lo),0),1);
    pred=w.*base.local+(1-w).*base.hybrid;
    name="confidence_blend_local_hybrid_"+sprintf('%03d',round(100*lo))+"_"+sprintf('%03d',round(100*hi));
    [cube,specs]=add_strategy(cube,specs,pred,name,hi,false,"blend");
    if lo==.5&&hi==.8, D.w_local=w; end
end

D.small_source=zeros(n,numel(CORR.ConfidenceThresholds),'uint8');
D.edge_source=zeros(n,numel(CORR.ConfidenceThresholds),'uint8');
D.combined_source=zeros(n,numel(CORR.ConfidenceThresholds),'uint8');
for ti=1:numel(CORR.ConfidenceThresholds)
    threshold=CORR.ConfidenceThresholds(ti); low=base.meta.confidence<threshold;
    tag="c"+sprintf('%03d',round(100*threshold));
    pred=base.hybrid; use=low&C.small_available&isfinite(C.small_local);
    pred(use)=C.small_local(use); D.small_source(use,ti)=1;
    [cube,specs]=add_strategy(cube,specs,pred, ...
        "small_window_low_confidence_"+tag,threshold,false,"small_window");

    pred_better=base.hybrid; accept=low&C.small_operational_accept;
    pred_better(accept)=C.small_local(accept);
    [cube,specs]=add_strategy(cube,specs,pred_better, ...
        "small_window_if_operationally_better_"+tag,threshold,false,"small_operational");

    [edge,edge_used]=edge_aware_correct(base.hybrid,base.local, ...
        base.meta.confidence,C.pseudo_region,base.meta,threshold,CORR);
    D.edge_source(edge_used,ti)=2;
    [cube,specs]=add_strategy(cube,specs,edge, ...
        "edge_aware_interpolation_low_confidence_"+tag,threshold,false,"edge_aware");

    combined=base.hybrid; combined(accept)=C.small_local(accept);
    unresolved=low&~accept;
    combined(unresolved&edge_used)=edge(unresolved&edge_used);
    D.combined_source(accept,ti)=1;
    D.combined_source(unresolved&edge_used,ti)=2;
    D.combined_source(unresolved&~edge_used,ti)=3;
    [cube,specs]=add_strategy(cube,specs,combined, ...
        "small_window_then_edge_aware_"+tag,threshold,false,"combined");
end

ti=find(CORR.ConfidenceThresholds==.8,1); operational=cube(:,specs.strategy_name==CORR.SelectedStrategy);
oracle=operational; pure=base.meta.patch_purity>=.90;
oracle(pure)=base.local(pure);
[cube,specs]=add_strategy(cube,specs,oracle,"oracle_region_reference",.8,true,"oracle");
D.pseudo_region=C.pseudo_region; D.small_operational_accept=C.small_operational_accept;
D.small_available=C.small_available; D.small_local=C.small_local;
D.small_width=C.small_width; D.base_width=C.base_width; D.small_tail=C.small_tail;
D.base_tail=C.base_tail; D.small_M=C.small_M;
D.selected_source=D.combined_source(:,ti);
D.reliability=base.meta.confidence;
end

function [cube,specs]=add_strategy(cube,specs,pred,name,threshold,diagnostic,family)
cube(:,end+1)=pred;
row=table(string(name),threshold,logical(diagnostic),string(family), ...
    'VariableNames',{'strategy_name','confidence_threshold','diagnostic_only','family'});
specs=[specs;row];
end

function [corrected,used]=edge_aware_correct(baseline,source,confidence,region,meta,threshold,CORR)
[B,iz,ix]=map_from_rows(meta,baseline); S=map_from_rows(meta,source);
C=map_from_rows(meta,confidence); R=map_from_rows(meta,region);
corrected_map=B; used_map=false(size(B)); low=C<threshold;
[rows,cols]=find(low);
for p=1:numel(rows)
    z=rows(p); x=cols(p); found=false;
    for radius=CORR.EdgeRadii
        zr=max(1,z-radius):min(size(B,1),z+radius);
        xr=max(1,x-radius):min(size(B,2),x+radius);
        [XX,ZZ]=meshgrid(xr,zr);
        valid=C(zr,xr)>=threshold&R(zr,xr)==R(z,x)&isfinite(S(zr,xr));
        if nnz(valid)>=CORR.EdgeMinimumNeighbors
            values=S(zr,xr); distances=hypot(XX-x,ZZ-z);
            corrected_map(z,x)=weighted_median(values(valid),1./(distances(valid)+.5));
            used_map(z,x)=true; found=true; break;
        end
    end
    if ~found, corrected_map(z,x)=B(z,x); end
end
corrected=corrected_map(sub2ind(size(B),iz,ix));
used=used_map(sub2ind(size(B),iz,ix));
end

function value=weighted_median(x,w)
[x,order]=sort(x(:)); w=w(order); value=x(find(cumsum(w)>=.5*sum(w),1));
end

function assert_operational_outputs(cube,specs,meta,CORR)
assert(size(cube,1)==height(meta)&&size(cube,2)==height(specs));
assert(all(isfinite(cube(:))&cube(:)>=CORR.PhysicalSwsRange(1)& ...
    cube(:)<=CORR.PhysicalSwsRange(2)),'Operational strategy produced invalid SWS.');
banned=["true_SWS","patch_purity","material_side","distance_to_interface", ...
    "q_true","q_theory"];
operational=specs(~specs.diagnostic_only,:);
assert(~any(contains(lower(operational.strategy_name),lower(banned))));
end

%% Compact outputs and statistics

function T=make_compact_table(meta,cube,specs,D,CORR)
T=meta;
for name=CORR.CompactStrategies
    i=find(specs.strategy_name==name,1);
    if ~isempty(i), T.(matlab.lang.makeValidName("sws_"+name))=cube(:,i); end
end
T.pseudo_region=D.pseudo_region; T.reliability=D.reliability;
T.small_window_available=D.small_available;
T.small_window_accepted=D.small_operational_accept;
T.small_window_M=D.small_M; T.small_window_sws=D.small_local;
T.small_window_spectral_width=D.small_width;
T.base_spectral_width=D.base_width;
T.small_window_highk_tail=D.small_tail; T.base_highk_tail=D.base_tail;
T.correction_source_c080=D.selected_source;
T.w_local_blend_050_080=D.w_local;
end

function T=condition_statistics(meta,cube,specs,baseline,~)
T=table();
for si=1:height(specs)
    pred=cube(:,si); R=make_eval_rows(meta,pred,baseline,specs(si,:));
    T=concat_tables(T,raw_stats(R,"strategy_name","overall"));
    T=concat_tables(T,raw_stats(R,["strategy_name","purity_bin"],"purity"));
    T=concat_tables(T,raw_stats(R,["strategy_name","distance_bin"],"distance"));
end
for name=["geometry","field_regime","f0","M","dx"]
    T.(name)=repmat(meta.(name)(1),height(T),1);
end
end

function R=make_eval_rows(meta,pred,baseline,spec)
R=meta; R.strategy_name=repmat(spec.strategy_name,height(meta),1);
R.diagnostic_only=repmat(spec.diagnostic_only,height(meta),1);
R.confidence_threshold=repmat(spec.confidence_threshold,height(meta),1);
R.error=100*(pred-R.true_SWS)./R.true_SWS; R.abs_error=abs(R.error);
R.high10=R.abs_error>10; R.high20=R.abs_error>20;
R.under=R.error<0; R.over=R.error>0;
R.modified=abs(pred-baseline)>1e-9;
base_error=abs(100*(baseline-R.true_SWS)./R.true_SWS);
R.improvement=base_error-R.abs_error;
if isfinite(spec.confidence_threshold), R.high_confidence=R.confidence>=spec.confidence_threshold;
else, R.high_confidence=true(height(R),1); end
end

function S=raw_stats(T,groups,scope)
G=findgroups(T(:,cellstr(groups))); ids=unique(G,'stable'); S=table();
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups)); row.stat_scope=string(scope);
    row.diagnostic_only=X.diagnostic_only(1); row.confidence_threshold=X.confidence_threshold(1);
    row.N=height(X); row.sum_abs=sum(X.abs_error); row.sum_signed=sum(X.error);
    row.condition_median_signed=median(X.error); row.under_count=sum(X.under);
    row.over_count=sum(X.over); row.high10_count=sum(X.high10); row.high20_count=sum(X.high20);
    row.modified_count=sum(X.modified); row.low_conf_count=sum(~X.high_confidence);
    row.modified_improvement_sum=sum(X.improvement(X.modified),'omitnan');
    row.high_conf_degradation_sum=sum(-X.improvement(X.high_confidence),'omitnan');
    row.high_conf_count=sum(X.high_confidence);
    for class=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
        idx=X.purity_bin==class; token=char(class);
        row.("N_"+token)=sum(idx); row.("sum_abs_"+token)=sum(X.abs_error(idx));
    end
    for side=["soft","hard"]
        idx=X.material_side==side; token=char(side);
        row.("N_"+token)=sum(idx); row.("sum_abs_"+token)=sum(X.abs_error(idx));
        row.("high20_"+token)=sum(X.high20(idx));
    end
    S=concat_tables(S,row);
end
end

function S=aggregate_stats(T,groups,scope)
T=T(T.stat_scope==scope,:); G=findgroups(T(:,cellstr(groups))); ids=unique(G,'stable'); S=table();
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups)); row.N=sum(X.N);
    row.MAPE=sum(X.sum_abs)/row.N; row.mean_signed_error=sum(X.sum_signed)/row.N;
    row.median_signed_error=weighted_median(X.condition_median_signed,X.N);
    row.underestimate_pct=100*sum(X.under_count)/row.N;
    row.overestimate_pct=100*sum(X.over_count)/row.N;
    row.high_error_10_pct=100*sum(X.high10_count)/row.N;
    row.high_error_20_pct=100*sum(X.high20_count)/row.N;
    row.modified_pct=100*sum(X.modified_count)/row.N;
    row.low_confidence_pct=100*sum(X.low_conf_count)/row.N;
    row.mean_improvement_modified=sum(X.modified_improvement_sum)/max(sum(X.modified_count),1);
    row.degradation_high_confidence=sum(X.high_conf_degradation_sum)/max(sum(X.high_conf_count),1);
    row.confidence_threshold=X.confidence_threshold(1);
    for class=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
        token=char(class); count=sum(X.("N_"+token));
        row.("MAPE_"+token)=sum(X.("sum_abs_"+token))/max(count,1);
        if count==0, row.("MAPE_"+token)=NaN; end
    end
    for side=["soft","hard"]
        token=char(side); count=sum(X.("N_"+token));
        row.("MAPE_"+token)=sum(X.("sum_abs_"+token))/max(count,1);
        row.("high_error_20_"+token+"_pct")=100*sum(X.("high20_"+token))/max(count,1);
        if count==0, row.("MAPE_"+token)=NaN; end
    end
    S=concat_tables(S,row);
end
end

function T=local_vs_hybrid_rows(meta,local,hybrid)
eL=abs(100*(local-meta.true_SWS)./meta.true_SWS);
eH=abs(100*(hybrid-meta.true_SWS)./meta.true_SWS);
T=meta(1,{'geometry','field_regime','f0','M','dx'});
T.N=height(meta); T.sum_abs_local=sum(eL); T.sum_abs_hybrid=sum(eH);
T.local_win_count=sum(eL<eH); T.disagreement_sum=sum(abs(local-hybrid));
for class=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
    idx=meta.purity_bin==class; token=char(class); T.("N_"+token)=sum(idx);
    T.("sum_local_"+token)=sum(eL(idx)); T.("sum_hybrid_"+token)=sum(eH(idx));
end
end

function S=aggregate_local_rows(T)
groups=["geometry","field_regime","f0","M"]; G=findgroups(T(:,cellstr(groups)));
ids=unique(G,'stable'); S=table();
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups)); row.N=sum(X.N);
    row.MAPE_local=sum(X.sum_abs_local)/row.N; row.MAPE_hybrid=sum(X.sum_abs_hybrid)/row.N;
    row.local_wins_pct=100*sum(X.local_win_count)/row.N;
    row.mean_local_hybrid_disagreement=sum(X.disagreement_sum)/row.N;
    for class=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
        token=char(class); n=sum(X.("N_"+token));
        row.("MAPE_local_"+token)=sum(X.("sum_local_"+token))/max(n,1);
        row.("MAPE_hybrid_"+token)=sum(X.("sum_hybrid_"+token))/max(n,1);
    end
    S=concat_tables(S,row);
end
end

function T=small_window_rows(meta,C,~,CORR)
T=table();
for ti=1:numel(CORR.ConfidenceThresholds)
    th=CORR.ConfidenceThresholds(ti); low=meta.confidence<th;
    row=meta(1,{'geometry','field_regime','f0','M','dx'});
    row.confidence_threshold=th; row.N=height(meta); row.low_count=sum(low);
    row.available_count=sum(low&C.small_available);
    row.accepted_count=sum(low&C.small_operational_accept);
    row.mean_base_width=mean(C.base_width(low),'omitnan');
    row.mean_small_width=mean(C.small_width(low),'omitnan');
    row.mean_base_tail=mean(C.base_tail(low),'omitnan');
    row.mean_small_tail=mean(C.small_tail(low),'omitnan');
    T=concat_tables(T,row);
end
end

function T=edge_aware_rows(meta,D,CORR)
T=table();
for ti=1:numel(CORR.ConfidenceThresholds)
    row=meta(1,{'geometry','field_regime','f0','M','dx'}); source=D.edge_source(:,ti);
    row.confidence_threshold=CORR.ConfidenceThresholds(ti); row.N=height(meta);
    row.low_count=sum(meta.confidence<CORR.ConfidenceThresholds(ti));
    row.edge_corrected_count=sum(source==2); row.fallback_count=row.low_count-row.edge_corrected_count;
    row.pseudo_region_2_pct=100*mean(D.pseudo_region==2);
    T=concat_tables(T,row);
end
end

function S=aggregate_simple_rows(T,groups)
G=findgroups(T(:,cellstr(groups))); ids=unique(G,'stable'); S=table();
numeric=setdiff(string(T.Properties.VariableNames),groups);
for i=1:numel(ids)
    X=T(G==ids(i),:); row=X(1,cellstr(groups));
    for name=numeric, row.(name)=sum(X.(name),'omitnan'); end
    if ismember('low_count',row.Properties.VariableNames)
        row.accepted_low_pct=100*getfield_default(row,'accepted_count',0)/max(row.low_count,1);
        row.edge_corrected_low_pct=100*getfield_default(row,'edge_corrected_count',0)/max(row.low_count,1);
    end
    S=concat_tables(S,row);
end
end

function value=getfield_default(S,name,default)
if ismember(name,S.Properties.VariableNames), value=S.(name); else, value=default; end
end

function B=best_strategy_candidates(overall,regime,frequency_M)
B=table(); B=concat_tables(B,best_rows(overall,strings(0,1),"overall"));
B=concat_tables(B,best_rows(regime,"field_regime","regime"));
B=concat_tables(B,best_rows(frequency_M,["f0","M"],"frequency_M"));
end

function B=best_rows(S,groups,scope)
S=S(~S.diagnostic_only,:); keys=string(groups); if isempty(keys), G=ones(height(S),1); else, G=findgroups(S(:,cellstr(keys))); end
B=table();
for id=unique(G,'stable')'
    X=S(G==id,:); [~,ig]=min(X.MAPE);
    mixed=mean([X.MAPE_moderately_mixed X.MAPE_strongly_mixed],2,'omitnan'); [~,im]=min(mixed);
    pure=mean([X.MAPE_pure X.MAPE_near_pure],2,'omitnan'); [~,ip]=min(pure);
    baseline=X(X.strategy_name=="hybrid_baseline",:);
    score=(X.MAPE-baseline.MAPE)+2*max(X.MAPE_homogeneous-baseline.MAPE_homogeneous-.75,0)+ ...
        max(pure-mean([baseline.MAPE_pure baseline.MAPE_near_pure],'omitnan'),0);
    [~,it]=min(score); row=X(1,cellstr(keys)); row.candidate_scope=string(scope);
    row.best_global_strategy=X.strategy_name(ig); row.best_mixed_strategy=X.strategy_name(im);
    row.best_pure_strategy=X.strategy_name(ip); row.best_tradeoff_strategy=X.strategy_name(it);
    row.best_tradeoff_score=score(it); B=concat_tables(B,row);
end
end

%% Figures

function plot_strategy_summary(T,CORR,OUT)
strategies=plot_strategies(CORR);
M=nan(numel(strategies),3); H=M;
for si=1:numel(strategies)
    for ri=1:3
        idx=T.purity_bin=="homogeneous";
        if ri==2, idx=ismember(T.purity_bin,["pure","near_pure"]);
        elseif ri==3, idx=ismember(T.purity_bin,["moderately_mixed","strongly_mixed"]); end
        pred=T.(pred_var(strategies(si))); err=abs(100*(pred-T.true_SWS)./T.true_SWS);
        M(si,ri)=mean(err(idx),'omitnan'); H(si,ri)=100*mean(err(idx)>20,'omitnan');
    end
end
fig=figure('Color','w','Position',[100 100 1250 470]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); bar(ax,M); title(ax,'MAPE by region'); ylabel(ax,'MAPE (%)'); style_axes(ax);
ax=nexttile(tl); bar(ax,H); title(ax,'Large-error rate'); ylabel(ax,'Error >20% (%)'); style_axes(ax);
for ax=findall(fig,'Type','axes')', set(ax,'XTick',1:numel(strategies),'XTickLabel',pretty_strategy(strategies),'XTickLabelRotation',18,'TickLabelInterpreter','none'); end
legend(["Homogeneous","Pure / near-pure","Mixed"],'Location','northoutside','Orientation','horizontal','Box','off');
export_fig(fig,OUT,'test27_strategy_summary.png');
end

function plot_error_vs_distance(T,CORR,OUT)
plot_binned_lines(T,CORR,OUT,"distance_bin",["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"], ...
    ["0–1","1–2","2–4","4–8",">8"],'Distance to interface (mm)','test27_error_vs_distance.png');
end

function plot_error_vs_purity(T,CORR,OUT)
plot_binned_lines(T,CORR,OUT,"purity_bin",["strongly_mixed","moderately_mixed","near_pure","pure"], ...
    ["Strongly mixed","Moderately mixed","Near-pure","Pure"],'Patch class','test27_error_vs_purity.png');
end

function plot_binned_lines(T,CORR,OUT,binvar,bins,labels,xlabel_text,file)
T=T(T.geometry_type~="homogeneous",:); strategies=plot_strategies(CORR);
fig=figure('Color','w','Position',[100 100 1150 480]); tl=tiledlayout(1,2,'TileSpacing','compact');
for gi=1:2
    geometry=["bilayer_2_3","circular_inclusion_2_3"]; ax=nexttile(tl); hold(ax,'on');
    for si=1:numel(strategies)
        pred=T.(pred_var(strategies(si))); err=abs(100*(pred-T.true_SWS)./T.true_SWS); y=nan(size(bins));
        for bi=1:numel(bins), idx=T.geometry==geometry(gi)&T.(binvar)==bins(bi); y(bi)=mean(err(idx),'omitnan'); end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.6,'DisplayName',pretty_strategy(strategies(si)));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',labels,'TickLabelInterpreter','none'); xlabel(ax,xlabel_text); ylabel(ax,'MAPE (%)');
    title(ax,ternary(gi==1,'Bilayer 2/3','Circular inclusion 2/3')); style_axes(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off'); export_fig(fig,OUT,file);
end

function plot_error_scatter(T,CORR,OUT)
baseline=T.sws_hybrid_baseline; corrected=T.(pred_var(CORR.SelectedStrategy));
e0=abs(100*(baseline-T.true_SWS)./T.true_SWS); e1=abs(100*(corrected-T.true_SWS)./T.true_SWS);
if height(T)>25000, rng(CORR.RandomSeed); idx=randperm(height(T),25000); else, idx=1:height(T); end
fig=figure('Color','w','Position',[100 100 1100 480]); tl=tiledlayout(1,2,'TileSpacing','compact');
for pi=1:2, ax=nexttile(tl); color=T.patch_purity; label='Patch purity'; if pi==2, color=T.distance_to_interface_mm; label='Distance (mm)'; end
    scatter(ax,e0(idx),e1(idx),9,color(idx),'filled','MarkerFaceAlpha',.22); hold(ax,'on'); lim=prctile([e0(idx);e1(idx)],99.5); plot(ax,[0 lim],[0 lim],'k--'); xlim(ax,[0 lim]); ylim(ax,[0 lim]); axis(ax,'square'); xlabel(ax,'Hybrid error (%)'); ylabel(ax,'Corrected error (%)'); title(ax,label); colorbar(ax); style_axes(ax); end
export_fig(fig,OUT,'test27_baseline_vs_corrected_scatter.png');
end

function plot_local_vs_hybrid(T,~,OUT)
eL=abs(100*(T.sws_local_baseline-T.true_SWS)./T.true_SWS); eH=abs(100*(T.sws_hybrid_baseline-T.true_SWS)./T.true_SWS);
classes=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"];
L=nan(size(classes)); H=L; W=L;
for i=1:numel(classes), idx=T.purity_bin==classes(i); L(i)=mean(eL(idx)); H(i)=mean(eH(idx)); W(i)=100*mean(eL(idx)<eH(idx)); end
fig=figure('Color','w','Position',[100 100 1050 450]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); bar(ax,[H' L']); ylabel(ax,'MAPE (%)'); legend(ax,["Hybrid","Local"],'Box','off'); style_axes(ax);
ax=nexttile(tl); bar(ax,W); ylabel(ax,'Pixels where Local wins (%)'); style_axes(ax);
for ax=findall(fig,'Type','axes')', set(ax,'XTick',1:numel(classes),'XTickLabel',replace(classes,"_"," "),'XTickLabelRotation',18); end
export_fig(fig,OUT,'test27_local_vs_hybrid.png');
end

function plot_small_window_diagnostics(T,~,OUT)
idx=T.small_window_available&isfinite(T.small_window_sws); base=T.sws_hybrid_baseline;
fig=figure('Color','w','Position',[100 100 1100 450]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); scatter(ax,base(idx),T.small_window_sws(idx),8,T.confidence(idx),'filled','MarkerFaceAlpha',.2); hold(ax,'on'); plot(ax,[.5 4],[.5 4],'k--'); xlabel(ax,'Hybrid SWS'); ylabel(ax,'Small-window SWS'); colorbar(ax); style_axes(ax);
ax=nexttile(tl); scatter(ax,T.base_spectral_width(idx),T.small_window_spectral_width(idx),8,T.confidence(idx),'filled','MarkerFaceAlpha',.2); hold(ax,'on'); lim=max([xlim(ax) ylim(ax)]); plot(ax,[0 lim],[0 lim],'k--'); xlabel(ax,'Baseline spectral width'); ylabel(ax,'Small-window spectral width'); colorbar(ax); style_axes(ax);
export_fig(fig,OUT,'test27_small_window_diagnostics.png');
end

function plot_edge_diagnostics(T,CORR,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
C=T(T.condition_key==keys(1),:); corrected=C.(pred_var(CORR.SelectedStrategy));
maps={map_from_rows(C,C.pseudo_region),map_from_rows(C,C.confidence<.8), ...
    map_from_rows(C,C.correction_source_c080),map_from_rows(C,corrected)};
fig=figure('Color','w','Position',[100 100 1100 300]); tl=tiledlayout(1,4,'TileSpacing','compact'); titles=["Pseudo-region","Low confidence","Correction source","Corrected SWS"];
for i=1:4, ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i)); colorbar(ax); end
export_fig(fig,OUT,'test27_edge_aware_diagnostics.png');
end

function plot_representative_maps(T,CORR,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
plot_map_panel(T(T.condition_key==keys(1),:),CORR,fullfile(OUT.figure_dir,'test27_representative_maps.png'));
end

function plot_all_maps(T,CORR,OUT)
root=fullfile(OUT.figure_dir,'maps_by_condition_dx200um'); keys=unique(T.condition_key,'stable');
for i=1:numel(keys)
    C=T(T.condition_key==keys(i),:); folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end
    file=fullfile(folder,"test27_maps__"+sanitize(keys(i))+".png");
    if exist(file,'file')~=2, plot_map_panel(C,CORR,file); end
    if mod(i,25)==0||i==numel(keys), fprintf('  Test 27 maps: %d/%d.\n',i,numel(keys)); end
end
end

function plot_map_panel(C,CORR,file)
C=sortrows(C,{'map_iz','map_ix'}); H=C.sws_hybrid_baseline; L=C.sws_local_baseline; Q=C.(pred_var(CORR.SelectedStrategy));
eH=abs(100*(H-C.true_SWS)./C.true_SWS); eL=abs(100*(L-C.true_SWS)./C.true_SWS); eQ=abs(100*(Q-C.true_SWS)./C.true_SWS);
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,H),map_from_rows(C,L),map_from_rows(C,Q), ...
    map_from_rows(C,eH),map_from_rows(C,eL),map_from_rows(C,eQ),map_from_rows(C,C.confidence), ...
    map_from_rows(C,C.confidence<.8),map_from_rows(C,C.pseudo_region),map_from_rows(C,Q-H)};
titles=["True SWS","Hybrid","Local","Selected correction","Hybrid error","Local error","Corrected error","Confidence","Confidence <0.8","Pseudo-region","Corrected - Hybrid"];
fig=figure('Color','w','Visible','off','Position',[30 60 1500 720]); tl=tiledlayout(2,6,'TileSpacing','compact','Padding','compact');
for i=1:11, ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',9); colorbar(ax); if i<=4, clim(ax,[1.5 3.5]); elseif i>=5&&i<=7, clim(ax,[0 50]); elseif i==8||i==9, clim(ax,[0 1]); end, end
title(tl,C.geometry(1)+" | "+C.field_regime(1)+" | f="+string(C.f0(1))+" | M="+string(C.M(1)),'Interpreter','none'); exportgraphics(fig,file,'Resolution',170); close(fig);
end

function strategies=plot_strategies(CORR)
strategies=["hybrid_baseline","local_baseline", ...
    "confidence_blend_local_hybrid_050_080", ...
    "small_window_if_operationally_better_c080", ...
    "edge_aware_interpolation_low_confidence_c080",CORR.SelectedStrategy];
end

function name=pred_var(strategy), name=matlab.lang.makeValidName("sws_"+strategy); end

function labels=pretty_strategy(x)
labels=replace(x,["hybrid_baseline","local_baseline","confidence_blend_local_hybrid_050_080", ...
    "small_window_if_operationally_better_c080","edge_aware_interpolation_low_confidence_c080", ...
    "small_window_then_edge_aware_c080"],["Hybrid","Local","Confidence blend","Small window","Edge-aware","Small + edge"]);
end

function print_interpretation(S,L,W,E,CORR)
S=S(~S.diagnostic_only,:); [~,ig]=min(S.MAPE); mixed=mean([S.MAPE_moderately_mixed S.MAPE_strongly_mixed],2,'omitnan'); [~,im]=min(mixed); pure=mean([S.MAPE_pure S.MAPE_near_pure],2,'omitnan'); [~,ip]=min(pure);
fprintf('\nInterpretation:\n  Best global: %s (%.2f%%).\n',S.strategy_name(ig),S.MAPE(ig));
fprintf('  Best mixed: %s (%.2f%%).\n',S.strategy_name(im),mixed(im));
fprintf('  Best pure/near-pure: %s (%.2f%%).\n',S.strategy_name(ip),pure(ip));
fprintf('  Local wins %.1f%% of evaluated pixels; mean disagreement %.3f m/s.\n',100*sum(L.local_wins_pct.*L.N)/sum(L.N)/100,sum(L.mean_local_hybrid_disagreement.*L.N)/sum(L.N));
fprintf('  Small-window accepted in %.1f%% of low-confidence opportunities.\n',100*sum(W.accepted_count)/max(sum(W.low_count),1));
fprintf('  Edge-aware corrected %.1f%% of low-confidence opportunities.\n',100*sum(E.edge_corrected_count)/max(sum(E.low_count),1));
base=S(S.strategy_name=="hybrid_baseline",:); selected=S(S.strategy_name==CORR.SelectedStrategy,:);
fprintf('  Selected homogeneous delta: %+.2f points.\n',selected.MAPE_homogeneous-base.MAPE_homogeneous);
fprintf('  Candidate output: %s + frozen confidence/reliability map.\n',CORR.SelectedStrategy);
end

%% Validation and generic helpers

function validate_test27(FILES,spectral_dir,CORR)
S=load(FILES(1).path,'RESULT'); R=S.RESULT;
required={'prediction_cube','strategy_specs','models','detectors','T_compact'};
assert(all(isfield(R,required)),'Test 26 checkpoint lacks required fields.');
base=load_test26_condition(FILES(1).path,CORR);
assert(height(base.meta)==numel(base.hybrid)&&numel(base.local)==numel(base.hybrid));
assert(all(isfinite(base.hybrid)&isfinite(base.local)));
labels=pseudo_segment(base.local,base.hybrid,base.meta,CORR);
[edge,~]=edge_aware_correct(base.hybrid,base.local,base.meta.confidence,labels,base.meta,.8,CORR);
assert(isequal(size(edge),size(base.hybrid))&&all(isfinite(edge)));
fallback=base.hybrid; bad=~isfinite(fallback)|fallback<CORR.PhysicalSwsRange(1)|fallback>CORR.PhysicalSwsRange(2);
fallback(bad)=base.hybrid(bad); assert(all(isfinite(fallback)));
sf=load_spectral_features(spectral_dir,FILES(1).key,base.meta);
assert(numel(sf.width)==height(base.meta));
small=struct('available',false,'meta',table(),'hybrid',[],'local',[], ...
    'width',[],'tail',[],'target_M',NaN);
context=build_operational_context(base,small,sf,CORR);
[cube,specs,diagnostics]=evaluate_strategies(base,context,CORR);
assert_operational_outputs(cube,specs,base.meta,CORR);
assert(height(specs)==19&&sum(specs.diagnostic_only)==1);
assert(all(isfinite(diagnostics.reliability))&& ...
    isequal(size(diagnostics.pseudo_region),size(base.hybrid)));
banned=["true_SWS","patch_purity","material_side","distance_to_interface","q_true","q_theory"];
operational_inputs=["hybrid","local","confidence","spectral_width","highk_tail","coordinates"];
assert(isempty(intersect(lower(banned),lower(operational_inputs))));
end

function tf=result_is_compatible(R,CORR)
tf=isfield(R,'config_signature')&&R.config_signature==config_signature(CORR)&& ...
    isfield(R,'T_compact')&&isfield(R,'T_stats');
end

function s=config_signature(CORR)
s=CORR.RunMode+"__dx"+string(round(1e6*CORR.Dx))+"__"+ ...
    strjoin(string(CORR.ConfidenceThresholds),"_")+"__v1";
end

function [A,iz,ix]=map_from_rows(T,v)
[~,~,iz]=unique(T.map_iz,'sorted'); [~,~,ix]=unique(T.map_ix,'sorted');
A=nan(max(iz),max(ix)); A(sub2ind(size(A),iz,ix))=v;
end

function style_axes(ax)
set(ax,'FontName','Helvetica','FontSize',9,'Box','off','TickDir','out'); grid(ax,'on'); ax.GridAlpha=.15;
end

function export_fig(fig,OUT,name)
for ax=findall(fig,'Type','axes')'
    try
        ax.Toolbar.Visible='off';
    catch
    end
end
exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',220); close(fig);
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
if isstring(example), x=repmat("",n,1); elseif islogical(example), x=false(n,1);
elseif isnumeric(example), x=nan(n,1); else, x=cell(n,1); end
end

function value=ternary(condition,a,b), if condition, value=a; else, value=b; end, end
function s=sanitize(x), s=regexprep(char(x),'[^A-Za-z0-9_-]','_'); end
