%% analyze_test_31_simple_confidence_interpolation.m
% Test 31: simple confidence-based q/SWS interpolation.
%
% This diagnostic compares simple low-confidence corrections that copy or
% interpolate either LocalOnly_T18 q or LocalOnly_T18 SWS from high-confidence
% pixels. It also compares against TheoryQDiscrete and the Test 30
% Theory-structure / Local-region-level result when available.
%
% Operational correction inputs are limited to frozen model outputs, frozen
% confidence, coordinates, and structures inferred from Local or Theory maps.
% true_SWS, material_side, patch_purity, and distance_to_interface are used
% only after prediction for evaluation.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST31_MODE          = quick | full
%   ADAPTIVE_REQ_TEST31_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST31_SAVE_ALL_MAPS = true | false

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.12, ...
    'defaultAxesLabelFontSizeMultiplier',1.04);

%% Configuration

CFG=struct();
mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST31_MODE'))));
mode = 'full';
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]),'ADAPTIVE_REQ_TEST31_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST31_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST31_SAVE_ALL_MAPS',true);
CFG.Version=2; CFG.Dx=0.2e-3; CFG.M=[2 3];
CFG.ConfidenceThreshold=env_double('ADAPTIVE_REQ_TEST31_CONFIDENCE_THRESHOLD',.8);
CFG.PhysicalRange=[.5 10]; CFG.MinimumRegionHighConf=6; CFG.NearestRadiusMaxPx=inf;
CFG.Geometries=["homogeneous_cs2","homogeneous_cs3","bilayer_2_3","circular_inclusion_2_3"];
if CFG.QuickMode
    CFG.Frequencies=500;
    CFG.Regimes=["directional_2D","diffuse_3D"];
else
    CFG.Frequencies=[300 400 500 600];
    CFG.Regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
end
CFG.Strategies=["local_baseline","theory_baseline", ...
    "q_nearest_highconf","sws_nearest_highconf", ...
    "q_mean_highconf_global","q_median_highconf_global", ...
    "q_median_highconf_by_region_local_structure", ...
    "sws_interp_by_region_local_structure", ...
    "test30_theory_region_levels", ...
    "edgeaware_sws_theory_structure","edgeaware_q_theory_structure", ...
    "edgeaware_sws_local_structure","edgeaware_q_local_structure"];
CFG.QStrategies=["q_nearest_highconf","q_mean_highconf_global", ...
    "q_median_highconf_global","q_median_highconf_by_region_local_structure", ...
    "edgeaware_q_theory_structure","edgeaware_q_local_structure"];
CFG.MainStrategies=["local_baseline","theory_baseline","q_nearest_highconf", ...
    "sws_nearest_highconf","q_median_highconf_by_region_local_structure", ...
    "sws_interp_by_region_local_structure","test30_theory_region_levels", ...
    "edgeaware_sws_theory_structure","edgeaware_q_theory_structure"];

OUT=make_output_dirs(root_dir,CFG);
write_config_json(CFG,fullfile(OUT.root_dir,'test31_configuration.json'));

SOURCE=resolve_test29_source(root_dir,CFG);
SPECTRAL_DIR=resolve_spectral_dir(root_dir,CFG);
T30_TABLE=resolve_test30_table(root_dir,CFG);

fprintf('\nTest 31: simple confidence-based q/SWS interpolation\n');
fprintf('Mode: %s | dx=%.1f mm | confidence threshold %.2f.\n',CFG.RunMode,1e3*CFG.Dx,CFG.ConfidenceThreshold);
fprintf('Source Test 29: %s\n',SOURCE);
fprintf('Spectral cache: %s\n',SPECTRAL_DIR);

T=load_source_patch_table(SOURCE,CFG);
T30=load_test30_table(T30_TABLE,CFG);
keys=unique(T.condition_key,'stable');
assert(~isempty(keys),'No Test 31 conditions matched the requested design.');
fprintf('Matched %d conditions and %d patches.\n',numel(keys),height(T));

if CFG.ValidateOnly
    validate_test31(T,T30,SPECTRAL_DIR,keys(1),CFG);
    fprintf('Test 31 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Evaluate each condition

parts=cell(numel(keys),1);
for ki=1:numel(keys)
    timer=tic; key=keys(ki);
    checkpoint=fullfile(OUT.condition_dir,"test31__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        S=load(checkpoint,'RESULT');
        if isfield(S,'RESULT')&&checkpoint_is_compatible(S.RESULT,CFG)
            parts{ki}=S.RESULT.T_patch;
            fprintf('[%d/%d] Reused %s (%d patches).\n',ki,numel(keys),key,height(parts{ki}));
            continue;
        end
    end
    X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
    SPEC=load_spectral_cache(SPECTRAL_DIR,key,X);
    R=evaluate_condition(X,SPEC,T30,CFG);
    RESULT=struct('version',CFG.Version,'cache_signature',cache_signature(CFG),'T_patch',R);
    save(checkpoint,'RESULT','-v7.3');
    parts{ki}=R;
    fprintf('[%d/%d] %s: low-conf %d/%d, q->SWS %.1f s total.\n', ...
        ki,numel(keys),key,sum(X.confidence<CFG.ConfidenceThreshold),height(X),toc(timer));
end
T_patch=vertcat(parts{:}); clear parts T T30;

%% Tables

T_overall=summarize_predictions(T_patch,"strategy_name",CFG);
T_by_M=summarize_predictions(T_patch,["strategy_name","M"],CFG);
T_by_frequency=summarize_predictions(T_patch,["strategy_name","f0"],CFG);
T_by_regime=summarize_predictions(T_patch,["strategy_name","field_regime"],CFG);
T_by_geometry=summarize_predictions(T_patch,["strategy_name","geometry"],CFG);
T_by_purity=summarize_predictions(T_patch,["strategy_name","geometry","purity_bin"],CFG);
T_by_distance=summarize_predictions(T_patch,["strategy_name","geometry","distance_bin"],CFG);
T_region_freq=region_sws_vs_frequency(T_patch,CFG);
T_roi_freq=roi_sws_vs_frequency(T_patch,CFG);
T_q=q_correction_summary(T_patch,CFG);
T_best=best_strategy_candidates(T_overall,T_by_purity,T_roi_freq,CFG);

writetable(T_patch,fullfile(OUT.table_dir,'test31_patch_level_results.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test31_strategy_summary_overall.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'test31_strategy_summary_by_M.csv'));
writetable(T_by_frequency,fullfile(OUT.table_dir,'test31_strategy_summary_by_frequency.csv'));
writetable(T_by_regime,fullfile(OUT.table_dir,'test31_strategy_summary_by_regime.csv'));
writetable(T_by_geometry,fullfile(OUT.table_dir,'test31_strategy_summary_by_geometry.csv'));
writetable(T_by_purity,fullfile(OUT.table_dir,'test31_strategy_summary_by_purity_bin.csv'));
writetable(T_by_distance,fullfile(OUT.table_dir,'test31_strategy_summary_by_distance_bin.csv'));
writetable(T_region_freq,fullfile(OUT.table_dir,'test31_region_sws_vs_frequency.csv'));
writetable(T_roi_freq,fullfile(OUT.table_dir,'test31_roi_sws_vs_frequency.csv'));
writetable(T_q,fullfile(OUT.table_dir,'test31_q_correction_summary.csv'));
writetable(T_best,fullfile(OUT.table_dir,'test31_best_strategy_candidates.csv'));
save(fullfile(OUT.data_dir,'test31_compact_results.mat'), ...
    'T_patch','T_overall','T_by_purity','T_by_distance','T_region_freq','T_roi_freq','T_q','CFG','-v7.3');

%% Figures

plot_strategy_ranking(T_patch,CFG,OUT);
plot_sws_region_vs_frequency(T_region_freq,CFG,OUT);
plot_roi_vs_frequency(T_roi_freq,CFG,OUT);
plot_error_vs_distance(T_patch,CFG,OUT);
plot_error_vs_confidence(T_patch,CFG,OUT);
plot_correction_vs_distance(T_patch,CFG,OUT);
plot_q_correction_summary(T_patch,CFG,OUT);
plot_local_corrected_scatter(T_patch,CFG,OUT);
plot_representative_maps(T_patch,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T_patch,CFG,OUT); end

print_interpretation(T_overall,T_by_distance,T_roi_freq,T_q,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 31 complete.\n',OUT.table_dir,OUT.figure_dir);

%% Input loading

function OUT=make_output_dirs(root,CFG)
OUT.root_dir=fullfile(root,'outputs','test_31_simple_confidence_interpolation');
if CFG.QuickMode, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables');
OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data');
OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
for d=string(struct2cell(OUT))', if exist(d,'dir')~=7, mkdir(d); end, end
end

function tf=checkpoint_is_compatible(RESULT,CFG)
tf=isfield(RESULT,'version')&&RESULT.version==CFG.Version&& ...
    isfield(RESULT,'cache_signature')&&isequal(RESULT.cache_signature,cache_signature(CFG));
end

function sig=cache_signature(CFG)
sig=struct();
sig.Version=CFG.Version;
sig.RunMode=CFG.RunMode;
sig.Dx=CFG.Dx;
sig.M=CFG.M;
sig.Frequencies=CFG.Frequencies;
sig.Regimes=CFG.Regimes;
sig.Geometries=CFG.Geometries;
sig.ConfidenceThreshold=CFG.ConfidenceThreshold;
sig.Strategies=CFG.Strategies;
sig.QStrategies=CFG.QStrategies;
end

function SOURCE=resolve_test29_source(root,CFG)
base=fullfile(root,'outputs','test_29_oracle_halfwindow_graph_tv');
if CFG.QuickMode
    candidates=[string(fullfile(base,'data','test29_compact_results.mat')), ...
        string(fullfile(base,'quick','data','test29_compact_results.mat'))];
else
    candidates=string(fullfile(base,'data','test29_compact_results.mat'));
end
SOURCE="";
for c=candidates
    if exist(c,'file')==2, SOURCE=c; break; end
end
assert(SOURCE~="",'Test 29 compact source is missing. Run Test 29 first.');
end

function folder=resolve_spectral_dir(root,CFG)
base=fullfile(root,'outputs','test_26_confidence_gated_corrections');
if CFG.QuickMode
    candidates=[string(fullfile(base,'data','spectral_checkpoints')), ...
        string(fullfile(base,'quick','data','spectral_checkpoints'))];
else
    candidates=string(fullfile(base,'data','spectral_checkpoints'));
end
folder="";
for c=candidates
    if exist(c,'dir')==7, folder=c; break; end
end
assert(folder~="",'Spectral cache missing. Run Test 26/Test 29 cache generation first.');
end

function file=resolve_test30_table(root,CFG)
base=fullfile(root,'outputs','test_30_theory_structure_local_m2');
if CFG.QuickMode
    candidates=[string(fullfile(base,'tables','test30_patch_level_results.csv')), ...
        string(fullfile(base,'quick','tables','test30_patch_level_results.csv'))];
else
    candidates=string(fullfile(base,'tables','test30_patch_level_results.csv'));
end
file="";
for c=candidates
    if exist(c,'file')==2, file=c; break; end
end
if file=="", warning('Test31:MissingTest30','Test 30 table not found; Test 31 will reconstruct theory structures but Test30 region levels become NaN.'); end
end

function T=load_source_patch_table(source,CFG)
S=load(source,'T_patch'); T=S.T_patch;
T=T(abs(T.dx-CFG.Dx)<1e-12&ismember(T.M,CFG.M)& ...
    ismember(T.geometry,CFG.Geometries)&ismember(T.f0,CFG.Frequencies)& ...
    ismember(T.field_regime,CFG.Regimes),:);
if ~ismember('sws_theory_discrete',T.Properties.VariableNames)
    error('Test31:MissingTheory','Test 29 table must include sws_theory_discrete. Re-run the updated Test 29 if needed.');
end
T=sortrows(T,{'condition_key','map_iz','map_ix'});
end

function T30=load_test30_table(file,CFG)
if file=="", T30=table(); return; end
T30=readtable(file,'TextType','string');
T30=T30(abs(T30.dx-CFG.Dx)<1e-12&ismember(T30.geometry,CFG.Geometries)& ...
    ismember(T30.f0,CFG.Frequencies)&ismember(T30.field_regime,CFG.Regimes),:);
end

function SPEC=load_spectral_cache(folder,key,X)
file=fullfile(folder,"spectra__"+sanitize(key)+".mat");
assert(exist(file,'file')==2,'Missing spectral cache: %s',file);
S=load(file,'SPEC'); SPEC=S.SPEC; assert(numel(SPEC.curves)==height(X));
cx=round(X.x/X.dx(1))+1; cz=round(X.z/X.dz(1))+1;
assert(isequal(double(SPEC.cx(:)),double(cx(:)))&&isequal(double(SPEC.cz(:)),double(cz(:))), ...
    'Spectral cache alignment failed for %s.',key);
end

%% Per-condition evaluation

function R=evaluate_condition(X,SPEC,T30,CFG)
X=sortrows(X,{'map_iz','map_ix'});
n=height(X); low=X.confidence<CFG.ConfidenceThreshold; high=~low;
local=double(X.sws_local_baseline); theory=double(X.sws_theory_discrete);
q_local=sws_to_q(SPEC,local,X.f0(1));

local_region=segment_map(X,local,"local",CFG);
theory_region=segment_map(X,theory,"theory",CFG);
[test30_sws,test30_region]=align_test30(X,T30,CFG);
if all(~isfinite(test30_sws))
    test30_sws=region_level_sws(local,X.confidence,theory_region,CFG);
    test30_region=theory_region;
end

[q_near,dist_q_near,src_q_near]=nearest_fill_map(X,q_local,high,ones(n,1),q_local);
[sws_near,dist_sws_near,src_sws_near]=nearest_fill_map(X,local,high,ones(n,1),local);

q_mean=q_local; if any(high), q_mean(low)=mean(q_local(high),'omitnan'); end
q_median=q_local; if any(high), q_median(low)=median(q_local(high),'omitnan'); end
q_region=region_stat_fill(q_local,high,low,local_region,"median");

[sws_region_local,dist_sws_region_local,src_sws_region_local]= ...
    nearest_fill_map(X,local,high,local_region,local);
[sws_edge_theory,dist_sws_edge_theory,src_sws_edge_theory]= ...
    nearest_fill_map(X,local,high,theory_region,local);
[q_edge_theory,dist_q_edge_theory,src_q_edge_theory]= ...
    nearest_fill_map(X,q_local,high,theory_region,q_local);
[sws_edge_local,dist_sws_edge_local,src_sws_edge_local]= ...
    nearest_fill_map(X,local,high,local_region,local);
[q_edge_local,dist_q_edge_local,src_q_edge_local]= ...
    nearest_fill_map(X,q_local,high,local_region,q_local);

pred=struct();
pred.local_baseline=local;
pred.theory_baseline=theory;
pred.q_nearest_highconf=q_to_sws(SPEC,q_near,X.f0(1),local,CFG);
pred.sws_nearest_highconf=sws_near;
pred.q_mean_highconf_global=q_to_sws(SPEC,q_mean,X.f0(1),local,CFG);
pred.q_median_highconf_global=q_to_sws(SPEC,q_median,X.f0(1),local,CFG);
pred.q_median_highconf_by_region_local_structure=q_to_sws(SPEC,q_region,X.f0(1),local,CFG);
pred.sws_interp_by_region_local_structure=sws_region_local;
pred.test30_theory_region_levels=test30_sws;
pred.edgeaware_sws_theory_structure=sws_edge_theory;
pred.edgeaware_q_theory_structure=q_to_sws(SPEC,q_edge_theory,X.f0(1),local,CFG);
pred.edgeaware_sws_local_structure=sws_edge_local;
pred.edgeaware_q_local_structure=q_to_sws(SPEC,q_edge_local,X.f0(1),local,CFG);

qcorr=struct('q_nearest_highconf',q_near,'q_mean_highconf_global',q_mean, ...
    'q_median_highconf_global',q_median, ...
    'q_median_highconf_by_region_local_structure',q_region, ...
    'edgeaware_q_theory_structure',q_edge_theory, ...
    'edgeaware_q_local_structure',q_edge_local);

assert_policy(pred,qcorr,X,CFG);
R=attach_outputs(X,pred,qcorr,q_local,CFG);
R.low_confidence=low; R.high_confidence=high;
R.local_structure=local_region; R.theory_structure=theory_region; R.test30_structure=test30_region;
R.donor_distance_q_nearest_highconf_mm=dist_q_near; R.correction_source_q_nearest_highconf=src_q_near;
R.donor_distance_sws_nearest_highconf_mm=dist_sws_near; R.correction_source_sws_nearest_highconf=src_sws_near;
R.donor_distance_sws_interp_by_region_local_structure_mm=dist_sws_region_local; R.correction_source_sws_interp_by_region_local_structure=src_sws_region_local;
R.donor_distance_edgeaware_sws_theory_structure_mm=dist_sws_edge_theory; R.correction_source_edgeaware_sws_theory_structure=src_sws_edge_theory;
R.donor_distance_edgeaware_q_theory_structure_mm=dist_q_edge_theory; R.correction_source_edgeaware_q_theory_structure=src_q_edge_theory;
R.donor_distance_edgeaware_sws_local_structure_mm=dist_sws_edge_local; R.correction_source_edgeaware_sws_local_structure=src_sws_edge_local;
R.donor_distance_edgeaware_q_local_structure_mm=dist_q_edge_local; R.correction_source_edgeaware_q_local_structure=src_q_edge_local;
end

function labels=segment_map(meta,value,kind,CFG)
[A,iz,ix]=map_from_rows(meta,value); A=medfilt2(A,[3 3],'symmetric');
if meta.geometry_type(1)=="homogeneous"
    L=ones(size(A)); labels=L(sub2ind(size(L),iz,ix)); return;
end
x=A(isfinite(A)); if isempty(x), L=ones(size(A)); labels=L(sub2ind(size(L),iz,ix)); return; end
c=[prctile(x,25) prctile(x,75)];
for iter=1:20
    first=abs(x-c(1))<=abs(x-c(2));
    if any(first), c(1)=median(x(first)); end
    if any(~first), c(2)=median(x(~first)); end
end
c=sort(c); threshold=mean(c);
switch meta.geometry_type(1)
    case "bilayer"
        col=median(A,1,'omitnan'); row=median(A,2,'omitnan');
        if range(col)>=range(row)
            split=col>threshold; L=repmat(split,size(A,1),1)+1;
        else
            split=row>threshold; L=repmat(split,1,size(A,2))+1;
        end
    case "inclusion"
        raw=A>threshold;
        if mean(raw(:),'omitnan')>.5, raw=~raw; end
        if kind=="local", raw=medfilt2(double(raw),[3 3],'symmetric')>.5;
        else, raw=medfilt2(double(raw),[5 5],'symmetric')>.5; end
        components=bwconncomp(raw,8);
        if components.NumObjects>0
            sizes=cellfun(@numel,components.PixelIdxList); [~,largest]=max(sizes);
            clean=false(size(raw)); clean(components.PixelIdxList{largest})=true;
            raw=imfill(clean,'holes');
        end
        L=ones(size(A)); L(raw)=2;
    otherwise
        L=ones(size(A)); L(A>threshold)=2;
end
labels=L(sub2ind(size(L),iz,ix));
% Keep region 2 as the faster/harder side when possible.
if numel(unique(labels))==2
    m1=median(value(labels==1),'omitnan'); m2=median(value(labels==2),'omitnan');
    if m1>m2, labels=3-labels; end
end
end

function [sws,region]=align_test30(X,T30,~)
sws=nan(height(X),1); region=nan(height(X),1);
if isempty(T30), return; end
base_key=regexprep(string(X.condition_key(1)),'__M[0-9]+__step[0-9]+px$','');
T30.base_key=regexprep(string(T30.condition_key),'__M[0-9]+__step[0-9]+px$','');
C=T30(T30.base_key==base_key,:);
if isempty(C), return; end
ok=isfinite(C.sws_userguess_region_levels);
F=scatteredInterpolant(C.x(ok),C.z(ok),C.sws_userguess_region_levels(ok),'nearest','nearest');
sws=F(X.x,X.z);
if ismember('userguess_region',C.Properties.VariableNames)
    G=scatteredInterpolant(C.x,C.z,double(C.userguess_region),'nearest','nearest');
    region=round(G(X.x,X.z));
end
end

function sws=region_level_sws(local,confidence,region,CFG)
sws=local;
for r=unique(region(:))'
    member=region==r; seed=member&confidence>=CFG.ConfidenceThreshold&isfinite(local);
    if nnz(seed)<CFG.MinimumRegionHighConf, seed=member&isfinite(local); end
    if any(seed), sws(member)=median(local(seed),'omitnan'); end
end
end

function [filled,dist_mm,source]=nearest_fill_map(meta,value,high,region,fallback)
[V,iz,ix]=map_from_rows(meta,value); [H,~,~]=map_from_rows(meta,high);
[R,~,~]=map_from_rows(meta,region); [F,~,~]=map_from_rows(meta,fallback);
out=V; dist_map=zeros(size(V)); source_map=strings(size(V)); source_map(:)="kept_high_conf";
low=~H;
dx_mm=1e3*mean([median(diff(unique(meta.x))) median(diff(unique(meta.z)))],'omitnan');
for rr=unique(R(isfinite(R)))'
    member=R==rr; donor=member&H&isfinite(V);
    target=member&low;
    if ~any(target(:)), continue; end
    if any(donor(:))
        [D,I]=bwdist(donor);
        out(target)=V(I(target)); dist_map(target)=D(target)*dx_mm;
        source_map(target)="nearest_same_region_high_conf";
    else
        med=median(V(member&isfinite(V)),'omitnan');
        if ~isfinite(med), med=median(F(member&isfinite(F)),'omitnan'); end
        if isfinite(med)
            out(target)=med; dist_map(target)=NaN; source_map(target)="regional_median_fallback";
        else
            out(target)=F(target); dist_map(target)=NaN; source_map(target)="baseline_fallback";
        end
    end
end
filled=out(sub2ind(size(out),iz,ix));
dist_mm=dist_map(sub2ind(size(out),iz,ix));
source=source_map(sub2ind(size(out),iz,ix));
end

function q=region_stat_fill(q_local,high,low,region,stat)
q=q_local;
for r=unique(region(:))'
    member=region==r; seed=member&high&isfinite(q_local);
    if ~any(seed), seed=member&isfinite(q_local); end
    if ~any(seed), continue; end
    switch stat
        case "mean", value=mean(q_local(seed),'omitnan');
        otherwise, value=median(q_local(seed),'omitnan');
    end
    q(member&low)=value;
end
q=clip_q(q);
end

function sws=q_to_sws(SPEC,q,f0,fallback,CFG)
q=clip_q(q); sws=fallback;
for i=1:numel(q)
    if isfinite(q(i))&&~isempty(SPEC.curves{i})
        value=adaptive_req.quantile.quantile_to_cs(SPEC.curves{i},q(i),f0);
        if isfinite(value)&&value>=CFG.PhysicalRange(1)&&value<=CFG.PhysicalRange(2)
            sws(i)=value;
        end
    end
end
end

function q=sws_to_q(SPEC,sws,f0)
q=nan(numel(sws),1);
for i=1:numel(sws)
    if ~isfinite(sws(i))||sws(i)<=0||isempty(SPEC.curves{i}), continue; end
    curve=SPEC.curves{i};
    k=double(curve.k_cent(:)); E=double(curve.Ecum(:));
    valid=isfinite(k)&isfinite(E);
    k=k(valid); E=E(valid);
    if numel(k)<2, continue; end
    [kuniq,ia]=unique(k,'stable'); Euniq=E(ia);
    target=2*pi*f0/sws(i);
    q(i)=interp1(kuniq,Euniq,target,'linear','extrap');
end
q=clip_q(q);
bad=~isfinite(q); q(bad)=0.5;
end

function q=clip_q(q)
q=min(max(q,0.001),0.999);
end

function assert_policy(pred,qcorr,X,CFG)
for s=CFG.Strategies
    y=pred.(s); assert(numel(y)==height(X)&&all(isfinite(y))&& ...
        all(y>=CFG.PhysicalRange(1)&y<=CFG.PhysicalRange(2)), ...
        'Invalid SWS output for %s.',s);
end
for s=CFG.QStrategies
    q=qcorr.(s); assert(numel(q)==height(X)&&all(isfinite(q))&&all(q>0&q<1), ...
        'Invalid q output for %s.',s);
end
operational_inputs=["LocalOnly_T18","TheoryQDiscrete","q_local","confidence","coordinates","estimated_structure"];
oracle=["true_SWS","patch_purity","material_side","distance_to_interface","q_true"];
assert(isempty(intersect(lower(operational_inputs),lower(oracle))));
end

function R=attach_outputs(X,pred,qcorr,q_local,CFG)
keep={'condition_key','geometry','geometry_type','field_regime','f0','M','dx','dz', ...
    'map_iz','map_ix','x','z','true_SWS','confidence','patch_purity','material_side', ...
    'purity_bin','distance_bin','distance_to_interface_mm'};
R=X(:,keep);
R.q_local=q_local;
for s=CFG.Strategies
    R.("sws_"+s)=pred.(s);
    R.("signed_error_"+s)=100*(pred.(s)-R.true_SWS)./R.true_SWS;
    R.("abs_error_"+s)=abs(R.("signed_error_"+s));
    R.("delta_vs_local_"+s)=pred.(s)-pred.local_baseline;
end
for s=CFG.QStrategies
    R.("q_corr_"+s)=qcorr.(s);
    R.("q_delta_"+s)=qcorr.(s)-R.q_local;
end
end

%% Summaries

function S=summarize_predictions(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for strategy=CFG.Strategies
    error=T.("signed_error_"+strategy);
    if isempty(basegroups), G=ones(height(T),1); else, G=findgroups(T(:,cellstr(basegroups))); end
    ids=unique(G,'stable')';
    for id=ids
        idx=G==id; row=table();
        for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy; row.N=sum(idx);
        row.MAPE=mean(abs(error(idx)),'omitnan');
        row.mean_signed_error=mean(error(idx),'omitnan');
        row.median_signed_error=median(error(idx),'omitnan');
        row.high_error_10_pct=100*mean(abs(error(idx))>10,'omitnan');
        row.high_error_20_pct=100*mean(abs(error(idx))>20,'omitnan');
        row.underestimate_pct=100*mean(error(idx)<0,'omitnan');
        row.overestimate_pct=100*mean(error(idx)>0,'omitnan');
        row.corrected_pct=100*mean(abs(T.("delta_vs_local_"+strategy)(idx))>1e-9,'omitnan');
        row.mean_delta_vs_local=mean(T.("delta_vs_local_"+strategy)(idx),'omitnan');
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            j=idx&T.purity_bin==cls; token=char(cls);
            row.("N_"+token)=sum(j);
            row.("MAPE_"+token)=mean(abs(error(j)),'omitnan');
            row.("high20_"+token)=100*mean(abs(error(j))>20,'omitnan');
        end
        for side=["soft","hard"]
            j=idx&T.material_side==side; token=char(side);
            row.("MAPE_"+token)=mean(abs(error(j)),'omitnan');
            row.("mean_SWS_"+token)=mean(T.("sws_"+strategy)(j),'omitnan');
            row.("std_SWS_"+token)=std(T.("sws_"+strategy)(j),'omitnan');
        end
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function T=region_sws_vs_frequency(P,CFG)
T=table(); strategies=CFG.MainStrategies;
for strategy=strategies
    for M=unique(P.M)'
        for f=unique(P.f0)'
            for region=["soft","hard"]
                idx=P.M==M&P.f0==f&P.material_side==region;
                if ~any(idx), continue; end
                row=table(strategy,M,f,region,sum(idx), ...
                    mean(P.("sws_"+strategy)(idx),'omitnan'), ...
                    std(P.("sws_"+strategy)(idx),'omitnan'), ...
                    iqr(P.("sws_"+strategy)(idx)), ...
                    'VariableNames',{'strategy_name','M','f0','material_side','N','mean_SWS','std_SWS','iqr_SWS'});
                T=concat_tables(T,row);
            end
        end
    end
end
end

function T=roi_sws_vs_frequency(P,CFG)
T=table(); strategies=CFG.MainStrategies;
for strategy=strategies
    for M=unique(P.M)'
        for f=unique(P.f0)'
            for geo=unique(P.geometry)'
                if startsWith(geo,"homogeneous"), continue; end
                for roi=["soft_core","hard_core","interface_0_2mm"]
                    switch roi
                        case "soft_core"
                            idx=P.geometry==geo&P.M==M&P.f0==f&P.material_side=="soft"&P.distance_to_interface_mm>8;
                        case "hard_core"
                            idx=P.geometry==geo&P.M==M&P.f0==f&P.material_side=="hard"&P.distance_to_interface_mm>4;
                        otherwise
                            idx=P.geometry==geo&P.M==M&P.f0==f&P.distance_to_interface_mm<=2;
                    end
                    if ~any(idx), continue; end
                    err=100*(P.("sws_"+strategy)(idx)-P.true_SWS(idx))./P.true_SWS(idx);
                    row=table(strategy,geo,M,f,roi,sum(idx), ...
                        mean(P.("sws_"+strategy)(idx),'omitnan'),std(P.("sws_"+strategy)(idx),'omitnan'), ...
                        mean(abs(err),'omitnan'),100*mean(abs(err)>20,'omitnan'), ...
                        'VariableNames',{'strategy_name','geometry','M','f0','roi','N','mean_SWS','std_SWS','MAPE','high20_pct'});
                    T=concat_tables(T,row);
                end
            end
        end
    end
end
end

function T=q_correction_summary(P,CFG)
T=table();
for strategy=CFG.QStrategies
    q=P.("q_corr_"+strategy); dq=P.("q_delta_"+strategy);
    for group=["all","low_confidence","high_confidence"]
        switch group
            case "low_confidence", idx=P.low_confidence;
            case "high_confidence", idx=~P.low_confidence;
            otherwise, idx=true(height(P),1);
        end
        row=table(strategy,group,sum(idx),mean(q(idx),'omitnan'),std(q(idx),'omitnan'), ...
            mean(abs(dq(idx)),'omitnan'),prctile(abs(dq(idx)),95), ...
            'VariableNames',{'strategy_name','confidence_group','N','mean_q_corr','std_q_corr','mean_abs_q_delta','p95_abs_q_delta'});
        T=concat_tables(T,row);
    end
end
end

function B=best_strategy_candidates(O,Purity,ROI,CFG)
B=table();
O=O(ismember(O.strategy_name,CFG.Strategies),:);
[~,i]=min(O.MAPE); B=concat_tables(B,table("global_MAE",O.strategy_name(i),O.MAPE(i), ...
    'VariableNames',{'criterion','strategy_name','value'}));
[~,i]=min(O.high_error_20_pct); B=concat_tables(B,table("global_high20",O.strategy_name(i),O.high_error_20_pct(i), ...
    'VariableNames',{'criterion','strategy_name','value'}));
P=Purity(Purity.purity_bin=="strongly_mixed",:); if ~isempty(P)
    [~,i]=min(P.MAPE); B=concat_tables(B,table("strongly_mixed_MAPE",P.strategy_name(i),P.MAPE(i), ...
        'VariableNames',{'criterion','strategy_name','value'}));
end
R=ROI(ROI.roi=="hard_core",:); if ~isempty(R)
    G=findgroups(R.strategy_name); S=splitapply(@mean,R.MAPE,G); names=splitapply(@(x)x(1),R.strategy_name,G);
    [~,i]=min(S); B=concat_tables(B,table("hard_core_MAPE",names(i),S(i), ...
        'VariableNames',{'criterion','strategy_name','value'}));
end
end

%% Figures

function plot_strategy_ranking(T,CFG,OUT)
strategies=CFG.MainStrategies; labels=["Homogeneous","Pure/near-pure","Mixed"];
M=nan(numel(strategies),3); H=M;
for si=1:numel(strategies)
    e=T.("abs_error_"+strategies(si));
    masks={T.purity_bin=="homogeneous",ismember(T.purity_bin,["pure","near_pure"]), ...
        ismember(T.purity_bin,["moderately_mixed","strongly_mixed"])};
    for j=1:3, M(si,j)=mean(e(masks{j}),'omitnan'); H(si,j)=100*mean(e(masks{j})>20,'omitnan'); end
end
fig=figure('Color','w','Position',[40 40 1400 780]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,M); xlabel(ax,'MAPE (%)'); title(ax,'Strategy ranking'); style_axes(ax);
ax=nexttile(tl); barh(ax,H); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style_axes(ax);
for ax=findall(fig,'Type','axes')', set(ax,'YTick',1:numel(strategies),'YTickLabel',pretty(strategies),'TickLabelInterpreter','none'); end
legend(labels,'Location','northoutside','Orientation','horizontal','Box','off');
export_fig(fig,OUT,'test31_strategy_ranking_heatmap.png');
end

function plot_sws_region_vs_frequency(T,CFG,OUT)
if isempty(T), return; end
for M=unique(T.M)'
    fig=figure('Color','w','Position',[80 80 1250 540]); tl=tiledlayout(1,2,'TileSpacing','compact');
    for side=["soft","hard"]
        ax=nexttile(tl); hold(ax,'on');
        for strategy=CFG.MainStrategies
            X=T(T.M==M&T.material_side==side&T.strategy_name==strategy,:);
            if isempty(X), continue; end
            X=sortrows(X,'f0');
            errorbar(ax,X.f0,X.mean_SWS,X.std_SWS,'-o','LineWidth',1.2,'DisplayName',pretty(strategy));
        end
        yline(ax,ternary(side=="soft",2,3),'k--','True');
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Mean SWS (m/s)'); title(ax,sprintf('M=%d, %s side',M,side)); style_axes(ax);
    end
    legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off');
    export_fig(fig,OUT,sprintf('test31_region_sws_vs_frequency_M%d.png',M));
end
end

function plot_roi_vs_frequency(T,CFG,OUT)
if isempty(T), return; end
for M=unique(T.M)'
    fig=figure('Color','w','Position',[70 70 1350 760]); tl=tiledlayout(2,2,'TileSpacing','compact');
    combos=unique(T(:,{'geometry','roi'}),'rows','stable');
    combos=combos(ismember(combos.roi,["soft_core","hard_core","interface_0_2mm"]),:);
    for ci=1:min(4,height(combos))
        ax=nexttile(tl); hold(ax,'on');
        for strategy=CFG.MainStrategies
            X=T(T.M==M&T.geometry==combos.geometry(ci)&T.roi==combos.roi(ci)&T.strategy_name==strategy,:);
            if isempty(X), continue; end
            X=sortrows(X,'f0'); errorbar(ax,X.f0,X.mean_SWS,X.std_SWS,'-o','LineWidth',1.1,'DisplayName',pretty(strategy));
        end
        title(ax,replace(combos.geometry(ci)+" / "+combos.roi(ci),"_"," "));
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'ROI SWS (m/s)'); style_axes(ax);
    end
    legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off');
    export_fig(fig,OUT,sprintf('test31_roi_sws_vs_frequency_M%d.png',M));
end
end

function plot_error_vs_distance(T,CFG,OUT)
bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"]; labels=["0-1","1-2","2-4","4-8",">8"];
fig=figure('Color','w','Position',[80 80 1300 520]); tl=tiledlayout(1,2,'TileSpacing','compact');
for geo=["bilayer_2_3","circular_inclusion_2_3"]
    ax=nexttile(tl); hold(ax,'on');
    for strategy=CFG.MainStrategies
        y=nan(size(bins));
        for bi=1:numel(bins)
            idx=T.geometry==geo&T.distance_bin==bins(bi);
            y(bi)=mean(T.("abs_error_"+strategy)(idx),'omitnan');
        end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.2,'DisplayName',pretty(strategy));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',labels); xlabel(ax,'Distance to interface (mm)');
    ylabel(ax,'MAPE (%)'); title(ax,replace(geo,"_"," ")); style_axes(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off');
export_fig(fig,OUT,'test31_error_vs_distance.png');
end

function plot_error_vs_confidence(T,CFG,OUT)
edges=0:.1:1; centers=(edges(1:end-1)+edges(2:end))/2;
fig=figure('Color','w','Position',[80 80 1250 500]); tl=tiledlayout(1,2,'TileSpacing','compact');
for metric=["MAPE","high20"]
    ax=nexttile(tl); hold(ax,'on');
    for strategy=CFG.MainStrategies
        y=nan(size(centers));
        for bi=1:numel(centers)
            idx=T.confidence>=edges(bi)&T.confidence<edges(bi+1);
            e=T.("abs_error_"+strategy);
            if metric=="MAPE", y(bi)=mean(e(idx),'omitnan'); else, y(bi)=100*mean(e(idx)>20,'omitnan'); end
        end
        plot(ax,centers,y,'-o','LineWidth',1.2,'DisplayName',pretty(strategy));
    end
    xline(ax,CFG.ConfidenceThreshold,'k--','threshold');
    xlabel(ax,'Confidence bin'); ylabel(ax,ternary(metric=="MAPE",'MAPE (%)','Error >20% (%)'));
    title(ax,"Error vs confidence: "+metric); style_axes(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off');
export_fig(fig,OUT,'test31_error_vs_confidence.png');
end

function plot_correction_vs_distance(T,CFG,OUT)
bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"]; labels=["0-1","1-2","2-4","4-8",">8"];
fig=figure('Color','w','Position',[100 100 1050 450]); ax=axes(fig); hold(ax,'on');
for strategy=CFG.MainStrategies(3:end)
    y=nan(size(bins));
    for bi=1:numel(bins)
        idx=T.distance_bin==bins(bi)&T.geometry_type~="homogeneous";
        y(bi)=mean(abs(T.("delta_vs_local_"+strategy)(idx)),'omitnan');
    end
    plot(ax,1:numel(bins),y,'-o','LineWidth',1.2,'DisplayName',pretty(strategy));
end
set(ax,'XTick',1:numel(bins),'XTickLabel',labels); xlabel(ax,'Distance to interface (mm)');
ylabel(ax,'|SWS corrected - Local| (m/s)'); title(ax,'Correction magnitude vs distance'); style_axes(ax);
legend(ax,'Location','eastoutside','Box','off');
export_fig(fig,OUT,'test31_correction_vs_distance.png');
end

function plot_q_correction_summary(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
C=sortrows(T(T.condition_key==keys(1),:),{'map_iz','map_ix'});
methods=["q_nearest_highconf","q_median_highconf_by_region_local_structure","edgeaware_q_theory_structure"];
fig=figure('Color','w','Position',[30 30 1350 820],'Visible','off'); tl=tiledlayout(numel(methods),3,'TileSpacing','compact');
for mi=1:numel(methods)
    q=C.("q_corr_"+methods(mi)); maps={map_from_rows(C,C.q_local),map_from_rows(C,q),map_from_rows(C,q-C.q_local)};
    titles=["q original",pretty(methods(mi))+" q","q corrected - original"];
    for j=1:3
        ax=nexttile(tl); imagesc(ax,maps{j}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
        title(ax,titles(j),'Interpreter','none'); colorbar(ax);
    end
end
title(tl,C.condition_key(1),'Interpreter','none');
exportgraphics(fig,fullfile(OUT.figure_dir,'test31_q_correction_maps.png'),'Resolution',180); close(fig);
end

function plot_local_corrected_scatter(T,CFG,OUT)
strategies=["q_median_highconf_by_region_local_structure","edgeaware_sws_theory_structure","edgeaware_q_theory_structure","test30_theory_region_levels"];
fig=figure('Color','w','Position',[80 80 1250 880]); tl=tiledlayout(2,2,'TileSpacing','compact');
for si=1:numel(strategies)
    ax=nexttile(tl); idx=randperm(height(T),min(height(T),70000));
    scatter(ax,T.sws_local_baseline(idx),T.("sws_"+strategies(si))(idx),8,T.confidence(idx),'filled','MarkerFaceAlpha',.18);
    hold(ax,'on'); plot(ax,[.5 4],[.5 4],'k--'); xlabel(ax,'Local SWS'); ylabel(ax,'Corrected SWS');
    title(ax,pretty(strategies(si))); colorbar(ax); style_axes(ax);
end
export_fig(fig,OUT,'test31_local_vs_corrected_scatter.png');
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable');
if isempty(keys), keys=unique(T.condition_key,'stable'); end
if isempty(keys), return; end
plot_case_maps(T(T.condition_key==keys(1),:),CFG,fullfile(OUT.figure_dir,'test31_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable'); root=fullfile(OUT.figure_dir,'maps_by_condition_dx200um');
for i=1:numel(keys)
    C=T(T.condition_key==keys(i),:);
    folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end
    file=fullfile(folder,"test31__"+sanitize(keys(i))+".png");
    plot_case_maps(C,CFG,file);
end
end

function plot_case_maps(C,~,file)
C=sortrows(C,{'map_iz','map_ix'});
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_theory_baseline), ...
    map_from_rows(C,C.sws_local_baseline),map_from_rows(C,C.q_local), ...
    map_from_rows(C,C.local_structure),map_from_rows(C,C.theory_structure), ...
    map_from_rows(C,C.confidence),map_from_rows(C,C.low_confidence), ...
    map_from_rows(C,C.sws_q_nearest_highconf),map_from_rows(C,C.sws_sws_nearest_highconf), ...
    map_from_rows(C,C.sws_q_median_highconf_global), ...
    map_from_rows(C,C.sws_q_median_highconf_by_region_local_structure), ...
    map_from_rows(C,C.sws_sws_interp_by_region_local_structure), ...
    map_from_rows(C,C.sws_test30_theory_region_levels), ...
    map_from_rows(C,C.sws_edgeaware_sws_theory_structure), ...
    map_from_rows(C,C.sws_edgeaware_q_theory_structure), ...
    map_from_rows(C,C.abs_error_local_baseline), ...
    map_from_rows(C,C.abs_error_q_nearest_highconf), ...
    map_from_rows(C,C.abs_error_sws_nearest_highconf), ...
    map_from_rows(C,C.abs_error_q_median_highconf_by_region_local_structure), ...
    map_from_rows(C,C.abs_error_test30_theory_region_levels), ...
    map_from_rows(C,C.abs_error_edgeaware_sws_theory_structure), ...
    map_from_rows(C,C.abs_error_edgeaware_q_theory_structure), ...
    map_from_rows(C,C.donor_distance_edgeaware_sws_theory_structure_mm)};
titles=["True SWS","TheoryQDiscrete SWS","LocalOnly SWS","LocalOnly q", ...
    "Local structure","Theory/Test30 structure","Confidence","Low-confidence", ...
    "q nearest SWS","SWS nearest","q median global","q median by local region", ...
    "SWS local-region interp","Test30 region levels","Edge-aware SWS theory", ...
    "Edge-aware q theory","Local abs error","q nearest abs error", ...
    "SWS nearest abs error","q region abs error","Test30 abs error", ...
    "Edge-aware SWS abs error","Edge-aware q abs error","Donor distance (mm)"];
fig=figure('Color','w','Visible','off','Position',[20 20 1700 1100]); tl=tiledlayout(4,6,'TileSpacing','compact');
for i=1:numel(maps)
    ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
    title(ax,titles(i),'Interpreter','none','FontSize',8); colorbar(ax);
end
title(tl,C.condition_key(1),'Interpreter','none');
exportgraphics(fig,file,'Resolution',170); close(fig);
end

%% Console and validation

function print_interpretation(O,D,R,Q,CFG)
oper=O(ismember(O.strategy_name,CFG.Strategies),:);
[~,i]=min(oper.MAPE); best_global=oper(i,:);
soft=R(R.roi=="soft_core",:); hard=R(R.roi=="hard_core",:);
near=D(D.distance_bin=="0_1mm",:);
fprintf('\nPreliminary interpretation:\n');
fprintf('  Best global strategy: %s (MAPE %.2f%%, >20 %.2f%%).\n', ...
    best_global.strategy_name,best_global.MAPE,best_global.high_error_20_pct);
if ~isempty(soft)
    S=groupsummary(soft,'strategy_name','mean','MAPE'); [~,j]=min(S.mean_MAPE);
    fprintf('  Best soft-core ROI: %s (MAPE %.2f%%).\n',S.strategy_name(j),S.mean_MAPE(j));
end
if ~isempty(hard)
    H=groupsummary(hard,'strategy_name','mean','MAPE'); [~,j]=min(H.mean_MAPE);
    fprintf('  Best hard-core ROI: %s (MAPE %.2f%%).\n',H.strategy_name(j),H.mean_MAPE(j));
end
if ~isempty(near)
    N=groupsummary(near,'strategy_name','mean','MAPE'); [~,j]=min(N.mean_MAPE);
    fprintf('  Best near-interface distance bin: %s (MAPE %.2f%%).\n',N.strategy_name(j),N.mean_MAPE(j));
end
M2=oper(oper.MAPE==oper.MAPE,:); %#ok<NASGU>
if any(O.strategy_name=="local_baseline")
    L=O(O.strategy_name=="local_baseline",:);
    fprintf('  Local baseline: MAPE %.2f%%, homogeneous %.2f%%, strong mixed %.2f%%.\n', ...
        L.MAPE,L.MAPE_homogeneous,L.MAPE_strongly_mixed);
end
qNames=Q.strategy_name(Q.confidence_group=="low_confidence");
if ~isempty(qNames)
    swsQ=O(ismember(O.strategy_name,unique(qNames)),:);
    swsS=O(ismember(O.strategy_name,["sws_nearest_highconf","sws_interp_by_region_local_structure","edgeaware_sws_theory_structure"]),:);
    fprintf('  Best q-interpolation MAPE: %.2f%%; best SWS-interpolation MAPE: %.2f%%.\n', ...
        min(swsQ.MAPE,[],'omitnan'),min(swsS.MAPE,[],'omitnan'));
end
E=O(O.strategy_name=="edgeaware_sws_theory_structure",:);
N=O(O.strategy_name=="sws_nearest_highconf",:);
if ~isempty(E)&&~isempty(N)
    fprintf('  Theory edge-aware SWS vs no-structure nearest: %.2f%% vs %.2f%% global MAPE.\n',E.MAPE,N.MAPE);
end
fprintf('  Recommendation: compare the winning simple method against Test30 region levels; keep the confidence map and donor-distance map as reliability outputs.\n');
end

function validate_test31(T,T30,spectral_dir,key,CFG)
X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
SPEC=load_spectral_cache(spectral_dir,key,X);
R=evaluate_condition(X,SPEC,T30,CFG);
assert(height(R)==height(X));
assert(all(isfinite(R.sws_edgeaware_q_theory_structure)));
assert(all(isfinite(R.q_corr_edgeaware_q_theory_structure)));
assert(all(R.low_confidence==(R.confidence<CFG.ConfidenceThreshold)));
fprintf('  q/SWS fill, q-to-SWS conversion, structures, and output dimensions passed for %s.\n',key);
end

%% Generic helpers

function tf=env_true(name,default)
value=strtrim(getenv(name));
if isempty(value), tf=default; else, tf=any(strcmpi(value,{'true','1','yes','on'})); end
end

function value=env_double(name,default)
raw=strtrim(getenv(name));
if isempty(raw)
    value=default;
else
    value=str2double(raw);
    assert(isfinite(value),'%s must be numeric.',name);
end
end

function write_config_json(CFG,file)
C=CFG; for f=string(fieldnames(C))', if isstring(C.(f)), C.(f)=cellstr(C.(f)); end, end
fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char'); clear cleanup;
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
missing=setdiff(A.Properties.VariableNames,B.Properties.VariableNames);
for m=missing, B.(m)=missing_value_like(A.(m),height(B)); end
extra=setdiff(B.Properties.VariableNames,A.Properties.VariableNames);
for e=extra, A.(e)=missing_value_like(B.(e),height(A)); end
B=B(:,A.Properties.VariableNames); T=[A;B];
end

function value=missing_value_like(example,n)
if isstring(example), value=strings(n,1);
elseif iscellstr(example), value=repmat({''},n,1);
elseif islogical(example), value=false(n,1);
else, value=nan(n,1);
end
end

function label=pretty(name)
label=replace(string(name),["_","q nearest highconf","sws nearest highconf"],[" ","q nearest","SWS nearest"]);
label=replace(label,["local baseline","theory baseline","q median highconf by region local structure", ...
    "sws interp by region local structure","test30 theory region levels", ...
    "edgeaware sws theory structure","edgeaware q theory structure", ...
    "edgeaware sws local structure","edgeaware q local structure", ...
    "q mean highconf global","q median highconf global"], ...
    ["Local","TheoryQDiscrete","q median by local region", ...
    "SWS interp local region","Test30 region levels", ...
    "Edge-aware SWS Theory","Edge-aware q Theory", ...
    "Edge-aware SWS Local","Edge-aware q Local", ...
    "q mean global","q median global"]);
end

function value=sanitize(value)
value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_');
end

function style_axes(ax)
grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9;
end

function export_fig(fig,OUT,name)
exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',210); close(fig);
end

function y=ternary(condition,a,b)
if condition, y=a; else, y=b; end
end
