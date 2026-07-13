%% analyze_test_28_two_component_edge_masked_req.m
% Test 28: confidence-gated two-component and edge-masked REQ.
%
% This test does not train or modify a model. It reuses frozen Test 27 SWS
% predictions, Test 22 confidence, Test 26 radial spectra, and Test 22 field
% caches. Operational corrections use only Local/Hybrid SWS, confidence,
% patch coordinates, and measured spectra. True SWS, material side, patch
% purity, and interface distance are attached only after inference.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST28_MODE          = quick | full
%   ADAPTIVE_REQ_TEST28_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST28_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST28_EDGE_MASKED   = true | false

clear; clc; close all;
format compact;

%% Configuration

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.15, ...
    'defaultAxesLabelFontSizeMultiplier',1.05);

CFG=struct();
mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST28_MODE'))));
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]), ...
    'ADAPTIVE_REQ_TEST28_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST28_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST28_SAVE_ALL_MAPS',true);
CFG.EnableEdgeMasked=env_true('ADAPTIVE_REQ_TEST28_EDGE_MASKED',true);
CFG.Version=2;
CFG.ConfidenceThreshold=0.80;
CFG.PhysicalSwsRange=[0.5 10];
CFG.PseudoRegionMinSeparation=0.25;
CFG.MinPeakHeightRatio=0.10;
CFG.MinPeakSeparationFraction=0.12;
CFG.MaxValleyRatio=0.82;
CFG.MinComponentWeight=0.08;
CFG.MinMixtureScore=0.004;
CFG.MaskTransitionSws=0.08;
CFG.MaskFloor=0.03;
CFG.MaskEffectiveFraction=[0.20 0.92];
CFG.MaskedWidthRatioMax=0.95;
CFG.MaskedScoreImprovementMin=0.005;
CFG.MaxRelativeCorrection=0.50;
CFG.RandomSeed=28001;
CFG.Dx=0.2e-3;
CFG.Geometries=["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"];
if CFG.QuickMode
    CFG.Frequencies=500; CFG.M=[2 3];
    CFG.Regimes=["directional_2D","diffuse_3D"];
else
    CFG.Frequencies=[300 400 500 600]; CFG.M=[2 3 4];
    CFG.Regimes=["directional_2D","diffuse_2D", ...
        "partial_3D","diffuse_3D"];
end
CFG.Strategies=["hybrid_baseline","local_baseline", ...
    "confidence_switch_c080","two_component_local_c080", ...
    "two_component_switch_c080","edge_masked_local_c080", ...
    "edge_masked_switch_c080","two_component_then_edge_masked_c080", ...
    "two_component_then_edge_masked_switch_c080"];
CFG.SelectedStrategy="two_component_then_edge_masked_switch_c080";

TEST27=fullfile(root_dir,'outputs','test_27_adaptive_window_edge_aware');
TEST26=fullfile(root_dir,'outputs','test_26_confidence_gated_corrections');
if CFG.QuickMode
    quick27=fullfile(TEST27,'quick'); quick26=fullfile(TEST26,'quick');
    using_quick27=false;
    if exist(fullfile(quick27,'tables','test27_patch_level_results.csv'),'file')==2
        TEST27=quick27; using_quick27=true;
    end
    if using_quick27&&exist(fullfile(quick26,'data','spectral_checkpoints'),'dir')==7
        TEST26=quick26;
    end
end
SOURCE_TABLE=fullfile(TEST27,'tables','test27_patch_level_results.csv');
SPECTRAL_DIR=fullfile(TEST26,'data','spectral_checkpoints');
FIELD_DIR=fullfile(root_dir,'outputs','test_22_confidence_external_validation', ...
    'analysis','data','field_cache');
OUT=make_output_dirs(root_dir,CFG);
write_config_json(CFG,fullfile(OUT.root_dir,'test28_configuration.json'));

fprintf('\nTest 28: confidence-gated two-component and edge-masked REQ\n');
fprintf('No training. Oracle variables are evaluation-only.\n');
fprintf('Mode: %s | edge-masked candidate: %d\n',CFG.RunMode,CFG.EnableEdgeMasked);
assert(exist(SOURCE_TABLE,'file')==2,'Test 27 patch table missing: %s',SOURCE_TABLE);
assert(exist(SPECTRAL_DIR,'dir')==7,'Test 26 spectral cache missing: %s',SPECTRAL_DIR);

T=readtable(SOURCE_TABLE,'TextType','string');
T=filter_design(T,CFG);
keys=unique(T.condition_key,'stable');
assert(~isempty(keys),'No Test 27 conditions match Test 28 configuration.');
fprintf('Matched %d cached conditions and %d patches.\n',numel(keys),height(T));

if CFG.ValidateOnly
    validate_test28(T,keys(1),SPECTRAL_DIR,FIELD_DIR,CFG);
    fprintf('Test 28 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Condition-level operational corrections

parts=cell(numel(keys),1);
for ci=1:numel(keys)
    timer=tic; key=keys(ci);
    checkpoint=fullfile(OUT.condition_dir,"test28__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        S=load(checkpoint,'RESULT');
        if isfield(S,'RESULT')&&isfield(S.RESULT,'version')&& ...
                S.RESULT.version==CFG.Version
            parts{ci}=S.RESULT.T_patch;
            fprintf('[%d/%d] Reused %s (%d patches).\n', ...
                ci,numel(keys),key,height(parts{ci}));
            continue;
        end
    end
    X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
    SPEC=load_spectral_cache(SPECTRAL_DIR,key,X);
    operational=build_operational_inputs(X,SPEC,CFG);
    field=struct();
    if CFG.EnableEdgeMasked
        field=load_field_cache(FIELD_DIR,X(1,:));
    end
    [pred,diag]=evaluate_condition(operational,SPEC,field,CFG);
    assert_operational_policy(pred,diag,operational,CFG);
    T_patch=attach_diagnostics(X,pred,diag,CFG);
    RESULT=struct('version',CFG.Version,'T_patch',T_patch);
    save(checkpoint,'RESULT','-v7.3');
    parts{ci}=T_patch;
    fprintf(['[%d/%d] %s: low confidence %d, mixture %d, ', ...
        'edge-mask %d in %.1f s.\n'],ci,numel(keys),key, ...
        sum(diag.low_confidence),sum(diag.mixture_accepted), ...
        sum(diag.edge_accepted),toc(timer));
end
T_patch=vertcat(parts{:}); clear parts T;

%% Tables

T_overall=summarize_predictions(T_patch,"strategy_name",CFG);
T_geometry=summarize_predictions(T_patch,["strategy_name","geometry"],CFG);
T_regime=summarize_predictions(T_patch,["strategy_name","field_regime"],CFG);
T_frequency_M=summarize_predictions(T_patch,["strategy_name","f0","M"],CFG);
T_purity=summarize_predictions(T_patch, ...
    ["strategy_name","geometry","purity_bin"],CFG);
T_distance=summarize_predictions(T_patch, ...
    ["strategy_name","geometry","distance_bin"],CFG);
T_acceptance=summarize_acceptance(T_patch);
T_best=best_candidates(T_overall);

writetable(T_patch,fullfile(OUT.table_dir,'test28_patch_level_results.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test28_strategy_summary_overall.csv'));
writetable(T_geometry,fullfile(OUT.table_dir,'test28_strategy_summary_by_geometry.csv'));
writetable(T_regime,fullfile(OUT.table_dir,'test28_strategy_summary_by_regime.csv'));
writetable(T_frequency_M,fullfile(OUT.table_dir,'test28_strategy_summary_by_frequency_M.csv'));
writetable(T_purity,fullfile(OUT.table_dir,'test28_strategy_summary_by_purity.csv'));
writetable(T_distance,fullfile(OUT.table_dir,'test28_strategy_summary_by_distance.csv'));
writetable(T_acceptance,fullfile(OUT.table_dir,'test28_correction_acceptance.csv'));
writetable(T_best,fullfile(OUT.table_dir,'test28_best_strategy_candidates.csv'));
save(fullfile(OUT.data_dir,'test28_compact_results.mat'),'T_patch','CFG','-v7.3');

%% Figures and summary

plot_strategy_summary(T_patch,CFG,OUT);
plot_error_by_purity(T_patch,CFG,OUT);
plot_error_by_distance(T_patch,CFG,OUT);
plot_component_diagnostics(T_patch,OUT);
plot_representative_spectra(T_patch,SPECTRAL_DIR,OUT);
plot_representative_maps(T_patch,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T_patch,CFG,OUT); end
print_interpretation(T_overall,T_acceptance,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 28 complete.\n', ...
    OUT.table_dir,OUT.figure_dir);

%% Configuration and cache helpers

function OUT=make_output_dirs(root_dir,CFG)
OUT.root_dir=fullfile(root_dir,'outputs', ...
    'test_28_two_component_edge_masked_req');
if CFG.QuickMode, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
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
if isempty(value), tf=default_value;
else, tf=any(strcmpi(value,{'true','1','yes','on'})); end
end

function write_config_json(CFG,file)
C=CFG;
for f=string(fieldnames(C))'
    if isstring(C.(f)), C.(f)=cellstr(C.(f)); end
end
fid=fopen(file,'w'); assert(fid>=0,'Cannot write %s.',file);
cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char');
end

function T=filter_design(T,CFG)
T=T(abs(T.dx-CFG.Dx)<1e-12 & ismember(T.geometry,CFG.Geometries) & ...
    ismember(T.f0,CFG.Frequencies) & ismember(T.M,CFG.M) & ...
    ismember(T.field_regime,CFG.Regimes),:);
end

function SPEC=load_spectral_cache(folder,key,X)
file=fullfile(folder,"spectra__"+sanitize(key)+".mat");
assert(exist(file,'file')==2,'Spectral cache missing: %s',file);
S=load(file,'SPEC'); SPEC=S.SPEC;
assert(numel(SPEC.curves)==height(X), ...
    'Spectral cache row count mismatch for %s.',key);
cx=round(X.x/X.dx(1))+1; cz=round(X.z/X.dz(1))+1;
assert(isequal(double(SPEC.cx(:)),double(cx(:)))&& ...
    isequal(double(SPEC.cz(:)),double(cz(:))), ...
    'Spectral cache coordinate mismatch for %s.',key);
end

function F=load_field_cache(folder,row)
key=sprintf('%s__f%g__%s__dx%gum',row.geometry,row.f0, ...
    lower(string(row.field_regime)),round(1e6*row.dx));
file=fullfile(folder,"field__"+sanitize(key)+".mat");
assert(exist(file,'file')==2,'Field cache missing: %s',file);
F=load(file,'sim','cfg_sim');
end

%% Operational inference

function O=build_operational_inputs(X,SPEC,CFG)
O=struct(); O.meta=X(:,{'condition_key','geometry','geometry_type', ...
    'field_regime','f0','M','dx','dz','map_iz','map_ix','x','z'});
O.local=double(X.sws_local_baseline);
O.hybrid=double(X.sws_hybrid_baseline);
O.switch=double(X.sws_confidence_switch_local_hybrid_c080);
O.confidence=double(X.confidence);
[O.pseudo_region,O.pseudo_threshold,O.pseudo_separation]= ...
    pseudo_segment(O.local,O.hybrid,X,CFG);
n=height(X); O.decomposition=repmat( ...
    adaptive_req.analysis.decompose_radial_spectrum(struct()),n,1);
f0=X.f0(1); krange=sort(2*pi*f0./CFG.PhysicalSwsRange);
for i=1:n
    O.decomposition(i)=adaptive_req.analysis.decompose_radial_spectrum( ...
        SPEC.curves{i},'KRange',krange, ...
        'MinPeakHeightRatio',CFG.MinPeakHeightRatio, ...
        'MinSeparationFraction',CFG.MinPeakSeparationFraction, ...
        'MaxValleyRatio',CFG.MaxValleyRatio, ...
        'MinComponentWeight',CFG.MinComponentWeight);
end
end

function [region,threshold,separation]=pseudo_segment(local,hybrid,meta,CFG)
seed=.75*local+.25*hybrid;
[A,iz,ix]=map_from_rows(meta,seed);
A=medfilt2(A,[3 3],'symmetric'); x=A(isfinite(A));
if isempty(x), region=ones(size(seed)); threshold=NaN; separation=0; return; end
c=[prctile(x,25) prctile(x,75)];
for iter=1:20
    first=abs(x-c(1))<=abs(x-c(2));
    if any(first), c(1)=median(x(first)); end
    if any(~first), c(2)=median(x(~first)); end
end
c=sort(c); separation=diff(c); threshold=mean(c);
L=ones(size(A));
if separation>=CFG.PseudoRegionMinSeparation, L(A>threshold)=2; end
region=L(sub2ind(size(L),iz,ix));
end

function [pred,D]=evaluate_condition(O,SPEC,F,CFG)
n=numel(O.local); low=O.confidence<CFG.ConfidenceThreshold;
mix=O.local; mix_switch=O.switch; mix_accept=false(n,1);
selected_k=nan(n,1); low_k=nan(n,1); high_k=nan(n,1);
score=nan(n,1); valley=nan(n,1); detected=false(n,1);
for i=1:n
    d=O.decomposition(i); low_k(i)=d.low_k; high_k(i)=d.high_k;
    score(i)=d.two_component_score; valley(i)=d.valley_ratio;
    detected(i)=d.two_component_detected;
    if ~(low(i)&&d.two_component_detected&& ...
            d.two_component_score>=CFG.MinMixtureScore&& ...
            O.pseudo_separation>=CFG.PseudoRegionMinSeparation)
        continue;
    end
    if O.pseudo_region(i)==1, k=d.high_k; else, k=d.low_k; end
    candidate=2*pi*O.meta.f0(i)/k;
    if operational_candidate_ok(candidate,O.local(i),CFG)
        mix(i)=candidate; mix_switch(i)=candidate;
        selected_k(i)=k; mix_accept(i)=true;
    end
end

edge=O.local; edge_switch=O.switch; edge_accept=false(n,1); edge_sws=nan(n,1);
edge_effective=nan(n,1); edge_width_ratio=nan(n,1);
edge_score_delta=nan(n,1);
if CFG.EnableEdgeMasked
    [edge_sws,edge_accept,edge_effective,edge_width_ratio,edge_score_delta]= ...
        edge_masked_candidates(O,SPEC,F,low,CFG);
    edge(edge_accept)=edge_sws(edge_accept);
    edge_switch(edge_accept)=edge_sws(edge_accept);
end
combined=O.local; source=zeros(n,1,'uint8');
combined(mix_accept)=mix(mix_accept); source(mix_accept)=1;
use_edge=edge_accept&~mix_accept;
combined(use_edge)=edge_sws(use_edge); source(use_edge)=2;
combined_switch=O.switch;
combined_switch(mix_accept)=mix_switch(mix_accept);
combined_switch(use_edge)=edge_sws(use_edge);

pred=struct('hybrid_baseline',O.hybrid,'local_baseline',O.local, ...
    'confidence_switch_c080',O.switch, ...
    'two_component_local_c080',mix, ...
    'two_component_switch_c080',mix_switch, ...
    'edge_masked_local_c080',edge, ...
    'edge_masked_switch_c080',edge_switch, ...
    'two_component_then_edge_masked_c080',combined, ...
    'two_component_then_edge_masked_switch_c080',combined_switch);
D=struct('low_confidence',low,'pseudo_region',O.pseudo_region, ...
    'pseudo_threshold',repmat(O.pseudo_threshold,n,1), ...
    'pseudo_separation',repmat(O.pseudo_separation,n,1), ...
    'two_component_detected',detected,'mixture_accepted',mix_accept, ...
    'edge_accepted',edge_accept,'selected_k',selected_k, ...
    'component_low_k',low_k,'component_high_k',high_k, ...
    'two_component_score',score,'valley_ratio',valley, ...
    'edge_candidate_sws',edge_sws,'edge_effective_fraction',edge_effective, ...
    'edge_width_ratio',edge_width_ratio,'edge_score_delta',edge_score_delta, ...
    'correction_source',source);
end

function [candidate,accepted,effective,width_ratio,score_delta]= ...
        edge_masked_candidates(O,SPEC,F,low,CFG)
n=numel(O.local); candidate=nan(n,1); accepted=false(n,1);
effective=nan(n,1); width_ratio=nan(n,1); score_delta=nan(n,1);
if O.pseudo_separation<CFG.PseudoRegionMinSeparation, return; end

seed=.75*O.local+.25*O.hybrid;
ok=isfinite(seed); prior=scatteredInterpolant(O.meta.x(ok),O.meta.z(ok), ...
    seed(ok),'linear','nearest');
feat=adaptive_req.config.default_feature_config('M',O.meta.M(1), ...
    'cs_guess',3,'gamma_win',1,'pad_factor',1);
[req,feat]=adaptive_req.config.default_req_config(F.cfg_sim,feat, ...
    'Nbins','auto','Nbins_auto_oversample',1,'Nbins_min',16, ...
    'smooth_sigma',1);
half=floor(feat.win_size/2); cfg_patch=F.cfg_sim; cfg_patch.UseParfor=false;
idx=find(low);
for jj=1:numel(idx)
    i=idx(jj); cx=SPEC.cx(i); cz=SPEC.cz(i);
    zi=(cz-half):(cz+half); xi=(cx-half):(cx+half);
    if min(zi)<1||max(zi)>size(F.sim.Uxz,1)|| ...
            min(xi)<1||max(xi)>size(F.sim.Uxz,2), continue; end
    [XX,ZZ]=meshgrid((xi-1)*F.cfg_sim.dx,(zi-1)*F.cfg_sim.dz);
    local_prior=prior(XX,ZZ); scale=max(CFG.MaskTransitionSws,eps);
    if O.pseudo_region(i)==1
        weight=1./(1+exp((local_prior-O.pseudo_threshold)/scale));
    else
        weight=1./(1+exp((O.pseudo_threshold-local_prior)/scale));
    end
    weight=CFG.MaskFloor+(1-CFG.MaskFloor)*weight;
    effective(i)=mean(weight(:));
    if effective(i)<CFG.MaskEffectiveFraction(1)|| ...
            effective(i)>CFG.MaskEffectiveFraction(2), continue; end
    patch=F.sim.Uxz(zi,xi).*sqrt(weight);
    try
        [~,masked_curve]=adaptive_req.quantile.compute_quantile_from_patch( ...
            patch,cfg_patch,req);
        q_local=q_at_sws(SPEC.curves{i},O.local(i),O.meta.f0(i));
        value=adaptive_req.quantile.quantile_to_cs( ...
            masked_curve,q_local,O.meta.f0(i));
        original=O.decomposition(i);
        krange=sort(2*pi*O.meta.f0(i)./CFG.PhysicalSwsRange);
        masked=adaptive_req.analysis.decompose_radial_spectrum(masked_curve, ...
            'KRange',krange,'MinPeakHeightRatio',CFG.MinPeakHeightRatio, ...
            'MinSeparationFraction',CFG.MinPeakSeparationFraction, ...
            'MaxValleyRatio',CFG.MaxValleyRatio, ...
            'MinComponentWeight',CFG.MinComponentWeight);
        width_ratio(i)=masked.radial_width_10_90/max(original.radial_width_10_90,eps);
        score_delta(i)=original.two_component_score-masked.two_component_score;
        quality=width_ratio(i)<=CFG.MaskedWidthRatioMax|| ...
            score_delta(i)>=CFG.MaskedScoreImprovementMin;
        if quality&&operational_candidate_ok(value,O.local(i),CFG)
            candidate(i)=value; accepted(i)=true;
        end
    catch
        % Invalid masked spectra fall back to Local without changing output.
    end
end
end

function q=q_at_sws(curve,sws,f0)
k=double(curve.k_cent(:)); E=double(curve.Ecum(:));
valid=isfinite(k)&isfinite(E); [k,ia]=unique(k(valid),'stable');
E=E(valid); E=E(ia); target=2*pi*f0/sws;
if numel(k)<2, q=NaN; else, q=interp1(k,E,target,'linear','extrap'); end
q=min(max(q,0),1);
end

function tf=operational_candidate_ok(value,fallback,CFG)
tf=isfinite(value)&&value>=CFG.PhysicalSwsRange(1)&& ...
    value<=CFG.PhysicalSwsRange(2)&& ...
    abs(value-fallback)/max(fallback,eps)<=CFG.MaxRelativeCorrection;
end

function assert_operational_policy(pred,D,O,CFG)
for name=CFG.Strategies
    value=pred.(name);
    assert(numel(value)==height(O.meta)&&all(isfinite(value))&& ...
        all(value>=CFG.PhysicalSwsRange(1)&value<=CFG.PhysicalSwsRange(2)), ...
        'Invalid output from %s.',name);
end
assert(~any(D.mixture_accepted&~D.low_confidence));
assert(~any(D.edge_accepted&~D.low_confidence));
operational_inputs=["local","hybrid","confidence","radial_spectrum", ...
    "coordinates","pseudo_region"];
oracle=["true_SWS","patch_purity","material_side", ...
    "distance_to_interface","q_true","q_theory"];
assert(isempty(intersect(lower(operational_inputs),lower(oracle))));
end

%% Evaluation-only outputs and summaries

function R=attach_diagnostics(X,pred,D,CFG)
keep={'condition_key','geometry','geometry_type','field_regime','f0','M', ...
    'dx','dz','map_iz','map_ix','x','z','true_SWS','confidence', ...
    'patch_purity','material_side','purity_bin','distance_bin', ...
    'distance_to_interface_mm'};
R=X(:,keep);
for name=CFG.Strategies, R.("sws_"+name)=pred.(name); end
for f=string(fieldnames(D))', R.(f)=D.(f); end
end

function S=summarize_predictions(T,groups,CFG)
S=table();
for strategy=CFG.Strategies
    pred=T.("sws_"+strategy); error=100*(pred-T.true_SWS)./T.true_SWS;
    X=T(:,cellstr(setdiff(groups,"strategy_name",'stable')));
    X.strategy_name=repmat(strategy,height(T),1); X.error=error;
    X.abs_error=abs(error); X.high10=X.abs_error>10; X.high20=X.abs_error>20;
    if isempty(groups)||isequal(groups,"strategy_name"), G=ones(height(X),1);
    else, G=findgroups(X(:,cellstr(setdiff(groups,"strategy_name",'stable')))); end
    for id=unique(G,'stable')'
        Y=X(G==id,:); row=Y(1,cellstr(groups)); row.N=height(Y);
        row.MAPE=mean(Y.abs_error); row.mean_signed_error=mean(Y.error);
        row.median_signed_error=median(Y.error);
        row.high_error_10_pct=100*mean(Y.high10);
        row.high_error_20_pct=100*mean(Y.high20);
        row.underestimate_pct=100*mean(Y.error<0);
        row.overestimate_pct=100*mean(Y.error>0);
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            idx=T.purity_bin==cls;
            if ~isempty(groups)&&~isequal(groups,"strategy_name")
                basevars=setdiff(groups,"strategy_name",'stable');
                for v=basevars, idx=idx&T.(v)==Y.(v)(1); end
            end
            token=char(cls); row.("N_"+token)=sum(idx);
            row.("MAPE_"+token)=mean(abs(100*(pred(idx)-T.true_SWS(idx))./T.true_SWS(idx)),'omitnan');
        end
        S=concat_tables(S,row);
    end
end
end

function S=summarize_acceptance(T)
groups=["geometry","field_regime","f0","M"];
[G,S]=findgroups(T(:,cellstr(groups)));
S.N=splitapply(@numel,T.confidence,G);
S.low_confidence_pct=100*splitapply(@mean,T.low_confidence,G);
S.two_component_detected_pct=100*splitapply(@mean,T.two_component_detected,G);
S.mixture_accepted_pct=100*splitapply(@mean,T.mixture_accepted,G);
S.edge_accepted_pct=100*splitapply(@mean,T.edge_accepted,G);
S.mean_two_component_score=splitapply(@(x)mean(x,'omitnan'),T.two_component_score,G);
S.mean_mask_effective_fraction=splitapply(@(x)mean(x,'omitnan'),T.edge_effective_fraction,G);
end

function B=best_candidates(S)
[~,ig]=min(S.MAPE); mixed=mean([S.MAPE_moderately_mixed S.MAPE_strongly_mixed],2,'omitnan');
[~,im]=min(mixed); pure=mean([S.MAPE_pure S.MAPE_near_pure],2,'omitnan'); [~,ip]=min(pure);
base=S(S.strategy_name=="local_baseline",:);
penalty=max(S.MAPE_homogeneous-base.MAPE_homogeneous-.5,0)*3+ ...
    max(pure-mean([base.MAPE_pure base.MAPE_near_pure],2,'omitnan'),0);
[~,it]=min(S.MAPE+penalty);
B=table(S.strategy_name(ig),S.strategy_name(im),S.strategy_name(ip), ...
    S.strategy_name(it),'VariableNames',{'best_global','best_mixed', ...
    'best_pure','best_tradeoff'});
end

%% Figures

function plot_strategy_summary(T,CFG,OUT)
M=nan(numel(CFG.Strategies),3); H=M;
for si=1:numel(CFG.Strategies)
    pred=T.("sws_"+CFG.Strategies(si)); err=abs(100*(pred-T.true_SWS)./T.true_SWS);
    masks={T.purity_bin=="homogeneous",ismember(T.purity_bin,["pure","near_pure"]), ...
        ismember(T.purity_bin,["moderately_mixed","strongly_mixed"])};
    for j=1:3, M(si,j)=mean(err(masks{j}),'omitnan'); H(si,j)=100*mean(err(masks{j})>20,'omitnan'); end
end
fig=figure('Color','w','Position',[80 60 1300 720]); tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); barh(ax,M); xlabel(ax,'MAPE (%)'); title(ax,'Mean absolute error'); style_axes(ax);
ax=nexttile(tl); barh(ax,H); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style_axes(ax);
for ax=findall(fig,'Type','axes')'
    set(ax,'YTick',1:numel(CFG.Strategies),'YTickLabel',pretty(CFG.Strategies), ...
        'TickLabelInterpreter','none','FontSize',8);
end
lg=legend(["Homogeneous","Pure / near-pure","Mixed"],'Location','northoutside', ...
    'Orientation','horizontal','Box','off'); lg.FontSize=8;
export_fig(fig,fullfile(OUT.figure_dir,'test28_strategy_summary.png'));
end

function plot_error_by_purity(T,CFG,OUT)
bins=["strongly_mixed","moderately_mixed","near_pure","pure"];
plot_binned_error(T,CFG,bins,"purity_bin",replace(bins,"_"," "), ...
    'Patch purity class',fullfile(OUT.figure_dir,'test28_error_by_purity.png'));
end

function plot_error_by_distance(T,CFG,OUT)
bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
plot_binned_error(T,CFG,bins,"distance_bin",["0-1","1-2","2-4","4-8",">8"], ...
    'Distance to interface (mm)',fullfile(OUT.figure_dir,'test28_error_by_distance.png'));
end

function plot_binned_error(T,CFG,bins,var,labels,xlabel_text,file)
T=T(T.geometry_type~="homogeneous",:); cases=["bilayer_2_3","circular_inclusion_2_3"];
selected=["local_baseline","confidence_switch_c080", ...
    "two_component_local_c080","edge_masked_switch_c080",CFG.SelectedStrategy];
fig=figure('Color','w','Position',[80 80 1180 480]); tl=tiledlayout(1,2,'TileSpacing','compact');
for gi=1:2
    ax=nexttile(tl); hold(ax,'on');
    for strategy=selected
        pred=T.("sws_"+strategy); err=abs(100*(pred-T.true_SWS)./T.true_SWS); y=nan(size(bins));
        for bi=1:numel(bins), idx=T.geometry==cases(gi)&T.(var)==bins(bi); y(bi)=mean(err(idx),'omitnan'); end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.5,'DisplayName',pretty(strategy));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',labels); xlabel(ax,xlabel_text); ylabel(ax,'MAPE (%)');
    title(ax,replace(cases(gi),"_"," ")); style_axes(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off');
export_fig(fig,file);
end

function plot_component_diagnostics(T,OUT)
fig=figure('Color','w','Position',[80 80 1150 470]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); scatter(ax,T.confidence,T.two_component_score,7,T.patch_purity, ...
    'filled','MarkerFaceAlpha',.18); xlabel(ax,'Confidence'); ylabel(ax,'Two-component score'); colorbar(ax); style_axes(ax);
classes=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"];
rate=nan(numel(classes),3);
for i=1:numel(classes)
    idx=T.purity_bin==classes(i); rate(i,:)=[100*mean(T.two_component_detected(idx)), ...
        100*mean(T.mixture_accepted(idx)),100*mean(T.edge_accepted(idx))];
end
ax=nexttile(tl); bar(ax,rate); set(ax,'XTick',1:numel(classes), ...
    'XTickLabel',replace(classes,"_"," "),'XTickLabelRotation',20); ylabel(ax,'Patches (%)');
legend(ax,["Two components","Mixture accepted","Edge mask accepted"], ...
    'Box','off','Location','northoutside','Orientation','horizontal','FontSize',8); style_axes(ax);
export_fig(fig,fullfile(OUT.figure_dir,'test28_component_diagnostics.png'));
end

function plot_representative_spectra(T,folder,OUT)
X=T(T.two_component_detected,:); if isempty(X), return; end
[~,order]=sort(X.two_component_score,'descend'); X=X(order(1:min(6,height(X))),:);
fig=figure('Color','w','Position',[50 50 1300 720]); tl=tiledlayout(2,3,'TileSpacing','compact');
for i=1:height(X)
    S=load(fullfile(folder,"spectra__"+sanitize(X.condition_key(i))+".mat"),'SPEC');
    key=[X.map_iz(i) X.map_ix(i)]; C=T(T.condition_key==X.condition_key(i),:);
    C=sortrows(C,{'map_iz','map_ix'}); j=find(C.map_iz==key(1)&C.map_ix==key(2),1);
    curve=S.SPEC.curves{j}; ax=nexttile(tl); plot(ax,curve.k_cent,curve.Srad/max(curve.Srad),'LineWidth',1.2); hold(ax,'on');
    xline(ax,X.component_low_k(i),'b--','low-k'); xline(ax,X.component_high_k(i),'r--','high-k');
    xlabel(ax,'k (rad/m)'); ylabel(ax,'Normalized radial power'); title(ax, ...
        sprintf('%s | purity %.2f | score %.3g',X.geometry(i),X.patch_purity(i),X.two_component_score(i)), ...
        'Interpreter','none','FontWeight','normal'); style_axes(ax);
end
export_fig(fig,fullfile(OUT.figure_dir,'test28_representative_two_component_spectra.png'));
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable');
if isempty(keys), return; end
plot_map_panel(T(T.condition_key==keys(1),:),CFG, ...
    fullfile(OUT.figure_dir,'test28_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable'); root=fullfile(OUT.figure_dir,'maps_by_condition_dx200um');
for i=1:numel(keys)
    C=T(T.condition_key==keys(i),:); folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end
    file=fullfile(folder,"test28_maps__"+sanitize(keys(i))+".png");
    if exist(file,'file')~=2, plot_map_panel(C,CFG,file); end
    if mod(i,25)==0||i==numel(keys), fprintf('  Test 28 maps: %d/%d.\n',i,numel(keys)); end
end
end

function plot_map_panel(C,CFG,file)
C=sortrows(C,{'map_iz','map_ix'}); selected=C.("sws_"+CFG.SelectedStrategy);
e0=abs(100*(C.sws_local_baseline-C.true_SWS)./C.true_SWS);
e1=abs(100*(selected-C.true_SWS)./C.true_SWS);
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_hybrid_baseline), ...
    map_from_rows(C,C.sws_local_baseline),map_from_rows(C,selected), ...
    map_from_rows(C,e0),map_from_rows(C,e1),map_from_rows(C,C.confidence), ...
    map_from_rows(C,C.two_component_score),map_from_rows(C,C.mixture_accepted), ...
    map_from_rows(C,C.edge_accepted),map_from_rows(C,C.correction_source), ...
    map_from_rows(C,selected-C.sws_local_baseline)};
titles=["True SWS","Hybrid","Local","Selected","Local error","Selected error", ...
    "Confidence","Two-component score","Mixture accepted","Edge accepted", ...
    "Correction source","Selected - Local"];
fig=figure('Color','w','Visible','off','Position',[20 40 1500 750]); tl=tiledlayout(2,6,'TileSpacing','compact','Padding','compact');
for i=1:12
    ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i)); colorbar(ax);
    if i<=4, clim(ax,[1.5 3.5]); elseif i==5||i==6, clim(ax,[0 50]); elseif i==7, clim(ax,[0 1]); end
end
title(tl,C.condition_key(1),'Interpreter','none'); exportgraphics(fig,file,'Resolution',170); close(fig);
end

function print_interpretation(S,A,CFG)
[~,ig]=min(S.MAPE); mixed=mean([S.MAPE_moderately_mixed S.MAPE_strongly_mixed],2,'omitnan'); [~,im]=min(mixed);
base=S(S.strategy_name=="local_baseline",:); selected=S(S.strategy_name==CFG.SelectedStrategy,:);
fprintf('\nInterpretation:\n');
fprintf('  Best global: %s (MAPE %.2f%%).\n',S.strategy_name(ig),S.MAPE(ig));
fprintf('  Best mixed: %s (MAPE %.2f%%).\n',S.strategy_name(im),mixed(im));
fprintf('  Selected homogeneous delta vs Local: %+.2f points.\n',selected.MAPE_homogeneous-base.MAPE_homogeneous);
fprintf('  Two components detected in %.1f%% of patches.\n',sum(A.N.*A.two_component_detected_pct)/sum(A.N));
fprintf('  Mixture accepted in %.1f%%; edge mask accepted in %.1f%%.\n', ...
    sum(A.N.*A.mixture_accepted_pct)/sum(A.N),sum(A.N.*A.edge_accepted_pct)/sum(A.N));
fprintf('  Candidate continuous output: %s + frozen confidence map.\n',CFG.SelectedStrategy);
end

%% Validation and generic helpers

function validate_test28(T,key,spectral_dir,field_dir,CFG)
k=linspace(100,2500,180)';
one=struct('k_cent',k,'Srad',exp(-0.5*((k-900)/90).^2));
two=struct('k_cent',k,'Srad',exp(-0.5*((k-750)/65).^2)+.7*exp(-0.5*((k-1250)/75).^2));
d1=adaptive_req.analysis.decompose_radial_spectrum(one);
d2=adaptive_req.analysis.decompose_radial_spectrum(two);
assert(d1.valid&&~d1.two_component_detected);
assert(d2.two_component_detected&&d2.low_k<d2.high_k);
X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
SPEC=load_spectral_cache(spectral_dir,key,X);
O=build_operational_inputs(X,SPEC,CFG);
F=load_field_cache(field_dir,X(1,:));
keep=false(height(X),1); low=find(O.confidence<CFG.ConfidenceThreshold,6,'first');
if isempty(low), low=1:min(6,height(X)); end
keep(low)=true; O.confidence(~keep)=1;
[pred,D]=evaluate_condition(O,SPEC,F,CFG);
assert_operational_policy(pred,D,O,CFG);
assert(all(isfinite(pred.(CFG.SelectedStrategy))));
fprintf('  Synthetic one/two-component decomposition passed.\n');
fprintf('  Cached-condition dimensions and masked fallback passed.\n');
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
B=B(:,A.Properties.VariableNames); T=[A;B];
end

function text=pretty(name)
text=replace(string(name),["hybrid_baseline","local_baseline","confidence_switch_c080", ...
    "two_component_local_c080","two_component_switch_c080","edge_masked_local_c080", ...
    "edge_masked_switch_c080","two_component_then_edge_masked_c080", ...
    "two_component_then_edge_masked_switch_c080"], ...
    ["Hybrid","Local","Switch","Two-component Local","Two-component Switch", ...
    "Edge-masked Local","Edge-masked Switch","Mixture + edge-mask Local", ...
    "Mixture + edge-mask Switch"]);
end

function style_axes(ax)
grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9;
end

function export_fig(fig,file)
exportgraphics(fig,file,'Resolution',200); close(fig);
end

function value=sanitize(value)
value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_');
end
