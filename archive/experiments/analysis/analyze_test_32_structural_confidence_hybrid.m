%% analyze_test_32_structural_confidence_hybrid.m
% Test 32: structural-confidence hybrid SWS maps.
%
% This analysis combines the less rigid Local/SWS-nearest maps from Test 31
% with the strong Theory/Test30 region-level estimate. The goal is to keep the
% robust region levels in trustworthy interiors, but relax toward local or
% high-confidence interpolation near estimated edges or when the structural
% prior looks less reliable. No true SWS, material labels, patch purity, or
% interface distance are used to decide the hybrid weights; those variables are
% used only for evaluation.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST32_MODE          = quick | full
%   ADAPTIVE_REQ_TEST32_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST32_SAVE_ALL_MAPS = true | false

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.12, ...
    'defaultAxesLabelFontSizeMultiplier',1.04);

%% Configuration

CFG=struct(); mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST32_MODE'))));
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]),'ADAPTIVE_REQ_TEST32_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST32_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST32_SAVE_ALL_MAPS',true);
CFG.Version=3; CFG.ConfidenceThreshold=env_double('ADAPTIVE_REQ_TEST32_CONFIDENCE_THRESHOLD',NaN);
CFG.BoundaryBandMm=env_double('ADAPTIVE_REQ_TEST32_BOUNDARY_BAND_MM',2.0);
CFG.InteriorBandMm=env_double('ADAPTIVE_REQ_TEST32_INTERIOR_BAND_MM',4.0);
CFG.DisagreementScale=env_double('ADAPTIVE_REQ_TEST32_DISAGREEMENT_SCALE',0.35);
CFG.EdgeBlendWidthMm=env_double('ADAPTIVE_REQ_TEST32_EDGE_BLEND_WIDTH_MM',4.0);
CFG.StructureDisagreementPenalty=env_double('ADAPTIVE_REQ_TEST32_STRUCTURE_DISAGREEMENT_PENALTY',0.35);
CFG.MinimumConfidence=0.05;
CFG.Strategies=["local_baseline","sws_nearest_highconf","test30_theory_region_levels", ...
    "hybrid_lowconf_region_else_local","hybrid_lowconf_region_else_sws_nearest", ...
    "hybrid_relaxed_region_blend","hybrid_boundary_protected_region", ...
    "hybrid_confidence_region_blend","hybrid_region_interior_sws_edge", ...
    "hybrid_conservative_structure","hybrid_adaptive_soft_region", ...
    "hybrid_final_candidate"];
CFG.MainStrategies=CFG.Strategies;

SOURCE=resolve_test31_source(root_dir,CFG);
OUT=make_output_dirs(root_dir,CFG); write_config_json(CFG,fullfile(OUT.root_dir,'test32_configuration.json'));

fprintf('\nTest 32: structural-confidence hybrid maps\n');
fprintf('Mode: %s | source: %s\n',CFG.RunMode,SOURCE);
S=load(SOURCE,'T_patch','CFG'); T=S.T_patch;
if ~isnan(CFG.ConfidenceThreshold)
    T.low_confidence=T.confidence<CFG.ConfidenceThreshold;
    T.high_confidence=~T.low_confidence;
end
fprintf('Loaded %d patches. Low-confidence fraction %.1f%%.\n',height(T),100*mean(T.low_confidence));

if CFG.ValidateOnly
    key=unique(T.condition_key,'stable'); validate_test32(T(T.condition_key==key(1),:),CFG);
    fprintf('Test 32 validation-only checks passed. No analysis was run.\n');
    return;
end

%% Evaluate condition maps

keys=unique(T.condition_key,'stable'); parts=cell(numel(keys),1); diag_parts=cell(numel(keys),1);
for ki=1:numel(keys)
    timer=tic; key=keys(ki); checkpoint=fullfile(OUT.condition_dir,"test32__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        C=load(checkpoint,'RESULT');
        if isfield(C,'RESULT')&&checkpoint_is_compatible(C.RESULT,CFG)
            parts{ki}=C.RESULT.T_patch; diag_parts{ki}=C.RESULT.T_diag;
            fprintf('[%d/%d] Reused %s.\n',ki,numel(keys),key); continue;
        end
    end
    X=sortrows(T(T.condition_key==key,:),{'map_iz','map_ix'});
    [R,D]=evaluate_condition(X,CFG);
    RESULT=struct('version',CFG.Version,'cache_signature',cache_signature(CFG),'T_patch',R,'T_diag',D); save(checkpoint,'RESULT','-v7.3');
    parts{ki}=R; diag_parts{ki}=D;
    fprintf('[%d/%d] %s: mean w_struct %.2f, edge %.1f%% in %.2f s.\n', ...
        ki,numel(keys),key,mean(R.hybrid_structural_weight,'omitnan'),100*mean(R.estimated_boundary_band),toc(timer));
end
T_patch=vertcat(parts{:}); T_diag=vertcat(diag_parts{:}); clear parts diag_parts T;

%% Summaries

T_overall=summarize_predictions(T_patch,"strategy_name",CFG);
T_by_M=summarize_predictions(T_patch,["strategy_name","M"],CFG);
T_by_geometry=summarize_predictions(T_patch,["strategy_name","geometry"],CFG);
T_by_regime=summarize_predictions(T_patch,["strategy_name","field_regime"],CFG);
T_by_frequency=summarize_predictions(T_patch,["strategy_name","f0"],CFG);
T_by_purity=summarize_predictions(T_patch,["strategy_name","geometry","purity_bin"],CFG);
T_by_distance=summarize_predictions(T_patch,["strategy_name","geometry","distance_bin"],CFG);
T_roi=roi_summary(T_patch,CFG);
T_best=best_strategy_candidates(T_overall,T_by_purity,T_roi,CFG);

writetable(T_patch,fullfile(OUT.table_dir,'test32_patch_level_results.csv'));
writetable(T_diag,fullfile(OUT.table_dir,'test32_condition_diagnostics.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test32_strategy_summary_overall.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'test32_strategy_summary_by_M.csv'));
writetable(T_by_geometry,fullfile(OUT.table_dir,'test32_strategy_summary_by_geometry.csv'));
writetable(T_by_regime,fullfile(OUT.table_dir,'test32_strategy_summary_by_regime.csv'));
writetable(T_by_frequency,fullfile(OUT.table_dir,'test32_strategy_summary_by_frequency.csv'));
writetable(T_by_purity,fullfile(OUT.table_dir,'test32_strategy_summary_by_purity_bin.csv'));
writetable(T_by_distance,fullfile(OUT.table_dir,'test32_strategy_summary_by_distance_bin.csv'));
writetable(T_roi,fullfile(OUT.table_dir,'test32_roi_summary.csv'));
writetable(T_best,fullfile(OUT.table_dir,'test32_best_strategy_candidates.csv'));
save(fullfile(OUT.data_dir,'test32_compact_results.mat'),'T_patch','T_overall','T_by_purity','T_by_distance','T_roi','T_diag','CFG','-v7.3');

plot_strategy_summary(T_patch,CFG,OUT);
plot_error_vs_distance(T_patch,CFG,OUT);
plot_weight_diagnostics(T_patch,CFG,OUT);
plot_representative_maps(T_patch,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T_patch,CFG,OUT); end
print_interpretation(T_overall,T_by_distance,T_roi,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 32 complete.\n',OUT.table_dir,OUT.figure_dir);

%% Setup

function SOURCE=resolve_test31_source(root,CFG)
base=fullfile(root,'outputs','test_31_simple_confidence_interpolation');
if CFG.QuickMode
    candidates=[string(fullfile(base,'quick','data','test31_compact_results.mat')), ...
        string(fullfile(base,'data','test31_compact_results.mat'))];
else
    candidates=string(fullfile(base,'data','test31_compact_results.mat'));
end
SOURCE="";
for c=candidates
    if exist(c,'file')==2, SOURCE=c; break; end
end
assert(SOURCE~="",'Missing Test 31 compact results. Run Test 31 first.');
end

function OUT=make_output_dirs(root,CFG)
OUT.root_dir=fullfile(root,'outputs','test_32_structural_confidence_hybrid');
if CFG.QuickMode, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables'); OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data'); OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
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
sig.ConfidenceThreshold=CFG.ConfidenceThreshold;
sig.BoundaryBandMm=CFG.BoundaryBandMm;
sig.InteriorBandMm=CFG.InteriorBandMm;
sig.DisagreementScale=CFG.DisagreementScale;
sig.EdgeBlendWidthMm=CFG.EdgeBlendWidthMm;
sig.StructureDisagreementPenalty=CFG.StructureDisagreementPenalty;
sig.Strategies=CFG.Strategies;
end

function write_config_json(CFG,file)
C=CFG; for f=string(fieldnames(C))', if isstring(C.(f)), C.(f)=cellstr(C.(f)); end, end
fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char'); clear cleanup;
end

%% Hybrid construction

function [R,D]=evaluate_condition(X,CFG)
X=sortrows(X,{'map_iz','map_ix'});
local=double(X.sws_local_baseline); nearest=double(X.sws_sws_nearest_highconf);
region=double(X.sws_test30_theory_region_levels); theory=double(X.sws_theory_baseline);
confidence=double(X.confidence); low=logical(X.low_confidence);
structure=choose_structure(X);
[boundary_mm,edge_band,interior_weight]=structure_distance(X,structure,CFG);
local_region=double(X.local_structure); theory_region=double(X.theory_structure);
structure_agree=double(local_region==theory_region | X.geometry_type=="homogeneous");
structure_agree_weight=CFG.StructureDisagreementPenalty+(1-CFG.StructureDisagreementPenalty)*structure_agree;
region_local_disagreement=abs(region-local)./max(abs(region),0.25);
disagreement_weight=exp(-region_local_disagreement/CFG.DisagreementScale);
low_weight=double(low);
if isfinite(CFG.ConfidenceThreshold)
    continuous_low_weight=min(max((CFG.ConfidenceThreshold-confidence)/max(CFG.ConfidenceThreshold,eps),0),1);
else
    continuous_low_weight=min(max((0.9-confidence)/0.9,0),1);
end
structural_weight=continuous_low_weight.*interior_weight.*structure_agree_weight.*max(disagreement_weight,0.35);
structural_weight(X.geometry_type=="homogeneous")=0;
structural_weight=min(max(structural_weight,0),1);
edge_weight=1-interior_weight;

pred=struct();
pred.local_baseline=local;
pred.sws_nearest_highconf=nearest;
pred.test30_theory_region_levels=region;
pred.hybrid_lowconf_region_else_local=local;
pred.hybrid_lowconf_region_else_local(low & X.geometry_type~="homogeneous")=region(low & X.geometry_type~="homogeneous");
pred.hybrid_lowconf_region_else_sws_nearest=nearest;
pred.hybrid_lowconf_region_else_sws_nearest(low & X.geometry_type~="homogeneous")=region(low & X.geometry_type~="homogeneous");
relaxed_weight=max(0.75*double(low),0.55*interior_weight.*structure_agree_weight);
relaxed_weight(X.geometry_type=="homogeneous")=0; relaxed_weight=min(max(relaxed_weight,0),1);
pred.hybrid_relaxed_region_blend=(1-relaxed_weight).*nearest + relaxed_weight.*region;
protected_weight=double(~edge_band).*structure_agree_weight;
protected_weight(X.geometry_type=="homogeneous")=0; protected_weight=min(max(protected_weight,0),1);
pred.hybrid_boundary_protected_region=(1-protected_weight).*nearest + protected_weight.*region;
pred.hybrid_confidence_region_blend=(1-structural_weight).*nearest + structural_weight.*region;
pred.hybrid_region_interior_sws_edge=(interior_weight.*region + edge_weight.*nearest);
pred.hybrid_region_interior_sws_edge(confidence>=max(.80,nanmedian(confidence)))=local(confidence>=max(.80,nanmedian(confidence)));
pred.hybrid_conservative_structure=nearest;
mask=(low & interior_weight>.5 & (structure_agree>0 | confidence<.35));
pred.hybrid_conservative_structure(mask)=0.65*region(mask)+0.35*nearest(mask);
pred.hybrid_adaptive_soft_region=(1-structural_weight).*local + structural_weight.*region;
edge_relax=min(max(boundary_mm/CFG.EdgeBlendWidthMm,0),1);
final_weight=max(structural_weight,0.75*edge_relax.*double(low).*structure_agree_weight);
final_weight(X.geometry_type=="homogeneous")=0;
pred.hybrid_final_candidate=(1-final_weight).*nearest + final_weight.*region;

R=attach_outputs(X,pred,CFG);
R.estimated_boundary_distance_mm=boundary_mm; R.estimated_boundary_band=edge_band;
R.hybrid_structural_weight=structural_weight; R.hybrid_final_weight=final_weight;
R.structure_agree=logical(structure_agree); R.region_local_disagreement=region_local_disagreement;
D=condition_diag(R,CFG);
end

function threshold=ConfidenceThreshold(CFG) %#ok<DEFNU>
threshold=CFG.ConfidenceThreshold;
end

function structure=choose_structure(X)
if ismember('test30_structure',X.Properties.VariableNames) && any(isfinite(X.test30_structure))
    structure=double(X.test30_structure);
else
    structure=double(X.theory_structure);
end
bad=~isfinite(structure); structure(bad)=1;
end

function [dist_mm,edge_band,interior_weight]=structure_distance(meta,labels,CFG)
[L,iz,ix]=map_from_rows(meta,labels); boundary=false(size(L));
boundary(:,1:end-1)=boundary(:,1:end-1)|(L(:,1:end-1)~=L(:,2:end));
boundary(:,2:end)=boundary(:,2:end)|(L(:,1:end-1)~=L(:,2:end));
boundary(1:end-1,:)=boundary(1:end-1,:)|(L(1:end-1,:)~=L(2:end,:));
boundary(2:end,:)=boundary(2:end,:)|(L(1:end-1,:)~=L(2:end,:));
dx_mm=1e3*mean([median(diff(unique(meta.x))) median(diff(unique(meta.z)))],'omitnan');
if any(boundary(:)), D=bwdist(boundary)*dx_mm; else, D=inf(size(L)); end
dist_mm=D(sub2ind(size(D),iz,ix)); edge_band=dist_mm<=CFG.BoundaryBandMm;
interior_weight=min(max((dist_mm-CFG.BoundaryBandMm)/(CFG.InteriorBandMm-CFG.BoundaryBandMm),0),1);
interior_weight(~isfinite(interior_weight))=1;
end

function R=attach_outputs(X,pred,CFG)
keep={'condition_key','geometry','geometry_type','field_regime','f0','M','dx','dz', ...
    'map_iz','map_ix','x','z','true_SWS','confidence','low_confidence','patch_purity', ...
    'material_side','purity_bin','distance_bin','distance_to_interface_mm'};
R=X(:,keep);
for s=CFG.Strategies
    y=pred.(s); R.("sws_"+s)=y;
    R.("signed_error_"+s)=100*(y-R.true_SWS)./R.true_SWS;
    R.("abs_error_"+s)=abs(R.("signed_error_"+s));
    R.("delta_vs_local_"+s)=y-pred.local_baseline;
end
end

function D=condition_diag(R,~)
D=R(1,{'condition_key','geometry','geometry_type','field_regime','f0','M','dx'});
D.N=height(R); D.low_confidence_pct=100*mean(R.low_confidence);
D.estimated_boundary_pct=100*mean(R.estimated_boundary_band);
D.mean_structural_weight=mean(R.hybrid_structural_weight,'omitnan');
D.mean_final_weight=mean(R.hybrid_final_weight,'omitnan');
D.structure_agreement_pct=100*mean(R.structure_agree,'omitnan');
D.mean_region_local_disagreement=mean(R.region_local_disagreement,'omitnan');
end

%% Summaries

function S=summarize_predictions(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for strategy=CFG.Strategies
    error=T.("signed_error_"+strategy);
    if isempty(basegroups), G=ones(height(T),1); else, G=findgroups(T(:,cellstr(basegroups))); end
    for id=unique(G,'stable')'
        idx=G==id; row=table(); for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy; row.N=sum(idx); row.MAPE=mean(abs(error(idx)),'omitnan');
        row.mean_signed_error=mean(error(idx),'omitnan'); row.median_signed_error=median(error(idx),'omitnan');
        row.high_error_10_pct=100*mean(abs(error(idx))>10,'omitnan');
        row.high_error_20_pct=100*mean(abs(error(idx))>20,'omitnan');
        row.underestimate_pct=100*mean(error(idx)<0,'omitnan'); row.corrected_pct=100*mean(abs(T.("delta_vs_local_"+strategy)(idx))>1e-9,'omitnan');
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            j=idx&T.purity_bin==cls; token=char(cls); row.("N_"+token)=sum(j);
            row.("MAPE_"+token)=mean(abs(error(j)),'omitnan'); row.("high20_"+token)=100*mean(abs(error(j))>20,'omitnan');
        end
        for side=["soft","hard"]
            j=idx&T.material_side==side; token=char(side); row.("MAPE_"+token)=mean(abs(error(j)),'omitnan');
            row.("mean_SWS_"+token)=mean(T.("sws_"+strategy)(j),'omitnan');
        end
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function T=roi_summary(P,CFG)
T=table();
for strategy=CFG.Strategies
    for geo=unique(P.geometry)'
        if startsWith(geo,"homogeneous"), continue; end
        for roi=["soft_core","hard_core","interface_0_2mm"]
            switch roi
                case "soft_core", idx=P.geometry==geo&P.material_side=="soft"&P.distance_to_interface_mm>8;
                case "hard_core", idx=P.geometry==geo&P.material_side=="hard"&P.distance_to_interface_mm>4;
                otherwise, idx=P.geometry==geo&P.distance_to_interface_mm<=2;
            end
            if ~any(idx), continue; end
            err=P.("signed_error_"+strategy);
            row=table(strategy,geo,roi,sum(idx),mean(P.("sws_"+strategy)(idx),'omitnan'), ...
                std(P.("sws_"+strategy)(idx),'omitnan'),mean(abs(err(idx)),'omitnan'),100*mean(abs(err(idx))>20,'omitnan'), ...
                'VariableNames',{'strategy_name','geometry','roi','N','mean_SWS','std_SWS','MAPE','high20_pct'});
            T=concat_tables(T,row);
        end
    end
end
end

function B=best_strategy_candidates(O,Purity,ROI,CFG)
B=table(); [~,i]=min(O.MAPE); B=concat_tables(B,table("global_MAPE",O.strategy_name(i),O.MAPE(i),'VariableNames',{'criterion','strategy_name','value'}));
[~,i]=min(O.high_error_20_pct); B=concat_tables(B,table("global_high20",O.strategy_name(i),O.high_error_20_pct(i),'VariableNames',{'criterion','strategy_name','value'}));
P=Purity(Purity.purity_bin=="strongly_mixed",:); if ~isempty(P), [~,i]=min(P.MAPE); B=concat_tables(B,table("strongly_mixed_MAPE",P.strategy_name(i),P.MAPE(i),'VariableNames',{'criterion','strategy_name','value'})); end
R=ROI(ROI.roi=="interface_0_2mm",:); if ~isempty(R), G=findgroups(R.strategy_name); m=splitapply(@mean,R.MAPE,G); names=splitapply(@(x)x(1),R.strategy_name,G); [~,i]=min(m); B=concat_tables(B,table("interface_ROI_MAPE",names(i),m(i),'VariableNames',{'criterion','strategy_name','value'})); end
end

%% Figures

function plot_strategy_summary(T,CFG,OUT)
strategies=CFG.MainStrategies; M=nan(numel(strategies),3); H=M;
for si=1:numel(strategies)
    e=T.("abs_error_"+strategies(si)); masks={T.purity_bin=="homogeneous",ismember(T.purity_bin,["pure","near_pure"]),ismember(T.purity_bin,["moderately_mixed","strongly_mixed"])};
    for j=1:3, M(si,j)=mean(e(masks{j}),'omitnan'); H(si,j)=100*mean(e(masks{j})>20,'omitnan'); end
end
fig=figure('Color','w','Position',[50 50 1350 720]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,M); xlabel(ax,'MAPE (%)'); title(ax,'Hybrid MAPE'); style_axes(ax);
ax=nexttile(tl); barh(ax,H); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style_axes(ax);
for ax=findall(fig,'Type','axes')', set(ax,'YTick',1:numel(strategies),'YTickLabel',pretty(strategies),'TickLabelInterpreter','none'); end
legend(["Homogeneous","Pure/near-pure","Mixed"],'Location','northoutside','Orientation','horizontal','Box','off'); export_fig(fig,OUT,'test32_strategy_summary.png');
end

function plot_error_vs_distance(T,CFG,OUT)
bins=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"]; labels=["0-1","1-2","2-4","4-8",">8"];
fig=figure('Color','w','Position',[80 80 1300 520]); tl=tiledlayout(1,2,'TileSpacing','compact');
for geo=["bilayer_2_3","circular_inclusion_2_3"]
    ax=nexttile(tl); hold(ax,'on');
    for strategy=CFG.MainStrategies
        y=nan(size(bins)); for bi=1:numel(bins), idx=T.geometry==geo&T.distance_bin==bins(bi); y(bi)=mean(T.("abs_error_"+strategy)(idx),'omitnan'); end
        plot(ax,1:numel(bins),y,'-o','LineWidth',1.2,'DisplayName',pretty(strategy));
    end
    set(ax,'XTick',1:numel(bins),'XTickLabel',labels); xlabel(ax,'True distance to interface (mm)'); ylabel(ax,'MAPE (%)'); title(ax,replace(geo,"_"," ")); style_axes(ax);
end
legend('Location','southoutside','Orientation','horizontal','NumColumns',3,'Box','off'); export_fig(fig,OUT,'test32_error_vs_distance.png');
end

function plot_weight_diagnostics(T,~,OUT)
fig=figure('Color','w','Position',[80 80 1100 460]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); scatter(ax,T.confidence,T.hybrid_final_weight,6,T.region_local_disagreement,'filled','MarkerFaceAlpha',.12); xlabel(ax,'Confidence'); ylabel(ax,'Final structural weight'); colorbar(ax); style_axes(ax);
ax=nexttile(tl); scatter(ax,T.estimated_boundary_distance_mm,T.hybrid_final_weight,6,T.confidence,'filled','MarkerFaceAlpha',.12); xlabel(ax,'Estimated boundary distance (mm)'); ylabel(ax,'Final structural weight'); colorbar(ax); style_axes(ax);
export_fig(fig,OUT,'test32_weight_diagnostics.png');
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
plot_case_maps(T(T.condition_key==keys(1),:),CFG,fullfile(OUT.figure_dir,'test32_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable'); root=fullfile(OUT.figure_dir,'maps_by_condition_dx200um');
for i=1:numel(keys)
    C=T(T.condition_key==keys(i),:); folder=fullfile(root,sanitize(C.geometry(1)),sanitize(C.field_regime(1))); if exist(folder,'dir')~=7, mkdir(folder); end
    plot_case_maps(C,CFG,fullfile(folder,"test32__"+sanitize(keys(i))+".png"));
end
end

function plot_case_maps(C,~,file)
C=sortrows(C,{'map_iz','map_ix'});
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_local_baseline),map_from_rows(C,C.sws_sws_nearest_highconf), ...
    map_from_rows(C,C.sws_test30_theory_region_levels),map_from_rows(C,C.sws_hybrid_confidence_region_blend), ...
    map_from_rows(C,C.sws_hybrid_lowconf_region_else_sws_nearest),map_from_rows(C,C.sws_hybrid_relaxed_region_blend), ...
    map_from_rows(C,C.sws_hybrid_region_interior_sws_edge),map_from_rows(C,C.sws_hybrid_conservative_structure), ...
    map_from_rows(C,C.sws_hybrid_final_candidate),map_from_rows(C,C.abs_error_local_baseline), ...
    map_from_rows(C,C.abs_error_sws_nearest_highconf), ...
    map_from_rows(C,C.abs_error_test30_theory_region_levels),map_from_rows(C,C.abs_error_hybrid_final_candidate), ...
    map_from_rows(C,C.confidence),map_from_rows(C,C.low_confidence),map_from_rows(C,C.estimated_boundary_distance_mm), ...
    map_from_rows(C,C.hybrid_final_weight),map_from_rows(C,C.region_local_disagreement)};
titles=["True","Local","SWS nearest","Test30 region","Conf region blend","Lowconf region / SWS", ...
    "Relaxed region blend","Region interior / SWS edge","Conservative structure","Final hybrid","Local error","SWS nearest error", ...
    "Test30 error","Final error","Confidence","Low confidence","Estimated boundary dist", ...
    "Final structural weight","Region-local disagreement"];
cols=5; rows=ceil(numel(maps)/cols);
fig=figure('Color','w','Visible','off','Position',[20 20 1700 1050]); tl=tiledlayout(rows,cols,'TileSpacing','compact');
for i=1:numel(maps), ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',8); colorbar(ax); end
title(tl,C.condition_key(1),'Interpreter','none'); exportgraphics(fig,file,'Resolution',170); close(fig);
end

function print_interpretation(O,D,R,~)
[~,i]=min(O.MAPE); fprintf('\nInterpretation:\n  Best global: %s (MAPE %.2f%%, >20 %.2f%%).\n',O.strategy_name(i),O.MAPE(i),O.high_error_20_pct(i));
near=D(D.distance_bin=="0_1mm",:); if ~isempty(near), G=groupsummary(near,'strategy_name','mean','MAPE'); [~,j]=min(G.mean_MAPE); fprintf('  Best 0-1 mm interface bin: %s (MAPE %.2f%%).\n',G.strategy_name(j),G.mean_MAPE(j)); end
iface=R(R.roi=="interface_0_2mm",:); if ~isempty(iface), G=groupsummary(iface,'strategy_name','mean','MAPE'); [~,j]=min(G.mean_MAPE); fprintf('  Best interface ROI: %s (MAPE %.2f%%).\n',G.strategy_name(j),G.mean_MAPE(j)); end
for s=["local_baseline","sws_nearest_highconf","test30_theory_region_levels","hybrid_final_candidate"]
    row=O(O.strategy_name==s,:); if ~isempty(row), fprintf('  %s: global %.2f%%, strong mixed %.2f%%, homogeneous %.2f%%.\n',pretty(s),row.MAPE,row.MAPE_strongly_mixed,row.MAPE_homogeneous); end
end
end

function validate_test32(C,CFG)
[R,D]=evaluate_condition(C,CFG); assert(height(R)==height(C)&&height(D)==1); assert(all(isfinite(R.sws_hybrid_final_candidate))); assert(all(R.hybrid_final_weight>=0&R.hybrid_final_weight<=1));
fprintf('  Hybrid maps, structural weights, and diagnostics passed for %s.\n',C.condition_key(1));
end

%% Generic helpers

function tf=env_true(name,default)
value=strtrim(getenv(name)); if isempty(value), tf=default; else, tf=any(strcmpi(value,{'true','1','yes','on'})); end
end
function value=env_double(name,default)
raw=strtrim(getenv(name)); if isempty(raw), value=default; else, value=str2double(raw); assert(isfinite(value)); end
end
function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux); A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end
function T=concat_tables(A,B)
if isempty(A), T=B; return; end
missing=setdiff(A.Properties.VariableNames,B.Properties.VariableNames); for m=missing, B.(m)=missing_value_like(A.(m),height(B)); end
extra=setdiff(B.Properties.VariableNames,A.Properties.VariableNames); for e=extra, A.(e)=missing_value_like(B.(e),height(A)); end
B=B(:,A.Properties.VariableNames); T=[A;B];
end
function value=missing_value_like(example,n)
if isstring(example), value=strings(n,1); elseif iscellstr(example), value=repmat({''},n,1); elseif islogical(example), value=false(n,1); else, value=nan(n,1); end
end
function label=pretty(name)
label=replace(string(name),["_","local baseline","sws nearest highconf","test30 theory region levels", ...
    "hybrid confidence region blend","hybrid region interior sws edge","hybrid conservative structure", ...
    "hybrid adaptive soft region","hybrid final candidate"], ...
    [" ","Local","SWS nearest","Test30 region","Confidence-region blend","Region interior / SWS edge", ...
    "Conservative structure","Adaptive Local-region","Final hybrid"]);
end
function value=sanitize(value), value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_'); end
function style_axes(ax), grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9; end
function export_fig(fig,OUT,name), exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',210); close(fig); end
