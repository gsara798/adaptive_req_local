%% analyze_test_34_kwave_mixedness_transfer.m
% Test 34: transfer Test 33 mixedness-aware corrections to k-Wave maps.
%
% This script does not retrain q models, confidence detectors, or Test 33
% post-processing models. It loads the cached Test 12 k-Wave REQ feature
% tables, applies frozen Test 19 q models, frozen Test 21 confidence detector
% ML_bagged_trees, and frozen Test 33 mixedness/log-k/q-selector models.
%
% k-Wave truth (2 m/s background, 3 m/s circular inclusion) is used only for
% evaluation, ROI summaries, and plots.

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9,'defaultLegendFontSize',8);

%% Configuration

CFG=struct();
CFG.Mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST34_MODE'))));
if CFG.Mode=="", CFG.Mode="full"; end
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST34_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST34_SAVE_ALL_MAPS',true);
CFG.PrimaryDetector="ML_bagged_trees";
CFG.FrequencyHz=500;
CFG.Geometry="kwave_circular_inclusion_2_3";
CFG.GeometryType="inclusion";
CFG.ConfidenceThreshold=0.8;
CFG.BoundaryBandMm=1.5;
CFG.InteriorBandMm=4.0;
CFG.StructureDisagreementPenalty=0.20;
CFG.DisagreementScale=0.25;
CFG.BaseStrategy="hybrid_lowconf_region_else_sws_nearest";
CFG.Strategies=["local_baseline","hybrid_baseline","theory_baseline", ...
    "sws_nearest_highconf","test30_theory_region_levels", ...
    "hybrid_lowconf_region_else_sws_nearest","hybrid_boundary_protected_region", ...
    "mixedness_logk_corrected","mixedness_q_candidate_selector"];

OUT=make_output_dirs(root_dir);
fprintf('\nTest 34: k-Wave transfer of mixedness-aware correction\n');
fprintf('Output: %s\n',OUT.root_dir);

KW_SRC=locate_kwave_source(root_dir);
T33_SRC=fullfile(root_dir,'outputs','test_33_mixedness_aware_q_correction','data','test33_models_and_predictions.mat');
assert(exist(KW_SRC,'file')==2,'Missing k-Wave source. Run Test 12 Level 06 first.');
assert(exist(T33_SRC,'file')==2,'Missing Test 33 full models. Run Test 33 full first.');

MODELS=load_all_models(root_dir,T33_SRC,CFG);
save(fullfile(OUT.model_dir,'test34_frozen_models_used.mat'),'-struct','MODELS','-v7.3');
write_config_json(CFG,fullfile(OUT.root_dir,'test34_configuration.json'));

if CFG.ValidateOnly
    validate_inputs(KW_SRC,MODELS,CFG);
    fprintf('Test 34 validation-only checks passed.\n');
    return;
end

S=load(KW_SRC,'CASE_OUT','KW');
CASE_OUT=S.CASE_OUT; KW=S.KW; clear S;

T_all=table(); T_roi=table(); T_side=table(); T_dist=table();
for ci=1:numel(CASE_OUT)
    fprintf('[%d/%d] %s M=%g\n',ci,numel(CASE_OUT),CASE_OUT(ci).case.case_name,CASE_OUT(ci).REQ_M);
    X=evaluate_kwave_condition(CASE_OUT(ci),KW,MODELS,CFG);
    T_all=concat_tables(T_all,X);
    T_roi=concat_tables(T_roi,summarize_roi(X,CFG));
    T_side=concat_tables(T_side,summarize_groups(X,["strategy_name","case_name","field_regime","M","material_side"],CFG));
    T_dist=concat_tables(T_dist,summarize_groups(X,["strategy_name","case_name","field_regime","M","distance_bin"],CFG));
    if CFG.SaveAllMaps
        plot_condition_maps(X,OUT);
    end
end

T_overall=summarize_groups(T_all,"strategy_name",CFG);
T_by_case=summarize_groups(T_all,["strategy_name","case_name","field_regime"],CFG);
T_by_M=summarize_groups(T_all,["strategy_name","M"],CFG);
T_by_region=summarize_groups(T_all,["strategy_name","region_label"],CFG);
T_by_confidence=summarize_confidence_bins(T_all,CFG);

writetable(T_all,fullfile(OUT.table_dir,'test34_kwave_patch_level_predictions.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test34_kwave_strategy_summary_overall.csv'));
writetable(T_by_case,fullfile(OUT.table_dir,'test34_kwave_strategy_summary_by_case.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'test34_kwave_strategy_summary_by_M.csv'));
writetable(T_roi,fullfile(OUT.table_dir,'test34_kwave_roi_summary.csv'));
writetable(T_side,fullfile(OUT.table_dir,'test34_kwave_soft_hard_summary.csv'));
writetable(T_dist,fullfile(OUT.table_dir,'test34_kwave_distance_summary.csv'));
writetable(T_by_region,fullfile(OUT.table_dir,'test34_kwave_region_summary.csv'));
writetable(T_by_confidence,fullfile(OUT.table_dir,'test34_kwave_confidence_reliability_summary.csv'));
save(fullfile(OUT.data_dir,'test34_kwave_transfer_results.mat'), ...
    'T_all','T_overall','T_by_case','T_by_M','T_roi','T_side','T_dist','T_by_region','T_by_confidence','CFG','KW','-v7.3');

plot_summary_bars(T_overall,OUT);
plot_roi_bars(T_roi,OUT);
plot_side_distance(T_side,T_dist,OUT);
plot_confidence_reliability(T_by_confidence,OUT);
write_summary_doc(root_dir,OUT);
print_summary(T_overall,T_roi,T_side);

fprintf('\nTables: %s\nFigures: %s\nModels copied: %s\nTest 34 complete.\n',OUT.table_dir,OUT.figure_dir,OUT.model_dir);

%% Main condition evaluation

function X=evaluate_kwave_condition(C,KW,MODELS,CFG)
F=C.T_feat;
F.dx=repmat(C.cfg.dx,height(F),1); F.dz=repmat(C.cfg.dz,height(F),1);
F.SIM_dx=F.dx; F.SIM_dz=F.dz;
F.cs_guess=F.REQ_cs_guess;
F.M=F.REQ_M;
F.f0=F.SIM_f0;
F.x=F.x_center_m; F.z=F.z_center_m;
F.geometry=repmat(CFG.Geometry,height(F),1);
F.geometry_type=repmat(CFG.GeometryType,height(F),1);
F.field_regime=repmat(map_kwave_regime(C.case.wave_model),height(F),1);
F.case_name=repmat(string(C.case.case_name),height(F),1);
F.condition_key=repmat(string(sprintf('kwave__%s__M%d',sanitize(C.case.case_name),C.REQ_M)),height(F),1);
F.true_SWS=F.cs_true;
F.material_side=repmat("soft",height(F),1); F.material_side(F.cs_true>2.5)="hard";
F.distance_to_interface_mm=abs(hypot(F.x-KW.inclusion_center_m(1),F.z-KW.inclusion_center_m(2))-KW.inclusion_radius_m)*1e3;
F.distance_bin=distance_bin(F.distance_to_interface_mm);
F.purity_bin=repmat("unknown_kWave",height(F),1);
F.region_label=F.material_side;

P=predict_q_models(F,MODELS.qmodels,C.cfg.f0);
F.sws_local_baseline=P.LocalOnly_T18.sws;
F.sws_hybrid_baseline=P.HybridLocalGlobal_T18_noUserRegime.sws;
q_theory=theory_q_kwave(F,C.case.wave_model,C.cfg.f0);
F.sws_theory_baseline=q_to_cs(q_theory,F.req_mapping,C.cfg.f0);
F.q_local_model=P.LocalOnly_T18.q;
F.q_hybrid_model=P.HybridLocalGlobal_T18_noUserRegime.q;
F.q_theory_model=q_theory;

F=add_detector_features(F);
F.confidence=apply_primary_detector(MODELS.detector,F);
F.low_confidence=F.confidence<CFG.ConfidenceThreshold;
F.sws_sws_nearest_highconf=nearest_highconf_by_region(F,F.sws_local_baseline);
F.sws_test30_theory_region_levels=region_levels_from_theory(F);
F=add_structural_hybrids(F,CFG);
F=apply_test33_models(F,MODELS.test33,CFG);
X=attach_strategy_errors(F,CFG);
end

function P=predict_q_models(F,qmodels,f0)
names=fieldnames(qmodels); P=struct();
for i=1:numel(names)
    Q=adaptive_req.analysis.predict_q_model_from_table(qmodels.(names{i}).model,F, ...
        'ModelType','bagged_trees','ModelName',string(names{i}));
    q=clamp_q(Q.q_pred);
    P.(names{i}).q=q;
    P.(names{i}).sws=q_to_cs(q,F.req_mapping,f0);
end
end

function F=add_detector_features(F)
F.q_pred_hybrid_T18=F.q_hybrid_model;
F.sws_pred_hybrid_T18=F.sws_hybrid_baseline;
[F.grad_q_hybrid,F.local_std_q_hybrid,F.local_range_q_hybrid]=map_local_stats(F,F.q_hybrid_model);
[F.grad_sws_hybrid,F.local_std_sws_hybrid,F.local_range_sws_hybrid]=map_local_stats(F,F.sws_hybrid_baseline);
F.abs_q_hybrid_minus_local=abs(F.q_hybrid_model-F.q_local_model);
F.abs_q_hybrid_minus_global=abs(F.q_hybrid_model-F.q_global_req);
F.abs_q_hybrid_minus_theory=abs(F.q_hybrid_model-F.q_theory_model);
F.abs_q_hybrid_minus_hybrid_guess=zeros(height(F),1);
F.abs_q_hybrid_T18_minus_old=zeros(height(F),1);
F.abs_sws_hybrid_minus_local=abs(F.sws_hybrid_baseline-F.sws_local_baseline);
F.abs_sws_hybrid_minus_global=abs(F.sws_hybrid_baseline-q_to_cs(clamp_q(F.q_global_req),F.req_mapping,F.SIM_f0(1)));
F.abs_sws_hybrid_minus_theory=abs(F.sws_hybrid_baseline-F.sws_theory_baseline);
F.abs_sws_hybrid_T18_minus_old=zeros(height(F),1);
end

function F=add_structural_hybrids(F,CFG)
structure=double(F.sws_test30_theory_region_levels>2.5);
[boundary_mm,edge_band,~]=structure_distance(F,structure,CFG);
F.estimated_boundary_distance_mm=boundary_mm;
F.estimated_boundary_band=edge_band;
F.structure_agree=true(height(F),1);
F.region_local_disagreement=abs(F.sws_test30_theory_region_levels-F.sws_local_baseline)./max(abs(F.sws_test30_theory_region_levels),0.25);
nearest=F.sws_sws_nearest_highconf; region=F.sws_test30_theory_region_levels;
F.sws_hybrid_lowconf_region_else_sws_nearest=nearest;
mask=F.low_confidence; F.sws_hybrid_lowconf_region_else_sws_nearest(mask)=region(mask);
protected_weight=double(~edge_band); F.hybrid_structural_weight=protected_weight;
F.hybrid_final_weight=protected_weight;
F.sws_hybrid_boundary_protected_region=(1-protected_weight).*nearest+protected_weight.*region;
end

function F=apply_test33_models(F,T33,CFG)
[Feat,~]=build_test33_features(F,CFG);
[p_mixed,p_strong,purity_pred]=predict_mixedness(T33.MIX,Feat);
F.predicted_mixed_probability=p_mixed;
F.predicted_strong_mixed_probability=p_strong;
F.predicted_patch_purity=min(max(purity_pred,0),1);
Fcorr=[Feat p_mixed p_strong F.predicted_patch_purity];
base=F.sws_hybrid_lowconf_region_else_sws_nearest;
base_k=2*pi*F.f0./base;
delta=predict(T33.CORR.model,Fcorr);
gate=max(p_mixed,0.35*p_strong);
F.mixedness_gate=gate; F.delta_logk_pred=delta;
logk=(1-gate).*log(base_k)+gate.*(log(base_k)+delta);
F.sws_mixedness_logk_corrected=clip_sws(2*pi*F.f0./exp(logk));
F.sws_mixedness_q_candidate_selector=F.sws_sws_nearest_highconf;
F.q_mixedness_q_candidate_selector=F.q_local_model;
F.q_candidate_selected=repmat("kwave_proxy_sws_nearest",height(F),1);
correction_mag=abs(F.sws_mixedness_logk_corrected-base)./max(base,0.25);
spread=estimate_spread(F,base,F.sws_mixedness_logk_corrected);
Freliab=[Fcorr correction_mag spread gate delta];
[~,score]=predict(T33.REL.model,Freliab);
F.posterior_high_error20_probability=positive_score(T33.REL.model,score);
F.posterior_reliability=1-F.posterior_high_error20_probability;
end

function X=attach_strategy_errors(F,CFG)
keep={'condition_key','case_name','field_regime','geometry','geometry_type','M','f0','dx','dz', ...
    'map_iz','map_ix','x','z','true_SWS','cs_true','roi_name','material_side','region_label', ...
    'distance_to_interface_mm','distance_bin','purity_bin','confidence','low_confidence', ...
    'posterior_reliability','predicted_mixed_probability','predicted_patch_purity'};
X=F(:,keep);
for s=CFG.Strategies
    sws=strategy_sws(F,s);
    X.("sws_"+s)=sws;
    err=100*(sws-F.true_SWS)./F.true_SWS;
    X.("signed_error_"+s)=err;
    X.("abs_error_"+s)=abs(err);
end
end

%% Loading

function MODELS=load_all_models(root,T33_SRC,CFG)
MODELS=struct();
MODELS.qmodels=load_q_models(root);
MODELS.detector=load_primary_detector(root,CFG.PrimaryDetector);
MODELS.test33=load(T33_SRC,'MIX','CORR','QSEL','REL','CFG','feature_names','f_corr_names','reliab_names');
end

function qmodels=load_q_models(root)
dir19=fullfile(root,'outputs','model_registry','test19_clean_field_regime');
requested=["LocalOnly_T18","HybridLocalGlobal_T18_noUserRegime"];
features=["CleanFieldRegime_noUser","CleanFieldRegime_noUser"];
qmodels=struct();
for i=1:numel(requested)
    [M,~,~]=adaptive_req.analysis.load_q_model_deployment(dir19, ...
        'ModelName',requested(i),'FeatureSet',features(i),'ModelType','bagged_trees');
    qmodels.(requested(i)).model=M;
end
end

function detector=load_primary_detector(root,name)
file=fullfile(root,'outputs','test_21_interface_confidence_detector','analysis','models','level21_detector_training_checkpoint.mat');
S=load(file,'DETECTORS'); D=S.DETECTORS; names=string({D.detector_name});
idx=find(names==name,1); assert(~isempty(idx),'Detector %s not found.',name);
detector=D(idx);
end

function src=locate_kwave_source(root)
src=fullfile(root,'outputs','test_12_cs_guess_window_sweep','test_12_cs_guess_window_sweep_2026-06-10_122903', ...
    'analysis','level_06_kwave_transfer','data','level12_level06_kwave_transfer.mat');
end

%% Summaries and plots

function T=summarize_groups(P,groups,CFG)
T=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for s=CFG.Strategies
    err=P.("signed_error_"+s);
    if isempty(basegroups), G=ones(height(P),1); else, G=findgroups(P(:,cellstr(basegroups))); end
    for id=unique(G,'stable')'
        idx=G==id; row=table(); for v=basegroups, row.(v)=P.(v)(find(idx,1)); end
        row.strategy_name=s; row.N=sum(idx);
        row.MAPE=mean(abs(err(idx)),'omitnan');
        row.mean_signed_error=mean(err(idx),'omitnan');
        row.high_error_20_pct=100*mean(abs(err(idx))>20,'omitnan');
        row.underestimate_pct=100*mean(err(idx)<0,'omitnan');
        T=concat_tables(T,row);
    end
end
T=movevars(T,'strategy_name','Before',1);
end

function T=summarize_roi(P,CFG)
T=summarize_groups(P(P.roi_name~="outside_roi",:),["strategy_name","case_name","field_regime","M","roi_name"],CFG);
end

function T=summarize_confidence_bins(P,CFG)
edges=0:0.1:1; T=table();
for s=CFG.Strategies
    err=P.("abs_error_"+s);
    for b=1:numel(edges)-1
        idx=P.confidence>=edges(b)&P.confidence<edges(b+1); if b==numel(edges)-1, idx=P.confidence>=edges(b)&P.confidence<=edges(b+1); end
        row=table(s,mean(edges(b:b+1)),sum(idx),mean(P.confidence(idx),'omitnan'),mean(P.posterior_reliability(idx),'omitnan'),mean(err(idx),'omitnan'),100*mean(err(idx)>20,'omitnan'), ...
            'VariableNames',{'strategy_name','confidence_bin_center','N','mean_confidence','mean_posterior_reliability','MAPE','high_error_20_pct'});
        T=concat_tables(T,row);
    end
end
end

function plot_summary_bars(T,OUT)
fig=figure('Color','w','Position',[80 80 1100 480]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,T.MAPE); set(ax,'YTick',1:height(T),'YTickLabel',pretty(T.strategy_name),'TickLabelInterpreter','none'); xlabel(ax,'MAPE (%)'); title(ax,'k-Wave overall MAPE'); style_axes(ax);
ax=nexttile(tl); barh(ax,T.high_error_20_pct); set(ax,'YTick',1:height(T),'YTickLabel',pretty(T.strategy_name),'TickLabelInterpreter','none'); xlabel(ax,'>20% error (%)'); title(ax,'k-Wave large-error rate'); style_axes(ax);
exportgraphics(fig,fullfile(OUT.figure_dir,'test34_kwave_strategy_summary.png'),'Resolution',220); close(fig);
end

function plot_roi_bars(T,OUT)
main=ismember(T.strategy_name,["local_baseline","theory_baseline","test30_theory_region_levels","hybrid_boundary_protected_region","mixedness_logk_corrected"]);
T=T(main,:);
fig=figure('Color','w','Position',[50 50 1450 650]); tl=tiledlayout(2,1,'TileSpacing','compact');
for metric=["MAPE","mean_signed_error"]
    ax=nexttile(tl); hold(ax,'on'); cats=strcat(string(T.case_name)," / M",string(T.M)," / ",string(T.roi_name));
    [G,labels]=findgroups(cats); strategies=unique(string(T.strategy_name),'stable'); width=.8/numel(strategies);
    for i=1:numel(strategies)
        idx=T.strategy_name==strategies(i); vals=splitapply(@(x)mean(x,'omitnan'),T.(metric)(idx),G(idx));
        bar(ax,(1:numel(labels))+(i-(numel(strategies)+1)/2)*width,vals,width,'DisplayName',pretty(strategies(i)));
    end
    set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',25,'TickLabelInterpreter','none');
    ylabel(ax,replace(metric,"_"," ")); legend(ax,'Location','bestoutside'); style_axes(ax);
end
exportgraphics(fig,fullfile(OUT.figure_dir,'test34_kwave_roi_error_bars.png'),'Resolution',220); close(fig);
end

function plot_side_distance(Tside,Tdist,OUT)
fig=figure('Color','w','Position',[60 60 1400 650]); tl=tiledlayout(1,2,'TileSpacing','compact');
main=["local_baseline","theory_baseline","test30_theory_region_levels","hybrid_boundary_protected_region","mixedness_logk_corrected"];
ax=nexttile(tl); hold(ax,'on'); T=Tside(ismember(Tside.strategy_name,main),:);
cats=strcat(string(T.material_side)," / M",string(T.M)); [G,labels]=findgroups(cats); width=.8/numel(main);
for i=1:numel(main), idx=T.strategy_name==main(i); vals=splitapply(@(x)mean(x,'omitnan'),T.MAPE(idx),G(idx)); bar(ax,(1:numel(labels))+(i-(numel(main)+1)/2)*width,vals,width,'DisplayName',pretty(main(i))); end
set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'TickLabelInterpreter','none'); ylabel(ax,'MAPE (%)'); title(ax,'Soft/hard summary'); legend(ax,'Location','bestoutside'); style_axes(ax);
ax=nexttile(tl); hold(ax,'on'); T=Tdist(ismember(Tdist.strategy_name,main)&Tdist.M==2,:);
order=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
for i=1:numel(main)
    vals=nan(size(order)); for b=1:numel(order), vals(b)=mean(T.MAPE(T.strategy_name==main(i)&T.distance_bin==order(b)),'omitnan'); end
    plot(ax,1:numel(order),vals,'-o','DisplayName',pretty(main(i)));
end
set(ax,'XTick',1:numel(order),'XTickLabel',replace(order,"_","-")); xlabel(ax,'Distance bin'); ylabel(ax,'MAPE (%)'); title(ax,'Distance to interface, M=2'); legend(ax,'Location','best'); style_axes(ax);
exportgraphics(fig,fullfile(OUT.figure_dir,'test34_kwave_soft_hard_distance.png'),'Resolution',220); close(fig);
end

function plot_confidence_reliability(T,OUT)
fig=figure('Color','w','Position',[80 80 800 430]); ax=axes(fig); hold(ax,'on');
for s=["mixedness_logk_corrected","hybrid_boundary_protected_region","local_baseline"]
    X=T(T.strategy_name==s,:); plot(ax,X.confidence_bin_center,X.MAPE,'-o','DisplayName',pretty(s));
end
xlabel(ax,'Confidence bin'); ylabel(ax,'MAPE (%)'); title(ax,'k-Wave error vs frozen confidence'); legend(ax,'Location','best'); style_axes(ax);
exportgraphics(fig,fullfile(OUT.figure_dir,'test34_kwave_error_vs_confidence.png'),'Resolution',220); close(fig);
end

function plot_condition_maps(X,OUT)
folder=fullfile(OUT.figure_dir,'maps_by_condition',sanitize(X.case_name(1))); if exist(folder,'dir')~=7, mkdir(folder); end
maps={map_from_rows(X,X.true_SWS),map_from_rows(X,X.sws_local_baseline),map_from_rows(X,X.sws_theory_baseline),map_from_rows(X,X.sws_test30_theory_region_levels), ...
    map_from_rows(X,X.sws_hybrid_boundary_protected_region),map_from_rows(X,X.sws_mixedness_logk_corrected),map_from_rows(X,X.abs_error_local_baseline), ...
    map_from_rows(X,X.abs_error_mixedness_logk_corrected),map_from_rows(X,X.confidence),map_from_rows(X,X.posterior_reliability),map_from_rows(X,X.predicted_mixed_probability),map_from_rows(X,X.distance_to_interface_mm)};
titles=["True SWS","Local","Theory","Region levels","Boundary hybrid","Mixedness log-k","Local error","Log-k error","Confidence","Posterior reliability","P(mixed)","Distance mm"];
fig=figure('Color','w','Visible','off','Position',[20 20 1500 850]); tl=tiledlayout(3,4,'TileSpacing','compact');
for i=1:numel(maps), ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',8); colorbar(ax); end
title(tl,sprintf('%s | M=%d',X.case_name(1),X.M(1)),'Interpreter','none');
exportgraphics(fig,fullfile(folder,sprintf('test34_kwave__%s__M%d.png',sanitize(X.case_name(1)),X.M(1))),'Resolution',180); close(fig);
end

%% Feature builders

function [F,names]=build_test33_features(T,~)
base=T.sws_hybrid_lowconf_region_else_sws_nearest; local=T.sws_local_baseline; nearest=T.sws_sws_nearest_highconf; region=T.sws_test30_theory_region_levels;
operational_spread=estimate_spread(T,base,base);
vars=[T.confidence,double(T.low_confidence),normalize_like(T.f0,500,129.1),normalize_like(T.M,2.5,.5), ...
    regime_code(T.field_regime),zeros(height(T),1),log(max(local,0.1)),log(max(nearest,0.1)),log(max(region,0.1)),log(max(base,0.1)), ...
    nearest-local,region-local,region-nearest,abs(nearest-local)./max(local,0.25),abs(region-local)./max(local,0.25), ...
    T.estimated_boundary_distance_mm,double(T.estimated_boundary_band),T.hybrid_structural_weight,T.hybrid_final_weight,double(T.structure_agree),T.region_local_disagreement,operational_spread];
names=["confidence","low_confidence","f0_norm","M_norm","field_regime_code","geometry_type_code","log_sws_local","log_sws_nearest","log_sws_region","log_sws_base","nearest_minus_local","region_minus_local","region_minus_nearest","rel_nearest_local","rel_region_local","estimated_boundary_distance_mm","estimated_boundary_band","hybrid_structural_weight","hybrid_final_weight","structure_agree","region_local_disagreement","ensemble_spread_operational"];
F=double(vars); F(~isfinite(F))=0;
end

function confidence=apply_primary_detector(det,F)
names=string(det.predictors); X=zeros(height(F),numel(names));
for j=1:numel(names)
    if ismember(names(j),string(F.Properties.VariableNames)), X(:,j)=double(F.(names(j))); else, X(:,j)=det.impute_median(j); end
end
for j=1:size(X,2), bad=~isfinite(X(:,j)); X(bad,j)=det.impute_median(j); end
X=(X-det.predictor_mu)./det.predictor_sigma;
[~,score]=predict(det.model,X); risk=positive_score(det.model,score); confidence=1-risk;
end

function [p_mixed,p_strong,purity_pred]=predict_mixedness(MIX,F)
[~,score]=predict(MIX.cls,F); p_mixed=positive_score(MIX.cls,score);
[~,score]=predict(MIX.strong,F); p_strong=positive_score(MIX.strong,score);
purity_pred=predict(MIX.purity,F);
end

%% Map operations

function y=nearest_highconf_by_region(T,values)
structure=double(T.sws_theory_baseline>2.5); high=T.confidence>=0.8; y=values;
for region=unique(structure)'
    idx=structure==region; donors=idx&high&isfinite(values);
    if ~any(donors), donors=idx&isfinite(values); end
    if any(donors)
        med=median(values(donors),'omitnan'); y(idx&~high)=med;
    end
end
end

function y=region_levels_from_theory(T)
structure=double(T.sws_theory_baseline>2.5); y=T.sws_theory_baseline;
for region=unique(structure)'
    idx=structure==region; donors=idx&T.confidence>=0.8&isfinite(T.sws_local_baseline);
    if ~any(donors), donors=idx&isfinite(T.sws_local_baseline); end
    if any(donors), y(idx)=median(T.sws_local_baseline(donors),'omitnan'); end
end
end

function [gradv,stdv,rangev]=map_local_stats(T,v)
[A,iz,ix]=map_from_rows(T,v); [gx,gz]=gradient(A); G=hypot(gx,gz);
gradv=G(sub2ind(size(A),iz,ix));
K=ones(3); valid=isfinite(A); A0=A; A0(~valid)=0; n=conv2(double(valid),K,'same');
mu=conv2(A0,K,'same')./max(n,1); mu2=conv2(A0.^2,K,'same')./max(n,1);
S=sqrt(max(mu2-mu.^2,0)); stdv=S(sub2ind(size(A),iz,ix));
mx=local_extreme(A,true); mn=local_extreme(A,false); R=mx-mn; rangev=R(sub2ind(size(A),iz,ix));
gradv(~isfinite(gradv))=0; stdv(~isfinite(stdv))=0; rangev(~isfinite(rangev))=0;
end

function E=local_extreme(A,ismax)
E=nan(size(A)); pad=nan(size(A)+2); pad(2:end-1,2:end-1)=A;
for i=1:size(A,1)
    for j=1:size(A,2)
        W=pad(i:i+2,j:j+2); if ismax, E(i,j)=max(W(:),[],'omitnan'); else, E(i,j)=min(W(:),[],'omitnan'); end
    end
end
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

%% Utilities

function sws=strategy_sws(T,s)
switch string(s)
    case "local_baseline", sws=T.sws_local_baseline;
    case "hybrid_baseline", sws=T.sws_hybrid_baseline;
    case "theory_baseline", sws=T.sws_theory_baseline;
    case "sws_nearest_highconf", sws=T.sws_sws_nearest_highconf;
    case "test30_theory_region_levels", sws=T.sws_test30_theory_region_levels;
    case "hybrid_lowconf_region_else_sws_nearest", sws=T.sws_hybrid_lowconf_region_else_sws_nearest;
    case "hybrid_boundary_protected_region", sws=T.sws_hybrid_boundary_protected_region;
    case "mixedness_logk_corrected", sws=T.sws_mixedness_logk_corrected;
    case "mixedness_q_candidate_selector", sws=T.sws_mixedness_q_candidate_selector;
    otherwise, error('Unknown strategy %s',s);
end
end

function cs=q_to_cs(q,mappings,f0)
cs=nan(size(q));
for i=1:numel(q)
    if isfinite(q(i)), cs(i)=adaptive_req.quantile.quantile_to_cs(mappings{i},q(i),f0); end
end
end

function q=clamp_q(q), q=min(max(q,0.001),0.999); end
function sws=clip_sws(sws), sws=min(max(sws,0.5),10); end
function q=theory_q_kwave(T,wave_model,f0)
switch string(wave_model)
    case "SingleWave", type="SingleWave";
    case "Diffuse2D", type="Diffuse2D";
    case "Diffuse3D", type="Diffuse3D";
    otherwise, type="Partial3D";
end
if type=="Partial3D"
    q=.5*(theory_one(T,"Diffuse2D",f0)+theory_one(T,"Diffuse3D",f0));
else
    q=theory_one(T,type,f0);
end
q=repmat(clamp_q(q),height(T),1);
end
function q=theory_one(T,type,f0)
o=adaptive_req.theory.q_theory_REQ_discrete_shearUZ(T.dx(1),T.dz(1), ...
    f0,T.REQ_cs_guess(1),'M',T.REQ_M(1),'Gamma',1, ...
    'PadFactor',1,'Nbins','auto','SmoothSigma',1,'TheoryMode','S2D', ...
    'FieldType',type,'Plot',false);
q=o.q_th;
end
function x=normalize_like(x,mu,sigma), x=(double(x)-mu)./max(sigma,eps); end
function c=regime_code(r), [~,~,k]=unique(string(r)); c=double(k)-mean(double(k)); end
function r=map_kwave_regime(w)
switch string(w)
    case "SingleWave", r="directional_2D";
    case "Diffuse2D", r="diffuse_2D";
    case "Diffuse3D", r="diffuse_3D";
    otherwise, r="partial_3D";
end
end

function b=distance_bin(d)
edges=[0 1 2 4 8 inf]; labs=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"]; b=strings(size(d));
for i=1:numel(labs), b(d>=edges(i)&d<edges(i+1))=labs(i); end
end

function spread=estimate_spread(T,base,corrected)
V=[T.sws_local_baseline T.sws_sws_nearest_highconf T.sws_test30_theory_region_levels T.sws_hybrid_lowconf_region_else_sws_nearest T.sws_hybrid_boundary_protected_region base corrected];
spread=std(V,0,2,'omitnan')./max(mean(abs(V),2,'omitnan'),0.25); spread(~isfinite(spread))=0;
end

function p=positive_score(model,score)
classes=string(model.ClassNames); pos=find(classes=="true"|classes=="1",1); if isempty(pos), pos=size(score,2); end
p=min(max(score(:,pos),0),1);
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
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
function tf=env_true(name,default), raw=strtrim(getenv(name)); if isempty(raw), tf=default; else, tf=any(strcmpi(raw,{'true','1','yes','on'})); end, end
function OUT=make_output_dirs(root)
OUT.root_dir=fullfile(root,'outputs','test_34_kwave_mixedness_transfer');
OUT.table_dir=fullfile(OUT.root_dir,'tables'); OUT.figure_dir=fullfile(OUT.root_dir,'figures'); OUT.data_dir=fullfile(OUT.root_dir,'data'); OUT.model_dir=fullfile(OUT.root_dir,'models');
for d=string(struct2cell(OUT))', if exist(d,'dir')~=7, mkdir(d); end, end
end
function write_config_json(CFG,file), fid=fopen(file,'w'); fwrite(fid,jsonencode(CFG,'PrettyPrint',true),'char'); fclose(fid); end
function s=sanitize(s), s=regexprep(char(string(s)),'[^A-Za-z0-9_-]','_'); end
function label=pretty(s), label=replace(string(s),"_"," "); end
function style_axes(ax), grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); end

function validate_inputs(src,MODELS,CFG)
S=load(src,'CASE_OUT'); assert(isfield(S,'CASE_OUT')&&~isempty(S.CASE_OUT));
assert(isfield(MODELS,'qmodels')&&isfield(MODELS.qmodels,'LocalOnly_T18'));
assert(isfield(MODELS,'detector')&&string(MODELS.detector.detector_name)==CFG.PrimaryDetector);
assert(isfield(MODELS,'test33')&&isfield(MODELS.test33,'MIX'));
end

function print_summary(T,Troi,Tside)
[~,i]=min(T.MAPE); fprintf('\nBest k-Wave global MAPE: %s %.2f%%.\n',T.strategy_name(i),T.MAPE(i));
disp(T(:,{'strategy_name','MAPE','mean_signed_error','high_error_20_pct'}));
fprintf('\nROI summary:\n'); disp(Troi(:,{'strategy_name','case_name','M','roi_name','MAPE','mean_signed_error','high_error_20_pct'}));
fprintf('\nSoft/hard summary:\n'); disp(Tside(:,{'strategy_name','case_name','M','material_side','MAPE','mean_signed_error','high_error_20_pct'}));
end

function write_summary_doc(root,OUT)
doc=fullfile(root,'docs','TEST34_KWAVE_MIXEDNESS_TRANSFER.md');
fid=fopen(doc,'w'); assert(fid>=0); c=onCleanup(@()fclose(fid));
fprintf(fid,'# Test 34: k-Wave Mixedness Transfer\n\n');
fprintf(fid,'Applies frozen Test 19 q models, frozen Test 21 confidence detector, and frozen Test 33 mixedness/log-k models to cached k-Wave inclusion simulations.\n\n');
fprintf(fid,'No model is retrained. k-Wave truth is used only for evaluation and plots.\n\n');
fprintf(fid,'Outputs: `%s`.\n\n',OUT.root_dir);
fprintf(fid,'Key tables include patch-level predictions, ROI summaries, soft/hard summaries, distance-to-interface summaries, and confidence/reliability summaries.\n');
end
