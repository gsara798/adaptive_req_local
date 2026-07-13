%% analyze_test_33_mixedness_aware_q_correction.m
% Test 33: mixedness-aware log-k correction and posterior reliability.
%
% This analysis trains three post-processing models using cached Test 32 maps:
%   1) an operational mixedness detector, trained against diagnostic
%      patch_purity labels;
%   2) a residual log-k corrector, conditioned on predicted mixedness;
%   3) a posterior reliability model for the corrected map.
%
% The correction and reliability features never include true SWS,
% patch_purity, true material side, or true distance-to-interface. Those are
% labels/evaluation variables only. The q model and confidence detectors remain
% frozen; this is a post-map corrector.
%
% Runtime controls:
%   ADAPTIVE_REQ_TEST33_MODE          = quick | full
%   ADAPTIVE_REQ_TEST33_VALIDATE_ONLY = true | false
%   ADAPTIVE_REQ_TEST33_SAVE_ALL_MAPS = true | false
%   ADAPTIVE_REQ_TEST33_BASE_STRATEGY = hybrid_lowconf_region_else_sws_nearest

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.12, ...
    'defaultAxesLabelFontSizeMultiplier',1.04);

%% Configuration
CFG=struct(); mode=lower(strtrim(string(getenv('ADAPTIVE_REQ_TEST33_MODE'))));
if mode=="", mode="quick"; end
assert(ismember(mode,["quick","full"]),'ADAPTIVE_REQ_TEST33_MODE must be quick or full.');
CFG.RunMode=mode; CFG.QuickMode=mode=="quick";
CFG.ValidateOnly=env_true('ADAPTIVE_REQ_TEST33_VALIDATE_ONLY',false);
CFG.SaveAllMaps=env_true('ADAPTIVE_REQ_TEST33_SAVE_ALL_MAPS',true);
CFG.Version=1; CFG.RandomSeed=33001;
CFG.BaseStrategy=string(getenv('ADAPTIVE_REQ_TEST33_BASE_STRATEGY'));
if CFG.BaseStrategy=="", CFG.BaseStrategy="hybrid_lowconf_region_else_sws_nearest"; end
CFG.MixedThreshold=0.95; CFG.StrongMixedThreshold=0.75;
CFG.TrainFraction=0.70; CFG.MaxTrainRows=env_double('ADAPTIVE_REQ_TEST33_MAX_TRAIN_ROWS',120000);
CFG.NumLearningCycles=env_double('ADAPTIVE_REQ_TEST33_NUM_TREES',120);
CFG.MinLeafSize=env_double('ADAPTIVE_REQ_TEST33_MIN_LEAF_SIZE',40);
CFG.Shrinkage=0.05; CFG.PhysicalRange=[0.5 10];
CFG.EvalStrategies=["local_baseline","sws_nearest_highconf","test30_theory_region_levels", ...
    "hybrid_lowconf_region_else_sws_nearest","hybrid_relaxed_region_blend", ...
    "hybrid_boundary_protected_region","mixedness_logk_corrected", ...
    "mixedness_q_candidate_selector"];

SOURCE=resolve_test32_source(root_dir,CFG);
OUT=make_output_dirs(root_dir,CFG); write_config_json(CFG,fullfile(OUT.root_dir,'test33_configuration.json'));

fprintf('\nTest 33: mixedness-aware log-k correction and posterior reliability\n');
fprintf('Mode: %s | source: %s\n',CFG.RunMode,SOURCE);
fprintf('Base strategy: %s\n',CFG.BaseStrategy);
S=load(SOURCE,'T_patch'); T=S.T_patch; clear S;
assert(ismember("sws_"+CFG.BaseStrategy,T.Properties.VariableNames),'Base strategy missing from Test 32 table.');
T=sortrows(T,{'condition_key','map_iz','map_ix'});
[T,q_available,q_meta]=attach_test31_q_candidates(T,root_dir,CFG);
fprintf('Loaded %d patches.\n',height(T));
if q_available
    fprintf('Attached Test 31 q candidates: %s.\n',strjoin(q_meta.strategy_names,", "));
else
    fprintf('Test 31 q candidates unavailable; q-selector strategy will fall back to LocalOnly.\n');
end

if CFG.ValidateOnly
    validate_test33(T(1:min(height(T),4096),:),CFG,q_available);
    fprintf('Test 33 validation-only checks passed. No analysis was run.\n');
    return;
end

rng(CFG.RandomSeed);
[train_mask,test_mask,split_table]=make_condition_split(T,CFG);
fprintf('Condition split: %d train rows, %d test rows.\n',sum(train_mask),sum(test_mask));

%% Feature matrix and targets

[F,feature_names]=build_features(T,CFG);
base_sws=T.("sws_"+CFG.BaseStrategy);
base_sws=clip_sws(base_sws,CFG);
true_k=2*pi*T.f0./T.true_SWS;
base_k=2*pi*T.f0./base_sws;
delta_logk_target=log(true_k)-log(base_k);
is_mixed=T.patch_purity<CFG.MixedThreshold;
is_strong_mixed=T.patch_purity<CFG.StrongMixedThreshold;
patch_purity_target=T.patch_purity;
high_error20_base=abs(100*(base_sws-T.true_SWS)./T.true_SWS)>20;

train_rows=find(train_mask & all(isfinite(F),2) & isfinite(delta_logk_target));
if numel(train_rows)>CFG.MaxTrainRows
    train_rows=train_rows(randperm(numel(train_rows),CFG.MaxTrainRows));
end
test_rows=find(test_mask & all(isfinite(F),2));
assert(~isempty(train_rows)&&~isempty(test_rows),'Empty train/test split.');

%% Train mixedness models

fprintf('Training mixedness detector on %d rows...\n',numel(train_rows));
tree_cls=templateTree('MinLeafSize',CFG.MinLeafSize,'MaxNumSplits',512);
tree_reg=templateTree('MinLeafSize',CFG.MinLeafSize,'MaxNumSplits',512);
MIX.cls=fitcensemble(F(train_rows,:),is_mixed(train_rows), ...
    'Method','Bag','Learners',tree_cls,'NumLearningCycles',CFG.NumLearningCycles, ...
    'ClassNames',[false true]);
MIX.strong=fitcensemble(F(train_rows,:),is_strong_mixed(train_rows), ...
    'Method','Bag','Learners',tree_cls,'NumLearningCycles',CFG.NumLearningCycles, ...
    'ClassNames',[false true]);
MIX.purity=fitrensemble(F(train_rows,:),patch_purity_target(train_rows), ...
    'Method','Bag','Learners',tree_reg,'NumLearningCycles',CFG.NumLearningCycles);

[p_mixed,p_strong,purity_pred]=predict_mixedness(MIX,F);
purity_pred=min(max(purity_pred,0),1);

%% Train mixedness-aware residual log-k corrector

Fcorr=[F p_mixed p_strong purity_pred];
f_corr_names=[feature_names "p_mixed" "p_strong_mixed" "predicted_patch_purity"];
train_corr=train_rows;
fprintf('Training log-k residual corrector...\n');
CORR.model=fitrensemble(Fcorr(train_corr,:),delta_logk_target(train_corr), ...
    'Method','LSBoost','Learners',tree_reg,'NumLearningCycles',CFG.NumLearningCycles, ...
    'LearnRate',CFG.Shrinkage);

delta_logk_pred=predict(CORR.model,Fcorr);
logk_corrected=log(base_k)+delta_logk_pred;
sws_corrected=2*pi*T.f0./exp(logk_corrected);
sws_corrected=clip_sws(sws_corrected,CFG);
% Blend correction by predicted mixedness: pure-like patches keep more base.
gate=max(p_mixed,0.35*p_strong);
gate(T.geometry_type=="homogeneous")=0;
logk_blend=(1-gate).*log(base_k)+gate.*log(2*pi*T.f0./sws_corrected);
sws_mixedness=clip_sws(2*pi*T.f0./exp(logk_blend),CFG);

%% Train mixedness-aware q-candidate selector

if q_available
    fprintf('Training q-candidate selector from Test 31 q maps...\n');
    [Q_sws,Q_q,q_strategy_names]=q_candidate_matrices(T);
    q_label=best_q_candidate_label(Q_sws,T.true_SWS);
    QSEL.model=fitcensemble(Fcorr(train_rows,:),categorical(q_label(train_rows)), ...
        'Method','Bag','Learners',tree_cls,'NumLearningCycles',CFG.NumLearningCycles);
    q_pred_label=double(string(predict(QSEL.model,Fcorr)));
    q_pred_label(~isfinite(q_pred_label)|q_pred_label<1|q_pred_label>numel(q_strategy_names))=1;
    [sws_q_selector,q_q_selector,q_selector_name]=select_q_candidate(Q_sws,Q_q,q_pred_label,q_strategy_names);
    use_q_selector=(p_mixed>=0.50 | p_strong>=0.25 | T.low_confidence) & T.geometry_type~="homogeneous";
    sws_q_selector(~use_q_selector)=T.sws_local_baseline(~use_q_selector);
    q_q_selector(~use_q_selector)=T.q_local(~use_q_selector);
    q_selector_name(~use_q_selector)="local_baseline";
else
    QSEL=struct();
    sws_q_selector=T.sws_local_baseline;
    q_q_selector=nan(height(T),1);
    q_selector_name=repmat("unavailable",height(T),1);
end
sws_q_selector=clip_sws(sws_q_selector,CFG);

%% Posterior reliability model

corrected_abs_error=abs(100*(sws_mixedness-T.true_SWS)./T.true_SWS);
high_error20_corrected=corrected_abs_error>20;
correction_mag=abs(sws_mixedness-base_sws)./max(base_sws,0.25);
ensemble_spread=estimate_spread(T,base_sws,sws_mixedness);
Freliab=[Fcorr correction_mag ensemble_spread gate delta_logk_pred];
reliab_names=[f_corr_names "relative_correction_magnitude" "ensemble_spread" "mixedness_gate" "delta_logk_pred"];
fprintf('Training posterior reliability detector...\n');
REL.model=fitcensemble(Freliab(train_rows,:),high_error20_corrected(train_rows), ...
    'Method','Bag','Learners',tree_cls,'NumLearningCycles',CFG.NumLearningCycles, ...
    'ClassNames',[false true]);
[~,score_rel]=predict(REL.model,Freliab);
p_high20_corrected=positive_score(REL.model,score_rel);
posterior_reliability=1-p_high20_corrected;

%% Attach outputs

T.predicted_patch_purity=purity_pred;
T.predicted_mixed_probability=p_mixed;
T.predicted_strong_mixed_probability=p_strong;
T.mixedness_gate=gate;
T.delta_logk_pred=delta_logk_pred;
T.sws_mixedness_logk_corrected=sws_mixedness;
T.signed_error_mixedness_logk_corrected=100*(sws_mixedness-T.true_SWS)./T.true_SWS;
T.abs_error_mixedness_logk_corrected=abs(T.signed_error_mixedness_logk_corrected);
T.sws_mixedness_q_candidate_selector=sws_q_selector;
T.q_mixedness_q_candidate_selector=q_q_selector;
T.q_candidate_selected=q_selector_name;
T.signed_error_mixedness_q_candidate_selector=100*(sws_q_selector-T.true_SWS)./T.true_SWS;
T.abs_error_mixedness_q_candidate_selector=abs(T.signed_error_mixedness_q_candidate_selector);
T.posterior_high_error20_probability=p_high20_corrected;
T.posterior_reliability=posterior_reliability;
T.split_role=repmat("unused",height(T),1); T.split_role(train_mask)="train"; T.split_role(test_mask)="test";

%% Evaluation tables

T_eval=T(test_mask,:);
T_mixedness=mixedness_metrics(T_eval,CFG);
T_overall=summarize_strategies(T_eval,"strategy_name",CFG);
T_by_geometry=summarize_strategies(T_eval,["strategy_name","geometry"],CFG);
T_by_frequency=summarize_strategies(T_eval,["strategy_name","geometry","f0"],CFG);
T_by_M=summarize_strategies(T_eval,["strategy_name","M"],CFG);
T_by_purity=summarize_strategies(T_eval,["strategy_name","geometry","purity_bin"],CFG);
T_by_distance=summarize_strategies(T_eval,["strategy_name","geometry","distance_bin"],CFG);
T_by_side=summarize_by_side(T_eval,["strategy_name","geometry","material_side"],CFG);
T_by_side_frequency=summarize_by_side(T_eval,["strategy_name","geometry","material_side","f0"],CFG);
T_roi=roi_summary(T_eval,CFG);
T_reliability=reliability_calibration(T_eval);
T_feature_importance=feature_importance_table(MIX,CORR,REL,feature_names,f_corr_names,reliab_names);

writetable(T,fullfile(OUT.table_dir,'test33_patch_level_predictions.csv'));
writetable(split_table,fullfile(OUT.table_dir,'test33_condition_split.csv'));
writetable(T_mixedness,fullfile(OUT.table_dir,'test33_mixedness_metrics.csv'));
writetable(T_overall,fullfile(OUT.table_dir,'test33_strategy_summary_overall_test.csv'));
writetable(T_by_geometry,fullfile(OUT.table_dir,'test33_strategy_summary_by_geometry_test.csv'));
writetable(T_by_frequency,fullfile(OUT.table_dir,'test33_strategy_summary_by_geometry_frequency_test.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'test33_strategy_summary_by_M_test.csv'));
writetable(T_by_purity,fullfile(OUT.table_dir,'test33_strategy_summary_by_purity_bin_test.csv'));
writetable(T_by_distance,fullfile(OUT.table_dir,'test33_strategy_summary_by_distance_bin_test.csv'));
writetable(T_by_side,fullfile(OUT.table_dir,'test33_strategy_summary_by_soft_hard_test.csv'));
writetable(T_by_side_frequency,fullfile(OUT.table_dir,'test33_strategy_summary_by_soft_hard_frequency_test.csv'));
writetable(T_roi,fullfile(OUT.table_dir,'test33_roi_summary_test.csv'));
writetable(T_reliability,fullfile(OUT.table_dir,'test33_posterior_reliability_calibration.csv'));
writetable(T_feature_importance,fullfile(OUT.table_dir,'test33_feature_importance.csv'));
save(fullfile(OUT.data_dir,'test33_models_and_predictions.mat'), ...
    'MIX','CORR','QSEL','REL','CFG','feature_names','f_corr_names','reliab_names', ...
    'T_overall','T_mixedness','q_available','q_meta','-v7.3');

plot_summary(T_eval,CFG,OUT);
plot_mixedness(T_eval,OUT);
plot_reliability(T_reliability,OUT);
plot_soft_hard_diagnostics(T_by_side,T_by_side_frequency,T_roi,OUT);
plot_distance_side_diagnostics(T_eval,CFG,OUT);
plot_representative_maps(T_eval,CFG,OUT);
if CFG.SaveAllMaps, plot_all_maps(T,CFG,OUT); end
print_interpretation(T_overall,T_mixedness,T_roi,T_reliability,CFG);
fprintf('\nTables: %s\nFigures: %s\nTest 33 complete.\n',OUT.table_dir,OUT.figure_dir);

%% Setup helpers

function SOURCE=resolve_test32_source(root,CFG)
base=fullfile(root,'outputs','test_32_structural_confidence_hybrid');
if CFG.QuickMode
    candidates=[string(fullfile(base,'quick','data','test32_compact_results.mat')), ...
        string(fullfile(base,'data','test32_compact_results.mat'))];
else
    candidates=string(fullfile(base,'data','test32_compact_results.mat'));
end
SOURCE="";
for c=candidates
    if exist(c,'file')==2, SOURCE=c; break; end
end
assert(SOURCE~="",'Missing Test 32 compact results. Run Test 32 first.');
end

function [T,q_available,q_meta]=attach_test31_q_candidates(T,root,CFG)
q_available=false; q_meta=struct('strategy_names',strings(1,0));
base=fullfile(root,'outputs','test_31_simple_confidence_interpolation');
if CFG.QuickMode
    candidates=[string(fullfile(base,'quick','data','test31_compact_results.mat')), ...
        string(fullfile(base,'data','test31_compact_results.mat'))];
else
    candidates=string(fullfile(base,'data','test31_compact_results.mat'));
end
source="";
for c=candidates
    if exist(c,'file')==2, source=c; break; end
end
if source=="", return; end
S=load(source,'T_patch'); Q=S.T_patch; clear S;
keys={'condition_key','map_iz','map_ix'};
needed=[keys {'q_local','q_corr_q_nearest_highconf','q_corr_q_median_highconf_global', ...
    'q_corr_q_median_highconf_by_region_local_structure','q_corr_edgeaware_q_theory_structure', ...
    'q_corr_edgeaware_q_local_structure','sws_q_nearest_highconf','sws_q_median_highconf_global', ...
    'sws_q_median_highconf_by_region_local_structure','sws_edgeaware_q_theory_structure', ...
    'sws_edgeaware_q_local_structure'}];
if ~all(ismember(needed,Q.Properties.VariableNames)), return; end
Q=Q(:,needed);
if height(T)==height(Q) && all(T.condition_key==Q.condition_key) && all(T.map_iz==Q.map_iz) && all(T.map_ix==Q.map_ix)
    for v=setdiff(string(Q.Properties.VariableNames),string(keys),'stable')
        T.(v)=Q.(v);
    end
else
    T=leftjoin(T,Q,'Keys',keys,'MergeKeys',true);
end
q_available=all(ismember(["q_local","q_corr_edgeaware_q_theory_structure"],T.Properties.VariableNames));
q_meta.source=source;
q_meta.strategy_names=["local_baseline","q_nearest_highconf","q_median_highconf_global", ...
    "q_median_highconf_by_region","edgeaware_q_theory_structure","edgeaware_q_local_structure"];
end

function OUT=make_output_dirs(root,CFG)
OUT.root_dir=fullfile(root,'outputs','test_33_mixedness_aware_q_correction');
if CFG.QuickMode, OUT.root_dir=fullfile(OUT.root_dir,'quick'); end
OUT.table_dir=fullfile(OUT.root_dir,'tables'); OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.data_dir=fullfile(OUT.root_dir,'data'); OUT.condition_dir=fullfile(OUT.figure_dir,'maps_by_condition');
for d=string(struct2cell(OUT))', if exist(d,'dir')~=7, mkdir(d); end, end
end

function write_config_json(CFG,file)
C=CFG; for f=string(fieldnames(C))', if isstring(C.(f)), C.(f)=cellstr(C.(f)); end, end
fid=fopen(file,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(C,'PrettyPrint',true),'char'); clear cleanup;
end

function [train_mask,test_mask,split_table]=make_condition_split(T,CFG)
conditions=unique(T(:,{'condition_key','geometry','field_regime','f0','M'}),'rows','stable');
conditions=sortrows(conditions,{'geometry','field_regime','f0','M'});
train_condition=false(height(conditions),1);
for geo=unique(conditions.geometry)'
    idx=find(conditions.geometry==geo);
    n=max(1,round(CFG.TrainFraction*numel(idx)));
    train_condition(idx(1:n))=true;
end
split_table=conditions; split_table.split_role=repmat("test",height(conditions),1); split_table.split_role(train_condition)="train";
[tf,loc]=ismember(T.condition_key,split_table.condition_key); assert(all(tf));
train_mask=split_table.split_role(loc)=="train"; test_mask=~train_mask;
end

%% Features and models

function [F,names]=build_features(T,CFG)
base=T.("sws_"+CFG.BaseStrategy);
local=T.sws_local_baseline; nearest=T.sws_sws_nearest_highconf; region=T.sws_test30_theory_region_levels;
operational_spread=estimate_spread(T,base,base);
vars=[ ...
    T.confidence, double(T.low_confidence), normalize_num(T.f0), normalize_num(T.M), ...
    encode_string(T.field_regime), encode_string(T.geometry_type), ...
    log(max(local,0.1)), log(max(nearest,0.1)), log(max(region,0.1)), log(max(base,0.1)), ...
    nearest-local, region-local, region-nearest, ...
    abs(nearest-local)./max(local,0.25), abs(region-local)./max(local,0.25), ...
    T.estimated_boundary_distance_mm, double(T.estimated_boundary_band), ...
    T.hybrid_structural_weight, T.hybrid_final_weight, double(T.structure_agree), T.region_local_disagreement, ...
    operational_spread ...
    ];
names=["confidence","low_confidence","f0_norm","M_norm","field_regime_code","geometry_type_code", ...
    "log_sws_local","log_sws_nearest","log_sws_region","log_sws_base", ...
    "nearest_minus_local","region_minus_local","region_minus_nearest", ...
    "rel_nearest_local","rel_region_local","estimated_boundary_distance_mm", ...
    "estimated_boundary_band","hybrid_structural_weight","hybrid_final_weight", ...
    "structure_agree","region_local_disagreement","ensemble_spread_operational"];
F=double(vars); F(~isfinite(F))=0;
end

function x=normalize_num(x)
x=double(x); x=(x-mean(x,'omitnan'))./max(std(x,'omitnan'),eps);
end

function code=encode_string(s)
[~,~,code]=unique(string(s)); code=normalize_num(double(code));
end

function spread=estimate_spread(T,base_sws,corrected)
V=[T.sws_local_baseline T.sws_sws_nearest_highconf T.sws_test30_theory_region_levels ...
    T.sws_hybrid_lowconf_region_else_sws_nearest T.sws_hybrid_boundary_protected_region base_sws corrected];
spread=std(V,0,2,'omitnan')./max(mean(abs(V),2,'omitnan'),0.25);
spread(~isfinite(spread))=0;
end

function [p_mixed,p_strong,purity_pred]=predict_mixedness(MIX,F)
[~,score]=predict(MIX.cls,F); p_mixed=positive_score(MIX.cls,score);
[~,score]=predict(MIX.strong,F); p_strong=positive_score(MIX.strong,score);
purity_pred=predict(MIX.purity,F);
end

function [Q_sws,Q_q,names]=q_candidate_matrices(T)
names=["local_baseline","q_nearest_highconf","q_median_highconf_global", ...
    "q_median_highconf_by_region","edgeaware_q_theory_structure","edgeaware_q_local_structure"];
Q_sws=[T.sws_local_baseline T.sws_q_nearest_highconf T.sws_q_median_highconf_global ...
    T.sws_q_median_highconf_by_region_local_structure T.sws_edgeaware_q_theory_structure ...
    T.sws_edgeaware_q_local_structure];
Q_q=[T.q_local T.q_corr_q_nearest_highconf T.q_corr_q_median_highconf_global ...
    T.q_corr_q_median_highconf_by_region_local_structure T.q_corr_edgeaware_q_theory_structure ...
    T.q_corr_edgeaware_q_local_structure];
Q_sws(~isfinite(Q_sws))=NaN; Q_q(~isfinite(Q_q))=NaN;
end

function label=best_q_candidate_label(Q_sws,true_sws)
err=abs(100*(Q_sws-true_sws)./true_sws);
err(~isfinite(err))=inf;
[~,label]=min(err,[],2);
label=double(label);
label(~isfinite(label))=1;
end

function [sws,q,choice]=select_q_candidate(Q_sws,Q_q,label,names)
n=size(Q_sws,1); label=max(1,min(numel(names),round(label(:))));
idx=sub2ind(size(Q_sws),(1:n)',label);
sws=Q_sws(idx); q=Q_q(idx); choice=names(label)';
bad=~isfinite(sws); sws(bad)=Q_sws(bad,1);
badq=~isfinite(q); q(badq)=Q_q(badq,1);
choice(bad)="local_baseline";
end

function p=positive_score(model,score)
classes=string(model.ClassNames); pos=find(classes=="true"|classes=="1",1);
if isempty(pos), pos=size(score,2); end
p=score(:,pos);
p=min(max(p,0),1);
end

function sws=clip_sws(sws,CFG)
sws=min(max(sws,CFG.PhysicalRange(1)),CFG.PhysicalRange(2));
end

%% Summaries

function Tm=mixedness_metrics(T,CFG)
y=T.patch_purity<CFG.MixedThreshold; ys=T.patch_purity<CFG.StrongMixedThreshold;
score=T.predicted_mixed_probability; scores=T.predicted_strong_mixed_probability;
Tm=table();
Tm.metric="mixed_lt_0p95"; Tm.N=height(T); Tm.roc_auc=auc_roc(y,score); Tm.pr_auc=auc_pr(y,score);
Tm.rmse_patch_purity=sqrt(mean((T.predicted_patch_purity-T.patch_purity).^2,'omitnan'));
Tm.mae_patch_purity=mean(abs(T.predicted_patch_purity-T.patch_purity),'omitnan');
for th=[0.25 0.5 0.75]
    pred=score>=th; tag=sprintf('%03d',round(100*th));
    Tm.("precision_"+tag)=sum(pred&y)/max(sum(pred),1);
    Tm.("recall_"+tag)=sum(pred&y)/max(sum(y),1);
end
row=Tm(1,:); row.metric="strong_mixed_lt_0p75"; row.roc_auc=auc_roc(ys,scores); row.pr_auc=auc_pr(ys,scores);
for th=[0.25 0.5 0.75]
    pred=scores>=th; tag=sprintf('%03d',round(100*th));
    row.("precision_"+tag)=sum(pred&ys)/max(sum(pred),1);
    row.("recall_"+tag)=sum(pred&ys)/max(sum(ys),1);
end
Tm=[Tm; row];
end

function S=summarize_strategies(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
for strategy=CFG.EvalStrategies
    sws=strategy_sws(T,strategy);
    error=100*(sws-T.true_SWS)./T.true_SWS;
    if isempty(basegroups), G=ones(height(T),1); else, G=findgroups(T(:,cellstr(basegroups))); end
    for id=unique(G,'stable')'
        idx=G==id; row=table(); for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy; row.N=sum(idx); row.MAPE=mean(abs(error(idx)),'omitnan');
        row.mean_signed_error=mean(error(idx),'omitnan'); row.high_error_20_pct=100*mean(abs(error(idx))>20,'omitnan');
        row.underestimate_pct=100*mean(error(idx)<0,'omitnan');
        for cls=["homogeneous","pure","near_pure","moderately_mixed","strongly_mixed"]
            j=idx&T.purity_bin==cls; token=char(cls); row.("MAPE_"+token)=mean(abs(error(j)),'omitnan'); row.("high20_"+token)=100*mean(abs(error(j))>20,'omitnan');
        end
        for side=["soft","hard"]
            j=idx&T.material_side==side; token=char(side); row.("MAPE_"+token)=mean(abs(error(j)),'omitnan');
        end
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function T=roi_summary(P,CFG)
T=table();
for strategy=CFG.EvalStrategies
    sws=strategy_sws(P,strategy);
    err=100*(sws-P.true_SWS)./P.true_SWS;
    for geo=["bilayer_2_3","circular_inclusion_2_3"]
        for roi=["soft_core","hard_core","interface_0_2mm"]
            switch roi
                case "soft_core", idx=P.geometry==geo&P.material_side=="soft"&P.distance_to_interface_mm>8;
                case "hard_core", idx=P.geometry==geo&P.material_side=="hard"&P.distance_to_interface_mm>4;
                otherwise, idx=P.geometry==geo&P.distance_to_interface_mm<=2;
            end
            row=table(strategy,geo,roi,sum(idx),mean(sws(idx),'omitnan'),std(sws(idx),'omitnan'),mean(abs(err(idx)),'omitnan'),100*mean(abs(err(idx))>20,'omitnan'), ...
                'VariableNames',{'strategy_name','geometry','roi','N','mean_SWS','std_SWS','MAPE','high20_pct'});
            T=concat_tables(T,row);
        end
    end
end
end

function sws=strategy_sws(T,strategy)
if strategy=="mixedness_logk_corrected"
    sws=T.sws_mixedness_logk_corrected;
elseif strategy=="mixedness_q_candidate_selector"
    sws=T.sws_mixedness_q_candidate_selector;
else
    sws=T.("sws_"+strategy);
end
end

function S=summarize_by_side(T,groups,CFG)
S=table(); basegroups=setdiff(groups,"strategy_name",'stable');
valid_side=ismember(string(T.material_side),["soft","hard"]);
T=T(valid_side,:);
if isempty(T), return; end
for strategy=CFG.EvalStrategies
    sws=strategy_sws(T,strategy);
    error=100*(sws-T.true_SWS)./T.true_SWS;
    G=findgroups(T(:,cellstr(basegroups)));
    for id=unique(G,'stable')'
        idx=G==id; row=table(); for v=basegroups, row.(v)=T.(v)(find(idx,1)); end
        row.strategy_name=strategy;
        row.N=sum(idx);
        row.mean_SWS=mean(sws(idx),'omitnan');
        row.std_SWS=std(sws(idx),'omitnan');
        row.MAPE=mean(abs(error(idx)),'omitnan');
        row.mean_signed_error=mean(error(idx),'omitnan');
        row.median_signed_error=median(error(idx),'omitnan');
        row.high_error_10_pct=100*mean(abs(error(idx))>10,'omitnan');
        row.high_error_20_pct=100*mean(abs(error(idx))>20,'omitnan');
        row.underestimate_pct=100*mean(error(idx)<0,'omitnan');
        S=concat_tables(S,row);
    end
end
S=movevars(S,'strategy_name','Before',1);
end

function T=reliability_calibration(P)
edges=0:0.1:1; T=table();
for b=1:numel(edges)-1
    idx=P.posterior_reliability>=edges(b)&P.posterior_reliability<edges(b+1);
    if b==numel(edges)-1, idx=P.posterior_reliability>=edges(b)&P.posterior_reliability<=edges(b+1); end
    err=P.abs_error_mixedness_logk_corrected;
    row=table(mean(edges(b:b+1)),sum(idx),mean(P.posterior_reliability(idx),'omitnan'),mean(err(idx),'omitnan'),100*mean(err(idx)>20,'omitnan'), ...
        'VariableNames',{'reliability_bin_center','N','mean_predicted_reliability','observed_MAPE','observed_high20_pct'});
    T=concat_tables(T,row);
end
end

function T=feature_importance_table(MIX,CORR,REL,feature_names,f_corr_names,reliab_names)
T=table();
T=concat_tables(T,importance_rows("mixed_classifier",feature_names,predictorImportance(MIX.cls)));
T=concat_tables(T,importance_rows("purity_regressor",feature_names,predictorImportance(MIX.purity)));
T=concat_tables(T,importance_rows("logk_corrector",f_corr_names,predictorImportance(CORR.model)));
T=concat_tables(T,importance_rows("posterior_reliability",reliab_names,predictorImportance(REL.model)));
end

function T=importance_rows(model_name,names,importance)
T=table(repmat(model_name,numel(names),1),names(:),importance(:),'VariableNames',{'model_name','feature_name','importance'});
T=sortrows(T,'importance','descend');
end

function a=auc_roc(y,score)
y=logical(y); [~,ord]=sort(score,'descend'); y=y(ord); P=sum(y); N=sum(~y);
if P==0||N==0, a=NaN; return; end
tp=cumsum(y)/P; fp=cumsum(~y)/N; a=trapz([0;fp;1],[0;tp;1]);
end

function a=auc_pr(y,score)
y=logical(y); [~,ord]=sort(score,'descend'); y=y(ord); P=sum(y);
if P==0, a=NaN; return; end
tp=cumsum(y); fp=cumsum(~y); recall=tp/P; precision=tp./max(tp+fp,1); a=trapz([0;recall],[1;precision]);
end

%% Figures

function plot_summary(T,CFG,OUT)
S=summarize_strategies(T,"strategy_name",CFG);
fig=figure('Color','w','Position',[70 70 1150 520]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); barh(ax,S.MAPE); set(ax,'YTick',1:height(S),'YTickLabel',pretty(S.strategy_name),'TickLabelInterpreter','none'); xlabel(ax,'MAPE (%)'); title(ax,'Held-out strategies'); style_axes(ax);
ax=nexttile(tl); barh(ax,S.high_error_20_pct); set(ax,'YTick',1:height(S),'YTickLabel',pretty(S.strategy_name),'TickLabelInterpreter','none'); xlabel(ax,'Error >20% (%)'); title(ax,'Large-error rate'); style_axes(ax);
export_fig(fig,OUT,'test33_strategy_summary_test.png');
end

function plot_mixedness(T,OUT)
fig=figure('Color','w','Position',[80 80 1150 480]); tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl); scatter(ax,T.patch_purity,T.predicted_patch_purity,6,T.confidence,'filled','MarkerFaceAlpha',.12); hold(ax,'on'); plot(ax,[0 1],[0 1],'k--'); xlabel(ax,'True patch purity'); ylabel(ax,'Predicted patch purity'); colorbar(ax); style_axes(ax);
ax=nexttile(tl); scatter(ax,T.patch_purity,T.predicted_mixed_probability,6,T.region_local_disagreement,'filled','MarkerFaceAlpha',.12); xlabel(ax,'True patch purity'); ylabel(ax,'Predicted P(mixed)'); colorbar(ax); style_axes(ax);
export_fig(fig,OUT,'test33_mixedness_predictions.png');
end

function plot_reliability(T,OUT)
fig=figure('Color','w','Position',[100 100 700 430]); ax=axes(fig);
bar(ax,T.reliability_bin_center,T.observed_high20_pct); xlabel(ax,'Predicted reliability bin'); ylabel(ax,'Observed error >20% (%)'); title(ax,'Posterior reliability calibration'); style_axes(ax);
export_fig(fig,OUT,'test33_posterior_reliability_calibration.png');
end

function plot_soft_hard_diagnostics(T_side,T_side_freq,T_roi,OUT)
if isempty(T_side), return; end
main=["local_baseline","test30_theory_region_levels","hybrid_lowconf_region_else_sws_nearest", ...
    "mixedness_logk_corrected","mixedness_q_candidate_selector"];
T_side=T_side(ismember(T_side.strategy_name,main),:);
T_side_freq=T_side_freq(ismember(T_side_freq.strategy_name,main),:);
fig=figure('Color','w','Position',[40 40 1450 650]); tl=tiledlayout(2,2,'TileSpacing','compact');
for metric=["MAPE","mean_signed_error"]
    ax=nexttile(tl); hold(ax,'on');
    cats=strcat(string(T_side.geometry)," / ",string(T_side.material_side));
    [G,labels]=findgroups(cats);
    width=0.8/numel(main);
    for s=1:numel(main)
        idx=T_side.strategy_name==main(s);
        vals=splitapply(@(x)mean(x,'omitnan'),T_side.(metric)(idx),G(idx));
        x=(1:numel(labels))+(s-(numel(main)+1)/2)*width;
        bar(ax,x,vals,width,'DisplayName',pretty(main(s)));
    end
    set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',25,'TickLabelInterpreter','none');
    ylabel(ax,replace(metric,"_"," ")); title(ax,replace(metric,"_"," ")+" by geometry/side"); legend(ax,'Location','bestoutside'); style_axes(ax);
end
ax=nexttile(tl); hold(ax,'on');
for s=1:numel(main)
    idx=T_roi.strategy_name==main(s)&ismember(T_roi.roi,["soft_core","hard_core","interface_0_2mm"]);
    cats=strcat(string(T_roi.geometry(idx))," / ",string(T_roi.roi(idx)));
    [G,labels]=findgroups(cats);
    vals=splitapply(@(x)mean(x,'omitnan'),T_roi.MAPE(idx),G);
    x=(1:numel(labels))+(s-(numel(main)+1)/2)*(0.8/numel(main));
    bar(ax,x,vals,0.8/numel(main),'DisplayName',pretty(main(s)));
end
set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',25,'TickLabelInterpreter','none');
ylabel(ax,'MAPE (%)'); title(ax,'Core/interface ROI MAPE'); legend(ax,'Location','bestoutside'); style_axes(ax);
ax=nexttile(tl); hold(ax,'on');
for s=1:numel(main)
    idx=T_side_freq.strategy_name==main(s)&T_side_freq.geometry=="circular_inclusion_2_3"&T_side_freq.material_side=="hard";
    plot(ax,T_side_freq.f0(idx),T_side_freq.MAPE(idx),'-o','DisplayName',pretty(main(s)));
end
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'MAPE (%)'); title(ax,'Inclusion hard side vs frequency'); legend(ax,'Location','best'); style_axes(ax);
export_fig(fig,OUT,'test33_soft_hard_error_diagnostics.png');
end

function plot_distance_side_diagnostics(T,~,OUT)
main=["local_baseline","test30_theory_region_levels","hybrid_lowconf_region_else_sws_nearest", ...
    "mixedness_logk_corrected","mixedness_q_candidate_selector"];
dist_order=["0_1mm","1_2mm","2_4mm","4_8mm","gt_8mm"];
fig=figure('Color','w','Position',[55 55 1450 720]); tl=tiledlayout(2,2,'TileSpacing','compact');
for geo=["bilayer_2_3","circular_inclusion_2_3"]
    for side=["soft","hard"]
        ax=nexttile(tl); hold(ax,'on');
        for strategy=main
            sws=strategy_sws(T,strategy);
            err=abs(100*(sws-T.true_SWS)./T.true_SWS);
            vals=nan(size(dist_order));
            for b=1:numel(dist_order)
                idx=T.geometry==geo&T.material_side==side&T.distance_bin==dist_order(b);
                vals(b)=mean(err(idx),'omitnan');
            end
            plot(ax,1:numel(dist_order),vals,'-o','DisplayName',pretty(strategy));
        end
        set(ax,'XTick',1:numel(dist_order),'XTickLabel',replace(dist_order,"_","-"),'TickLabelInterpreter','none');
        xlabel(ax,'Distance to interface bin'); ylabel(ax,'MAPE (%)');
        title(ax,replace(geo,"_"," ")+" / "+side); legend(ax,'Location','best'); style_axes(ax);
    end
end
export_fig(fig,OUT,'test33_error_vs_distance_soft_hard.png');
end

function plot_representative_maps(T,CFG,OUT)
keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable'); if isempty(keys), return; end
plot_case_maps(T(T.condition_key==keys(1),:),CFG,fullfile(OUT.figure_dir,'test33_representative_maps.png'));
end

function plot_all_maps(T,CFG,OUT)
keys=unique(T.condition_key,'stable');
for i=1:numel(keys)
    C=T(T.condition_key==keys(i),:); folder=fullfile(OUT.condition_dir,sanitize(C.geometry(1)),sanitize(C.field_regime(1))); if exist(folder,'dir')~=7, mkdir(folder); end
    plot_case_maps(C,CFG,fullfile(folder,"test33__"+sanitize(keys(i))+".png"));
end
end

function plot_case_maps(C,~,file)
C=sortrows(C,{'map_iz','map_ix'});
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,C.sws_local_baseline),map_from_rows(C,C.("sws_"+"hybrid_lowconf_region_else_sws_nearest")), ...
    map_from_rows(C,C.sws_mixedness_logk_corrected),map_from_rows(C,C.abs_error_local_baseline), ...
    map_from_rows(C,C.abs_error_hybrid_lowconf_region_else_sws_nearest),map_from_rows(C,C.abs_error_mixedness_logk_corrected), ...
    map_from_rows(C,C.sws_mixedness_q_candidate_selector),map_from_rows(C,C.abs_error_mixedness_q_candidate_selector), ...
    map_from_rows(C,C.predicted_mixed_probability),map_from_rows(C,C.predicted_patch_purity),map_from_rows(C,C.mixedness_gate), ...
    map_from_rows(C,C.delta_logk_pred),map_from_rows(C,C.posterior_reliability)};
titles=["True SWS","Local","Base hybrid","Mixedness log-k","Local error","Base error","Log-k error", ...
    "Mixedness q selector","q-selector error","P(mixed)","Predicted purity","Mixedness gate","Delta log-k","Posterior reliability"];
fig=figure('Color','w','Visible','off','Position',[20 20 1600 900]); tl=tiledlayout(4,4,'TileSpacing','compact');
for i=1:numel(maps), ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',8); colorbar(ax); end
title(tl,C.condition_key(1),'Interpreter','none'); exportgraphics(fig,file,'Resolution',170); close(fig);
end

function print_interpretation(S,M,R,C,CFG)
[~,i]=min(S.MAPE); fprintf('\nInterpretation on held-out conditions:\n');
fprintf('  Best global strategy: %s (MAPE %.2f%%, >20 %.2f%%).\n',S.strategy_name(i),S.MAPE(i),S.high_error_20_pct(i));
row=S(S.strategy_name=="mixedness_logk_corrected",:);
base=S(S.strategy_name==CFG.BaseStrategy,:);
if ~isempty(row)&&~isempty(base)
    fprintf('  Mixedness log-k corrected: MAPE %.2f%% vs base %.2f%%; strong-mixed %.2f%% vs %.2f%%.\n', ...
        row.MAPE,base.MAPE,row.MAPE_strongly_mixed,base.MAPE_strongly_mixed);
end
fprintf('  Mixedness detector ROC AUC %.3f, PR AUC %.3f, purity MAE %.3f.\n',M.roc_auc(1),M.pr_auc(1),M.mae_patch_purity(1));
iface=R(R.roi=="interface_0_2mm"&R.strategy_name=="mixedness_logk_corrected",:);
if ~isempty(iface), fprintf('  Corrected interface ROI mean MAPE %.2f%%.\n',mean(iface.MAPE,'omitnan')); end
valid=C(C.N>0,:); if ~isempty(valid), fprintf('  Reliability bins show high20 from %.1f%% to %.1f%% across occupied bins.\n',min(valid.observed_high20_pct),max(valid.observed_high20_pct)); end
end

function validate_test33(T,CFG,q_available)
[F,names]=build_features(T,CFG); assert(size(F,1)==height(T)&&numel(names)==size(F,2));
assert(all(isfinite(F(:)))); fprintf('  Feature construction passed with %d features.\n',numel(names));
if q_available
    [Q_sws,Q_q,q_names]=q_candidate_matrices(T);
    assert(size(Q_sws,1)==height(T)&&size(Q_sws,2)==numel(q_names));
    assert(size(Q_q,1)==height(T)&&size(Q_q,2)==numel(q_names));
    fprintf('  q-candidate construction passed with %d candidates.\n',numel(q_names));
end
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
function label=pretty(name), label=replace(string(name),"_"," "); end
function value=sanitize(value), value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_'); end
function style_axes(ax), grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off'); ax.FontSize=9; end
function export_fig(fig,OUT,name), exportgraphics(fig,fullfile(OUT.figure_dir,name),'Resolution',210); close(fig); end
