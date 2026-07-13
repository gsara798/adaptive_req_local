%% analyze_test_21_interface_confidence_detector.m
% Test 21: interface confidence detector for adaptive REQ.
%
% This analysis does not modify or retrain any q-estimation model. It learns
% whether the primary Hybrid T18 SWS estimate is likely to exceed 20% error,
% using only operational predictors available at inference time.
%
% Local map statistics use a 3-by-3 neighborhood. mixed_window_fraction,
% distance to interface, and true SWS are diagnostic-only variables.

clear; clc; close all;
format compact;

%% Setup

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontName','Helvetica', ...
    'defaultTextFontName','Helvetica', ...
    'defaultAxesFontSize',8,'defaultTextFontSize',8, ...
    'defaultLegendFontSize',7);

OUT=make_output_dirs(root_dir);
SOURCE=fullfile(root_dir,'outputs', ...
    'test_20_external_validation_and_aperture_q_tracking','analysis', ...
    'tables','level20_external_predictions.csv');
assert(exist(SOURCE,'file')==2,'Test 20 prediction table not found:\n%s',SOURCE);

PRIMARY_MODEL="HybridLocalGlobal_T18_noUserRegime";
RANDOM_SEED=21001;
TEST_FRACTION=0.25;
LOCAL_NEIGHBORHOOD=3;
THRESHOLDS=[0.10 0.20 0.30 0.50];

fprintf('\nTest 21: Interface confidence detector\n');
fprintf('Source: %s\n',SOURCE);
fprintf('Primary q/SWS model: %s\n',PRIMARY_MODEL);
fprintf('Local map neighborhood: %dx%d\n',LOCAL_NEIGHBORHOOD,LOCAL_NEIGHBORHOOD);

%% Build or reuse the one-row-per-window confidence dataset

dataset_mat=fullfile(OUT.data_dir,'level21_confidence_dataset.mat');
if exist(dataset_mat,'file')==2
    fprintf('\nLoading cached confidence dataset:\n%s\n',dataset_mat);
    S=load(dataset_mat,'T_confidence','PREDICTOR_INFO');
    T_confidence=S.T_confidence;
    PREDICTOR_INFO=S.PREDICTOR_INFO;
else
    fprintf('\nReading selected Test 20 columns...\n');
    T_long=read_test20_predictions(SOURCE);
    fprintf('Long predictions: %d rows.\n',height(T_long));

    [T_confidence,PREDICTOR_INFO]=build_confidence_dataset( ...
        T_long,root_dir,PRIMARY_MODEL,LOCAL_NEIGHBORHOOD);
    clear T_long;
    save(dataset_mat,'T_confidence','PREDICTOR_INFO','SOURCE', ...
        'PRIMARY_MODEL','LOCAL_NEIGHBORHOOD','-v7.3');
end

assert(~isempty(T_confidence),'The Test 21 confidence dataset is empty.');
fprintf('Confidence rows: %d | columns: %d\n',height(T_confidence),width(T_confidence));
fprintf('High-error >20%% prevalence: %.4f%% (%d/%d)\n', ...
    100*mean(T_confidence.high_error_gt20),sum(T_confidence.high_error_gt20), ...
    height(T_confidence));
fprintf('High-error >10%% prevalence: %.4f%%\n', ...
    100*mean(T_confidence.high_error_gt10));

%% Strict anti-leakage predictor policy

map_features=["q_pred_hybrid_T18";"sws_pred_hybrid_T18"; ...
    "grad_q_hybrid";"grad_sws_hybrid"; ...
    "local_std_q_hybrid";"local_std_sws_hybrid"; ...
    "local_range_q_hybrid";"local_range_sws_hybrid"];
disagreement_features=[ ...
    "abs_q_hybrid_minus_local";"abs_q_hybrid_minus_global"; ...
    "abs_q_hybrid_minus_theory";"abs_q_hybrid_minus_hybrid_guess"; ...
    "abs_q_hybrid_T18_minus_old"; ...
    "abs_sws_hybrid_minus_local";"abs_sws_hybrid_minus_global"; ...
    "abs_sws_hybrid_minus_theory";"abs_sws_hybrid_T18_minus_old"];
operational_metadata=["REQ_M";"SIM_f0";"cs_guess";"dx";"dz"];
ML_PREDICTORS=unique([PREDICTOR_INFO.spectral_predictors(:); ...
    operational_metadata;map_features;disagreement_features],'stable');
ML_PREDICTORS=ML_PREDICTORS(ismember(ML_PREDICTORS, ...
    string(T_confidence.Properties.VariableNames)));
assert_no_leakage(ML_PREDICTORS);

PREDICTOR_POLICY=table(ML_PREDICTORS,repmat("operational_predictor", ...
    numel(ML_PREDICTORS),1),'VariableNames',{'variable_name','policy'});
writetable(PREDICTOR_POLICY,fullfile(OUT.table_dir, ...
    'level21_predictor_policy.csv'));
fprintf('Operational ML predictors: %d\n',numel(ML_PREDICTORS));

%% Grouped train/test split

[train_mask,test_mask,group_key]=grouped_condition_split( ...
    T_confidence,TEST_FRACTION,RANDOM_SEED);
T_confidence.group_key=group_key;
T_confidence.split=repmat("train",height(T_confidence),1);
T_confidence.split(test_mask)="test";
assert(isempty(intersect(unique(group_key(train_mask)), ...
    unique(group_key(test_mask)))),'Grouped split leakage detected.');

fprintf('\nGrouped split: %d train rows, %d test rows.\n', ...
    sum(train_mask),sum(test_mask));
fprintf('Physical condition groups: %d train, %d test.\n', ...
    numel(unique(group_key(train_mask))),numel(unique(group_key(test_mask))));
fprintf('Positive prevalence: train %.4f%% | test %.4f%%.\n', ...
    100*mean(T_confidence.high_error_gt20(train_mask)), ...
    100*mean(T_confidence.high_error_gt20(test_mask)));

%% Predictor matrix and balanced weights

X=double(T_confidence{:,cellstr(ML_PREDICTORS)});
y=logical(T_confidence.high_error_gt20);
[X,impute_median,predictor_mu,predictor_sigma]=preprocess_from_training(X,train_mask);
weights=balanced_class_weights(y(train_mask));

%% A/B: calibrated rule-based detectors

detector_checkpoint=fullfile(OUT.model_dir,'level21_detector_training_checkpoint.mat');
if exist(detector_checkpoint,'file')==2
    Sdet=load(detector_checkpoint,'DETECTORS','checkpoint_predictors');
    assert(isequal(string(Sdet.checkpoint_predictors(:)),ML_PREDICTORS(:)), ...
        'Detector checkpoint predictors do not match the current analysis.');
    DETECTORS=Sdet.DETECTORS;
    fprintf('Reusing trained detector checkpoint:\n%s\n',detector_checkpoint);
else
RULES=[ ...
    struct('name',"Rule_q_gradient",'feature',"grad_q_hybrid",'family',"rule_gradient")
    struct('name',"Rule_q_hybrid_minus_local", ...
        'feature',"abs_q_hybrid_minus_local",'family',"rule_disagreement")
    struct('name',"Rule_q_hybrid_minus_theory", ...
        'feature',"abs_q_hybrid_minus_theory",'family',"rule_disagreement")
    struct('name',"Rule_max_local_theory_disagreement", ...
        'feature',"max_local_theory_disagreement",'family',"rule_disagreement")];

T_confidence.max_local_theory_disagreement=max( ...
    T_confidence.abs_q_hybrid_minus_local, ...
    T_confidence.abs_q_hybrid_minus_theory);

DETECTORS=struct('detector_name',{},'detector_family',{}, ...
    'model_type',{},'model',{},'predictors',{},'impute_median',{}, ...
    'predictor_mu',{},'predictor_sigma',{},'rule_feature',{});
for i=1:numel(RULES)
    raw=double(T_confidence.(RULES(i).feature));
    raw_med=median(raw(train_mask),'omitnan');
    raw(~isfinite(raw))=raw_med;
    calibrator=fitclinear(raw(train_mask),y(train_mask), ...
        'Learner','logistic','Regularization','ridge','Lambda',1e-4, ...
        'Weights',weights,'ClassNames',[false true]);
    DETECTORS(end+1)=make_detector(RULES(i).name,RULES(i).family, ...
        "calibrated_rule",calibrator,RULES(i).feature,raw_med,0,1, ...
        RULES(i).feature); %#ok<SAGROW>
end

%% C: ML confidence detectors

opts=statset('UseParallel',true);
logistic=fitclinear(X(train_mask,:),y(train_mask), ...
    'Learner','logistic','Regularization','ridge','Lambda',1e-4, ...
    'Weights',weights,'ClassNames',[false true]);
DETECTORS(end+1)=make_detector("ML_logistic_regression","ml", ...
    "logistic_regression",logistic,ML_PREDICTORS,impute_median, ...
    predictor_mu,predictor_sigma,"");

tree_bag=templateTree('MinLeafSize',20,'MaxNumSplits',512, ...
    'Reproducible',true);
bagged=fitcensemble(X(train_mask,:),y(train_mask), ...
    'Method','Bag','NumLearningCycles',150,'Learners',tree_bag, ...
    'Weights',weights,'ClassNames',[false true],'Options',opts);
DETECTORS(end+1)=make_detector("ML_bagged_trees","ml", ...
    "bagged_trees",bagged,ML_PREDICTORS,impute_median, ...
    predictor_mu,predictor_sigma,"");

tree_boost=templateTree('MinLeafSize',20,'MaxNumSplits',128);
boosted=fitcensemble(X(train_mask,:),y(train_mask), ...
    'Method','LogitBoost','NumLearningCycles',150,'LearnRate',0.05, ...
    'Learners',tree_boost,'Weights',weights, ...
    'ClassNames',[false true]);
DETECTORS(end+1)=make_detector("ML_boosted_trees","ml", ...
    "boosted_trees",boosted,ML_PREDICTORS,impute_median, ...
    predictor_mu,predictor_sigma,"");

checkpoint_predictors=ML_PREDICTORS;
save(detector_checkpoint,'DETECTORS','checkpoint_predictors', ...
    'RANDOM_SEED','-v7.3');
end

%% Held-out predictions and metrics

T_predictions=table();
T_metrics=table();
T_threshold=table();
test_rows=find(test_mask);
for i=1:numel(DETECTORS)
    p=predict_detector(DETECTORS(i),T_confidence,X,test_mask);
    Pi=prediction_table(T_confidence(test_mask,:),p,DETECTORS(i));
    T_predictions=concat_tables(T_predictions,Pi);

    Mi=evaluate_detector(Pi.high_error_gt20,Pi.risk_probability,0.20);
    Mi.detector_name=DETECTORS(i).detector_name;
    Mi.detector_family=DETECTORS(i).detector_family;
    Mi.model_type=DETECTORS(i).model_type;
    T_metrics=concat_tables(T_metrics,Mi);

    for ti=1:numel(THRESHOLDS)
        Mt=evaluate_detector(Pi.high_error_gt20,Pi.risk_probability,THRESHOLDS(ti));
        Mt.detector_name=DETECTORS(i).detector_name;
        Mt.detector_family=DETECTORS(i).detector_family;
        Mt.model_type=DETECTORS(i).model_type;
        T_threshold=concat_tables(T_threshold,Mt);
    end
end
T_metrics=movevars(T_metrics,{'detector_name','detector_family','model_type'}, ...
    'Before',1);
T_threshold=movevars(T_threshold, ...
    {'detector_name','detector_family','model_type'},'Before',1);

T_by_case=evaluate_by_group(T_predictions,'external_case',0.20);
T_by_regime=evaluate_by_group(T_predictions,'field_regime',0.20);
T_by_M=evaluate_by_group(T_predictions,'REQ_M',0.20);

%% Feature importance and best deployable ML detector

T_importance=table();
for i=1:numel(DETECTORS)
    T_importance=concat_tables(T_importance,detector_importance(DETECTORS(i)));
end

ml_idx=find([DETECTORS.detector_family]=="ml");
ml_metrics=T_metrics(ismember(T_metrics.detector_name, ...
    string({DETECTORS(ml_idx).detector_name})),:);
ml_metrics=sortrows(ml_metrics,{'Recall','AUC_PR'},{'descend','descend'});
best_name=ml_metrics.detector_name(1);
best_idx=find([DETECTORS.detector_name]==best_name,1);
BEST=DETECTORS(best_idx);

best_all_probability=predict_detector(BEST,T_confidence,X,true(height(T_confidence),1));
T_confidence.best_detector_probability=best_all_probability;
T_confidence.probability_high_error=best_all_probability;
T_confidence.best_detector_confidence=1-best_all_probability;

%% Save tables and trained model

writetable(T_confidence,fullfile(OUT.table_dir,'level21_confidence_dataset.csv'));
writetable(T_metrics,fullfile(OUT.table_dir,'level21_detector_metrics.csv'));
writetable(T_by_case,fullfile(OUT.table_dir,'level21_detector_metrics_by_case.csv'));
writetable(T_by_regime,fullfile(OUT.table_dir,'level21_detector_metrics_by_regime.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'level21_detector_metrics_by_M.csv'));
writetable(T_predictions,fullfile(OUT.table_dir,'level21_predictions.csv'));
writetable(T_threshold,fullfile(OUT.table_dir,'level21_threshold_sweep.csv'));
writetable(T_importance,fullfile(OUT.table_dir,'level21_feature_importance.csv'));

CONFIDENCE_MODEL=struct();
CONFIDENCE_MODEL.model_name="ConfidenceDetector_T21_HybridT18";
CONFIDENCE_MODEL.detector_name=BEST.detector_name;
CONFIDENCE_MODEL.model_type=BEST.model_type;
CONFIDENCE_MODEL.target="high_error_gt20";
CONFIDENCE_MODEL.base_q_model=PRIMARY_MODEL;
CONFIDENCE_MODEL.predictors=BEST.predictors;
CONFIDENCE_MODEL.impute_median=BEST.impute_median;
CONFIDENCE_MODEL.predictor_mu=BEST.predictor_mu;
CONFIDENCE_MODEL.predictor_sigma=BEST.predictor_sigma;
CONFIDENCE_MODEL.classifier=BEST.model;
CONFIDENCE_MODEL.risk_definition="abs_sws_error_pct > 20";
CONFIDENCE_MODEL.recommended_threshold=recommend_threshold( ...
    T_threshold(T_threshold.detector_name==best_name,:));
CONFIDENCE_MODEL.local_neighborhood=LOCAL_NEIGHBORHOOD;
CONFIDENCE_MODEL.anti_leakage_policy=PREDICTOR_POLICY;
model_file=fullfile(OUT.model_dir,'ConfidenceDetector_T21_HybridT18.mat');
save(model_file,'CONFIDENCE_MODEL','T_metrics','T_threshold', ...
    'PREDICTOR_INFO','-v7.3');

best_metric=T_metrics(T_metrics.detector_name==best_name,:);
notes=sprintf(['Grouped-condition external validation; AUC PR %.4g; ', ...
    'recall@0.20 %.4g; precision@0.20 %.4g.'], ...
    best_metric.AUC_PR(1),best_metric.Recall(1),best_metric.Precision(1));
REGISTRY_ENTRY=adaptive_req.analysis.register_trained_model( ...
    'RootDir',root_dir,'SourceModelFile',model_file, ...
    'ModelId',"test21__ConfidenceDetector_T21_HybridT18__"+BEST.model_type, ...
    'RegistrySubdir','test21_interface_confidence_detector', ...
    'TestName','test_21_interface_confidence_detector', ...
    'AnalysisLevel','analysis','ModelName','ConfidenceDetector_T21_HybridT18', ...
    'FeatureSet','operational_spectral_disagreement_map', ...
    'ModelType',BEST.model_type,'ModelRole','confidence_detector', ...
    'TrainingDataset','test_20_external_validation', ...
    'Target','high_error_gt20', ...
    'PredictorSummary',strjoin(BEST.predictors,', '), ...
    'MetricsFile',fullfile(OUT.table_dir,'level21_detector_metrics.csv'), ...
    'SplitType','grouped_external_condition_75_25', ...
    'PerformanceSummary',notes, ...
    'Notes','Predicts risk only; does not modify or predict q.');
writetable(REGISTRY_ENTRY,fullfile(OUT.table_dir, ...
    'level21_model_registry_entry.csv'));

save(fullfile(OUT.data_dir,'level21_analysis_results.mat'), ...
    'DETECTORS','BEST','T_metrics','T_threshold','T_by_case', ...
    'T_by_regime','T_by_M','T_importance','REGISTRY_ENTRY','-v7.3');

%% Figures

plot_roc_curves(T_predictions,OUT);
plot_pr_curves(T_predictions,OUT);
plot_confidence_calibration(T_predictions,OUT);
plot_error_vs_confidence(T_predictions,OUT);
plot_feature_importance(T_importance,best_name,OUT);
plot_rule_ml_comparison(T_metrics,OUT);
plot_representative_confidence_maps(T_confidence,best_name,OUT);
plot_interface_confidence_maps(T_confidence,best_name,OUT);

%% Scientific summary

print_summary(T_confidence,T_predictions,T_metrics,T_threshold, ...
    T_by_case,T_by_regime,T_by_M,best_name);
fprintf('\nTest 21 complete. Analysis folder:\n%s\n',OUT.analysis_dir);

%% Local functions

function OUT=make_output_dirs(root_dir)
OUT.analysis_dir=fullfile(root_dir,'outputs', ...
    'test_21_interface_confidence_detector','analysis');
OUT.table_dir=fullfile(OUT.analysis_dir,'tables');
OUT.figure_dir=fullfile(OUT.analysis_dir,'figures');
OUT.data_dir=fullfile(OUT.analysis_dir,'data');
OUT.model_dir=fullfile(OUT.analysis_dir,'models');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs), if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end, end
end

function T=read_test20_predictions(file)
keep=["condition_id","condition_label","geometry_id","geometry_type", ...
    "field_regime_label","field_regime_variant","SIM_f0","SIM_dx","SIM_dz", ...
    "REQ_M","REQ_cs_guess","realization_idx","map_iz","map_ix", ...
    "x_center_m","z_center_m","distance_to_interface_m","cs_true", ...
    "model_name","q_pred","cs_pred","sws_error_pct","abs_sws_error_pct"];
opts=detectImportOptions(file,'TextType','string','VariableNamingRule','preserve');
opts.SelectedVariableNames=cellstr(keep);
string_vars=intersect(["condition_label","geometry_id","geometry_type", ...
    "field_regime_label","field_regime_variant","model_name"],keep);
opts=setvartype(opts,cellstr(string_vars),'string');
T=readtable(file,opts);
expected=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess", ...
    "LocalOnly_old","GlobalOnly_old","HybridLocalGlobal_old", ...
    "TheoryQDiscrete"];
assert(all(ismember(expected,unique(T.model_name))), ...
    'Test 20 table is missing one or more required models.');
end

function [D,INFO]=build_confidence_dataset(T,root_dir,primary,neighborhood)
model_names=["LocalOnly_T18","GlobalOnly_T18",primary, ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","LocalOnly_old", ...
    "GlobalOnly_old","HybridLocalGlobal_old","TheoryQDiscrete"];
q_names=["q_pred_local_T18","q_pred_global_T18","q_pred_hybrid_T18", ...
    "q_pred_hybrid_T18_with_guess","q_pred_local_old", ...
    "q_pred_global_old","q_pred_hybrid_old","q_theory_discrete"];
sws_names=["sws_pred_local_T18","sws_pred_global_T18","sws_pred_hybrid_T18", ...
    "sws_pred_hybrid_T18_with_guess","sws_pred_local_old", ...
    "sws_pred_global_old","sws_pred_hybrid_old","sws_theory_discrete"];

H=sortrows(T(T.model_name==primary,:),{'condition_id','map_iz','map_ix'});
D=table();
D.condition_id=H.condition_id;
D.external_case=H.geometry_id;
D.external_case(D.external_case=="inclusion_2_3")="circular_inclusion_2_3";
D.field_regime=H.field_regime_label;
D.regime_variant=H.field_regime_variant;
D.REQ_M=H.REQ_M; D.SIM_f0=H.SIM_f0;
D.cs_true=H.cs_true; D.cs_guess=H.REQ_cs_guess;
D.dx=H.SIM_dx; D.dz=H.SIM_dz;
D.x=H.x_center_m; D.z=H.z_center_m;
D.map_iz=H.map_iz; D.map_ix=H.map_ix;
D.realization_idx=H.realization_idx;
D.true_cs=H.cs_true;
D.sws_error_pct=H.sws_error_pct;
D.abs_sws_error_pct=H.abs_sws_error_pct;
D.high_error_gt20=D.abs_sws_error_pct>20;
D.high_error_gt10=D.abs_sws_error_pct>10;
D.distance_to_interface_mm=1e3*H.distance_to_interface_m;
D.near_interface=isfinite(D.distance_to_interface_mm) & ...
    D.distance_to_interface_mm<=2;

key_ref=H(:,{'condition_id','map_iz','map_ix'});
for i=1:numel(model_names)
    Ti=sortrows(T(T.model_name==model_names(i),:), ...
        {'condition_id','map_iz','map_ix'});
    assert(height(Ti)==height(H) && isequal(Ti(:,{'condition_id','map_iz','map_ix'}),key_ref), ...
        'Pixel alignment failed for model %s.',model_names(i));
    D.(q_names(i))=Ti.q_pred;
    D.(sws_names(i))=Ti.cs_pred;
end
clear T H;

D.abs_q_hybrid_minus_local=abs(D.q_pred_hybrid_T18-D.q_pred_local_T18);
D.abs_q_hybrid_minus_global=abs(D.q_pred_hybrid_T18-D.q_pred_global_T18);
D.abs_q_hybrid_minus_theory=abs(D.q_pred_hybrid_T18-D.q_theory_discrete);
D.abs_q_hybrid_minus_hybrid_guess=abs( ...
    D.q_pred_hybrid_T18-D.q_pred_hybrid_T18_with_guess);
D.abs_q_hybrid_T18_minus_old=abs(D.q_pred_hybrid_T18-D.q_pred_hybrid_old);
D.abs_sws_hybrid_minus_local=abs(D.sws_pred_hybrid_T18-D.sws_pred_local_T18);
D.abs_sws_hybrid_minus_global=abs(D.sws_pred_hybrid_T18-D.sws_pred_global_T18);
D.abs_sws_hybrid_minus_theory=abs(D.sws_pred_hybrid_T18-D.sws_theory_discrete);
D.abs_sws_hybrid_T18_minus_old=abs(D.sws_pred_hybrid_T18-D.sws_pred_hybrid_old);

D=add_map_features(D,neighborhood);
[spectral,mixed,INFO]=load_spectral_and_mixed_features(D,root_dir);
for i=1:numel(INFO.spectral_predictors)
    name=INFO.spectral_predictors(i);
    if ~ismember(name,string(D.Properties.VariableNames))
        D.(name)=spectral.(name);
    end
end
D.mixed_window_fraction=mixed;
end

function D=add_map_features(D,nhood)
D.grad_q_hybrid=nan(height(D),1); D.grad_sws_hybrid=nan(height(D),1);
D.local_std_q_hybrid=nan(height(D),1); D.local_std_sws_hybrid=nan(height(D),1);
D.local_range_q_hybrid=nan(height(D),1); D.local_range_sws_hybrid=nan(height(D),1);
conditions=unique(D.condition_id,'stable');
for ci=1:numel(conditions)
    idx=find(D.condition_id==conditions(ci)); Ti=D(idx,:);
    Q=map_from_rows(Ti,Ti.q_pred_hybrid_T18);
    C=map_from_rows(Ti,Ti.sws_pred_hybrid_T18);
    [dqz,dqx]=gradient(Q,Ti.dz(1),Ti.dx(1));
    [dcz,dcx]=gradient(C,Ti.dz(1),Ti.dx(1));
    [qstd,qrange]=local_stats(Q,nhood);
    [cstd,crange]=local_stats(C,nhood);
    lin=sub2ind(size(Q),Ti.map_iz,Ti.map_ix);
    qgrad=hypot(dqx,dqz); cgrad=hypot(dcx,dcz);
    D.grad_q_hybrid(idx)=qgrad(lin); D.grad_sws_hybrid(idx)=cgrad(lin);
    D.local_std_q_hybrid(idx)=qstd(lin); D.local_std_sws_hybrid(idx)=cstd(lin);
    D.local_range_q_hybrid(idx)=qrange(lin); D.local_range_sws_hybrid(idx)=crange(lin);
end
end

function [S,R]=local_stats(A,nhood)
r=floor(nhood/2); K=ones(nhood);
valid=isfinite(A); A0=A; A0(~valid)=0;
N=conv2(double(valid),K,'same');
mu=conv2(A0,K,'same')./max(N,1);
mu2=conv2(A0.^2,K,'same')./max(N,1);
S=sqrt(max(mu2-mu.^2,0));
Amax=movmax(movmax(A,[r r],1,'omitnan'),[r r],2,'omitnan');
Amin=movmin(movmin(A,[r r],1,'omitnan'),[r r],2,'omitnan');
R=Amax-Amin; S(~valid)=NaN; R(~valid)=NaN;
end

function A=map_from_rows(T,v)
A=nan(max(T.map_iz),max(T.map_ix));
A(sub2ind(size(A),T.map_iz,T.map_ix))=v;
end

function [F,mixed,INFO]=load_spectral_and_mixed_features(D,root_dir)
model_dir=fullfile(root_dir,'outputs','model_registry','test19_clean_field_regime');
[MODEL,~,~]=adaptive_req.analysis.load_q_model_deployment(model_dir, ...
    'ModelName','HybridLocalGlobal_T18_noUserRegime', ...
    'FeatureSet','CleanFieldRegime_noUser','ModelType','bagged_trees');
predictors=string({MODEL.encoder.entries.name})';
assert_no_leakage(predictors);
INFO=struct('spectral_predictors',predictors, ...
    'base_model','HybridLocalGlobal_T18_noUserRegime');
F=table();
for i=1:numel(predictors), F.(predictors(i))=nan(height(D),1); end
mixed=nan(height(D),1);

conditions=unique(D.condition_id,'stable');
base_dir=fullfile(root_dir,'outputs', ...
    'test_17_model_comparison_heterogeneous_cases','analysis','data','conditions');
for ci=1:numel(conditions)
    idx=find(D.condition_id==conditions(ci)); Ti=D(idx,:);
    token=regime_file_token(Ti.field_regime(1));
    file=fullfile(base_dir,sprintf('level17_compare_%s_%s_M%g.mat', ...
        case_file_token(Ti.external_case(1)),token,Ti.REQ_M(1)));
    assert(exist(file,'file')==2,'Test 17 condition missing: %s',file);
    S=load(file,'T_feat','sim');
    Tf=sortrows(S.T_feat,{'map_iz','map_ix'});
    [~,ord]=sortrows([Ti.map_iz Ti.map_ix],[1 2]);
    assert(height(Tf)==numel(idx) && ...
        isequal([Tf.map_iz Tf.map_ix],[Ti.map_iz(ord) Ti.map_ix(ord)]), ...
        'Feature alignment failed for condition %d.',conditions(ci));
    target_idx=idx(ord);
    for pi=1:numel(predictors)
        name=predictors(pi);
        assert(ismember(name,string(Tf.Properties.VariableNames)), ...
            'Required spectral predictor %s missing from %s.',name,file);
        F.(name)(target_idx)=double(Tf.(name));
    end
    mixed(target_idx)=mixed_fraction_from_true_map(Tf,S.sim.cs_map, ...
        Ti.REQ_M(1),Ti.cs_guess(1),Ti.SIM_f0(1),Ti.dx(1));
end
end

function token=regime_file_token(regime)
switch string(regime)
    case "directional_2D", token='directional_2d';
    case "diffuse_2D", token='diffuse_2d';
    case "partial_3D", token='partial_diffuse_3d';
    case "diffuse_3D", token='diffuse_3d';
    otherwise, error('Unknown field regime: %s',regime);
end
end

function token=case_file_token(case_name)
token=char(string(case_name));
if string(case_name)=="circular_inclusion_2_3", token='inclusion_2_3'; end
end

function frac=mixed_fraction_from_true_map(T,cs_map,M,cs_guess,f0,dx)
win=round(M*(cs_guess/f0)/dx); if mod(win,2)==0, win=win+1; end
h=floor(win/2); frac=zeros(height(T),1);
if max(cs_map(:))-min(cs_map(:))<1e-9, return; end
threshold=0.5*(min(cs_map(:))+max(cs_map(:)));
for i=1:height(T)
    zi=(T.cz(i)-h):(T.cz(i)+h); xi=(T.cx(i)-h):(T.cx(i)+h);
    high=cs_map(zi,xi)>threshold; p=mean(high(:)); frac(i)=min(p,1-p);
end
end

function assert_no_leakage(predictors)
p=lower(string(predictors(:)));
banned=lower(["q_theory";"q_true";"q_global_theory";"q_local_minus_global"; ...
    "abs_q_local_minus_global";"m_eff_true_diag";"cs_true";"true_cs"; ...
    "cs_pred_error";"sws_error";"abs_sws_error";"sws_error_pct"; ...
    "abs_sws_error_pct";"high_error_gt20";"high_error_gt10"; ...
    "distance_to_interface_mm";"near_interface";"mixed_window_fraction"; ...
    "aperture_weight";"solid_angle_weight";"true_aperture_weight"]);
bad=intersect(p,banned);
assert(isempty(bad),'Leakage predictors detected: %s',strjoin(bad,', '));
end

function [train_mask,test_mask,key]=grouped_condition_split(T,fraction,seed)
key=T.external_case+"__"+T.field_regime+"__M"+string(T.REQ_M)+ ...
    "__f"+string(T.SIM_f0)+"__dx"+string(T.dx)+ ...
    "__condition"+string(T.condition_id);
groups=unique(key,'stable'); ntest=max(1,round(fraction*numel(groups)));
for attempt=0:99
    rng(seed+attempt,'twister'); order=randperm(numel(groups));
    test_groups=groups(order(1:ntest));
    test_mask=ismember(key,test_groups); train_mask=~test_mask;
    if any(T.high_error_gt20(train_mask)) && any(~T.high_error_gt20(train_mask)) && ...
            any(T.high_error_gt20(test_mask)) && any(~T.high_error_gt20(test_mask))
        return;
    end
end
error('Unable to create a grouped split containing both classes.');
end

function [X,med,mu,sigma]=preprocess_from_training(X,train_mask)
med=median(X(train_mask,:),1,'omitnan'); med(~isfinite(med))=0;
for j=1:size(X,2), bad=~isfinite(X(:,j)); X(bad,j)=med(j); end
mu=mean(X(train_mask,:),1); sigma=std(X(train_mask,:),0,1);
sigma(~isfinite(sigma)|sigma<eps)=1;
X=(X-mu)./sigma;
end

function w=balanced_class_weights(y)
n=numel(y); np=sum(y); nn=sum(~y);
assert(np>0 && nn>0,'Both classes are required for weighted training.');
w=zeros(n,1); w(y)=n/(2*np); w(~y)=n/(2*nn);
end

function D=make_detector(name,family,type,model,predictors,med,mu,sigma,rule)
D=struct('detector_name',string(name),'detector_family',string(family), ...
    'model_type',string(type),'model',model,'predictors',string(predictors(:)), ...
    'impute_median',double(med(:))','predictor_mu',double(mu(:))', ...
    'predictor_sigma',double(sigma(:))','rule_feature',string(rule));
end

function p=predict_detector(D,T,X,mask)
if islogical(mask), rows=find(mask); else, rows=mask; end
if D.detector_family=="ml"
    Xin=X(rows,:);
else
    raw=double(T.(D.rule_feature)(rows)); raw(~isfinite(raw))=D.impute_median;
    Xin=raw;
end
[~,scores]=predict(D.model,Xin); p=positive_score(D.model,scores);
end

function p=positive_score(model,scores)
classes=string(model.ClassNames); idx=find(classes=="1"|classes=="true",1);
if isempty(idx), idx=size(scores,2); end
p=scores(:,idx); p=min(max(double(p),0),1);
end

function P=prediction_table(T,p,D)
vars=["condition_id","external_case","field_regime","regime_variant", ...
    "REQ_M","SIM_f0","dx","dz","x","z","map_iz","map_ix", ...
    "true_cs","abs_sws_error_pct","sws_error_pct","high_error_gt20", ...
    "high_error_gt10","distance_to_interface_mm","near_interface", ...
    "mixed_window_fraction"];
P=T(:,cellstr(vars));
P.detector_name=repmat(D.detector_name,height(P),1);
P.detector_family=repmat(D.detector_family,height(P),1);
P.model_type=repmat(D.model_type,height(P),1);
P.risk_probability=p; P.confidence=1-p;
P.probability_high_error=p;
P.predicted_high_error_020=p>=0.20;
P=movevars(P,{'detector_name','detector_family','model_type'},'Before',1);
end

function M=evaluate_detector(y,p,threshold)
y=logical(y); valid=isfinite(p); y=y(valid); p=p(valid);
if any(y) && any(~y)
    [~,~,~,~,auc_roc,auc_pr]=binary_curves(y,p);
else
    auc_roc=NaN; auc_pr=NaN;
end
pred=p>=threshold; tp=sum(pred&y); fp=sum(pred&~y);
tn=sum(~pred&~y); fn=sum(~pred&y);
M=table(); M.Threshold=threshold; M.N=numel(y);
M.N_positive=sum(y); M.PositiveFraction=mean(y);
M.AUC_ROC=auc_roc; M.AUC_PR=auc_pr;
M.Precision=safe_div(tp,tp+fp); M.Recall=safe_div(tp,tp+fn);
M.F1=safe_div(2*M.Precision*M.Recall,M.Precision+M.Recall);
M.Specificity=safe_div(tn,tn+fp); M.FalseNegativeRate=safe_div(fn,fn+tp);
M.FalsePositiveRate=safe_div(fp,fp+tn);
M.TP=tp; M.FP=fp; M.TN=tn; M.FN=fn;
end

function [fpr,tpr,recall,precision,auc_roc,auc_pr]=binary_curves(y,p)
[~,ord]=sort(p,'descend'); y=double(y(ord));
tp=cumsum(y); fp=cumsum(1-y); P=sum(y); N=sum(1-y);
tpr=[0;tp/max(P,1);1]; fpr=[0;fp/max(N,1);1];
recall=[0;tp/max(P,1)]; precision=[1;tp./max(tp+fp,1)];
auc_roc=trapz(fpr,tpr); auc_pr=trapz(recall,precision);
end

function v=safe_div(a,b)
if b==0, v=NaN; else, v=a/b; end
end

function Tout=evaluate_by_group(T,var,threshold)
Tout=table(); detectors=unique(T.detector_name,'stable'); values=unique(T.(var),'stable');
for di=1:numel(detectors)
    for vi=1:numel(values)
        idx=T.detector_name==detectors(di) & T.(var)==values(vi);
        if ~any(idx), continue; end
        M=evaluate_detector(T.high_error_gt20(idx),T.risk_probability(idx),threshold);
        M.detector_name=detectors(di); M.detector_family=T.detector_family(find(idx,1));
        M.(var)=values(vi); Tout=concat_tables(Tout,M);
    end
end
Tout=movevars(Tout,{'detector_name','detector_family',var},'Before',1);
end

function T=detector_importance(D)
if D.detector_family~="ml"
    imp=1; names=D.rule_feature;
elseif D.model_type=="logistic_regression"
    imp=abs(double(D.model.Beta(:))); names=D.predictors;
else
    imp=double(predictorImportance(D.model))'; names=D.predictors;
end
imp=imp(:); names=names(:); s=sum(imp,'omitnan'); if s>0, imp=imp/s; end
T=table(repmat(D.detector_name,numel(names),1), ...
    repmat(D.detector_family,numel(names),1),names,imp, ...
    'VariableNames',{'detector_name','detector_family','feature_name','importance'});
T=sortrows(T,'importance','descend'); T.rank=(1:height(T))';
end

function threshold=recommend_threshold(T)
eligible=T(T.Recall>=0.90,:);
if ~isempty(eligible)
    eligible=sortrows(eligible,{'Precision','Threshold'},{'descend','descend'});
    threshold=eligible.Threshold(1);
else
    T=sortrows(T,{'Recall','Precision'},{'descend','descend'});
    threshold=T.Threshold(1);
end
end

function plot_roc_curves(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 12]); ax=axes(fig); hold(ax,'on');
names=unique(T.detector_name,'stable');
for i=1:numel(names)
    Ti=T(T.detector_name==names(i),:); [fpr,tpr,~,~,auc,~]= ...
        binary_curves(Ti.high_error_gt20,Ti.risk_probability);
    plot(ax,fpr,tpr,'LineWidth',1.3,'DisplayName',short_detector(names(i))+ ...
        sprintf(' (%.3f)',auc));
end
plot(ax,[0 1],[0 1],'k:'); xlabel(ax,'False positive rate'); ylabel(ax,'Recall');
title(ax,'Held-out grouped ROC curves','FontWeight','normal'); axis(ax,'square'); grid(ax,'on');
legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'level21_roc_curves.png'));
end

function plot_pr_curves(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 12]); ax=axes(fig); hold(ax,'on');
names=unique(T.detector_name,'stable'); prevalence=mean(T.high_error_gt20(T.detector_name==names(1)));
for i=1:numel(names)
    Ti=T(T.detector_name==names(i),:); [~,~,r,p,~,auc]= ...
        binary_curves(Ti.high_error_gt20,Ti.risk_probability);
    plot(ax,r,p,'LineWidth',1.3,'DisplayName',short_detector(names(i))+ ...
        sprintf(' (%.3f)',auc));
end
yline(ax,prevalence,'k:','Prevalence'); xlabel(ax,'Recall'); ylabel(ax,'Precision');
title(ax,'Held-out grouped precision-recall curves','FontWeight','normal');
xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,'on');
legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'level21_precision_recall_curves.png'));
end

function plot_confidence_calibration(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 19 12]); ax=axes(fig); hold(ax,'on');
names=unique(T.detector_name,'stable'); edges=0:0.1:1;
for i=1:numel(names)
    Ti=T(T.detector_name==names(i),:); b=discretize(Ti.confidence,edges);
    x=nan(10,1); y=nan(10,1);
    for j=1:10, idx=b==j; x(j)=mean(Ti.confidence(idx),'omitnan'); y(j)=100*mean(Ti.high_error_gt20(idx),'omitnan'); end
    plot(ax,x,y,'-o','MarkerSize',4,'DisplayName',short_detector(names(i)));
end
xlabel(ax,'Confidence (1 - predicted risk)'); ylabel(ax,'Observed high-error rate (%)');
title(ax,'High-error rate by confidence bin','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'level21_high_error_rate_by_confidence_bin.png'));
end

function plot_error_vs_confidence(T,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 19 12]); ax=axes(fig); hold(ax,'on');
names=unique(T.detector_name,'stable'); edges=0:0.1:1;
for i=1:numel(names)
    Ti=T(T.detector_name==names(i),:); b=discretize(Ti.confidence,edges);
    x=nan(10,1); y=nan(10,1);
    for j=1:10, idx=b==j; x(j)=mean(Ti.confidence(idx),'omitnan'); y(j)=mean(Ti.abs_sws_error_pct(idx),'omitnan'); end
    plot(ax,x,y,'-o','MarkerSize',4,'DisplayName',short_detector(names(i)));
end
xlabel(ax,'Confidence (1 - predicted risk)'); ylabel(ax,'Mean absolute SWS error (%)');
title(ax,'SWS error versus confidence','FontWeight','normal'); grid(ax,'on');
legend(ax,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'level21_error_vs_confidence.png'));
end

function plot_feature_importance(T,best,OUT)
Ti=T(T.detector_name==best,:); Ti=sortrows(Ti,'importance','descend'); Ti=Ti(1:min(20,height(Ti)),:);
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 13]); ax=axes(fig);
barh(ax,flip(Ti.importance)); yticks(ax,1:height(Ti));
yticklabels(ax,flip(strrep(Ti.feature_name,'_',' '))); ax.TickLabelInterpreter='none';
xlabel(ax,'Normalized importance'); title(ax,'Top confidence-detector features','FontWeight','normal');
grid(ax,'on'); export_fig(fig,fullfile(OUT.figure_dir,'level21_feature_importance.png'));
end

function plot_rule_ml_comparison(T,OUT)
T=sortrows(T,'AUC_PR','descend');
fig=figure('Color','w','Units','centimeters','Position',[2 2 20 12]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
vars=["AUC_PR","Recall"]; labels=["AUC PR","Recall at risk threshold 0.20"];
for j=1:2
    ax=nexttile(tl); barh(ax,flip(T.(vars(j)))); yticks(ax,1:height(T));
    yticklabels(ax,flip(short_detector(T.detector_name))); ax.TickLabelInterpreter='none';
    xlabel(ax,labels(j)); xlim(ax,[0 1]); grid(ax,'on');
end
title(tl,'Rule-based versus ML confidence detectors','FontWeight','normal');
export_fig(fig,fullfile(OUT.figure_dir,'level21_rule_vs_ml_detector_comparison.png'));
end

function plot_representative_confidence_maps(T,best,OUT)
cases=["homogeneous_cs2","homogeneous_cs3","bilayer_2_3","circular_inclusion_2_3"];
regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 30 19]);
tl=tiledlayout(4,6,'TileSpacing','tight','Padding','compact');
for i=1:4
    Ti=select_map(T,cases(i),regimes(i),2);
    values={Ti.true_cs,Ti.sws_pred_hybrid_T18,Ti.abs_sws_error_pct, ...
        Ti.best_detector_confidence,Ti.distance_to_interface_mm,Ti.grad_q_hybrid};
    names={'true SWS','predicted SWS','absolute error (%)', ...
        'confidence','distance (mm)','|grad q|'};
    for j=1:6
        ax=nexttile(tl); A=map_from_rows(Ti,values{j}); imagesc(ax,A); axis(ax,'image','off');
        colorbar(ax,'FontSize',5); if i==1, title(ax,names{j},'FontSize',7,'FontWeight','normal'); end
        if j==1
            text(ax,-0.05,0.5,short_case(cases(i))+" / "+short_regime(regimes(i)), ...
                'Units','normalized','HorizontalAlignment','right','FontSize',6, ...
                'Interpreter','none','Clipping','off');
        end
    end
end
title(tl,"Representative confidence maps: "+short_detector(best), ...
    'FontSize',9,'FontWeight','normal');
export_fig(fig,fullfile(OUT.figure_dir,'level21_confidence_maps_representative.png'));
end

function plot_interface_confidence_maps(T,best,OUT)
cases=["bilayer_2_3","circular_inclusion_2_3"];
regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 24 27]);
tl=tiledlayout(8,4,'TileSpacing','tight','Padding','compact'); row=0;
for ci=1:numel(cases)
    for ri=1:numel(regimes)
        row=row+1; Ti=select_map(T,cases(ci),regimes(ri),2);
        values={Ti.abs_sws_error_pct,Ti.best_detector_confidence, ...
            Ti.distance_to_interface_mm,Ti.grad_q_hybrid};
        names={'absolute error (%)','confidence','distance (mm)','|grad q|'};
        for j=1:4
            ax=nexttile(tl); imagesc(ax,map_from_rows(Ti,values{j})); axis(ax,'image','off');
            colorbar(ax,'FontSize',5); if row==1, title(ax,names{j},'FontSize',7,'FontWeight','normal'); end
            if j==1
                text(ax,-0.05,0.5,short_case(cases(ci))+" / "+short_regime(regimes(ri)), ...
                    'Units','normalized','HorizontalAlignment','right','FontSize',5.5, ...
                    'Interpreter','none','Clipping','off');
            end
        end
    end
end
title(tl,"Interface error and confidence: "+short_detector(best), ...
    'FontSize',9,'FontWeight','normal');
export_fig(fig,fullfile(OUT.figure_dir,'level21_error_confidence_interface_maps.png'));
end

function Ti=select_map(T,case_name,regime,M)
idx=T.external_case==case_name & T.field_regime==regime & T.REQ_M==M;
assert(any(idx),'Representative map unavailable: %s / %s / M=%g',case_name,regime,M);
Ti=T(idx,:);
end

function print_summary(T,P,M,Tsweep,Tcase,Tregime,TM,best)
Mb=M(M.detector_name==best,:); threshold=recommend_threshold(Tsweep(Tsweep.detector_name==best,:));
qgrad=M(M.detector_name=="Rule_q_gradient",:);
dis=M(M.detector_name=="Rule_max_local_theory_disagreement",:);
ranked=sortrows(M,{'Recall','AUC_PR'},{'descend','descend'});
best_recall=ranked(1,:);
Pb=P(P.detector_name==best,:); fn=Pb.high_error_gt20 & Pb.risk_probability<0.20;
near=mean(Pb.confidence(Pb.near_interface),'omitnan');
hom=mean(Pb.confidence(contains(Pb.external_case,"homogeneous")),'omitnan');
prevalence=mean(T.high_error_gt20);
fprintf('\n================ Test 21 summary ================\n');
fprintf('Best deployable ML detector: %s\n',best);
fprintf('AUC PR: %.4f\n',Mb.AUC_PR);
fprintf('AUC ROC: %.4f\n',Mb.AUC_ROC);
fprintf('Recall at threshold 0.20: %.4f\n',Mb.Recall);
fprintf('Precision at threshold 0.20: %.4f\n',Mb.Precision);
fprintf('False negative rate: %.4f\n',Mb.FalseNegativeRate);
fprintf('High-error prevalence: %.4f%%\n',100*prevalence);
fprintf('Recommended high-error risk threshold: %.2f\n',threshold);
fprintf('Equivalent minimum confidence: %.2f\n',1-threshold);
fprintf('q-gradient AUC PR: %.4f | max disagreement AUC PR: %.4f\n', ...
    qgrad.AUC_PR,dis.AUC_PR);
fprintf('Mean confidence near interface: %.4f | homogeneous: %.4f\n',near,hom);
fprintf('Highest recall at threshold 0.20: %s (%.4f).\n', ...
    best_recall.detector_name,best_recall.Recall);
if Mb.AUC_PR>prevalence
    fprintf('Operational detection signal: YES; AUC PR exceeds prevalence baseline.\n');
else
    fprintf('Operational detection signal: WEAK; AUC PR does not exceed prevalence baseline.\n');
end
if dis.AUC_PR>qgrad.AUC_PR
    fprintf('Model disagreement adds signal beyond q-gradient alone.\n');
else
    fprintf('q-gradient matches or exceeds the tested disagreement rule.\n');
end
if near<hom
    fprintf('Confidence decreases near interfaces without using interface distance as input.\n');
else
    fprintf('Confidence did not decrease near interfaces in the held-out set.\n');
end

if any(fn)
    F=Pb(fn,:);
    [mode_case,~]=mode(categorical(F.external_case));
    [mode_regime,~]=mode(categorical(F.field_regime));
    [mode_M,~]=mode(F.REQ_M);
    fprintf('Main false-negative mode: %s | %s | M=%g\n', ...
        string(mode_case),string(mode_regime),mode_M);
else
    fprintf('Main false-negative mode: none at threshold 0.20.\n');
end
fprintf(['Operational detector uses no interface distance, true SWS, mixed-window ', ...
    'fraction, or error label as an input.\n']);
fprintf('Detailed failure tables: %d case rows, %d regime rows, %d M rows.\n', ...
    height(Tcase),height(Tregime),height(TM));
fprintf('=================================================\n');
end

function label=short_detector(name)
name=string(name); label=name;
label=replace(label,"Rule_q_gradient","q-gradient rule");
label=replace(label,"Rule_q_hybrid_minus_local","hybrid-local rule");
label=replace(label,"Rule_q_hybrid_minus_theory","hybrid-theory rule");
label=replace(label,"Rule_max_local_theory_disagreement","max-disagreement rule");
label=replace(label,"ML_logistic_regression","logistic ML");
label=replace(label,"ML_bagged_trees","bagged trees ML");
label=replace(label,"ML_boosted_trees","boosted trees ML");
end

function label=short_case(x)
label=replace(string(x),["homogeneous_cs2","homogeneous_cs3", ...
    "bilayer_2_3","circular_inclusion_2_3"], ...
    ["homogeneous 2","homogeneous 3","bilayer 2/3","inclusion 2/3"]);
end

function label=short_regime(x)
label=replace(string(x),["directional_2D","diffuse_2D","partial_3D","diffuse_3D"], ...
    ["directional 2D","diffuse 2D","partial 3D","diffuse 3D"]);
end

function export_fig(fig,file)
exportgraphics(fig,file,'Resolution',220); close(fig);
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
if isempty(B), T=A; return; end
vars=unique([string(A.Properties.VariableNames),string(B.Properties.VariableNames)],'stable');
A=add_missing(A,vars); B=add_missing(B,vars); T=[A(:,cellstr(vars));B(:,cellstr(vars))];
end

function T=add_missing(T,vars)
for i=1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)), continue; end
    T.(vars(i))=nan(height(T),1);
end
end
