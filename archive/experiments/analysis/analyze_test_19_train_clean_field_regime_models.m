%% analyze_test_19_train_clean_field_regime_models.m
% Test 19: train adaptive-q models on the clean Test 18 field regimes.
%
% Primary models infer field regime from operational spectral features.
% One additional Hybrid model accepts a categorical user field-regime guess.
% Old Test 12 deployments and discrete theory-q are evaluation baselines.

clear; clc; close all;
format compact;

set(groot,'defaultAxesFontSize',11);
set(groot,'defaultTextFontSize',11);
set(groot,'defaultLegendFontSize',10);

%% Runtime profile

PROFILE = lower(string(getenv('ADAPTIVE_REQ_TEST19_PROFILE')));
if PROFILE == "", PROFILE = "full"; end
switch PROFILE
    case "full"
        DATASET_FOLDER = "dataset";
        NUM_TREES_PRIMARY = 150;
        NUM_TREES_GENERALIZATION = 100;
        MIN_LEAF_SIZE = 8;
    case "fast"
        DATASET_FOLDER = "dataset";
        NUM_TREES_PRIMARY = 80;
        NUM_TREES_GENERALIZATION = 40;
        MIN_LEAF_SIZE = 10;
    case "smoke"
        DATASET_FOLDER = "smoke_test";
        NUM_TREES_PRIMARY = 10;
        NUM_TREES_GENERALIZATION = 5;
        MIN_LEAF_SIZE = 2;
    otherwise
        error('Unknown Test 19 profile: %s',PROFILE);
end
USE_PARALLEL = true;
TRAIN_FRACTION = 0.75;
RANDOM_SEED = 19001;

%% Project and data

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

dataset_dir = fullfile(root_dir,'outputs', ...
    'test_18_clean_field_regime_training_dataset',DATASET_FOLDER);
dataset_file = fullfile(dataset_dir,'data', ...
    'test18_clean_field_regime_dataset.mat');
assert(exist(dataset_file,'file')==2, ...
    ['Test 18 dataset not found:\n%s\nRun ', ...
    'experiments/run_test_18_clean_field_regime_training_dataset.m first.'], ...
    dataset_file);
S18 = load(dataset_file,'T_dataset','MC','CFG','PREDICTOR_POLICY');
T = S18.T_dataset;
MC18 = S18.MC;

OUT = make_output_dirs(root_dir,PROFILE);
fprintf('\nTest 19 profile: %s\nDataset: %s\nRows: %d\n', ...
    PROFILE,dataset_file,height(T));

%% Prepare table and enforce leakage policy

required = ["q_theory","req_mapping","condition_id","realization_idx", ...
    "patch_idx","REQ_M","SIM_f0","SIM_cs_bg","SIM_dx", ...
    "field_regime_label","field_regime_variant"];
require_vars(T,required,'Test 19 input');

T.step_idx = T.REQ_M;
T.row_key = make_test18_row_key(T);
T.user_field_guess = map_user_field_guess(T.field_regime_label);
if ~ismember('global_REQ_Nbins_effective',T.Properties.VariableNames)
    T.global_REQ_Nbins_effective = cellfun(@mapping_nbins,T.global_req_mapping);
end
T = adaptive_req.analysis.Test12Analysis.prepareFeatureTable(T);
T = adaptive_req.analysis.Test12Analysis.addBins(T);

SPECS = build_test19_specs(T);
for i=1:numel(SPECS)
    assert_operational_predictors(SPECS(i).predictors,SPECS(i).model_name, ...
        SPECS(i).allow_user_guess);
end

%% Strict grouped primary split

% condition_id contains regime variant, f0, true SWS, and dx=dz. Holding out
% complete condition IDs is stricter than splitting patches or realizations.
[train_mask,test_mask] = grouped_condition_split(T,TRAIN_FRACTION,RANDOM_SEED);
assert(~any(train_mask & test_mask) && any(train_mask) && any(test_mask));
assert(isempty(intersect(unique(T.condition_id(train_mask)), ...
    unique(T.condition_id(test_mask)))), ...
    'Physical condition leakage detected in primary split.');

%% Train primary Test 18 models

MODELS = struct([]);
T_new_all = table();
T_q_metrics_new = table();
if PROFILE=="full"
    registry_dir = fullfile(root_dir,'outputs','model_registry', ...
        'test19_clean_field_regime');
else
    registry_dir = fullfile(OUT.model_dir,'smoke_deployment_validation');
end
if exist(registry_dir,'dir')~=7, mkdir(registry_dir); end

for i=1:numel(SPECS)
    spec = SPECS(i);
    primary_file=fullfile(OUT.model_dir, ...
        "primary_checkpoint__"+sanitize_filename(spec.model_name)+".mat");
    if exist(primary_file,'file')==2
        fprintf('\nReusing primary model checkpoint: %s\n',spec.model_name);
        Sp=load(primary_file,'MODEL_i','T_pred_i','T_q_i');
        MODEL_i=Sp.MODEL_i; T_pred_i=Sp.T_pred_i; T_q_i=Sp.T_q_i;
    else
        fprintf('\n=== Primary training: %s ===\n',spec.model_name);
        [MODEL_i,T_pred_i,T_q_i] = adaptive_req.analysis.train_q_model_fixed_split( ...
            T,spec.predictors,train_mask,test_mask, ...
            'QVar','q_theory', ...
            'ModelName',spec.model_name, ...
            'ModelRole','operational', ...
            'ModelTypes',"bagged_trees", ...
            'NumLearningCycles',NUM_TREES_PRIMARY, ...
            'MinLeafSize',MIN_LEAF_SIZE, ...
            'UseParallel',USE_PARALLEL, ...
            'Verbose',true);

        T_pred_i = enrich_predictions(T_pred_i,T,spec);
        T_pred_i = compact_prediction_table(T_pred_i);
        T_q_i.feature_set = repmat(spec.feature_set,height(T_q_i),1);
        T_q_i.model_role = repmat(spec.model_role,height(T_q_i),1);
        save(primary_file,'MODEL_i','T_pred_i','T_q_i','spec','-v7.3');
    end

    local_deploy = adaptive_req.analysis.save_q_model_deployment( ...
        MODEL_i,OUT.model_dir,'ModelName',spec.model_name, ...
        'FeatureSet',spec.feature_set,'ModelRole',spec.model_role, ...
        'ModelTypes',"bagged_trees",'Overwrite',true);
    registry_deploy = adaptive_req.analysis.save_q_model_deployment( ...
        MODEL_i,registry_dir,'ModelName',spec.model_name, ...
        'FeatureSet',spec.feature_set,'ModelRole',spec.model_role, ...
        'ModelTypes',"bagged_trees",'Overwrite',true);

    MODELS(i).model_name = spec.model_name;
    MODELS(i).feature_set = spec.feature_set;
    MODELS(i).model_role = spec.model_role;
    MODELS(i).predictors = spec.predictors;
    MODELS(i).model = MODEL_i;
    MODELS(i).local_deployment = local_deploy;
    MODELS(i).registry_deployment = registry_deploy;

    T_new_all = concat_tables(T_new_all,T_pred_i);
    T_q_metrics_new = concat_tables(T_q_metrics_new,T_q_i);
end

%% Existing model and theory baselines on the same rows

OLD_SPECS = load_old_models(root_dir);
T_old_all = table();
for i=1:numel(OLD_SPECS)
    assert_operational_predictors( ...
        string({OLD_SPECS(i).model.encoder.entries.name}), ...
        OLD_SPECS(i).model_name,false);
    Tq = adaptive_req.analysis.predict_q_model_from_table( ...
        OLD_SPECS(i).model,T,'ModelType','bagged_trees', ...
        'ModelName',OLD_SPECS(i).model_name);
    spec = OLD_SPECS(i);
    Tq.split = repmat("external",height(Tq),1);
    Tq = enrich_predictions(Tq,T,spec);
    Tq = compact_prediction_table(Tq);
    T_old_all = concat_tables(T_old_all,Tq);
end

q_discrete = compute_discrete_theory_q(T);
theory_spec = struct('model_name',"TheoryQDiscrete", ...
    'feature_set',"field_matched_discrete_theory", ...
    'model_role',"diagnostic_only",'allow_user_guess',false);
T_theory_all = baseline_prediction_table(T,q_discrete,theory_spec,"theory_no_ml");

%% Primary test predictions and metrics

test_keys = T.row_key(test_mask);
T_new_test = T_new_all(T_new_all.split=="test",:);
T_old_test = T_old_all(ismember(T_old_all.row_key,test_keys),:);
T_theory_test = T_theory_all(ismember(T_theory_all.row_key,test_keys),:);
T_predictions = concat_tables(concat_tables(T_new_test,T_old_test),T_theory_test);

T_q_metrics = summarize_q(T_predictions, ...
    ["model_name","feature_set","model_role","model_type"]);
T_sws_metrics = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type"]);
T_by_regime = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","field_regime_label"]);
T_by_variant = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","field_regime_variant"]);
T_by_dx = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","SIM_dx"]);
T_by_M = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","REQ_M"]);
T_by_cs = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","SIM_cs_bg"]);
T_by_frequency = summarize_sws(T_predictions, ...
    ["model_name","feature_set","model_role","model_type","SIM_f0"]);
T_by_Meff = table();
if ismember('M_eff_guess_bin',T_predictions.Properties.VariableNames)
    T_by_Meff = summarize_sws(T_predictions, ...
        ["model_name","feature_set","model_role","model_type","M_eff_guess_bin"]);
end

%% Leave-one-group-out generalization

split_specs = [ ...
    struct('name',"leave_one_field_regime_out",'var',"field_regime_label")
    struct('name',"leave_one_dx_out",'var',"SIM_dx")
    struct('name',"leave_one_frequency_out",'var',"SIM_f0")
    struct('name',"leave_one_cs_out",'var',"SIM_cs_bg")
    struct('name',"leave_one_M_out",'var',"REQ_M")];

T_generalization = table();
for si=1:numel(split_specs)
    split_name = split_specs(si).name;
    heldout_var = split_specs(si).var;
    values = unique(T.(char(heldout_var)),'stable');
    if numel(values)<2
        warning('Skipping %s: only one value in %s.',split_name,heldout_var);
        continue;
    end

    for vi=1:numel(values)
        heldout_value = values(vi);
        [train_i,test_i] = make_generalization_split( ...
            T,heldout_var,heldout_value,RANDOM_SEED+100*si+vi);
        assert_grouped_split(T,train_i,test_i,heldout_var,heldout_value);

        % Old deployments and theory require no retraining.
        baseline_all = concat_tables(T_old_all,T_theory_all);
        heldout_test_keys=T.row_key(test_i);
        B = baseline_all(ismember(baseline_all.row_key,heldout_test_keys),:);
        T_generalization = concat_tables(T_generalization, ...
            generalization_summary(B,split_name,heldout_var,heldout_value, ...
            NaN,sum(test_i)));

        for mi=1:numel(SPECS)
            spec = SPECS(mi);
            job_key = sanitize_filename("splitv2__"+split_name+"__"+heldout_var+"__"+ ...
                value_string(heldout_value)+"__"+spec.model_name);
            job_file = fullfile(OUT.job_dir,job_key+".mat");

            if exist(job_file,'file')==2
                fprintf('Reusing generalization job: %s\n',job_key);
                Sj = load(job_file,'T_job');
                T_job = Sj.T_job;
            else
                fprintf('\n=== %s | %s=%s | %s ===\n', ...
                    split_name,heldout_var,value_string(heldout_value),spec.model_name);
                [~,Tj,~] = adaptive_req.analysis.train_q_model_fixed_split( ...
                    T,spec.predictors,train_i,test_i, ...
                    'QVar','q_theory','ModelName',spec.model_name, ...
                    'ModelRole','operational','ModelTypes',"bagged_trees", ...
                    'NumLearningCycles',NUM_TREES_GENERALIZATION, ...
                    'MinLeafSize',MIN_LEAF_SIZE,'UseParallel',USE_PARALLEL, ...
                    'Verbose',false);
                Tj = enrich_predictions(Tj,T,spec);
                T_job = compact_prediction_table(Tj(Tj.split=="test",:));
                save(job_file,'T_job','spec','split_name','heldout_var', ...
                    'heldout_value','-v7.3');
            end

            T_generalization = concat_tables(T_generalization, ...
                generalization_summary(T_job,split_name,heldout_var, ...
                heldout_value,sum(train_i),sum(test_i)));
        end
    end
end

%% Save tables and primary model bundle

writetable(remove_cell_columns(T_predictions),fullfile(OUT.table_dir, ...
    'level19_predictions.csv'));
writetable(T_q_metrics,fullfile(OUT.table_dir,'level19_q_metrics.csv'));
writetable(T_sws_metrics,fullfile(OUT.table_dir,'level19_sws_metrics.csv'));
writetable(T_by_regime,fullfile(OUT.table_dir,'level19_metrics_by_regime.csv'));
writetable(T_by_variant,fullfile(OUT.table_dir, ...
    'level19_metrics_by_regime_variant.csv'));
writetable(T_by_dx,fullfile(OUT.table_dir,'level19_metrics_by_dx.csv'));
writetable(T_by_M,fullfile(OUT.table_dir,'level19_metrics_by_M.csv'));
writetable(T_by_cs,fullfile(OUT.table_dir,'level19_metrics_by_cs.csv'));
writetable(T_by_frequency,fullfile(OUT.table_dir, ...
    'level19_metrics_by_frequency.csv'));
writetable(T_by_Meff,fullfile(OUT.table_dir, ...
    'level19_metrics_by_M_eff_guess.csv'));
writetable(T_generalization,fullfile(OUT.table_dir, ...
    'level19_generalization_tests.csv'));

save(fullfile(OUT.model_dir,'level19_clean_field_regime_models.mat'), ...
    'MODELS','SPECS','OLD_SPECS','T_q_metrics','T_sws_metrics', ...
    'T_generalization','MC18','PROFILE','-v7.3');

%% Registry manifest

T_manifest_update = table();
metrics_file = fullfile(OUT.table_dir,'level19_sws_metrics.csv');
for i=1:numel(MODELS)
    metric_i = T_sws_metrics(T_sws_metrics.model_name==MODELS(i).model_name,:);
    notes = sprintf(['Grouped condition split; Test MAPE %.4g%%; ', ...
        'HighError>20 %.4g%%.'],metric_i.MAPE_pct(1), ...
        metric_i.HighError_gt20_pct(1));
    if PROFILE=="full"
        entry = adaptive_req.analysis.register_trained_model( ...
            'RootDir',root_dir, ...
            'SourceModelFile',MODELS(i).registry_deployment(1), ...
            'ModelId',"test19__"+MODELS(i).model_name+"__bagged_trees", ...
            'RegistrySubdir','test19_clean_field_regime', ...
            'TestName','test_19_clean_field_regime_model_training', ...
            'AnalysisLevel','analysis', ...
            'ModelName',MODELS(i).model_name, ...
            'FeatureSet',MODELS(i).feature_set, ...
            'ModelType','bagged_trees', ...
            'ModelRole',MODELS(i).model_role, ...
            'TrainingDataset','test_18_clean_field_regime_training_dataset', ...
            'Target','q_theory', ...
            'PredictorSummary',strjoin(MODELS(i).predictors,', '), ...
            'MetricsFile',metrics_file, ...
            'SplitType','grouped_condition_75_25', ...
            'PerformanceSummary',notes, ...
            'Notes',notes);
    else
        entry = manifest_preview_row(MODELS(i),metrics_file,notes);
    end
    entry.split_type = "grouped_condition_75_25";
    entry.performance_summary = string(notes);
    T_manifest_update = concat_tables(T_manifest_update,entry);
end
writetable(T_manifest_update,fullfile(OUT.table_dir, ...
    'level19_model_manifest_update.csv'));

%% Figures

plot_q_scatter(T_predictions,OUT.fig_dir);
plot_metric_bar(T_sws_metrics,"MAPE_pct",'SWS MAPE (%)', ...
    'level19_mape_by_model.png',OUT.fig_dir);
plot_metric_bar(T_sws_metrics,"HighError_gt20_pct",'High-error >20% (%)', ...
    'level19_high_error_gt20_by_model.png',OUT.fig_dir);
plot_group_metric(T_by_regime,"field_regime_label", ...
    'MAPE by field regime','level19_mape_by_field_regime.png',OUT.fig_dir);
plot_group_metric(T_by_variant,"field_regime_variant", ...
    'MAPE by regime variant','level19_mape_by_regime_variant.png',OUT.fig_dir);
plot_group_metric(T_by_dx,"SIM_dx", ...
    'MAPE by spatial resolution','level19_mape_by_dx.png',OUT.fig_dir);
plot_group_metric(T_by_M,"REQ_M", ...
    'MAPE by REQ M','level19_mape_by_REQ_M.png',OUT.fig_dir);
plot_old_new_hybrid(T_sws_metrics,OUT.fig_dir);
plot_generalization(T_generalization,OUT.fig_dir);

%% Console summary

fprintf('\nTest 19 complete. Analysis folder:\n%s\n',OUT.analysis_dir);
fprintf('\nPrimary grouped-test performance:\n');
disp(sortrows(T_sws_metrics(:,{'model_name','model_role','MAPE_pct', ...
    'HighError_gt10_pct','HighError_gt20_pct'}),'MAPE_pct'));
fprintf(['\nMain operational model: HybridLocalGlobal_T18_noUserRegime.\n', ...
    'The withUserRegimeGuess model is user-informed and is reported separately.\n', ...
    'field_regime_label and SIM_Nwaves were not used by no-user models.\n']);

%% Local functions

function OUT = make_output_dirs(root_dir,profile)
OUT.analysis_dir = fullfile(root_dir,'outputs', ...
    'test_19_clean_field_regime_model_training','analysis');
if profile=="smoke", OUT.analysis_dir=fullfile(OUT.analysis_dir,'smoke_test'); end
OUT.fig_dir = fullfile(OUT.analysis_dir,'figures');
OUT.table_dir = fullfile(OUT.analysis_dir,'tables');
OUT.model_dir = fullfile(OUT.analysis_dir,'models');
OUT.job_dir = fullfile(OUT.model_dir,'generalization_jobs');
dirs={OUT.analysis_dir,OUT.fig_dir,OUT.table_dir,OUT.model_dir,OUT.job_dir};
for i=1:numel(dirs), if exist(dirs{i},'dir')~=7, mkdir(dirs{i}); end, end
end

function SPECS = build_test19_specs(T)
base = adaptive_req.analysis.Test12Analysis.buildModelSpecs(T);
local = pick_spec(base,"LocalOnly","NoCsGuess");
global_i = pick_spec(base,"GlobalOnly","NoCsGuess");
hybrid = pick_spec(base,"HybridLocalGlobal","WithCsGuess");
SPECS = [ ...
    make_spec("LocalOnly_T18","CleanFieldRegime_noUser", ...
        "operational",local.predictors,false)
    make_spec("GlobalOnly_T18","CleanFieldRegime_noUser", ...
        "operational",global_i.predictors,false)
    make_spec("HybridLocalGlobal_T18_noUserRegime", ...
        "CleanFieldRegime_noUser","operational",hybrid.predictors,false)
    make_spec("HybridLocalGlobal_T18_withUserRegimeGuess", ...
        "CleanFieldRegime_withUserGuess","user_informed", ...
        unique([hybrid.predictors;"user_field_guess"],'stable'),true)];
end

function out = pick_spec(specs,name,feature)
idx=find([specs.model_name]==name & [specs.feature_set]==feature,1);
assert(~isempty(idx),'Missing predictor spec %s | %s.',name,feature);
out=specs(idx);
end

function S=make_spec(name,feature,role,predictors,allow_user)
S=struct('model_name',string(name),'feature_set',string(feature), ...
    'model_role',string(role),'predictors',string(predictors(:)), ...
    'allow_user_guess',logical(allow_user));
end

function OLD = load_old_models(root_dir)
model_dir=fullfile(root_dir,'outputs','model_registry','test12_hybrid_baseline');
requests=[ ...
    struct('stored',"LocalOnly",'feature',"NoCsGuess",'report',"LocalOnly_old")
    struct('stored',"GlobalOnly",'feature',"NoCsGuess",'report',"GlobalOnly_old")
    struct('stored',"HybridLocalGlobal",'feature',"WithCsGuess",'report',"HybridLocalGlobal_old")];
OLD=struct([]);
for i=1:numel(requests)
    [M,I,F]=adaptive_req.analysis.load_q_model_deployment(model_dir, ...
        'ModelName',requests(i).stored,'FeatureSet',requests(i).feature, ...
        'ModelType','bagged_trees');
    OLD(i).model_name=requests(i).report;
    OLD(i).feature_set="Test12_"+string(I.feature_set);
    OLD(i).model_type="bagged_trees";
    OLD(i).model_role="external_old_model";
    OLD(i).allow_user_guess=false;
    OLD(i).model_file=string(F);
    OLD(i).model=M;
end
end

function assert_operational_predictors(predictors,model_name,allow_user)
p=lower(string(predictors(:)));
banned=lower(["q_theory","q_true","cs_true","cs_pred","sws_error", ...
    "abs_sws_error","M_eff_true_diag","aperture_weight", ...
    "solid_angle_weight","true_aperture_weight","SIM_Nwaves", ...
    "SIM_Is2D","SIM_ForceInPlaneWave","field_regime_variant"]);
banned=banned(:);
if ~allow_user, banned=[banned;"field_regime_label";"user_field_guess"]; end
bad=intersect(p,banned);
bad_pattern=contains(p,"error")|contains(p,"residual")|contains(p,"target");
assert(isempty(bad) && ~any(bad_pattern), ...
    'Leakage/diagnostic predictors in %s: %s',model_name, ...
    strjoin(unique([bad;p(bad_pattern)]),', '));
end

function labels=map_user_field_guess(regime)
regime=string(regime); labels=strings(size(regime));
labels(regime=="directional_2D")="directional_like";
labels(regime=="partial_3D")="partially_diffuse";
labels(regime=="diffuse_2D"|regime=="diffuse_3D")="diffuse_like";
labels(labels=="")="unknown";
labels=categorical(labels,["directional_like","partially_diffuse", ...
    "diffuse_like","unknown"]);
end

function n=mapping_nbins(mapping)
if isempty(mapping)
    n=NaN;
elseif isfield(mapping,'Nbins_effective')
    n=double(mapping.Nbins_effective);
elseif isfield(mapping,'k_cent')
    n=numel(mapping.k_cent);
else
    n=NaN;
end
end

function key=make_test18_row_key(T)
key=string(T.condition_id)+"|M"+string(T.REQ_M)+"|R"+ ...
    string(T.realization_idx)+"|P"+string(T.patch_idx);
end

function [train_mask,test_mask]=grouped_condition_split(T,fraction,seed)
rng(seed); groups=unique(T.condition_id,'stable');
order=groups(randperm(numel(groups)));
n=max(1,min(numel(groups)-1,round(fraction*numel(groups))));
train_mask=ismember(T.condition_id,order(1:n)); test_mask=~train_mask;
end

function Tpred=enrich_predictions(Tpred,Tref,spec)
Tpred.row_key=make_test18_row_key(Tpred);
[tf,loc]=ismember(Tpred.row_key,Tref.row_key);
assert(all(tf),'Could not map prediction rows to Test 18 rows.');
meta=["field_regime_label","field_regime_variant","SIM_Nwaves", ...
    "SIM_Is2D","SIM_ForceInPlaneWave","SIM_dx","SIM_dz", ...
    "REQ_cs_guess","M_eff_guess","M_eff_guess_bin"];
meta=meta(ismember(meta,string(Tref.Properties.VariableNames)));
for meta_idx=1:numel(meta)
    v=meta(meta_idx);
    source_values=Tref.(char(v));
    Tpred.(char(v))=source_values(loc);
end
Tpred.feature_set=repmat(spec.feature_set,height(Tpred),1);
Tpred.model_role=repmat(spec.model_role,height(Tpred),1);
Tpred.model_type=repmat("bagged_trees",height(Tpred),1);
Tpred.cs_true=Tref.SIM_cs_bg(loc);
Tpred.cs_pred=adaptive_req.analysis.Test12Analysis.qToCs( ...
    Tpred.q_pred,Tref.req_mapping(loc),Tref.SIM_f0(loc));
Tpred.sws_error_pct=100*(Tpred.cs_pred-Tpred.cs_true)./Tpred.cs_true;
Tpred.abs_sws_error_pct=abs(Tpred.sws_error_pct);
end

function Tpred=baseline_prediction_table(T,q,spec,model_type)
keep=["condition_id","realization_idx","patch_idx","REQ_M", ...
    "SIM_f0","SIM_cs_bg","SIM_dx","SIM_dz","field_regime_label", ...
    "field_regime_variant","SIM_Nwaves","M_eff_guess","M_eff_guess_bin"];
keep=keep(ismember(keep,string(T.Properties.VariableNames)));
Tpred=T(:,cellstr(keep)); Tpred.row_key=T.row_key;
Tpred.model_name=repmat(spec.model_name,height(T),1);
Tpred.feature_set=repmat(spec.feature_set,height(T),1);
Tpred.model_role=repmat(spec.model_role,height(T),1);
Tpred.model_type=repmat(string(model_type),height(T),1);
Tpred.split=repmat("external",height(T),1);
Tpred.q_true=T.q_theory; Tpred.q_pred_raw=q; Tpred.q_pred=q;
Tpred.residual=q-T.q_theory; Tpred.abs_error=abs(Tpred.residual);
Tpred.cs_true=T.SIM_cs_bg;
Tpred.cs_pred=adaptive_req.analysis.Test12Analysis.qToCs(q,T.req_mapping,T.SIM_f0);
Tpred.sws_error_pct=100*(Tpred.cs_pred-Tpred.cs_true)./Tpred.cs_true;
Tpred.abs_sws_error_pct=abs(Tpred.sws_error_pct);
end

function T=compact_prediction_table(T)
keep=["condition_id","realization_idx","patch_idx","REQ_M", ...
    "SIM_f0","SIM_cs_bg","SIM_dx","SIM_dz","REQ_cs_guess", ...
    "M_eff_guess","M_eff_guess_bin","field_regime_label", ...
    "field_regime_variant","SIM_Nwaves","row_key","model_name", ...
    "feature_set","model_role","model_type","split","q_true", ...
    "q_pred_raw","q_pred","residual","abs_error","cs_true", ...
    "cs_pred","sws_error_pct","abs_sws_error_pct"];
keep=keep(ismember(keep,string(T.Properties.VariableNames)));
T=T(:,cellstr(keep));
end

function q=compute_discrete_theory_q(T)
q=nan(height(T),1); cache=containers.Map('KeyType','char','ValueType','double');
for i=1:height(T)
    label=string(T.field_regime_label(i));
    key=sprintf('%s_dx%.12g_f%g_M%g',label,T.SIM_dx(i),T.SIM_f0(i),T.REQ_M(i));
    if isKey(cache,key), q(i)=cache(key); continue; end
    if label=="partial_3D"
        q_i=.5*(theory_one(T,i,"Diffuse2D")+theory_one(T,i,"Diffuse3D"));
    elseif label=="directional_2D"
        q_i=theory_one(T,i,"SingleWave");
    elseif label=="diffuse_2D"
        q_i=theory_one(T,i,"Diffuse2D");
    else
        q_i=theory_one(T,i,"Diffuse3D");
    end
    cache(key)=q_i; q(i)=q_i;
end
end

function q=theory_one(T,i,field_type)
out=adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    T.SIM_dx(i),T.SIM_dz(i),T.SIM_f0(i),T.REQ_cs_guess(i), ...
    'M',T.REQ_M(i),'Gamma',1,'PadFactor',1,'Nbins','auto', ...
    'SmoothSigma',1,'TheoryMode','S2D','FieldType',field_type,'Plot',false);
q=out.q_th;
end

function Tsum=summarize_q(T,groups)
[G,Tsum]=findgroups(T(:,cellstr(groups)));
Tsum.N=splitapply(@numel,T.q_true,G);
Tsum.MAE_q=splitapply(@(a,b)mean(abs(b-a),'omitnan'),T.q_true,T.q_pred,G);
Tsum.RMSE_q=splitapply(@(a,b)sqrt(mean((b-a).^2,'omitnan')),T.q_true,T.q_pred,G);
Tsum.MAPE_q_pct=splitapply(@(a,b)100*mean(abs((b-a)./a),'omitnan'),T.q_true,T.q_pred,G);
Tsum.bias_q=splitapply(@(a,b)mean(b-a,'omitnan'),T.q_true,T.q_pred,G);
end

function Tsum=summarize_sws(T,groups)
[G,Tsum]=findgroups(T(:,cellstr(groups)));
Tsum.N=splitapply(@numel,T.abs_sws_error_pct,G);
Tsum.MAPE_pct=splitapply(@(x)mean(x,'omitnan'),T.abs_sws_error_pct,G);
Tsum.RMSE_pct=splitapply(@(x)sqrt(mean(x.^2,'omitnan')),T.sws_error_pct,G);
Tsum.MedAE_pct=splitapply(@(x)median(x,'omitnan'),T.abs_sws_error_pct,G);
Tsum.bias_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_error_pct,G);
Tsum.HighError_gt10_pct=splitapply(@(x)100*mean(x>10,'omitnan'),T.abs_sws_error_pct,G);
Tsum.HighError_gt20_pct=splitapply(@(x)100*mean(x>20,'omitnan'),T.abs_sws_error_pct,G);
Tsum=sortrows(Tsum,'MAPE_pct');
end

function Tout=generalization_summary(T,split_name,heldout_var,heldout_value,ntrain,ntest)
Tout=summarize_sws(T,["model_name","feature_set","model_role","model_type"]);
Q=summarize_q(T,["model_name","feature_set","model_role","model_type"]);
Tout.MAE_q=Q.MAE_q; Tout.RMSE_q=Q.RMSE_q; Tout.MAPE_q_pct=Q.MAPE_q_pct;
Tout.generalization_test=repmat(string(split_name),height(Tout),1);
Tout.heldout_var=repmat(string(heldout_var),height(Tout),1);
Tout.heldout_value=repmat(value_string(heldout_value),height(Tout),1);
Tout.N_train=ntrain*ones(height(Tout),1); Tout.N_test=ntest*ones(height(Tout),1);
Tout=movevars(Tout,{'generalization_test','heldout_var','heldout_value'},'Before',1);
end

function tf=is_group_value(x,value)
if isstring(x)||iscategorical(x)||iscellstr(x), tf=string(x)==string(value);
else, tf=x==value; end
tf=tf(:);
end

function assert_grouped_split(T,train_mask,test_mask,var,value)
assert(any(train_mask)&&any(test_mask)&&~any(train_mask&test_mask));
x=T.(char(var));
assert(~any(is_group_value(x(train_mask),value)) && ...
    all(is_group_value(x(test_mask),value)), ...
    'Invalid grouped split for %s=%s.',var,value_string(value));
assert(isempty(intersect(unique(T.condition_id(train_mask)), ...
    unique(T.condition_id(test_mask)))), ...
    'Physical condition IDs overlap in grouped split %s=%s.', ...
    var,value_string(value));
end

function [train_mask,test_mask]=make_generalization_split(T,var,value,seed)
is_holdout=is_group_value(T.(char(var)),value);
if var~="REQ_M"
    train_mask=~is_holdout;
    test_mask=is_holdout;
    return;
end

% Every Test 18 physical condition contains all M values. Split physical
% conditions first, then use non-heldout M only in training and heldout M
% only in test. This prevents repeated global features from crossing splits.
rng(seed);
conditions=unique(T.condition_id,'stable');
conditions=conditions(randperm(numel(conditions)));
n_train=max(1,min(numel(conditions)-1,round(.75*numel(conditions))));
train_conditions=conditions(1:n_train);
test_conditions=conditions(n_train+1:end);
train_mask=ismember(T.condition_id,train_conditions) & ~is_holdout;
test_mask=ismember(T.condition_id,test_conditions) & is_holdout;
end

function s=value_string(x)
if isnumeric(x), s=string(sprintf('%.12g',x)); else, s=string(x); end
end

function require_vars(T,vars,context)
missing=setdiff(string(vars),string(T.Properties.VariableNames));
assert(isempty(missing),'%s missing: %s',context,strjoin(missing,', '));
end

function T=concat_tables(A,B)
if isempty(A),T=B;return;end
if isempty(B),T=A;return;end
vars=unique([string(A.Properties.VariableNames),string(B.Properties.VariableNames)],'stable');
A=add_missing_from(A,B,vars); B=add_missing_from(B,A,vars);
T=[A(:,cellstr(vars));B(:,cellstr(vars))];
end

function T=add_missing_from(T,prototype,vars)
for i=1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)),continue;end
    source=prototype.(char(vars(i)));
    if iscategorical(source)
        T.(char(vars(i)))=categorical(strings(height(T),1),categories(source));
    elseif isstring(source)
        T.(char(vars(i)))=strings(height(T),1);
    elseif islogical(source)
        T.(char(vars(i)))=false(height(T),1);
    elseif iscell(source)
        T.(char(vars(i)))=cell(height(T),1);
    else
        T.(char(vars(i)))=nan(height(T),1);
    end
end
end

function T=remove_cell_columns(T)
drop=false(1,width(T));
for i=1:width(T), drop(i)=iscell(T.(T.Properties.VariableNames{i}))||isstruct(T.(T.Properties.VariableNames{i})); end
T(:,drop)=[];
end

function plot_q_scatter(T,fig_dir)
models=unique(T.model_name,'stable');
fig=figure('Color','w','Units','centimeters','Position',[2 2 25 16]);
tl=tiledlayout(ceil(numel(models)/4),4,'TileSpacing','compact','Padding','compact');
for i=1:numel(models)
    ax=nexttile(tl); Ti=T(T.model_name==models(i),:);
    scatter(ax,Ti.q_true,Ti.q_pred,6,'filled','MarkerFaceAlpha',.2); hold(ax,'on');
    plot(ax,[0 1],[0 1],'k--'); axis(ax,'equal'); xlim(ax,[0 1]);ylim(ax,[0 1]);grid(ax,'on');
    title(ax,models(i),'Interpreter','none','FontWeight','normal','FontSize',9);
    xlabel(ax,'q true');ylabel(ax,'q predicted');
end
title(tl,'Test 19 q true vs predicted','FontWeight','normal');
export_clean(fig,fullfile(fig_dir,'level19_q_true_vs_q_pred_by_model.png'));
end

function plot_metric_bar(T,var,ylabel_text,file,fig_dir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 22 10]);
bar(categorical(T.model_name),T.(char(var))); grid on; xtickangle(30);
ylabel(ylabel_text); title('Test 19 primary grouped test','FontWeight','normal');
export_clean(fig,fullfile(fig_dir,file));
end

function plot_group_metric(T,xvar,title_text,file,fig_dir)
x=string(T.(char(xvar))); models=unique(T.model_name,'stable'); xv=unique(x,'stable');
Y=nan(numel(xv),numel(models));
for i=1:numel(xv)
    for j=1:numel(models)
        idx=x==xv(i)&T.model_name==models(j);
        if any(idx)
            Y(i,j)=mean(T.MAPE_pct(idx),'omitnan');
        end
    end
end
fig=figure('Color','w','Units','centimeters','Position',[2 2 24 11]);
bar(categorical(xv),Y);grid on;ylabel('SWS MAPE (%)');xtickangle(25);
title(title_text,'Interpreter','none','FontWeight','normal');
legend(models,'Location','eastoutside','Interpreter','none');
export_clean(fig,fullfile(fig_dir,file));
end

function plot_old_new_hybrid(T,fig_dir)
keep=ismember(T.model_name,["HybridLocalGlobal_old", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess"]);
plot_metric_bar(T(keep,:),"MAPE_pct",'SWS MAPE (%)', ...
    'level19_old_vs_new_hybrid.png',fig_dir);
end

function plot_generalization(T,fig_dir)
models=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess"];
T=T(ismember(T.model_name,models) & ...
    T.generalization_test=="leave_one_field_regime_out",:);
plot_group_metric(T,"heldout_value", ...
    'Leave-one-group-out generalization','level19_leave_one_regime_out_summary.png',fig_dir);
end

function entry=manifest_preview_row(MODEL,metrics_file,notes)
entry=table();
entry.model_id="smoke__test19__"+MODEL.model_name+"__bagged_trees";
entry.test_name="test_19_smoke_validation";
entry.analysis_level="smoke_test";
entry.model_name=MODEL.model_name;
entry.feature_set=MODEL.feature_set;
entry.model_type="bagged_trees";
entry.model_role=MODEL.model_role;
entry.training_dataset="test_18_smoke_test";
entry.target="q_theory";
entry.predictor_summary=strjoin(MODEL.predictors,', ');
entry.model_file=MODEL.local_deployment(1);
entry.metrics_file=string(metrics_file);
entry.split_type="grouped_condition_75_25";
entry.performance_summary=string(notes);
entry.created_datetime=string(datetime('now'));
entry.notes="Smoke validation only; not registered as a deployable final model.";
end

function export_clean(fig,path_i)
axs=findall(fig,'Type','axes');
for i=1:numel(axs)
    try
        axs(i).Toolbar.Visible='off';
    catch
    end
end
drawnow;exportgraphics(fig,path_i,'Resolution',240,'BackgroundColor','white');close(fig);
end

function s=sanitize_filename(x)
s=regexprep(string(x),'[^A-Za-z0-9_=-]+','_');
end
