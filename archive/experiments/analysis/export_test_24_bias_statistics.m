%% export_test_24_bias_statistics.m
% Export homogeneous and heterogeneous non-mixed SWS bias statistics.
%
% Homogeneous rows come from frozen Test 20 predictions. Heterogeneous rows
% come from the Test 23 joined diagnostic cache and exclude mixed patches:
% patch_soft_fraction > 0.05 AND patch_hard_fraction > 0.05.
% No model is trained or run.

clear; clc;

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();

SOURCE20=fullfile(root_dir,'outputs', ...
    'test_20_external_validation_and_aperture_q_tracking','analysis', ...
    'tables','level20_external_predictions.csv');
SOURCE23=fullfile(root_dir,'outputs','test_23_interface_patch_contamination', ...
    'analysis','level_23_interface_patch_contamination_step01','data', ...
    'test23_joined_predictions.mat');
OUT_DIR=fullfile(root_dir,'outputs','test_24_interface_spectral_failure_modes','tables');
assert(exist(SOURCE20,'file')==2,'Missing Test 20 predictions: %s',SOURCE20);
assert(exist(SOURCE23,'file')==2,'Missing Test 23 cache: %s',SOURCE23);
if exist(OUT_DIR,'dir')~=7, mkdir(OUT_DIR); end

models=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","TheoryQDiscrete", ...
    "LocalOnly_old","GlobalOnly_old","HybridLocalGlobal_old"];

%% Homogeneous frozen predictions

H=read_homogeneous_predictions(SOURCE20,models);
H.region_label=H.geometry_id;
H.material_side=repmat("not_applicable",height(H),1);
H.sws_signed_error_pct=100*(H.cs_pred-H.cs_true)./H.cs_true;
T_hom=summarize_bias(H,["model_name","geometry_id","region_label"]);
T_hom.context=repmat("homogeneous",height(T_hom),1);
T_hom.material_side=repmat("not_applicable",height(T_hom),1);
T_hom=movevars(T_hom,["context","material_side"],'Before',1);

%% Heterogeneous non-mixed predictions

S=load(SOURCE23,'T_joined'); X=S.T_joined; clear S;
X=X(ismember(X.model_name,models),:);
X.patch_is_mixed=X.patch_soft_fraction>0.05 & X.patch_hard_fraction>0.05;
X=X(~X.patch_is_mixed,:);
X.material_side=repmat("soft",height(X),1);
X.material_side(X.sws_true>2.5)="hard";
X.sws_signed_error_pct=100*(X.sws_pred-X.sws_true)./X.sws_true;

T_het=summarize_bias(X,["model_name","material_side"]);
T_het.context=repmat("heterogeneous_nonmixed",height(T_het),1);
T_het.region_label="heterogeneous_"+T_het.material_side;
T_het=movevars(T_het,["context","region_label"],'Before',1);

T_het_case=summarize_bias(X,["model_name","geometry_case","material_side"]);
T_het_case.context=repmat("heterogeneous_nonmixed",height(T_het_case),1);
T_het_case.region_label=T_het_case.geometry_case+"_"+T_het_case.material_side;
T_het_case=movevars(T_het_case,["context","region_label"],'Before',1);

%% Unified apples-to-apples table

C_hom=standardize_comparison(T_hom);
C_het=standardize_comparison(T_het);
C_case=standardize_comparison(T_het_case);
T_comparison=[C_hom;C_het;C_case];
T_comparison=sortrows(T_comparison,{'model_name','context','region_label'});

%% Write CSVs

files=[ ...
    string(fullfile(OUT_DIR,'test24_homogeneous_bias_statistics.csv'))
    string(fullfile(OUT_DIR,'test24_heterogeneous_nonmixed_bias_statistics.csv'))
    string(fullfile(OUT_DIR,'test24_heterogeneous_nonmixed_bias_by_case.csv'))
    string(fullfile(OUT_DIR,'test24_homogeneous_vs_heterogeneous_bias_comparison.csv'))];
writetable(T_hom,files(1));
writetable(T_het,files(2));
writetable(T_het_case,files(3));
writetable(T_comparison,files(4));

fprintf('\nBias statistics exported without training or inference:\n');
for i=1:numel(files), fprintf('%s\n',files(i)); end

%% Local functions

function T=read_homogeneous_predictions(file,models)
keep={'geometry_id','model_name','cs_true','cs_pred'};
ds=tabularTextDatastore(file,'TextType','string');
ds.SelectedVariableNames=keep;
parts=cell(0,1); n=0;
while hasdata(ds)
    A=read(ds); A.geometry_id=string(A.geometry_id); A.model_name=string(A.model_name);
    idx=ismember(A.geometry_id,["homogeneous_cs2","homogeneous_cs3"]) & ...
        ismember(A.model_name,models);
    if any(idx), n=n+1; parts{n,1}=A(idx,:); end
end
assert(~isempty(parts),'No homogeneous predictions found.');
T=vertcat(parts{:});
end

function S=summarize_bias(T,groups)
[G,S]=findgroups(T(:,cellstr(groups)));
e=T.sws_signed_error_pct;
S.N=splitapply(@numel,e,G);
S.mean_signed_error_pct=splitapply(@(x)mean(x,'omitnan'),e,G);
S.median_signed_error_pct=splitapply(@(x)median(x,'omitnan'),e,G);
S.underestimate_pct=100*splitapply(@(x)mean(x<0,'omitnan'),e,G);
S.overestimate_pct=100*splitapply(@(x)mean(x>0,'omitnan'),e,G);
S.mean_abs_error_pct=splitapply(@(x)mean(abs(x),'omitnan'),e,G);
S.high_error_gt10_pct=100*splitapply(@(x)mean(abs(x)>10,'omitnan'),e,G);
S.high_error_gt20_pct=100*splitapply(@(x)mean(abs(x)>20,'omitnan'),e,G);
end

function C=standardize_comparison(T)
vars=["context","region_label","model_name","material_side","N", ...
    "mean_signed_error_pct","median_signed_error_pct", ...
    "underestimate_pct","overestimate_pct","mean_abs_error_pct", ...
    "high_error_gt10_pct","high_error_gt20_pct"];
if ~ismember('material_side',T.Properties.VariableNames)
    T.material_side=repmat("not_applicable",height(T),1);
end
C=T(:,cellstr(vars));
end
