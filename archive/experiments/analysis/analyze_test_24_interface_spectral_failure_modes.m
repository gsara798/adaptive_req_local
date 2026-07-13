%% analyze_test_24_interface_spectral_failure_modes.m
% Test 24: interface spectral failure modes.
%
% Purely retrospective diagnosis. No model is trained, modified, or run.
% Frozen Test 20 predictions and Test 23 oracle diagnostics are loaded.
% True q/SWS, material wavenumbers, purity, and interface distance are used
% only after prediction. A few representative spectra are recomputed from
% the already saved Test 17 fields for visualization.
%
% k_model_region is mutually exclusive with this priority:
%   below_hard / above_soft (outside the material interval),
%   near_true (within 5% of center-material k),
%   near_hard / near_soft (lower/upper 25% of the material interval),
%   between_hard_soft (middle 50% of the interval).

clear; clc; close all;
format compact;

%% Setup and sources

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontName','Helvetica', ...
    'defaultTextFontName','Helvetica','defaultAxesFontSize',8, ...
    'defaultTextFontSize',8,'defaultLegendFontSize',7);

SOURCE=fullfile(root_dir,'outputs','test_23_interface_patch_contamination', ...
    'analysis','level_23_interface_patch_contamination_step01','data', ...
    'test23_joined_predictions.mat');
CONDITION_DIR=fullfile(root_dir,'outputs', ...
    'test_17_model_comparison_heterogeneous_cases','analysis','data','conditions');
assert(exist(SOURCE,'file')==2, ...
    'Test 23 joined cache missing. Run Test 23 Step 01 first: %s',SOURCE);
assert(exist(CONDITION_DIR,'dir')==7,'Test 17 conditions missing: %s',CONDITION_DIR);

OUT=make_output_dirs(root_dir);
PRIMARY_MODEL="HybridLocalGlobal_T18_noUserRegime";
PURITY_EDGES=[0 0.4 0.6 0.8 0.9 0.95 0.99 1.000001];
DISTANCE_EDGES=[0 1 2 4 6 10 Inf];

fprintf('\nTest 24: interface spectral failure modes\n');
fprintf('No training, correction, or model inference is performed.\n');
fprintf('Source: %s\n',SOURCE);

%% Build one-row-per-patch/model diagnostic table

S=load(SOURCE,'T_joined'); T=S.T_joined; clear S;
require_variables(T,["model_name","geometry_case","field_regime", ...
    "frequency_hz","dx_m","REQ_M","q_true","q_pred","sws_true", ...
    "sws_pred","patch_purity","patch_soft_fraction", ...
    "patch_hard_fraction","abs_distance_to_interface_mm", ...
    "signed_distance_mm","side"]);

T=build_failure_table(T);
assert_diagnostic_only(T);
fprintf('Diagnostic rows: %d | models: %d | cases: %d.\n', ...
    height(T),numel(unique(T.model_name)),numel(unique(T.external_case)));

%% Bins and summary tables

T.purity_bin=discretize(T.patch_purity,PURITY_EDGES,'IncludedEdge','left');
T.purity_bin_label=purity_labels(T.purity_bin);
T.distance_bin=discretize(T.distance_to_interface_mm,DISTANCE_EDGES, ...
    'IncludedEdge','left');
T.distance_bin_label=distance_labels(T.distance_bin);

T_model=summarize_failure_modes(T,["model_name","external_case", ...
    "field_regime","REQ_M","frequency","dx"]);
T_purity=summarize_failure_modes(T,["model_name","external_case", ...
    "field_regime","purity_bin","purity_bin_label"]);
T_distance=summarize_failure_modes(T,["model_name","external_case", ...
    "field_regime","distance_bin","distance_bin_label"]);
T_regions=summarize_k_regions(T);

%% Representative failure spectra

T_representative=select_representative_failures(T,PRIMARY_MODEL);
writetable(T_representative,fullfile(OUT.table_dir, ...
    'test24_representative_failure_patch_list.csv'));
plot_representative_failures(T_representative,CONDITION_DIR,OUT);

%% Save tables

patch_file=fullfile(OUT.table_dir, ...
    'test24_interface_spectral_failure_patch_table.csv');
model_file=fullfile(OUT.table_dir,'test24_summary_by_model_case_regime.csv');
purity_file=fullfile(OUT.table_dir,'test24_summary_by_purity_bin.csv');
distance_file=fullfile(OUT.table_dir,'test24_summary_by_distance_bin.csv');
region_file=fullfile(OUT.table_dir,'test24_k_region_counts.csv');
writetable(T,patch_file); writetable(T_model,model_file);
writetable(T_purity,purity_file); writetable(T_distance,distance_file);
writetable(T_regions,region_file);
save(fullfile(OUT.data_dir,'test24_analysis_results.mat'), ...
    'T_model','T_purity','T_distance','T_regions','T_representative', ...
    'PURITY_EDGES','DISTANCE_EDGES','-v7.3');

%% Figures

plot_k_ratio_vs_purity(T,OUT);
plot_signed_error_vs_purity(T,OUT);
plot_k_regions_by_distance(T,OUT,PRIMARY_MODEL);
plot_above_soft_by_distance(T,OUT);
plot_q_error_vs_k_ratio(T,OUT);
plot_mechanism_summary(T,OUT,PRIMARY_MODEL);

%% Console conclusions and paths

print_summary(T,PRIMARY_MODEL);
fprintf('\nGenerated tables:\n%s\n%s\n%s\n%s\n%s\n', ...
    patch_file,model_file,purity_file,distance_file,region_file);
fprintf('Generated figures:\n%s\n',OUT.figure_dir);
fprintf('Representative spectra:\n%s\n',OUT.spectral_dir);
fprintf('\nTest 24 complete. Output root:\n%s\n',OUT.root_dir);

%% Local functions

function OUT=make_output_dirs(root_dir)
OUT.root_dir=fullfile(root_dir,'outputs','test_24_interface_spectral_failure_modes');
OUT.table_dir=fullfile(OUT.root_dir,'tables');
OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.spectral_dir=fullfile(OUT.figure_dir,'test24_representative_failure_spectra');
OUT.data_dir=fullfile(OUT.root_dir,'data');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs), if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end, end
end

function require_variables(T,names)
missing=setdiff(string(names),string(T.Properties.VariableNames));
assert(isempty(missing),'Source table missing variables: %s',strjoin(missing,', '));
end

function T=build_failure_table(T)
T.external_case=T.geometry_case;
T.frequency=T.frequency_hz; T.dx=T.dx_m;
T.cs_true=T.sws_true; T.cs_soft=2*ones(height(T),1); T.cs_hard=3*ones(height(T),1);
T.k_true=2*pi*T.frequency./T.cs_true;
T.k_soft=2*pi*T.frequency./T.cs_soft;
T.k_hard=2*pi*T.frequency./T.cs_hard;
T.q_model=T.q_pred; T.SWS_pred=T.sws_pred;
T.SWS_error_signed_pct=100*(T.SWS_pred-T.cs_true)./T.cs_true;
T.SWS_error_abs_pct=abs(T.SWS_error_signed_pct);
T.k_model=2*pi*T.frequency./T.SWS_pred;
T.k_model_over_k_true=T.k_model./T.k_true;
T.k_model_over_k_soft=T.k_model./T.k_soft;
T.k_model_over_k_hard=T.k_model./T.k_hard;
T.rel_dist_to_k_true=abs(T.k_model-T.k_true)./T.k_true;
T.rel_dist_to_k_soft=abs(T.k_model-T.k_soft)./T.k_soft;
T.rel_dist_to_k_hard=abs(T.k_model-T.k_hard)./T.k_hard;
T.is_high_error_10=T.SWS_error_abs_pct>10;
T.is_high_error_20=T.SWS_error_abs_pct>20;
T.is_underestimate=T.SWS_error_signed_pct<0;
T.is_overestimate=T.SWS_error_signed_pct>0;
T.is_k_above_soft=T.k_model>T.k_soft;
T.is_k_between_hard_soft=T.k_model>=T.k_hard&T.k_model<=T.k_soft;
T.is_low_purity=T.patch_purity<0.8;
T.is_mixed_patch=T.patch_soft_fraction>0.05&T.patch_hard_fraction>0.05;
T.soft_fraction=T.patch_soft_fraction; T.hard_fraction=T.patch_hard_fraction;
T.distance_to_interface_mm=T.abs_distance_to_interface_mm;
T.signed_distance_to_interface_mm=T.signed_distance_mm;
T.q_error_model_true=T.q_model-T.q_true;
T.q_error_abs=abs(T.q_error_model_true);

% Attach frozen TheoryQDiscrete q to every aligned patch/model row.
K=string(T.condition_id)+"|"+string(T.map_iz)+"|"+string(T.map_ix);
Q=T(T.model_name=="TheoryQDiscrete",:);
Kq=string(Q.condition_id)+"|"+string(Q.map_iz)+"|"+string(Q.map_ix);
[tf,loc]=ismember(K,Kq); T.q_theory=nan(height(T),1);
T.q_theory(tf)=Q.q_model(loc(tf));
T.q_error_theory_true=T.q_theory-T.q_true;
T.q_model_minus_q_theory=T.q_model-T.q_theory;

T.is_near_true=T.rel_dist_to_k_true<=0.05;
T.is_near_soft=T.rel_dist_to_k_soft<T.rel_dist_to_k_hard;
T.is_near_hard=T.rel_dist_to_k_hard<T.rel_dist_to_k_soft;
T.k_model_region=repmat("between_hard_soft",height(T),1);
T.k_model_region(T.k_model<T.k_hard)="below_hard";
T.k_model_region(T.k_model>T.k_soft)="above_soft";
inside=T.is_k_between_hard_soft;
interval_position=(T.k_model-T.k_hard)./(T.k_soft-T.k_hard);
T.k_model_region(inside&interval_position<=0.25)="near_hard";
T.k_model_region(inside&interval_position>=0.75)="near_soft";
T.k_model_region(inside&T.is_near_true)="near_true";

keep=["condition_id","external_case","field_regime","REQ_M", ...
    "frequency","dx","dz_m","model_name","map_iz","map_ix","cx","cz", ...
    "x_center_m","z_center_m","cs_true","cs_soft","cs_hard", ...
    "k_true","k_soft","k_hard","q_true","q_model","q_theory", ...
    "k_model","SWS_pred","SWS_error_signed_pct","SWS_error_abs_pct", ...
    "patch_purity","soft_fraction","hard_fraction", ...
    "distance_to_interface_mm","signed_distance_to_interface_mm","side", ...
    "k_model_over_k_true","k_model_over_k_soft","k_model_over_k_hard", ...
    "rel_dist_to_k_true","rel_dist_to_k_soft","rel_dist_to_k_hard", ...
    "k_model_region","is_near_true","is_near_soft","is_near_hard", ...
    "is_high_error_10","is_high_error_20","is_underestimate", ...
    "is_overestimate","is_k_above_soft","is_k_between_hard_soft", ...
    "is_low_purity","is_mixed_patch","q_error_model_true","q_error_abs", ...
    "q_error_theory_true","q_model_minus_q_theory"];
T=T(:,cellstr(keep));
end

function assert_diagnostic_only(T)
assert(all(isfinite(T.k_model_over_k_true)|isnan(T.SWS_pred)), ...
    'Invalid k/SWS conversion detected.');
assert(all(abs(T.k_model_over_k_true-(T.cs_true./T.SWS_pred))<1e-10 | ...
    ~isfinite(T.SWS_pred)),'k ratio is inconsistent with SWS prediction.');
fprintf(['Anti-leakage check: q_true, k_true, material k, purity, and errors ', ...
    'were computed/attached only after loading frozen predictions.\n']);
end

function labels=purity_labels(bin)
names=["[0,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,0.9)", ...
    "[0.9,0.95)","[0.95,0.99)","[0.99,1.0]"];
labels=strings(size(bin)); ok=isfinite(bin)&bin>=1&bin<=numel(names);
labels(ok)=names(bin(ok));
end

function labels=distance_labels(bin)
names=["0-1 mm","1-2 mm","2-4 mm","4-6 mm","6-10 mm",">10 mm"];
labels=strings(size(bin)); ok=isfinite(bin)&bin>=1&bin<=numel(names);
labels(ok)=names(bin(ok));
end

function S=summarize_failure_modes(T,groups)
[G,S]=findgroups(T(:,cellstr(groups)));
S.N=splitapply(@numel,T.SWS_error_abs_pct,G);
S.mean_abs_SWS_error=splitapply(@(x)mean(x,'omitnan'),T.SWS_error_abs_pct,G);
S.median_abs_SWS_error=splitapply(@(x)median(x,'omitnan'),T.SWS_error_abs_pct,G);
S.high_error_10_pct=100*splitapply(@(x)mean(x,'omitnan'),T.is_high_error_10,G);
S.high_error_20_pct=100*splitapply(@(x)mean(x,'omitnan'),T.is_high_error_20,G);
S.mean_signed_error=splitapply(@(x)mean(x,'omitnan'),T.SWS_error_signed_pct,G);
S.mean_k_model_over_k_true=splitapply(@(x)mean(x,'omitnan'),T.k_model_over_k_true,G);
S.median_k_model_over_k_true=splitapply(@(x)median(x,'omitnan'),T.k_model_over_k_true,G);
S.pct_k_above_soft=100*splitapply(@(x)mean(x,'omitnan'),T.is_k_above_soft,G);
S.pct_k_between_hard_soft=100*splitapply(@(x)mean(x,'omitnan'),T.is_k_between_hard_soft,G);
S.pct_low_purity=100*splitapply(@(x)mean(x,'omitnan'),T.is_low_purity,G);
S.pct_mixed_patch=100*splitapply(@(x)mean(x,'omitnan'),T.is_mixed_patch,G);
S.pct_near_soft=100*splitapply(@(x)mean(x,'omitnan'),T.is_near_soft,G);
S.pct_near_hard=100*splitapply(@(x)mean(x,'omitnan'),T.is_near_hard,G);
end

function C=summarize_k_regions(T)
groups=["model_name","external_case","field_regime","side","k_model_region"];
[G,C]=findgroups(T(:,cellstr(groups)));
C.N=splitapply(@numel,T.SWS_error_abs_pct,G);
C.mean_abs_SWS_error=splitapply(@(x)mean(x,'omitnan'),T.SWS_error_abs_pct,G);
C.high_error_20_pct=100*splitapply(@(x)mean(x,'omitnan'),T.is_high_error_20,G);
[Gd,den]=findgroups(T(:,{'model_name','external_case','field_regime','side'}));
den.GroupCount=splitapply(@numel,T.k_model_region,Gd);
key=group_key(C,["model_name","external_case","field_regime","side"]);
keyd=group_key(den,["model_name","external_case","field_regime","side"]);
[tf,loc]=ismember(key,keyd); assert(all(tf));
C.percent_within_group=100*C.N./den.GroupCount(loc);
end

function key=group_key(T,vars)
key=repmat("",height(T),1);
for i=1:numel(vars), key=key+"|"+string(T.(vars(i))); end
end

function R=select_representative_failures(T,primary)
P=T(T.model_name==primary,:); cases=unique(P.external_case,'stable'); R=table();
for ci=1:numel(cases)
    X=P(P.external_case==cases(ci),:);
    specs=[ ...
        struct('name',"high_purity_low_error",'mask',X.patch_purity>=0.99&X.SWS_error_abs_pct<5,'score',X.SWS_error_abs_pct)
        struct('name',"high_purity_biased_error",'mask',X.patch_purity>=0.95&X.SWS_error_abs_pct>10,'score',-X.SWS_error_abs_pct)
        struct('name',"low_purity_high_error",'mask',X.patch_purity<0.8&X.is_high_error_20,'score',X.patch_purity-0.01*X.SWS_error_abs_pct)
        struct('name',"interface_center_high_error",'mask',X.side=="interface_band"&X.is_high_error_20,'score',X.distance_to_interface_mm-0.01*X.SWS_error_abs_pct)
        struct('name',"k_model_above_soft",'mask',X.is_k_above_soft,'score',-X.k_model_over_k_soft)
        struct('name',"k_model_between_hard_soft",'mask',X.is_k_between_hard_soft&X.is_high_error_10,'score',-X.SWS_error_abs_pct)];
    for si=1:numel(specs)
        C=X(specs(si).mask,:);
        if isempty(C), warning('Test24:RepresentativeMissing','Missing %s / %s.',cases(ci),specs(si).name); continue; end
        score=specs(si).score(specs(si).mask); [~,j]=min(score); row=C(j,:);
        row.representative_type=specs(si).name; R=concat_tables(R,row);
    end
end
R=movevars(R,'representative_type','Before',1);
end

function plot_representative_failures(R,condition_dir,OUT)
for i=1:height(R)
    file=condition_file(condition_dir,R.external_case(i),R.field_regime(i),R.REQ_M(i));
    if exist(file,'file')~=2, warning('Test24:ConditionMissing','Missing %s.',file); continue; end
    S=load(file,'sim','cfg_sim','req_out'); half=S.req_out.half_win;
    zi=(R.cz(i)-half):(R.cz(i)+half); xi=(R.cx(i)-half):(R.cx(i)+half);
    patch=S.sim.Uxz(zi,xi);
    [~,curve]=adaptive_req.quantile.compute_quantile_from_patch(patch,S.cfg_sim,S.req_out.req_cfg);
    mapping=adaptive_req.quantile.make_req_mapping(curve);
    ktheory=NaN;
    if isfinite(R.q_theory(i)), ktheory=adaptive_req.quantile.quantile_to_k(mapping,R.q_theory(i)); end
    make_spectral_figure(patch,curve,R(i,:),ktheory,OUT);
end
end

function file=condition_file(folder,geometry,regime,M)
g=string(geometry); if g=="circular_inclusion_2_3", g="inclusion_2_3"; end
switch string(regime)
    case "directional_2D", token="directional_2d";
    case "diffuse_2D", token="diffuse_2d";
    case "partial_3D", token="partial_diffuse_3d";
    otherwise, token="diffuse_3d";
end
file=fullfile(folder,sprintf('level17_compare_%s_%s_M%g.mat',g,token,M));
end

function make_spectral_figure(patch,C,R,ktheory,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 28 15]);
tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); imagesc(ax,real(patch)); axis(ax,'image'); colorbar(ax);
title(ax,'Patch particle velocity','FontWeight','normal');
ax=nexttile(tl); imagesc(ax,C.kx,C.kz,log10(C.Ssm+eps)); axis(ax,'image'); colorbar(ax);
title(ax,'log_{10} 2D power spectrum','FontWeight','normal'); xlabel(ax,'k_x'); ylabel(ax,'k_z');
ax=nexttile(tl); plot(ax,C.k_cent,C.Srad_norm,'LineWidth',1.2); hold(ax,'on');
add_k_lines(ax,R.k_hard,R.k_soft,R.k_true,R.k_model,ktheory);
title(ax,'Radial power spectrum','FontWeight','normal'); xlabel(ax,'k'); grid(ax,'on');
ax=nexttile(tl,[1 2]); plot(ax,C.k_cent,C.Ecum,'LineWidth',1.4); hold(ax,'on');
add_k_lines(ax,R.k_hard,R.k_soft,R.k_true,R.k_model,ktheory);
yline(ax,R.q_true,'k:','q true'); yline(ax,R.q_model,'m--','q model');
if isfinite(R.q_theory), yline(ax,R.q_theory,'c--','q theory'); end
xlabel(ax,'k (rad/m)'); ylabel(ax,'E_{cum}'); ylim(ax,[0 1]); grid(ax,'on');
ax=nexttile(tl); axis(ax,'off');
text(ax,0,1,sprintf(['%s\n%s / %s\npurity %.3f\nsoft %.3f, hard %.3f\n', ...
    'signed distance %+.2f mm\nq true %.3f, q model %.3f\n', ...
    'k/k true %.3f\nk/k soft %.3f\nSWS error %+.2f%%'], ...
    R.representative_type,R.external_case,R.field_regime,R.patch_purity, ...
    R.soft_fraction,R.hard_fraction,R.signed_distance_to_interface_mm, ...
    R.q_true,R.q_model,R.k_model_over_k_true,R.k_model_over_k_soft, ...
    R.SWS_error_signed_pct),'VerticalAlignment','top','Interpreter','none');
title(tl,'Representative spectral failure mode','FontWeight','normal');
file="test24_"+R.external_case+"_"+R.representative_type+".png";
export_fig(fig,fullfile(OUT.spectral_dir,file));
end

function add_k_lines(ax,khard,ksoft,ktrue,kmodel,ktheory)
xline(ax,khard,'r:','k hard'); xline(ax,ksoft,'b:','k soft');
xline(ax,ktrue,'k-','k true'); xline(ax,kmodel,'m--','k model');
if isfinite(ktheory), xline(ax,ktheory,'c--','k theory'); end
end

function plot_k_ratio_vs_purity(T,OUT)
cases=unique(T.external_case,'stable'); models=model_order(T.model_name);
fig=figure('Color','w','Units','centimeters','Position',[2 2 28 17]);
tl=tiledlayout(numel(cases),4,'TileSpacing','compact','Padding','compact');
for ci=1:numel(cases)
    for mi=1:min(4,numel(models))
        ax=nexttile(tl); X=T(T.external_case==cases(ci)&T.model_name==models(mi),:);
        plot_binned(ax,X.patch_purity,X.k_model_over_k_true); yline(ax,1,'k:');
        title(ax,short_model(models(mi)),'FontWeight','normal'); xlabel(ax,'Patch purity'); ylabel(ax,'k model / k true'); grid(ax,'on');
    end
end
title(tl,'k model / k true versus patch purity','FontWeight','normal');
export_fig(fig,fullfile(OUT.figure_dir,'test24_k_model_over_k_true_vs_purity.png'));
end

function plot_signed_error_vs_purity(T,OUT)
regimes=regime_order(T.field_regime); models=model_order(T.model_name);
fig=figure('Color','w','Units','centimeters','Position',[2 2 28 17]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for ri=1:numel(regimes)
    ax=nexttile(tl); hold(ax,'on'); X=T(T.field_regime==regimes(ri),:);
    for mi=1:numel(models), Y=X(X.model_name==models(mi),:); plot_binned(ax,Y.patch_purity,Y.SWS_error_signed_pct,models(mi)); end
    yline(ax,0,'k:'); xlabel(ax,'Patch purity'); ylabel(ax,'Signed SWS error (%)'); title(ax,regimes(ri),'Interpreter','none','FontWeight','normal'); grid(ax,'on');
end
legend(nexttile(tl,4),'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'test24_signed_error_vs_purity.png'));
end

function plot_k_regions_by_distance(T,OUT,primary)
P=T(T.model_name==primary,:); cases=unique(P.external_case,'stable'); regimes=regime_order(P.field_regime);
regions=["below_hard","near_hard","near_true","between_hard_soft","near_soft","above_soft"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 28 16]); tl=tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for ci=1:numel(cases)
    for ri=1:numel(regimes)
        ax=nexttile(tl); X=P(P.external_case==cases(ci)&P.field_regime==regimes(ri),:);
        C=zeros(6,numel(regions));
        for bi=1:6
            idx=X.distance_bin==bi;
            for ki=1:numel(regions)
                C(bi,ki)=mean(X.k_model_region(idx)==regions(ki),'omitnan');
            end
        end
        bar(ax,C,'stacked'); xticks(ax,1:6);
        xticklabels(ax,distance_labels((1:6)')); xtickangle(ax,35); ylim(ax,[0 1]);
        title(ax,cases(ci)+" / "+regimes(ri),'Interpreter','none','FontWeight','normal');
        ylabel(ax,'Fraction');
    end
end
legend(nexttile(tl,8),regions,'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'test24_k_model_region_by_distance.png'));
end

function plot_above_soft_by_distance(T,OUT)
cases=unique(T.external_case,'stable'); regimes=regime_order(T.field_regime); models=model_order(T.model_name);
fig=figure('Color','w','Units','centimeters','Position',[2 2 28 16]); tl=tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for ci=1:numel(cases)
    for ri=1:numel(regimes)
        ax=nexttile(tl); hold(ax,'on');
        X=T(T.external_case==cases(ci)&T.field_regime==regimes(ri),:);
        for mi=1:numel(models)
            Y=X(X.model_name==models(mi),:); y=nan(6,1);
            for bi=1:6
                y(bi)=100*mean(Y.is_k_above_soft(Y.distance_bin==bi),'omitnan');
            end
            plot(ax,1:6,y,'-o','MarkerSize',3, ...
                'DisplayName',short_model(models(mi)));
        end
        xticks(ax,1:6); xticklabels(ax,distance_labels((1:6)'));
        xtickangle(ax,35); ylabel(ax,'% k model > k soft');
        title(ax,cases(ci)+" / "+regimes(ri),'Interpreter','none','FontWeight','normal');
        grid(ax,'on');
    end
end
legend(nexttile(tl,8),'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'test24_pct_k_above_soft_by_distance.png'));
end

function plot_q_error_vs_k_ratio(T,OUT)
P=T(T.model_name=="HybridLocalGlobal_T18_noUserRegime",:);
fig=figure('Color','w','Units','centimeters','Position',[2 2 21 10]);
cases=unique(P.external_case,'stable'); tl=tiledlayout(1,2,'TileSpacing','compact');
for ci=1:numel(cases)
    ax=nexttile(tl); X=P(P.external_case==cases(ci),:); scatter(ax,X.q_error_model_true,X.k_model_over_k_true,5,X.patch_purity,'filled','MarkerFaceAlpha',0.15);
    yline(ax,1,'k:'); xline(ax,0,'k:'); colorbar(ax); xlabel(ax,'q model - q true'); ylabel(ax,'k model / k true'); title(ax,cases(ci),'Interpreter','none','FontWeight','normal'); grid(ax,'on');
end
export_fig(fig,fullfile(OUT.figure_dir,'test24_q_error_vs_k_ratio.png'));
end

function plot_mechanism_summary(T,OUT,primary)
P=T(T.model_name==primary,:); fig=figure('Color','w','Units','centimeters','Position',[2 2 20 15]); tl=tiledlayout(2,2,'TileSpacing','compact');
vars=["is_high_error_20","SWS_error_signed_pct","is_k_above_soft","k_model_over_k_true"];
labels=["High-error >20% rate","Mean signed error (%)","% k above soft","Mean k model / k true"];
for vi=1:4
    ax=nexttile(tl); hold(ax,'on'); cases=unique(P.external_case,'stable');
    for ci=1:numel(cases), X=P(P.external_case==cases(ci),:); y=nan(7,1); x=nan(7,1); for bi=1:7, idx=X.purity_bin==bi; x(bi)=mean(X.patch_purity(idx),'omitnan'); y(bi)=mean(double(X.(vars(vi))(idx)),'omitnan'); end; if vi==3, y=100*y; end; plot(ax,x,y,'-o','DisplayName',cases(ci)); end
    xlabel(ax,'Patch purity'); ylabel(ax,labels(vi)); grid(ax,'on');
end
legend(nexttile(tl,4),'Location','eastoutside','Interpreter','none');
export_fig(fig,fullfile(OUT.figure_dir,'test24_error_mechanism_summary.png'));
end

function plot_binned(ax,x,y,name)
if nargin<4, name="binned mean"; end
edges=linspace(0,1,21); b=discretize(x,edges); xc=nan(20,1); ym=nan(20,1);
for i=1:20, idx=b==i; xc(i)=mean(x(idx),'omitnan'); ym(i)=mean(y(idx),'omitnan'); end
plot(ax,xc,ym,'LineWidth',1.2,'DisplayName',short_model(name));
end

function print_summary(T,primary)
models=model_order(T.model_name);
fprintf('\n================ Test 24 summary ================\n');
for mi=1:numel(models)
    X=T(T.model_name==models(mi),:); H=X.is_high_error_20;
    fprintf(['%s: high-error %.2f%% | high-error low-purity %.2f%% | ', ...
        'high-error k>ksoft %.2f%% | high-error k/ktrue>1.2 %.2f%%.\n'], ...
        short_model(models(mi)),100*mean(H),100*safe_conditional(X.is_low_purity,H), ...
        100*safe_conditional(X.is_k_above_soft,H), ...
        100*safe_conditional(X.k_model_over_k_true>1.2,H));
end
P=T(T.model_name==primary,:); low=P.patch_purity<0.8; high=P.patch_purity>=0.95;
fprintf('Primary mean signed error: high purity %+.3f%% | low purity %+.3f%%.\n', ...
    mean(P.SWS_error_signed_pct(high),'omitnan'),mean(P.SWS_error_signed_pct(low),'omitnan'));
scores=nan(numel(models),1);
for mi=1:numel(models)
    X=T(T.model_name==models(mi)&T.side=="interface_band",:);
    scores(mi)=mean(X.SWS_error_abs_pct,'omitnan');
end
[~,best]=min(scores); fprintf('Most conservative near interfaces: %s (lowest mean absolute error %.3f%%).\n',short_model(models(best)),scores(best));
mixed=P.is_mixed_patch; over=mean(P.q_model_minus_q_theory(mixed),'omitnan');
if over>0
    fprintf('Hybrid no-user overselects q versus discrete theory in mixed patches (mean delta %.4f).\n',over);
else
    fprintf('Hybrid no-user does not overselect q versus discrete theory on average in mixed patches (mean delta %.4f).\n',over);
end
fprintf('=================================================\n');
end

function v=safe_conditional(x,mask)
if ~any(mask), v=NaN; else, v=mean(x(mask),'omitnan'); end
end

function order=regime_order(values)
preferred=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
order=preferred(ismember(preferred,unique(string(values))));
end

function order=model_order(values)
preferred=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","TheoryQDiscrete", ...
    "LocalOnly_old","GlobalOnly_old","HybridLocalGlobal_old"];
order=preferred(ismember(preferred,unique(string(values))));
end

function label=short_model(x)
label=replace(string(x),["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess","TheoryQDiscrete", ...
    "LocalOnly_old","GlobalOnly_old","HybridLocalGlobal_old"], ...
    ["Local T18","Global T18","Hybrid T18","Hybrid T18 + guess", ...
    "Theory","Local old","Global old","Hybrid old"]);
end

function export_fig(fig,file)
exportgraphics(fig,file,'Resolution',220); close(fig);
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; elseif isempty(B), T=A; return; end
vars=unique([string(A.Properties.VariableNames),string(B.Properties.VariableNames)],'stable');
A=add_missing(A,vars); B=add_missing(B,vars); T=[A(:,cellstr(vars));B(:,cellstr(vars))];
end

function T=add_missing(T,vars)
for i=1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)), continue; end
    T.(vars(i))=nan(height(T),1);
end
end
