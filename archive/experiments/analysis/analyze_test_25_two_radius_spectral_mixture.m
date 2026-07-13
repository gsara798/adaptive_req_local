%% analyze_test_25_two_radius_spectral_mixture.m
% Test 25: quantitative two-radius spectral-mixture diagnosis.
%
% Measures radial energy near the expected hard/soft wavenumbers in pure and
% mixed heterogeneous patches, then normalizes the apparent opposite-material
% energy against matched homogeneous fields. No model is trained or modified.
% Frozen Test 23 predictions are used only for retrospective error association.

clear; clc; close all;
format compact;

%% Setup

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
assert(exist(SOURCE,'file')==2,'Test 23 cache missing: %s',SOURCE);
assert(exist(CONDITION_DIR,'dir')==7,'Test 17 conditions missing: %s',CONDITION_DIR);

OUT=make_output_dirs(root_dir);
PRIMARY_MODEL="HybridLocalGlobal_T18_noUserRegime";
MAX_PER_STRATUM=30;
MAX_HOMOGENEOUS_PER_CONDITION=60;
RANDOM_SEED=25001;

fprintf('\nTest 25: quantitative two-radius spectral mixture\n');
fprintf('No training or model inference. Stratified spectral sample only.\n');

%% Frozen heterogeneous patch metadata

S=load(SOURCE,'T_joined'); P=S.T_joined; clear S;
P=P(P.model_name==PRIMARY_MODEL,:);
P.external_case=P.geometry_case;
P.material_side=repmat("soft",height(P),1);
P.material_side(P.sws_true>2.5)="hard";
P.purity_class=purity_class(P.patch_purity);
P.sws_signed_error_pct=100*(P.sws_pred-P.sws_true)./P.sws_true;
P.sws_abs_error_pct=abs(P.sws_signed_error_pct);
P.k_model_over_k_true=P.sws_true./P.sws_pred;

%% Extract spectra by condition with resumable checkpoints

T_spectral=table();
conditions=unique(P(:,{'geometry_case','field_regime','REQ_M'}),'rows','stable');
for ci=1:height(conditions)
    C=conditions(ci,:);
    key=C.geometry_case+"__"+C.field_regime+"__M"+string(C.REQ_M);
    checkpoint=fullfile(OUT.condition_dir,"heterogeneous__"+sanitize(key)+".mat");
    if exist(checkpoint,'file')==2
        S=load(checkpoint,'T_condition'); T_condition=S.T_condition;
        T_condition.external_case=repmat(C.geometry_case,height(T_condition),1);
        T_condition.field_regime=repmat(C.field_regime,height(T_condition),1);
        T_condition.REQ_M=C.REQ_M*ones(height(T_condition),1);
        save(checkpoint,'T_condition','-v7.3');
        T_spectral=concat_tables(T_spectral,T_condition);
        fprintf('Reused heterogeneous %s.\n',key); continue;
    end
    X=P(P.geometry_case==C.geometry_case & P.field_regime==C.field_regime & ...
        P.REQ_M==C.REQ_M,:);
    X=sample_stratified(X,MAX_PER_STRATUM,RANDOM_SEED+ci);
    file=condition_file(CONDITION_DIR,C.geometry_case,C.field_regime,C.REQ_M);
    T_condition=extract_condition_spectra(X,file,"heterogeneous");
    save(checkpoint,'T_condition','-v7.3');
    T_spectral=concat_tables(T_spectral,T_condition);
    fprintf('Heterogeneous %s: %d spectra.\n',key,height(T_condition));
end

%% Matched homogeneous spectral baselines

regimes=regime_order(P.field_regime); M_values=unique(P.REQ_M,'stable');
for cs=[2 3]
    for ri=1:numel(regimes)
        for mi=1:numel(M_values)
            M=M_values(mi); case_name="homogeneous_cs"+string(cs);
            key=case_name+"__"+regimes(ri)+"__M"+string(M);
            checkpoint=fullfile(OUT.condition_dir,"reference__"+sanitize(key)+".mat");
            if exist(checkpoint,'file')==2
                S=load(checkpoint,'T_condition'); T_spectral=concat_tables(T_spectral,S.T_condition);
                fprintf('Reused reference %s.\n',key); continue;
            end
            file=condition_file(CONDITION_DIR,case_name,regimes(ri),M);
            X=make_homogeneous_sample(file,case_name,regimes(ri),M,cs, ...
                MAX_HOMOGENEOUS_PER_CONDITION,RANDOM_SEED+1000*cs+100*ri+mi);
            T_condition=extract_condition_spectra(X,file,"homogeneous_reference");
            save(checkpoint,'T_condition','-v7.3');
            T_spectral=concat_tables(T_spectral,T_condition);
            fprintf('Reference %s: %d spectra.\n',key,height(T_condition));
        end
    end
end

%% Normalize against matched homogeneous material baselines

[T_reference,T_spectral]=add_homogeneous_normalization(T_spectral);
assert_spectral_policy(T_spectral);

%% Summaries

T_hetero=T_spectral(T_spectral.sample_type=="heterogeneous",:);
T_summary=summarize_spectral(T_hetero,["external_case","field_regime", ...
    "REQ_M","material_side","purity_class"]);
T_summary_overall=summarize_spectral(T_hetero,["external_case", ...
    "material_side","purity_class"]);
T_reference_summary=summarize_spectral(T_reference, ...
    ["external_case","field_regime","REQ_M","material_side"]);

writetable(T_spectral,fullfile(OUT.table_dir,'test25_patch_spectral_mixture_metrics.csv'));
writetable(T_reference_summary,fullfile(OUT.table_dir,'test25_homogeneous_spectral_baseline.csv'));
writetable(T_summary,fullfile(OUT.table_dir,'test25_spectral_mixture_summary.csv'));
writetable(T_summary_overall,fullfile(OUT.table_dir,'test25_spectral_mixture_summary_overall.csv'));
save(fullfile(OUT.data_dir,'test25_spectral_mixture_results.mat'), ...
    'T_spectral','T_reference','T_summary','T_summary_overall','-v7.3');

%% Diagnostic figures

plot_mixture_vs_purity(T_spectral,OUT);
plot_excess_energy_vs_purity(T_spectral,OUT);
plot_two_radius_rate(T_spectral,OUT);
plot_mixture_vs_error(T_hetero,OUT);
plot_distance_vs_mixture(T_hetero,OUT);
plot_pure_vs_homogeneous(T_spectral,OUT);

T_representative=select_representative_spectra(T_spectral);
writetable(T_representative,fullfile(OUT.table_dir, ...
    'test25_representative_spectral_patch_list.csv'));
plot_representative_spectra(T_representative,CONDITION_DIR,OUT);

%% Interpretation

print_summary(T_spectral);
fprintf('\nTables: %s\nFigures: %s\n',OUT.table_dir,OUT.figure_dir);
fprintf('Test 25 complete: %s\n',OUT.root_dir);

%% Local functions

function OUT=make_output_dirs(root_dir)
OUT.root_dir=fullfile(root_dir,'outputs','test_25_two_radius_spectral_mixture');
OUT.table_dir=fullfile(OUT.root_dir,'tables');
OUT.figure_dir=fullfile(OUT.root_dir,'figures');
OUT.spectral_dir=fullfile(OUT.figure_dir,'representative_spectra');
OUT.data_dir=fullfile(OUT.root_dir,'data');
OUT.condition_dir=fullfile(OUT.data_dir,'condition_checkpoints');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs), if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end, end
end

function labels=purity_class(p)
labels=repmat("strongly_mixed",size(p));
labels(p>=0.80)="moderately_mixed";
labels(p>=0.95)="near_pure";
labels(p>=0.99)="pure";
end

function X=sample_stratified(X,max_n,seed)
rng(seed,'twister'); classes=unique(X.purity_class,'stable'); sides=unique(X.material_side,'stable');
keep=false(height(X),1);
for ci=1:numel(classes)
    for si=1:numel(sides)
        idx=find(X.purity_class==classes(ci)&X.material_side==sides(si));
        if isempty(idx), continue; end
        n=min(max_n,numel(idx)); idx=idx(randperm(numel(idx),n)); keep(idx)=true;
    end
end
X=X(keep,:);
end

function X=make_homogeneous_sample(file,case_name,regime,M,cs,nmax,seed)
S=load(file,'sim','cfg_sim'); cfg=S.cfg_sim;
[~,feat]=req_settings(cfg,M); half=floor(feat.win_size/2);
x=(1+half):(size(S.sim.Uxz,2)-half); z=(1+half):(size(S.sim.Uxz,1)-half);
[CX,CZ]=meshgrid(x,z); cx=CX(:); cz=CZ(:); rng(seed,'twister');
n=min(nmax,numel(cx)); take=randperm(numel(cx),n); cx=cx(take); cz=cz(take);
X=table(); X.external_case=repmat(case_name,n,1); X.field_regime=repmat(regime,n,1);
X.REQ_M=M*ones(n,1); X.frequency_hz=cfg.f0*ones(n,1); X.dx=cfg.dx*ones(n,1);
X.cx=cx; X.cz=cz; X.map_iz=(1:n)'; X.map_ix=ones(n,1);
X.material_side=repmat(ternary(cs==2,"soft","hard"),n,1);
X.patch_purity=ones(n,1); X.purity_class=repmat("homogeneous_reference",n,1);
X.sws_signed_error_pct=nan(n,1); X.sws_abs_error_pct=nan(n,1);
X.k_model_over_k_true=nan(n,1); X.abs_distance_to_interface_mm=nan(n,1);
X.signed_distance_mm=nan(n,1); X.condition_id=zeros(n,1);
end

function T=extract_condition_spectra(X,file,sample_type)
S=load(file,'sim','cfg_sim'); cfg=S.cfg_sim; M=X.REQ_M(1);
[req,feat]=req_settings(cfg,M); half=floor(feat.win_size/2);
n=height(X); rows=cell(n,1);
for i=1:n
    zi=(X.cz(i)-half):(X.cz(i)+half); xi=(X.cx(i)-half):(X.cx(i)+half);
    patch=S.sim.Uxz(zi,xi);
    [~,C]=adaptive_req.quantile.compute_quantile_from_patch(patch,cfg,req);
    rows{i}=spectral_metrics(C,cfg.f0);
end
Mtab=struct2table([rows{:}]');
vars=["condition_id","external_case","field_regime","REQ_M", ...
    "frequency_hz","dx","cx","cz","map_iz","map_ix", ...
    "material_side","patch_purity","purity_class", ...
    "sws_signed_error_pct","sws_abs_error_pct","k_model_over_k_true", ...
    "abs_distance_to_interface_mm","signed_distance_mm"];
vars=vars(ismember(vars,string(X.Properties.VariableNames)));
T=[X(:,cellstr(vars)) Mtab]; T.sample_type=repmat(string(sample_type),n,1);
T=movevars(T,'sample_type','Before',1);
end

function [req,feat]=req_settings(cfg,M)
feat=adaptive_req.config.default_feature_config('M',M,'cs_guess',3, ...
    'gamma_win',1,'pad_factor',1);
[req,feat]=adaptive_req.config.default_req_config(cfg,feat, ...
    'Nbins','auto','Nbins_auto_oversample',1,'Nbins_min',16,'smooth_sigma',1);
end

function R=spectral_metrics(C,f0)
ksoft=2*pi*f0/2; khard=2*pi*f0/3; gap=ksoft-khard;
k=C.k_cent(:); s=max(double(C.Srad(:)),0); s(~isfinite(s))=0;
if sum(s)>0, s=s/sum(s); end
dk=median(diff(k),'omitnan'); half_band=min(0.45*gap,max(0.10*gap,1.5*dk));
hard_band=abs(k-khard)<=half_band; soft_band=abs(k-ksoft)<=half_band;
Ehard=sum(s(hard_band)); Esoft=sum(s(soft_band));
mix=2*min(Ehard,Esoft)/(Ehard+Esoft+eps);
peaks=find(s(2:end-1)>=s(1:end-2)&s(2:end-1)>s(3:end)& ...
    s(2:end-1)>=0.08*max(s))+1;
hard_peaks=peaks(abs(k(peaks)-khard)<=half_band);
soft_peaks=peaks(abs(k(peaks)-ksoft)<=half_band);
R=struct(); R.k_hard=khard; R.k_soft=ksoft; R.band_half_width=half_band;
R.hard_band_energy=Ehard; R.soft_band_energy=Esoft;
R.spectral_mixture_index=mix;
R.has_hard_radius_peak=~isempty(hard_peaks); R.has_soft_radius_peak=~isempty(soft_peaks);
R.two_radius_detected=R.has_hard_radius_peak&&R.has_soft_radius_peak;
R.num_radial_peaks=numel(peaks);
R.hard_peak_k=nearest_peak(k,s,hard_peaks); R.soft_peak_k=nearest_peak(k,s,soft_peaks);
R.radial_width_10_90=quantile_k(C,0.90)-quantile_k(C,0.10);
R.radial_width_25_75=quantile_k(C,0.75)-quantile_k(C,0.25);
R.radial_entropy=-sum(s(s>0).*log(s(s>0)))/log(max(numel(s),2));
end

function kp=nearest_peak(k,s,idx)
if isempty(idx), kp=NaN; else, [~,j]=max(s(idx)); kp=k(idx(j)); end
end

function kv=quantile_k(C,q)
E=C.Ecum(:); k=C.k_cent(:); ok=isfinite(E)&isfinite(k);
if nnz(ok)<2, kv=NaN; else, [Eu,ia]=unique(E(ok),'stable'); ku=k(ok); ku=ku(ia); kv=interp1(Eu,ku,q,'linear','extrap'); end
end

function [R,T]=add_homogeneous_normalization(T)
R0=T(T.sample_type=="homogeneous_reference",:);
groups=["field_regime","REQ_M","material_side"];
[G,B]=findgroups(R0(:,cellstr(groups)));
B.baseline_opposite_energy=splitapply(@(e,s)median(opposite(e,s),'omitnan'), ...
    [R0.hard_band_energy R0.soft_band_energy],R0.material_side,G);
B.baseline_mixture_index=splitapply(@(x)median(x,'omitnan'),R0.spectral_mixture_index,G);
B.baseline_two_radius_rate=splitapply(@mean,R0.two_radius_detected,G);
key=group_key(T,groups); keyb=group_key(B,groups); [tf,loc]=ismember(key,keyb);
assert(all(tf),'Homogeneous spectral baseline join failed.');
opp=opposite([T.hard_band_energy T.soft_band_energy],T.material_side);
T.opposite_material_energy=opp;
T.baseline_opposite_energy=B.baseline_opposite_energy(loc);
T.excess_opposite_energy=opp-T.baseline_opposite_energy;
T.opposite_energy_ratio_to_homogeneous=opp./max(T.baseline_opposite_energy,eps);
T.baseline_mixture_index=B.baseline_mixture_index(loc);
T.excess_spectral_mixture_index=T.spectral_mixture_index-T.baseline_mixture_index;
R=T(T.sample_type=="homogeneous_reference",:);
end

function v=opposite(energy,side)
v=energy(:,2); v(side=="soft")=energy(side=="soft",1);
end

function assert_spectral_policy(T)
assert(all(T.spectral_mixture_index>=0&T.spectral_mixture_index<=1+eps));
fprintf(['Diagnostic policy: material labels and errors were attached only ', ...
    'after frozen fields/predictions; no predictor or model was changed.\n']);
end

function S=summarize_spectral(T,groups)
[G,S]=findgroups(T(:,cellstr(groups)));
S.N=splitapply(@numel,T.spectral_mixture_index,G);
S.mean_mixture_index=splitapply(@(x)mean(x,'omitnan'),T.spectral_mixture_index,G);
S.median_mixture_index=splitapply(@(x)median(x,'omitnan'),T.spectral_mixture_index,G);
S.mean_excess_mixture=splitapply(@(x)mean(x,'omitnan'),T.excess_spectral_mixture_index,G);
S.mean_opposite_energy=splitapply(@(x)mean(x,'omitnan'),T.opposite_material_energy,G);
S.mean_excess_opposite_energy=splitapply(@(x)mean(x,'omitnan'),T.excess_opposite_energy,G);
S.two_radius_detected_pct=100*splitapply(@mean,T.two_radius_detected,G);
S.mean_radial_width_10_90=splitapply(@(x)mean(x,'omitnan'),T.radial_width_10_90,G);
S.mean_abs_sws_error_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_abs_error_pct,G);
S.mean_signed_sws_error_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_signed_error_pct,G);
end

function plot_mixture_vs_purity(T,OUT)
H=T(T.sample_type=="heterogeneous",:); cases=unique(H.external_case,'stable'); sides=["soft","hard"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 22 16]); tl=tiledlayout(2,2,'TileSpacing','compact');
for ci=1:2
    for si=1:2
        ax=nexttile(tl);
        X=H(H.external_case==cases(ci)&H.material_side==sides(si),:);
        scatter(ax,X.patch_purity,X.spectral_mixture_index,10,X.REQ_M, ...
            'filled','MarkerFaceAlpha',0.25); hold(ax,'on');
        plot_binned(ax,X.patch_purity,X.spectral_mixture_index);
        yline(ax,median(X.baseline_mixture_index,'omitnan'), ...
            'k:','homogeneous median');
        xlabel(ax,'Patch purity'); ylabel(ax,'Spectral mixture index');
        title(ax,cases(ci)+" / "+sides(si),'Interpreter','none','FontWeight','normal');
        colorbar(ax); grid(ax,'on');
    end
end
export_fig(fig,fullfile(OUT.figure_dir,'test25_spectral_mixture_index_vs_purity.png'));
end

function plot_excess_energy_vs_purity(T,OUT)
H=T(T.sample_type=="heterogeneous",:); cases=unique(H.external_case,'stable');
fig=figure('Color','w','Units','centimeters','Position',[2 2 21 10]); tl=tiledlayout(1,2,'TileSpacing','compact');
for ci=1:2
    ax=nexttile(tl); X=H(H.external_case==cases(ci),:); scatter(ax,X.patch_purity,X.excess_opposite_energy,10,X.sws_abs_error_pct,'filled','MarkerFaceAlpha',0.25); hold(ax,'on');
    plot_binned(ax,X.patch_purity,X.excess_opposite_energy); yline(ax,0,'k:'); xlabel(ax,'Patch purity'); ylabel(ax,'Opposite-band energy minus homogeneous baseline'); title(ax,cases(ci),'Interpreter','none','FontWeight','normal'); colorbar(ax); grid(ax,'on');
end
export_fig(fig,fullfile(OUT.figure_dir,'test25_excess_opposite_energy_vs_purity.png'));
end

function plot_two_radius_rate(T,OUT)
H=T(T.sample_type=="heterogeneous",:); classes=["pure","near_pure","moderately_mixed","strongly_mixed"];
cases=unique(H.external_case,'stable'); fig=figure('Color','w','Units','centimeters','Position',[2 2 20 11]); ax=axes(fig);
C=nan(numel(classes),numel(cases));
for ci=1:numel(cases), for pi=1:numel(classes), idx=H.external_case==cases(ci)&H.purity_class==classes(pi); C(pi,ci)=100*mean(H.two_radius_detected(idx),'omitnan'); end, end
bar(ax,C); xticks(ax,1:numel(classes)); xticklabels(ax,classes); xtickangle(ax,25); ylabel(ax,'Two-radius detection (%)'); legend(ax,cases,'Interpreter','none'); grid(ax,'on');
export_fig(fig,fullfile(OUT.figure_dir,'test25_two_radius_detection_by_purity.png'));
end

function plot_mixture_vs_error(H,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 21 10]); cases=unique(H.external_case,'stable'); tl=tiledlayout(1,2,'TileSpacing','compact');
for ci=1:2
    ax=nexttile(tl); X=H(H.external_case==cases(ci),:); scatter(ax,X.excess_spectral_mixture_index,X.sws_signed_error_pct,10,X.patch_purity,'filled','MarkerFaceAlpha',0.3); yline(ax,0,'k:'); xlabel(ax,'Excess mixture index vs homogeneous'); ylabel(ax,'Signed SWS error (%)'); title(ax,cases(ci),'Interpreter','none','FontWeight','normal'); colorbar(ax); grid(ax,'on');
end
export_fig(fig,fullfile(OUT.figure_dir,'test25_spectral_mixture_vs_sws_error.png'));
end

function plot_distance_vs_mixture(H,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 21 10]); cases=unique(H.external_case,'stable'); tl=tiledlayout(1,2,'TileSpacing','compact');
for ci=1:2
    ax=nexttile(tl); X=H(H.external_case==cases(ci),:); scatter(ax,X.abs_distance_to_interface_mm,X.excess_opposite_energy,10,X.patch_purity,'filled','MarkerFaceAlpha',0.3); yline(ax,0,'k:'); xlabel(ax,'Absolute distance to interface (mm)'); ylabel(ax,'Excess opposite-material energy'); title(ax,cases(ci),'Interpreter','none','FontWeight','normal'); colorbar(ax); grid(ax,'on');
end
export_fig(fig,fullfile(OUT.figure_dir,'test25_excess_opposite_energy_vs_distance.png'));
end

function plot_pure_vs_homogeneous(T,OUT)
classes=["homogeneous_reference","pure","near_pure","moderately_mixed","strongly_mixed"];
sides=["soft","hard"]; fig=figure('Color','w','Units','centimeters','Position',[2 2 22 11]); tl=tiledlayout(1,2,'TileSpacing','compact');
for si=1:2
    ax=nexttile(tl); data=[]; group=strings(0,1);
    for ci=1:numel(classes)
        if classes(ci)=="homogeneous_reference", idx=T.sample_type=="homogeneous_reference"&T.material_side==sides(si); else, idx=T.sample_type=="heterogeneous"&T.material_side==sides(si)&T.purity_class==classes(ci); end
        data=[data;T.excess_opposite_energy(idx)]; group=[group;repmat(classes(ci),sum(idx),1)]; %#ok<AGROW>
    end
    boxchart(ax,categorical(group,classes),data); yline(ax,0,'k:');
    ylabel(ax,'Excess opposite-band energy vs matched homogeneous');
    title(ax,sides(si),'FontWeight','normal'); grid(ax,'on');
end
export_fig(fig,fullfile(OUT.figure_dir,'test25_pure_heterogeneous_vs_homogeneous.png'));
end

function R=select_representative_spectra(T)
R=table(); sides=["soft","hard"];
for si=1:2
    X=T(T.material_side==sides(si),:);
    specs=[ ...
        struct('name',"homogeneous_reference",'mask',X.sample_type=="homogeneous_reference",'score',X.spectral_mixture_index)
        struct('name',"pure_low_excess",'mask',X.sample_type=="heterogeneous"&X.purity_class=="pure",'score',abs(X.excess_opposite_energy))
        struct('name',"pure_high_excess",'mask',X.sample_type=="heterogeneous"&X.purity_class=="pure",'score',-X.excess_opposite_energy)
        struct('name',"mixed_high_mixture",'mask',X.sample_type=="heterogeneous"&X.purity_class=="strongly_mixed",'score',-X.spectral_mixture_index)];
    for qi=1:numel(specs)
        C=X(specs(qi).mask,:); if isempty(C), continue; end
        score=specs(qi).score(specs(qi).mask); [~,j]=min(score); row=C(j,:);
        row.representative_type=specs(qi).name; R=concat_tables(R,row);
    end
end
R=movevars(R,'representative_type','Before',1);
end

function plot_representative_spectra(R,folder,OUT)
for i=1:height(R)
    Ri=R(i,:);
    file=condition_file(folder,Ri.external_case,Ri.field_regime,Ri.REQ_M);
    S=load(file,'sim','cfg_sim'); [req,feat]=req_settings(S.cfg_sim,Ri.REQ_M); half=floor(feat.win_size/2);
    patch=S.sim.Uxz((Ri.cz-half):(Ri.cz+half),(Ri.cx-half):(Ri.cx+half));
    [~,C]=adaptive_req.quantile.compute_quantile_from_patch(patch,S.cfg_sim,req);
    fig=figure('Color','w','Units','centimeters','Position',[2 2 27 9]); tl=tiledlayout(1,3,'TileSpacing','compact');
    ax=nexttile(tl); imagesc(ax,real(patch)); axis(ax,'image'); colorbar(ax); title(ax,'Particle velocity');
    ax=nexttile(tl); imagesc(ax,C.kx,C.kz,log10(C.Ssm+eps)); axis(ax,'image'); colorbar(ax); title(ax,'2D power spectrum');
    ax=nexttile(tl); plot(ax,C.k_cent,C.Srad_norm,'LineWidth',1.3); xline(ax,Ri.k_hard,'r:','k hard'); xline(ax,Ri.k_soft,'b:','k soft'); title(ax,sprintf('%s | mix %.3f | excess %.3g',Ri.representative_type,Ri.spectral_mixture_index,Ri.excess_opposite_energy),'Interpreter','none','FontWeight','normal'); xlabel(ax,'k'); grid(ax,'on');
    title(tl,Ri.material_side+" / "+Ri.external_case+" / "+Ri.purity_class,'Interpreter','none');
    export_fig(fig,fullfile(OUT.spectral_dir,"test25_"+Ri.material_side+"_"+Ri.representative_type+".png"));
end
end

function print_summary(T)
H=T(T.sample_type=="heterogeneous",:); fprintf('\n================ Test 25 summary ================\n');
for side=["soft","hard"]
    X=H(H.material_side==side,:); base=T(T.sample_type=="homogeneous_reference"&T.material_side==side,:);
    fprintf('%s homogeneous opposite-energy median: %.4f.\n',side,median(base.opposite_material_energy,'omitnan'));
    for cls=["pure","near_pure","moderately_mixed","strongly_mixed"]
        Y=X(X.purity_class==cls,:); fprintf('  %s: N=%d, mixture %.3f, excess opposite %.4f, two-radius %.1f%%, signed error %+.2f%%.\n',cls,height(Y),mean(Y.spectral_mixture_index,'omitnan'),mean(Y.excess_opposite_energy,'omitnan'),100*mean(Y.two_radius_detected,'omitnan'),mean(Y.sws_signed_error_pct,'omitnan'));
    end
end
fprintf('=================================================\n');
end

function plot_binned(ax,x,y)
edges=linspace(0.5,1,16); b=discretize(x,edges); xc=nan(15,1); ym=nan(15,1);
for i=1:15, idx=b==i; xc(i)=mean(x(idx),'omitnan'); ym(i)=mean(y(idx),'omitnan'); end
plot(ax,xc,ym,'k-','LineWidth',1.4);
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
assert(exist(file,'file')==2,'Condition file missing: %s',file);
end

function key=group_key(T,vars)
key=repmat("",height(T),1);
for i=1:numel(vars), key=key+"|"+string(T.(vars(i))); end
end

function order=regime_order(values)
preferred=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
order=preferred(ismember(preferred,unique(string(values))));
end

function out=ternary(test,a,b)
if test, out=a; else, out=b; end
end

function s=sanitize(x)
s=regexprep(char(x),'[^A-Za-z0-9_-]','_');
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
