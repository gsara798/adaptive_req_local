%% analyze_test_23_interface_patch_contamination_step01.m
% Test 23, Step 01: interface-patch contamination diagnosis.
%
% This script trains and predicts nothing. It reads frozen external
% predictions from Test 20, adds oracle geometry variables only for
% retrospective diagnosis, and recomputes REQ spectra only for a small set
% of representative patches. Diagnostic variables never enter a model.
%
% Signed-distance convention:
%   bilayer:  x-interface_x; negative=soft, positive=hard.
%   inclusion: radius-distance_to_center; negative=soft background,
%              positive=hard inclusion.

% Test 20 is used because it contains all T18, old, and TheoryQDiscrete
% predictions on the same Test 17 fields. The underlying individual Test 17
% condition files provide the frozen wavefield for spectral illustrations.

clear; clc; close all;
format compact;

%% Setup and immutable diagnostic policy

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontName','Helvetica', ...
    'defaultTextFontName','Helvetica','defaultAxesFontSize',8, ...
    'defaultTextFontSize',8,'defaultLegendFontSize',7);

SOURCE=fullfile(root_dir,'outputs', ...
    'test_20_external_validation_and_aperture_q_tracking','analysis', ...
    'tables','level20_external_predictions.csv');
CONDITION_DIR=fullfile(root_dir,'outputs', ...
    'test_17_model_comparison_heterogeneous_cases','analysis','data','conditions');
assert(exist(SOURCE,'file')==2,'Test 20 predictions not found: %s',SOURCE);
assert(exist(CONDITION_DIR,'dir')==7,'Test 17 condition folder not found: %s',CONDITION_DIR);

OUT=make_output_dirs(root_dir);
PRIMARY_MODEL="HybridLocalGlobal_T18_noUserRegime";
MODEL_ORDER=["LocalOnly_T18","GlobalOnly_T18", ...
    "HybridLocalGlobal_T18_noUserRegime", ...
    "HybridLocalGlobal_T18_withUserRegimeGuess", ...
    "TheoryQDiscrete","LocalOnly_old","GlobalOnly_old", ...
    "HybridLocalGlobal_old"];
PURITY_EDGES=[0 0.5 0.75 0.9 0.95 1.000001];
DISTANCE_EDGES_MM=[-Inf -10 -5 -2 -1 0 1 2 5 10 Inf];

fprintf('\nTest 23, Step 01: interface-patch contamination diagnosis\n');
fprintf('No model training or model inference is performed.\n');
fprintf('Predictions: %s\n',SOURCE);
fprintf(['Diagnostic-only variables: signed distance, material mask, ', ...
    'patch purity, and error labels.\n']);

%% Load frozen predictions and add diagnostic geometry

cache_file=fullfile(OUT.data_dir,'test23_joined_predictions.mat');
if exist(cache_file,'file')==2
    S=load(cache_file,'T_joined'); T_joined=S.T_joined;
    fprintf('Reused joined diagnostic cache: %d rows.\n',height(T_joined));
else
    T=read_test20_heterogeneous_predictions(SOURCE,MODEL_ORDER);
    assert(~isempty(T),'No heterogeneous Test 20 predictions were loaded.');
    assert(all(ismember(["bilayer_2_3","inclusion_2_3"],unique(T.geometry_id))), ...
        'Both bilayer and inclusion cases are required.');
    T_joined=add_patch_geometry_diagnostics(T);
    clear T;
    save(cache_file,'T_joined','SOURCE','-v7.3');
end

assert_no_inference_leakage(T_joined);
assert_side_coverage(T_joined);
report_contamination_counts(T_joined,PRIMARY_MODEL);

%% Grouped diagnostic tables

base_groups=["geometry_case","field_regime","frequency_hz","dx_m", ...
    "dz_m","REQ_M","model_name"];
T_signed=summarize_errors(T_joined,[base_groups "side"]);
T_purity=summarize_errors(T_joined,[base_groups "side" "patch_is_mixed" ...
    "patch_is_strongly_mixed"]);

T_bins=T_joined;
T_bins.patch_purity_bin=discretize(T_bins.patch_purity,PURITY_EDGES, ...
    'IncludedEdge','left');
T_bins.patch_purity_bin_label=purity_bin_labels(T_bins.patch_purity_bin);
T_purity_bin=summarize_errors(T_bins, ...
    [base_groups "side" "patch_purity_bin" "patch_purity_bin_label"]);

T_bins.signed_distance_bin=discretize(T_bins.signed_distance_mm, ...
    DISTANCE_EDGES_MM,'IncludedEdge','right');
T_bins.signed_distance_bin_label=distance_bin_labels( ...
    T_bins.signed_distance_bin,DISTANCE_EDGES_MM);
T_distance_bin=summarize_errors(T_bins, ...
    [base_groups "side" "signed_distance_bin" "signed_distance_bin_label"]);

%% Representative patches and spectral diagnostics

T_representative=select_representative_patches(T_joined,PRIMARY_MODEL);
writetable(T_representative,fullfile(OUT.table_dir, ...
    'test23_representative_patch_list.csv'));
plot_representative_spectra(T_representative,T_joined,CONDITION_DIR,OUT);

%% Save requested tables

writetable(T_joined,fullfile(OUT.table_dir, ...
    'test23_patch_contamination_joined_predictions.csv'));
writetable(T_signed,fullfile(OUT.table_dir, ...
    'test23_signed_distance_summary.csv'));
writetable(T_purity,fullfile(OUT.table_dir, ...
    'test23_patch_purity_summary.csv'));
writetable(T_purity_bin,fullfile(OUT.table_dir, ...
    'test23_error_by_patch_purity_bin.csv'));
writetable(T_distance_bin,fullfile(OUT.table_dir, ...
    'test23_error_by_side_and_distance_bin.csv'));
save(fullfile(OUT.data_dir,'test23_analysis_results.mat'), ...
    'T_signed','T_purity','T_purity_bin','T_distance_bin', ...
    'T_representative','PURITY_EDGES','DISTANCE_EDGES_MM','-v7.3');

%% Main diagnostic figures

plot_distance_relationship(T_joined,OUT,'sws_signed_error_pct', ...
    'Signed SWS error (%)','signed', ...
    'test23_signed_sws_error_vs_signed_distance');
plot_distance_relationship(T_joined,OUT,'sws_abs_error_pct', ...
    'Absolute SWS error (%)','absolute', ...
    'test23_absolute_sws_error_vs_abs_distance');
plot_distance_relationship(T_joined,OUT,'q_signed_error', ...
    'Signed q error','signed','test23_q_error_vs_signed_distance');
plot_error_vs_purity(T_bins,OUT);
plot_high_error_vs_purity(T_bins,OUT);

%% Console interpretation

print_summary(T_joined,PRIMARY_MODEL);
fprintf('\nTest 23, Step 01 complete. Analysis folder:\n%s\n',OUT.analysis_dir);

%% Local functions

function OUT=make_output_dirs(root_dir)
OUT.analysis_dir=fullfile(root_dir,'outputs', ...
    'test_23_interface_patch_contamination','analysis', ...
    'level_23_interface_patch_contamination_step01');
OUT.table_dir=fullfile(OUT.analysis_dir,'tables');
OUT.figure_dir=fullfile(OUT.analysis_dir,'figures');
OUT.spectral_dir=fullfile(OUT.figure_dir,'representative_spectra');
OUT.data_dir=fullfile(OUT.analysis_dir,'data');
dirs=string(struct2cell(OUT));
for i=1:numel(dirs), if exist(dirs(i),'dir')~=7, mkdir(dirs(i)); end, end
end

function T=read_test20_heterogeneous_predictions(file,model_order)
keep=["condition_id","geometry_id","geometry_type", ...
    "field_regime_label","field_regime_variant","SIM_f0", ...
    "SIM_dx","SIM_dz","REQ_M","REQ_cs_guess","map_iz","map_ix", ...
    "cx","cz","x_center_m","z_center_m","model_name", ...
    "q_true","q_pred","cs_true","cs_pred"];
ds=tabularTextDatastore(file,'TextType','string');
available=string(ds.VariableNames);
missing=setdiff(keep,available);
assert(isempty(missing),'Test 20 table missing: %s',strjoin(missing,', '));
ds.SelectedVariableNames=cellstr(keep);
parts=cell(0,1); n=0;
while hasdata(ds)
    X=read(ds);
    idx=ismember(string(X.geometry_id),["bilayer_2_3","inclusion_2_3"]) & ...
        ismember(string(X.model_name),model_order);
    if any(idx), n=n+1; parts{n,1}=X(idx,:); end
end
assert(~isempty(parts),'No requested rows found in Test 20 CSV.');
T=vertcat(parts{:}); clear parts;
string_vars=["geometry_id","geometry_type","field_regime_label", ...
    "field_regime_variant","model_name"];
for i=1:numel(string_vars), T.(string_vars(i))=string(T.(string_vars(i))); end
T=sortrows(T,{'condition_id','map_iz','map_ix','model_name'});
end

function T=add_patch_geometry_diagnostics(T)
T.geometry_case=T.geometry_id;
T.geometry_case(T.geometry_id=="inclusion_2_3")="circular_inclusion_2_3";
T.field_regime=T.field_regime_label;
T.frequency_hz=T.SIM_f0; T.dx_m=T.SIM_dx; T.dz_m=T.SIM_dz;
T.sws_true=T.cs_true; T.sws_pred=T.cs_pred;
T.sws_signed_error_pct=100*(T.sws_pred-T.sws_true)./T.sws_true;
T.sws_abs_error_pct=abs(T.sws_signed_error_pct);
T.q_signed_error=T.q_pred-T.q_true; T.q_abs_error=abs(T.q_signed_error);
T.high_error_gt10=T.sws_abs_error_pct>10;
T.high_error_gt20=T.sws_abs_error_pct>20;

keys=T(:,{'condition_id','map_iz','map_ix','cx','cz','x_center_m', ...
    'z_center_m','geometry_id','SIM_f0','SIM_dx','SIM_dz','REQ_M', ...
    'REQ_cs_guess','cs_true'});
[~,ia]=unique(keys(:,{'condition_id','map_iz','map_ix'}),'rows','stable');
B=keys(ia,:);
B=compute_patch_purity(B);

key_all=string(T.condition_id)+"|"+string(T.map_iz)+"|"+string(T.map_ix);
key_base=string(B.condition_id)+"|"+string(B.map_iz)+"|"+string(B.map_ix);
[tf,loc]=ismember(key_all,key_base); assert(all(tf),'Patch diagnostic join failed.');
vars=["signed_distance_mm","abs_distance_to_interface_mm", ...
    "patch_win_size_px","patch_radius_mm","patch_soft_fraction", ...
    "patch_hard_fraction","patch_purity","patch_is_mixed", ...
    "patch_is_strongly_mixed","side","signed_distance_convention"];
for i=1:numel(vars)
    values=B.(vars(i));
    T.(vars(i))=values(loc,:);
end
end

function B=compute_patch_purity(B)
n=height(B);
B.signed_distance_mm=nan(n,1); B.abs_distance_to_interface_mm=nan(n,1);
B.patch_win_size_px=nan(n,1); B.patch_radius_mm=nan(n,1);
B.patch_soft_fraction=nan(n,1); B.patch_hard_fraction=nan(n,1);
B.patch_purity=nan(n,1); B.patch_is_mixed=false(n,1);
B.patch_is_strongly_mixed=false(n,1); B.side=strings(n,1);
B.signed_distance_convention=strings(n,1);
groups=unique(B.condition_id,'stable');
for gi=1:numel(groups)
    idx=find(B.condition_id==groups(gi)); X=B(idx,:);
    dx=X.SIM_dx(1); dz=X.SIM_dz(1); f0=X.SIM_f0(1);
    M=X.REQ_M(1); cs_guess=X.REQ_cs_guess(1);
    win=round(M*(cs_guess/f0)/dx); win=max(win,3);
    if mod(win,2)==0, win=win+1; end
    half=floor(win/2); nx=round(0.05/dx)+1; nz=round(0.05/dz)+1;
    xv=(0:nx-1)*dx; zv=(0:nz-1)*dz; [XX,ZZ]=meshgrid(xv,zv);
    if X.geometry_id(1)=="bilayer_2_3"
        hard=XX>=0.025;
        sd=1e3*(X.x_center_m-0.025);
        convention="x-25mm; negative soft, positive hard";
    else
        hard=hypot(XX-0.025,ZZ-0.025)<=0.010;
        sd=1e3*(0.010-hypot(X.x_center_m-0.025,X.z_center_m-0.025));
        convention="10mm-radius; negative soft background, positive hard inclusion";
    end
    hard_count=conv2(double(hard),ones(win),'same');
    lin=sub2ind(size(hard_count),X.cz,X.cx);
    hard_fraction=hard_count(lin)/(win^2); soft_fraction=1-hard_fraction;
    center_hard=X.cs_true>2.5;
    purity=soft_fraction; purity(center_hard)=hard_fraction(center_hard);
    radius_mm=1e3*half*max(dx,dz);
    side=repmat("soft",numel(idx),1); side(sd>0)="hard";
    side(abs(sd)<=radius_mm)="interface_band";
    B.signed_distance_mm(idx)=sd;
    B.abs_distance_to_interface_mm(idx)=abs(sd);
    B.patch_win_size_px(idx)=win;
    B.patch_radius_mm(idx)=radius_mm;
    B.patch_soft_fraction(idx)=soft_fraction;
    B.patch_hard_fraction(idx)=hard_fraction;
    B.patch_purity(idx)=purity;
    B.patch_is_mixed(idx)=purity<0.95;
    B.patch_is_strongly_mixed(idx)=purity<0.75;
    B.side(idx)=side;
    B.signed_distance_convention(idx)=convention;
end
end

function assert_no_inference_leakage(T)
assert(all(isfinite(T.sws_signed_error_pct)|isnan(T.sws_pred)), ...
    'Unexpected nonfinite diagnostic errors.');
assert(~any(contains(lower(string(T.Properties.VariableNames)),"trained")), ...
    'Unexpected training output found.');
fprintf(['Anti-leakage check: predictions were loaded before signed distance, ', ...
    'purity, and error diagnostics were computed.\n']);
end

function assert_side_coverage(T)
cases=unique(T.geometry_case,'stable');
for i=1:numel(cases)
    sides=unique(T.side(T.geometry_case==cases(i)));
    required=["soft","hard","interface_band"];
    missing=setdiff(required,sides);
    if isempty(missing)
        fprintf('%s contains soft, hard, and interface-band points.\n',cases(i));
    else
        warning('Test23:SideCoverage','%s missing sides: %s',cases(i),strjoin(missing,', '));
    end
end
end

function report_contamination_counts(T,primary)
P=T(T.model_name==primary,:);
fprintf('Primary-model diagnostic points: %d.\n',height(P));
fprintf('patch_purity <0.95: %d; <0.75: %d.\n', ...
    sum(P.patch_purity<0.95),sum(P.patch_purity<0.75));
fprintf('Mixed and high-error >20%%: %d; strongly mixed and high-error: %d.\n', ...
    sum(P.patch_purity<0.95&P.high_error_gt20), ...
    sum(P.patch_purity<0.75&P.high_error_gt20));
end

function S=summarize_errors(T,groups)
[G,S]=findgroups(T(:,cellstr(groups)));
S.N=splitapply(@numel,T.sws_signed_error_pct,G);
S.mean_signed_sws_error_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_signed_error_pct,G);
S.median_signed_sws_error_pct=splitapply(@(x)median(x,'omitnan'),T.sws_signed_error_pct,G);
S.mean_abs_sws_error_pct=splitapply(@(x)mean(x,'omitnan'),T.sws_abs_error_pct,G);
S.high_error_gt10_rate=splitapply(@(x)mean(x,'omitnan'),T.high_error_gt10,G);
S.high_error_gt20_rate=splitapply(@(x)mean(x,'omitnan'),T.high_error_gt20,G);
S.mean_q_signed_error=splitapply(@(x)mean(x,'omitnan'),T.q_signed_error,G);
S.median_q_signed_error=splitapply(@(x)median(x,'omitnan'),T.q_signed_error,G);
S.mean_q_abs_error=splitapply(@(x)mean(x,'omitnan'),T.q_abs_error,G);
S.mean_patch_purity=splitapply(@(x)mean(x,'omitnan'),T.patch_purity,G);
S.mean_signed_distance_mm=splitapply(@(x)mean(x,'omitnan'),T.signed_distance_mm,G);
end

function labels=purity_bin_labels(bin)
names=["[0,0.5)","[0.5,0.75)","[0.75,0.9)", ...
    "[0.9,0.95)","[0.95,1.0]"];
labels=strings(size(bin)); ok=isfinite(bin)&bin>=1&bin<=numel(names);
labels(ok)=names(bin(ok));
end

function labels=distance_bin_labels(bin,edges)
labels=strings(size(bin));
for i=1:numel(edges)-1
    if isinf(edges(i)), left="-Inf"; else, left=string(edges(i)); end
    if isinf(edges(i+1)), right="Inf"; else, right=string(edges(i+1)); end
    labels(bin==i)="("+left+","+right+"]";
end
end

function R=select_representative_patches(T,primary)
P=T(T.model_name==primary,:);
cases=unique(P.geometry_case,'stable'); R=table();
for ci=1:numel(cases)
    X=P(P.geometry_case==cases(ci),:);
    specs=[ ...
        struct('label',"hard_far_high_purity",'side',"hard",'near',false,'low',false)
        struct('label',"hard_near_low_purity",'side',"hard",'near',true,'low',true)
        struct('label',"soft_far_high_purity",'side',"soft",'near',false,'low',false)
        struct('label',"soft_near_low_purity",'side',"soft",'near',true,'low',true)
        struct('label',"interface_center_low_purity",'side',"interface_band",'near',true,'low',true)];
    for si=1:numel(specs)
        s=specs(si);
        if s.side=="hard"
            idx=X.signed_distance_mm>0;
        elseif s.side=="soft"
            idx=X.signed_distance_mm<0;
        else
            idx=X.side=="interface_band";
        end
        if s.near, idx=idx&X.abs_distance_to_interface_mm<=X.patch_radius_mm; end
        % Far means the complete rectangular REQ window stays on one side.
        if ~s.near, idx=idx&X.abs_distance_to_interface_mm>X.patch_radius_mm; end
        if s.low, idx=idx&X.patch_purity<0.95; else, idx=idx&X.patch_purity>=0.95; end
        C=X(idx,:); if isempty(C), warning('Test23:RepresentativeMissing','Missing %s / %s.',cases(ci),s.label); continue; end
        if s.low
            score=C.patch_purity-0.001*C.sws_abs_error_pct;
        else
            score=-C.patch_purity-0.001*C.abs_distance_to_interface_mm;
        end
        [~,j]=min(score); row=C(j,:);
        Ri=row(:,{'geometry_case','field_regime','frequency_hz','dx_m','dz_m', ...
            'REQ_M','condition_id','map_iz','map_ix','cx','cz','x_center_m', ...
            'z_center_m','side','signed_distance_mm','patch_purity', ...
            'patch_soft_fraction','patch_hard_fraction','q_true','q_pred', ...
            'sws_true','sws_pred','sws_signed_error_pct','sws_abs_error_pct'});
        Ri.representative_type=s.label;
        R=concat_tables(R,Ri);
    end
end
R=movevars(R,'representative_type','Before',1);
end

function plot_representative_spectra(R,T,condition_dir,OUT)
for i=1:height(R)
    file=condition_file(condition_dir,R.geometry_case(i),R.field_regime(i),R.REQ_M(i));
    if exist(file,'file')~=2
        warning('Test23:ConditionMissing','Cannot plot spectrum; missing %s.',file); continue;
    end
    S=load(file,'sim','cfg_sim','req_out');
    half=S.req_out.half_win;
    zi=(R.cz(i)-half):(R.cz(i)+half); xi=(R.cx(i)-half):(R.cx(i)+half);
    patch=S.sim.Uxz(zi,xi);
    [~,curve]=adaptive_req.quantile.compute_quantile_from_patch( ...
        patch,S.cfg_sim,S.req_out.req_cfg);
    mapping=adaptive_req.quantile.make_req_mapping(curve);
    q_model=R.q_pred(i); k_model=adaptive_req.quantile.quantile_to_k(mapping,q_model);
    rows=T(T.condition_id==R.condition_id(i)&T.map_iz==R.map_iz(i)& ...
        T.map_ix==R.map_ix(i),:);
    theory=rows(rows.model_name=="TheoryQDiscrete",:);
    q_theory=NaN; k_theory=NaN;
    if ~isempty(theory)
        q_theory=theory.q_pred(1);
        k_theory=adaptive_req.quantile.quantile_to_k(mapping,q_theory);
    end
    k_true=2*pi*R.frequency_hz(i)/R.sws_true(i);
    k_soft=2*pi*R.frequency_hz(i)/2; k_hard=2*pi*R.frequency_hz(i)/3;
    make_spectral_figure(patch,curve,R(i,:),k_true,k_soft,k_hard, ...
        k_model,k_theory,q_model,q_theory,OUT);
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

function make_spectral_figure(patch,C,R,ktrue,ksoft,khard,kmodel,ktheory,qmodel,qtheory,OUT)
fig=figure('Color','w','Units','centimeters','Position',[2 2 27 15]);
tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); imagesc(ax,real(patch)); axis(ax,'image'); colorbar(ax);
title(ax,'Local particle velocity','FontWeight','normal'); xlabel(ax,'x pixel'); ylabel(ax,'z pixel');
ax=nexttile(tl); imagesc(ax,C.kx,C.kz,log10(C.Ssm+eps)); axis(ax,'image'); colorbar(ax);
title(ax,'log_{10} 2D power spectrum','FontWeight','normal'); xlabel(ax,'k_x'); ylabel(ax,'k_z');
ax=nexttile(tl); plot(ax,C.k_cent,C.Srad_norm,'LineWidth',1.2); hold(ax,'on');
add_k_lines(ax,ktrue,ksoft,khard,kmodel,ktheory); title(ax,'Radial power spectrum','FontWeight','normal');
xlabel(ax,'k (rad/m)'); ylabel(ax,'Normalized power'); grid(ax,'on');
ax=nexttile(tl,[1 2]); plot(ax,C.k_cent,C.Ecum,'LineWidth',1.4); hold(ax,'on');
add_k_lines(ax,ktrue,ksoft,khard,kmodel,ktheory); yline(ax,qmodel,'m--','q model');
if isfinite(qtheory), yline(ax,qtheory,'c--','q theory'); end
yline(ax,R.q_true,'k:','q true'); xlabel(ax,'k (rad/m)'); ylabel(ax,'E_{cum}(k)');
title(ax,'Cumulative radial energy','FontWeight','normal'); grid(ax,'on'); ylim(ax,[0 1]);
ax=nexttile(tl); axis(ax,'off');
text(ax,0,1,sprintf(['%s\n%s / %s\nM=%g, dx=%.1f mm\nside=%s\n', ...
    'signed distance=%.2f mm\npurity=%.3f\nsoft=%.3f, hard=%.3f\n', ...
    'SWS error=%+.2f%%'],R.representative_type,R.geometry_case, ...
    R.field_regime,R.REQ_M,1e3*R.dx_m,R.side,R.signed_distance_mm, ...
    R.patch_purity,R.patch_soft_fraction,R.patch_hard_fraction, ...
    R.sws_signed_error_pct),'VerticalAlignment','top','Interpreter','none');
title(tl,'Representative interface-patch spectral diagnosis','FontWeight','normal');
file=sprintf('test23_spectrum_%s_%s.png',R.geometry_case,R.representative_type);
export_fig(fig,fullfile(OUT.spectral_dir,file));
end

function add_k_lines(ax,ktrue,ksoft,khard,kmodel,ktheory)
xline(ax,ksoft,'b:','k soft'); xline(ax,khard,'r:','k hard');
xline(ax,ktrue,'k-','k true'); xline(ax,kmodel,'m--','k model');
if isfinite(ktheory), xline(ax,ktheory,'c--','k theory'); end
end

function plot_distance_relationship(T,OUT,var,ylab,mode,prefix)
cases=unique(T.geometry_case,'stable'); regimes=regime_order(T.field_regime);
models=model_order(T.model_name);
for ci=1:numel(cases)
    fig=figure('Color','w','Units','centimeters','Position',[2 2 28 15]);
    tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for ri=1:numel(regimes)
        ax=nexttile(tl); hold(ax,'on'); X=T(T.geometry_case==cases(ci)&T.field_regime==regimes(ri),:);
        for mi=1:numel(models)
            Y=X(X.model_name==models(mi),:);
            if mode=="absolute", x=Y.abs_distance_to_interface_mm; else, x=Y.signed_distance_mm; end
            plot_binned_line(ax,x,Y.(var),models(mi));
        end
        if mode=="signed", xline(ax,0,'k-'); end
        title(ax,regimes(ri),'Interpreter','none','FontWeight','normal');
        xlabel(ax,replace(mode,["signed","absolute"],["Signed distance (mm)","Absolute distance (mm)"]));
        ylabel(ax,ylab); grid(ax,'on');
    end
    legend(nexttile(tl,4),'Location','eastoutside','Interpreter','none');
    title(tl,cases(ci),'Interpreter','none','FontWeight','normal');
    export_fig(fig,fullfile(OUT.figure_dir,prefix+"_"+cases(ci)+".png"));
end
end

function plot_binned_line(ax,x,y,name)
ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok); if isempty(x), return; end
edges=linspace(min(x),max(x),31); b=discretize(x,edges);
xc=nan(30,1); ym=nan(30,1);
for i=1:30, idx=b==i; xc(i)=mean(x(idx),'omitnan'); ym(i)=mean(y(idx),'omitnan'); end
plot(ax,xc,ym,'LineWidth',1.1,'DisplayName',short_model(name));
end

function plot_error_vs_purity(T,OUT)
cases=unique(T.geometry_case,'stable'); sides=["soft","hard"];
models=model_order(T.model_name);
for ci=1:numel(cases)
    fig=figure('Color','w','Units','centimeters','Position',[2 2 22 10]);
    tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    for si=1:2
        ax=nexttile(tl); hold(ax,'on'); X=T(T.geometry_case==cases(ci)&T.side==sides(si),:);
        for mi=1:numel(models)
            Y=X(X.model_name==models(mi),:);
            plot_binned_line(ax,Y.patch_purity,Y.sws_signed_error_pct,models(mi));
        end
        yline(ax,0,'k:'); xlabel(ax,'Patch purity'); ylabel(ax,'Signed SWS error (%)');
        title(ax,sides(si),'FontWeight','normal'); grid(ax,'on');
    end
    legend(nexttile(tl,2),'Location','eastoutside','Interpreter','none');
    title(tl,cases(ci),'Interpreter','none','FontWeight','normal');
    export_fig(fig,fullfile(OUT.figure_dir, ...
        "test23_signed_sws_error_vs_patch_purity_"+cases(ci)+".png"));
end
end

function plot_high_error_vs_purity(T,OUT)
cases=unique(T.geometry_case,'stable'); models=model_order(T.model_name);
for ci=1:numel(cases)
    fig=figure('Color','w','Units','centimeters','Position',[2 2 21 12]); ax=axes(fig); hold(ax,'on');
    X=T(T.geometry_case==cases(ci),:);
    for mi=1:numel(models)
        Y=X(X.model_name==models(mi),:); x=nan(5,1); y=nan(5,1);
        for bi=1:5
            idx=Y.patch_purity_bin==bi;
            x(bi)=mean(Y.patch_purity(idx),'omitnan');
            y(bi)=mean(Y.high_error_gt20(idx),'omitnan');
        end
        plot(ax,x,y,'-o', ...
            'MarkerSize',3,'DisplayName',short_model(models(mi)));
    end
    xlabel(ax,'Mean patch purity bin'); ylabel(ax,'High-error >20% rate');
    grid(ax,'on'); legend(ax,'Location','eastoutside','Interpreter','none');
    title(ax,cases(ci),'Interpreter','none','FontWeight','normal');
    export_fig(fig,fullfile(OUT.figure_dir, ...
        "test23_high_error_vs_patch_purity_"+cases(ci)+".png"));
end
end

function print_summary(T,primary)
P=T(T.model_name==primary,:); cases=unique(P.geometry_case,'stable');
fprintf('\n================ Test 23 summary ================\n');
for ci=1:numel(cases)
    X=P(P.geometry_case==cases(ci),:);
    near=X.abs_distance_to_interface_mm<=X.patch_radius_mm;
    far=X.abs_distance_to_interface_mm>2*X.patch_radius_mm;
    soft=X.side=="soft"; hard=X.side=="hard";
    rho=corr(X.patch_purity,X.sws_abs_error_pct,'Rows','complete','Type','Spearman');
    fprintf('\n%s\n',cases(ci));
    fprintf('High-error >20%%: far %.3f | near/interface %.3f.\n', ...
        mean(X.high_error_gt20(far),'omitnan'),mean(X.high_error_gt20(near),'omitnan'));
    fprintf('Mean signed SWS error: soft %+.3f%% | hard %+.3f%%.\n', ...
        mean(X.sws_signed_error_pct(soft),'omitnan'), ...
        mean(X.sws_signed_error_pct(hard),'omitnan'));
    fprintf('Spearman(patch purity, absolute error) = %.3f.\n',rho);
    hard_near=X.side=="interface_band"&X.signed_distance_mm>0;
    soft_near=X.side=="interface_band"&X.signed_distance_mm<0;
    hard_bias=mean(X.sws_signed_error_pct(hard_near),'omitnan');
    soft_bias=mean(X.sws_signed_error_pct(soft_near),'omitnan');
    if hard_bias<0
        fprintf(['Hard side shows SWS underestimation near mixed patches, ', ...
            'consistent with soft-region spectral contamination increasing ', ...
            'the effective selected k.\n']);
    else
        fprintf('Hard-side near-interface underestimation is not observed clearly.\n');
    end
    if soft_bias>0
        fprintf(['Soft side shows SWS overestimation near mixed patches, ', ...
            'consistent with hard-region spectral contamination decreasing ', ...
            'the effective selected k.\n']);
    else
        fprintf(['The soft-side response is asymmetric/no clear symmetric bias ', ...
            'and needs further inspection.\n']);
    end
end
fprintf('=================================================\n');
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
