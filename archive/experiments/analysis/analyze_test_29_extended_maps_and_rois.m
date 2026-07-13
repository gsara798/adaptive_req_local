%% analyze_test_29_extended_maps_and_rois.m
% Extended Test 29 maps and the same four 6x6 mm ROIs used after Test 27.
% ROI error bars are mean +/- one standard deviation across independent
% frequency/regime/M conditions after averaging pixels within each ROI.

clear; clc; close all;
format compact;

this_file=mfilename('fullpath'); root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir); root_dir=setup_adaptive_req(); adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.15, ...
    'defaultAxesLabelFontSizeMultiplier',1.05);

SOURCE=fullfile(root_dir,'outputs','test_29_oracle_halfwindow_graph_tv', ...
    'tables','test29_patch_level_results.csv');
OUT_ROOT=fullfile(root_dir,'outputs','test_29_oracle_halfwindow_graph_tv');
TABLE_DIR=fullfile(OUT_ROOT,'tables');
FIGURE_DIR=fullfile(OUT_ROOT,'figures');
MAP_DIR=fullfile(FIGURE_DIR,'maps_heterogeneous_all_with_theory');
REP_DIR=fullfile(FIGURE_DIR,'maps_representative_extended_with_theory');
ROI_DIR=fullfile(FIGURE_DIR,'roi_error_bars');
for folder=string({MAP_DIR,REP_DIR,ROI_DIR})
    if exist(folder,'dir')~=7, mkdir(folder); end
end
assert(exist(SOURCE,'file')==2,'Full Test 29 patch table missing: %s',SOURCE);

STRATEGIES=[ ...
    make_strategy("hybrid_baseline","Hybrid"); ...
    make_strategy("local_baseline","Local"); ...
    make_strategy("theory_discrete","Theory discrete"); ...
    make_strategy("confidence_switch_c080","Switch c0.8"); ...
    make_strategy("halfwindow_switch_c080","Half-window Switch"); ...
    make_strategy("graph_halfwindow_switch_l075","Graph + half-window")];

% Same physical ROIs as analyze_test_27_roi_error_bars.m.
ROIS=[ ...
    make_roi("bilayer_soft","Bilayer soft","bilayer_2_3",15,25,6,2); ...
    make_roi("bilayer_hard","Bilayer hard","bilayer_2_3",35,25,6,3); ...
    make_roi("inclusion_background","Inclusion background", ...
        "circular_inclusion_2_3",12,12,6,2); ...
    make_roi("inclusion_core","Inclusion hard core", ...
        "circular_inclusion_2_3",25,25,6,3)];

fprintf('\nTest 29 extended maps and ROI analysis\n');
T=readtable(SOURCE,'TextType','string'); T.x_mm=1e3*T.x; T.z_mm=1e3*T.z;
assert(all(ismember("sws_"+string({STRATEGIES.name}),string(T.Properties.VariableNames))), ...
    'Test 29 patch table lacks one or more requested strategies.');

%% Export every heterogeneous condition map and a 24-map representative set

hetero_keys=unique(T.condition_key(T.geometry_type~="homogeneous"),'stable');
representative_keys=select_representative_keys(T);
fprintf('Saving %d heterogeneous condition maps; %d are copied to the representative set.\n', ...
    numel(hetero_keys),numel(representative_keys));
for i=1:numel(hetero_keys)
    C=T(T.condition_key==hetero_keys(i),:);
    folder=fullfile(MAP_DIR,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
    if exist(folder,'dir')~=7, mkdir(folder); end
    file=fullfile(folder,"test29_extended__"+sanitize(hetero_keys(i))+".png");
    if exist(file,'file')~=2, plot_condition_map(C,file); end
    if ismember(hetero_keys(i),representative_keys)
        rep_folder=fullfile(REP_DIR,sanitize(C.geometry(1)),sanitize(C.field_regime(1)));
        if exist(rep_folder,'dir')~=7, mkdir(rep_folder); end
        rep_file=fullfile(rep_folder,"test29_representative__"+sanitize(hetero_keys(i))+".png");
        if exist(rep_file,'file')~=2, copyfile(file,rep_file); end
    end
    if mod(i,12)==0||i==numel(hetero_keys), fprintf('  Maps: %d/%d.\n',i,numel(hetero_keys)); end
end

%% ROI condition means and summaries

condition_vars={'condition_key','field_regime','f0','M','dx'};
T_condition=table();
for ri=1:numel(ROIS)
    roi=ROIS(ri); half=roi.width_mm/2;
    in_roi=T.geometry==roi.geometry&abs(T.x_mm-roi.center_x_mm)<=half& ...
        abs(T.z_mm-roi.center_z_mm)<=half;
    X=T(in_roi,:); assert(~isempty(X),'ROI %s contains no patches.',roi.roi_id);
    keys=unique(X(:,condition_vars),'rows','stable');
    for ci=1:height(keys)
        idx=true(height(X),1);
        for v=string(condition_vars), idx=idx&string(X.(v))==string(keys.(v)(ci)); end
        C=X(idx,:);
        for si=1:numel(STRATEGIES)
            pred=C.("sws_"+STRATEGIES(si).name);
            error=100*(pred-C.true_SWS)./C.true_SWS;
            row=keys(ci,:); row.roi_id=roi.roi_id; row.roi_label=roi.label;
            row.geometry=roi.geometry; row.roi_true_sws=roi.true_sws;
            row.center_x_mm=roi.center_x_mm; row.center_z_mm=roi.center_z_mm;
            row.roi_width_mm=roi.width_mm; row.strategy_name=STRATEGIES(si).name;
            row.strategy_label=STRATEGIES(si).label; row.N_pixels=height(C);
            row.mean_sws_pred=mean(pred,'omitnan');
            row.mean_signed_error_pct=mean(error,'omitnan');
            row.MAPE_pct=mean(abs(error),'omitnan');
            row.high_error_10_pct=100*mean(abs(error)>10,'omitnan');
            row.high_error_20_pct=100*mean(abs(error)>20,'omitnan');
            row.mean_patch_purity=mean(C.patch_purity,'omitnan');
            row.low_confidence_pct=100*mean(C.confidence<.8,'omitnan');
            row.halfwindow_accepted_pct=100*mean(C.halfwindow_accepted,'omitnan');
            T_condition=concat_tables(T_condition,row);
        end
    end
end
T_summary=summarize_roi_conditions(T_condition);
writetable(T_condition,fullfile(TABLE_DIR,'test29_roi_condition_metrics.csv'));
writetable(T_summary,fullfile(TABLE_DIR,'test29_roi_summary.csv'));

plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_MAPE_pct','sd_MAPE_pct', ...
    'MAPE (%)','test29_roi_mape_errorbars.png',ROI_DIR,false);
plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_sws_pred', ...
    'sd_sws_pred_across_conditions','Predicted SWS (m/s)', ...
    'test29_roi_predicted_sws_errorbars.png',ROI_DIR,true);
plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_signed_error_pct', ...
    'sd_signed_error_pct','Signed SWS error (%)', ...
    'test29_roi_signed_error_errorbars.png',ROI_DIR,false);
plot_roi_mape_heatmap(T_summary,ROIS,STRATEGIES,ROI_DIR);
print_roi_summary(T_summary,ROIS,STRATEGIES);
fprintf('\nROI tables: %s\nAll maps: %s\nRepresentative maps: %s\n', ...
    TABLE_DIR,MAP_DIR,REP_DIR);

%% Local functions

function s=make_strategy(name,label)
s=struct('name',string(name),'label',string(label));
end

function roi=make_roi(id,label,geometry,cx,cz,width,true_sws)
roi=struct('roi_id',string(id),'label',string(label),'geometry',string(geometry), ...
    'center_x_mm',cx,'center_z_mm',cz,'width_mm',width,'true_sws',true_sws);
end

function keys=select_representative_keys(T)
keys=strings(0,1); cases=["bilayer_2_3","circular_inclusion_2_3"];
regimes=["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
design=[300 2;500 3;600 4];
for geometry=cases
    for regime=regimes
        for di=1:size(design,1)
            idx=T.geometry==geometry&T.field_regime==regime& ...
                T.f0==design(di,1)&T.M==design(di,2);
            candidate=unique(T.condition_key(idx),'stable');
            if ~isempty(candidate), keys(end+1)=candidate(1); end %#ok<AGROW>
        end
    end
end
keys=unique(keys,'stable');
end

function plot_condition_map(C,file)
C=sortrows(C,{'map_iz','map_ix'});
S=C.sws_confidence_switch_c080; H=C.sws_halfwindow_switch_c080;
G=C.sws_graph_halfwindow_switch_l075; D=C.sws_theory_discrete;
eS=abs(100*(S-C.true_SWS)./C.true_SWS); eH=abs(100*(H-C.true_SWS)./C.true_SWS);
eG=abs(100*(G-C.true_SWS)./C.true_SWS); eD=abs(100*(D-C.true_SWS)./C.true_SWS);
maps={map_from_rows(C,C.true_SWS),map_from_rows(C,S),map_from_rows(C,H), ...
    map_from_rows(C,G),map_from_rows(C,D),map_from_rows(C,eS),map_from_rows(C,eH), ...
    map_from_rows(C,eG),map_from_rows(C,eD), ...
    map_from_rows(C,C.confidence),map_from_rows(C,C.halfwindow_accepted), ...
    map_from_rows(C,C.halfwindow_angle_rad),map_from_rows(C,C.halfwindow_width_ratio), ...
    map_from_rows(C,G-S)};
titles=["True SWS","Switch","Half-window Switch","Graph + half-window","Theory discrete", ...
    "Switch error (%)","Half-window error (%)","Graph error (%)","Theory error (%)","Confidence", ...
    "Half-window accepted","Selected angle (rad)","Spectral width ratio","Graph - Switch"];
fig=figure('Color','w','Visible','off','Position',[20 30 1500 850]);
tl=tiledlayout(fig,3,5,'TileSpacing','compact','Padding','compact');
for i=1:numel(maps)
    ax=nexttile(tl); imagesc(ax,maps{i}); axis(ax,'image');
    set(ax,'XTick',[],'YTick',[]); title(ax,titles(i),'FontSize',9); colorbar(ax);
    if i<=5, clim(ax,[1.5 3.5]); elseif i>=6&&i<=9, clim(ax,[0 50]);
    elseif i==10||i==11, clim(ax,[0 1]); elseif i==14, clim(ax,[-.5 .5]); end
end
title(tl,C.condition_key(1),'Interpreter','none','FontSize',11);
exportgraphics(fig,file,'Resolution',180); close(fig);
end

function S=summarize_roi_conditions(T)
groups={'roi_id','roi_label','geometry','roi_true_sws','center_x_mm', ...
    'center_z_mm','roi_width_mm','strategy_name','strategy_label'};
[G,S]=findgroups(T(:,groups));
S.N_conditions=splitapply(@numel,T.MAPE_pct,G);
S.total_pixels=splitapply(@sum,T.N_pixels,G);
S.mean_sws_pred=splitapply(@mean,T.mean_sws_pred,G);
S.sd_sws_pred_across_conditions=splitapply(@std,T.mean_sws_pred,G);
S.mean_MAPE_pct=splitapply(@mean,T.MAPE_pct,G);
S.sd_MAPE_pct=splitapply(@std,T.MAPE_pct,G);
S.mean_signed_error_pct=splitapply(@mean,T.mean_signed_error_pct,G);
S.sd_signed_error_pct=splitapply(@std,T.mean_signed_error_pct,G);
S.mean_high_error_10_pct=splitapply(@mean,T.high_error_10_pct,G);
S.mean_high_error_20_pct=splitapply(@mean,T.high_error_20_pct,G);
S.mean_patch_purity=splitapply(@mean,T.mean_patch_purity,G);
S.mean_low_confidence_pct=splitapply(@mean,T.low_confidence_pct,G);
S.mean_halfwindow_accepted_pct=splitapply(@mean,T.halfwindow_accepted_pct,G);
end

function plot_roi_metric(S,ROIS,STRATEGIES,value_var,error_var,ylab,file,folder,truth)
fig=figure('Color','w','Position',[70 60 1300 720]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact'); colors=lines(numel(STRATEGIES));
for ri=1:numel(ROIS)
    ax=nexttile(tl); hold(ax,'on'); X=S(S.roi_id==ROIS(ri).roi_id,:);
    values=nan(numel(STRATEGIES),1); errors=values;
    for si=1:numel(STRATEGIES)
        row=X(X.strategy_name==STRATEGIES(si).name,:); values(si)=row.(value_var); errors(si)=row.(error_var);
    end
    b=bar(ax,1:numel(STRATEGIES),values,.72,'FaceColor','flat'); b.CData=colors;
    errorbar(ax,1:numel(STRATEGIES),values,errors,'k.','LineWidth',1.1,'CapSize',7);
    if truth, yline(ax,ROIS(ri).true_sws,'k--','True SWS','LineWidth',1.1); end
    if string(value_var)=="mean_signed_error_pct", yline(ax,0,'k--'); end
    set(ax,'XTick',1:numel(STRATEGIES),'XTickLabel',string({STRATEGIES.label}), ...
        'XTickLabelRotation',18,'TickLabelInterpreter','none','FontSize',8);
    title(ax,ROIS(ri).label,'Interpreter','none'); ylabel(ax,ylab); grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off');
end
title(tl,'ROI means ± SD across frequency / regime / M conditions','FontWeight','bold');
exportgraphics(fig,fullfile(folder,file),'Resolution',220); close(fig);
end

function plot_roi_mape_heatmap(S,ROIS,STRATEGIES,folder)
M=nan(numel(ROIS),numel(STRATEGIES));
for ri=1:numel(ROIS)
    for si=1:numel(STRATEGIES)
        row=S(S.roi_id==ROIS(ri).roi_id&S.strategy_name==STRATEGIES(si).name,:);
        M(ri,si)=row.mean_MAPE_pct;
    end
end
fig=figure('Color','w','Position',[100 100 950 420]); ax=axes(fig); imagesc(ax,M); colorbar(ax);
set(ax,'YTick',1:numel(ROIS),'YTickLabel',string({ROIS.label}), ...
    'XTick',1:numel(STRATEGIES),'XTickLabel',string({STRATEGIES.label}), ...
    'XTickLabelRotation',18,'TickLabelInterpreter','none'); title(ax,'ROI MAPE (%)');
for r=1:size(M,1)
    for c=1:size(M,2)
        text(ax,c,r,sprintf('%.2f',M(r,c)),'HorizontalAlignment','center', ...
            'Color',contrast_color(M(r,c),M));
    end
end
exportgraphics(fig,fullfile(folder,'test29_roi_mape_heatmap.png'),'Resolution',220); close(fig);
end

function color=contrast_color(value,M)
if value>mean(M(:),'omitnan'), color='w'; else, color='k'; end
end

function print_roi_summary(S,ROIS,STRATEGIES)
fprintf('\nROI summary: MAPE mean +/- SD across conditions\n');
for ri=1:numel(ROIS)
    fprintf('  %s:\n',ROIS(ri).label); X=S(S.roi_id==ROIS(ri).roi_id,:);
    for si=1:numel(STRATEGIES)
        R=X(X.strategy_name==STRATEGIES(si).name,:);
        fprintf('    %-21s %6.2f +/- %5.2f%% | signed %+6.2f%% | Ncond %d\n', ...
            STRATEGIES(si).label,R.mean_MAPE_pct,R.sd_MAPE_pct, ...
            R.mean_signed_error_pct,R.N_conditions);
    end
end
end

function [A,iz,ix]=map_from_rows(T,values)
uz=unique(T.map_iz); ux=unique(T.map_ix); [~,iz]=ismember(T.map_iz,uz); [~,ix]=ismember(T.map_ix,ux);
A=nan(numel(uz),numel(ux)); A(sub2ind(size(A),iz,ix))=values;
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end; B=B(:,A.Properties.VariableNames); T=[A;B];
end

function value=sanitize(value), value=regexprep(char(string(value)),'[^A-Za-z0-9_-]','_'); end
