%% analyze_test_27_roi_error_bars.m
% Test 27 ROI analysis: interior soft/hard regions for inclusion and bilayer.
% Error bars are mean +/- one standard deviation across independent
% frequency/regime/M conditions, after averaging pixels within each ROI.

clear; clc; close all;
format compact;

this_file=mfilename('fullpath');
root_dir=fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir=setup_adaptive_req();

SOURCE=fullfile(root_dir,'outputs','test_27_adaptive_window_edge_aware', ...
    'tables','test27_patch_level_results.csv');
OUT_ROOT=fullfile(root_dir,'outputs','test_27_adaptive_window_edge_aware');
TABLE_DIR=fullfile(OUT_ROOT,'tables');
FIGURE_DIR=fullfile(OUT_ROOT,'figures','roi_error_bars');
if exist(FIGURE_DIR,'dir')~=7, mkdir(FIGURE_DIR); end
assert(exist(SOURCE,'file')==2,'Test 27 full patch table missing: %s',SOURCE);

STRATEGIES=[ ...
    struct('name',"hybrid_baseline",'label',"Hybrid"); ...
    struct('name',"local_baseline",'label',"Local"); ...
    struct('name',"confidence_switch_local_hybrid_c080",'label',"Switch c0.8"); ...
    struct('name',"confidence_blend_local_hybrid_050_080",'label',"Blend"); ...
    struct('name',"edge_aware_interpolation_low_confidence_c080",'label',"Edge-aware"); ...
    struct('name',"small_window_then_edge_aware_c080",'label',"Small + edge")];

% All squares are specified by center and full width in physical millimeters.
ROIS=[ ...
    make_roi("bilayer_soft","Bilayer soft","bilayer_2_3",15,25,6,2); ...
    make_roi("bilayer_hard","Bilayer hard","bilayer_2_3",35,25,6,3); ...
    make_roi("inclusion_background","Inclusion background", ...
        "circular_inclusion_2_3",12,12,6,2); ...
    make_roi("inclusion_core","Inclusion hard core", ...
        "circular_inclusion_2_3",25,25,6,3)];

fprintf('\nTest 27 ROI error-bar analysis\n');
T=readtable(SOURCE,'TextType','string');
T.x_mm=1e3*T.x; T.z_mm=1e3*T.z;
condition_vars={'condition_key','field_regime','f0','M','dx'};
T_condition=table();

for ri=1:numel(ROIS)
    roi=ROIS(ri); half=roi.width_mm/2;
    in_roi=T.geometry==roi.geometry & ...
        abs(T.x_mm-roi.center_x_mm)<=half & ...
        abs(T.z_mm-roi.center_z_mm)<=half;
    X=T(in_roi,:);
    assert(~isempty(X),'ROI %s contains no patches.',roi.roi_id);
    keys=unique(X(:,condition_vars),'rows','stable');
    for ci=1:height(keys)
        idx=true(height(X),1);
        for v=string(condition_vars)
            idx=idx & string(X.(v))==string(keys.(v)(ci));
        end
        C=X(idx,:);
        for si=1:numel(STRATEGIES)
            variable=matlab.lang.makeValidName("sws_"+STRATEGIES(si).name);
            assert(ismember(variable,C.Properties.VariableNames), ...
                'Missing strategy column %s.',variable);
            pred=C.(variable);
            signed_error=100*(pred-C.true_SWS)./C.true_SWS;
            row=keys(ci,:);
            row.roi_id=roi.roi_id; row.roi_label=roi.label;
            row.geometry=roi.geometry; row.roi_true_sws=roi.true_sws;
            row.center_x_mm=roi.center_x_mm; row.center_z_mm=roi.center_z_mm;
            row.roi_width_mm=roi.width_mm;
            row.strategy_name=STRATEGIES(si).name;
            row.strategy_label=STRATEGIES(si).label;
            row.N_pixels=height(C);
            row.mean_sws_pred=mean(pred,'omitnan');
            row.pixel_sd_sws_pred=std(pred,'omitnan');
            row.mean_signed_error_pct=mean(signed_error,'omitnan');
            row.MAPE_pct=mean(abs(signed_error),'omitnan');
            row.high_error_10_pct=100*mean(abs(signed_error)>10,'omitnan');
            row.high_error_20_pct=100*mean(abs(signed_error)>20,'omitnan');
            row.mean_patch_purity=mean(C.patch_purity,'omitnan');
            row.mixed_patch_pct=100*mean(C.patch_purity<0.95,'omitnan');
            T_condition=concat_tables(T_condition,row);
        end
    end
end

T_summary=summarize_conditions(T_condition);
writetable(T_condition,fullfile(TABLE_DIR,'test27_roi_condition_metrics.csv'));
writetable(T_summary,fullfile(TABLE_DIR,'test27_roi_summary.csv'));

plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_sws_pred', ...
    'sd_sws_pred_across_conditions','Predicted SWS (m/s)', ...
    'test27_roi_predicted_sws_errorbars.png',FIGURE_DIR,true);
plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_MAPE_pct', ...
    'sd_MAPE_pct','MAPE (%)','test27_roi_mape_errorbars.png',FIGURE_DIR,false);
plot_roi_metric(T_summary,ROIS,STRATEGIES,'mean_signed_error_pct', ...
    'sd_signed_error_pct','Signed SWS error (%)', ...
    'test27_roi_signed_error_errorbars.png',FIGURE_DIR,false);

print_roi_summary(T_summary,ROIS,STRATEGIES);
fprintf('\nROI tables: %s\nROI figures: %s\n',TABLE_DIR,FIGURE_DIR);

function roi=make_roi(id,label,geometry,cx,cz,width,true_sws)
roi=struct('roi_id',string(id),'label',string(label),'geometry', ...
    string(geometry),'center_x_mm',cx,'center_z_mm',cz, ...
    'width_mm',width,'true_sws',true_sws);
end

function S=summarize_conditions(T)
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
S.mean_mixed_patch_pct=splitapply(@mean,T.mixed_patch_pct,G);
end

function plot_roi_metric(S,ROIS,STRATEGIES,value_var,error_var,ylab,file,folder,show_truth)
fig=figure('Color','w','Position',[80 80 1250 720]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
colors=lines(numel(STRATEGIES));
for ri=1:numel(ROIS)
    ax=nexttile(tl); hold(ax,'on'); roi=ROIS(ri);
    X=S(S.roi_id==roi.roi_id,:);
    values=nan(numel(STRATEGIES),1); errors=values;
    for si=1:numel(STRATEGIES)
        row=X(X.strategy_name==STRATEGIES(si).name,:);
        values(si)=row.(value_var); errors(si)=row.(error_var);
    end
    b=bar(ax,1:numel(STRATEGIES),values,.72,'FaceColor','flat');
    b.CData=colors; errorbar(ax,1:numel(STRATEGIES),values,errors, ...
        'k.','LineWidth',1.1,'CapSize',7);
    if show_truth, yline(ax,roi.true_sws,'k--','True SWS','LineWidth',1.2); end
    if string(value_var)=="mean_signed_error_pct", yline(ax,0,'k--'); end
    set(ax,'XTick',1:numel(STRATEGIES),'XTickLabel', ...
        string({STRATEGIES.label}),'XTickLabelRotation',18, ...
        'TickLabelInterpreter','none','FontSize',9);
    title(ax,roi.label,'Interpreter','none','FontSize',11);
    ylabel(ax,ylab); grid(ax,'on'); ax.GridAlpha=.15; box(ax,'off');
end
title(tl,"ROI means ± SD across frequency / regime / M conditions", ...
    'FontSize',13,'FontWeight','bold');
for ax=findall(fig,'Type','axes')'
    try
        ax.Toolbar.Visible='off';
    catch
        % Older MATLAB graphics objects may not expose a toolbar.
    end
end
exportgraphics(fig,fullfile(folder,file),'Resolution',220); close(fig);
end

function print_roi_summary(S,ROIS,STRATEGIES)
fprintf('\nROI summary: mean across conditions (MAPE / signed error)\n');
for ri=1:numel(ROIS)
    fprintf('  %s:\n',ROIS(ri).label);
    X=S(S.roi_id==ROIS(ri).roi_id,:);
    for si=1:numel(STRATEGIES)
        R=X(X.strategy_name==STRATEGIES(si).name,:);
        fprintf('    %-16s %6.2f%% / %+6.2f%% | purity %.3f | mixed %.1f%%\n', ...
            STRATEGIES(si).label,R.mean_MAPE_pct,R.mean_signed_error_pct, ...
            R.mean_patch_purity,R.mean_mixed_patch_pct);
    end
end
end

function T=concat_tables(A,B)
if isempty(A), T=B; return; end
B=B(:,A.Properties.VariableNames); T=[A;B];
end
