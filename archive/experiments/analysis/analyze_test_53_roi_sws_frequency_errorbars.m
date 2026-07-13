%% analyze_test_53_roi_sws_frequency_errorbars.m
% Test 53 post-hoc ROI frequency error bars.
%
% This script does not train and does not run REQ. It uses the saved Test 53
% held-out predictions and builds material-specific ROIs for M=2.
% ROIs use true material/distance only for evaluation and visualization.
%
% Output:
%   outputs/test_53_paper_final_clean_q_training/figures_no_theory/roi_frequency_m2/

clear; clc; close all;
format compact;

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();
set(groot,'defaultAxesFontSize',9,'defaultTextFontSize',9, ...
    'defaultLegendFontSize',8,'defaultAxesTitleFontSizeMultiplier',1.05);

CFG = struct();
CFG.SourceRoot = fullfile(root_dir, 'outputs', 'test_53_paper_final_clean_q_training');
CFG.ResultsFile = fullfile(CFG.SourceRoot, 'data', 'test38_results.mat');
CFG.ModelName = string(getenv_default('ADAPTIVE_REQ_TEST53_ROI_MODEL', ...
    'q_spectrum_plus_composition'));
CFG.M = env_number('ADAPTIVE_REQ_TEST53_ROI_M', 2);
CFG.RoiSideMm = env_number('ADAPTIVE_REQ_TEST53_ROI_SIDE_MM', 8);
CFG.ThinLayerHardWidthFraction = env_number('ADAPTIVE_REQ_TEST53_THIN_HARD_WIDTH_FRACTION', 0.80);
CFG.SoftCoreMinDistanceMm = env_number('ADAPTIVE_REQ_TEST53_ROI_SOFT_DIST_MM', 8);
CFG.HardCoreMinDistanceMm = env_number('ADAPTIVE_REQ_TEST53_ROI_HARD_DIST_MM', 4);
CFG.MinRoiPoints = env_number('ADAPTIVE_REQ_TEST53_ROI_MIN_POINTS', 3);

OUT = struct();
OUT.root_dir = fullfile(CFG.SourceRoot, 'figures_no_theory', 'roi_frequency_m2');
OUT.table_dir = fullfile(OUT.root_dir, 'tables');
OUT.figure_dir = fullfile(OUT.root_dir, 'figures');
OUT.overlay_dir = fullfile(OUT.root_dir, 'roi_overlays');
for d = string(struct2cell(OUT))'
    if exist(d,'dir') ~= 7, mkdir(d); end
end

fprintf('\nTest 53 ROI SWS-vs-frequency error bars\n');
fprintf('Source: %s\n', CFG.ResultsFile);
fprintf('Model: %s | M=%.1f | default ROI side %.1f mm\n', ...
    CFG.ModelName, CFG.M, CFG.RoiSideMm);
assert(exist(CFG.ResultsFile,'file') == 2, 'Missing Test53 results file: %s', CFG.ResultsFile);

S = load(CFG.ResultsFile, 'T_held');
T = S.T_held;
clear S;

T.model_name = string(T.model_name);
T.case_id = string(T.case_id);
T.case_family = string(T.case_family);
T.field_regime_ood = string(T.field_regime_ood);
T.condition_key = string(T.condition_key);
T = T(T.model_name == CFG.ModelName & abs(T.M - CFG.M) < 1e-9, :);
assert(~isempty(T), 'No rows found for model %s and M=%.1f.', CFG.ModelName, CFG.M);

fprintf('Loaded %d rows across %d conditions and %d geometries.\n', ...
    height(T), numel(unique(T.condition_key)), numel(unique(T.case_id)));

T_roi = build_condition_roi_table(T, CFG);
T_freq = summarize_frequency_roi(T_roi);
T_field = summarize_frequency_roi_by_field(T_roi);

writetable(T_roi, fullfile(OUT.table_dir, 'test53_roi_square_m2_condition_summary.csv'));
writetable(T_freq, fullfile(OUT.table_dir, 'test53_roi_square_m2_frequency_summary.csv'));
writetable(T_field, fullfile(OUT.table_dir, 'test53_roi_square_m2_frequency_by_field_type.csv'));

plot_all_cases(T_freq, T_field, OUT, CFG);
plot_roi_overlay_examples(T, T_roi, OUT, CFG);
audit_roi_purity(T_roi, OUT);

fprintf('\nSaved %d condition-ROI rows and %d frequency-ROI rows.\n', ...
    height(T_roi), height(T_freq));
fprintf('Tables: %s\nFigures: %s\nDone.\n', OUT.table_dir, OUT.root_dir);

%% Local functions

function R = build_condition_roi_table(T, CFG)
keys = unique(T.condition_key, 'stable');
rows = {};
for ki = 1:numel(keys)
    X = T(T.condition_key == keys(ki), :);
    specs = define_square_rois(X, CFG);
    for ri = 1:numel(specs)
        mask = specs(ri).mask;
        if nnz(mask) < CFG.MinRoiPoints
            continue;
        end
        Y = X(mask,:);
        rows{end+1,1} = table( ...
            X.condition_key(1), X.case_id(1), X.case_family(1), ...
            X.field_regime_ood(1), field_type_from_regime(X.field_regime_ood(1)), ...
            X.f0(1), X.M(1), ...
            string(specs(ri).roi_name), string(specs(ri).material_label), ...
            specs(ri).roi_center_x_mm, specs(ri).roi_center_z_mm, ...
            specs(ri).roi_width_mm, specs(ri).roi_height_mm, ...
            specs(ri).N_roi_rect_total, specs(ri).N_wrong_material, ...
            specs(ri).roi_material_fraction, height(Y), ...
            mean(Y.true_SWS,'omitnan'), std_omitnan(Y.true_SWS), ...
            mean(Y.sws_pred,'omitnan'), std_omitnan(Y.sws_pred), ...
            mean(Y.sws_abs_error_pct,'omitnan'), ...
            mean(Y.sws_signed_error_pct,'omitnan'), ...
            100*mean(Y.high_error20,'omitnan'), ...
            mean(Y.patch_purity,'omitnan'), ...
            mean(Y.distance_to_boundary_mm,'omitnan'), ...
            'VariableNames', {'condition_key','case_id','case_family', ...
            'field_regime_ood','field_type','f0','M','roi_name','material_label', ...
            'roi_center_x_mm','roi_center_z_mm','roi_width_mm','roi_height_mm', ...
            'N_roi_rect_total','N_wrong_material','roi_material_fraction','N', ...
            'true_sws_mean','true_sws_std','pred_sws_mean','pred_sws_std', ...
            'MAPE_pct','mean_signed_error_pct','high_error20_pct', ...
            'mean_patch_purity','mean_distance_to_boundary_mm'});
    end
end
if isempty(rows)
    R = table();
else
    R = vertcat(rows{:});
end
end

function specs = define_square_rois(X, CFG)
cs = X.true_SWS;
cmin = min(cs, [], 'omitnan');
cmax = max(cs, [], 'omitnan');
crange = cmax - cmin;
half = 0.5*CFG.RoiSideMm;
specs = struct('roi_name',{},'material_label',{},'roi_center_x_mm',{}, ...
    'roi_center_z_mm',{},'roi_width_mm',{},'roi_height_mm',{}, ...
    'mask',{},'square_mask',{},'N_roi_rect_total',{}, ...
    'N_wrong_material',{},'roi_material_fraction',{});

xmm = 1e3*X.x_center_m;
zmm = 1e3*X.z_center_m;
dist = X.distance_to_boundary_mm;
fam = string(X.case_family(1));

if ~isfinite(crange) || crange < 0.05 || fam == "homogeneous"
    mat = repmat("homogeneous", height(X), 1);
    desired = [median(xmm,'omitnan'), median(zmm,'omitnan')];
    [cx,cz,roi_mask,square_mask,frac,nwrong] = choose_rect_roi( ...
        xmm, zmm, mat, dist, "homogeneous", true(height(X),1), ...
        CFG.RoiSideMm, CFG.RoiSideMm, desired, CFG);
    specs(end+1) = make_spec("homogeneous center", "homogeneous", cx, cz, ...
        CFG.RoiSideMm, CFG.RoiSideMm, roi_mask, square_mask, frac, nwrong);
    return;
end

mat = classify_material(cs, fam);
roi_jobs = build_roi_jobs(mat, xmm, zmm, fam, CFG);

for li = 1:numel(roi_jobs)
    lab = roi_jobs(li).material_label;
    roi_name = roi_jobs(li).roi_name;
    base = mat == lab;
    if isfield(roi_jobs, 'component_mask') && ~isempty(roi_jobs(li).component_mask)
        base = base & roi_jobs(li).component_mask;
    end
    if ~any(base), continue; end
    if contains(lab, "soft")
        far = base & dist >= CFG.SoftCoreMinDistanceMm;
    elseif contains(lab, "hard")
        far = base & dist >= CFG.HardCoreMinDistanceMm;
    else
        far = base & dist >= CFG.HardCoreMinDistanceMm;
    end
    if nnz(far) < CFG.MinRoiPoints
        far = base & dist >= prctile(dist(base & isfinite(dist)), 70);
    end
    if nnz(far) < CFG.MinRoiPoints
        far = base;
    end

    desired = roi_jobs(li).desired_center_mm;
    width_mm = roi_jobs(li).width_mm;
    height_mm = roi_jobs(li).height_mm;
    [cx,cz,roi_mask,square_mask,frac,nwrong] = choose_rect_roi( ...
        xmm, zmm, mat, dist, lab, far, width_mm, height_mm, desired, CFG);
    if nnz(roi_mask) < CFG.MinRoiPoints
        [cx,cz,roi_mask,square_mask,frac,nwrong] = choose_rect_roi( ...
            xmm, zmm, mat, dist, lab, base, width_mm, height_mm, desired, CFG);
    end
    specs(end+1) = make_spec(roi_name, lab, cx, cz, width_mm, height_mm, ...
        roi_mask, square_mask, frac, nwrong);
end
end

function jobs = build_roi_jobs(mat, xmm, zmm, fam, CFG)
jobs = struct('roi_name',{},'material_label',{},'desired_center_mm',{}, ...
    'width_mm',{},'height_mm',{},'component_mask',{});
labels = unique(mat, 'stable');
order = ["soft background","soft","mid","hard","hard inclusion"];
labels = order(ismember(order, labels));
for li = 1:numel(labels)
    lab = labels(li);
    base = mat == lab;
    if ~any(base), continue; end
    width = CFG.RoiSideMm;
    height = CFG.RoiSideMm;
    if fam == "thin_layer" && lab == "hard"
        xr = max(xmm(base)) - min(xmm(base));
        width = max(min(CFG.RoiSideMm, CFG.ThinLayerHardWidthFraction*xr), min(xr, 1.0));
        height = CFG.RoiSideMm;
        jobs(end+1) = make_job("hard thin-layer rectangular ROI", lab, base, xmm, zmm, width, height); %#ok<AGROW>
        continue;
    end
    if fam == "two_inclusions" && lab == "hard inclusion"
        xmed = median(xmm(base),'omitnan');
        comp_masks = {base & xmm <= xmed, base & xmm > xmed};
        for ci = 1:2
            cbase = comp_masks{ci};
            if nnz(cbase) < 3, continue; end
            name = sprintf('hard/inclusion %d rectangular ROI', ci);
            cxr = max(xmm(cbase)) - min(xmm(cbase));
            czr = max(zmm(cbase)) - min(zmm(cbase));
            cwidth = max(min(CFG.RoiSideMm, 0.50*cxr), min(cxr, 2.0));
            cheight = max(min(CFG.RoiSideMm, 0.50*czr), min(czr, 2.0));
            jobs(end+1) = make_job(name, lab, cbase, xmm, zmm, cwidth, cheight); %#ok<AGROW>
        end
    else
        jobs(end+1) = make_job(pretty_roi_name(lab), lab, base, xmm, zmm, width, height); %#ok<AGROW>
    end
end
end

function job = make_job(name, lab, base, xmm, zmm, width, height)
job = struct();
job.roi_name = string(name);
job.material_label = string(lab);
job.desired_center_mm = [0.5*(min(xmm(base))+max(xmm(base))), ...
                         0.5*(min(zmm(base))+max(zmm(base)))];
job.width_mm = width;
job.height_mm = height;
job.component_mask = base;
end

function mat = classify_material(cs, fam)
cmin = min(cs, [], 'omitnan');
cmax = max(cs, [], 'omitnan');
crange = cmax - cmin;
mat = repmat("soft", numel(cs), 1);
if fam == "three_material"
    mat(cs >= cmin + crange/3 & cs < cmin + 2*crange/3) = "mid";
    mat(cs >= cmin + 2*crange/3) = "hard";
elseif contains(fam, ["inclusion","ellipse","two_inclusions"])
    mat(:) = "soft background";
    mat(cs >= cmin + 0.5*crange) = "hard inclusion";
else
    mat(cs >= cmin + 0.5*crange) = "hard";
end
end

function [cx,cz,roi_mask,rect_mask,best_frac,nwrong] = choose_rect_roi( ...
    xmm, zmm, mat, dist, target_label, candidate_mask, width_mm, height_mm, desired_center_mm, CFG)
idx = find(candidate_mask & isfinite(xmm) & isfinite(zmm));
if isempty(idx)
    cx = desired_center_mm(1);
    cz = desired_center_mm(2);
    rect_mask = abs(xmm - cx) <= 0.5*width_mm & abs(zmm - cz) <= 0.5*height_mm;
    roi_mask = rect_mask & mat == target_label;
    best_frac = nnz(roi_mask) / max(nnz(rect_mask),1);
    nwrong = nnz(rect_mask & mat ~= target_label);
    return;
end

best_score = -Inf;
cx = xmm(idx(1)); cz = zmm(idx(1));
rect_mask = false(size(xmm));
roi_mask = false(size(xmm));
best_frac = NaN;
nwrong = Inf;
for ii = 1:numel(idx)
    j = idx(ii);
    inside_domain = xmm(j) - 0.5*width_mm >= min(xmm) & ...
        xmm(j) + 0.5*width_mm <= max(xmm) & ...
        zmm(j) - 0.5*height_mm >= min(zmm) & ...
        zmm(j) + 0.5*height_mm <= max(zmm);
    rect = abs(xmm - xmm(j)) <= 0.5*width_mm & abs(zmm - zmm(j)) <= 0.5*height_mm;
    n = nnz(rect);
    if n < CFG.MinRoiPoints, continue; end
    same = rect & mat == target_label;
    wrong = rect & mat ~= target_label;
    frac = nnz(same) / max(n,1);
    md = mean(dist(same), 'omitnan');
    if ~isfinite(md), md = 0; end
    dc = hypot(xmm(j)-desired_center_mm(1), zmm(j)-desired_center_mm(2));
    % First maximize material purity, then center the ROI within the region,
    % then keep comparable support and distance from the interface.
    score = 1e9*frac - 1e5*dc + 1e2*min(n,200) + md;
    if ~inside_domain
        score = score - 5e8;
    end
    if score > best_score
        best_score = score;
        cx = xmm(j);
        cz = zmm(j);
        rect_mask = rect;
        roi_mask = same;
        best_frac = frac;
        nwrong = nnz(wrong);
    end
end
if ~isfinite(best_score)
    cx = desired_center_mm(1);
    cz = desired_center_mm(2);
    rect_mask = abs(xmm - cx) <= 0.5*width_mm & abs(zmm - cz) <= 0.5*height_mm;
    roi_mask = rect_mask & mat == target_label;
    best_frac = nnz(roi_mask) / max(nnz(rect_mask),1);
    nwrong = nnz(rect_mask & mat ~= target_label);
end
end

function s = make_spec(name, label, cx, cz, width_mm, height_mm, roi_mask, square_mask, frac, nwrong)
s = struct();
s.roi_name = name;
s.material_label = label;
s.roi_center_x_mm = cx;
s.roi_center_z_mm = cz;
s.roi_width_mm = width_mm;
s.roi_height_mm = height_mm;
s.mask = roi_mask;
s.square_mask = square_mask;
s.N_roi_rect_total = nnz(square_mask);
s.N_wrong_material = nwrong;
s.roi_material_fraction = frac;
end

function name = pretty_roi_name(label)
switch string(label)
    case "soft background"
        name = "soft/background ROI";
    case "hard inclusion"
        name = "hard/inclusion ROI";
    case "soft"
        name = "soft ROI";
    case "hard"
        name = "hard ROI";
    case "mid"
        name = "mid-material ROI";
    otherwise
        name = string(label) + " ROI";
end
end

function F = summarize_frequency_roi(R)
if isempty(R)
    F = table();
    return;
end
[G, groups] = findgroups(R(:, {'case_id','case_family','f0','M','roi_name','material_label'}));
F = groups;
F.N_conditions = splitapply(@numel, R.pred_sws_mean, G);
F.N_points_total = splitapply(@sum, R.N, G);
F.true_sws_mean = splitapply(@(x,w) weighted_mean(x,w), R.true_sws_mean, R.N, G);
F.pred_sws_mean = splitapply(@mean_omitnan, R.pred_sws_mean, G);
F.pred_sws_sd_across_conditions = splitapply(@std_omitnan, R.pred_sws_mean, G);
F.pred_sws_sem_across_conditions = F.pred_sws_sd_across_conditions ./ sqrt(max(F.N_conditions,1));
F.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), R.MAPE_pct, R.N, G);
F.mean_signed_error_pct = splitapply(@(x,w) weighted_mean(x,w), R.mean_signed_error_pct, R.N, G);
F.high_error20_pct = splitapply(@(x,w) weighted_mean(x,w), R.high_error20_pct, R.N, G);
F.mean_patch_purity = splitapply(@(x,w) weighted_mean(x,w), R.mean_patch_purity, R.N, G);
F.mean_roi_material_fraction = splitapply(@(x,w) weighted_mean(x,w), ...
    R.roi_material_fraction, R.N_roi_rect_total, G);
end

function F = summarize_frequency_roi_by_field(R)
if isempty(R)
    F = table();
    return;
end
[G, groups] = findgroups(R(:, {'case_id','case_family','field_type','f0','M','roi_name','material_label'}));
F = groups;
F.N_conditions = splitapply(@numel, R.pred_sws_mean, G);
F.N_points_total = splitapply(@sum, R.N, G);
F.true_sws_mean = splitapply(@(x,w) weighted_mean(x,w), R.true_sws_mean, R.N, G);
F.pred_sws_mean = splitapply(@mean_omitnan, R.pred_sws_mean, G);
F.pred_sws_sd_across_conditions = splitapply(@std_omitnan, R.pred_sws_mean, G);
F.pred_sws_sem_across_conditions = F.pred_sws_sd_across_conditions ./ sqrt(max(F.N_conditions,1));
F.MAPE_pct = splitapply(@(x,w) weighted_mean(x,w), R.MAPE_pct, R.N, G);
F.mean_signed_error_pct = splitapply(@(x,w) weighted_mean(x,w), R.mean_signed_error_pct, R.N, G);
F.high_error20_pct = splitapply(@(x,w) weighted_mean(x,w), R.high_error20_pct, R.N, G);
F.mean_patch_purity = splitapply(@(x,w) weighted_mean(x,w), R.mean_patch_purity, R.N, G);
F.mean_roi_material_fraction = splitapply(@(x,w) weighted_mean(x,w), ...
    R.roi_material_fraction, R.N_roi_rect_total, G);
end

function plot_all_cases(F, Ffield, OUT, CFG)
if isempty(F), return; end
cases = unique(F.case_id, 'stable');
for ci = 1:numel(cases)
    X = F(F.case_id == cases(ci), :);
    plot_case_frequency_errorbars(X, OUT, CFG);
    plot_case_frequency_by_field(Ffield(Ffield.case_id == cases(ci), :), OUT, CFG);
end
plot_family_overview(F, OUT, CFG);
end

function plot_case_frequency_errorbars(X, OUT, CFG)
roi_names = unique(X.roi_name, 'stable');
colors = roi_colors(roi_names);
fig = figure('Color','w','Units','centimeters','Position',[1 1 24 11]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile; hold(ax,'on');
for ri = 1:numel(roi_names)
    R = sortrows(X(X.roi_name == roi_names(ri),:), 'f0');
    errorbar(ax, R.f0, R.pred_sws_mean, R.pred_sws_sd_across_conditions, ...
        '-o', 'LineWidth',1.3, 'Color',colors(ri,:), ...
        'MarkerFaceColor', colors(ri,:), 'DisplayName', pretty_label(roi_names(ri)));
    plot(ax, R.f0, R.true_sws_mean, '--', 'LineWidth',1.0, ...
        'Color', colors(ri,:), 'HandleVisibility','off');
end
xlabel(ax,'Frequency (Hz)');
ylabel(ax,'SWS in ROI (m/s)');
title(ax, sprintf('%s: predicted SWS by ROI, M=%.0f', pretty_label(X.case_id(1)), CFG.M), ...
    'FontWeight','normal');
grid(ax,'on');
legend(ax,'Location','bestoutside','Interpreter','none');

ax = nexttile; hold(ax,'on');
for ri = 1:numel(roi_names)
    R = sortrows(X(X.roi_name == roi_names(ri),:), 'f0');
    plot(ax, R.f0, R.MAPE_pct, '-o', 'LineWidth',1.3, ...
        'Color',colors(ri,:), 'MarkerFaceColor', colors(ri,:), ...
        'DisplayName', pretty_label(roi_names(ri)));
end
xlabel(ax,'Frequency (Hz)');
ylabel('MAPE in ROI (%)');
title(ax,'ROI error versus frequency','FontWeight','normal');
grid(ax,'on');
legend(ax,'Location','bestoutside','Interpreter','none');

out_case = fullfile(OUT.figure_dir, sanitize(X.case_id(1)));
if exist(out_case,'dir') ~= 7, mkdir(out_case); end
exportgraphics(fig, fullfile(out_case, sprintf('test53_roi_sws_vs_frequency__%s__M%.0f.png', ...
    sanitize(X.case_id(1)), CFG.M)), 'Resolution', 220);
close(fig);
end

function plot_case_frequency_by_field(X, OUT, CFG)
if isempty(X), return; end
roi_names = unique(X.roi_name, 'stable');
field_types = ["directional 2D","diffuse 2D","partial 3D","diffuse 3D"];
field_types = field_types(ismember(field_types, unique(X.field_type)));
if isempty(field_types), return; end
colors = lines(numel(field_types));

fig = figure('Color','w','Units','centimeters', ...
    'Position',[1 1 max(22,10*numel(roi_names)) 14]);
tl = tiledlayout(fig, 2, numel(roi_names), 'TileSpacing','compact','Padding','compact');

for ri = 1:numel(roi_names)
    Rroi = X(X.roi_name == roi_names(ri),:);
    ax = nexttile(tl, ri); hold(ax,'on');
    true_by_f = groupsummary(Rroi, 'f0', 'mean', 'true_sws_mean');
    plot(ax, true_by_f.f0, true_by_f.mean_true_sws_mean, 'k--', ...
        'LineWidth',1.0, 'DisplayName','true SWS');
    for fi = 1:numel(field_types)
        R = sortrows(Rroi(Rroi.field_type == field_types(fi),:), 'f0');
        if isempty(R), continue; end
        errorbar(ax, R.f0, R.pred_sws_mean, R.pred_sws_sd_across_conditions, ...
            '-o', 'LineWidth',1.1, 'Color',colors(fi,:), ...
            'MarkerFaceColor', colors(fi,:), 'DisplayName', field_types(fi));
    end
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'SWS in ROI (m/s)');
    title(ax, pretty_label(roi_names(ri)), 'FontWeight','normal');
    grid(ax,'on');
    if ri == numel(roi_names)
        legend(ax,'Location','bestoutside','Interpreter','none');
    end

    ax = nexttile(tl, numel(roi_names)+ri); hold(ax,'on');
    for fi = 1:numel(field_types)
        R = sortrows(Rroi(Rroi.field_type == field_types(fi),:), 'f0');
        if isempty(R), continue; end
        plot(ax, R.f0, R.MAPE_pct, '-o', 'LineWidth',1.1, ...
            'Color',colors(fi,:), 'MarkerFaceColor', colors(fi,:), ...
            'DisplayName', field_types(fi));
    end
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'MAPE in ROI (%)');
    title(ax,'Error by field type','FontWeight','normal');
    grid(ax,'on');
end
title(tl, sprintf('%s: ROI curves separated by field type, M=%.0f', ...
    pretty_label(X.case_id(1)), CFG.M), 'FontWeight','normal');

out_case = fullfile(OUT.figure_dir, sanitize(X.case_id(1)));
if exist(out_case,'dir') ~= 7, mkdir(out_case); end
exportgraphics(fig, fullfile(out_case, sprintf('test53_roi_sws_vs_frequency_by_field__%s__M%.0f.png', ...
    sanitize(X.case_id(1)), CFG.M)), 'Resolution', 220);
close(fig);
end

function plot_family_overview(F, OUT, CFG)
families = unique(F.case_family, 'stable');
for fi = 1:numel(families)
    X = F(F.case_family == families(fi), :);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 24 13]);
    ax = axes(fig); hold(ax,'on');
    roi_names = unique(X.roi_name, 'stable');
    colors = roi_colors(roi_names);
    for ri = 1:numel(roi_names)
        R = X(X.roi_name == roi_names(ri),:);
        [G,g] = findgroups(R.f0);
        y = splitapply(@mean_omitnan, R.MAPE_pct, G);
        e = splitapply(@std_omitnan, R.MAPE_pct, G);
        errorbar(ax, g, y, e, '-o', 'LineWidth',1.3, ...
            'Color',colors(ri,:), 'MarkerFaceColor', colors(ri,:), ...
            'DisplayName', pretty_label(roi_names(ri)));
    end
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'MAPE across cases/regimes (%)');
    title(ax, sprintf('%s: ROI error summary, M=%.0f', pretty_label(families(fi)), CFG.M), ...
        'FontWeight','normal');
    grid(ax,'on');
    legend(ax,'Location','bestoutside','Interpreter','none');
    exportgraphics(fig, fullfile(OUT.figure_dir, sprintf('test53_roi_family_mape_vs_frequency__%s__M%.0f.png', ...
        sanitize(families(fi)), CFG.M)), 'Resolution', 220);
    close(fig);
end
end

function plot_roi_overlay_examples(T, R, OUT, CFG)
if isempty(R), return; end
cases = unique(R.case_id, 'stable');
for ci = 1:numel(cases)
    Rc = R(R.case_id == cases(ci),:);
    if isempty(Rc), continue; end
    key = Rc.condition_key(1);
    X = T(T.condition_key == key,:);
    plot_single_overlay(X, Rc(Rc.condition_key == key,:), OUT, CFG);
end
end

function plot_single_overlay(X, R, OUT, CFG)
[A,nz,nx,xu,zu] = rows_to_grid(X, X.true_SWS);
fig = figure('Color','w','Units','centimeters','Position',[1 1 14 12]);
ax = axes(fig);
imagesc(ax, xu, zu, A); axis(ax,'image'); set(ax,'YDir','normal');
c = colorbar(ax); ylabel(c,'True SWS (m/s)');
title(ax, sprintf('%s ROI placement, M=%.0f', pretty_label(X.case_id(1)), CFG.M), ...
    'FontWeight','normal');
xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
hold(ax,'on');
for i = 1:height(R)
    col = roi_colors(R.roi_name(i));
    x0 = R.roi_center_x_mm(i) - 0.5*R.roi_width_mm(i);
    z0 = R.roi_center_z_mm(i) - 0.5*R.roi_height_mm(i);
    rectangle(ax,'Position',[x0 z0 R.roi_width_mm(i) R.roi_height_mm(i)], ...
        'EdgeColor',col(1,:), 'LineWidth',1.8, 'LineStyle','--');
    text(ax, R.roi_center_x_mm(i), R.roi_center_z_mm(i), pretty_short_roi(R.roi_name(i)), ...
        'Color','w','FontWeight','bold','HorizontalAlignment','center', ...
        'BackgroundColor',[0 0 0 0.35], 'Interpreter','none');
end
out_case = fullfile(OUT.overlay_dir, sanitize(X.case_id(1)));
if exist(out_case,'dir') ~= 7, mkdir(out_case); end
exportgraphics(fig, fullfile(out_case, sprintf('test53_roi_overlay__%s.png', sanitize(X.condition_key(1)))), ...
    'Resolution', 220);
close(fig);
end

function audit_roi_purity(R, OUT)
if isempty(R), return; end
Audit = R(:, {'condition_key','case_id','case_family','field_regime_ood','field_type', ...
    'f0','M','roi_name','material_label','roi_width_mm','roi_height_mm','N_roi_rect_total', ...
    'N_wrong_material','roi_material_fraction','N','mean_patch_purity'});
Audit = sortrows(Audit, {'roi_material_fraction','N_roi_rect_total'}, {'ascend','descend'});
writetable(Audit, fullfile(OUT.table_dir, 'test53_roi_square_m2_purity_audit.csv'));

Bad = Audit(Audit.roi_material_fraction < 0.999, :);
fprintf('\nROI purity audit: %d/%d condition-ROIs have ROI material fraction < 0.999.\n', ...
    height(Bad), height(Audit));
if ~isempty(Bad)
    fprintf('  Worst ROI rectangle: %s | %s | f%.0f | %s | material fraction %.3f | wrong points %d/%d.\n', ...
        Bad.case_id(1), Bad.roi_name(1), Bad.f0(1), Bad.field_regime_ood(1), ...
        Bad.roi_material_fraction(1), Bad.N_wrong_material(1), Bad.N_roi_rect_total(1));
end
end

function [G,nz,nx,xu,zu] = rows_to_grid(T, values)
x = 1e3*double(T.x_center_m(:));
z = 1e3*double(T.z_center_m(:));
v = double(values(:));
xu = unique(x); zu = unique(z);
nx = numel(xu); nz = numel(zu);
G = nan(nz,nx);
[~,ix] = ismember(x, xu);
[~,iz] = ismember(z, zu);
ok = ix > 0 & iz > 0 & isfinite(v);
G = accumarray([iz(ok) ix(ok)], v(ok), [nz nx], @mean, NaN);
end

function t = field_type_from_regime(regime)
r = lower(string(regime));
if startsWith(r, "directional")
    t = "directional 2D";
elseif startsWith(r, "diffuse_2d")
    t = "diffuse 2D";
elseif startsWith(r, "partial_3d")
    t = "partial 3D";
elseif startsWith(r, "diffuse_3d")
    t = "diffuse 3D";
else
    t = "unknown";
end
end

function C = roi_colors(names)
names = string(names);
C = zeros(numel(names),3);
for i = 1:numel(names)
    n = lower(names(i));
    if contains(n, "soft") || contains(n, "background") || contains(n, "homogeneous")
        C(i,:) = [0.10 0.35 0.95];
    elseif contains(n, "mid")
        C(i,:) = [0.10 0.60 0.20];
    elseif contains(n, "hard") || contains(n, "inclusion")
        C(i,:) = [0.85 0.10 0.10];
    else
        C(i,:) = [0.20 0.20 0.20];
    end
end
end

function s = pretty_label(x)
s = strrep(string(x), "_", " ");
s = strrep(s, "p", ".");
end

function s = pretty_short_roi(x)
x = string(x);
if contains(x, "soft") || contains(x, "background")
    s = "soft ROI";
elseif contains(x, "hard") || contains(x, "inclusion")
    s = "hard ROI";
elseif contains(x, "mid")
    s = "mid ROI";
elseif contains(x, "homogeneous")
    s = "center ROI";
else
    s = x;
end
end

function y = mean_omitnan(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(x); end
end

function y = std_omitnan(x)
x = x(isfinite(x));
if numel(x) <= 1, y = NaN; else, y = std(x); end
end

function y = weighted_mean(x,w)
ok = isfinite(x) & isfinite(w) & w > 0;
if ~any(ok), y = NaN; else, y = sum(x(ok).*w(ok))/sum(w(ok)); end
end

function x = env_number(name, default)
v = strtrim(string(getenv(name)));
if v == "", x = default; else, x = str2double(v); end
if ~isfinite(x), x = default; end
end

function v = getenv_default(name, default)
v = string(getenv(name));
if strlength(strtrim(v)) == 0, v = string(default); end
end

function s = sanitize(s)
s = regexprep(char(string(s)), '[^A-Za-z0-9_-]', '_');
end
