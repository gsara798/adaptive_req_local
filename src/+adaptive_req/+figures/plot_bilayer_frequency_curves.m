function fig = plot_bilayer_frequency_curves(T_roi, model_name, roi_specs, varargin)
%PLOT_BILAYER_FREQUENCY_CURVES Plot predicted c_s vs frequency for one model.

p = inputParser;
p.FunctionName = 'adaptive_req.figures.plot_bilayer_frequency_curves';
addRequired(p, 'T_roi', @istable);
addRequired(p, 'model_name', @(x) ischar(x) || isstring(x));
addRequired(p, 'roi_specs');
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'FigureSize', [18 13], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FontSize', 13, @(x) isnumeric(x) && isscalar(x));
parse(p, T_roi, model_name, roi_specs, varargin{:});

model_name = string(model_name);
T = T_roi(string(T_roi.model_name) == model_name & ...
    T_roi.roi_name ~= "outside_roi" & T_roi.roi_name ~= "all", :);

fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');

colors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
markers = {'o', 's'};
handles = gobjects(0);
labels = strings(0);

for i = 1:numel(roi_specs)
    roi_name = string(roi_specs(i).name);
    Ti = T(string(T.roi_name) == roi_name, :);
    if isempty(Ti)
        continue;
    end
    Ti = sortrows(Ti, 'frequency_hz');
    c = colors(min(i, size(colors, 1)), :);
    h = errorbar(ax, Ti.frequency_hz, Ti.cs_pred_mean, Ti.cs_pred_std, ...
        [markers{min(i, numel(markers))} '-'], ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', c, ...
        'Color', c, ...
        'DisplayName', char(roi_specs(i).label));
    yline(ax, roi_specs(i).true_cs, '--', ...
        'Color', c, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    handles(end + 1) = h; %#ok<AGROW>
    labels(end + 1) = string(roi_specs(i).label); %#ok<AGROW>
end

grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'c_s (m/s)', 'Interpreter', 'tex');
if strlength(string(p.Results.Title)) > 0
    title(ax, string(p.Results.Title), 'Interpreter', 'none');
else
    title(ax, readable_model_name(model_name), 'Interpreter', 'none');
end
legend(ax, handles, cellstr(labels), 'Location', 'best');

if ismember('frequency_hz', string(T.Properties.VariableNames)) && ~isempty(T)
    freqs = unique(T.frequency_hz);
    xlim(ax, [min(freqs) - 40, max(freqs) + 40]);
    xticks(ax, freqs);
end

y_vals = [roi_specs.true_cs];
if ismember('cs_pred_mean', string(T.Properties.VariableNames)) && ~isempty(T)
    y_vals = [y_vals, T.cs_pred_mean.'];
end
ylim(ax, [min(y_vals) - 0.35, max(y_vals) + 0.35]);

adaptive_req.templates.apply_paper_style(fig, ax, ...
    'Times New Roman', p.Results.FontSize);
set(fig, 'Units', 'centimeters');
fig.Position(3:4) = p.Results.FigureSize;

end

function out = readable_model_name(model_name)

out = string(model_name);
out = replace(out, "TheoryDiffuse3D", "Theory diffuse 3D");
out = replace(out, "GlobalQSingleModel", "Global q model");
out = replace(out, "LocalOnly", "Local only");
out = replace(out, "HybridLocalGlobal", "Hybrid");

end
