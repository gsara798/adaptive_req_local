%% run_test_18_clean_field_regime_training_dataset.m
% Generate the clean Test 18 field-regime training dataset.
%
% No model is trained here. The field-regime labels and wave counts are
% retained for validation only and are explicitly excluded from the future
% primary operational predictor sets.

clear; clc; close all;
format compact;

%% Project setup

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(this_file));
addpath(root_dir);
root_dir = setup_adaptive_req();
adaptive_req.templates.setup_style();

CFG = adaptive_req.config.load_profile_config( ...
    'test_18_clean_field_regime_training_dataset', ...
    'RootDir', root_dir);

% Optional short validation:
% setenv('ADAPTIVE_REQ_TEST18_MODE','smoke') before running this script.
SMOKE_TEST = strcmpi(getenv('ADAPTIVE_REQ_TEST18_MODE'), 'smoke');
if SMOKE_TEST
    CFG.FIELD_REGIMES = CFG.FIELD_REGIMES([1 4 7 10]);
    CFG.GRID.f0 = 500;
    CFG.GRID.cs_bg = [2 3];
    CFG.GRID.dx_dz = 0.5e-3;
    CFG.GRID.REQ_M = [2 3];
    CFG.EXP.num_realizations = 1;
    CFG.EXP.num_patches = 3;
end

%% Stable, resumable output folder

SAVE.root_dir = fullfile(root_dir, 'outputs', CFG.EXP.name);
if SMOKE_TEST
    SAVE.run_name = 'smoke_test';
else
    SAVE.run_name = 'dataset';
end
SAVE.output_dir = fullfile(SAVE.root_dir, SAVE.run_name);
SAVE.condition_dir = fullfile(SAVE.output_dir, 'conditions');
SAVE.data_dir = fullfile(SAVE.output_dir, 'data');
SAVE.table_dir = fullfile(SAVE.output_dir, 'tables');
SAVE.figure_dir = fullfile(SAVE.output_dir, 'figures');
make_dirs(SAVE);

%% Build linked-regime condition table

T_conditions = build_condition_table(CFG);
writetable(T_conditions, fullfile(SAVE.table_dir, ...
    'test18_condition_design.csv'));

fprintf('\nTest 18 clean field-regime dataset.\n');
fprintf('Physical conditions : %d\n', height(T_conditions));
fprintf('Realizations        : %d\n', CFG.EXP.num_realizations);
fprintf('Patches per M       : %d\n', CFG.EXP.num_patches);
fprintf('REQ M values        : %s\n', mat2str(CFG.GRID.REQ_M));
fprintf('Expected rows       : %d\n', height(T_conditions) * ...
    CFG.EXP.num_realizations * CFG.EXP.num_patches * numel(CFG.GRID.REQ_M));
fprintf('Output folder:\n%s\n', SAVE.output_dir);

%% Generate or reuse each physical condition

T_status = table();
for ci = 1:height(T_conditions)
    cond = T_conditions(ci, :);
    condition_file = fullfile(SAVE.condition_dir, sprintf( ...
        'condition_%04d_%s.mat', cond.condition_id, ...
        sanitize_filename(cond.field_regime_variant)));

    fprintf('\n=== Condition %d / %d ===\n', ci, height(T_conditions));
    fprintf('%s | N=%d | f=%g Hz | cs=%g m/s | dx=dz=%.3f mm\n', ...
        cond.field_regime_variant, cond.SIM_Nwaves, cond.SIM_f0, ...
        cond.SIM_cs_bg, cond.SIM_dx*1e3);

    if exist(condition_file, 'file') == 2
        fprintf('Reusing completed condition file.\n');
        S = load(condition_file, 'META');
        meta = S.META;
    else
        regime = CFG.FIELD_REGIMES(cond.regime_idx);
        [T_condition, REP, META] = ...
            adaptive_req.studies.run_clean_field_regime_condition( ...
            CFG, regime, cond.SIM_f0, cond.SIM_cs_bg, cond.SIM_dx, ...
            cond.condition_id, 'Verbose', true);
        save(condition_file, 'T_condition', 'REP', 'META', '-v7.3');
        meta = META;
    end

    status_i = struct2table(meta, 'AsArray', true);
    T_status = concat_tables(T_status, status_i);
    writetable(T_status, fullfile(SAVE.table_dir, ...
        'test18_generation_status.csv'));
end

%% Aggregate condition files

fprintf('\nAggregating condition files...\n');
T_dataset = table();
representative_cells = {};
for ci = 1:height(T_conditions)
    cond = T_conditions(ci, :);
    condition_file = fullfile(SAVE.condition_dir, sprintf( ...
        'condition_%04d_%s.mat', cond.condition_id, ...
        sanitize_filename(cond.field_regime_variant)));
    S = load(condition_file, 'T_condition', 'REP');
    T_dataset = concat_tables(T_dataset, S.T_condition);
    if ~isempty(S.REP)
        representative_cells{end+1,1} = S.REP; %#ok<SAGROW>
    end
end
if isempty(representative_cells)
    REPRESENTATIVES = struct([]);
else
    REPRESENTATIVES = vertcat(representative_cells{:});
end

%% Validation and derived metadata

required = [ ...
    "field_regime_label", "field_regime_variant", "SIM_Nwaves", ...
    "SIM_Is2D", "SIM_ForceInPlaneWave", "SIM_f0", "SIM_cs_bg", ...
    "SIM_dx", "SIM_dz", "REQ_M", "q_theory", "req_mapping", ...
    "global_req_mapping", "global_radial_entropy", "global_ang_entropy"];
missing = setdiff(required, string(T_dataset.Properties.VariableNames));
assert(isempty(missing), 'Missing Test 18 variables: %s', strjoin(missing, ', '));
assert(~ismember('req_curve', T_dataset.Properties.VariableNames), ...
    'Test 18 unexpectedly contains the heavy req_curve.');
assert(all(T_dataset.SIM_dx == T_dataset.SIM_dz), ...
    'Test 18 currently requires dx=dz.');
assert(all(ismember(unique(T_dataset.field_regime_label), ...
    ["directional_2D","diffuse_2D","partial_3D","diffuse_3D"])), ...
    'Unexpected field regime label.');

expected_rows = height(T_conditions) * CFG.EXP.num_realizations * ...
    CFG.EXP.num_patches * numel(CFG.GRID.REQ_M);
assert(height(T_dataset) == expected_rows, ...
    'Generated rows (%d) do not match expected rows (%d).', ...
    height(T_dataset), expected_rows);

PREDICTOR_POLICY = table( ...
    CFG.PREDICTOR_POLICY.diagnostic_only(:), ...
    repmat("diagnostic_only_do_not_use_in_primary_operational_model", ...
    numel(CFG.PREDICTOR_POLICY.diagnostic_only), 1), ...
    'VariableNames', {'variable_name','policy'});

%% Save final dataset

MC = struct();
MC.CFG = CFG;
MC.SAVE = SAVE;
MC.T_conditions = T_conditions;
MC.T_status = T_status;
MC.expected_rows = expected_rows;
MC.actual_rows = height(T_dataset);
MC.predictor_policy = PREDICTOR_POLICY;

save(fullfile(SAVE.data_dir, 'test18_clean_field_regime_dataset.mat'), ...
    'T_dataset', 'MC', 'CFG', 'REPRESENTATIVES', 'PREDICTOR_POLICY', '-v7.3');
writetable(remove_heavy_columns(T_dataset), fullfile(SAVE.table_dir, ...
    'test18_clean_field_regime_dataset.csv'));
writetable(PREDICTOR_POLICY, fullfile(SAVE.table_dir, ...
    'test18_predictor_policy.csv'));

%% Sanity-check figures

plot_representative_velocity(REPRESENTATIVES, SAVE.figure_dir, CFG);
plot_representative_spectra(REPRESENTATIVES, SAVE.figure_dir, CFG);
plot_q_distribution(T_dataset, SAVE.figure_dir, CFG);
plot_dataset_occupancy(T_dataset, SAVE.figure_dir, CFG);

fprintf('\nTest 18 completed successfully.\n');
fprintf('Rows generated: %d\n', height(T_dataset));
fprintf('Regimes       : %s\n', strjoin(unique(T_dataset.field_regime_label)', ', '));
fprintf('Frequencies   : %s Hz\n', mat2str(unique(T_dataset.SIM_f0)'));
fprintf('SWS values    : %s m/s\n', mat2str(unique(T_dataset.SIM_cs_bg)'));
fprintf('dx=dz values  : %s mm\n', mat2str(unique(T_dataset.SIM_dx)'*1e3));
fprintf('REQ M values  : %s\n', mat2str(unique(T_dataset.REQ_M)'));
fprintf(['No baseline model was trained in Test 18. Training and old-vs-new ', ...
    'comparisons are intentionally deferred to Test 19.\n']);
fprintf('Dataset folder:\n%s\n', SAVE.output_dir);

%% Local functions

function T = build_condition_table(CFG)
rows = struct([]);
idx = 0;
for ri = 1:numel(CFG.FIELD_REGIMES)
    regime = CFG.FIELD_REGIMES(ri);
    for di = 1:numel(CFG.GRID.dx_dz)
        for ci = 1:numel(CFG.GRID.cs_bg)
            for fi = 1:numel(CFG.GRID.f0)
                idx = idx+1;
                rows(idx).condition_id = idx; %#ok<AGROW>
                rows(idx).regime_idx = ri;
                rows(idx).field_regime_label = regime.field_regime_label;
                rows(idx).field_regime_variant = regime.field_regime_variant;
                rows(idx).SIM_Nwaves = regime.Nwaves;
                rows(idx).SIM_Is2D = regime.Is2D;
                rows(idx).SIM_ForceInPlaneWave = regime.ForceInPlaneWave;
                rows(idx).SIM_f0 = CFG.GRID.f0(fi);
                rows(idx).SIM_cs_bg = CFG.GRID.cs_bg(ci);
                rows(idx).SIM_dx = CFG.GRID.dx_dz(di);
                rows(idx).SIM_dz = CFG.GRID.dx_dz(di);
            end
        end
    end
end
T = struct2table(rows);
end

function plot_representative_velocity(REP, fig_dir, CFG)
REP = one_per_regime(REP);
fig = figure('Color','w','Units','centimeters','Position',[2 2 27 7]);
tl = tiledlayout(1,numel(REP),'TileSpacing','compact','Padding','compact');
for i = 1:numel(REP)
    ax = nexttile(tl);
    imagesc(ax,REP(i).x*100,REP(i).z*100,real(REP(i).Uxz));
    axis(ax,'image'); set(ax,'YDir','normal','FontSize',9);
    colorbar(ax); colormap(ax,parula);
    title(ax,sprintf('%s | N=%d',REP(i).field_regime_label,REP(i).Nwaves), ...
        'Interpreter','none','FontWeight','normal');
    xlabel(ax,'x (cm)'); ylabel(ax,'z (cm)');
end
title(tl,'Representative clean particle-velocity fields','FontWeight','normal');
export_clean(fig,fullfile(fig_dir,'test18_representative_particle_velocity.png'),CFG);
end

function plot_representative_spectra(REP, fig_dir, CFG)
REP = one_per_regime(REP);
fig = figure('Color','w','Units','centimeters','Position',[2 2 27 7]);
tl = tiledlayout(1,numel(REP),'TileSpacing','compact','Padding','compact');
for i = 1:numel(REP)
    ax = nexttile(tl);
    imagesc(ax,REP(i).central_power_spectrum);
    axis(ax,'image'); set(ax,'YDir','normal','FontSize',9);
    colorbar(ax); colormap(ax,turbo);
    title(ax,sprintf('%s | N=%d',REP(i).field_regime_label,REP(i).Nwaves), ...
        'Interpreter','none','FontWeight','normal');
    xlabel(ax,'k_x bin'); ylabel(ax,'k_z bin');
end
title(tl,'Central-patch spectra, f=500 Hz, c_s=3 m/s, M=3, dx=0.5 mm', ...
    'FontWeight','normal');
export_clean(fig,fullfile(fig_dir,'test18_representative_power_spectra.png'),CFG);
end

function plot_q_distribution(T, fig_dir, CFG)
fig = figure('Color','w','Units','centimeters','Position',[2 2 19 11]);
boxchart(categorical(T.field_regime_label),T.q_theory, ...
    'GroupByColor',categorical(T.REQ_M));
grid on; ylim([0 1]);
xlabel('field regime'); ylabel('q target');
title('Test 18 q-target distribution by regime and REQ M','FontWeight','normal');
legend(compose('M=%g',unique(T.REQ_M)),'Location','best');
export_clean(fig,fullfile(fig_dir,'test18_q_distribution_by_regime.png'),CFG);
end

function plot_dataset_occupancy(T, fig_dir, CFG)
[G,S] = findgroups(T(:,{'field_regime_label','field_regime_variant'}));
S.N = splitapply(@numel,T.q_theory,G);
fig = figure('Color','w','Units','centimeters','Position',[2 2 20 11]);
bar(categorical(S.field_regime_variant),S.N);
grid on; ylabel('dataset rows'); xlabel('regime variant');
title('Test 18 dataset occupancy','FontWeight','normal');
xtickangle(35);
export_clean(fig,fullfile(fig_dir,'test18_dataset_occupancy.png'),CFG);
end

function out = one_per_regime(REP)
labels = ["directional_2D","diffuse_2D","partial_3D","diffuse_3D"];
indices = nan(numel(labels),1);
for i = 1:numel(labels)
    idx = find(arrayfun(@(x) x.field_regime_label == labels(i),REP),1);
    assert(~isempty(idx),'Missing representative for %s.',labels(i));
    indices(i) = idx;
end
out = REP(indices);
end

function T = remove_heavy_columns(T)
drop = false(1,width(T));
for i = 1:width(T)
    drop(i) = iscell(T.(T.Properties.VariableNames{i})) || ...
        isstruct(T.(T.Properties.VariableNames{i}));
end
T(:,drop) = [];
end

function T = concat_tables(A,B)
if isempty(A), T=B; return; end
if isempty(B), T=A; return; end
vars = unique([string(A.Properties.VariableNames), ...
    string(B.Properties.VariableNames)],'stable');
A = add_missing(A,vars); B = add_missing(B,vars);
T = [A(:,cellstr(vars)); B(:,cellstr(vars))];
end

function T = add_missing(T,vars)
for i = 1:numel(vars)
    if ismember(vars(i),string(T.Properties.VariableNames)), continue; end
    if any(contains(vars(i),["label","name","type","variant"]))
        T.(char(vars(i))) = strings(height(T),1);
    else
        T.(char(vars(i))) = nan(height(T),1);
    end
end
end

function make_dirs(SAVE)
dirs = string({SAVE.root_dir, SAVE.output_dir, SAVE.condition_dir, ...
    SAVE.data_dir, SAVE.table_dir, SAVE.figure_dir});
for i = 1:numel(dirs)
    if exist(dirs(i),'dir') ~= 7, mkdir(dirs(i)); end
end
end

function name = sanitize_filename(x)
name = regexprep(char(string(x)),'[^A-Za-z0-9_=-]+','_');
end

function export_clean(fig,path_i,CFG)
axs = findall(fig,'Type','axes');
for i = 1:numel(axs)
    try, axs(i).Toolbar.Visible='off'; catch, end
end
drawnow;
exportgraphics(fig,path_i,'Resolution',CFG.SAVE.png_resolution, ...
    'BackgroundColor','white');
close(fig);
end
