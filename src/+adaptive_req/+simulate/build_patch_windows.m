function patch_pack = build_patch_windows(cfg, feat_cfg, varargin)
%BUILD_PATCH_WINDOWS Build patch index lists for local spectral analysis.
%
% Usage:
%   patch_pack = adaptive_req.simulate.build_patch_windows(cfg, feat_cfg);
%
%   patch_pack = adaptive_req.simulate.build_patch_windows(cfg, feat_cfg, ...
%       'NumPatches', 9);
%
%   patch_pack = adaptive_req.simulate.build_patch_windows(cfg, feat_cfg, ...
%       'GridSize', [5 5]);
%
%   patch_pack = adaptive_req.simulate.build_patch_windows(cfg, feat_cfg, ...
%       'Pattern', 'center');
%
%   patch_pack = adaptive_req.simulate.build_patch_windows(cfg, feat_cfg, ...
%       'Pattern', 'manual', ...
%       'ManualCenters', [250 250; 200 250; 300 250]);
%
% Output convention:
%   Uxz(z,x)
%   x_idx_list contains x indices
%   z_idx_list contains z indices

p = inputParser;

addRequired(p, 'cfg', @isstruct);
addRequired(p, 'feat_cfg', @isstruct);

addParameter(p, 'Pattern', 'grid', @(x) ischar(x) || isstring(x));
addParameter(p, 'NumPatches', 9, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'GridSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'CoverageFraction', 0.65, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'ManualCenters', [], @(x) isempty(x) || isnumeric(x));

parse(p, cfg, feat_cfg, varargin{:});

pattern = lower(char(p.Results.Pattern));
num_patches_requested = round(p.Results.NumPatches);
grid_size = p.Results.GridSize;
coverage_fraction = p.Results.CoverageFraction;
manual_centers = p.Results.ManualCenters;

validate_required_fields(cfg, {'Lx', 'Lz', 'dx', 'dz'});
validate_required_fields(feat_cfg, {'win_size'});

Nx = round(cfg.Lx / cfg.dx) + 1;
Nz = round(cfg.Lz / cfg.dz) + 1;

win_size = round(feat_cfg.win_size);

if mod(win_size, 2) == 0
    win_size = win_size + 1;
end

half_win = floor(win_size / 2);

x_min_valid = 1 + half_win;
x_max_valid = Nx - half_win;

z_min_valid = 1 + half_win;
z_max_valid = Nz - half_win;

if x_min_valid > x_max_valid || z_min_valid > z_max_valid
    error('Window size is too large for the simulation domain.');
end

switch pattern

    case 'center'
        cx_list = round(Nx / 2);
        cz_list = round(Nz / 2);
        patch_labels = "patch_001";

    case {'grid', 'auto'}
        if isempty(grid_size)
            [n_rows, n_cols] = choose_grid_size(num_patches_requested);
        else
            n_rows = round(grid_size(1));
            n_cols = round(grid_size(2));

            if n_rows < 1 || n_cols < 1
                error('GridSize must contain positive integers.');
            end

            num_patches_requested = n_rows * n_cols;
        end

        [cx_grid, cz_grid] = make_center_grid( ...
            x_min_valid, x_max_valid, ...
            z_min_valid, z_max_valid, ...
            n_rows, n_cols, ...
            coverage_fraction);

        cx_all = cx_grid(:);
        cz_all = cz_grid(:);

        if numel(cx_all) > num_patches_requested
            center_x = round(Nx / 2);
            center_z = round(Nz / 2);

            dist2 = (cx_all - center_x).^2 + (cz_all - center_z).^2;
            [~, order] = sort(dist2, 'ascend');

            keep = sort(order(1:num_patches_requested), 'ascend');

            cx_all = cx_all(keep);
            cz_all = cz_all(keep);
        end

        cx_list = cx_all(:).';
        cz_list = cz_all(:).';

        patch_labels = strings(1, numel(cx_list));
        for i = 1:numel(cx_list)
            patch_labels(i) = "patch_" + sprintf('%03d', i);
        end

    case 'manual'
        if isempty(manual_centers)
            error('ManualCenters must be provided when Pattern is manual.');
        end

        if size(manual_centers, 2) ~= 2
            error('ManualCenters must be an N by 2 matrix: [cx, cz].');
        end

        cx_list = round(manual_centers(:, 1).');
        cz_list = round(manual_centers(:, 2).');

        patch_labels = strings(1, numel(cx_list));
        for i = 1:numel(cx_list)
            patch_labels(i) = "manual_" + sprintf('%03d', i);
        end

    case 'nine'
        warning(['Pattern = ''nine'' is deprecated. ', ...
                 'Use Pattern = ''grid'', NumPatches = 9 instead.']);

        [cx_grid, cz_grid] = make_center_grid( ...
            x_min_valid, x_max_valid, ...
            z_min_valid, z_max_valid, ...
            3, 3, ...
            coverage_fraction);

        cx_list = cx_grid(:).';
        cz_list = cz_grid(:).';

        patch_labels = strings(1, numel(cx_list));
        for i = 1:numel(cx_list)
            patch_labels(i) = "patch_" + sprintf('%03d', i);
        end

    otherwise
        error('Unknown patch pattern: %s', pattern);
end

n_patches = numel(cx_list);

x_idx_list = cell(1, n_patches);
z_idx_list = cell(1, n_patches);

for pidx = 1:n_patches

    x_idx = (cx_list(pidx) - half_win):(cx_list(pidx) + half_win);
    z_idx = (cz_list(pidx) - half_win):(cz_list(pidx) + half_win);

    if min(x_idx) < 1 || max(x_idx) > Nx
        error('Patch %d exceeds x bounds. Reduce win_size or CoverageFraction.', pidx);
    end

    if min(z_idx) < 1 || max(z_idx) > Nz
        error('Patch %d exceeds z bounds. Reduce win_size or CoverageFraction.', pidx);
    end

    x_idx_list{pidx} = x_idx;
    z_idx_list{pidx} = z_idx;
end

patch_pack = struct();

patch_pack.Nx = Nx;
patch_pack.Nz = Nz;

patch_pack.win_size = win_size;
patch_pack.half_win = half_win;

patch_pack.cx_list = cx_list;
patch_pack.cz_list = cz_list;

patch_pack.x_idx_list = x_idx_list;
patch_pack.z_idx_list = z_idx_list;

patch_pack.patch_labels = patch_labels;
patch_pack.n_patches = n_patches;

patch_pack.pattern = pattern;
patch_pack.num_patches_requested = num_patches_requested;
patch_pack.coverage_fraction = coverage_fraction;

patch_pack.x_min_valid = x_min_valid;
patch_pack.x_max_valid = x_max_valid;
patch_pack.z_min_valid = z_min_valid;
patch_pack.z_max_valid = z_max_valid;

end

function [n_rows, n_cols] = choose_grid_size(num_patches)

n_cols = ceil(sqrt(num_patches));
n_rows = ceil(num_patches / n_cols);

end

function [cx_grid, cz_grid] = make_center_grid( ...
    x_min_valid, x_max_valid, ...
    z_min_valid, z_max_valid, ...
    n_rows, n_cols, ...
    coverage_fraction)

x_center = 0.5 * (x_min_valid + x_max_valid);
z_center = 0.5 * (z_min_valid + z_max_valid);

x_half_span = 0.5 * (x_max_valid - x_min_valid) * coverage_fraction;
z_half_span = 0.5 * (z_max_valid - z_min_valid) * coverage_fraction;

if n_cols == 1
    cx_vals = round(x_center);
else
    cx_vals = round(linspace( ...
        x_center - x_half_span, ...
        x_center + x_half_span, ...
        n_cols));
end

if n_rows == 1
    cz_vals = round(z_center);
else
    cz_vals = round(linspace( ...
        z_center - z_half_span, ...
        z_center + z_half_span, ...
        n_rows));
end

[cx_grid, cz_grid] = meshgrid(cx_vals, cz_vals);

end

function validate_required_fields(s, names)

for i = 1:numel(names)
    if ~isfield(s, names{i})
        error('Missing required field: %s', names{i});
    end
end

end