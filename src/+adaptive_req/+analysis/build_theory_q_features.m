function [T_out, T_cache] = build_theory_q_features(T, varargin)
%BUILD_THEORY_Q_FEATURES Add discrete theory-q candidate predictors.
%
% The calculation is cached by the REQ/simulation geometry so repeated
% windows from the same condition do not recompute the same theory curves.

p = inputParser;
p.FunctionName = 'adaptive_req.analysis.build_theory_q_features';
addRequired(p, 'T', @istable);
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, T, varargin{:});

verbose = logical(p.Results.Verbose);
T_out = T;

required = ["SIM_dx", "SIM_dz", "SIM_f0", "REQ_cs_guess", ...
    "REQ_M", "REQ_gamma_win", "REQ_pad_factor", "REQ_smooth_sigma"];
missing = required(~ismember(required, string(T.Properties.VariableNames)));
assert(isempty(missing), 'Missing variables for theory q features: %s', ...
    strjoin(missing, ', '));

q_names = [
    "q_theory_dir2D"
    "q_theory_diffuse2D"
    "q_theory_projected3D"
    "q_theory_mean_dir2D_projected3D"
    "q_theory_mean_all"];

for i = 1:numel(q_names)
    if ~ismember(q_names(i), string(T_out.Properties.VariableNames))
        T_out.(char(q_names(i))) = nan(height(T_out), 1);
    end
end

key_vars = ["SIM_dx", "SIM_dz", "SIM_f0", "REQ_cs_guess", ...
    "REQ_M", "REQ_gamma_win", "REQ_pad_factor", "REQ_smooth_sigma"];
if ismember("REQ_Nbins_requested", string(T.Properties.VariableNames))
    key_vars(end + 1) = "REQ_Nbins_requested";
end

[G, T_cache] = findgroups(T(:, cellstr(key_vars)));
T_cache.q_theory_dir2D = nan(height(T_cache), 1);
T_cache.q_theory_diffuse2D = nan(height(T_cache), 1);
T_cache.q_theory_projected3D = nan(height(T_cache), 1);
T_cache.q_theory_mean_dir2D_projected3D = nan(height(T_cache), 1);
T_cache.q_theory_mean_all = nan(height(T_cache), 1);

if verbose
    fprintf('Computing discrete theory-q candidates for %d unique REQ geometries...\n', ...
        height(T_cache));
end

for gi = 1:height(T_cache)
    rows_i = find(G == gi);
    r = rows_i(1);

    nbins = 'auto';
    if ismember("REQ_Nbins_requested", string(T.Properties.VariableNames))
        nbins_i = T.REQ_Nbins_requested(r);
        if isnumeric(nbins_i) && isfinite(nbins_i)
            nbins = nbins_i;
        end
    end

    common = { ...
        'M', T.REQ_M(r), ...
        'Gamma', T.REQ_gamma_win(r), ...
        'PadFactor', T.REQ_pad_factor(r), ...
        'Nbins', nbins, ...
        'SmoothSigma', T.REQ_smooth_sigma(r), ...
        'TheoryMode', 'S2D', ...
        'Plot', false, ...
        'UseDonutFilter', false};

    q_dir = compute_one(T.SIM_dx(r), T.SIM_dz(r), T.SIM_f0(r), ...
        T.REQ_cs_guess(r), 'SingleWave', common);
    q_diff2 = compute_one(T.SIM_dx(r), T.SIM_dz(r), T.SIM_f0(r), ...
        T.REQ_cs_guess(r), 'Diffuse2D', common);
    q_proj3 = compute_one(T.SIM_dx(r), T.SIM_dz(r), T.SIM_f0(r), ...
        T.REQ_cs_guess(r), 'Diffuse3D', common);

    T_cache.q_theory_dir2D(gi) = q_dir;
    T_cache.q_theory_diffuse2D(gi) = q_diff2;
    T_cache.q_theory_projected3D(gi) = q_proj3;
    T_cache.q_theory_mean_dir2D_projected3D(gi) = mean([q_dir q_proj3], 'omitnan');
    T_cache.q_theory_mean_all(gi) = mean([q_dir q_diff2 q_proj3], 'omitnan');

    T_out.q_theory_dir2D(rows_i) = T_cache.q_theory_dir2D(gi);
    T_out.q_theory_diffuse2D(rows_i) = T_cache.q_theory_diffuse2D(gi);
    T_out.q_theory_projected3D(rows_i) = T_cache.q_theory_projected3D(gi);
    T_out.q_theory_mean_dir2D_projected3D(rows_i) = ...
        T_cache.q_theory_mean_dir2D_projected3D(gi);
    T_out.q_theory_mean_all(rows_i) = T_cache.q_theory_mean_all(gi);
end

end

function q = compute_one(dx, dz, f0, cs_guess, field_type, common)

out = adaptive_req.theory.q_theory_REQ_discrete_shearUZ( ...
    dx, dz, f0, cs_guess, common{:}, 'FieldType', field_type);
q = out.q_th;

end
