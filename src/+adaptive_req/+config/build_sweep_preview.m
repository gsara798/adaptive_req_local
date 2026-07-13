function [T_axes, T_conditions, CFG_base] = build_sweep_preview(CFG)
%BUILD_SWEEP_PREVIEW Build a preview of the MC sweep configuration.
%
% This function does not run simulations.
%
% Inputs:
%   CFG:
%       Raw configuration struct loaded with load_profile_config.
%
% Outputs:
%   T_axes:
%       Table with the active sweep axes.
%
%   T_conditions:
%       Cartesian parameter matrix. Each row is one MC condition.
%
%   CFG_base:
%       Scalar baseline configuration created by taking the first value of
%       each swept parameter.
%
% Supported sweep value types:
%   numeric vectors
%   logical vectors
%   string arrays
%   character vectors
%   cell arrays of character vectors

if ~isfield(CFG, 'SWEEP') || ~isfield(CFG.SWEEP, 'paths')
    error('CFG must contain CFG.SWEEP.paths.');
end

if ~isfield(CFG.SWEEP, 'enabled') || ~logical(CFG.SWEEP.enabled)
    paths = strings(0, 1);
else
    paths = string(CFG.SWEEP.paths(:));
end

n_axes = numel(paths);

axis_values = cell(n_axes, 1);
axis_types = strings(n_axes, 1);
axis_names = strings(n_axes, 1);
n_values = zeros(n_axes, 1);
values_text = strings(n_axes, 1);

for i = 1:n_axes

    path_i = paths(i);
    raw_value_i = get_by_path(CFG, path_i);

    validate_sweep_values(path_i, raw_value_i);

    [values_i, type_i] = normalize_sweep_values(raw_value_i);

    axis_values{i} = values_i;
    axis_types(i) = type_i;
    axis_names(i) = path_to_variable_name(path_i);
    n_values(i) = numel(values_i);
    values_text(i) = values_to_text(values_i);

end

T_axes = table();

T_axes.axis_id = (1:n_axes).';
T_axes.path = paths;
T_axes.variable_name = axis_names;
T_axes.type = axis_types;
T_axes.n_values = n_values;
T_axes.values = values_text;

CFG_base = CFG;

for i = 1:n_axes
    first_value_i = axis_values{i}(1);
    CFG_base = set_by_path(CFG_base, paths(i), first_value_i);
end

if n_axes == 0

    T_conditions = table();
    T_conditions.condition_id = 1;
    return;

end

index_vectors = cell(n_axes, 1);

for i = 1:n_axes
    index_vectors{i} = 1:n_values(i);
end

index_grids = cell(n_axes, 1);
[index_grids{:}] = ndgrid(index_vectors{:});

n_conditions = numel(index_grids{1});

T_conditions = table();
T_conditions.condition_id = (1:n_conditions).';

for i = 1:n_axes

    var_name_i = char(axis_names(i));
    idx_i = index_grids{i}(:);
    values_i = axis_values{i};

    switch axis_types(i)

        case "numeric"
            T_conditions.(var_name_i) = values_i(idx_i);

        case "logical"
            T_conditions.(var_name_i) = values_i(idx_i);

        case "string"
            T_conditions.(var_name_i) = values_i(idx_i);

        otherwise
            error('Unsupported sweep axis type: %s', axis_types(i));

    end
end

end

% =========================================================================
% Local helper functions
% =========================================================================

function value = get_by_path(S, path)

parts = cellstr(split(string(path), '.'));
value = S;

for i = 1:numel(parts)

    field_i = parts{i};

    if ~isstruct(value) || ~isfield(value, field_i)
        error('Invalid configuration path: %s', string(path));
    end

    value = value.(field_i);

end

end

function S = set_by_path(S, path, value)

parts = cellstr(split(string(path), '.'));
S = set_by_parts(S, parts, value);

end

function S = set_by_parts(S, parts, value)

field_i = parts{1};

if numel(parts) == 1
    S.(field_i) = value;
else
    S.(field_i) = set_by_parts(S.(field_i), parts(2:end), value);
end

end

function validate_sweep_values(path, value)

is_valid_type = ...
    isnumeric(value) || ...
    islogical(value) || ...
    isstring(value) || ...
    ischar(value) || ...
    iscellstr(value);

if ~is_valid_type
    error(['Sweep path %s must contain numeric, logical, string, ', ...
           'char, or cellstr values.'], string(path));
end

if isempty(value)
    error('Sweep path %s contains an empty value.', string(path));
end

if ischar(value)
    return;
end

if ~isvector(value)
    error('Sweep path %s must contain a vector, not a matrix.', string(path));
end

end

function [values, value_type] = normalize_sweep_values(value)

if isnumeric(value)

    values = value(:);
    value_type = "numeric";

elseif islogical(value)

    values = value(:);
    value_type = "logical";

elseif isstring(value)

    values = value(:);
    value_type = "string";

elseif ischar(value)

    values = string(value);
    values = values(:);
    value_type = "string";

elseif iscellstr(value)

    values = string(value(:));
    value_type = "string";

else

    error('Unsupported sweep value type.');

end

end

function name = path_to_variable_name(path)

name = matlab.lang.makeValidName(strrep(char(path), '.', '_'));
name = string(name);

end

function txt = values_to_text(values)

if isnumeric(values) || islogical(values)

    txt = string(mat2str(values(:).'));

elseif isstring(values)

    txt = strjoin(values(:).', ", ");

else

    txt = string(values);

end

end