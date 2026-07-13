function CFG = load_profile_config(profile_name, varargin)
%LOAD_PROFILE_CONFIG Load a raw configuration profile.
%
% This function loads the default high-level configuration and then applies
% a user profile from the configs/ folder.
%
% Unlike load_run_config, this function does not resolve the configuration
% into a single run. Therefore, it can safely load sweep configurations
% containing vector-valued fields.

p = inputParser;
p.FunctionName = 'adaptive_req.config.load_profile_config';

addRequired(p, 'profile_name', @(x) ischar(x) || isstring(x));

addParameter(p, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'ConfigDir', '', @(x) ischar(x) || isstring(x));

parse(p, profile_name, varargin{:});

profile_name = char(p.Results.profile_name);
root_dir = char(p.Results.RootDir);
config_dir = char(p.Results.ConfigDir);

if isempty(config_dir)
    config_dir = fullfile(root_dir, 'configs');
end

if exist(config_dir, 'dir')
    addpath(config_dir);
else
    error('Configuration folder not found: %s', config_dir);
end

CFG = adaptive_req.config.default_run_config();

if ~strcmpi(profile_name, 'default')

    if exist(profile_name, 'file') ~= 2
        error('Configuration profile not found: %s', profile_name);
    end

    profile_fun = str2func(profile_name);
    CFG = profile_fun(CFG);
end

CFG.EXP.profile_name = profile_name;
CFG.EXP.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));

end