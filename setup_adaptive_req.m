function root_dir = setup_adaptive_req()
%SETUP_ADAPTIVE_REQ Add the adaptive_req source folder to the MATLAB path.
%
% Usage:
%   root_dir = setup_adaptive_req();

root_dir = fileparts(mfilename('fullpath'));

src_dir = fullfile(root_dir, 'src');

if ~exist(src_dir, 'dir')
    error('Source folder not found: %s', src_dir);
end

addpath(src_dir);
rehash toolboxcache;

fprintf('adaptive_req setup completed.\n');
fprintf('Project root: %s\n', root_dir);
fprintf('Source path : %s\n', src_dir);

end