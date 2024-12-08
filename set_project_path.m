[PROJ_ROOT,~,~] = fileparts(mfilename('fullpath'));

% --- Set path ---
% First check if path is already set
s       = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr, [s, MODEL_SWEEP_DIR, s], 'IgnoreCase', ispc);

% Set if needed
if ~onPath
    fprintf('Setting up path with subfolders:\n\t%s\n',PROJ_ROOT)
    addpath(genpath(PROJ_ROOT))
end