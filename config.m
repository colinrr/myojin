% Configuring matlab directories, path, other useful variables

% THIS NEEDS CONFIGURING ON THE LOCAL SYSTEM
%  - Directory containing the large model parameter sweep data files 
%   (not included in git repository)
MAIN_DATA_DIR = '/Users/crrowell/Kahuna/data/myojin/';


% ---- set other directories ----
[PROJ_ROOT,~,~] = fileparts(mfilename('fullpath'));

DATA_DIR        = fullfile(PROJ_ROOT,'data');
FIGURES_DIR     = fullfile(PROJ_ROOT,'figures');
SCRIPTS_DIR     = fullfile(PROJ_ROOT,'sandbox_scripts');
MODEL_SWEEP_DIR = fullfile(PROJ_ROOT,'conduit_sweep_scripts');


