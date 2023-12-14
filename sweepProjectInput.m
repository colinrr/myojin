%% ======================================================================
%            SETUP FOR MJYOJIN KNOLL SWEEPS/MONTE CARLO
% =======================================================================
% C Rowell, Sep 2023

clear all; close all
% Uses Hajimirza conduit model, V6
codeDir = '~/code/research-projects/hydroVolc/';
addpath(genpath(codeDir))

%% THE PLAN

% Manual baselines?
% Zw = either 0 or min water depth?
% N0, phi0, dP = 0

% parameter sweep sets:

% Do for: Q = [1e8 1e9]
% Do for: N0 = [1e14 1e15] - but see con model papers for BSD's and approp. vals
% `-> Gives minimum 4 run sets, plus some baseline scenarios with any
% of N0, phi0, dP = 0;
%
% For each set, vary:
%   dP: 0 - ~40 MPa? catch errors in sat pressure, maybe calc upper lim?
%   n_ex0: 0-50% excess exsolved gas mass, use to calc phi_range
%   a: always allow to vary as needed for solution
%
% THEN runs for Zw = min to max by 10 to 50 m step size
%
% So 3D solution search in dP, n_ex0/phi0, Zw
%   - may need to consider how strongly dP is f'n of n_ex0, Zw - plot P components vs Zw
%   - Zw may actually map onto dP in a simple way to reduce dimensionality?
%
% Need:
%   Output behavior interpreter

%% Parameter input

% conIn.conduit_radius = [];          % Unsure - use extrap val and run fit sweep

% Shallow chamber depth given as 1.4-2.2 km bsl (Tsuru et al 2008)
%   '-> They further suggest this would result from volumetric expansion at
%   shallow depth, presumably if, say, the magma stalled there for some short time
%   '-> Follow by explosive eruption
%   '-> They further suggest the caldera is explosive in origin (NOT
%   collapse)
% Fiske et al 2001 gives estimated pre-caldera summit depth of 600 m
%  '-> but shoaled to within about 300m bsl during eruption early stages

max_chamber_depth           = 2200; % Tsuru ea argue this is chamber BOTTOM
current_caldera_floor_depth = 1400;  % -'-> and that this is potential chamber top
min_vent_depth              = 300;
max_vent_depth              = 1300; % Seems like a reasonable absolute max
% So MAX Z0 is 1900, ranging down to just a few hundred (assuming
% progressive vent deepening...)

% SET UP INITIAL COARSE RUNS -------
% -> Do these plus two controls runs with chamber at 6 km
Zw = [0 min_vent_depth:200:max_vent_depth]';
chamber_depths = [current_caldera_floor_depth max_chamber_depth];

% Start sweep list with 2 control runs:
sweepList.Zw = [0; 300; Zw; Zw];
sweepList.Z0 = [6e3; 6e3-300; current_caldera_floor_depth-Zw; max_chamber_depth-Zw];

% For each run, sweep these parameters:
sweepParams.dP.range = [0 4e6]; % 0 - 40 MPa overpressures
sweepParams.dP.n     = 21;

sweepParams.n0_excess.range = [0 0.7];
sweepParams.n0_excess.n     = 16;

% Do for:
Q = [1e8 1e9];
N0 = [1e14 1e15]; % Fix these initial BND's roughly to Q based on early conduit results
                   % '-> (see N_vs_phi_checks.m and bnd_h20_MLRFits.m)
%% MYOJIN COMMON PARAMS ------
conIn.Z0             = 2200;        % Below vent or sea level? 
conIn.Zw             = 0;         % (Still processing quench depths, but possibly up to 1km?) can run a range
% conIn.T              = 900+273.15;  % Assumed similar to Sumisu (Shukuno et al)
conIn.phi_frag       = 0.75;
conIn.rho_melt       = 2350;
conIn.vh0            = -conIn.Zw; % Vent below sea level
conIn.conduit_radius = 39.66; % based on initial runs;

conIn.N0             = 1.2e14;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn.phi0           = 0.3;           % unsure
% conIn.dP             = 0;           % Chamber overpressure

% Conduit run failure thresholds: DEFAULTS are fine for now

conIn.composition = struct(... % Weight fractions
    'SiO2', 74.53e-2, ...
    'TiO2', .32e-2, ...
    'Al2O3', 13.50e-2, ...
    'FeO',  1.94e-2, ...    % combine w/ Fe2O3? Should check Hui and Zhang
    'Fe2O3', 0.94e-2, ...
    'MnO',  .11e-2, ...
    'MgO',  .61e-2, ...
    'CaO',  2.74e-2, ...
    'Na2O', 4.41e-2, ...
    'K2O',  0.81e-2 ...
    );

% ------- SWEEP PARAMS -------

sweepParams.phi0 = []; % Crudely - get this range from manual tests?
sweepParams.dP = [];
% sweepParams.N0 = f(phi0)

% At 2.2 km conduit length, phi0=31.4 and N0=1.198e14 are good starting
% values that correspond to the default run...huh

% Other:
% vent altitude to be set at -Zw (below sea level)
% conduit depth to be set at (abs chamber depth (bsl) - vent depth (bsl))

% - default atmosphere (only for small pressure contribution of atmosphere)
% - default conduit friction coefficient
 

% ----------- Search params -------------
% Define conduit range params to search for shooting solutions
aRange = [0.5 1.5]; % Fractional range of conduit radius to search?

% Search Tolerances
dRminScale  = 0.2e-3;
dQminScale  = 0.2e-3;
maxIter     = 10;

% Original workflow switches - m
getBaseRuns = false; % Run conduit with no water depth to get base R?
getQfromR   = true;  % Run with adjusted MER?

%%


%% OLD CRAP: Start with a single run, put into a function later
% conIn.vh0 = -conIn.Zw; % Vent below sea level
% conIn.conduit_radius = extrapVentRadius(conIn.Q);
% Zfailthresh = ZfailScale*conIn.conduit_radius;
% 
% % Do a first test run
% conIn = getConduitSource(conIn);
% conOut = Conduit_flow_with_nucleation_V6(conIn);
% 
% % Check for a valid result
% [~,~,~,~,~,~,valid] = checkConduitResult(conOut,Zfailthresh,Mfailthresh,Pfailthresh);
% if valid
%     Rvalid = conIn.conduit_radius;
% else
%     Rvalid = [];
% end    
% 
% [Rrange,cIi,conOut,success_base] = conduitRadiusFromQ(conIn,conIn.conduit_radius.*aRange,...
%     'Rvalid',Rvalid,'maxIter',maxIter,...
%     'Zfailthresh',Zfailthresh,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
%     'dRminScale',dRminScale,'verbose',true,'output',true);
% 
% 
% D.cI = cIi; D.cO = conOut;
% plotConduitOutput(D)
