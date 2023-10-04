%% ======================================================================
%            SETUP FOR MJYOJIN KNOLL SWEEPS/MONTE CARLO
% =======================================================================
% C Rowell, Sep 2023

clear all; close all
% Uses Hajimirza conduit model, V6
codeDir = '~/code/research-projects/hydroVolc/';
addpath(genpath(codeDir))

%% Parameter input

conIn.Q = 1e8;
conIn.Zw = 100;
conIn.Z0 = 2200;

% MYOJIN PARAMS
% conIn.Q              = 1e8;          % >108? range:
% conIn.conduit_radius = [];          % Unsure - use extrap val and run fit sweep
% conIn.Z0             = 2200;        % Below vent or sea level? (Tsuru et al 2008)
% conIn.T              = 900+273.15;  % Assumed similar to Sumisu (Shukuno et al)
% conIn.phi_frag       = 0.75;
% conIn.rho_melt       = 2350;
% conIn.Zw             = 100;         % (Still processing quench depths, but possibly up to 1km?) can run a range
conIn.N0             = 1.2e14;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn.phi0           = 0.3;           % unsure
% conIn.dP             = 0;           % Chamber overpressure

% At 2.2 km conduit length, phi0=31.4 and N0=1.198e14 are good starting
% values that correspond to the default run...huh

% Other:
% vent altitude to be set at -Zw (below sea level)

% - default atmosphere (only for small pressure contribution of atmosphere)
% - default conduit friction coefficient
 

% COMPOSITION - NOT YET APPLIED
mjoyin.composition = Composition(...
    'SiO2', 74.53e-2, ...
    'TiO2', .32e-2, ...
    'Al2O3', 13.50e-2, ...
    'FeO',  1.94e-2, ...    % combine?
    'MnO',  .11e-2, ...
    'MgO',  .61e-2, ...
    'CaO',  2.74e-2, ...
    'Na2O', 4.41e-2, ...
    'K2O',  0.81e-2 ...
    );

% ----------- Search params -------------
% Define conduit range params to search for shooting solutions
aRange = [0.5 1.5]; % Fractional range of conduit radius to search

% Conduit run failure thresholds
ZfailScale  = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .04;  % Over-/underpressure threshold

% Tolerances
dRminScale  = 0.2e-3;
dQminScale  = 0.2e-3;
maxIter     = 10;

% Original workflow switches - m
getBaseRuns = false; % Run conduit with no water depth to get base R?
getQfromR   = true;  % Run with adjusted MER?

%% Start with a single run, put into a function later
conIn.vh0 = -conIn.Zw; % Vent below sea level
conIn.conduit_radius = extrapVentRadius(conIn.Q);
Zfailthresh = ZfailScale*conIn.conduit_radius;

% Do a first test run
conIn = getConduitSource(conIn);
conOut = Conduit_flow_with_nucleation_V6(conIn);

% Check for a valid result
[~,~,~,~,~,~,valid] = checkConduitResult(conOut,Zfailthresh,Mfailthresh,Pfailthresh);
if valid
    Rvalid = conIn.conduit_radius;
else
    Rvalid = [];
end    

[Rrange,cIi,conOut,success_base] = conduitRadiusFromQ(conIn,conIn.conduit_radius.*aRange,...
    'Rvalid',Rvalid,'maxIter',maxIter,...
    'Zfailthresh',Zfailthresh,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
    'dRminScale',dRminScale,'verbose',true,'output',true);


D.cI = cIi; D.cO = conOut;
plotConduitOutput(D)
