%% ======================================================================
%       MANUAL TEST RANGES FOR MJYOJIN KNOLL SWEEPS/MONTE CARLO
% =======================================================================
% C Rowell, Sep 2023

clear all; close all
% Uses Hajimirza conduit model, V6
codeDir = '~/code/research-projects/hydroVolc/';
addpath(genpath(codeDir))

%% Parameter input

% ------------------------------------------------------------------------
% 1st pair of input params to match the default as closely as possible
%   '-> Z0 is shortened to 2200, but N0+phi0 are adjusted to match output
%       of default run
% ------------------------------------------------------------------------
% MYOJIN PARAMS
conIn.Q              = 1e8;          % >108? range:
% conIn.conduit_radius = [];          % Unsure - use extrap val and run fit sweep
conIn.Z0             = 2200;        % Below vent or sea level? (Tsuru et al 2008)
% conIn.T              = 900+273.15;  % Assumed similar to Sumisu (Shukuno et al)
conIn.phi_frag       = 0.75;
conIn.rho_melt       = 2350;
conIn.Zw             = 0;         % (Still processing quench depths, but possibly up to 1km?) can run a range
conIn.N0             = 1.2e14;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn.phi0           = 0.3;           % unsure
% conIn.dP             = 0;           % Chamber overpressure

% At 2.2 km conduit length, phi0=31.4 and N0=1.198e14 are good starting
% values that correspond to the default conduit run...huh
%  For a Q range of say [1e8 to 1e9]:
    % 1) get a sweep range of params for dP, phi0, N0 @ Ze = 0 - adjusting radius as needed
        % Explore co-variance in the 3 params?
    % 2) explore ranges for Ze up to...800m?
    
    % Key outputs: 
        % fragmentation depth ranges (allow for adjustments to phi_frag to explore bubble OP peaks/dropoffs?)
        % total exsolved gas, estimated BSD (fixed at frag, yes? (though BSD likely unrealistic?)

% Other:
% vent altitude to be set at -Zw (below sea level)

% - default atmosphere (only for small pressure contribution of atmosphere)
% - default conduit friction coefficient

conIn.vh0 = -conIn.Zw; % Vent below sea level
conIn.conduit_radius = extrapVentRadius(conIn.Q)*1.041;

conIn2 = conIn;
conIn2.conduit_radius = conIn.conduit_radius*0.837; % Tweaked manually to be about right for choking
conIn2.composition = struct(... % Weight fractions
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

% ========== Single run tests ================
% ------------------------------------------------------------------------
% Try with no initial exsolved gas/bubbles but enhanced chamber
% overpressure
conIn3 = conIn2;
conIn3.N0             = 0;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn3.phi0           = 0;           % unsure
conIn3.dP             = 3.38e7;           % Chamber overpressure
% ------------------------------------------------------------------------
conIn4 = conIn2;
conIn4.N0             = 0;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn4.phi0           = 0;           % unsure
conIn4.Q              = 1e9;          % >108? range:
conIn4.conduit_radius = extrapVentRadius(conIn4.Q)*1.4;  % Unsure - use extrap val and run fit sweep

% ========= Initial sweep params ============
% Starting from reasonable core scenario of conIn2
%   - sweep dP and (covary?) phi0,N0?
%   - eval fragmenting and non-fragmenting valid results

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

%% 

D(1).cI = getConduitSource(conIn);
D(2).cI = getConduitSource(conIn2);
D(3).cI = getConduitSource(conIn3);
D(4).cI = getConduitSource(conIn4);

% Run a single test pair for composition - works?
% D(1).cO = Conduit_flow_with_nucleation_V6(D(1).cI); % Run v6 control
% D(2).cO = Conduit_flow_with_nucleation_V7(D(2).cI); % Run v7 with standard composition
% D(3).cO = Conduit_flow_with_nucleation_V7(D(3).cI); % Run v7 with standard composition
% D(4).cO = Conduit_flow_with_nucleation_V7(D(4).cI); % Run v7 with standard composition

% plotConduitOutput(D)
% [ax,th] = plotConduitObjective(D);
%% Run some simple sweep tests, output initial plots of result + validity

radius = conIn3.conduit_radius .* (0.96:0.002:0.992);
% radius = conIn3.conduit_radius .* (0.98:0.001:0.992);

ax = [];
for ri = length(radius):-1:1
    dd(ri).cI = conIn3;
    dd(ri).cI.conduit_radius = radius(ri);
    dd(ri).cI = getConduitSource(dd(ri).cI);
    dd(ri).cO = Conduit_flow_with_nucleation_V7(dd(ri).cI);
    
%     [ax] = plotConduitObjective(dd(ri),[],ax);
end
[ax] = plotConduitObjective(dd,'conduit_radius');
% plotConduitOutput(dd(10:13))