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


% MYOJIN PARAMS
conIn.Q              = 1e8;          % >108? range:
% conIn.conduit_radius = [];          % Unsure - use extrap val and run fit sweep
conIn.Z0             = 2200;        % Below vent or sea level? (Tsuru et al 2008)
% conIn.T              = 900+273.15;  % Assumed similar to Sumisu (Shukuno et al)
conIn.phi_frag       = 0.75;
conIn.rho_melt       = 2350;
conIn.Zw             = 0;         % (Still processing quench depths, but possibly up to 1km?) can run a range
conIn.vh0            = -conIn.Zw; % Vent below sea level
conIn.conduit_radius = 39.66; % based on initial runs;

conIn.N0             = 1.2e14;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
conIn.phi0           = 0.3;           % unsure
% conIn.dP             = 0;           % Chamber overpressure

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
% ----------- Search params -------------
% Define conduit range params to search for shooting solutions
qRange = [1e8 1e9];
% aRange = extrapVentRadius([1e8 1e9]) .* 0.85; % Absolute radii initial guess for each Q0
aVals = [38.7 70];
% aRange = [0.836 0.8642]; % Fractional range of conduit radius to search based on initial 1e8 runs
aRange = [0.97 1.03];
na = 21;
nq = length(qRange);
aScale = linspace(aRange(1),aRange(2),na); % These may need to change for each run set

% Conduit run failure tolerances
conIn.zFailTol  = 2;   % Z (depth) threshold in conduit radii - RELAXED FOR NOW to accommodate near-frag edge cases
% conIn.mFailTol = 0.95; % Mach number threshold
% conIn.pFailTol = .001;  % Over-/underpressure relative tolerance

% Search Step Tolerances - NOT YET ON USE
% dRminScale  = 0.2e-3;
% dQminScale  = 0.2e-3;
% maxIter     = 10;

% ALTERNATE FRAG'N THRESHOLDS
dPg = 5e6; % Estimated bubble overpressure threshold, Pa
epsilon = []; % Extensional strain rate threshold?
% Shear strain rate threshold?


% Run set, for each: run Q = [1e8, 1e9]
run_set(1) = false;    % 1) Base set, no excess P or gas
run_set(2) = true;    % 2) Excess P to reach surface
run_set(3) = true;    % 3) Excess gas up to ~50-65%, using 2.P0
run_set(4) = true;    % 4) Vary P0 from best 3.N0 result?
run_set(5) = false;    % 5) Given success values in 4, run some reprentative Ze ranges
run_set(6) = true;     % 6) Run radius variation tests (good for testing conduitRadiusFromQ)

plotObj = true;

% Original workflow switches - not in use yet but something similar by
% variable may happen
% getBaseRuns = false; % Run conduit with no water depth to get base R?
% getQfromR   = true;  % Run with adjusted MER?

% STAGED ADJUSTMENTS
% conIn2 = conIn;
% conIn2.conduit_radius = conIn.conduit_radius*0.837; % Tweaked manually to be about right for choking

% ========== Short sweep tests ================
% 1) Base set, Myojin composition, no excess P or gas -----------------------------------------
cI(1) = getConduitSource(conIn);
cI(1).N0             = 0;           
cI(1).phi0           = 0;          
cI(1).dP             = 0;    % Chamber overpressure
desc1 = 'Base set, Myojin composition, no excess P or gas';
% plotsets1 = []; % Plot all
for qi = nq:-1:1
    dat1(qi).cI = cI(1);
    dat1(qi).cI.Q = qRange(qi);
%     dat1(qi).cI.conduit_radius = extrapVentRadius(qRange(qi)).*mean(aRange); % Just 1 each for now
    dat1(qi).cI.conduit_radius = aVals(qi).*mean(aRange); % Just 1 each for now
    dat1(qi).cI = getConduitSource(dat1(qi).cI);
end
 % -->OUTCOME: Buoyancy alone is not enough to overcome pressure losses,
 % results in no frag'n and stall/intrusion somewhere around 1km depth

% 2) Excess P to reach surface --------------------------------------------
% Try with no initial exsolved gas/bubbles but enhanced chamber
% overpressure
desc2  = 'Excess P to reach surface, plus Effusive R search test';
plotsets2 = [3 6 9 12 15 ; 7 10 13 16 19]; % Center index is "best" result for now
cI(2) = getConduitSource(conIn);
cI(2).N0             = 0;           
cI(2).phi0           = 0;          
cI(2).dP             = 3.38e7;    % Chamber overpressure
% radius = cI(2).conduit_radius .* (0.96:0.002:0.992); % Choose param range
for ri = na:-1:1
    for qi = nq:-1:1
        dat2(ri,qi).cI                  = cI(2);
        dat2(ri,qi).cI.Q                = qRange(qi);
%         dat2(ri,qi).cI.conduit_radius = extrapVentRadius(qRange(qi)).*aScale(ri);
        dat2(ri,qi).cI.conduit_radius   = aVals(qi).*aScale(ri);
        dat2(ri,qi).cI                  = getConduitSource(dat2(ri,qi).cI);
    end
end

 % -->OUTCOME: This set results in runs with fragmentation condition very close to
  % Z=0, resulting in somewhat ambiguous conditions for "valid" runs
    % -> 75% porosity occurs within a radius or two of surface, velocity is
    % so pressures are "somewhat" unbalanced and whether or not the model
    % should proceed to a very brief fragmentation step is not clear.
    % -> Could be considered valid and effusive if Pm is within a couple
    % bars of atmospheric, based on Sahand's original conditions for frag'n
    % -> Overpressure and porosity are very high, however, so frag is very
    % possible at or near vent?
    
% 3) Excess gas up to ~50-65%, using 2.P0 ---------------------------------
desc3  = 'Excess gas up to ~50-65%, using 2.P0';
plotsets3 = [8 10 12 14 16 ; 7 9 11 13 15]; % Center index is "best" result for now
cI(3) = getConduitSource(conIn);        
conduit_radius = [38.35 70]; %[39 85]; % [35.77 78];
% cI(3).dP             = 1e7; %3.38e7;    % Chamber overpressure
dP                  = [1e7 2e7];
N0_range            = [1e14 1e15];  % Largely MER dependent in initial tests, so tie to Q for now
excess_gas_range    = [0.0 0.65]; % Fraction initial exsolved volatiles in EXCESS of solubility
excess_gas_vals     = linspace(excess_gas_range(1),excess_gas_range(2),na);
for pi = na:-1:1
    for qi = nq:-1:1
        dat3(pi,qi).cI                  = cI(3);
        dat3(pi,qi).cI.Q                = qRange(qi);
        dat3(pi,qi).cI.conduit_radius   = conduit_radius(qi);
        if excess_gas_vals(pi) ==0
            dat3(pi,qi).cI.N0               = 0;
        else
            dat3(pi,qi).cI.N0               = N0_range(qi);
        end
        dat3(pi,qi).cI.dP               = dP(qi);
        dat3(pi,qi).cI.phi0             = h2oExcessMass2ReservoirPorosity(excess_gas_vals(pi),dat3(pi,qi).cI);
        dat3(pi,qi).cI                  = getConduitSource(dat3(pi,qi).cI);
    end
end

% 4) Vary P0 from best 3.N0 result  --------------------------------------
desc4  = 'Vary P0 from best 3.N0 result';
plotsets4 = [7:2:15 ; 7:2:15]; % Center index is "best" result for now
cI(4) = getConduitSource(dat3(12,1).cI); % Best 10^8 result from 4       
conduit_radius      = [39 85];
% cI(3).dP             = 1e7; %3.38e7;   
N0_range            = [1e14 1e15];  % Largely MER dependent in initial tests, so tie to Q for now
excess_gas_range    = [0.3665 0.3350]; % EXCESS exsolved volatiles, best from (3)
dP_range             = [5e6 2.5e7];  % Chamber overpressure range
dP                  = [1e7 2e7];
dP_mult             = linspace(0.25,1.25,na)';
dP_vals             = dP.*dP_mult;  %linspace(dP_range(1),dP_range(2),na);
for pi = na:-1:1
    for qi = nq:-1:1
        dat4(pi,qi).cI                  = cI(4);
        dat4(pi,qi).cI.Q                = qRange(qi);
        dat4(pi,qi).cI.conduit_radius   = conduit_radius(qi);
        dat4(pi,qi).cI.N0               = N0_range(qi);
        dat4(pi,qi).cI.dP               = dP_vals(pi,qi);
        dat4(pi,qi).cI.phi0             = h2oExcessMass2ReservoirPorosity(excess_gas_range(qi),dat4(pi,qi).cI);
        dat4(pi,qi).cI                  = getConduitSource(dat4(pi,qi).cI);
    end
end

% 5) Run some reprentative Ze ranges  --------------------------------------
desc5  = 'Vary P0 from best 3.N0 result';
plotsets5 = [10 13 17 19 21 ; 7 10 13 16 19]; % Center index is "best" result for now
cI(5) = getConduitSource(dat3(12,1).cI); % Best 10^8 result from 3       
conduit_radius      = [39 85];
% cI(3).dP             = 1e7; %3.38e7;   
N0_range            = [1e14 1e15];  % Largely MER dependent in initial tests, so tie to Q for now
excess_gas_range    = [0.1663 0.1451]; % EXCESS exsolved volatiles, best from (3)
dP                  = [1e7 2e7];  % Chamber overpressure range
Zw_range            = [0 800];
Zw_vals             = linspace(Zw_range(1),Zw_range(2),na);
for pi = na:-1:1
    for qi = nq:-1:1
        dat5(pi,qi).cI                  = cI(5);
        dat5(pi,qi).cI.Q                = qRange(qi);
        dat5(pi,qi).cI.conduit_radius   = conduit_radius(qi);
        dat5(pi,qi).cI.N0               = N0_range(qi);
        dat5(pi,qi).cI.dP               = dP(qi);
        dat5(pi,qi).cI.phi0             = h2oExcessMass2ReservoirPorosity(excess_gas_range(qi),dat5(pi,qi).cI);
        dat5(pi,qi).cI.Zw               = Zw_vals(pi);
        dat5(pi,qi).cI                  = getConduitSource(dat5(pi,qi).cI);
    end
end

% 6) Test R searches  --------------------------------------
desc6  = 'Test radius variation range using runSet 4 example';
plotsets6 = [3 7 9 11 12 15 17]; % Center index is "best" result for now
% cI(6) = getConduitSource(dat4(15,2).cI); % OPTION 1: Best 10^9 result from 4       
cI(6) = getConduitSource(dat4(15,2).cI); % OPTION 2: Failing 10^9 result from 4     
conduit_radius      = 85;
radius_range        = [0.8 1.1];
radii = conduit_radius.*linspace(radius_range(1),radius_range(2),na);

for pi = na:-1:1
        dat6(pi).cI                  = cI(6);
        dat6(pi).cI.conduit_radius   = radii(pi);
        dat6(pi).cI                  = getConduitSource(dat6(pi).cI);
end
dat6 = dat6';



% ========= Initial sweep params ============
% Starting from reasonable core scenario of conIn2
%   - sweep dP and (covary?) phi0,N0?
%   - eval fragmenting and non-fragmenting valid results





%% Single runs

% D(1).cI = getConduitSource(conIn);
% D(2).cI = getConduitSource(conIn2);
% D(3).cI = getConduitSource(conIn3);
% D(4).cI = getConduitSource(conIn4);

% Run a single test pair for composition - works?
% D(1).cO = Conduit_flow_with_nucleation_V6(D(1).cI); % Run v6 control
% D(2).cO = Conduit_flow_with_nucleation_V7(D(2).cI); % Run v7 with standard composition
% D(3).cO = Conduit_flow_with_nucleation_V7(D(3).cI); % Run v7 with standard composition
% D(4).cO = Conduit_flow_with_nucleation_V7(D(4).cI); % Run v7 with standard composition

% plotConduitOutput(D)
% [ax,th] = plotConduitObjective(D);
%% Run some simple sweep tests, output initial plots of result + validity

if run_set(1)  % Base set
   [dat1,ax1] = getQuickConduitSweep(dat1,false,'conduit_radius',desc1);
%     for di = size(dat1,1):-1:1
%         for dj = size(dat1,2):-1:1
%             dat1(ri,qi).cO = Conduit_flow_with_nucleation_V7(dat1(ri).cI);
%         end
%     end
    plotConduitOutput(dat1)
    set(gcf,'Name',desc1)
end

if run_set(2)  % Excess P to reach surface
   [dat2,ax2] = getQuickConduitSweep(dat2,plotObj,'conduit_radius',desc2);
   
   plotConduitOutput(dat2(plotsets2(1,:),1))
   set(gcf,'Name',[desc2 sprintf(', Q = %.0e kg/s',dat2(1,1).cI.Q) ])
   plotConduitOutput(dat2(plotsets2(2,:),2))
   set(gcf,'Name',[desc2 sprintf(', Q = %.0e kg/s',dat2(1,2).cI.Q) ])

   % Test radius search - unbounded
   [Rlims,cIo,cOo,success] = conduitRadiusFromQ(cI(2),[],'verbose',true);
%    % Test radius search - bracketed but no known valid cases
   [Rlims,cIo,cOo,success] = conduitRadiusFromQ(cI(2),aRange.*aVals(1),'verbose',true);
%    % Test radius search, known valid cases
%    [Rlims2,cIo2,cOo2,success2] = conduitRadiusFromQ(cI(2),aRange.*aVals(1),...
%        'Rvalid',aVals(1).*aScale(8),'verbose',true);

end

if run_set(3)  % Excess gas up to ~50-65%, using 2.P0
   [dat3,ax3] = getQuickConduitSweep(dat3,plotObj,'phi0',desc3);
   
   plotConduitOutput(dat3(plotsets3(1,:),1))
   set(gcf,'Name',[desc3 sprintf(', Q = %.0e kg/s',dat3(1,1).cI.Q) ])
   plotConduitOutput(dat3(plotsets3(2,:),2))
   set(gcf,'Name',[desc3 sprintf(', Q = %.0e kg/s',dat3(1,2).cI.Q) ])

end

if run_set(4)  % Vary P0 from best 3.N0 result
   [dat4,ax4] = getQuickConduitSweep(dat4,plotObj,'dP',desc4);
   
   plotConduitOutput(dat4(plotsets4(1,:),1))
   set(gcf,'Name',[desc4 sprintf(', Q = %.0e kg/s',dat4(1,1).cI.Q) ])
   plotConduitOutput(dat4(plotsets4(2,:),2))
   set(gcf,'Name',[desc4 sprintf(', Q = %.0e kg/s',dat4(1,2).cI.Q) ])

end

if run_set(5)  % Vary P0 from best 3.N0 result
   [dat5,ax5] = getQuickConduitSweep(dat5,plotObj,'dP',desc5);
   
   plotConduitOutput(dat5(plotsets5(1,:),1))
   set(gcf,'Name',[desc5 sprintf(', Q = %.0e kg/s',dat5(1,1).cI.Q) ])
   plotConduitOutput(dat5(plotsets5(2,:),2))
   set(gcf,'Name',[desc5 sprintf(', Q = %.0e kg/s',dat5(1,2).cI.Q) ])

end

if run_set(6)  % Testing radii
   [dat6,ax6] = getQuickConduitSweep(dat6,plotObj,'conduit_radius',desc6);
   
   plotConduitOutput(dat6(plotsets6(1,:),1))
   set(gcf,'Name',[desc6 sprintf(', Q = %.0e kg/s',dat6(1,1).cI.Q) ])
   
%    % Test radius search - no known valid cases
%    [Rlims,cIo,cOo,success] = conduitRadiusFromQ(cI(6),radius_range.*cI(6).conduit_radius,'verbose',true);
%    % Test radius search, known valid cases
%    [Rlims2,cIo2,cOo2,success2] = conduitRadiusFromQ(cI(6),radius_range.*cI(6).conduit_radius,...
%        'Rvalid',cI(6).conduit_radius,'verbose',true);
end


%% Get me some estimates for reasonably universal radii search bounds/steps
if all(run_set([2 3 4 6]))
    di = [dat2(:); dat3(:); dat4(:); dat6(:)];
    nn = numel(di);
%     nn = sum([numel(dat2) numel(dat3) numel(dat4) numel(dat6)]);
    Rtest.conduit_radius     = nan(nn,1);
    Rtest.dP                 = nan(nn,1);
    Rtest.phi0               = nan(nn,1);
    Rtest.Q                  = nan(nn,1);

    ff = fieldnames(Rtest);
    
    for ii=1:nn
        if di(ii).cO.Outcome.Valid
            for fi = 1:length(ff)
                Rtest.(ff{fi})(ii) = di(ii).cI.(ff{fi});
            end
        end
    end

    % Other reference sweep data
    dataDir = '~/Kahuna/data/glaciovolc/coupledSweeps/manuscript_v1_submitted_data';
    sweepFile = 'coupledSweep_2021-09-02_MWIv2_hiLat_n1785_35q_0zw.mat';
    summFile = ['outputSummary_' sweepFile];
    load(fullfile(dataDir,summFile),'qA')
    
    % The original R search data for Rowell et al:
    dataDir = '/Users/crrowell/Kahuna/data/glaciovolc/conduitSweeps/';
    interpFile = fullfile(dataDir,'conduitV6_fineSweep_n64821_21-05-16_compressedV4_fineV2_rInterpV1.mat');
    load(interpFile,'Rmax','Rmin','MER','pf')
    MER = repmat(MER,[1,size(Rmax,2)]);
    pf  = repmat(pf,[size(Rmax,1),1]);
    
    figure
    ax1 = subplot(1,1,1);
    scatter(ax1,log10(MER(:)), Rmax(:)./extrapVentRadius(MER(:)),...
        50,Rmax(:).*0,'Marker','+')
    hold on
    scatter(ax1,log10(MER(:)), Rmin(:)./extrapVentRadius(MER(:)),...
        50,Rmin(:).*0,'Marker','x')    
    scatter(ax1,log10(qA.cI.Q(:)),qA.cI.conduit_radius(:)./extrapVentRadius(qA.cI.Q(:)),...
        50,qA.cI.pf(:).*0+1e7,'Marker','o')
    hold on
    scatter(ax1,log10(Rtest.Q),Rtest.conduit_radius./extrapVentRadius(Rtest.Q),...
                50,Rtest.dP,'filled','Marker','d')

    xlabel('log10(Q)')
    ylabel('Valid R / predicted R')
%     scatter(ax1,Rtest.dP,Rtest.conduit_radius./extrapVentRadius(Rtest.Q),...
%                 50,Rtest.phi0,'filled')
%     hold on
%     scatter(ax1,qA.cI.pf(:).*0,qA.cI.conduit_radius(:)./extrapVentRadius(qA.cI.Q(:)),...
%         50,Rtest.phi0,'filled','Marker','*')
end

% radius = conIn3.conduit_radius .* (0.98:0.001:0.992);

% ax = [];
% for ri = length(radius):-1:1
%     dd(ri).cI = cI;
%     dd(ri).cI.conduit_radius = radius(ri);
%     dd(ri).cI = getConduitSource(dd(ri).cI);
%     dd(ri).cO = Conduit_flow_with_nucleation_V7(dd(ri).cI);
%     
%     [ax] = plotConduitObjective(dd(ri),[],ax);
% end
% [ax] = plotConduitObjective(dd,'conduit_radius');
% plotConduitOutput(dd(10:13))

%% Define quick sweep function
function [dat,ax] = getQuickConduitSweep(dat,plotObj,varname,desc)
    for dj = size(dat,2):-1:1
        for di = size(dat,1):-1:1
            [dat(di,dj).cO,dat(di,dj).cI] = conduitFlowRun(dat(di,dj).cI);
%             try
%                 dat(di,dj).cO = Conduit_flow_with_nucleation_V7(dat(di,dj).cI);
%             catch ME
%                 dat(di,dj).cI
%                 dat(di,dj).cO = ConduitOutcome.getErrorOutcomeFields;
%                 dat(di,dj).cO.Outcome = ConduitOutcome(ME);
% %                 dat(di,dj).cO.msg = ME.message;
% %                 dat(di,dj).cO.Except = ME;
%             end
        end
        
        if plotObj
            [fig,ax,~,~] = plotConduitObjective(dat(:,dj),varname);
            set(fig,'Name',[desc sprintf(', Q = %.0e kg/s',dat(1,dj).cI.Q) ])
        else
            ax = [];
        end
    end
end
