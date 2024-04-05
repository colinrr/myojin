c%% ======================================================================
%            SETUP FOR MJYOJIN KNOLL SWEEPS/MONTE CARLO
% =======================================================================
% C Rowell, Sep 2023

clear all; close all
% Uses Hajimirza conduit model, V7

% Mac
codeDir = '~/code/research-projects/hydroVolc/';
outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep1/';

% SJ
% codeDir = 'C:\Users\crowell\Documents\GitHub\hydroVolc\';
% outDir = 'C:\Users\crowell\Kahuna\data\myojin\mainSweep1';


addpath(genpath(codeDir))

run_test_sweep      = false;
run_full_sweep      = false;
run_sweep_summary   = true;
plot_sweep_summary  = true;
    simplifyCodes   = true;
    printfigs       = true;
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
sweepVars.dP.range = [0 4e7]; % 0 - 40 MPa overpressures
sweepVars.dP.n     = 21;

% Additional KEY note here: DOMINANT range of quartz-hosted melt inclusion
% H20 is ~5.5-6.5 wt% H20. Model assumes something SOLUBILITY based, so we
% should highlight the set of results that correspond to about 6 wt % total
% for n0 + n0_excess.
sweepVars.n0_excess.range = [0 0.8];
sweepVars.n0_excess.n     = 17;

% Do for:
Qset = [1e8 1e9];
N0set = [1e14 1e15]; % Fix these initial BND's roughly to Q based on early conduit results
                   % '-> (see N_vs_phi_checks.m and bnd_h20_MLRFits.m)
 
% ----------- Search params -------------
sweepParams.verbose    = false;
sweepParams.cores      = 8;
sweepParams.fixed      = 'R';
descriptor = 'myojin';

% Define conduit range params to search for shooting solutions
% aRange = [0.6 1.4]; % Fractional range of conduit radius to search?

% Search Tolerances
% dRminScale  = 0.2e-3;
% dQminScale  = 0.2e-3;
% maxIter     = 10;

% OTHER NOTES
% At 2.2 km conduit length, phi0=31.4 and N0=1.198e14 are good starting
% values that correspond to the default run...huh

% Other:
% vent altitude to be set at -Zw (below sea level)
% conduit depth to be set at (abs chamber depth (bsl) - vent depth (bsl))

% - default atmosphere (only for small pressure contribution of atmosphere)
% - default conduit friction coefficient

% ------ TEST SET -----
testSweep.dP.range = [1e6 1e7]; % 0 - 40 MPa overpressures
testSweep.dP.n     = 4;
testSweep.n0_excess.range = [0.1 0.3];
testSweep.n0_excess.n     = 2;

testPars.verbose    = false;
testPars.cores      = 2;
testPars.fixed      = 'R';
testPars.descriptor = 'myojinSweepTest';
%% MYOJIN COMMON PARAMS ------
% conIn.Z0             = 2200;        % Below vent or sea level? 
% conIn.Zw             = 0;         % (Still processing quench depths, but possibly up to 1km?) can run a range
% conIn.T              = 900+273.15;  % Assumed similar to Sumisu (Shukuno et al)
conIn.phi_frag       = 0.75;
conIn.rho_melt       = 2350;

% conIn.vh0            = -conIn.Zw; % Vent below sea level
% conIn.conduit_radius = 39.66; % based on initial runs;

% conIn.N0             = 1.2e14;           % most juvenile (A13) = 1e9 cm-3 = 1e15 m^-3
% conIn.phi0           = 0;           % unsure, but 0.3 is about the value 
                                        % of the juvenile population with spherical bubbles
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




%% Run TEST set

if run_test_sweep
    cI_test = conIn;
    cI_test.Zw = min_vent_depth;
    cI_test.Z0 = current_caldera_floor_depth-min_vent_depth;

    [dat,Rmin,Rmax,outcomeCodes] = conduitParameterSweep(cI_test,testSweep,outDir,testPars);
end

%% Full run set
if run_full_sweep
    nS = length(Qset).*length(sweepList.Zw);
    fprintf('Total sweep set: %i\n',nS);
    for qi = 1:length(Qset)
        for si = 1:length(sweepList.Zw)
            
            cI = conIn;
            cI.Q  = Qset(qi);
            cI.N0 = N0set(qi);
            
            cI.Zw = sweepList.Zw(si);
            cI.Z0 = sweepList.Z0(si);

            sweepParams.descriptor = sprintf('%s_Q%i_Z0%i_Zw%i',...
                descriptor,log10(cI.Q),cI.Z0,cI.Zw);

            fprintf('Sweep %i/%i: %s\n', qi.*si, nS, sweepParams.descriptor)
            [dat,Rmin,Rmax,outcomeCodes] = conduitParameterSweep(cI,sweepVars,outDir,sweepParams);
        end
    end
end

%% RUN SWEEP SUMMARY DATA GATHERING AND PLOTS

knownFailsToMap = {'Check Pg, T'}; % update code assignment below if needed

if run_sweep_summary
    % Get list of sweep files
    fileList = dir(fullfile(outDir, ['*' descriptor '*']));
    fileNames = {fileList.name}';
    
    % 32 sweeps in N0, dP
    %  - [2 chamber depths x 7 water depths plus 1 control depths x 2 water
    %  depths] x 2 MER
    %  - so control runs gives n = 4 [2x2] plot
    %  - plus 4 x 7 sweeps over water depth
    dP = linspace(sweepVars.dP.range(1),sweepVars.dP.range(2),sweepVars.dP.n);
    n0_excess = linspace(sweepVars.n0_excess.range(1),sweepVars.n0_excess.range(2),sweepVars.n0_excess.n);
    % -> transform this to porosity?
    
    % ----------- FIGURE SETUP ----------
    if plot_sweep_summary
        

        % Control sweep figures
        nr = 2;
        nc = 2;
        dx = 0.01;
        dy = 0.09;
        ppads = [0.08 0.2 0.11 0.05];
        cbposCTL = [0.82 0.025];
        fs = 6;
        fscb = 6;
        
        ctrlFig = figure('position',[10 200 1100 800],'Name','Control runs');
        for ii=1:4 % 2 Zw x 2 Q
            ctlax(ii) = tightSubplot(2,2,ii, dx, dy, ppads);

            if ~or(ii==1,ii==8)
                set(ctlax(ii),'YTickLabel',[])
            end
            if ii<8
                set(ctlax(ii),'XTickLabel',[])
            end
        end

        % Main sweep figures
        nr = 2;
        nc = 7;
        dx = 0.01;
        dy = 0.11;
        ppads = [0.05 0.15 0.07 0.05];
        cbpos = [0.86 0.025];

        for qi = 1:2 % 2 Q's

            mjFig(qi) = figure('position',[10 200 1500 600],'Name',sprintf('Q_0 = 10^%.1f kg/s',log10(Qset(qi)))); 
            for ii=1:14 % 2 Z0 x 7 Zw
                ax(ii,qi) = tightSubplot(nr,nc,ii, dx, dy, ppads);

                if ~or(ii==1,ii==8)
                    set(ax(ii),'YTickLabel',[])
                end
                if ii<8
                    set(ax(ii),'XTickLabel',[])
                end
            end
        end
    end
    % -------------------------------------
    
    for qi = 1:length(Qset)
        for si = 1:length(sweepList.Zw)
            thisDescriptor = sprintf('%s_Q%i_Z0%i_Zw%i',...
                descriptor,log10(Qset(qi)),sweepList.Z0(si),sweepList.Zw(si));
            
            fileIdx = find(contains(fileNames,thisDescriptor));
            
            load(fullfile(outDir,fileNames{fileIdx}))
            nSuccess = cellfun(@(x) sum(x>0),outcomeCodes.all); % Check for more than 1 success type
            if any(nSuccess(:) > 1)
                pause(0.1)  %Hold up
            end
            
            % FIX UP CODES
            % Some temporary stuff to retrieve most important code results
            %   - since maxR is spitting out zero codes instead of negatives
            maxCode = cellfun(@max,outcomeCodes.all);
            % Account for new to-be-mapped errors
            unmappedFails = (maxCode == -22);
            umi = find(unmappedFails);
            for um = 1:sum(unmappedFails(:))
                if ismember(outcomeCodes.failMsg(umi(um)).ME.message,knownFailsToMap)
                    maxCode(umi(um)) = -21; % update as needed if knownFailsToMap is updated
                end
            end
            % Update
            codeZero = or(outcomeCodes.maxR < 1,isnan(outcomeCodes.maxR));
            plotCodes = outcomeCodes.maxR;
            plotCodes(codeZero) = maxCode(codeZero);
            
            % getSweepSummaryHERE when needed
            
            if plot_sweep_summary

                % For now just plot outcome codes
                
                % Plot CONTROL RUNS
                if si < 3
                    axi = (qi-1)*2 + si;
                    axes(ctlax(axi))
                    [cmap,cax,cticks,clabels,outcomeIndex,~] = outcomeColorMap(plotCodes,simplifyCodes, false);
                    imagesc(n0_excess*100,dP/1e6, outcomeIndex )
                    
                    colormap(gca,cmap)
 
                    set(gca,'YDir','normal','FontSize',fs)
                    caxis(cax)
                    
                    if qi==1
                        title(sprintf('Z_e = %.0fm,   Z_0 = %.1fkm',sweepList.Zw(si),sweepList.Z0(si)/1e3))
                    else
                        xlabel('Excess n_0 (%)')
                    end
                    if si==1
                        ylabel({sprintf('\\bfQ_0 = 10^{%.1f} kg/s\\rm',log10(Qset(qi))),'','\Delta P (MPa)'})
                    else
                        set(gca,'YTickLabel',[])
                    end
                    if qi==1 && si==1
                        ctlcb = colorbar(gca,'location','eastoutside');
                        caxis(cax);
                        ctlcb.Ticks = cticks;
                        ctlcb.TickLabels = clabels;

                        axpos1 = get(ctlax(1),'Position');
                        axpos2 = get(ctlax(3),'Position');
                        ctlcb.Position = [cbposCTL(1) axpos2(2) cbposCTL(2) sum(axpos1([2 4]))-axpos2(2)];
                        ctlcb.FontSize = fscb;
                    end
                 
                % Plot MYOJIN RUNS
                else
                    axi = si-2;
                    axes(ax(axi,qi))
                    [cmap,cax,cticks,clabels,outcomeIndex,~] = outcomeColorMap(plotCodes,simplifyCodes, false);
                    imagesc(n0_excess*100,dP/1e6, outcomeIndex )

                    colormap(gca,cmap)
 
                    set(gca,'YDir','normal','FontSize',fs)
                    caxis(cax)
                    
                    title(sprintf('Z_e = %.0fm,   Z_0 = %.1fkm',sweepList.Zw(si),sweepList.Z0(si)/1e3))
                    if si>=10
                        xlabel('Excess n_0 (%)')
                    end
                    if si==3
                        ylabel({sprintf('\\bfChamber Depth = %.1f\\rm',chamber_depths(1)),'','\Delta P (MPa)'})
                    elseif si ==10
                        ylabel({sprintf('\\bfChamber Depth = %.1f\\rm',chamber_depths(2)),'','\Delta P (MPa)'})
                    else
                        set(gca,'YTickLabel',[])
                    end
                    if si==3
                        ctlcb = colorbar(gca,'location','eastoutside');
                        caxis(cax);
                        ctlcb.Ticks = cticks;
                        ctlcb.TickLabels = clabels;

                        axpos1 = get(ax(1),'Position');
                        axpos2 = get(ax(end),'Position');
                        ctlcb.Position = [cbpos(1) axpos2(2) cbpos(2) sum(axpos1([2 4]))-axpos2(2)];
                        ctlcb.FontSize = fscb;
                    end

                end
            end
        end
    end
    

    
end

if plot_sweep_summary && printfigs
    figure(ctrlFig)
    printpdf('Control_runs_OutcomeCodes',outDir,[16.5 12])
    figure(mjFig(1))
    printpdf('Q_1e8_OutcomeCodes',outDir,[25 12])
    figure(mjFig(2))
    printpdf('Q_1e9_OutcomeCodes',outDir,[25 12])
    
end
%% PLOT SUMMARY


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
