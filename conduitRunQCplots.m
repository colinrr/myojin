% QC plots for some of the "invalid effusive" cases

outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/';
foi = '2024-04-26_myojin_Q8_Z05300_Zw700_357n_dP_21_n0_excess_17.mat' ;

    simplifyCodes   = true;

% index pairs in 'dat' of Runs to plot
ip = [8 9; 8 10; 8 11;];
%%

load(fullfile(outDir,foi))

% Get the run sweep vectors
dP = linspace(sweepParams.dP.range(1),sweepParams.dP.range(2),sweepParams.dP.n);
n0_excess = linspace(sweepParams.n0_excess.range(1),sweepParams.n0_excess.range(2),sweepParams.n0_excess.n);

plotCodes = getCodeSummaryArray(outcomeCodes); % Get corrected code array

% Get plottable codes
[cmap,cax,cticks,clabels,outcomeIndex,~] = outcomeColorMap(plotCodes,simplifyCodes, false);

% Make the plot
figure('name',foi)
imagesc(n0_excess*100,dP/1e6, outcomeIndex )
colormap(gca,cmap)
set(gca,'YDir','normal','FontSize',fs)
caxis(cax)

%% Test runs

% This run with no specified R bounds produces no good result - as in the
% first parameter sweep
[Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(1,1),ip(1,2)).cI,[],'verbose',true,'output',true);

% THIS run with Rbounds based on adjacent successful runs finds a solution
% just fine. So (a) the radius search is flawed when unbouned (which I
% knew) and (b) refining the results for better successes may be a
% matter of doing a simple R search and interpolation