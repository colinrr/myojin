% QC plots for some of the "invalid effusive" cases

% 1 -- These were to span an explosive-p balanced-effusive transition with
% missing solutions in the first, unrefined search ---
% outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/';
% foi = '2024-04-26_myojin_Q8_Z05300_Zw700_357n_dP_21_n0_excess_17.mat' ;
% --> index pairs in 'dat' of Runs to plot
% ip = [8 9; 8 10; 8 11;];
% ----------

% 2 --- These are to QC a valid effusive - invalid effusive transition in the
% refined sweep
% outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
% foi = '2024-08-26_myojin_Q9_Z0300_Zw1100_357n_dP_21_n0_excess_17.mat' ;
% ip = [6 5; 7 5; 12 5;];
% -------------

% 3 --- These are to QC a valid effusive - invalid effusive -explosive transition in the
% refined sweep
% outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
% foi = '2024-08-26_myojin_Q9_Z02200_Zw0_357n_dP_21_n0_excess_17.mat' ;
% ip = [5 2; 5 4; 5 5;];
% -------------

% 4 --- These are to QC invalid explosive - valid explosive transition in the
% refined sweep
outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
% outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/pass2/';
% foi = '2024-04-26_myojin_Q8_Z05100_Zw900_357n_dP_21_n0_excess_17.mat' ;
foi = '2024-04-26_myojin_Q8_Z04900_Zw1100_357n_dP_21_n0_excess_17.mat' ;
% foi = '2024-04-26_myojin_Q8_Z01500_Zw700_357n_dP_21_n0_excess_17.mat' ;

ip = [1 16; 3 16; 4 16;];
% -------------

simplifyCodes   = true;


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
hold on
scatter(n0_excess(ip(:,2))*100, dP(ip(:,1))/1e6)

for ii=1:size(ip,1)
    plotDat(ii) = dat(ip(ii,1),ip(ii,2));
end

plotConduitOutput(plotDat)

%% Test runs for initial refined sweep

% This run with no specified R bounds produces no good result - as in the
% first parameter sweep
% [Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(1,1),ip(1,2)).cI,[],'verbose',true,'output',true);

% THIS run with Rbounds based on adjacent successful runs finds a solution
% just fine. So (a) the radius search is flawed when unbouned (which I
% knew) and (b) refining the results for better successes may be a
% matter of doing a simple R search and interpolation

%% QC checks for refined sweep output

% 2 --> Valid/Invalud effusive transition
% %   -> show radius range for the relevent transect
%   [varMin(1:12,5)  varMax(1:12,5) plotCodes(1:12,5)]
% 
%     newDat = plotDat;
%     [Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(3,1),ip(3,2)).cI,[50 70],'verbose',true,'output',true);
%     %    --> This easily produced a wide range of valid effusive solutions with
%     %        R ~ 51 - 59
%     newDat(2).cI = cI; newDat(2).cO = cO;
%     [Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(3,1),ip(3,2)).cI,[40 70],'verbose',true,'output',true);
%     newDat(3).cI = cI; newDat(3).cO = cO;
%     plotConduitOutput(newDat)

% 3 --> Valid/Invalud effusive/Explosive transition
%   -> show radius range for the relevent transect
%   [varMin(5,1:7)'  varMax(5,1:7)' plotCodes(5,1:7)']
% 
%     newDat = plotDat;
%     [Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(2,1),ip(2,2)).cI,[85 100],'verbose',true,'output',true);
%     %    --> For this case, it's still impossible to find a valid
%     solution, but the transition is narrow and the fragmenting solution
%     occurs with frag depth JUST at the surface, so arguably impossible to
%     fragmenting solutions to occur in the invalid zone

% 3 --> Invalid explosive/Explosive transition
%   -> show radius range for the relevent transect
  [varMin(1:9,16)  varMax(1:9,16) plotCodes(1:9,16)]
    [Rlims,cI,cO,validCodes,allCodes] = conduitRadiusFromQ(dat(ip(2,1),ip(2,2)).cI,[27 30],'verbose',true,'output',true);
    % --> ALSO true for invalid explosive results (codes -3 to -6) -> 2nd
    % quick refinement is warranted where we can't obviously infill
    
    % SO FOR NOW, FROM BOTH TESTS ABOVE, for the purpose of plotting outcome codes and regime
    % transitions, I will treat invalid effusive as valid effusive
        % --> ie outcomeCode = -1 --> +1

 