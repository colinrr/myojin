dataDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
codesFile = 'outcomeCodeSummary.mat'

% Test file
% foi = '2024-04-26_myojin_Q8_Z0900_Zw500_357n_dP_21_n0_excess_17.mat';
foi = 'myojin_Q8_Z0900_Zw500';


% Test run on multiples
% files = {'2024-04-26_myojin_Q8_Z0900_Zw500_357n_dP_21_n0_excess_17.mat',
%         } ;

    
    
%% load 
load(fullfile(dataDir,codesFile))

% 1 Cut to a specified max pressure
% 2 replace invalid effusive with effusive
% 3 Nearest neighbour interp to replace errored results?
% 4 Conncomp to get remaining chunks
%   -> get outlines
%   -> Any gaps to fill?
%   -> smooth outlines
%   -> exclude regions
