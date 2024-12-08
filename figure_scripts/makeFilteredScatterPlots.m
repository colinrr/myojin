%% Make filtered scatterplots

clear all; %close all
config

% How to make this plot...
% -> For each conduit length, pick a single fixed excess pressure
%      (gas depletion would make this go down, vent deepening would make it
%      go up, so we'll ignore variation for now)
% -> Find the n0_excess which most closely corresponds to a total gas
% content.
% --> From the resulting filtered subsets, plot the key variables, e.g.:
% --- Psat vs n0 ---- 
% --- Psat vs LAST Z_nucl ----
% --- Psat vs dPdt ----
% --- Psat vs phi ---- ???


% Set up data filter scenarios
filterScenarios.Z_total = [1400 2200 6000];
filterScenarios.dP = [20000000 16000000 10000000];

runQ = [1e9]; % 

n0_target = 0.06;
%% 

load(fullfile(DATA_DIR,'ConduitSweepsDataTableB.mat'),'T')
allZw = unique(T.Zw);

for qi = 1:length(runQ)
    Tsub = T((T.Q==runQ(qi)),:);
    
    for si = 1:length(filterScenarios.Z_total)
        
        
        for zwi = 1:length(allZw)
            subsetIdx = (Tsub.Z_total == filterScenarios.Z_total(si) ...
                & Tsub.dP == filterScenarios.dP(si) ...
                & Tsub.Zw == allZw(zwi));

            n0t_subset = Tsub(subsetIdx,:).n0_total;
            [n0_vals(si,zwi),~] = closest(n0_target,n0t_subset);
        end
        

    end
    % General filter
        subsetIdx = (Tsub.Z_total == filterScenarios.Z_total(si) & Tsub.dP == filterScenarios.dP(si));
    
end