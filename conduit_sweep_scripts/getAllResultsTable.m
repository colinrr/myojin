% Build a tabular dataset of sweep output
%  GET:
%   -> MER
%   -> n0
%   -> Ztotal, Zw, Z0 (conduit length)
%   -> # nucleation events and their (onset) depths
%   -> Ascent time
%   -> Avg. decompression rate over ascent time
%   -> Outcome Code - full and simple
%   -> BND at fragmentation or surface - whichever is first
%   -> Pressure, depth, and dissolved H20 at saturation
%       -> For dissolved H20, also get equilibrium value


outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';

outcomeFile = 'outcomeCodeSummary.mat';


%% Do the thing

load(fullfile(outDir,outcomeFile))
% Sweep names
fn = fieldnames(allOutcomeCodes);
fn = fn(cellfun(@(x) contains(x,'myojin'),fn));

% File names
ff = dir(outDir);
ff = {ff.name}';
ff = ff(cellfun(@(x) contains(x,'myojin'),ff));

% loop through all sweep files and pull parameters
for fi = 1 %:length(fn)
    sweep = fn{fi};
    whichFile = contains(sweep, ff);
    sweepFile  = ff{whichFile};
    
    load(fullfile(outDir,sweepFile))
end