
% Control script to run refined conduit model radius searches

outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/';
refinedDir = fullfile(outDir,'refinedSweep');
testFile = '2024-04-26_myojin_Q8_Z05300_Zw700_357n_dP_21_n0_excess_17.mat' ;

% For loading files
descriptor = 'myojin';

varTol = 0.06;
nCores = 4;

singleRunTest = true;
runAll = false;
%%
fileList = dir(fullfile(outDir, ['*' descriptor '*']));
fileNames = {fileList.name}';

%% Single run test
if singleRunTest
    pool = parpool(nCores);
    [dat,varMin,varMax,outcomeCodes] = refineConduitSweep(fullfile(outDir,testFile),refinedDir,varTol,pool);
    delete(pool)
end

%% Run all
if runAll
    pool = parpool(nCores);
    for fi = 1:length(fileNames)
        fprintf('%i/%i ...',fi,length(fileNames))
        [dat,varMin,varMax,outcomeCodes] = refineConduitSweep(fullfile(outDir,testFile),refinedDir,varTol,pool);
        fprintf('Done\n')
    end
    delete(pool)
end