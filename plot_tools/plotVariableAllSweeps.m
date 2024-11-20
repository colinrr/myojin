function plotVariableAllSweeps(varName)
% Easy wrapper to plot a variable field across the various parameter sweeps
% 
if nargin==0
    dispVars = true;
else
    dispVars = false;
end

DeltaP_max = 2e7; % Maximum chamber overpressure to keep

try
    run ../config
catch
    error('Could not find config.m. Check working directory!')
end


%% sweep data directory and file names
outputTableName = fullfile(DATA_DIR,'ConduitSweepsDataTableB.mat');
load(outputTableName)

if dispVars
    disp('Argument should be one of the following variable string names:')
    disp(string(T.Properties.VariableNames'))
    return
else
    assert(ismember(varName,string(T.Properties.VariableNames')),'Variable not found in table!')
end

% outDir = fullfile(MAIN_DATA_DIR,'mainSweep2/refinedSweep/');
% sweepFiles = dir(outDir);
% sweepFiles = {sweepFiles.name}';
% sweepFiles = sweepFiles(cellfun(@(x) contains(x,'myojin'),sweepFiles));

%% Select the subset and x/y vars
n0 = unique(T.n0_excess);
dP = unique(T.dP);
Q  = unique(T.Q);
Zw = unique(T.Zw);
Ztot = unique(T.Z_total);

Pscale = 1/1e6;

% Main sweep figure params
nr = 3;
nc = 7;
dx = 0.015;
dy = 0.11;
ppads = [0.05 0.15 0.07 0.05];
msz = 90;
cbX = 0.86;
cbdX = 0.025;
cbY = ppads(3);
cbdY = 1-sum(ppads(3:4));




% Outcome alpha mapping
sCodes = unique(T.simpleCode);
codeAlpha = (sCodes - min(sCodes)) ./ (max(sCodes) - min(sCodes));
codeMarkers = {'o','s','h'};
cax = [];

for qi = 1:length(Q) % 1 figure per Q
    figure('position',[10 200 1500 600],'Name',sprintf('Q_0 = 10^%.1f kg/s',log10(Q(qi))))
    
    count = 0;
    for zti = 1:length(Ztot) % 1 row per chamber depth bsl
        for zwi = 1:length(Zw) % 1 column per zw
            
            subsetIdx = (T.Q==Q(qi) & T.Zw==Zw(zwi) & T.Z_total ==Ztot(zti));
            Tsub = T(subsetIdx,:);
            cax = [cax; min(Tsub.(varName)) max(Tsub.(varName))];
            
            count = count+1;
            ax(zti,zwi) = tightSubplot(length(Ztot),length(Zw),count, dx, dy, ppads);
            
            for ci = 1:length(sCodes)
                codeSub = Tsub.simpleCode==sCodes(ci);
                scatter(Tsub.n0_excess(codeSub),Tsub.dP(codeSub).*Pscale,msz,Tsub.(varName)(codeSub),...
                    'filled',codeMarkers{ci},...
                    'MarkerEdgeColor','k',...
                    'MarkerEdgeAlpha',codeAlpha(ci))
                hold on
            end
            if zwi==1
                ylabel('{\Delta}P (MPa)')
            end
            if zti==length(Ztot)
                xlabel('n_0 excess (wt.%)')
            end
            title(sprintf('Z_{tot} = %.1f km, Z_w = %.0f m',Ztot(zti)/1e3,Zw(zwi)))
%             caxis([0 3])
        end
    end
    cb = colorbar(gca,'location','eastoutside');
    cax = [min(cax(:,1)) max(cax(:,2))];
%     caxis(cax);
    cb.Label.String = strrep(varName,'_','\_');
    cb.Position = [cbX cbY cbdX cbdY];
    for ai = 1:numel(ax)
        caxis(ax(ai),[cax])
    end
end


%%
end