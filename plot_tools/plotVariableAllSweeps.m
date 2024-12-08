function varargout = plotVariableAllSweeps(varName, customVar)
% Easy wrapper to plot a variable field across the various parameter sweeps
% Can also add 'customVar' to plot your own, just must match the length of
% the table.
% 
% Optionally output the plot axes
%
if nargin==0
    dispVars = true;
else
    dispVars = false;
end
if nargin<2
    customVar = [];
    useCustomVar = false;
else
    useCustomVar = true;
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

if useCustomVar
    assert(and(isvector(customVar),length(customVar)==size(T,1)),"'customVar' must be a vector with length equal to the number of rows in the table.")
else
    if dispVars
        disp('Argument should be one of the following variable string names:')
        disp(string(T.Properties.VariableNames'))
        return
    else
        assert(ismember(varName,string(T.Properties.VariableNames')),'Variable not found in table!')
    end
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
ppads = [0.05 0.15 0.07 0.09];
msz = 90;
cbX = 0.86;
cbdX = 0.025;
cbY = ppads(3);
cbdY = 1-sum(ppads(3:4));
lpos = [0.5 1-ppads(4)];



% Outcome alpha mapping
sCodes = unique(T.simpleCode);
sCodeDescriptions = ["Intrusive","Effusive","Fragmenting"];
codeAlpha = (sCodes - min(sCodes)) ./ (max(sCodes) - min(sCodes));
codeMarkers = {'o','s','^'};
cax = [];

axout = {};
for qi = 1:length(Q) % 1 figure per Q
    figure('position',[10 200 1500 600],'Name',sprintf('Q_0 = 10^%.1f kg/s, %s',log10(Q(qi)),varName))
    
    count = 0;
    for zti = 1:length(Ztot) % 1 row per chamber depth bsl
        for zwi = 1:length(Zw) % 1 column per zw
            
            subsetIdx = (T.Q==Q(qi) & T.Zw==Zw(zwi) & T.Z_total ==Ztot(zti));
            Tsub = T(subsetIdx,:);
            if useCustomVar
                cVar = customVar(subsetIdx);
            else
                cVar = Tsub.(varName);
            end
            
            cax = [cax; min(cVar) max(cVar)];
            
            count = count+1;
            ax(zti,zwi) = tightSubplot(length(Ztot),length(Zw),count, dx, dy, ppads);
            
            for ci = 1:length(sCodes)
                codeSub = Tsub.simpleCode==sCodes(ci);
                scatter(Tsub.n0_excess(codeSub),Tsub.dP(codeSub).*Pscale,msz,cVar(codeSub),...
                    'filled',codeMarkers{ci},...
                    'MarkerEdgeColor','k',...
                    'MarkerEdgeAlpha',codeAlpha(ci))
                hold on
                lh(ci) = scatter(nan,nan,msz,[0.5 0.5 0.5],...
                    'filled',codeMarkers{ci},...
                    'MarkerEdgeColor','k',...
                    'MarkerEdgeAlpha',codeAlpha(ci),...
                    'DisplayName','Intrusive');                
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
    for ci = 1:length(sCodes)
        lh(ci) = scatter(nan,nan,msz,[0.5 0.5 0.5],...
        'filled',codeMarkers{ci},...
        'MarkerEdgeColor','k',...
        'MarkerEdgeAlpha',codeAlpha(ci),...
        'DisplayName',sCodeDescriptions(ci));                
    end
    la = legend(lh,'Location','southoutside','Orientation','horizontal');
    set(la,'Position',[0.5-la.Position(3)/2 1-ppads(4)+0.06  la.Position(3) la.Position(4)])
    if nargout == 1
        axout{qi} = ax; 
    end
end


%% output axes handles?
if nargout==1
    varargout{1} = axout;
end
end