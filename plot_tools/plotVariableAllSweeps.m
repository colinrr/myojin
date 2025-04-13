function varargout = plotVariableAllSweeps(varName, customVar, logCax, forPrint)
% Easy wrapper to plot a variable field across the various parameter sweeps
% Can also add 'customVar' to plot your own, just must match the length of
% the table.
%  - call with no args to show the variable options for plotting
%
% Optionally output the plot axes, legend and colorbar handles
if nargin==0
    dispVars = true;
else
    dispVars = false;
end
if nargin<2 || isempty(customVar)
    customVar = [];
    useCustomVar = false;
else
    useCustomVar = true;
end
if nargin<3 || isempty(logCax)
    logCax = false;
end
if nargin<4
    forPrint = false;
end

DeltaP_max = 2e7; % Maximum chamber overpressure to keep

try
    run config
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
xl = [-0.05 0.85];
yl = [-.5 20.5];

% Main sweep figure params
nr = 3;
nc = 7;
dx = 0.015;
dy = 0.11;

if forPrint
    msz = 20;
    ppads = [0.05 0.12 0.07 0.1];
    cbX = 0.92;
    cbdX = 0.015;
    cbY = ppads(3);
    cbdY = 1-sum(ppads(3:4));
else
    msz = 90;
    ppads = [0.05 0.15 0.07 0.1];
    cbX = 0.86;
    cbdX = 0.025;
    cbY = ppads(3);
    cbdY = 1-sum(ppads(3:4));
end

% Outcome alpha mapping
sCodes = unique(T.simpleCode);
codeAlpha = (sCodes - min(sCodes)) ./ (max(sCodes) - min(sCodes));
cax = [];

axout = {};
for qi = 1:length(Q) % 1 figure per Q
    figs{qi} = figure('position',[10 200 1500 600],'Name',sprintf('Q_0 = 10^%.1f kg/s, %s',log10(Q(qi)),varName));
    
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
                [codeMarker, ~, cAlpha] = getCodeMarker(sCodes(ci));
                
                scatter(Tsub.n0_excess(codeSub),Tsub.dP(codeSub).*Pscale,msz,cVar(codeSub),...
                    'filled',codeMarker,...
                    'MarkerEdgeColor','k',...
                    'MarkerEdgeAlpha',cAlpha)
                hold on
%                 lh(ci) = scatter(nan,nan,msz,[0.5 0.5 0.5],...
%                     'filled',codeMarker,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerEdgeAlpha',codeAlpha(ci),...
%                     'DisplayName','Intrusive'); 

            end
            
            xlim(xl)
            ylim(yl)
            if forPrint
                if zwi==1
                    ylabel('${\Delta}P$ (MPa)', 'Interpreter', 'latex')
                end
                if zti==length(Ztot)
                    xlabel('$n_{0e}$ (wt.\%)', 'Interpreter', 'latex')
                end
                title(sprintf('$Z_{tot}$ = %.1f km\n $Z_w$ = %.0f m', Ztot(zti)/1e3,Zw(zwi)), 'Interpreter','latex')
                if logCax
                    set(gca,'colorscale','log')
                end
            else
                if zwi==1
                    ylabel('{\Delta}P (MPa)')
                end
                if zti==length(Ztot)
                    xlabel('Excess H_2O n_{0e} (wt.%)')
                end
                title(sprintf('Z_{tot} = %.1f km, Z_w = %.0f m',Ztot(zti)/1e3,Zw(zwi)))
                if logCax
                    set(gca,'colorscale','log')
                end
            end
        end
    end
    cb{qi} = colorbar(gca,'location','eastoutside');
    cax = [min(cax(:,1)) max(cax(:,2))];
%     caxis(cax);
    cb{qi}.Label.String = strrep(varName,'_','\_');
    cb{qi}.Position = [cbX cbY cbdX cbdY];
    for ai = 1:numel(ax)
        caxis(ax(ai),[cax])
    end

    for ci = [1, 2, 3]  % 1:length(sCodes)
        [cm, cd, ca] = getCodeMarker(ci-1);
        lh(ci) = scatter(nan,nan,msz,[0.5 0.5 0.5],...
        'filled', cm,...
        'MarkerEdgeColor','k',...
        'MarkerEdgeAlpha',ca,...
        'DisplayName',cd);     
    end
    try
        la{qi} = legend(lh,'Location','southoutside','Orientation','horizontal');
    catch ME
        pause(0.5)
        warning('Error creating legend')
        
    end
    set(la{qi},'Position',[0.5-la{qi}.Position(3)/2 1-ppads(4)+0.07  la{qi}.Position(3) la{qi}.Position(4)])
    if nargout >= 1
        axout{qi} = ax; 
    end

end


%% output axes handles?
if nargout>=1
    varargout{1} = axout;
end
if nargout >= 2
    varargout{2} = la;
end
if nargout >= 3
    varargout{3} = cb;
end
if nargout >= 4
    varargout{4} = figs;
end
end
% 
% function [codeMarker, codeDesc] = getCodeMarker(codeSub)
% % codeMarkers = {'o','s','^'};
% % sCodeDescriptions = ["Intrusive","Effusive","Fragmenting"];
% 
%     if codeSub==3 || codeSub==2 % P-bal and explosive frag
%         codeDesc = "Fragmenting";
%         codeMarker = '^';
%     elseif codeSub==1 || codeSub==-1 % Taking all effusive as valid (see Obsidian notes)
%         codeDesc = "Effusive";
%         codeMarker = 's';
%     elseif codeSub==0
%         codeDesc = "Intrusive";
%         codeMarker = 'o';
%     else
%         codeDesc = "Invalid";
%         codeMarker = '.';
%     end
% end