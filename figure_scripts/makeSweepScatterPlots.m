%% Build some plots from Myojin conduit sweeps data table

clear all; close all
run ../config % Assumes run from local path

outputTableFile = fullfile(DATA_DIR,'ConduitSweepsDataTable.mat');

printFigs = true;
%% Load the table
disp('Loading data table...')
load(outputTableFile)

disp('Making plots...')
%% Pair plots
fname = 'Sweep variable pair plots';

% Continous variables
scatter_vars = { 
    'Zf','Frag. Depth','Z_{fr} (m)'
    'Cm_0','Ini. Dissolved H_2O', 'H_2O wt.%'
    'Psat','Supersat. Pressure','P_{sat} (Pa)'
    'dP_dt_bar','Decomp. Rate','dP/dt (Pa/s)'
    'BND_f','Bubble # Dens.','BND_f'
    'phi_0','Ini. Porosity','\phi_0'
    'phi_f','Fin. Porosity','\phi_f'
    't_ascent','Ascent Time','t (s)'
    };

% cVar = T.n0_excess;
cVars = {'simpleCode',[0 2]
         'n0_excess',[0 0.8]
         'dP',[0 2e7]
         'Zw',[0 1300]
         'Z0',[100 6000]
         'N_nucl',[0 2]};

dx = 0.01;
dy = 0.01;
ppads = [0.05 0.06 0.06 0.05];
msz = 20;
cbW = 0.015;
cbX = 0.95;
cbH = 1-sum(ppads(3:4));

nHist = 50;
nS = length(scatter_vars);

for ci = 1:length(cVars)
    cVar = T.(cVars{ci,1});
    figure('position',[50 50 1500 1200],'name',sprintf('%s: %s',fname,cVars{ci,1}))
    count = 0;
    for ii = 1:nS
        yVar = T.(scatter_vars{ii,1});
        yLab = scatter_vars{ii,3};     

        for jj = 1:nS
            count = count+1;
            tightSubplot(nS,nS,count,dx,dy,ppads)
            xVar = T.(scatter_vars{jj,1});
            tLab = scatter_vars{jj,2};
            xLab = scatter_vars{jj,3};        
            if ii==jj
                yLab = 'Counts';
                histogram(xVar,nHist)
            else
                scatter(xVar,yVar,msz,cVar,'.')
            end
            if ii==1
                title(tLab)
            end
            if ii==nS
                xlabel(xLab)
                axp = get(gca,'Position');
                cbY = axp(2);
            else
                set(gca,'XTickLabel',[])
            end
            if jj==1
                ylabel(yLab)
            else
                set(gca,'YTickLabel',[])
            end
            if strcmp(scatter_vars{jj,1},'BND_f')
                set(gca,'XScale','log')
                xlim([1e13 1e18])
            end
            if strcmp(scatter_vars{ii,1},'BND_f')
                set(gca,'YScale','log')
                ylim([1e13 1e18])
            end
            box on
            caxis(cVars{ci,2})

        end
    end
    cb = colorbar;
    cb.Position = [cbX cbY cbW cbH];
    cb.Label.String = cVars{ci,1};
end
% Discrete Vars
%     'C0'
%     'N_nucl','# Nucl. Events'      
%     'Z_nucl','Nucl. Depths'
%     'dZ_nucl' 
%% Decompression rates vs Chamber overpressure
fname = 'Chamber overpressure vs Decompression rates';

% Symbol, code, edge alpha
code_symbols = {'o',0,0.3;
                's',1,0.5;
                '^',2,0.8};
% null_symbol = ''

% Sort by MER for better marker size-plotting
T = sortrows(T,'Q','descend');
msz = [25 35]; % Marker sized for MER

Zt = [6000 2200 1400]; % Break up plots by conduit length

figure('position',[50 100 1500 600],'name',fname)
for zi = 1:length(Zt)
    ax(zi) = subplot(1,length(Zt),zi);
    Tsub = T(T.Z_total==Zt(zi),:);
    for si=1:size(code_symbols,1)
        symb = code_symbols{si,1};
        code = code_symbols{si,2};
        edgealph = code_symbols{si,3};

        % Set marker sizes based on MER
        getCode = Tsub.simpleCode==code;
        size_map = zeros(sum(getCode),1);
        size_map(Tsub.Q(getCode)==1e8) = msz(1);
        size_map(Tsub.Q(getCode)==1e9) = msz(2);
        
        if code==0
            scatter(Tsub.dP(getCode)/1e6,Tsub.dP_dt_bar(getCode)/1e6,size_map,...
                Tsub.Zw(getCode),symb,...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        else
            scatter(Tsub.dP(getCode)/1e6,Tsub.dP_dt_bar(getCode)/1e6,size_map,...
                Tsub.Zw(getCode),symb,'filled',...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        end
        hold on
    end
%     set(gca,'YScale','log')
%     colormap(pasteljet(length(unique(T.Zw))))
    colormap(flipud(pasteljet))
    if zi==length(Zt)
        cb = colorbar;
        cb.Label.String = 'Water Depth (m)';
    end
    ylabel('dP/dt (MPa s^{-1})')
    xlabel('Chamber Overpressure (MPa)')
    grid on
    title(sprintf('Chamber depth: %.0f m',Zt(zi)))
end
linkaxes(ax,'xy')

% ylim([1e13 1e17])
% xlim([0 7.5])
%% Decompression rates vs bubble number densities
fname = 'Bubble number densities vs Decompression rate';

% Symbol, code, edge alpha
code_symbols = {'o',0,0.3;
                's',1,0.5;
                '^',2,0.8};
% null_symbol = ''

% Sort by MER for better marker size-plotting
T = sortrows(T,'Q','descend');
msz = [25 45]; % Marker sized for MER

Zt = [6000 2200 1400]; % Break up plots by conduit length

figure('position',[50 100 1500 600],'name',fname)
for zi = 1:length(Zt)
    ax(zi) = subplot(1,length(Zt),zi);
    Tsub = T(T.Z_total==Zt(zi),:);
    for si=1:size(code_symbols,1)
        symb = code_symbols{si,1};
        code = code_symbols{si,2};
        edgealph = code_symbols{si,3};

        % Set marker sizes based on MER
        getCode = Tsub.simpleCode==code;
        size_map = zeros(sum(getCode),1);
        size_map(Tsub.Q(getCode)==1e8) = msz(1);
        size_map(Tsub.Q(getCode)==1e9) = msz(2);
        
        if code==0
            scatter(Tsub.dP_dt_bar(getCode)/1e6,Tsub.BND_f(getCode),size_map,...
                Tsub.Zw(getCode),symb,...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        else
            scatter(Tsub.dP_dt_bar(getCode)/1e6,Tsub.BND_f(getCode),size_map,...
                Tsub.Zw(getCode),symb,'filled',...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        end
        hold on
    end
    set(gca,'YScale','log')
%     colormap(pasteljet(length(unique(T.Zw))))
    colormap(flipud(pasteljet))
    if zi==length(Zt)
        cb = colorbar;
        cb.Label.String = 'Water Depth (m)';
    end
    xlabel('dP/dt (MPa s^{-1})')
    ylabel('BND (m^{-3})')
    grid on
    title(sprintf('Chamber depth: %.0f m',Zt(zi)))
end
linkaxes(ax,'xy')

ylim([1e13 1e17])
xlim([0 7.5])

%% Final magma/pyroclast porosity vs Decompression rate
fname = 'Final magma/pyroclast porosity vs Decompression rate';
% Symbol, code, edge alpha
code_symbols = {'o',0,0.3;
                's',1,0.5;
                '^',2,0.8};
% null_symbol = ''

% Sort by MER for better marker size-plotting
T = sortrows(T,'Q','descend');
msz = [25 45]; % Marker sized for MER

Zt = [6000 2200 1400]; % Break up plots by conduit length

figure('position',[50 100 1500 600],'name',fname)
for zi = 1:length(Zt)
    ax(zi) = subplot(1,length(Zt),zi);
    Tsub = T(T.Z_total==Zt(zi),:);
    for si=1:size(code_symbols,1)
        symb = code_symbols{si,1};
        code = code_symbols{si,2};
        edgealph = code_symbols{si,3};

        % Set marker sizes based on MER
        getCode = Tsub.simpleCode==code;
        size_map = zeros(sum(getCode),1);
        size_map(Tsub.Q(getCode)==1e8) = msz(1);
        size_map(Tsub.Q(getCode)==1e9) = msz(2);
        
        if code==0
            scatter(Tsub.dP_dt_bar(getCode)/1e6,Tsub.phi_f(getCode),size_map,...
                Tsub.Zw(getCode),symb,...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        else
            scatter(Tsub.dP_dt_bar(getCode)/1e6,Tsub.phi_f(getCode),size_map,...
                Tsub.Zw(getCode),symb,'filled',...
                'MarkerFaceAlpha',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerEdgeAlpha',edgealph)
        end
        hold on
    end
%     set(gca,'YScale','log')
%     colormap(pasteljet(length(unique(T.Zw))))
    colormap(flipud(pasteljet))
    if zi==length(Zt)
        cb = colorbar;
        cb.Label.String = 'Water Depth (m)';
    end
    xlabel('dP/dt (MPa s^{-1})')
    ylabel('Porosity')
    grid on
    title(sprintf('Chamber depth: %.0f m',Zt(zi)))
end
linkaxes(ax,'xy')

% ylim([1e13 1e17])
% xlim([0 7.5])

%% Number, depth, etc associated with nucleation events

if printFigs
    % Figure print dimensions
%     fs = 8;
%     lfs = 7;
%     ppads = [0.9 0.105 0.16 0.06];
%     Qsz = [25 50];
%     oDims = [18 10];
    
% Temp for legibility while drafting
    fs = 8;
    lfs = 7;
    ppads = [0.09 0.105 0.16 0.06];
    Qsz = [25 50];
    oDims = [30 15];
else
    fs = 12;
    lfs = 11;
    ppads = [0.05 0.1 0.12 0.05];
    Qsz = [35 90];
end


nc = 3;
nr = 1;
dx = 0.05;
dy = 0;
cbpos = [0.015+(1-ppads(2)) ppads(3) 0.02 1-sum(ppads(3:4)) ];

Pscale = 1/1e6;

Q  = unique(T.Q);
Zw = unique(T.Zw);
Ztot = unique(T.Z_total);
sCodes = unique(T.simpleCode);

codeMarkers = {'o','s','h'};
N_nucl_edge_colors = {'k','k','r'};
codeAlpha = [.5 1 1];
Nuc_alpha = [0.5 1.0];

figure('position',[50 100 1600 600],'name','Nucleation events, explosivity, and water depth')

for ai = 1:nc
    ax(ai) = tightSubplot(nr,nc,ai,dx,dy,ppads);
    
    subsetIdx = (T.Z_total ==Ztot(ai));
    Tsub = sortrows(T(subsetIdx,:),'Q','descend');
    
    % Size by W
    msz = zeros(size(Tsub,1),1);
    msz(Tsub.Q==1e8) = Qsz(1);
    msz(Tsub.Q==1e9) = Qsz(2);
    
    
% --- Psat vs phi ----
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & Tsub.N_nucl<2);
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.phi_f(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci))
%         hold on
%     end
%     
%     % Do 2nd nucleation events 2nd
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & Tsub.N_nucl==2);
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.phi_f(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerEdgeColor','r',...
%             'MarkerEdgeAlpha',codeAlpha(ci))
%     end
    
% --- Psat vs Zw ----
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & Tsub.N_nucl<2);
%     scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zw(codeSub),msz(codeSub),Tsub.phi_f(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerEdgeColor','k',...
%         'MarkerEdgeAlpha',codeAlpha(ci))
%     hold on
% end
% 
% % Do 2nd nucleation events 2nd
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & Tsub.N_nucl==2);
%     scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zw(codeSub),msz(codeSub),Tsub.phi_f(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerEdgeColor','r',...
%         'MarkerEdgeAlpha',codeAlpha(ci))
% end

% --- Psat vs n0 ---- MONEY
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerFaceAlpha',Nuc_alpha(1))
%         %             'MarkerEdgeColor','k',...
% %             'MarkerEdgeAlpha',codeAlpha(ci),...
%         hold on
%     end
%     
%     % Do 2nd nucleation events 2nd
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci),...
%             'MarkerFaceAlpha',Nuc_alpha(2))
%     end
    
% --- Psat vs Zf ----
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zf(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerFaceAlpha',Nuc_alpha(1))
%         %             'MarkerEdgeColor','k',...
% %             'MarkerEdgeAlpha',codeAlpha(ci),...
%         hold on
%     end
%     
%     % Do 2nd nucleation events 2nd
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zf(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarkers{ci},...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci),...
%             'MarkerFaceAlpha',Nuc_alpha(2))
%     end
    
% --- Psat vs dPdt ----
yname = 'Avg. Decompression Rate (MPa s^{-1})';
xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
oname = 'sweepScatter_dPdt_vs_Psat';
legax = 1;
for ci = 1:length(sCodes)
    codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
    scatter(Tsub.Psat(codeSub).*Pscale,Tsub.dP_dt_bar(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
        'filled',codeMarkers{ci},...
        'MarkerFaceAlpha',Nuc_alpha(1))
    %             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci),...
    hold on
end
% Do 2nd nucleation events 2nd
for ci = 1:length(sCodes)
    codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
    scatter(Tsub.Psat(codeSub).*Pscale,Tsub.dP_dt_bar(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
        'filled',codeMarkers{ci},...
        'MarkerEdgeColor',[0.3 0.3 0.3],...
        'MarkerEdgeAlpha',codeAlpha(ci),...
        'MarkerFaceAlpha',Nuc_alpha(2))
end

% --- Psat vs LAST Z_nucl ----
% yname = 'Final nucleation depth (m)';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_Z_nucl_vs_Psat';
% legax = 1;
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%     Z_nucl = cellfun(@min,Tsub.Z_nucl);
%     
%     scatter(Tsub.Psat(codeSub).*Pscale,Z_nucl(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerFaceAlpha',Nuc_alpha(1))
%     %             'MarkerEdgeColor','k',...
% %             'MarkerEdgeAlpha',codeAlpha(ci),...
%     hold on
% end
% % Do 2nd nucleation events 2nd
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%     Z_nucl = cellfun(@min,Tsub.Z_nucl);
%     
%     scatter(Tsub.Psat(codeSub).*Pscale,Z_nucl(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerEdgeColor','k',...
%         'MarkerEdgeAlpha',codeAlpha(ci),...
%         'MarkerFaceAlpha',Nuc_alpha(2))
% end
% axis tight

    if ai==1
        ylabel(yname)
    end
    title(sprintf('Conduit depth (km b.s.l.) = %.1f',Ztot(ai)))
    xlabel(xname)
    set(gca,'FontSize',fs)
end
cb = colorbar(gca,'location','eastoutside');
cax = [min(Zw) max(Zw)];
caxis(cax);
cb.Label.String = 'Water Depth, Z_w (m)';
cb.Position = cbpos;
colormap(flipud(pasteljet))

% Build a custom legend
lh(1) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{1},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Intrusive');
lh(2) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Effusive');
lh(3) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{3},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Fragmenting');
lh(4) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Faded : 1 Nucleation Event');
lh(5) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(2),'MarkerEdgeColor','k','DisplayName','Bold/border : 2 Nucleation Events');
lh(6) = scatter(ax(legax),nan,nan,min(Qsz),codeMarkers{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Small : Q = 10^8 kg/s');
lh(7) = scatter(ax(legax),nan,nan,max(Qsz),codeMarkers{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(2),'MarkerEdgeColor','k','DisplayName','Large : Q = 10^9 kg/s');
legend(ax(legax),lh,'location','northwest','FontSize',lfs)

if printFigs
    printpdf(oname,FIGURES_DIR,oDims)
end
