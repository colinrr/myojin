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

printFigs = false;
%%  Get table subset for plotting

load(fullfile(DATA_DIR,'ConduitSweepsDataTableB.mat'),'T')
allZw = unique(T.Zw);
all_n0e = unique(T.n0_excess);

rows_to_get = false(size(T,1),1);
for qi = 1:length(runQ)
    Tsub = T((T.Q==runQ(qi)),:);
    
    for si = 1:length(filterScenarios.Z_total)
        
        
        for zwi = 1:length(allZw)
            subsetIdx = (T.Q == runQ(qi) ...
                & T.Z_total == filterScenarios.Z_total(si) ...
                & T.dP == filterScenarios.dP(si) ...
                & T.Zw == allZw(zwi));

            n0t_subset = T(subsetIdx,:).n0_total;
            [n0_vals(si,zwi),n0e_idx(si,zwi)] = closest(n0_target,n0t_subset);
            % --> Now use the n0 index to get which n0_excess value to
            % pull. This has to be sample by sample, so initialized values
            % above
            this_row = subsetIdx;% & T.n0_excess == all_n0e(n0e_idx(si,zwi));
            rows_to_get = rows_to_get | this_row;
        end
        

    end
    % General filter
%     subsetIdx = (Tsub.Z_total == filterScenarios.Z_total(si) & Tsub.dP == filterScenarios.dP(si));
%     pause(0.1)
end

T = T(rows_to_get,:);

%% Begin the plots

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
    
    % Secondary filter to get just a few key model results
    
    
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
% yname = 'Total gass mass fraction';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_n_tot_vs_Psat';
% legax = 1;
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
% yname = 'Fragmentation Depth (m)';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_Zf_vs_Psat';
% legax = 1;
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         [codeMarker,~] = getCodeMarker(sCodes(ci));
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zf(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarker,...
%             'MarkerFaceAlpha',Nuc_alpha(1))
%         %             'MarkerEdgeColor','k',...
% %             'MarkerEdgeAlpha',codeAlpha(ci),...
%         hold on
%     end
%     
%     % Do 2nd nucleation events 2nd
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         [codeMarker,~] = getCodeMarker(sCodes(ci));
%         scatter(Tsub.Psat(codeSub).*Pscale,Tsub.Zf(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarker,...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci),...
%             'MarkerFaceAlpha',Nuc_alpha(2))
%     end
%     
% --- Psat vs dPdt ----
% yname = 'Avg. Decompression Rate (MPa s^{-1})';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_dPdt_vs_Psat';
% legax = 1;
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%     scatter(Tsub.Psat(codeSub).*Pscale,Tsub.dP_dt_bar(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerFaceAlpha',Nuc_alpha(1))
%     %             'MarkerEdgeColor','k',...
% %             'MarkerEdgeAlpha',codeAlpha(ci),...
%     hold on
% end
% % Do 2nd nucleation events 2nd
% for ci = 1:length(sCodes)
%     codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%     scatter(Tsub.Psat(codeSub).*Pscale,Tsub.dP_dt_bar(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%         'filled',codeMarkers{ci},...
%         'MarkerEdgeColor',[0.3 0.3 0.3],...
%         'MarkerEdgeAlpha',codeAlpha(ci),...
%         'MarkerFaceAlpha',Nuc_alpha(2))
% end

% --- Psat vs LAST Z_nucl ----
% yname = 'Final nucleation depth (m)';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% % oname = 'sweepScatter_Z_nucl_vs_Psat';
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

% --- Psat vs n0 ---- MONEY - tight filter for just a few points
% n_total = (T.n0_excess(codeSub)+1).*T.Cm_0(codeSub);
% subsetIdx = (T.Z_total ==Ztot(ai) & ...
%             and(n_total>=0.55,n_total<=0.065) & ...
%             T.dP );

yname = 'Total gass mass fraction';
xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
oname = 'sweepScatter_n_tot_vs_Psat_filtered';
legax = 1;
    for ci = 1:length(sCodes)
        codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
        [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
        scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
            'filled',codeMarker,...
            'MarkerFaceAlpha',Nuc_alpha(1))
        %             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',codeAlpha(ci),...
        hold on
    end
    
    % Do 2nd nucleation events 2nd
    for ci = 1:length(sCodes)
        codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
        [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
        scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
            'filled',codeMarker,...
            'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',ca,...
            'MarkerFaceAlpha',Nuc_alpha(2))
    end

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
lh(1) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(1),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Intrusive');
lh(2) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(2),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Effusive');
lh(3) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(3),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Fragmenting');
lh(4) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(2),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Faded : 1 Nucleation Event');
lh(5) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(2),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(2),'MarkerEdgeColor','k','DisplayName','Bold/border : 2 Nucleation Events');
lh(6) = scatter(ax(legax),nan,nan,min(Qsz),getCodeMarker(2),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(1),'MarkerEdgeAlpha',0,'DisplayName','Small : Q = 10^8 kg/s');
lh(7) = scatter(ax(legax),nan,nan,max(Qsz),getCodeMarker(2),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',Nuc_alpha(2),'MarkerEdgeColor','k','DisplayName','Large : Q = 10^9 kg/s');
legend(ax(legax),lh,'location','northwest','FontSize',lfs)

if printFigs
    printpdf(oname,FIGURES_DIR,oDims)
end
