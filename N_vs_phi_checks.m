% Load original hydroVolc manuscript conduit results to look at BND/phi
% relationships to get a sense for how to vary inputs appropriately
        % -> immediately after nucleation
        % -> over height z or for varying pressures
clear all; % close all

% dataDir = '~/Kahuna/data/glaciovolc/coupledSweeps/Rowell_etal_hydrovolcanism_scenario_data_2022-02-24';
% sweepFile = 'coupledSweep_2021-09-02_MWIv2_hiLat_n1785_35q_0zw_Reference.mat';

dataDir = '~/Kahuna/data/glaciovolc/coupledSweeps/manuscript_v1_submitted_data';
sweepFile = 'coupledSweep_2021-09-02_MWIv2_hiLat_n1785_35q_0zw.mat';

summFile = ['outputSummary_' sweepFile];

%% Do the thing

load(fullfile(dataDir,sweepFile))
load(fullfile(dataDir,summFile))

%%
% Find end depths of first nucleation event
Q0 = qA.Q0;
Zw = qA.Zw;
logQ0 = log10(Q0);



zNuc = NaN(size(qSet));
pmNuc = zNuc; % pressure after 1st nucleation
ndNuc = zNuc; % BND after 1st nucleation
phiNuc = zNuc; % Porosity after 1st nucleation
h20_0   = zNuc; % Initial dissolved h20
h20_Nuc = zNuc; % Dissolved H20 after 1st nucleation
h20_ex  = zNuc; % Exsolved h20 after 1st nucleation
CsolNuc = zNuc; % Solubility after 1st nucleation

% Properties at depth z=2200 m for comparison
z2.z  = 2200;
z2.pm = NaN(size(qSet));
z2.nd = zNuc;
z2.phi = zNuc;
z2.h20_ex = zNuc; % Total exsolved
z2.h20 = zNuc; % Remaining dissolved
z2.Csol = zNuc; % Solubility
z2.frag = zNuc;

for ii=1:numel(qSet)
    if ~isempty(qSet(ii).cO)
        nucIdx      = find(and(-gradient(qSet(ii).cO.J,qSet(ii).cO.Z)<0,qSet(ii).cO.J<1e-9),1,'first');
        zNuc(ii)    = qSet(ii).cO.Z(nucIdx);
        pmNuc(ii)   = qSet(ii).cO.pm(nucIdx);
        ndNuc(ii)   = qSet(ii).cO.M0(nucIdx); % get BND, phi - then run variations on conduit depth, dP
        phiNuc(ii)  = qSet(ii).cO.porosity(nucIdx);
        h20_0(ii)   = qSet(ii).cO.Cm(1);
        h20_Nuc(ii) = qSet(ii).cO.Cm(nucIdx);
        h20_ex(ii)  = qSet(ii).cO.Cm(1) - qSet(ii).cO.Cm(nucIdx);
        CsolNuc(ii) = qSet(ii).cO.Csol(nucIdx);
%         faafo % Reconstruct BSD here

        % Do a second set at the depth of 2200m?
        [~,z2Idx] = min(abs(qSet(ii).cO.Z-z2.z));
        z2.pm(ii)   = qSet(ii).cO.pm(z2Idx);
        z2.nd(ii)   = qSet(ii).cO.M0(z2Idx); % get BND, phi - then run variations on conduit depth, dP
        z2.phi(ii)  = qSet(ii).cO.porosity(z2Idx);
        z2.h20(ii)  = qSet(ii).cO.Cm(z2Idx);
        z2.h20_ex(ii)  = qSet(ii).cO.Cm(1) - qSet(ii).cO.Cm(z2Idx);
        z2.Csol(ii) = qSet(ii).cO.Csol(z2Idx);
        z2.frag(ii) = qSet(ii).cO.Par.Zf>z2.z;
    end
end

% Calc solubility vs pressure,depth (adjusting for water depth?)
% `-> Use this to get phi0

% Checking
% [axC,lh] = plotConduitOutput([qSet(26,1) qSet(26,10)]);
% plot(axC(5),[1e-10 1e16],zNuc(26,1).*[1 1]/1e3,'--k')

%% ----- Plot results -----
fs = 13;
fscb = 12;
dx = 0.18;
dy = 0.07;
ppads = [.08 .09 0.06 0.02];
nr = 3;
nc = 2;
lw = 2;
lw2 = 1.0;

oDims  = [18 20];
oUnits = 'centimeters';
dpi    = 300;

cfig = figure('position',[50 50 900 1000]);
for ii=1:nr*nc
    cax(ii) = tightSubplot(nr,nc,ii,dx,dy,ppads);
end


% ---------  1st Nuc depth for all? --------- 
axes(cax(1))
imagesc(gca,Q0,Zw,zNuc'/1e3); %shading flat
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Control MER, $Q_0$ (kg/s)','Interpreter','Latex')
ylabel('Water Depth, $Z_e$ (m)','Interpreter','Latex')
colormap(gca,plasma)
cb = colorbar;
cb.Label.String = '1st Nucleation Depth (km)';
cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
apos = get(gca,'Position');
cb.Position = [0.42 apos(2) .016 apos(4)];

%  --------- Initial dissolved h20 --------- 
axes(cax(2))
imagesc(gca,Q0,Zw,h20_0'*1e2); %shading flat
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Control MER, $Q_0$ (kg/s)','Interpreter','Latex')
ylabel('Water Depth, $Z_e$ (m)','Interpreter','Latex')
colormap(gca,plasma)
cb = colorbar;
cb.Label.String = 'Initial dissolved H$_2$O ($wt.\%$)';
cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
apos = get(gca,'Position');
cb.Position = [0.92 apos(2) .016 apos(4)];

%  --------- BND after 1st nuc vs Frac exsolved after 1st nucl --------- 
axes(cax(3))

% -- by MER ---
scatter(ndNuc(:),h20_ex(:)./h20_0(:),30,log10(qA.cI.Q(:)),'filled')
cb = colorbar;
cbl = cb.TickLabels;
for cbi=1:length(cbl)
    cbl{cbi} = sprintf('10^{%s}',cbl{cbi});
end
cb.TickLabels = cbl;
cb.Label.String = 'Adjusted MER $q_c$ (kg/s)';

% -- by Zw --
% scatter(ndNuc(:),phiNuc(:),30,qA.cI.Zw(:),'filled')
% cb = colorbar;
% cb.Label.String = '$Z_w$ (m)';

colormap(gca,plasma)
set(gca,'FontSize',fs)
set(gca,'Xscale','log') %,'XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('BND after 1st nucleation (m$^{-3}$)','Interpreter','Latex')
ylabel({'Fraction exsolved H$_2$0', 'after 1st nucleation'},'Interpreter','Latex')

cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
hold on
grid on
apos = get(gca,'Position');
cb.Position = [0.42 apos(2) .016 apos(4)];

% --------- Fraction exsolved H20 vs Por after 1st nuc --------- 
axes(cax(4))

% -- by MER ---
scatter(h20_ex(:)./h20_0(:),phiNuc(:),30,log10(qA.cI.Q(:)),'filled')
cb = colorbar;
cbl = cb.TickLabels;
for cbi=1:length(cbl)
    cbl{cbi} = sprintf('10^{%s}',cbl{cbi});
end
cb.TickLabels = cbl;
cb.Label.String = 'Adjusted MER $q_c$ (kg/s)';

% -- by Zw --
% scatter(ND_Nuc(:),h20_ex(:)./h20_0(:),30,qA.cI.Zw(:),'filled')
% cb = colorbar;
% cb.Label.String = '$Z_w$ (m)';

colormap(gca,plasma)
set(gca,'FontSize',fs)
set(gca,'FontSize',fs)
ylabel('$\phi$ after 1st nucleation','Interpreter','Latex')
xlabel({'Fraction exsolved H$_2$0', 'after 1st nucleation'},'Interpreter','Latex')

cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
hold on
grid on
apos = get(gca,'Position');
cb.Position = [0.92 apos(2) .016 apos(4)];

% ---------  BND vs pm --------- 
axes(cax(5))

% -- by MER ---
scatter(pmNuc(:),ndNuc(:),30,log10(qA.cI.Q(:)),'filled')
cb = colorbar;
cbl = cb.TickLabels;
for cbi=1:length(cbl)
    cbl{cbi} = sprintf('10^{%s}',cbl{cbi});
end
cb.TickLabels = cbl;
cb.Label.String = 'Adjusted MER $q_c$ (kg/s)';

% -- by Zw --
% scatter(ND_Nuc(:),h20_ex(:)./h20_0(:),30,qA.cI.Zw(:),'filled')
% cb = colorbar;
% cb.Label.String = '$Z_w$ (m)';

colormap(gca,plasma)
set(gca,'FontSize',fs)
set(gca,'FontSize',fs)
xlabel('$P_m$ after 1st nucleation','Interpreter','Latex')
ylabel('BND after 1st nucleation (m$^{-3}$)','Interpreter','Latex')
set(gca,'Yscale','log') %,'XTick',10.^round(min(logQ0):max(logQ0)))

cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
hold on
grid on
apos = get(gca,'Position');
cb.Position = [0.42 apos(2) .016 apos(4)];

% ---------  BND vs Supersaturation --------- 
axes(cax(6))

imagesc(gca,Q0,Zw,(h20_Nuc./CsolNuc)'); %shading flat
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Control MER, $Q_0$ (kg/s)','Interpreter','Latex')
ylabel('Water Depth, $Z_e$ (m)','Interpreter','Latex')
colormap(gca,plasma)
cb = colorbar;
cb.Label.String = 'Supersaturation after 1st nucleation';
cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
apos = get(gca,'Position');
cb.Position = [0.92 apos(2) .016 apos(4)];

%% Properties at 2200m

nr = 2;
nc = 2;
cfig2 = figure('position',[50 50 900 800]);
for ii=1:nr*nc
    cax(ii) = tightSubplot(nr,nc,ii,dx,dy,ppads);
end

%  --------- BND vs Por --------- 
axes(cax(1))

% -- by MER ---
scatter(z2.nd(:),z2.phi(:),30,log10(qA.cI.Q(:)),'filled')
cb = colorbar;
cbl = cb.TickLabels;
for cbi=1:length(cbl)
    cbl{cbi} = sprintf('10^{%s}',cbl{cbi});
end
cb.TickLabels = cbl;
cb.Label.String = 'Adjusted MER $q_c$ (kg/s)';

% -- by Zw --
% scatter(ndNuc(:),phiNuc(:),30,qA.cI.Zw(:),'filled')
% cb = colorbar;
% cb.Label.String = '$Z_w$ (m)';

colormap(gca,plasma)
set(gca,'FontSize',fs)
set(gca,'Xscale','log') %,'XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('BND @ 2.2km (m$^{-3}$)','Interpreter','Latex')
ylabel('$\phi$ @ 2.2km','Interpreter','Latex')

cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
hold on
grid on
apos = get(gca,'Position');
cb.Position = [0.42 apos(2) .016 apos(4)];

% ---------  Frag yet? --------- 
axes(cax(2))

imagesc(gca,Q0,Zw,z2.frag'); %shading flat
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Control MER, $Q_0$ (kg/s)','Interpreter','Latex')
ylabel('Water Depth, $Z_e$ (m)','Interpreter','Latex')
colormap(gca,plasma)
cb = colorbar;
cb.Label.String = 'Frag below 2.2km?';
cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
apos = get(gca,'Position');
cb.Position = [0.92 apos(2) .016 apos(4)];

% ---------  BND vs pm --------- 
axes(cax(3))

% -- by MER ---
scatter(z2.pm(:),z2.nd(:),30,log10(qA.cI.Q(:)),'filled')
cb = colorbar;
cbl = cb.TickLabels;
for cbi=1:length(cbl)
    cbl{cbi} = sprintf('10^{%s}',cbl{cbi});
end
cb.TickLabels = cbl;
cb.Label.String = 'Adjusted MER $q_c$ (kg/s)';

% -- by Zw --
% scatter(ND_Nuc(:),h20_ex(:)./h20_0(:),30,qA.cI.Zw(:),'filled')
% cb = colorbar;
% cb.Label.String = '$Z_w$ (m)';

colormap(gca,plasma)
set(gca,'FontSize',fs)
set(gca,'FontSize',fs)
xlabel('$P_m$ @ 2.2km','Interpreter','Latex')
ylabel('BND 2.2km (m$^{-3}$)','Interpreter','Latex')
set(gca,'Yscale','log') %,'XTick',10.^round(min(logQ0):max(logQ0)))

cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
hold on
grid on
apos = get(gca,'Position');
cb.Position = [0.42 apos(2) .016 apos(4)];

% ---------  Fraction exsolved volatiles --------- 
axes(cax(4))

z2.frag(isnan(z2.frag)) = false;
ex_frac = (z2.h20_ex./h20_0);
ex_frac(logical(z2.frag)) = NaN;

imagesc(gca,Q0,Zw,ex_frac'); %shading flat
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Control MER, $Q_0$ (kg/s)','Interpreter','Latex')
ylabel('Water Depth, $Z_e$ (m)','Interpreter','Latex')
colormap(gca,plasma)
cb = colorbar;
cb.Label.String = 'Fraction exsolved H$_2$O @ 2.2km';
cb.Label.Interpreter = 'Latex';
cb.FontSize = fscb;
apos = get(gca,'Position');
cb.Position = [0.92 apos(2) .016 apos(4)];
% caxis([1.0 1.2])

%% Doing a check on excess volatiles and solubility
qi = 26;
P0 = zeros(length(Zw));
for zi=1:length(Zw)
    P0 = qSet(qi,zi).cO.pm(1);
end
sol0 = meltH2Osolubility(P0,0,qSet(1,1).cO.Par.T);
sol2_2 = meltH2Osolubility(z2.pm(qi,:),0,qSet(1,1).cO.Par.T);

% Predicted h2o exsolution assuming exquilibrium
h2o_ex_pred = (sol0-sol2_2)./sol0;

% Modeled exsolved h2o and supersaturated component
h2o_ex_calc = ex_frac(qi,:);
h2o_ex_sup  = (z2.h20(qi,:)-z2.Csol(qi,:))./h20_0(qi,:);

P_range = (1e5:1e5:1.5e8)'; 
h2oSol = [meltH2Osolubility(P_range,0,qSet(1,1).cO.Par.T) ...
    meltH2Osolubility(P_range,0.025,qSet(1,1).cO.Par.T)...
    meltH2Osolubility(P_range,0.05,qSet(1,1).cO.Par.T)];

figure
subplot(1,2,1)
plot(P_range/1e6,h2oSol)
xlabel('Pressure (MPa)','Interpreter','Latex')
ylabel('H$_2$O Solubility ($wt.\%$)','Interpreter','Latex')
lg = legend({'0.0','0.025','0.05'},'location','northwest');
lg.Title.String = 'CO_2 volatile fraction';
set(gca,'FontSize',fs)
hold on
plot((2400*9.81*2200)/1e6*[1 1],[0 max(h2oSol(:))],'--k')
plot((2400*9.81*6000)/1e6*[1 1],[0 max(h2oSol(:))],'--k')

subplot(1,2,2)
plot([0.35 0.5],[0.35 0.5],'--k')
set(gca,'ColorOrderIndex',1)
hold on
plot(h2o_ex_pred,h2o_ex_calc,'o')
plot(h2o_ex_pred,h2o_ex_calc+h2o_ex_sup,'o')
xlabel('Predicted exolved fraction')
ylabel('Modeled exsolved fraction')
legend({'1-1','Exsolved','Exsolved + supersat'},'location','northwest')
set(gca,'FontSize',fs)



