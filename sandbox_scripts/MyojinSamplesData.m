%% Make Myojoin sample data table
run ../config

clear my
% Myojin decompression rate/textural data
my.ID = ["B17","F19","B2","C16","B17","A24","C18","G22","E18","D10","A13"]';
my.h2o_wtp = [5.9 6.5 3.0 6.6 6.5 6.4 6.4 6.4 6.4 6.4 6.4]';
my.nucl_P_Pa  = [191 234 49 243 238 226 225 228 226 227 228]' * 1e6; % MPa to Pa
my.dP_dt_Pas  = [1.3 4.4 3.9 3.1 4.5 6.9 5.8 7 7.2 8.5 4.6]' * 1e6; % MPa to Pa
my.Nv = [2.07e4 2.03e5 7.51e4 1.14e5 1.87e5 7.28e5 3.46e5 9.15e5 1.18e6 8.72e5 7.91e5]' *1e9; % mm^-3 to m^-3
my.Nv_melt_corr = [1.23e8 1.09e9 3.4e8 4.51e8 6.58e8 1.99e9 8.63e8 1.85e9 3.18e9 1.38e9 1.02e9]' *1e6; % cm^-3 to m^-3
my.Nv_lt_10um = [1.7e4 6.7e5 3.2e5 4.0e5 6.8e5 1.3e6 1.0e6 1.3e6 1.4e6 1.8e6 7.1e5]' *1e9; % mm^-3?? to m^-3
my.phi = [.836 .817 .784 .753 .721 .641 .529 .514 .401 .38 .238]';
my.depth_2600 = my.nucl_P_Pa ./ (9.81 * 2600);

my = struct2table(my);
my.Properties.VariableUnits = ["","wt.%","Pa","Pa s^{-1}","m^{-3}","m^{-3}","m^{-3}","","m"];
my.Properties.VariableDescriptions = [
    "Sample ID"
    "Water content (final)"
    "Nucleation Pressure - Water content (final)"
    "Decompression rate (calculated)"
    "Total bubble number density"
    "Bubble number density (melt corrected)"
    "Vesicles < 10 microns, corrected for vesicularity"
    "Vesicularity fraction"
    "Nucleation depth for rock density = 2600 kg/m^3"
    ];

%% Display the table and properties
disp(my)
varDescriptions = table(string(my.Properties.VariableNames'), ...
    string(my.Properties.VariableUnits'), ...
    string(my.Properties.VariableDescriptions'),...
    'VariableNames',["Variable","Unit","Description"]);
disp(varDescriptions)

save(fullfile(DATA_DIR,'MyojinSampleDataTable.mat'),'my','varDescriptions')
%% Solubility model check
% Horizontal black lines show Rebecca's reported melt inclusion H2O range
% Vertical black lines are the 3-5 km depth range



% Myojin samples: CO2 suggested to have negligible influence on H2O solubility, so ignoring that for now
co2_ppm = 0; 

rho_rock = [2400 2600 2800];

lstyles = {'-' , '--'};

co2_mole_frac = co2_ppm ./ 1e6;
T_mag = 850+273.15;

nz = 1201;
Z = linspace(0,12000,nz)'; % Depth range to 10 km

P = 9.81 .* rho_rock .* Z;


%% Plot solubility vs pressure depth curves for myojin samples
sol = zeros(nz,length(rho_rock),length(co2_mole_frac));
for co = 1:length(co2_mole_frac)
    sol(:,:,co) = meltH2Osolubility(P,co2_mole_frac(co),T_mag);
end


figure
subplot(2,1,1)
lh1 = plot(P(:,1)/1e6,sol(:,1,1)*100);
xlabel('P (MPa)')
ylabel('H_2O Solubility (wt.%)')
hold on
lh2 =  plot(xlim'.*[1 1],[4.2 6.5].*[1;1],'--k');
lh3 = plot(xlim,6.*[1;1],'k');
sh = scatter(my.nucl_P_Pa/1e6, my.h2o_wtp, 'sk');
legend([lh1 lh2(1) lh3 sh],["Liu et al. (2005) equation","Melt inclusion range","Melt inclusion mean","Myojin detailed samples"],'location','southeast')


subplot(2,1,2)
for co = 1:length(co2_mole_frac)
    ch{co} = plot(Z,sol(:,:,co)*100,lstyles{co});
    hold on
    set(gca,'ColorOrderIndex',1)
end
xlabel('Depth (m)')
ylabel('H_2O Solubility (wt.%)')
% plot([3e3 5e3].*[1;1],ylim'.*[1 1],'--k')
grid on
lh2 =  plot(xlim'.*[1 1],[4.2 6.5].*[1;1],'--k');
lh3 = plot(xlim,6.*[1;1],'k');
sh = scatter(my.depth_2600, my.h2o_wtp, 'sk');
legend([ch{1}; lh2(1); lh3],{'\rho_{rock} = 2400 kg/m^3','\rho_{rock} = 2600 kg/m^3','\rho_{rock} = 2800 kg/m^3',"Melt inclusion range","Melt inclusion mean","Myojin samples, \rho_{rock} = 2600 kg/m^3"},'location','southeast')

