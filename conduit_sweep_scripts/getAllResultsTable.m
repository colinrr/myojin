% Build a tabular dataset of sweep output
%  GET:
%   -> cI: Q, n0, Zw, Z0 (conduit length)
%   -> Ztotal
%   -> # nucleaton events and their (onset) depths
%   -> Ascent time
%   -> Time avg. decompression rate over ascent time
%   -> Outcome Code - full and simple
%   -> BND at fragmentation or surface - whichever is first
%   -> Pressure, depth, dissolved H20, BND at {nucleation? fragmentation? values matching myojin?}
%       -> For dissolved H20, also get equilibrium (solubility) value? (depends on depth choice)

clear all; close all;

outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
dataSaveDir = '/Users/crrowell/code/research-projects/myojin_knoll/data';

outcomeFile = 'outcomeCodeSummary.mat';

outcomes = {'OutcomeCode','simpleCode'};                            % Conduit OutcomeCodes
cI_pars = {'Q','Z0','Zw','dP','n0_excess','N0','conduit_radius'};   % Scalar input parameters
cO_pars = {'Zf','C0'};                                              % Scalar output parameters
output_pars = {'Cm_0','Psat','t_ascent','dP_dt_bar',...    % Output fields needing specific retreival methods
    'BND_f','N_nucl','Z_nucl','dZ_nucl'};
calc_pars = {'Z_total'};                                            % Post-processing calculations

n_total = 41*17*21; % All simulations


DeltaP_max = 2e7; % Maximum chamber overpressure to keep
J_threshold = 1e6; % Minimum nucleation rate threshold for tracking nucleation events


%% Do the thing

load(fullfile(outDir,outcomeFile))
% Sweep names
fn = fieldnames(allOutcomeCodes);
fn = fn(cellfun(@(x) contains(x,'myojin'),fn));

% File names
ff = dir(outDir);
ff = {ff.name}';
ff = ff(cellfun(@(x) contains(x,'myojin'),ff));

%% Make some quick sample plots for the first run
 % referring to Hajimirza ea 2021 figures/details as well

load(fullfile(outDir,ff{1}))
P = linspace(sweepParams.dP.range(1),sweepParams.dP.range(2),sweepParams.dP.n);
P_keep = P<=DeltaP_max;
% --> Clip:
dat = dat(P_keep,:);

dd = dat(177);

titles = {'Ascent Velocity','Decomp. Rate','Dissolved H$_2$O','Solubility','Supersaturation','Max. Supersat. Press.','Nucleation Rate','Bubble \# Density','Ascent Time'};
xlabs = {'U (m/s)','dP/dt (MPa s$^{-1}$)','Cm ($wt.\%$)','Csol ($wt.\%$)','Cm-Csol ($wt.\%$)','Sat. P (MPa)','J (m$^{-3}$ s$^{-1}$)','BND (m$^{-3}$)','t (s)'};
exprs = {'dd.cO.U',...
         'dd.cO.dPdt/1e6',...
         'dd.cO.Cm',...
         'dd.cO.Csol',...
         'dd.cO.Cm - dd.cO.Csol',...
         '(dd.cO.psat-dd.cI.pf)/1e6',...
         'dd.cO.J',...
         'dd.cO.M0',...
         'dd.cO.t'};
%          'dd.cO.Csol',...

ne = length(exprs);
figure('position',[50 200 1300 600])
for ai = 1:ne
    ax(ai) = tightSubplot(1,ne,ai);
    plot(eval(exprs{ai}),dd.cO.Z);
    set(gca,'ydir','reverse'); 
    grid on; 
%     axis tight;
    xlabel(xlabs{ai})
    title(titles{ai})
    if ai==1
        ylabel('Z (m)')
    else
        set(gca,'YTickLabel',[])
    end
    if contains(exprs{ai},'M0')
        set(gca,'XScale','log')
        xlim([1e5 1e17])
    elseif contains(exprs{ai},'J')
        set(gca,'XScale','log')
        xlim([1e5 1e20])
    end
    ylim([0 dd.cI.Z0])
end
linkaxes(ax,'y')
%% Build data table
clear temps 

% Initialize table
blank_series = []; % zeros(n_total,1);
all_pars = [outcomes cI_pars cO_pars output_pars calc_pars];
for fi = 1:length(all_pars)
    T.(all_pars{fi}) = blank_series;
end

% loop through all sweep files and pull parameters
for fi = 1:length(fn)
    sweep = fn{fi};
    whichFile = contains(ff,sweep);
    sweepFile  = ff{whichFile};
    
    load(fullfile(outDir,sweepFile))
    [a,b] = regexp(sweepFile,'myojin_Q\d+_Z0\d+_Zw\d+');
    run_name = sweepFile(a:b);
    
    % Truncate to desired max Delta P
    P = linspace(sweepParams.dP.range(1),sweepParams.dP.range(2),sweepParams.dP.n);
    P_keep = P<=DeltaP_max;
    % --> Clip:
    dat = dat(P_keep,:);
    modelFail = outcomeCodes.modelFail(P_keep,:);
    
    % ---- Get outcome codes ---- 
    oc = allOutcomeCodes.(run_name)(P_keep,:);
    sc = simplePlotIndex.(run_name)(P_keep,:);
    T.OutcomeCode = [T.OutcomeCode; oc(:)];
    T.simpleCode  = [T.simpleCode; sc(:)];
    
    % ---- Input parameters ---- 
    cI_temp = [dat(:).cI]';
    for ii = 1:length(cI_pars)
        T.(cI_pars{ii}) = [T.(cI_pars{ii}); [cI_temp.(cI_pars{ii})]'];
    end
    
    % ---- Output parameters ---- 
    f_blank = nan(numel(dat),1);
    for ii = 1:length(cO_pars)
        % Many of these fail, so loop with brute force
        f_temp = f_blank;
        for jj=1:numel(dat)
            if dat(jj).cO.Outcome.Code > -20 % Get for non-failed runs
                f_temp(jj) = dat(jj).cO.Par.(cO_pars{ii});
            end
        end
        T.(cO_pars{ii}) = [T.(cO_pars{ii}); f_temp];
    end
    
    % ---- Calculated output fields ---- 
    %     Z_total       : Conduit length plus water depth
    %     Csol          : NOT CURRENTLY USED. Can revisit if solublity info
    %                       is needed to compare with melt inclusion data
    %     Cm_0          : INITIAL dissolved H20 wt%
    %     Psat          : Max. supersaturation pressure = initial saturation pressure - p_final
    %     t_ascent      : Total modeled ascent time
    %     dP_dt_bar     : AVERAGE decompression rate over modeled rise time
    %     BND_f         : FINAL bubble number density
    %     N_nucl        : Number of discrete nucleation events
    %     Z_nucl        : Onset depth of discrete nucleation events
    %     dZ_nucl       : Difference between fragmentation depth and nucleation depth (if fragmentation occurred)
    [temps{1:length(output_pars)}] = deal(f_blank);
    temps = cell2struct(temps',output_pars);
    temps.Z_nucl = num2cell(f_blank);
    temps.dZ_nucl = num2cell(f_blank);
    for jj=1:numel(dat)
        if dat(jj).cO.Outcome.Code > -20 % Get for non-failed runs

            temps.BND_f(jj)     = dat(jj).cO.M0(end);
            temps.dP_dt_bar(jj) = (dat(jj).cO.pm(1) - dat(jj).cO.pm(1))/dat(jj).cO.t(end);
            temps.Cm_0(jj)      = dat(jj).cO.Cm(1);
            temps.Psat(jj)      = dat(jj).cO.psat(1) - dat(jj).cI.pf;
            temps.t_ascent(jj)  = dat(jj).cO.t(end);

            % Nucleation events
            nuc = dat(jj).cO.J >= J_threshold;
            bw = bwconncomp(nuc);
            temps.N_nucl(jj)    = bw.NumObjects;
            if bw.NumObjects > 0
                for ni = 1:bw.NumObjects
                    z_nucl(ni) = dat(jj).cO.Z(bw.PixelIdxList{ni}(1));
        %             dz_nucl(ni) = z_nucl(ni) - f_temp.Zf(jj);
                end
                temps.Z_nucl{jj} = z_nucl;

            end
    %         temps.dz_nucl{jj} = dz_nucl;
            clear z_nucl
        end
    end
    for ii=1:length(output_pars)
        T.(output_pars{ii}) = [T.(output_pars{ii}); temps.(output_pars{ii})];
    end
    
    clear temps
    
    % VALIDITY CHECK - only works where dat is from the first file (ff{1})
    if fi==1
        TT = struct2table(T);
        icheck = 177;
        dd = dat(icheck);
        bw = bwconncomp(dd.cO.J>J_threshold);
        Z_nucl_test = dd.cO.Z(bw.PixelIdxList{1}(1));
        fprintf('Nucleation depth from raw:\t %.2f m\nNucleation depth from table:\t %.2f m\n',Z_nucl_test,TT.Z_nucl{icheck})
    end
end

T.Z_total = T.Zw + T.Z0;
for ii = 1:length(T.Z_nucl)
    T.dZ_nucl{ii} = T.Z_nucl{ii} - T.Zf(ii);
end

%% Quick QC p to show nucleation depths
figure
for jj=1:numel(dat)
    if dat(jj).cO.Outcome.Code>-20
        nuc = dat(jj).cO.J>=1e6;
        zz = dat(jj).cO.Z;
        zz(~nuc) = nan;
        plot(ones(size(nuc)).*jj, zz, '-b')
        hold on
    end
end
plot(1:numel(dat), T.Zf, 'xr')
ll(1) = plot(nan,nan,'-b');
ll(2) = plot(nan,nan,'xr');
legend(ll,{'Nucleation Events','Fragmentation Depth'})
xlabel('Simulation \#')
ylabel('Depth (m)')

set(gca,'YDir','reverse')
title('Nucleation event and fragmentation depths')

%% Clean simplified outcome codes
% -> using same methods as in process_conduit_outcomes.process_outcome_codes (py)

%% To table and save

%% Build