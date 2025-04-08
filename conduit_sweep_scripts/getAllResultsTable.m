% Build the tabular dataset of parameter sweep outputs
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
%
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
%     phi_0         : INITIAL porosity
%     phi_f         : FINAL porosity (at surfac/end)
%
% C Rowell, Nov 2024

clear all; %close all;
run ../config % Assumes run from local path

% outDir  = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/';
outDir = fullfile(MAIN_DATA_DIR,'mainSweep2/refinedSweep/');

outcomeFile = 'outcomeCodeSummary.mat';


outcomes = {'OutcomeCode','simpleCode'};                            % Conduit OutcomeCodes
cI_pars = {'Q','Z0','Zw','dP','n0_excess','N0','conduit_radius'};   % Scalar input parameters
cO_pars = {'Zf','C0'};                                              % Scalar output parameters
output_pars = {'Cm_0','Psat','t_ascent','dP_dt_bar', 'dP_dt_max',...    % Output fields needing specific retreival methods
    'BND_f','N_nucl','Z_nucl','dZ_nucl','phi_0','phi_f','n0_total'};
calc_pars = {'Z_total'};                                            % Post-processing calculations

n_total = 41*17*21; % All simulations

DeltaP_max = 2e7; % Maximum chamber overpressure to keep
J_threshold = 1e6; % Minimum nucleation rate threshold for tracking nucleation events

% Make QC plots?
qcplot = true;

% ---- save output table? ----
saveTable = true;
% dataSaveDir = '/Users/crrowell/code/research-projects/myojin_knoll/data/';
outputTableName = 'ConduitSweepsDataTableB.mat';

%% Do the thing

load(fullfile(outDir,outcomeFile))
% Sweep names
fn = fieldnames(allOutcomeCodes);
fn = fn(cellfun(@(x) contains(x,'myojin'),fn));

% File names
sweepFiles = dir(outDir);
sweepFiles = {sweepFiles.name}';
sweepFiles = sweepFiles(cellfun(@(x) contains(x,'myojin'),sweepFiles));

%% Make some quick sample plots for the first run
 % referring to Hajimirza ea 2021 figures/details as well

% load(fullfile(outDir,sweepFiles{39}))
% P = linspace(sweepParams.dP.range(1),sweepParams.dP.range(2),sweepParams.dP.n);
% P_keep = P<=DeltaP_max;
% % --> Clip:
% dat = dat(P_keep,:);
% dd = dat(177);

% 2nd nucleation at effusive/explosive transition
% load(fullfile(outDir,sweepFiles{39}))
% dd = dat(7,9);

% 3rd nucleation??
load(fullfile(outDir,sweepFiles{29}))
dd = dat(8,2);

titles = {'Ascent Velocity','Decomp. Rate','Dissolved H_2O','Solubility','Supersaturation','Max. Supersat. Press.','Nucleation Rate','Bubble \# Density','Porosity','Ascent Time'};
xlabs = {'U (m/s)','dP/dt (MPa s^{-1})','Cm (wt.\%)','Csol (wt.\%)','Cm-Csol (wt.\%)','Sat. P (MPa)','J (m^{-3} s^{-1})','BND (m^{-3})','\phi (%)','t (s)'};
exprs = {'dd.cO.U',...
         'dd.cO.dPdt/1e6',...
         'dd.cO.Cm',...
         'dd.cO.Csol',...
         'dd.cO.Cm - dd.cO.Csol',...
         '(dd.cO.psat-dd.cI.pf)/1e6',...
         'dd.cO.J',...
         'dd.cO.M0',...
         'internalPorosity(dd.cO)',...
         'dd.cO.t'};
%          'dd.cO.Csol',...
%         'dd.cO.porosity',...


if qcplot
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
end
%% Build data table
clear temps 

% Initialize table
blank_series = []; % zeros(n_total,1);
all_pars = [outcomes cI_pars cO_pars output_pars calc_pars];
for fi = 1:length(all_pars)
    T.(all_pars{fi}) = blank_series;
end

% loop through all sweep files and pull parameters
disp('Processing sweep files...')
for fi = 1:length(fn)
    sweep = fn{fi};
    whichFile = contains(sweepFiles,sweep);
    sweepFile  = sweepFiles{whichFile};
    fprintf('\t%i/%i %s...\n',fi,length(fn),sweepFile)
    
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
    %     dP_dt_max     : MAX decompression rate
    %     BND_f         : FINAL bubble number density
    %     N_nucl        : Number of discrete nucleation events
    %     Z_nucl        : Onset depth of discrete nucleation events
    %     dZ_nucl       : Difference between fragmentation depth and nucleation depth (if fragmentation occurred)
    %     phi_0         : INITIAL porosity
    %     phi_f         : FINAL porosity (at surface)
    %     n0_total      : Excess gas plus initial solubility = bulk H20 content
    [temps{1:length(output_pars)}] = deal(f_blank);
    temps = cell2struct(temps',output_pars);
    temps.Z_nucl = num2cell(f_blank);
    temps.dZ_nucl = num2cell(f_blank);
    for jj=1:numel(dat)
        if dat(jj).cO.Outcome.Code > -20  && ~any(dat(jj).cO.Z<-10) % Get for non-failed runs

            temps.BND_f(jj)     = dat(jj).cO.M0(end);
            temps.dP_dt_bar(jj) = (dat(jj).cO.pm(1) - dat(jj).cO.pm(end))/dat(jj).cO.t(end); % Taken as the negative
            temps.dP_dt_max(jj) = max(abs(dat(jj).cO.dPdt));
            temps.Cm_0(jj)      = dat(jj).cO.Cm(1);
            temps.Psat(jj)      = dat(jj).cO.psat(1) - dat(jj).cI.pf;
            temps.t_ascent(jj)  = dat(jj).cO.t(end);
            temps.phi_0(jj)     = dat(jj).cO.porosity(1);
            temps.n0_total(jj)  = (dat(jj).cI.n0_excess+1).*dat(jj).cO.Cm(1);
            
            try
                porosity_in = internalPorosity(dat(jj).cO);
            catch ME
                keyboard
            end
            temps.phi_f(jj) = porosity_in(end);

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
            % Try:
            % - case where initial non-zero BND was set, and J is not
            % already above threshold at onset : TREAT THIS AS ADDITIONAL
            % (PRESCRIBED) NUCLEATION EVENT
            if ~nuc(1) && dat(jj).cI.n0_excess>0
                temps.N_nucl(jj) = temps.N_nucl(jj) + 1;
                temps.Z_nucl{jj} = [temps.Z_nucl{jj} NaN];
                temps.dZ_nucl{jj} = [temps.dZ_nucl{jj} NaN];
            end
            %
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
        TT = struct2table(rmfield(T,'Z_total'));
        icheck = 177;
        dd = dat(icheck);
        bw = bwconncomp(dd.cO.J>J_threshold);
        Z_nucl_test = dd.cO.Z(bw.PixelIdxList{1}(1));
        fprintf('\nVALIDATION CHECK:\n\tNucleation depth from raw:\t %.2f m\n\tNucleation depth from table:\t %.2f m\n\n',Z_nucl_test,TT.Z_nucl{icheck})
        
        if qcplot
            % Quick QC plot to show nucleation depths for a single sweep
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
            title(sprintf("Nucleation event and fragmentation depths: '%s'",strrep(run_name,'_','\_')))
        end
    end
end

T.Z_total = T.Zw + T.Z0;
for ii = 1:length(T.Z_nucl)
    T.dZ_nucl{ii} = T.Z_nucl{ii} - T.Zf(ii);
end

T = struct2table(T);
%% 


%% Clean simplified outcome codes
% -> using same methods as in process_conduit_outcomes.process_outcome_codes (py)



% Filter out rows with error results since we can't use these
removeRows = T.simpleCode <= -2;
numErrorRowsRemoved = sum(removeRows);
l_before = size(T,1);
T(removeRows,:) = [];
l_after = size(T,1);
fprintf('Removing error codes: num table rows =  %i --> %i\n',l_before,l_after)
%% To table and save
if saveTable
    ofile = fullfile(DATA_DIR, outputTableName);
    fprintf('Saving output file:\n\t%s\n',ofile)
    save(ofile,'T','DeltaP_max','J_threshold','numErrorRowsRemoved','sweepFiles')
end
%% Function to calculate pyroclast internal porosity after fragmentation
function porosity = internalPorosity(out)
% out = conduit output struct (ie dat.cO)
%
    porosity = zeros(size(out.porosity));
    
    por_non_zero = out.porosity~=0;

    vg_in = 4/3 * pi * out.M3(por_non_zero);
    [rho_in,~] = EoS_H2O_2(out.pg(por_non_zero),out.Par.T);
    mg_in = rho_in.*vg_in;
    
    mg_out = out.mg(por_non_zero) - mg_in;
    pg_out = out.pm(por_non_zero);

    [rho_out,~] = EoS_H2O_2(pg_out,out.Par.T);
    vg_out = mg_out ./ rho_out;

    % THESE 'porosities' measure gas volume fraction of the total volume,
    % which includes escaped gas after fragmentation
    porosity_in = (vg_in) ./ (1+vg_in+vg_out);
    porosity_out = (vg_out) ./ (1+vg_in+vg_out);
    
    % '-> Therefore renormalize to pyroclast + internal porosity volume
    porosity(por_non_zero) = porosity_in ./ (1-porosity_out);

end