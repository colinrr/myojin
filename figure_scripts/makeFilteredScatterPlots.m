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
filterScenarios.Z_total = [1400 2200 6000 1400 2200 6000];
filterScenarios.dP = [20 16 10 20 16 10] * 1e6;
filterScenarios.Q = [1e8 1e8 1e8 1e9 1e9 1e9];

runQ = [1e9]; % 

n0_targets = [0.04 0.05 0.06];
n0_tol = [0.002]; % Discard values more than this distance from target

xvar = 'Psat';
yvar = 'n0_total';
common_lines_var = 'Zw';

printFigs = false;
%%  Get table subset for plotting

load(fullfile(DATA_DIR,'ConduitSweepsDataTableB.mat'),'T')


%% Filter v1
% rows_to_get = false(size(T,1),1);
% for qi = 1:length(runQ)
%     Tsub = T((T.Q==runQ(qi)),:);
%     
%     for si = 1:length(filterScenarios.Z_total)
%         
%         
%         for zwi = 1:length(allZw)
%             subsetIdx = (T.Q == runQ(qi) ...
%                 & T.Z_total == filterScenarios.Z_total(si) ...
%                 ... & ismember(T.dP, filterScenarios.dP) ... 
%                 & T.dP == filterScenarios.dP(si) ...
%                 & T.Zw == allZw(zwi));
% 
%             n0t_subset = T(subsetIdx,:).n0_total;
%             [n0_vals(si,zwi),n0e_idx(si,zwi)] = closest(n0_target,n0t_subset);
%             % --> Now use the n0 index to get which n0_excess value to
%             % pull. This has to be sample by sample, so initialized values
%             % above
%             this_row = subsetIdx; % & T.n0_excess == all_n0e(n0e_idx(si,zwi));
%             rows_to_get = rows_to_get | this_row;
%         end
%         
% 
%     end
%     % General filter
% %     subsetIdx = (Tsub.Z_total == filterScenarios.Z_total(si) & Tsub.dP == filterScenarios.dP(si));
% %     pause(0.1)
% end
% 
% T = T(rows_to_get,:);

%% Filter v2
allZw = unique(T.Zw);
all_n0e = unique(T.n0_excess);
all_dP = unique(T.dP);
all_Q = unique(T.Q);
all_Zt = unique(T.Z_total);

nf = length(filterScenarios.Q);
fn = fieldnames(filterScenarios);
nn0 = length(n0_targets);
nzw = length(allZw);

% dat_vecs = zeroes([nf nn0]);

% [ Zw x filterScenarios x n0 targets]
% -> filterScenarios = conduit lengths x MERs = 6
n0_vals  = zeros([nzw nn0 nf]);
n0e_idx   = n0_vals;

rows_to_get = false(size(T,1),1);
masks = zeros([size(T,1), nf]);

% First apply flag the grab n0_total
% -> this SHOULD get me:
%     3 (Z_tot) * 7 (Zw) * 11 (dP) * 2 * (Q) * num(n0_targets) 
%          = 462 * num(n0_targets) values
% -> because the actual n0 target values (e.g. 0.06) may not be represented
%    in all runs, the actual number will be smaller
% -> now can grab these when further subsetting by the filterScenarios
T.is_n0_target = false(size(T,1),1);
count = 0;
for qi = 1:length(all_Q)
    for zti = 1:length(all_Zt)
        for pi = 1:length(all_dP)
            for zwi = 1:nzw
                mask = (T.dP==all_dP(pi)) & (T.Zw==allZw(zwi)) & (T.Q==all_Q(qi)) & (T.Z_total==all_Zt(zti));
                for ni = 1:length(n0_targets)
                    count = count + 1;
                    [n0_val, n0e_idx] = closest(n0_targets, T{mask, 'n0_total'});
                    
                    % Verify that something was retrieved
                    assert(any(mask & ismember(T.n0_total, n0_val)))
                    
                    % Verify the right shape for retrieve values and
                    % targets
                    assert(all(size(n0_val - n0_targets')==[length(n0_targets),1]))
                    
                    % Filter out values that are beyond tolerance
                    n0_val = n0_val(abs(n0_val - n0_targets') <= n0_tol);
                    
                    % Flag remaining
                    T.is_n0_target = T.is_n0_target | (mask & ismember(T.n0_total, n0_val));
                end
            end
        end
    end
end

Faafo here - need to change the filter below to use the new is_n0_target flag above
for ii = 1:nf
    mask = true(size(T,1),1);
    
    % Get mask matching all in filterScenario
    for fi = 1:length(fn)
        this_var = (T.(fn{fi}) == filterScenarios.(fn{fi})(ii));
        assert(sum(this_var)>0, sprintf('None found: %2 = %f', fn{fi}, filterScenarios.(fn{fi})(ii)))
        mask = mask & this_var;
    end
    
    % For each zw, get the n0e value closest to each n0 target
    n0_masks_ii = zeros([size(T,1) nzw]);
    for zwi = 1:nzw

        water_mask = mask & (T.Zw == allZw(zwi));
        n0t_subset = T(water_mask,:).n0_total;
        for ni = 1:nn0
            [n0_vals(zwi,ni,ii), n0e_idx(zwi,ni,ii)] = closest(n0_targets(ni), n0t_subset);
%             if abs(n0_vals(ii,ni) - n0_targets(ni)) > n0_tol
%                 n0_vals(ii,ni) = NaN;
%             end

        end        
        
%         [n0_vals(si,zwi),n0e_idx(si,zwi)] = closest(n0_target,n0t_subset);
        % --> Now use the n0 index to get which n0_excess value to
        % pull. This has to be sample by sample, so initialized values
        % above
        n0_masks_ii(:,zwi) = water_mask & ismember(T.n0_total, n0_vals(zwi,:,ii));
%         this_row = subsetIdx; % & T.n0_excess == all_n0e(n0e_idx(si,zwi));
%         rows_to_get = rows_to_get | this_row;
    end
    n0_mask = any(n0_masks_ii, 2);
    
    mask = mask & n0_mask;
%     n0t_subset = T(mask,:).n0_total;
%     for ni = 1:nn0
%         [n0_vals(ii,ni), n0e_idx(ii,ni)] = closest(n0_targets(ni), n0t_subset);
% %         if abs(n0_vals(ii,ni) - n0_targets(ni)) > n0_tol
% %             n0_vals(ii,ni) = NaN;
% %         end
%         
%     end
%     mask = mask & ismember(T.n0_total, n0_vals(ii,:));
    masks(:,ii) = mask;

%     q = filterScenarios.Q(ii);
%     dP = filterScenarios.dP(ii);
%     Zt = filterScenarios.Z_total(ii);
end
all_mask = any(masks, 2);
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
    
    
    subsetIdx = (T.Z_total == Ztot(ai)) & all_mask;
    [Tsub, sortIdx] = sortrows(T(subsetIdx,:),'Q','descend');
    subMasks = masks(subsetIdx, :);
    subMasks = subMasks(sortIdx, :);
    
    
    % Size by W
    msz = zeros(size(Tsub,1),1);
    msz(Tsub.Q==1e8) = Qsz(1);
    msz(Tsub.Q==1e9) = Qsz(2);
    
    % Secondary filter to get just a few key model results
    unique_lines = unique(Tsub.(common_lines_var));
    
    for ci = 1:length(sCodes)
        codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));

        [x, xlab] = getVar(Tsub, xvar, codeSub);
        [y, ylab] = getVar(Tsub, yvar, codeSub);
        [c, clab] = getVar(Tsub, xvar, mask(subsetIdx));
        legax = 1;
        scatter(x, y, msz)
    end
    
    % Do 2nd nucleation events 2nd
    for ci = 1:length(sCodes)
        codeSub = (Tsub.simpleCode==sCodes(ci) & (Tsub.N_nucl==2 & Tsub.simpleCode>=1));
        scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
            'filled',codeMarkers{ci},...
            'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',codeAlpha(ci),...
            'MarkerFaceAlpha',Nuc_alpha(2))
    end
    
    % --- Psat vs n0 ---- MONEY - tight filter for just a few points
% n_total = (T.n0_excess(codeSub)+1).*T.Cm_0(codeSub);
% subsetIdx = (T.Z_total ==Ztot(ai) & ...
%             and(n_total>=0.55,n_total<=0.065) & ...
%             T.dP );
% 
% yname = 'Total gass mass fraction';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_n_tot_vs_Psat_filtered';
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


    if ai==1
        ylabel(yname)
    end
    title(sprintf('Conduit depth (km b.s.l.) = %.1f',Ztot(ai)))
    xlabel(xname)
    set(gca,'FontSize',fs)
end

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

% yname = 'Total gass mass fraction';
% xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
% oname = 'sweepScatter_n_tot_vs_Psat_filtered';
% xv = 'Psat';
% yv = 'n0_excess';
% 
% legax = 1;
%     for li = 1:length(unique_lines)
%         tsub_sub_mask = (Tsub.(common_lines_var) == unique_lines(li));
%         xxl = Tsub{tsub_sub_mask, xv}.*Pscale;
%         yyl = (Tsub.n0_excess(tsub_sub_mask)+1).*Tsub.Cm_0(tsub_sub_mask);
%         plot(xxl, yyl, 'k')
%         hold on
%     end
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
%         xx = Tsub{codeSub, xv}.*Pscale;
%         yy = (Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub);
%         
%         scatter(xx, yy, msz(codeSub),Tsub.Zw(codeSub),...
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
%         [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
%         scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarker,...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',ca,...
%             'MarkerFaceAlpha',Nuc_alpha(2))
%     end

%     if ai==1
%         ylabel(yname)
%     end
%     title(sprintf('Conduit depth (km b.s.l.) = %.1f',Ztot(ai)))
%     xlabel(xname)
%     set(gca,'FontSize',fs)
% end
%     for ci = 1:length(sCodes)
%         codeSub = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
%         [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
%         scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
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
%         [codeMarker,~, ca] = getCodeMarker(sCodes(ci));
%         scatter(Tsub.Psat(codeSub).*Pscale,(Tsub.n0_excess(codeSub)+1).*Tsub.Cm_0(codeSub),msz(codeSub),Tsub.Zw(codeSub),...
%             'filled',codeMarker,...
%             'MarkerEdgeColor','k',...
%             'MarkerEdgeAlpha',ca,...
%             'MarkerFaceAlpha',Nuc_alpha(2))
%     end
% 
%     if ai==1
%         ylabel(yname)
%     end
%     title(sprintf('Conduit depth (km b.s.l.) = %.1f',Ztot(ai)))
%     xlabel(xname)
%     set(gca,'FontSize',fs)
% end
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

function [x, xname] = getVar(Tsub, thisVar, mask)

    if nargin<3
        mask = true(size(Tsub,1),1);
    end
    switch thisVar
        case 'Psat'  %, 'dP_dT_max', 'dP_dT_bar'}
            x = Tsub{mask, thisVar}.*Pscale;
            xname = {'Supersaturation Pressure','P_{H_2O} - P_f (MPa)'};
        case 'dP_dt_max'
            x = Tsub{mask, thisVar}.*Pscale;
            xname = 'Max. Decompression Rate (MPa s^{-1})';
        case 'dP_dt_bar'
            x = Tsub{mask, thisVar}.*Pscale;
            xname = 'Mean Decompression Rate (MPa s^{-1})';
        case 'phi_f'
            xname = 'Vesicularity, $\phi$';
        case 'Zw'
            xname = 'Water depth (m)';
        case 'Zf'
            xname = 'Fragmentation Depth (m)';
%             x = Tsub{mask, thisVar}; %./1e3;
        case 'Z_nucl'
            xname = 'Final nucleation depth (m)';
            %     mask = (Tsub.simpleCode==sCodes(ci) & ~(Tsub.N_nucl==2 & Tsub.simpleCode>=1));
            x = cellfun(@min,Tsub{mask, thisVar});
            
        case 'n0_excess'
            xname = 'Excess gas mass fraction, $n_{0e}$';
        case 'n0_total'
            xname = 'Total gas mass fraction, $n_{0t}$';

%         otherwise
%             xScale = 1;
            
    end
    if ~exist('x','var')
        x = Tsub{mask, thisVar};
    end
end
