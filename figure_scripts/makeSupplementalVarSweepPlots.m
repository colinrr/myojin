%% Make supp figs with all results plotted

run config
close all

printFigs = true;
dpi = 500;

vars = {'n0_total', 'dP_dt_max', 'N_nucl'};
% vars = {'N_nucl'};

logQs = [8 9];

% Figure options
figdims = [24 12];
fs = 6;
lfs = 8;
tfs = 8;

n0_cax = [0.01 0.11];

for vi = 1:length(vars)
    vari = vars{vi};
    fprintf('%s:\n', vari)
    switch vari
        case 'dP_dt_max'
            log_plot = true;
            cbl = 'Max. Decompression Rate (MPa/s)';
        case 'n0_total'
            log_plot = false;
            cbl = 'Total Volatile Mass Fraction, $n_{0t}$ (wt.\%)';
        case 'N_nucl'
            log_plot = false;
            cbl = 'Number of discreet nucleation events';
    end
    
    [axes, lh, cbh, figs] = plotVariableAllSweeps(vari, [], log_plot, true);
    

    
    for ii = 1:2
        for ai = 1:numel(axes{ii})
            
            set(axes{ii}(ai), 'FontSize', fs)  
            xl = get(axes{ii}(ai), 'xlabel');
            yl = get(axes{ii}(ai), 'ylabel');
            tl = get(axes{ii}(ai), 'title');
            
            switch vari
                case 'dP_dt_max'
%                     cbh.Label.String = 'Max. Decompression Rate (MPa/s)$';

                case 'n0_total'
                    colormap(axes{ii}(ai), redblue); 
                    caxis(n0_cax); 
%                     cbh.Label.String = 'Total Volatile Mass Fraction, $n_{0t}$ (wt.\%)';

                case 'N_nucl'
%                     cbh.Label.String = 'Number of discreet nucleation events';
                    co = get(0,'DefaultAxesColorOrder');
                    nc = 4;
                    caxis(axes{ii}(ai),[-.5 nc-.5])
                    colormap([0.5 0.5 0.5; co([1 2 3],:)])
                    cbh{ii}.Ticks = [0:nc];
            end
        end
        set(lh{ii}, 'FontSize', lfs)
        set(cbh{ii}, 'FontSize', lfs)
        cbh{ii}.Label.Interpreter = 'latex';
        cbh{ii}.Label.String = cbl;
    end
    
    if printFigs
        for qi = 1:length(logQs)
            % Save both pdf and png for now
            fig = figs{qi};
            figure(fig)
            fprintf('\tQ %i\n', qi)
            disp(' --> pdf...')
            figname = sprintf('SuppFig_%s_AllSweeps_Q%i', vari, logQs(qi));
            printpdf(figname, FIGURES_DIR, figdims)
            
            disp(' --> png...')
            fullpath = fullfile(FIGURES_DIR, [figname, '.png']);
            set(fig,'paperunits', 'centimeters' ,'paperposition',[0 0 figdims],'paperpositionmode','manual')
            set(fig,'papersize', figdims)
            saveas(fig, fullpath, 'png')
%             exportgraphics(fig, fullpath, 'Resolution', dpi);
        end
    end
end