function tbs_rseeg_plotica2(comps, suggested_comps)

%% Parameters
FREQ            = [0.5 70];
NDISP           = 5; % Increase if you want more components per page
SIG_THRESHOLD   = 2.5;
NCOMP           = 25;
WIDTH           = 4;
BACKGROUND      = [0.9 0.9 0.9];
SPREAD          = 0.95;
TOP             = 1.05;
SUB_HEIGHT      = 0.1;
SUB_WIDTH       = 0.2;

%% Screen parameters
figure
pos = get(groot, 'Screensize');
pos(3) = pos(4) / NDISP * WIDTH;
set(gcf, 'position', pos)

%% PSD Configuration
cfg             = [];
cfg.output      = 'pow';
cfg.channel     = 1:64;
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';
cfg.foilim      = FREQ;

psd_hann = ft_freqanalysis(cfg, comps);

%% Topology z-scores 
ztopology   = zscore(comps.topo, 0, 1);

%% RS traces
traces = cat(3, comps.trial{:});

%% N-back plotting data
for t = 1:length(comps.time)
    temp_time{t} = comps.cfg.previous.hdr.trialindex(t) * ones(1, length(comps.time{t}));
end

cnttime     = horzcat(temp_time{:});
cntdata     = horzcat(comps.trial{:});
plottime    = 1:length(cnttime);

%% Plot ICA Components
for n = 1:ceil(NCOMP / NDISP)
    for m = 1:NDISP
        c = (n - 1) * NDISP + m;
        if c > NCOMP
            break
        end

        % Time-series plot
        subplot('position', [0.07, TOP - m / NDISP * SPREAD, SUB_WIDTH, SUB_HEIGHT])
        
        for a = unique(cnttime)
            trl_ind = cnttime == a;
            plot(plottime(trl_ind), cntdata(c, trl_ind), 'linewidth', 1)
            hold on
            axis tight
        end

        if ismember(c, suggested_comps)
            set(gca, 'color', BACKGROUND)
        end

        box off
        ylabel(num2str(c), 'fontsize', 10);

        if m == 1
            title('Time-series', 'fontsize', 10);
        end

        if m == NDISP
            xlabel('time [s]')
            legend(comps.cfg.previous.hdr.triallist, 'fontsize', 6, ...
                'location', 'none', 'position', [0.1, 0.77, 0.1, 0.07])
        end

        % PSD plot
        subplot('position', [0.31, TOP - m / NDISP * SPREAD, SUB_WIDTH, SUB_HEIGHT])
        semilogy(repmat(psd_hann.freq, NCOMP, 1)', psd_hann.powspctrm(1:NCOMP, :)', ...
            'linewidth', 0.05, 'color', [0.5 0.5 0.5])
        hold on
        semilogy(psd_hann.freq, psd_hann.powspctrm(c, :), 'linewidth', 2.5)
        axis tight
        box off
        ylim([min(psd_hann.powspctrm(:)) max(psd_hann.powspctrm(:))])
        xlim(FREQ)
        set(gca, 'fontsize', 8)

        if m == 1
            title('PSDs', 'fontsize', 10);
        end

        if m == NDISP
            xlabel('freq [Hz]', 'fontsize', 9)
        end

        if ismember(c, suggested_comps)
            set(gca, 'color', BACKGROUND)
        end

        % Topoplot
        subplot('position', [0.54, (TOP - 0.02) - m / NDISP * SPREAD, 0.20, 0.113])
        cfg             = [];
        cfg.layout      = 'easycapM1.mat';
        cfg.component   = c;
        cfg.comment     = 'no';
        cfg.feedback    = 'no';
        cfg.colorbar    = 'yes';

        sigchannels             = find(abs(ztopology(:, c)) > SIG_THRESHOLD);
        cfg.highlight           = 'on';
        cfg.highlightchannel    = sigchannels;
        cfg.highlightsymbol     = '*';
        cfg.highlightcolor      = 'red';
        cfg.highlightsize       = 6;
        cfg.highlightfontsize   = 6;

        ft_topoplotIC(cfg, comps);
        set(gca, 'fontsize', 8)

        if m == 1
            title('Topoplot', 'fontsize', 10);
        end

        % N-back trials plot
        subplot('position', [0.77, TOP - m / NDISP * SPREAD, SUB_WIDTH, SUB_HEIGHT])
        traces_comp = squeeze(traces(c, :, :));
        imagesc(traces_comp')
        box off

        if m == 1
            title('N-back Trials', 'fontsize', 10);
        end
    end

    % Control Button for moving to the next set of components
    drawnow
    h = uicontrol('Position', [pos(3)/2-80, 10, 150, 30], 'String', 'Next', ...
        'Callback', 'uiresume(gcbf)', 'FontWeight', 'bold', 'FontSize', 12, ...
        'ForegroundColor', [0 0 1], 'BackgroundColor', [1 1 1]);
    uiwait(gcf);
end
