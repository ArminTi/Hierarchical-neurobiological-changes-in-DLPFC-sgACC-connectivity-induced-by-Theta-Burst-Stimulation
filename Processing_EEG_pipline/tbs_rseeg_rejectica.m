function postica = tbs_rseeg_rejectica(interp,comps)

for m = 1:length(interp)
        
    % Remove non-EEG channels
    cfg             = [];
    cfg.channel     = {'all','-ECG','-EOG_VERT','-EOG_LAT','-Marker'};
    interp(m)       = ft_selectdata(cfg,interp(m));

    % Create ind component time series using original data
    cfg             = [];
    cfg.unmixing    = comps.unmixing; % NxN unmixing matrix
    cfg.topolabel   = comps.topolabel; % Nx1 cell-array with the channel labels
    cfg.channel     = 1:64;
    comp_orig       = ft_componentanalysis(cfg,interp(m));

    % Original data reconstructed excluding rejected components
    cfg             = [];
    cfg.component   = comps.rejected;
    cfg.channel     = 1:64;
    postica(m)      = ft_rejectcomponent(cfg,comp_orig,interp(m));
end

% for a = 1:length(postica)
%     
%     % Quick visual inspection of all trials - identify possible BAD trials
%     cfg             = [];
%     cfg.channel     = 1:64;
%     cfg.viewmode    = 'butterfly'; % butterfly/vertical
%     cfg.alim        = [-100 100];
%     cfg.blocksize   = round(postica(a).time{1}(end) - postica(a).time{1}(1));
%     cfg.continuous  = 'no';
%     cfg.colorgroups = 'sequential'; 
%     cfg.layout      = 'easycapM1.mat';
%     
%     % Note down potential trials for rejection in separate spreadsheet
%     ft_databrowser(cfg,postica(a))
% end
end