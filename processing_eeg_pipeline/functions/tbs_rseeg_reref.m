function reref = tbs_rseeg_reref(postica)

% Rereference
cfg             = [];
%cfg.channel     = 'all';
cfg.channel     = 'EEG'; %{'all', '-ECG','-EOG_LAT', '-Marker'}';
cfg.refchannel  = 'all';
cfg.reref       = 'yes';
cfg.refmethod   = 'avg';
cfg.trials      = 'all';

for n = 1:length(postica)
    reref(n) = ft_preprocessing(cfg,postica(n));
end

% Add event and block-name structures to file
for o = 1:length(postica)
    reref(o).event = postica(o).event;
    reref(o).block = postica(o).block;
end

fprintf('\nINSPECT DATA AND NOTE ITEMS FOR REJECTION IN SEPARATE SPREADSHEET\n\n')

% inspect all trials (without rejecting yet) - this is doe to get an
% overall feel for the quality of the final data
for a = 1:length(reref)
    
    % Quick visual inspection of all trials - identify possible BAD trials
    cfg             = [];
    cfg.channel     = 1:64;
    cfg.viewmode    = 'butterfly'; % butterfly/vertical
    cfg.alim        = [-100 100];
    cfg.blocksize   = round(reref(a).time{1}(end) - reref(a).time{1}(1));
    cfg.continuous  = 'no';
    cfg.colorgroups = 'sequential'; 
    cfg.layout      = 'easycapM1.mat';
    
    % Note down potential trials for rejection in separate spreadsheet
    ft_databrowser(cfg,reref(a))
end

end