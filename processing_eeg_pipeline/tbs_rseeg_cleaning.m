function allrej = tbs_rseeg_cleaning(input)

% channel type and unit correction
for i = 1:length(input)
    input(i).hdr.chantype(63) = {'EEG'};
    input(i).hdr.chantype(64) = {'EEG'};
    input(i).hdr.chantype(65) = {'ECG'};
    input(i).hdr.chantype(66) = {'EOG'};
    input(i).hdr.chantype(67) = {'EOG'};
    input(i).hdr.chantype(68) = {'Marker'};
end
for i = 1:length(input)
    input(i).hdr.chanunit(67) = {'uV'};
    input(i).hdr.chanunit(68) = {'V'};
end
fprintf('\nAutomatic\n\n')
% Automatic rejection of trials    
for row = 1:length(input)
    autorej(row) = tbs_rseeg_trialrejection(input(row));
    
end

% instructions
fprintf('\nINSPECT DATA AND NOTE ITEMS FOR REJECTION IN SEPARATE SPREADSHEET\n\n')

% inspect all trials (without rejecting yet) - this is doe to get an
% overall feel for the quality of the data
for a = 1:length(autorej)
    
    % Quick visual inspection of all trials - identify possible BAD trials
    cfg             = [];
    cfg.channel     = 1:64;
    cfg.viewmode    = 'butterfly'; % butterfly/vertical
    cfg.alim        = [-100 100];
    cfg.blocksize   = round(autorej(a).time{1}(end) - autorej(a).time{1}(1));
    cfg.continuous  = 'no';
    cfg.colorgroups = 'sequential'; 
    cfg.layout      = 'easycapM1.mat';
    
    % Note down potential trials for rejection in separate spreadsheet
    ft_databrowser(cfg,autorej(a))
end

% Select and remove BAD trials and channels using visual inspection 
for a = 1:length(autorej)
    cfg                 = [];
    cfg.method          = 'trial';
    cfg.alim            = 50;
    cfg.eogscale        = 0.1;
    cfg.ecgscale        = 0.05;
    cfg.preproc.detrend = 'yes';
    cfg.preproc.demean  = 'yes';
    visrej(a)   = ft_rejectvisual(cfg,autorej(a));
end

% Select and remove BAD trials and channels using summary statistics
for a = 1:length(visrej)
    cfg         = [];
    cfg.channel	= {'all','-ECG','-EOG_VERT','-EOG_LAT','-Marker'};
    cfg.method  = 'summary';
    cfg.metric  = 'var'; % or: 'zvalue';
    allrej(a)   = ft_rejectvisual(cfg,visrej(a));
end

end