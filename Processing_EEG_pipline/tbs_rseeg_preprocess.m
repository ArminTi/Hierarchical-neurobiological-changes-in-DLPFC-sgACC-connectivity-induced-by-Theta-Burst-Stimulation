function data = tbs_rseeg_preprocess(data)

% Parameters
LINE_NOISE  = 50;

% Hide non-EEG channels from the preprocessing step
EEG_chans   = strcmp(data.hdr.chantype,'EEG'); % ICA pre_processing
EEG_chans   = ~strcmp(data.label,'Marker'); % Gratton pre_processing
extra_dat   = data.trial{1}(~EEG_chans,:);
extra_lab   = data.label(~EEG_chans);

% Remove line noise
[b,a] = butter(2, ([LINE_NOISE*0.98 LINE_NOISE*1.02])/(data.fsample * 0.5),'stop');
data.trial{1}(EEG_chans,:) = filtfilt(b,a,data.trial{1}(EEG_chans,:)')';

% Filtering
cfg             = [];
cfg.bpfilter    = 'yes';
cfg.bpfiltord   = 2;
cfg.bpfreq      = [1 50];
cfg.bsfilter    = 'yes';
cfg.bsfiltord   = 2;
cfg.bsfreq      = [LINE_NOISE*0.98 LINE_NOISE*1.02]; % Second line-noise removal
cfg.demean      = 'yes';
cfg.detrend     = 'yes';
cfg.channel     = {'all','-Marker'}'; %{'all', '-ECG','-HEOG', '-VEOG','-Marker'}';
data            = ft_preprocessing(cfg,data);

% Add extra channels back to data
data.trial{1}(~EEG_chans,:) = extra_dat;
data.label(~EEG_chans)      = extra_lab;

% Downsample 

cfg = [];
cfg.resamplefs = 256;
cfg.method = 'downsample';
data = ft_resampledata(cfg, data);
targetFs = 256;
data.hdr.Fs = targetFs;
data.hdr.nSamples = size(data.time{1}, 2);



