
function [data] = tbs_rseeg_resample(data)

% Initialize configuration for resampling
cfg = [];
cfg.resamplefs = 256;
cfg.method = 'downsample';
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 128;
data = ft_resampledata(cfg, data);
targetFs = 256;
data.hdr.Fs = targetFs;
data.hdr.nSamples = size(data.time{1}, 2);
end

