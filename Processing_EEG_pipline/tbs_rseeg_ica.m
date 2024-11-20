function [comps] = tbs_rseeg_ica(interp)

% Concatenate data structure to allow ICA
ICA_data    = tbs_rseeg_icadata(interp);

% ICA
cfg         = [];
cfg.method  = 'runica';
cfg.channel = 1:64;

comps = ft_componentanalysis(cfg,ICA_data);
end