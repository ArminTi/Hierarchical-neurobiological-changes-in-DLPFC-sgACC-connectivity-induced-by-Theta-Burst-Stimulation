function ICA_data = tbs_rseeg_icadata(interp)
% Takes input and creates concatenated file

cfg = [];

if length(interp) == 4

    ICA_data    = ft_appenddata(cfg, interp(1),...
        interp(2),...
        interp(3),...
        interp(4));
    
elseif length(interp) == 5
    
    ICA_data    = ft_appenddata(cfg, interp(1),...
        interp(2),...
        interp(3),...
        interp(4),...
        interp(5));
    
elseif length(interp) == 3
    
    ICA_data    = ft_appenddata(cfg, interp(1),...
        interp(2),...
        interp(3));
end

% Create event list to aid in plotting ICA components
for a = 1:length(interp)
    temp(a).trials = a*ones(1,length(interp(a).trial));
end

ICA_data.cfg.hdr.trialindex = horzcat(temp.trials);
ICA_data.cfg.hdr.triallist  = {interp.block};
ICA_data.label = interp(1).label;
end
