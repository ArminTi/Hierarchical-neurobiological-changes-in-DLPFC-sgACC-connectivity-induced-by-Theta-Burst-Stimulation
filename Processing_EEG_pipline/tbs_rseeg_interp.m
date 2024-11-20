function interp = tbs_rseeg_interp(reject,trials)

% load necessary files
load('easycapM1_neighb.mat'); % electrode neighbours
elec = ft_read_sens('easycap-M1.txt'); % 3D electrode layout

for a = 1:length(reject)

    % Store list of removed trials and channels
    bad_trial{a} = ismember(trials(a).sampleinfo(:,1), ...
        setdiff(trials(a).sampleinfo(:,1), reject(a).sampleinfo(:,1)));
    bad_chans{a} = setdiff(trials(a).label(1:64),reject(a).label);

    % Interpolate data for any removed channels
    cfg                 =[];
    cfg.method          = 'weighted'; 
    cfg.trials          = 'all';
    cfg.neighbours      = neighbours;
    cfg.elec            = elec;
    cfg.lambda          = 1e-5;
    cfg.order           = 4;
    cfg.missingchannel  = bad_chans{a};
    cfg.senstype        = 'EEG';
    interp(a)           = ft_channelrepair(cfg, reject(a));

    
    interp(a).hdr.rmv_trls = bad_trial{a};
    interp(a).hdr.rmv_chns = bad_chans{a};
    

    % Sort channels to consistent order
    b = zeros(size(trials(a).label(1:64)));
    % Loop through each label
    for i = 1:length(trials(a).label(1:64))
        % Find the index
        b(i) = find(strcmp(interp(a).label, trials(a).label{i}));
    end
    % Reorder

%     % Sort channels to consistent order
%     flipped = reverse(interp(a).label);
%     [~, b]  = sort(flipped);
    
    
    interp(a).label(1:length(b))        = interp(a).label(b);
    interp(a).hdr.label                 = interp(a).label(b);

    % Sort channels in trials
    for t = 1:length(interp(a).trial)
        interp(a).trial{t}(1:length(b),:)   = interp(a).trial{t}(b,:);
    end

    % Add back ECG | EOG | Marker from trials - less bad trials
    interp(a).label(65:68) = trials(a).label(65:68);     
    good_trial = find(~bad_trial{a});

    for trl = 1:length(good_trial)
        interp(a).trial{trl}(65:68,:) = trials(a).trial{good_trial(trl)}(65:68,:);
    end
end

end
