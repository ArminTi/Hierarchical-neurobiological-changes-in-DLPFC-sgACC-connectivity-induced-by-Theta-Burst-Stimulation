function result = tbs_rseeg_trialrejection(input)
% Removes trials that have a maxzvalue > THRESHOLD
% Rejection Criteria/Threshold
    % Range     > 200
    % Var       > 3000
    % Max z     > 12
    % Range z   > 3 
    
THRESHOLD   = 12; % maxzscore allowed = 12
ZTHRESH     = 3;

count = 0;

for row = 1:length(input)
    data = input(row);
    
    if ~isempty(data.trial)
        
        index = [];
        count = count + 1;
        
        % Identify key data features
        EEG         = strcmp(data.hdr.chantype,'EEG');
        SAMPLES     = length(data.time{1});
        TRIALS      = length(data.trial);
        
        % Range-based bad trial selection process
        for t = 1:length(data.trial)
            % Find trials with largest spread across all EEG channels
            trial_range = range(data.trial{t}(EEG,:));
            temp(t)     = sum(trial_range);
        end
        temp            = zscore(temp);
        index_bad_trls  = temp > ZTHRESH;

        % Generate concatenated data set for calculations
        catdat           = horzcat(data.trial{:});

        % Remove non-EEG channels from zscore calculations
        %catdat(~EEG,:) = [];
        CHANNELS = size(catdat,1);

        % Convert data into zscores
        catdat = zscore(catdat,0,2);

        % Put data back into fieldtrip cell format
        trials = mat2cell(catdat,CHANNELS,repmat(SAMPLES,TRIALS,1));

        % Find the max zscore for channels and trials
        maxzscore = zeros(CHANNELS,TRIALS);
        for a = 1:TRIALS
            for b = 1:CHANNELS
                maxzscore(b,a) = max(trials{a}(b,:));
            end
        end

        % Find max zscore for trials
        maxtrials = max(maxzscore,[],1);

        % Find trials with max zscores above threshold
        index_max = maxtrials > THRESHOLD;
        
        index = logical(index_max + index_bad_trls);

        % Remove trials from data
        data.trial(index)        = [];
        data.time(index)         = [];
        data.event(index)        = [];
        data.sampleinfo(index,:) = [];
        data.trialinfo.rejtrls   = index;
        result(count)            = data;
        
    end
end
end