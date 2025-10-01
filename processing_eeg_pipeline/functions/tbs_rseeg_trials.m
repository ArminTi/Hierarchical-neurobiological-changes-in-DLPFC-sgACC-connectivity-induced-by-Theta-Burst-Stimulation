function trials = tbs_rseeg_trials(events)
% Segments the data into trials

% Find all 'RS-EEG' events
EEG = find(strcmp({events.event.type}, 'RS-EEG'));
EEG_LABELS = {'RSpre1', 'RSpre2', 'RSpost1', 'RSpost2', 'RSpost3'};
chs = 1:length(events.label);

% Define segment length (2 seconds)
segment_length = 2;

for s = 1:length(EEG)
    tempdata = events;
    EEG_start = tempdata.event(EEG(s)).sample;
    EEG_end = tempdata.event(EEG(s)).sample + tempdata.event(EEG(s)).sampdur - 1;
    index = EEG_start:EEG_end;
   
    % Calculate the number of full 2-second segments
    total_secs = floor(length(index) / (tempdata.fsample * segment_length));
    remainder = length(index) - total_secs * tempdata.fsample * segment_length;
   
    % Ignore first <2s remainder of data and keep the rest
    A = tempdata.trial{1}(chs, index((remainder + 1):end));
    B = mat2cell(A, length(chs), repmat(tempdata.fsample * segment_length, total_secs, 1));
    tempdata.trial = B;
   
    % Extract corresponding time data and segment into 2-second chunks
    C = tempdata.time{1}(1, index((remainder + 1):end));
    D = mat2cell(C, 1, repmat(tempdata.fsample * segment_length, total_secs, 1));
    tempdata.time = D;
   
    % Create event structure for each 2-second segment
    for t = 1:length(tempdata.trial)
        tempevent(t).type = 'eeg';
        tempevent(t).value = 254;
        tempevent(t).sample = index(remainder + 1) + ((t - 1) * tempdata.fsample * segment_length);
        tempevent(t).timestamp = tempdata.time{t}(1, 1);
        tempevent(t).sampdur = size(tempdata.trial{t}, 2);
        tempevent(t).duration = tempevent(t).sampdur / tempdata.fsample;
    end
   
    % Update sampleinfo for each segment
    tempdata.sampleinfo = [];
    tempdata.sampleinfo(:, 1) = index(remainder + 1):tempdata.fsample * segment_length:EEG_end;
    tempdata.sampleinfo(:, 2) = [(tempdata.sampleinfo(2:end, 1) - 1); EEG_end];
   
    % Assign the events and block label to the segmented data
    tempdata.event = tempevent;
    tempdata.block = EEG_LABELS{s};
   
    % Store the segmented data in the output structure
    trials(s) = tempdata;    
end
end
