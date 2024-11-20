function data = tbs_rseeg_events(data)

MARKER      = strcmp(data.label,'Marker');

RS_EEG      = 245;

%% Find trigger values
i               = [];
trigger         = data.trial{1}(MARKER,:);
trigger(trigger == 0) = 1;
trigger(1)      = 100;                      % set initial trigger to value 100
trigger(2)      = 101;
trigger(end)    = 100;                      % set final trigger to end value 100   
dif_trigger     = [0,diff(trigger)~=0];     % is 1 when value changes
trigger         = trigger.*dif_trigger;     % set trigger to 0 if it doesn't change
i               = find(trigger);            % indices of triggers

event_data      = [];
event_data(:,1) = trigger(i);                    % Value of the amplifier (252:255)
event_data(:,2) = i;                             % Sample location of the events
event_data(:,3) = i/data.fsample;                % Time of event in seconds
event_data(:,4) = [diff(event_data(:,2)); 0];    % Duration of events in samples
event_data(:,5) = event_data(:,4)/data.fsample;  % Duration of events in seconds

short_trls      = event_data(:,5) < 1e-03;
short_trls(1)   = 0;
event_data(short_trls,:) = [];

%% Categorise events 
% Alternative
% temp = [event_data(:,1)];
% hits = strfind(temp,[254,252]);
% event_data(hits,6) = 1;
% miss = strfind(temp,[254,255]);
% event_data(miss,6) = 2;
% fa   = strfind(temp,[255,253,252]);
% event_data(fa,6)   = 3;
% nost = strfind(temp,[255,253,255]);
% event_data(nost,6) = 4;

% RESTING STATE
EEG_triggers = find(event_data(:,5) > 230 & event_data(:,5) < 250 & event_data(:,1) == RS_EEG);
event_data(EEG_triggers,6) = 1; % Resting state             


%% Create Events list
TRIGGERS = 1;
TYPE = 'RS-EEG';

event_index = find(event_data(:,6) == 1);
for n = 1:length(event_index)
    [~,j] = intersect(TRIGGERS,event_data(event_index(n),6));
    data.event(n).type      = TYPE;                       % event type (hit etc...)
    data.event(n).value     = event_data(event_index(n),1);  % trigger value
    data.event(n).sample    = event_data(event_index(n),2);  % sample index
    data.event(n).timestamp = event_data(event_index(n),3);  % time of event
    data.event(n).sampdur   = event_data(event_index(n),4);  % duration (samples)
    data.event(n).duration  = event_data(event_index(n),5);  % duration (seconds)
end

end