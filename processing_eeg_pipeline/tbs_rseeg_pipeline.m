

%% Add paths to toolboxes/dependencies
restoredefaultpath
addpath C:/Users/ASUS/Desktop/Apps/fieldtrip-20231220
addpath C:/Users/ASUS/Desktop/Apps/auto-ica-checking
addpath C:/Users/ASUS/Desktop/Apps/spm12


ft_defaults

%% Loop through subjets and process data one at a time 
% Modify directory below to where you've saved the open tbs data
DIRECTORY = 'D:/Research Topic/DR nasrabadi/2. dataset';

for subject = 7
    for session = 1        
        % generate filename
        filename = fullfile(DIRECTORY,...
            ['sub-',sprintf('%02d',subject),'_',...
            'ses-',sprintf('%02d',session),'_',...
            'sourcedata.mat']);
        
        % load source data
        load(filename);

        % line noise removal and filtering  
        %preproc_resam = tbs_rseeg_preprocess(data); %its already done if you download
        %data babak sent to you
        
        % generate resting state events
        events  = tbs_rseeg_events(preproc_resam);

        % segment the data into trials
        trials  = tbs_rseeg_trials(events);
        
        % visual inspection and trial/channel rejection
        reject  = tbs_rseeg_cleaning(trials);
        
        % interpolate removed channels
        interp  = tbs_rseeg_interp(reject,trials);
        
        % ica cleaning
        comps   = tbs_rseeg_ica(interp);
        
        % inspect ica components for rejection - addpath for icacheck
        icacheck   = tbs_rseeg_inspectica(comps);

        % reject ica components
        postica = tbs_rseeg_rejectica(interp,icacheck);
        
        % rereferencing to common average and final checking
        reref   = tbs_rseeg_reref(postica);

        %save reref according to subject, session and TBS name
        % Save 'reref' with a filename based on subject and session
        save_filename = fullfile(DIRECTORY, ...
            ['sub-', sprintf('%02d', subject), '_', ...
            'ses-', sprintf('%02d', session), '_', ...
            'preprocessed.mat']);
        
        save(save_filename, 'reref');

        % convert to spm data
        spm_data = tbs_rseeg_spmload(subject, session, reref);
    end
end




