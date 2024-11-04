%%%%%%% scripts for pre-processing physionet EEG data (see manuscript for the data source)
%%% eeglab 13.3.2b was used for these purposes.
%%%% note that IC(s were manually removed from the data after running this script. 
%% clear
clc;
clear all
close all
restoredefaultpath;
addpath c:\__matlab_lic_\fieldtrip-20231220\
addpath c:\__matlab_lic_\eeglab2024.0\
%% set directory
eeglab
%hamed% cd /home/frederik/data/EEG_rest/raw/
path_01 = 'c:\_WORKS_\Moffa_data_VanDeSteen_codes\data\EEG_rest\raw\';
path_02 = 'c:\_WORKS_\Moffa_data_VanDeSteen_codes\data\EEG_rest\prep\';

cd(path_01);

%% pre-process both EO and EC data
fold = dir(fullfile([pwd '/sub*.set']))
tel = 0;
%%
for i = 1:length(fold)
    try
       
        %% load data
        
        EEG = pop_fileio([path_01,fold(i).name]);

        %% get channel locations
       
        %hamed% for 2022_Moffa Data
        DatasetChanlocs = ["FT9","F7","FT7","T7","AF7","F5","FC5","C5","FP1", ...
        "AF3","F3","FC3","F1","FC1","Fz","FCz","TP7","P7","PO9","CP5","P5", ...
        "PO7","C3","CP3","O1","Oz","P3","PO3","POz","P1","Pz","CP1","FT10", ...
        "F8","FT8","T8","AF8","F6","FC6","C6","FP2","AF4","F4","FC4","C4", ...
        "F2","FC2","Cz","TP8","P8","PO10","CP6","P6","PO8","O2","PO4","CP4", ...
        "P4","C2","CP2","P2","CPz","C1","AFz"];
        channel_loc_data = readtable([path_01 fold(i).name(1:end-4) '.csv']);
        % Extract the relevant columns
        channel_labels = channel_loc_data.name(1:65);
        channel_x = channel_loc_data.x(1:65);
        channel_y = channel_loc_data.y(1:65);
        channel_z = channel_loc_data.z(1:65);
        for j = 1:64
            EEG.chanlocs(j).labels =DatasetChanlocs{j};
            index = find(strcmp(channel_labels, EEG.chanlocs(j).labels));
            EEG.chanlocs(j).X = channel_x(index);
            EEG.chanlocs(j).Y = channel_y(index);
            EEG.chanlocs(j).Z = channel_z(index);
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',fold(i).name(1:end-4),'gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_raw.set'],'filepath',path_02);
        %% Selecting 238 s from ~ 240 s
        
        EEG = pop_select( EEG,'time',[0 238] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
        EEG = eeg_checkset( EEG );
         
        %% Downsampling to 256 Hz
        
        EEG = pop_resample(EEG,256);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp.set'],'filepath',path_02);

        %% Filter: highpass (1 Hz) & lowpass (60 Hz)

        EEG = pop_eegfiltnew(EEG, 1.0, 60,[], false, [], 0);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed.set'],'filepath',path_02);

        %% Using Zapline to remove line noise
        
        EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('noisefreqs','line','plotResults',false))
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln.set'],'filepath',path_02);

        %% Clean rawdata and remove bad channels

        % Make a copy of the original channel locations (with 64 channels)
        originalChanlocs = EEG.chanlocs;

        [EEG,~,~,removed_channels] = clean_artifacts(EEG, ...
            'FlatlineCriterion', 5, ...          % Identifies flat channels
            'ChannelCriterion', 0.8, ...         % Correlation threshold for bad channels
            'RemoveFlatChannels', 1, ...         % remove flat channels
            'Remove', 1, ...                     % remove bad channels
            'LineNoiseCriterion', 4, ...         % Threshold for line noise
            'Highpass', 'off', ...               % Keep high-pass filter off
            'BurstCriterion', 20, ...            % Artifact burst criterion
            'WindowCriterion', 0.25, ...         % Identifies bad windows
            'BurstRejection', 'off', ...         % Prevent rejection of bursts
            'CleanArtifacts', 0, ...             % Prevent removal of any flagged artifacts
            'Distance', 'Euclidian', ...
            'WindowCriterionTolerances', [-Inf +Inf]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv.set'],'filepath',path_02);

        %% Interpolating removed channels
       
        % Interpolate missing channels based on original channel locations
        EEG = pop_interp(EEG, originalChanlocs, 'spherical');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv_interp.set'],'filepath',path_02);

        %% Re-Referencing All Channels to common-average

        % Re-reference to common-average
        EEG = pop_reref(EEG, []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv_interp_avrgd.set'],'filepath',path_02);

        %% running ICA
        
        EEG = pop_runica(EEG, 'extended',1,'interupt','on');
        EEG = eeg_checkset(EEG, 'ica');  % Ensures the EEG structure has ICA fields calculated
        EEG = pop_iclabel(EEG, 'default');  % Classify ICA components with ICLabel
        % Extract the classification matrix from ICLabel
        classifications = EEG.etc.ic_classification.ICLabel.classifications;
        % The first column corresponds to the probability of being "brain"
        brain_probabilities = classifications(:, 1);
        % Find components where the probability of being "brain" is greater than 70%
        brain_components = find(brain_probabilities > 0.7);
        
        % Load standard head model and set up DIPFIT model
        EEG = pop_dipfit_settings(EEG, ...
            'hdmfile', 'c:\__matlab_lic_\eeglab2024.0\plugins\dipfit\standard_BEM\standard_vol.mat', ...  % Standard BEM volume conductor model
            'coordformat', 'MNI', ...
            'mrifile', 'c:\__matlab_lic_\eeglab2024.0\plugins\dipfit\standard_BEM\standard_mri.mat', ...  % Standard MRI
            'chanfile', 'C:\__matlab_lic_\eeglab2024.0\plugins\dipfit\standard_BEM\elec\standard_1005.elc', ...
            'coord_transform', [1.04721     -15.3702     -67.1375    -0.327775    0.0218727     -1.57753      1.12798      1.12798      1.12798], ...  % Transformation (example values)
            'chansel', 1:EEG.nbchan);  % Use all channels for the fit
        % Fit dipoles to all ICA components
        EEG = pop_dipfit_batch(EEG, 1:size(EEG.icaweights, 1));
        % Extract the residual variance values
        rv_values = [EEG.dipfit.model.rv];  % Assuming dipfit is used and initialized
        % Find the indices of components with RV < 0.15
        low_rv_components = find(rv_values < 0.15);
        % Display the indices of components that meet the threshold
        disp('Components with >70% probability of being brain:');
        disp(brain_components');
        disp('Components with RV < 0.15:');
        disp(low_rv_components);
        acceptable_components = union(brain_components,low_rv_components);
        disp('Acceptable Components:');
        disp(acceptable_components');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv_interp_avrgd_ica.set'],'filepath',path_02);
        
        %% Removing ICA components
        
        components_to_remove = setxor([1:size(EEG.icaweights, 1)], acceptable_components);
        EEG = pop_subcomp(EEG, components_to_remove, 0);  % Remove specified components and reconstruct EEG
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui', 'off');
        EEG = eeg_checkset(EEG);
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv_interp_avrgd_ica_cleaned.set'],'filepath',path_02);

        %% epoching

        epoch_length = 2;  % 2 seconds
        samples_per_epoch = epoch_length * EEG.srate;  % Number of samples in each epoch
        % Initialize the event structure
        EEG.event = [];
        % Get the total number of samples in the dataset
        total_samples = EEG.pnts;
        % Create events at every 2 seconds (samples_per_epoch)
        num_epochs = floor(total_samples / samples_per_epoch);  % Total number of full epochs
        for iii = 1:num_epochs
            EEG.event(iii).type = 'epoch';  % Label the event as 'epoch'
            EEG.event(iii).latency = (iii - 1) * samples_per_epoch + 1;  % Latency in samples
        end
        EEG = pop_epoch(EEG, {'epoch'}, [0 epoch_length]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        %EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_resamp_h_l_passed_zapln_badchrmv_interp_avrgd_ica_cleaned_epoched.set'],'filepath',path_02);

        %% STFT
        
        % Define parameters for STFT
        % Frequency range in Hz
        fmin = 1;
        fmax = 50;
        epoch_length = 2;  % in seconds
        window_length = EEG.srate * epoch_length;  % Window length for each epoch (no overlap)
        % Initialize cell array to store STFT results
        stft_results = cell(EEG.nbchan, EEG.trials);
        % Initialize array to store average power values for each epoch
        epoch_powers = zeros(1, EEG.trials);
        % Loop through each channel and each epoch to compute the STFT
        for epoch_idx = 1:EEG.trials
            % Initialize power accumulator for the current epoch
            total_power = 0;
            for chan_idx = 1:EEG.nbchan
                % Extract the data for this epoch and channel
                epoch_data = EEG.data(chan_idx, :, epoch_idx);
                [s, f] = stft(epoch_data, EEG.srate,...
                    'Window', hann(window_length),...
                    'OverlapLength', 0, ...
                    'FrequencyRange', 'onesided');
                % Limit frequencies to range of interest (1-50 Hz)
                freq_indices = find(f >= fmin & f <= fmax);
                power = abs(s(freq_indices, :)).^2;
                % Compute mean power for the frequency range across time and frequency bins
                mean_power = mean(power(:));
                % Accumulate mean power for all channels in this epoch
                total_power = total_power + mean_power;
            end
            % Calculate average power across channels for this epoch
            epoch_powers(epoch_idx) = total_power / EEG.nbchan;
        end
        % Calculate Z-scores of the power values across epochs
        z_scores = (epoch_powers - mean(epoch_powers)) / std(epoch_powers);
        % Identify epochs where Z-score exceeds ±3
        bad_epochs = find(abs(z_scores) > 3);
        disp('bad_epochs:');
        disp(bad_epochs);
        % Remove bad epochs from the EEG dataset
        EEG = pop_rejepoch(EEG, bad_epochs, 0);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename', [fold(i).name(1:end-4) '_final.set'],'filepath',path_02);


        %%
        STUDY = []; 
        CURRENTSTUDY = 0; 
        ALLEEG = []; 
        EEG=[]; 
        CURRENTSET=[];
        eeglab redraw
    catch
        tel = tel+1
        fprintf('an error occured for ')
        error_id{tel,1} = fold(i).name;
        error_id{tel,2} = i;
    end
end
try
    save('/home/frederik/data/EEG_rest/errors_subjects','error_id')
    save('c:\_WORKS_\Moffa_data_VanDeSteen_codes\data\EEG_rest\errors_subjects','error_id')
catch
    disp('no error occured')
end