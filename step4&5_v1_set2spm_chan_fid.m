%Electrodes were grounded to Fpz
%%% this script converts preprocessed ICA EEGLAB .set data into SPM12 M/EEG objacts
%%% both pre (last 10secs) and post (first 10secs) data stored in a single M/EEG object
clear all;
close all;
clc;
%run feildtrip (it is recomended to not be directly add to perminant paths)
original_path = path;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/fieldtrip-20240620/fieldtrip-20240620');
ft_defaults;
% starting spm and eeglab
spm('defaults', 'eeg');
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui');

% Define paths
base_path = 'C:/Users/Growth fire/Documents/research/signal/convert/ICA clean'; % Directory for preprocessed ICA EEG.set data 
pre_files = dir(fullfile(base_path, 'sub*_pre_*_ICA_clean.set')); % Directory to pre TBS data
post_files = dir(fullfile(base_path, 'sub*_post_*_ICA_clean.set')); % Directory to post TBS data
output_folder = 'C:/Users/Growth fire/Documents/research/signal/convert/new out';
Elc_path = 'C:/Users/Growth fire/Documents/research/signal/convert/electrodes'; % Directory for CSV files storing EEG electrodes and fiducial data

% Defining the channel labels (if the letters were small or capital they wouldn't be recognized as eeg by SPM)
chan64 = {'AFz','Fp1','Fp2','AF3','AF4','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8',...
          'CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','PO3','PO4','O1','O2','AF7','AF8','F5','F1','F2','F6',...
          'FC3','FCz','FC4','FT8','C5','C1','C2','C6','TP7','FT7','CP3','CPz','CP4','TP8','P5','P1','P2','P6',...
          'PO7','POz','PO8','Oz','P9','P10','TP9','FP1','FZ','PZ','FP2','CPZ','AFZ','PO9','FT10','PO10','FT9'};
original_labels = {'FP1', 'FZ', 'PZ', 'FP2', 'CPZ', 'AFZ'};
corrected_labels = {'Fp1', 'Fz', 'Pz', 'Fp2', 'CPz', 'AFz'};

% Iterate through pre and post files for each subject and condition 
for i = 1:length(pre_files)
    % Load pre EEGLAB data convert to Feildtrip Convert to SPM
    EEG_pre = pop_loadset('filename', pre_files(i).name, 'filepath', base_path);
    EEG_pre.data = EEG_pre.data(:,end-2499:end);  % the last 10s of data 
    data_pre_ft = eeglab2fieldtrip(EEG_pre, 'preprocessing', 'none');
    % Loop through the original labels and replace them with corrected labels
    for s = 1:numel(original_labels)
        idx = find(strcmp(data_pre_ft.label, original_labels{s}));
        if ~isempty(idx)
            data_pre_ft.label{idx} = corrected_labels{s};
        end
    end
    output_path_pre = fullfile(output_folder, ['_SPM_', pre_files(i).name]);
    D_pre = spm_eeg_ft2spm(data_pre_ft, output_path_pre);
    D_pre = meeg(D_pre);  % Ensure it's a MEEG object

    % Load post EEGLAB data
    EEG_post = pop_loadset('filename', post_files(i).name, 'filepath', base_path);
    EEG_post.data = EEG_post.data(:,1:2500);  % the first 10s of data 
    data_post_ft = eeglab2fieldtrip(EEG_post, 'preprocessing', 'none');
    % Loop through the original labels and replace them with corrected labels
    for o = 1:numel(original_labels)
        idx = find(strcmp( data_post_ft.label, original_labels{o}));
        if ~isempty(idx)
             data_post_ft.label{idx} = corrected_labels{o};
        end
    end
    output_path_post = fullfile(output_folder, ['_SPM_', post_files(i).name]);
    D_post = spm_eeg_ft2spm(data_post_ft, output_path_post);
    D_post = meeg(D_post);  % Ensure it's a MEEG object

    % seting units and data type
    chanidx_pre = ismember(D_pre.chanlabels, chan64);
    D_pre = units(D_pre, find(chanidx_pre), 'uV');
    D_pre = type(D_pre, 'single');
    chanidx_post = ismember(D_post.chanlabels, chan64);
    D_post = units(D_post, find(chanidx_post), 'uV');
    D_post = type(D_post, 'single');

    % Extracting subject and condition info to read csv
    tokens = regexp(pre_files(i).name, 'sub(\d+)_pre_(\w+)_ICA_clean.set', 'tokens');
    subject_num = tokens{1}{1};
    condition = tokens{1}{2};
    csv_filename = fullfile(Elc_path, ['sub', subject_num, '_', condition, '_electrodes.csv']);
    electrode_data = readtable(csv_filename);

    % Converting electrode labels to match expected format
    electrode_data.name = strrep(electrode_data.name, 'FP1', 'Fp1');
    electrode_data.name = strrep(electrode_data.name, 'FP2', 'Fp2');
    electrode_data.name = strrep(electrode_data.name, 'FZ', 'Fz');
    electrode_data.name = strrep(electrode_data.name, 'PZ', 'Pz');
    electrode_data.name = strrep(electrode_data.name, 'CPZ', 'CPz');
    electrode_data.name = strrep(electrode_data.name, 'AFZ', 'AFz');
    
    % Extracting and seting fiducials
    fiducial_names = {'Nasion', 'LPAP', 'RPAP'};
    fiducial_short_labels = {'nas', 'lpa', 'rpa'};
    fiducial_coords = zeros(3,3); 
    labels = cell(3,1);  
    for j = 1:length(fiducial_names)
        idx = strcmp(electrode_data.name, fiducial_names{j});
        fiducial_coords(j, :) = [electrode_data.x(idx), electrode_data.y(idx), electrode_data.z(idx)];
        labels{j} = fiducial_short_labels{j};
    end
    fid_struct = struct('label', {labels}, 'pnt', fiducial_coords);
    newfiducials = struct('fid', fid_struct, 'unit', 'mm', 'pnt', fiducial_coords);

    D_pre = fiducials(D_pre, newfiducials);
    D_post = fiducials(D_post, newfiducials);

    % Extracting and seting sensor positions
    sensor_labels = electrode_data.name(ismember(electrode_data.name, chan64));
    sensor_positions = electrode_data{ismember(electrode_data.name, chan64), {'x', 'y', 'z'}};
    new_sensors = struct('label', {sensor_labels}, 'pos', sensor_positions, 'type', 'eeg1010');

    D_pre = sensors(D_pre, 'EEG', new_sensors);
    D_post = sensors(D_post, 'EEG', new_sensors);

    % new MEEG object with the extra dimentions to combine pre and post 
    D_combined = clone(D_pre, ['_SPM_', 'sub', subject_num, '_', condition, '_combined'], [D_pre.nchannels, D_pre.nsamples, 2]);
    % Combining D_pre and D_post
    combined_data = cat(3, D_pre(:,:,:), D_post(:,:,:));
    D_combined(:,:,:) = combined_data;
    
    % Updating the trial info
    trial_labels = {'pre', 'post'};
    D_combined = conditions(D_combined, 1:2, trial_labels);
    
    % Save the combined MEEG object
    save(D_combined);
    save(D_pre);
    save(D_post);

%-----------------------------------these steps are not necessary
    % .sfp electrode and fiducial file
    combined_sfp_filename = fullfile(Elc_path, ['sub', subject_num, '_', condition, '_electrodes_fiducials_combined.sfp']);
    combined_file = fopen(combined_sfp_filename, 'w');
    % sensor data
    for j = 1:length(sensor_labels)
        fprintf(combined_file, '%-6s %10.3f %10.3f %10.3f\n', sensor_labels{j}, sensor_positions(j, 1), sensor_positions(j, 2), sensor_positions(j, 3));
    end   
    % fiducial data
    for j = 1:length(labels)
        fprintf(combined_file, '%-6s %10.3f %10.3f %10.3f\n', labels{j}, fiducial_coords(j, 1), fiducial_coords(j, 2), fiducial_coords(j, 3));
    end
%-------------------------------the previous step is not necessary
 
end
restoredefaultpath;
path(original_path);
