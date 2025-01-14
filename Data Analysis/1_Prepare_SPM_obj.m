clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
% since SPM12 and FieldTrip have conflicting functions we want to restore
% the paths to its original after running our program so that they are not
% both in the general paths at the same time 
original_path = path;
restoredefaultpath;

% Add FieldTrip 
addpath('C:/Users/Growth fire/Programs/Matlab plugins/fieldtrip-20240731');
ft_defaults; % Initialize FieldTrip

% Add SPM12
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg'); % Start SPM with EEG defaults

% --------------------------- Define Paths ---------------------------

% Define the root directory containing the subfolders
root_dir = 'D:/signal_data/new_preprocessing_rawdata/doing';
file_list = dir(fullfile(root_dir, '**', '*.mat'));
% if you want to specifiy a single file you can use this instead
%file_list = dir(fullfile(root_dir, '**', 'sub-01_ses-01_epoch-01_converted_data.mat'));
length(file_list)

if isempty(file_list)
    error('No .mat files found in the specified directory or its subdirectories.');
end

fprintf('Found %d .mat files:\n', length(file_list));
for i = 1:length(file_list)
    fprintf('%s\n', fullfile(file_list(i).folder, file_list(i).name));
end

elec_path = 'D:/signal_data/raw_blocked_data/19297701';
elec_file = 'sub-*_*-*_electrodes.csv';
elecs = dir(fullfile(elec_path, elec_file));

output_folder = 'D:/signal_data/all_spm_obj/doing';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Load the data that shows participants recived which modality in that session
participants_csv = 'D:/signal_data/new_preprocessing_rawdata/participants.csv'; 
participants_data = readtable(participants_csv);

% Defining channel labels for consistency with SPM format
%'FPz' is the reference channel
chan64 = {'AFz','Fp1','Fp2','AF3','AF4','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8',...
          'CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','PO3','PO4','O1','O2','AF7','AF8','F5','F1','F2','F6',...
          'FC3','FCz','FC4','FT8','C5','C1','C2','C6','TP7','FT7','CP3','CPz','CP4','TP8','P5','P1','P2','P6',...
          'PO7','POz','PO8','Oz','P9','P10','TP9','FP1','FZ','PZ','FP2','CPZ','AFZ','FT9','FT10','PO9','PO10'};
original_labels = {'FP1', 'FZ', 'PZ', 'FP2', 'CPZ', 'AFZ'};
corrected_labels = {'Fp1', 'Fz', 'Pz', 'Fp2', 'CPz', 'Afz'}; % SPM is sensetive to name capitilization

%% ================== adding channel locations fiducials and session modality ==============
% Loops through all found .mat files
for i = 1:length(file_list)
    mat_file_path = fullfile(file_list(i).folder, file_list(i).name);
    fprintf('Processing file: %s\n', mat_file_path);
    D = spm_eeg_load(mat_file_path);

    % Sets units and data type
    chanidx = ismember(D.chanlabels, chan64); % Matchs channel indices
    D = units(D, find(chanidx), 'uV');        % Sets units to microvolts
    D = type(D, 'single');                    % Sets datatype

    %% Extracting subject csv for Electrode and Fid info
    tokens = regexp(file_list(i).name, 'sub-(\d+)_ses-(\d+)_epoch-(\d+)_converted_data\.mat', 'tokens');
    subject_num = tokens{1}{1}; % The subject number (e.g., '01')
    session_num = tokens{1}{2}; % The session number (e.g., '03')
    epoch_num = tokens{1}{3};   % The epoch number (e.g., '04')
    subj_num = str2double(tokens{1}{1}); % Convert subject number to numeric
    ses_num = str2double(tokens{1}{2}); % Convert session number to numeric

    csv_pattern = fullfile(elec_path, ['sub-', subject_num], ['ses-', session_num], ...
    ['sub-', subject_num, '_ses-', session_num, '_electrodes.csv']);
    csv_files = dir(csv_pattern);
    if ~isempty(csv_files)
        csv_filename = fullfile(csv_files(1).folder, csv_files(1).name); % Full path to the CSV file
        fprintf('Found CSV file: %s\n', csv_filename);
    else
        csv_pattern
        error('No electrode CSV file found for subject %s, session %s.\n', subject_num, session_num);
    end;
    csv_filename = fullfile(csv_files.folder, csv_files.name);
    %% Adding Fiducial locations 
    electrode_data = readtable(csv_pattern);
    fiducial_names = {'Nasion', 'LPAP', 'RPAP'};
    fiducial_short_labels = {'nas', 'lpa', 'rpa'};% spm is sensetive to capitalization
    fiducial_coords = zeros(3,3); 
    labels = cell(3,1);  
    
    for j = 1:length(fiducial_names)
        idx = strcmp(electrode_data.name, fiducial_names{j});
        fiducial_coords(j, :) = [electrode_data.x(idx), electrode_data.y(idx), electrode_data.z(idx)];
        labels{j} = fiducial_short_labels{j};
    end

    fid_struct = struct('label', {labels}, 'pnt', fiducial_coords);
    newfiducials = struct('fid', fid_struct, 'unit', 'mm', 'pnt', fiducial_coords);
    D = fiducials(D, newfiducials);
    %% Extracting and setting sensor locations from CSV
    % Load electrode data
    electrode_data = readtable(csv_pattern);

    % Standardize electrode labels to match expected format (spm is
    % sensetive to capitilization)
    electrode_data.name = strrep(electrode_data.name, 'FP1', 'Fp1');
    electrode_data.name = strrep(electrode_data.name, 'FP2', 'Fp2');
    electrode_data.name = strrep(electrode_data.name, 'FZ', 'Fz');
    electrode_data.name = strrep(electrode_data.name, 'PZ', 'Pz');
    electrode_data.name = strrep(electrode_data.name, 'CPZ', 'CPz');
    electrode_data.name = strrep(electrode_data.name, 'AFZ', 'AFz');

    % Extract the channel labels in the exact order from D.channels
    channel_order = D.chanlabels; % Retrieves the channel labels as a cell array

    % Initialize new sensor positions and labels arrays
    ordered_sensor_positions = nan(length(channel_order), 3);
    ordered_sensor_labels = cell(length(channel_order), 1);

    % Loop through the desired channel order and find the corresponding positions
    for k = 1:length(channel_order)
        channel_label = channel_order{k};
        idx = strcmp(electrode_data.name, channel_label); % Find matching label in CSV data
        if any(idx)
            ordered_sensor_positions(k, :) = [electrode_data.x(idx), electrode_data.y(idx), electrode_data.z(idx)];
            ordered_sensor_labels{k} = channel_label;
        else
            warning('Channel %s not found in the electrode CSV file.', channel_label);
        end
    end

    % Check for empty cells and print a message if any are found
    if any(cellfun(@isempty, ordered_sensor_labels))
        fprintf('Empty cells found in ordered_sensor_labels.\n');
    end

    % Remove empty rows in case some labels were not found
    valid_idx = ~cellfun(@isempty, ordered_sensor_labels); % Find non-empty labels
    ordered_sensor_positions = ordered_sensor_positions(valid_idx, :);
    ordered_sensor_labels = ordered_sensor_labels(valid_idx);


    % Create the new sensors structure and add the locations
    new_sensors = struct('label', {ordered_sensor_labels}, 'pos', ordered_sensor_positions, 'type', 'eeg1010');
    D = sensors(D, 'EEG', new_sensors);

    %% Find the modality based on the session
    session_column = sprintf('s%d_cond', ses_num); 
    modality_idx = find(participants_data.participant == subj_num); 

    if isempty(modality_idx)
        error('Subject %d not found in the participants table.', subj_num);
    end
    
    if ismember(session_column, participants_data.Properties.VariableNames)
        modality = participants_data.(session_column){modality_idx}; % Gets the modality (e.g., 'iTBS', 'cTBS', or 'sham')
    else
        error('Session column %s not found in the participants table.', ses_column);
    end

    fprintf('Subject: %d, Session: %d, Modality: %s\n', subj_num, ses_num, modality);
    fprintf('Saving file to: %s\n', output_folder );

    % Prepare the output filename
    [~, name, ~] = fileparts(file_list(i).name); % Extracts file name without extension
    output_filename = sprintf('SPM_%s_%s.mat', name, modality); % Adds modality to the filename
    output_path = fullfile(output_folder, output_filename); % Full path to the file

    % Copys the SPM object to the new location
    D_copy = D.copy(output_path);

    fprintf('Saved processed file: %s\n', output_path);

end

% Restores the original path (we don't want spm and FieldTrip functions to conflict)
restoredefaultpath;
path(original_path);
fprintf('Conversion to SPM object completed for all files.\n');

