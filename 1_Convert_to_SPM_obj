clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
%run feildtrip (it is recomended to not be directly add to perminant paths)
original_path = path;
restoredefaultpath;

% Add FieldTrip 
addpath('C:/Users/Growth fire/Programs/Matlab plugins/fieldtrip-20240731');
ft_defaults; % Initialize FieldTrip

% Add SPM12
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg'); % Start SPM with EEG defaults

% Add EEGLAB 
addpath('C:/Users/Growth fire/Programs/Matlab plugins/eeglab2024.0');
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui'); % Start EEGLAB without GUI
%==========================================================================
%-------------------------Define file paths--------------------------------
file_path = 'E:/signal_data/set_data';
file_name = 'sub-*_*_*-rest_run-*_final.set';  
files = dir(fullfile(file_path, file_name));

elec_path = 'E:/signal_data/electrode_fid_table';
elec_file = 'sub-*_*-*_electrodes.csv';
elecs = dir(fullfile(elec_path, elec_file));

output_folder = 'E:/signal_data/all_spm_obj';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define channel labels for consistency with SPM format
chan64 = {'AFz','Fp1','Fp2','AF3','AF4','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8',...
          'CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','PO3','PO4','O1','O2','AF7','AF8','F5','F1','F2','F6',...
          'FC3','FCz','FC4','FT8','C5','C1','C2','C6','TP7','FT7','CP3','CPz','CP4','TP8','P5','P1','P2','P6',...
          'PO7','POz','PO8','Oz','P9','P10','TP9','FP1','FZ','PZ','FP2','CPZ','AFZ'};
original_labels = {'FP1', 'FZ', 'PZ', 'FP2', 'CPZ', 'AFZ'};
corrected_labels = {'Fp1', 'Fz', 'Pz', 'Fp2', 'CPz', 'AFz'};

%---------------------------------Convert-----------------------------------
for i = 1:length(files)
    % Load EEGLAB data 
    EEG = pop_loadset('filename', files(i).name, 'filepath', file_path);
    % Convert to FieldTrip structure
    data_ft = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

    %% Correct channel labels if necessary
    for s = 1:numel(original_labels)
        idx = find(strcmp(data_ft.label, original_labels{s}));
        if ~isempty(idx)
            data_ft.label{idx} = corrected_labels{s};
        end
    end

    %% Extract sensor information from FieldTrip structure before conversion
    sensor_labels = data_ft.label;               % Sensor labels
    sensor_positions = data_ft.elec.elecpos  ;     % Sensor positions (electrode locations)

    % Convert FieldTrip data to SPM MEEG object
    output_path = fullfile(output_folder, ['SPM_', files(i).name]);    %%saving in the output location
    D = spm_eeg_ft2spm(data_ft, output_path);
    D = meeg(D);  % Ensure it's a MEEG object

    %% Extracting subject and condition info to read csv
    tokens = regexp(files(i).name, 'sub-(\d+)_(\w+)_\w+-rest_run-\d+_final.set', 'tokens');  
    subject_num = tokens{1}{1};
    modality = tokens{1}{2};  
    csv_pattern = fullfile(elec_path, ['sub-', subject_num, '_', modality, '*_electrodes.csv']);    
    csv_files = dir(csv_pattern);
    csv_filename = fullfile(elec_path, csv_files(1).name);

    electrode_data = readtable(csv_filename);
    fiducial_names = {'Nasion', 'LPAP', 'RPAP'};
    fiducial_short_labels = {'nas', 'lpa', 'rpa'};
    fiducial_coords = zeros(3,3); 
    labels = cell(3,1);  
    
    for j = 1:length(fiducial_names)
        idx = strcmp(electrode_data.name, fiducial_names{j});
        fiducial_coords(j, :) = [electrode_data.x(idx), electrode_data.y(idx), electrode_data.z(idx)];
        labels{j} = fiducial_short_labels{j};
    end

    %% Create fiducial structure
    fid_struct = struct('label', {labels}, 'pnt', fiducial_coords);

    %% Add fiducials to the MEEG object
    newfiducials = struct('fid', fid_struct, 'unit', 'mm', 'pnt', fiducial_coords);
    D = fiducials(D, newfiducials);

    %% Add sensor information to the MEEG object (using sensor labels and positions from FieldTrip)
    D = sensors(D, 'EEG', struct('label', {sensor_labels}, 'pos', sensor_positions, 'type', 'eeg1010'));

    %% Set the units and datatype of the MEEG object
    chanidx = ismember(D.chanlabels, chan64);  % Find matching channel indices
    D = units(D, find(chanidx), 'uV');  % Set units to microvolts
    D = type(D, 'single');  % Set datatype

    % Save the updated MEEG object
    save(D);
end
% Restore the original path
restoredefaultpath;
path(original_path);
fprintf('convasion to SPM object completed for all files.\n');

