clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
original_path = path;
restoredefaultpath;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg');
% ----------------------------- path ---------------------------------
file_path = 'E:/signal_data/DCM_fitted';
file_name = 'fitted_DCM_SPM_sub-*_*_*-rest_run-*_final.mat';  
files = dir(fullfile(file_path, file_name));
cd(file_path)

output_folder = 'E:/signal_data/PEB_of_PEB_results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% ---------------------- Organize Files in GCMs --------------------------
files_pre_sham = {};
files_post_sham = {};
files_pre_iTBS = {};
files_post_iTBS = {};
files_pre_cTBS = {};
files_post_cTBS = {};

for i = 1:length(files)
    file_name = files(i).name;  % Get the current file name
    tokens = regexp(file_name, 'fitted_DCM_SPM_sub-(\d+)_(\w+)_(pre|post)-rest_run-\d+_final.mat', 'tokens');  
    subject_num = tokens{1}{1};  % Subject number
    modality = tokens{1}{2};     % Modality
    timing = tokens{1}{3};   
    
    if modality == 'sham'
        % Sham condition
        if strcmp(timing, 'pre')
            files_pre_sham{end+1,1} = fullfile(file_path, file_name);  % Pre-sham
        else
            files_post_sham{end+1,1} = fullfile(file_path, file_name);  % Post-sham
        end
    elseif modality == 'iTBS'
        % iTBS condition
        if strcmp(timing, 'pre')
            files_pre_iTBS{end+1,1} = fullfile(file_path, file_name);  % Pre-iTBS
        else 
            files_post_iTBS{end+1,1} = fullfile(file_path, file_name);  % Post-iTBS
        end
    elseif modality == 'cTBS'
        % cTBS condition
        if strcmp(timing, 'pre')
            files_pre_cTBS{end+1,1} = fullfile(file_path, file_name);  % Pre-cTBS
        else 
            files_post_cTBS{end+1,1} = fullfile(file_path, file_name);  % Post-cTBS
        end
    end
end

%% ------------------------ subject level PEB --------------------------

Ns = length(files_pre_sham);  
num_conditions = 2;  
col1 = ones(Ns * num_conditions, 1);
col2 = [ones(Ns, 1); -ones(Ns, 1)];
X = [col1, col2];

GCMiTBS = vertcat(files_pre_iTBS, files_post_iTBS);
PEBiTBS = spm_dcm_peb(GCMiTBS, X);

GCMcTBS = vertcat(files_pre_cTBS, files_post_cTBS);
PEBcTBS = spm_dcm_peb(GCMcTBS, X);

GCMsham = vertcat(files_pre_sham, files_post_sham);
PEBsham = spm_dcm_peb(GCMsham, X);
%% ------------------------ group level PEB of PEB --------------------------
X2 =  [1 -1  0  1; 
       1  1 -1  0;
       1  0  1 -1];
PEBs = {PEBiTBS; PEBcTBS;PEBsham};
PEB3 = spm_dcm_peb(PEBs,X2);
[BMA,BMR] = spm_dcm_peb_bmc(PEB3);
spm_dcm_peb_review(BMA,PEBs)

new_filename = ['PEB_of_Peb_result'];
save(fullfile(output_folder, new_filename), 'PEB3');
