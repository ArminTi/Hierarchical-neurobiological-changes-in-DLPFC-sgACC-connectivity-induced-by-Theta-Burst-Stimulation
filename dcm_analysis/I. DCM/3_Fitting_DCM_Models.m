clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
restoredefaultpath;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg')
% --------------------------- Define Paths ---------------------------
file_path = 'D:/signal_data/DCM_Specified/doing';
file_name = 'DCM_SPM_sub-*_ses-*_epoch-*.mat';  
% Or choose a specific file to be processed
%file_name = 'DCM_SPM_sub-09_ses-02_epoch-03_converted_data_itbs.mat';
files = dir(fullfile(file_path, file_name));
cd(file_path)

addpath('D:/signal_data/all_spm_obj/doing');
addpath('D:/signal_data/code_looped'); % function of save_dcm_original

output_folder = 'D:/signal_data/DCM_fitted/doing';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
tic

%% ----------------------------- Fitting the DCM to data  --------------------------------
% paralell loop for fitting 
parfor idcm = 1:length(files)
    files(idcm)
    tmp = load(files(idcm).name, 'DCM');
    disp('data loaded')
    DCM = tmp.DCM;
    disp('DCM extracted')
    DCM = spm_dcm_csd(DCM);
    save_dcm_original(output_folder, ['fitted_' files(idcm).name(1:end-4)], DCM);
end
toc

disp('fitting for all subjects are done')


