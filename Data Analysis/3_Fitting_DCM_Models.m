clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
restoredefaultpath;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg')
% ------------------------------- path ----------------------------------
file_path = 'E:/signal_data/DCM_Specified';
file_name = 'DCM_SPM_sub-*_*_*-rest_run-*_final.mat'; 
files = dir(fullfile(file_path, file_name));
cd(file_path)
addpath('E:/signal_data/all_spm_obj');
addpath('E:/signal_data/code_looped'); % where save_dcm function is saved

output_folder = 'E:/signal_data/DCM_fitted';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% -------------- paralell loop for fitting the DCM to data -----------------
tic
parfor idcm = 1:length(files)
        tmp = load(files(idcm).name, 'DCM');  
        DCM = tmp.DCM;
        DCM = spm_dcm_csd(DCM);
        save_dcm(output_folder, ['fitted_' files(idcm).name(1:end-4)], DCM);
end
toc
disp('fitting for all subjects are done')

