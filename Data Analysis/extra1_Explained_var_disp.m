clear all;
close all;
clc;
% ------------------------------- path ----------------------------------
file_path = 'E:/signal_data/DCM_fitted';
file_name = 'fitted_DCM_SPM_sub-*_*_*-rest_run-*_final.mat';  % Specific file to be processed
files = dir(fullfile(file_path, file_name));
cd(file_path)
addpath('E:/signal_data/code_looped');

% -------------------- visualize explaied varience ----------------------
for i = 1:length(files)
    disp(files(i).name)
    tmp = load(files(i).name, 'DCM');  
    DCM = tmp.DCM;
    x = explained_var_dcm(DCM);
    disp(x);
end
