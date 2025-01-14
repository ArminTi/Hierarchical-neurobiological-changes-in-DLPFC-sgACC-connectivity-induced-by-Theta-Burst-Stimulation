clear all;
close all;
clc;
% There is a manual code to check files one by one to make sure the head
% models are being efined correctly the indevidual fiducials and MRI data
% are working correctly once we determined this for each session of each
% subject we used the dcmspecify_looped_auto to automaticaly go through all
% the subjects 
% ----------------------------- Setup --------------------------------
original_path = path;
restoredefaultpath;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg'); % Start SPM with EEG defaults
% --------------------------- Define Paths ---------------------------
file_path = 'D:/signal_data/all_spm_obj/doing';
% if you want to choose multiple files you can use wild card
%file_name = 'SPM_sub-*_ses-*_epoch-*_converted_data_*.mat';  
file_name = 'SPM_sub-22_ses-04_epoch-03_converted_data_ctbs.mat';
files = dir(fullfile(file_path, file_name));
cd(file_path)

output_folder = 'D:/signal_data/DCM_Specified/doing';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% ------------------------- Define DCM Model ---------------------------
%#####################################################################
spm('defaults','EEG');
DCM.options.analysis = 'CSD'; % Observation modesl (cross spectral density) for resting EEG
DCM.options.model    = 'CMM_NMDA'; % Neuronal model
% Conductance based version of canonical microcircuates (CMC) with NMDA receptors (4 neural population)
% data and design 
DCM.options.trials   = 1; 
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 2000;  % end of peri-stimulus time to be modelled    % [start endtime] 
DCM.options.Fdcm(1)  = 1;
DCM.options.Fdcm(2)  = 40;
DCM.options.Nmodes   = 8;     %> nr of modes for data selection 
DCM.options.h        = 1;     % nr of DCT components
DCM.options.D        = 1;     % downsampling
DCM.options.spatial  = 'IMG'; % spatial model
%#####################################################################
%% ============== looping through all the selected subjects ================
for i = 1:length(files)
    disp(files(i).name)
    DCM.xY.Dfile =  [files(i).name(1:end-4)];
    DCM.options.mod_int = 1;
    % ------------------------- connections ---------------------------
    %############################################################
    DCM.Lpos  = [[-46; -66; 30] [49; -63; 33] [0 ;-58; 0] [-1 ;54; 27] ];
    DCM.Sname = {'lLP', 'rLP', 'PC', 'MPC'};
    %'left LP', 'right LP', 'posterior cingulate', 'medial prefrontal cortex'
    Nareas    = size(DCM.Lpos,2);
    DCM.A{1} = [0 1 0 0;
                1 0 0 0;
                1 1 0 0;
                1 1 1 0]; %Forward model

    DCM.A{2} = [0 1 1 1;
                1 0 1 1;
                0 0 0 1;
                0 0 0 0]; %Backward model

    DCM.A{3} = [0 0 0 0;
                0 0 0 0;
                0 0 0 0;
                0 0 0 0]; %Lateral
    DCM.B = cell(0, 0); % rest data
    DCM.C = []; % rest data
    %############################################################
    %design matrix
    xU.name = '';
    xU.X = [];    %'design matrix' 
    DCM.xU = xU; 
    % -------------------------  Loop --------------------------

    fprintf('Processing file %d of %d: %s\n', i, length(files), files(i).name);
    % ------------------- Specify DCM Data File ------------------------

    DCM  = spm_dcm_erp_data(DCM); %% get the data

    % ------------------- Specify Neuronal Model -----------------------
    DCM.M.dipfit.model = 'CMM_NMDA'; %% specific the 'CMM models'
    DCM.M.dipfit.type  = 'IMG'; %%% 'IMG' forward model 
    DCM.M.nograph=1;%% supress visual feedback during model fitting
    % -------------------- Specify Priors -----------------------------
    %--- neural priors
    [pE,pC] = spm_cmm_NMDA_priors(DCM.A,DCM.B,DCM.C);
    DCM.M.pE = pE; % prior
    DCM.M.pC = pC; % prior covarience
    DCM.M.pE.G = sparse(length(DCM.Sname),10);
    DCM.M.pC.G = sparse(length(DCM.Sname),10)+1/8;
    % Set state and observation equations
    DCM = spm_dcm_erp_dipfit(DCM, 1 );%% get lead field etc. for forward models
    % CSD data
    DCM = spm_dcm_csd_data(DCM); % Estimating cross-spectral data
    %--- priors on the forward model 
    [DCM.M.pE,DCM.M.pC] = spm_L_priors(DCM.M.dipfit,pE,pC);
    DCM.M.pE.G = sparse(length(DCM.Sname),10);
    DCM.M.pC.G = sparse(length(DCM.Sname),10)+1/8;
    %--- max number of VB iterations
    name = []; 
    DCM.M.Nmax    = 256;
    DCM.name = [files(i).name];
    DCM.M.nograph = 1;
    % -------------------- save specified DCM file -----------------------------
    [~, original_filename, ~] = fileparts(files(i).name);
    new_filename = ['DCM_' original_filename];
    save(fullfile(output_folder, new_filename), 'DCM');
end

restoredefaultpath;
fprintf('DCM specification completed for all files.\n');
