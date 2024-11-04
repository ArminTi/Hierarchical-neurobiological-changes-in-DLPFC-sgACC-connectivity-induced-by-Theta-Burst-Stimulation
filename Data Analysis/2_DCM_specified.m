clear all;
close all;
clc;
% ----------------------------- Setup --------------------------------
original_path = path;
restoredefaultpath;
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg'); % Start SPM with EEG defaults
% --------------------------- Define Paths ---------------------------
file_path = 'E:/signal_data/all_spm_obj';
file_name = 'SPM_sub-*_*_*-rest_run-*_final.mat';  
files = dir(fullfile(file_path, file_name));
cd(file_path)

output_folder = 'E:/signal_data/DCM_Specified';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% ------------------------- Define DCM Model ---------------------------
%#####################################################################
spm('defaults','EEG');
DCM.options.analysis = 'CSD'; % step 1: Observation modesl (cross spectral density) for resting EEG
DCM.options.model    = 'CMM_NMDA'; % neuronal model
% Conductance based version of canonical microcircuates (CMC) with NMDA receptors (4 neural population)
DCM.options.trials   = 1;
% data and design step5
DCM.options.Tdcm(1)  = 0;   % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 2000;   % end of peri-stimulus time to be modelled
DCM.options.Fdcm(1)  = 1;
DCM.options.Fdcm(2)  = 40;
DCM.options.Nmodes   = 8;     %> nr of modes for data selection 
DCM.options.h        = 1;     % nr of DCT components
DCM.options.D        = 1;     % downsampling
DCM.options.spatial  = 'IMG'; % spatial model
%#####################################################################
% -----------------------------  Loop --------------------------------
for i = 1:length(files)
    disp(files(i).name)
    DCM.xY.Dfile =  [files(i).name(1:end-4)];
    DCM.options.mod_int = 1;
    % ------------------------- connections ---------------------------
    %#################################################################
    DCM.Lpos  = [[-46; -66; 30] [49; -63; 33] [0 ;-58; 0] [-1 ;54; 27] ];
    DCM.Sname = {'lLP', 'rLP', 'PC', 'MPC'};
    %'left LP', 'right LP', 'posterior cingulate', 'medial prefrontal cortex'
    Nareas    = size(DCM.Lpos,2);
    DCM.A{1} = [1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1]; %Forward model

    DCM.A{2} = [1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1]; %Backward model

    DCM.A{3} = [1 1 1 1;
                1 1 1 1;
                1 1 1 1;
                1 1 1 1]; %Lateral
    DCM.B = cell(0, 0);
    DCM.C = [];
    %#################################################################
    %design matrix
    xU.name = '';
    xU.X = [];   
    DCM.xU = xU; 
    fprintf('Processing file %d of %d: %s\n', i, length(files), files(i).name);
    % ------------------- Specify DCM Data File ------------------------
    DCM  = spm_dcm_erp_data(DCM);
    % ------------------- Specify Neuronal Model -----------------------
    DCM.M.dipfit.model = 'CMM_NMDA'; 
    DCM.M.dipfit.type  = 'IMG';
    DCM.M.nograph=1;
    % -------------------- Specify Priors -----------------------------
    % neural priors
    [pE,pC] = spm_cmm_NMDA_priors(DCM.A,DCM.B,DCM.C);
    DCM.M.pE = pE; % prior
    DCM.M.pC = pC; % prior covarience
    % state and observation equations
    DCM = spm_dcm_erp_dipfit(DCM, 1 );
    % CSD data
    DCM = spm_dcm_csd_data(DCM);
    
    %% priors on the forward model 
    [DCM.M.pE,DCM.M.pC] = spm_L_priors(DCM.M.dipfit,pE,pC);
    %% max number of VB iterations
    name = []; 
    DCM.M.Nmax    = 256;
    DCM.name = [files(i).name];
    DCM.M.nograph = 1;
    %% save specified DCM.
    [~, original_filename, ~] = fileparts(files(i).name);
    new_filename = ['DCM_' original_filename];
    save(fullfile(output_folder, new_filename), 'DCM');
end
restoredefaultpath;
fprintf('DCM specification completed for all files.\n');

