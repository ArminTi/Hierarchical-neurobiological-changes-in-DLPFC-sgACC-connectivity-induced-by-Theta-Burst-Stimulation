%%%%%% specifiy a DCM for each subject
% clear all;clc;close all
%% go to folder were the MEEG objects are stored
clear all;
close all;
clc;
dir_feed = 'C:/Users/Growth fire/Documents/research/signal/updated data/out all/'
cd(dir_feed)
fold = dir(fullfile([dir_feed '*_combined.mat']));
clear DCM %--------------------------------------------------------
%#####################################################################
%%%% load a template DCM options
%load('/Users/fvdsteen/data/EEG_rest/DCM_options.mat'); %%%see github
%%%
spm('defaults','EEG');
DCM.options.analysis = 'CSD'; % analyze evoked responses -> (cross spectral density) neuronal innovations (rest EEG)
DCM.options.model    = 'NMM'; %> ERP model -> neuronal mass model conductance based
%Conductance-based dynamic causal modeling: A mathematical review of its
%application to cross-power spectral densities
DCM.options.spatial  = 'IMG'; %> spatial model -> local feild potential
DCM.options.trials   = [1 2]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled-> time window start
DCM.options.Tdcm(2)  = 2500;   % end of peri-stimulus time to be modelled-> time window end
%DCM.options.Nmodes   = 8;     %> nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
%DCM.options.onset    = 60;    % selection of onset (prior mean)??? 
DCM.options.D        = 1;     % downsampling
%#####################################################################

for i = 1:length(fold)
   
    disp(fold(i).name)
    
    DCM.xY.Dfile =  [fold(i).name(1:end-4)]
    %DCM.options=options;% set the correct DCM options################
    
    DCM.options.mod_int = 1;%%%% model intrinsic connectivity changes
    %################################################################
    %DCM.Lpos =[0 -88 4; -42 -72 0; 44 -68 0]'; %%% put the sources positions (in MNI)
    %DCM.Sname = {'V1' 'lV5' 'rV5'}%% names sources
    DCM.Lpos  = [[-46; -66; 30] [49; -63; 33] [0 ;-58; 0] [-1 ;54; 27] ];
    DCM.Sname = {'left LP', 'right LP', 'posterior cingulate', 'medial prefrontal cortex'};

    %DCM.A{1} = [1 0 0;1 0 0;1 0 0]%% set forward connections
    %DCM.A{2} = [0 1 1; 0 0 0;0 0 0]%% set backward connections
    %DCM.A{3} = zeros(3,3)%% this is currently not used but needs to be set for the code to properly run
    %DCM.B{1} =ones(3,3);%% allow extrinsic connectivity modulations
    %DCM.C = sparse(3,0);%%% no external input
    Nareas    = size(DCM.Lpos,2);
    DCM.A{1} = [0 0 1 1;
                0 0 1 1;
                0 0 0 1;
                0 0 0 0]; %Forward model

    DCM.A{2} = [0 0 0 0;
                0 0 0 0;
                1 1 0 0;
                1 1 1 0]; %Backward model

    DCM.A{3} = [0 1 0 0;
                1 0 0 0;
                0 0 0 0;
                0 0 0 0]; %Lateral
    DCM.B{1} = DCM.A{1} + DCM.A{2}; % input on node interactions model
    DCM.B{1}(:) = 0;
    DCM.C = [0; 0; 0; 0]; % input directly on node model
    %###############################################################
  
    xU.X = [0;1]% 'design matrix' voor connectivity modulations
    DCM.xU = xU;
    DCM  = spm_dcm_erp_data(DCM); %% get the data  (see spm_dcm_csd)????????????????????? WHY ERP MODEL FOR REST
    DCM.M.dipfit.model = 'CMC';%% specific the 'CMC models'?????????????????????????????????????????????????????
    DCM.M.dipfit.type  = 'IMG'; %%% 'IMG' forward model (i.e. patch)????????????????????????????????????????????
    DCM.options.Nmodes = 4;%% number of data modes??????????????????????????????????????? why 4
    DCM = spm_dcm_erp_dipfit(DCM,1);%% get lead field etc. for forward models?????????????? WHY ERP MODEL FOR REST
    DCM = spm_dcm_csd_data(DCM);%% extract cross spectral densities as data features
    DCM.options.DATA = 0;%% put this option to 0, because we extracted the data feature outside this function
    DCM.M.nograph=0;%% supress visual feedback during model fitting----------------------------------------------
    %%% specify priors

    % neural priors
    [pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    %%% allow all intrisic connections to be estimated
    DCM.M.pE.G = sparse(length(DCM.Sname),10);
    DCM.M.pC.G = sparse(length(DCM.Sname),10)+1/8;
    %% priors on the forward model
    [DCM.M.pE,DCM.M.pC] = spm_L_priors(DCM.M.dipfit,pE,pC);
    DCM.M.pE.G = sparse(length(DCM.Sname),10);
    DCM.M.pC.G = sparse(length(DCM.Sname),10)+1/8;
    DCM.M.pE.L = zeros(6,3);
    DCM.M.pE.L(3,1) = 0.5;
    DCM.M.pE.L(1,2) = 0.5;
    DCM.M.pE.L(1,3) = 0.5;
   
    %% max number of VB iterations
    name = []; %-------------------------------------------------------------------------
    DCM.M.Nmax    = 256;
    %DCM.name = [name 'sub' num2str(i,'%.3u')];
    DCM.M.nograph = 1;
    %% save specified DCM.
    [~, original_filename, ~] = fileparts(fold(i).name);
    new_filename = ['DCM_' original_filename];
    save(['C:/Users/Growth fire/Documents/research/signal/updated data/DCM/' new_filename],'DCM')
end

% spm_dcm_csd(DCM)