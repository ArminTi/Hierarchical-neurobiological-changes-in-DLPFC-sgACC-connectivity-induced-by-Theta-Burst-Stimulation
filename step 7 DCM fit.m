clear all;
close all;
clc;
cd  'C:/Users/Growth fire/Documents/research/signal/updated data/DCM/'
dir_feed = 'C:/Users/Growth fire/Documents/research/signal/updated data/DCM/'
fold = dir(fullfile([dir_feed '*.mat']))
%% paralell loop for fitting 
spm('defaults','EEG');

parfor idcm =1:length(fold)
   
        tmp = load(fold(idcm).name,'DCM');
        DCM = tmp.DCM;
        DCM = spm_dcm_csd(DCM);
        save_dcm([fold(idcm).name(1:end-4) '_results'],DCM);
    
end

%function [] = save_dcm(d_name,DCM)

%save(d_name,'DCM');

%end
