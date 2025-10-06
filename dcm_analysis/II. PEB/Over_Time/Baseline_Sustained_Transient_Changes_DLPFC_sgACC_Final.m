%% Loading GCMs
clear all;
clc;
addpath('C:/Users/Documents/Programs/spm_25.01.02/spm');
spm('defaults', 'eeg'); 
load('GCM_ctbs.mat');
load('GCM_itbs.mat');
load("GCM_sham.mat");

%% Peb First Level

fields = {'A', 'AN', 'H','T','CV'}; % Parameters to include in PEB
sub = 22;
Peb_itbs = cell(sub,1);
Peb_ctbs = cell(sub,1);
Peb_sham = cell(sub,1);
pattern_baseline =  [1;  1;  1;  1];
pattern_sustained = [-1; 1; 1;  1];
pattern_transient = [-1; 1; 1; -1];

M = struct();
M.X = [pattern_baseline, pattern_sustained, pattern_transient];

for e = 1:sub
    DCM_itbs_epoch = GLM_itbs_complete(e, :)';
    DCM_ctbs_epoch = GLM_ctbs_complete(e, :)';
    DCM_sham_epoch = GLM_sham_complete(e, :)';
    PEB_itbs_epoch = spm_dcm_peb(DCM_itbs_epoch, M, fields);
    PEB_ctbs_epoch = spm_dcm_peb(DCM_ctbs_epoch, M, fields);
    PEB_sham_epoch = spm_dcm_peb(DCM_sham_epoch, M, fields);
    Peb_itbs{e} = PEB_itbs_epoch;
    Peb_ctbs{e} = PEB_ctbs_epoch;
    Peb_sham{e} = PEB_sham_epoch;
   
end

nSham = length(Peb_sham);
ncTBS = length(Peb_ctbs);
niTBS = length(Peb_itbs);

save('Peb_sham_30min', 'Peb_sham')
save('Peb_ctbs_30min', "Peb_ctbs")
save('Peb_itbs_30min', "Peb_itbs")
Peb_cTBSvsiTBS = [Peb_ctbs; Peb_itbs];
Peb_cTBSvsSham = [Peb_sham; Peb_ctbs];
Peb_iTBSvsSham = [Peb_sham; Peb_itbs];

%% peb second level: itbs vs ctbs
clear M

M = struct();
M.X = [ones(ncTBS+niTBS,1),[repmat(-1,[ncTBS,1]); repmat(+1,[niTBS,1])]];

cTBSvsiTBS_final = spm_dcm_peb(Peb_cTBSvsiTBS, M, fields);
BMA_cTBSvsiTBS = spm_dcm_peb_bmc(cTBSvsiTBS_final)
save( "BMA_cTBSvsiTBS_BST_-1c_1i",'BMA_cTBSvsiTBS')
%spm_dcm_peb_review(BMA_cTBSvsiTBS);

%% peb second level: itbs
clear M

M = struct();
M.X = [ones(ncTBS,1)];

iTBS_final = spm_dcm_peb(Peb_itbs, M, fields);
BMA_itbs_final = spm_dcm_peb_bmc(iTBS_final);
save( 'BMA_itbs_final_BST_n1i','BMA_itbs_final')
%spm_dcm_peb_review(BMA_itbs_final);
%% peb second level: ctbs

cTBS_final = spm_dcm_peb(Peb_ctbs, M, fields);
BMA_ctbs_final = spm_dcm_peb_bmc(cTBS_final);
save( 'BMA_ctbs_final_BST_n1c', 'BMA_ctbs_final')
%spm_dcm_peb_review(BMA_ctbs_final);

%% peb second level: sham

sham_final = spm_dcm_peb(Peb_sham, M, fields);
BMA_sham_final = spm_dcm_peb_bmc(sham_final);
save( 'BMA_sham_final_BST_n1s', 'BMA_sham_final')
%spm_dcm_peb_review(BMA_sham_final);

%% peb second level: itbs vs sham
clear M

M = struct();
M.X = [ones(niTBS+nSham,1),[repmat(-1,[niTBS,1]); repmat(+1,[niTBS,1])]];

iTBSvsSham_final = spm_dcm_peb(Peb_iTBSvsSham, M, fields);
BMA_itbs = spm_dcm_peb_bmc(iTBSvsSham_final);
save( 'BMA_itbs_BST_-1s_1i','BMA_itbs')


%% peb second level: ctbs vs sham
clear M

M = struct();
M.X = [ones(ncTBS+nSham,1),[repmat(-1,[ncTBS,1]); repmat(+1,[ncTBS,1])]];

cTBSvsSham_final = spm_dcm_peb(Peb_cTBSvsSham, M, fields);
BMA_ctbs = spm_dcm_peb_bmc(cTBSvsSham_final);
save('BMA_ctbs_BST_-1s_1c','BMA_ctbs')

%% Review

load("BMA_cTBS_winning_model.mat")
load("BMA_iTBS_winning_model.mat")

spm_dcm_peb_review(BMA_itbs);
input("enter to continue")
spm_dcm_peb_review(BMA_ctbs);

