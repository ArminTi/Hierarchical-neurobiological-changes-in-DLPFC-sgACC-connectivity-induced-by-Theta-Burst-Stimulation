%% Loading GCMs
clear all;
clc;
addpath('C:/Users/ASUS/Desktop/Apps/spm12');
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
pattern_baseline =  [1;  1];
pattern_change = [-1; 1];

M = struct();
%M.X = [pattern_baseline, pattern_sustained, pattern_transient];
M.X = [pattern_baseline, pattern_change];
for e = 1:sub
    DCM_itbs_epoch = GLM_itbs_complete(e, 1:2)';
    DCM_ctbs_epoch = GLM_ctbs_complete(e, 1:2)';
    DCM_sham_epoch = GLM_sham_complete(e, 1:2)';
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

Peb_cTBSvsSham = [Peb_sham; Peb_ctbs];
Peb_iTBSvsSham = [Peb_sham; Peb_itbs];

save('Peb_sham', "Peb_sham")
save('Peb_ctbs', "Peb_ctbs")
save('Peb_itbs', "Peb_itbs")
%% peb second; effect of sham
clear M
M = struct();
M.X = [ones(ncTBS,1)];
[Sham_final, Peb_sham_update] = spm_dcm_peb(Peb_sham,M, fields);
BMA_Sham_average = spm_dcm_peb_bmc(Sham_final);
save('BMA_Sham_average', "BMA_Sham_average")
save('Peb_sham_update', "Peb_sham_update")

%% peb second; effect of cTBS

[cTBS_final,Peb_ctbs_update]  = spm_dcm_peb(Peb_ctbs,M, fields);
BMA_cTBS_average = spm_dcm_peb_bmc(cTBS_final);
save('BMA_cTBS_average', "BMA_cTBS_average")
save('Peb_ctbs_update', "Peb_ctbs_update")

%% peb second; effect of iTBS

[iTBS_final,Peb_itbs_update] = spm_dcm_peb(Peb_itbs,M, fields);
BMA_iTBS_average = spm_dcm_peb_bmc(iTBS_final);
save('BMA_iTBS_average', "BMA_iTBS_average")
save('Peb_itbs_update', "Peb_itbs_update")

%% Review

load("BMA_cTBS_average.mat")
load("BMA_iTBS_average.mat")
load("BMA_Sham_average.mat")
spm_dcm_peb_review(BMA_cTBS_average);
input("enter to continue")
spm_dcm_peb_review(BMA_iTBS_average);
input("enter to continue")
spm_dcm_peb_review(BMA_Sham_average);

%% Individual difference in H(4,2,1) cTBS, and H(3,4,1) iTBS

load("Peb_sham.mat")
load("Peb_ctbs.mat")
load("Peb_itbs.mat")
%% Individual difference in H(4,2,1) cTBS, and H(3,4,1) iTBS after updating

load("Peb_sham_update.mat")
load("Peb_ctbs_update.mat")
load("Peb_itbs_update.mat")

%% H(4,2,1) for cTBS

paramName = 'H(4,2,1)';
numSubs = length(Peb_ctbs_update);


all_Ep = nan(numSubs,1);
all_Pp = nan(numSubs,1);

for s = 1:numSubs
    peb = Peb_ctbs_update{s};
    idx = find(strcmp(peb.Pnames, paramName));
    thisEp = peb.Ep(idx, 2);
    diagCp = diag(peb.Cp);
    thisCp = diagCp(idx*2);
    T = 0;
    thisPp = 1 - spm_Ncdf(T, abs(thisEp), thisCp);
    all_Ep(s) = thisEp;
    all_Pp(s) = thisPp;
end

% Display results for all subjects
fprintf('Parameter %s (Covariate 2):\n', paramName);
disp(table((1:numSubs)', all_Ep, all_Pp, 'VariableNames', {'Subject','Ep','Pp'}));


% Plotting


posCount = sum(all_Ep > 0);
negCount = sum(all_Ep < 0);
total = posCount + negCount;
percentPos = (posCount / total) * 100;
percentNeg = (negCount / total) * 100;
labels = {sprintf('Positive (%.1f%%)', percentPos), sprintf('Negative (%.1f%%)', percentNeg)};
colors = [1 0 0; 0 0 1]; % Blue, Red

% Plot pie chart
figure;
pie([posCount, negCount]);
colormap(colors); % Apply custom colors
legend(labels, 'Location', 'best');
title('Distribution of Ep Values');
%% H(3,4,1) iTBS
clear all_Pp all_Ep posCount negCount total percentPos percentNeg
paramName = 'H(3,4,1)';
numSubs = length(Peb_itbs_update);

% Preallocate arrays to store Ep and Pp for each subject
all_Ep = nan(numSubs,1);
all_Pp = nan(numSubs,1);

for s = 1:numSubs
    peb = Peb_itbs_update{s};
    idx = find(strcmp(peb.Pnames, paramName));
    thisEp = peb.Ep(idx, 2);
    diagCp = diag(peb.Cp);
    thisCp = diagCp(idx*2);
    T = 0;
    thisPp = 1 - spm_Ncdf(T, abs(thisEp), thisCp);
    all_Ep(s) = thisEp;
    all_Pp(s) = thisPp;
end


% Plotting


posCount = sum(all_Ep > 0);
negCount = sum(all_Ep < 0);
total = posCount + negCount;
percentPos = (posCount / total) * 100;
percentNeg = (negCount / total) * 100;
labels = {sprintf('Positive (%.1f%%)', percentPos), sprintf('Negative (%.1f%%)', percentNeg)};
colors = [1 0 0; 0 0 1]; % Blue, Red

% Plot pie chart
figure;
pie([posCount, negCount]);
colormap(colors); % Apply custom colors
legend(labels, 'Location', 'best');
title('Distribution of Ep Values');

