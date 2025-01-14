clear all;
close all;
clc;

% ----------------------------- Setup --------------------------------
original_path = path;
restoredefaultpath;

% Add SPM12
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm12/spm12');
spm('defaults', 'eeg'); % Start SPM with EEG defaults

% --------------------------- Define Paths ---------------------------
file_path = 'D:/signal_data/DCM_fitted/doing';
file_pattern = 'fitted_DCM_SPM_sub-*_ses-*_epoch-*_converted_data_*.mat';
files = dir(fullfile(file_path, file_pattern));
cd(file_path);

addpath('D:/signal_data/code_looped');% function of explained_var_dcm

% ------------------------ Organizing subjects ------------------------

% Identify subjects
subjectList = [];
for f = 1:length(files)
    fname = files(f).name;
    subj_expr = 'sub-(\d+)';
    subj_token = regexp(fname, subj_expr, 'tokens', 'once');
    if ~isempty(subj_token)
        subj_num = str2double(subj_token{1});
        subjectList(end+1) = subj_num;
    end
end

uniqueSubjects = unique(subjectList);
numSubjects = length(uniqueSubjects);

% we want epochs 2 to 5 (4 epochs total) to have sessiones after TMS
numEpochs = 4;
epochIndices = 2:5;

% Initialize GLM cell arrays (GCM for the subject level was named GLM)
GLM_sham = cell(numSubjects, numEpochs);
GLM_ctbs = cell(numSubjects, numEpochs);
GLM_itbs = cell(numSubjects, numEpochs);

% Create a mapping from subject numbers to row indices
subjectIndexMap = containers.Map(uniqueSubjects, 1:numSubjects);

% Fill structures
for f = 1:length(files)
    fname = files(f).name;
    % Parse subject number
    subj_expr = 'sub-(\d+)';
    subj_token = regexp(fname, subj_expr, 'tokens', 'once');
    if isempty(subj_token), continue; end
    subj_num = str2double(subj_token{1});
    subj_idx = subjectIndexMap(subj_num);

    % Parse epoch number
    epoch_expr = 'epoch-(\d+)';
    epoch_token = regexp(fname, epoch_expr, 'tokens', 'once');
    if isempty(epoch_token), continue; end
    epoch_num = str2double(epoch_token{1});
    if epoch_num < 2 || epoch_num > 5, continue; end

    % Identify modality
    modality = '';
    if contains(fname, 'sham.mat')
        modality = 'sham';
    elseif contains(fname, 'ctbs.mat')
        modality = 'ctbs';
    elseif contains(fname, 'itbs.mat')
        modality = 'itbs';
    else
        continue;
    end

    % Load the DCM file
    data = load(fullfile(file_path, fname));
    if ~isfield(data, 'DCM')
        warning(['No DCM field found in ' fname]);
        continue;
    end
    DCM = data.DCM;

    % Column index for epoch
    col_idx = find(epochIndices == epoch_num);

    % Assign to GLM
    switch modality
        case 'sham'
            GLM_sham{subj_idx, col_idx} = DCM;
        case 'ctbs'
            GLM_ctbs{subj_idx, col_idx} = DCM;
        case 'itbs'
            GLM_itbs{subj_idx, col_idx} = DCM;
    end
end

% -----------------------------------------------------------------------
%            Check explained varience of data before proceeding
% -----------------------------------------------------------------------

threshold_explained_variance = 80; % Set the threshold

% Initialize variables to track explained variance filtering
low_explained_variance = false(numSubjects, numEpochs, 3); % 3 modalities: sham, itbs, ctbs

for iSub = 1:numSubjects
    for e = 1:numEpochs
        % Checking sham
        if ~isempty(GLM_sham{iSub, e})
            ev = explained_var_dcm(GLM_sham{iSub, e});
            if ev < threshold_explained_variance
                fprintf('Subject %d, sham epoch %d excluded due to low explained variance (%.2f%%)\n', uniqueSubjects(iSub), e+1, ev);
                GLM_sham{iSub, e} = [];
                low_explained_variance(iSub, e, 1) = true;
            end
        end
        
        % Checking itbs
        if ~isempty(GLM_itbs{iSub, e})
            ev = explained_var_dcm(GLM_itbs{iSub, e});
            if ev < threshold_explained_variance
                fprintf('Subject %d, itbs epoch %d excluded due to low explained variance (%.2f%%)\n', uniqueSubjects(iSub), e+1, ev);
                GLM_itbs{iSub, e} = [];
                low_explained_variance(iSub, e, 2) = true;
            end
        end
        
        % Checking ctbs
        if ~isempty(GLM_ctbs{iSub, e})
            ev = explained_var_dcm(GLM_ctbs{iSub, e});
            if ev < threshold_explained_variance
                fprintf('Subject %d, ctbs epoch %d excluded due to low explained variance (%.2f%%)\n', uniqueSubjects(iSub), e+1, ev);
                GLM_ctbs{iSub, e} = [];
                low_explained_variance(iSub, e, 3) = true;
            end
        end
    end
end

% Recompute missing data flags after filtering
missing_sham = any(cellfun('isempty', GLM_sham), 2);
missing_itbs = any(cellfun('isempty', GLM_itbs), 2);
missing_ctbs = any(cellfun('isempty', GLM_ctbs), 2);

missing_any = missing_sham | missing_itbs | missing_ctbs;

% Print out details of remaining missing subjects and epochs
for iSub = find(missing_any)'
    subjNum = uniqueSubjects(iSub);

    sham_missing_epochs = find(cellfun('isempty', GLM_sham(iSub,:)));
    for ep = sham_missing_epochs
        fprintf('Subject %d missing sham epoch %d\n', subjNum, ep+1);
    end

    itbs_missing_epochs = find(cellfun('isempty', GLM_itbs(iSub,:)));
    for ep = itbs_missing_epochs
        fprintf('Subject %d missing itbs epoch %d\n', subjNum, ep+1);
    end

    ctbs_missing_epochs = find(cellfun('isempty', GLM_ctbs(iSub,:)));
    for ep = ctbs_missing_epochs
        fprintf('Subject %d missing ctbs epoch %d\n', subjNum, ep+1);
    end
end

% Keep only subjects with complete data
GLM_sham_complete = GLM_sham(~missing_any, :);
GLM_itbs_complete = GLM_itbs(~missing_any, :);
GLM_ctbs_complete = GLM_ctbs(~missing_any, :);

numCompleteSubjects = size(GLM_sham_complete, 1);


% -----------------------------------------------------------------------
%                           Subject level PEB     
%                 Run PEB between groups for each epoch
% -----------------------------------------------------------------------
% itbs vs sham, and ctbs vs sham, separately, for each epoch.

fields = {'A','H','AN','T','CV'}; % Parameters to include in PEB

% Preallocate cell arrays to store results for each epoch
PEB_itbs_vs_sham_all = cell(numEpochs,1);
BMA_itbs_vs_sham_all = cell(numEpochs,1);
PEB_ctbs_vs_sham_all = cell(numEpochs,1);
BMA_ctbs_vs_sham_all = cell(numEpochs,1);

for e = 1:numEpochs
    % ---------------------------
    % itbs vs sham for epoch e
    % ---------------------------
    DCM_sham_epoch = GLM_sham_complete(:, e); 
    DCM_itbs_epoch = GLM_itbs_complete(:, e);

    % Combine into one GCM (just this epoch)
    GCM_itbs_vs_sham_epoch = [DCM_sham_epoch; DCM_itbs_epoch];

    nSham = length(DCM_sham_epoch);
    nItbs = length(DCM_itbs_epoch);

    % Design matrix: first column = intercept, second column = group difference
    M = struct();
    M.X = [ones(nSham+nItbs,1), [repmat(-1,[nSham,1]); repmat(+1,[nItbs,1])]];

    PEB_itbs_vs_sham_epoch = spm_dcm_peb(GCM_itbs_vs_sham_epoch, M, fields);
    BMA_itbs_vs_sham_epoch = spm_dcm_peb_bmc(PEB_itbs_vs_sham_epoch);

    % Store results in cell arrays
    PEB_itbs_vs_sham_all{e} = PEB_itbs_vs_sham_epoch;
    BMA_itbs_vs_sham_all{e} = BMA_itbs_vs_sham_epoch;

    
    fprintf('Results for itbs vs sham, Epoch %d:\n', epochIndices(e));
    %spm_dcm_peb_review(BMA_itbs_vs_sham_epoch, GCM_itbs_vs_sham_epoch, C);

    % ---------------------------
    % ctbs vs sham for epoch e
    % ---------------------------
    DCM_ctbs_epoch = GLM_ctbs_complete(:, e); 
    GCM_ctbs_vs_sham_epoch = [DCM_sham_epoch; DCM_ctbs_epoch];

    nCtbs = length(DCM_ctbs_epoch);

    M = struct();
    M.X = [ones(nSham+nCtbs,1), [repmat(-1,[nSham,1]); repmat(+1,[nCtbs,1])]];

    PEB_ctbs_vs_sham_epoch = spm_dcm_peb(GCM_ctbs_vs_sham_epoch, M, fields);
    BMA_ctbs_vs_sham_epoch = spm_dcm_peb_bmc(PEB_ctbs_vs_sham_epoch);

    % Store results for ctbs vs sham
    PEB_ctbs_vs_sham_all{e} = PEB_ctbs_vs_sham_epoch;
    BMA_ctbs_vs_sham_all{e} = BMA_ctbs_vs_sham_epoch;

    fprintf('Results for ctbs vs sham, Epoch %d:\n', epochIndices(e));
    % spm_dcm_peb_review(BMA_ctbs_vs_sham_epoch, GCM_ctbs_vs_sham_epoch, C);
end


%% Remove the intercept (i.e., the first column) for both itbs and ctbs PEBs
for e = 1:numEpochs
    
    % ================ PEB_itbs_vs_sham_all ================
    PEB = PEB_itbs_vs_sham_all{e};
    
    % Some epochs might be empty if no valid data, so check:
    if ~isempty(PEB)
        % 1) Keep only the second column in M.X (the group difference).
        %    Originally M.X was [N × 2]. Now it becomes [N × 1].
        PEB.M.X = PEB.M.X(:,2);

        % 2) Ep was [88 × 2]. Keep only the second column → [88 × 1].
        PEB.Ep = PEB.Ep(:,2);

        % 3) Cp was [176 × 176].  The block (89:176, 89:176) 
        %    corresponds to the second column’s parameters.
        PEB.Cp = PEB.Cp(89:176, 89:176);

        % 4) M.pE was [176 × 1].  Keep (89:176) → [88 × 1].
        PEB.M.pE = PEB.M.pE(89:176);

        % 5) M.pC was [176 × 176]. Keep the same block (89:176, 89:176).
        PEB.M.pC = PEB.M.pC(89:176, 89:176);

        % 6) (Optional but recommended) Fix Xnames to show only the second covariate
        %    If your original was something like {'Covariate 1','Covariate 2'},
        %    now you only want {'Covariate 2'}.
        if isfield(PEB,'Xnames') && numel(PEB.Xnames) >= 2
            PEB.Xnames = {PEB.Xnames{2}};
        end

        % 7) Store the updated PEB back into your cell array
        PEB_itbs_vs_sham_all{e} = PEB;
    end
    
    % ================ PEB_ctbs_vs_sham_all ================
    PEB = PEB_ctbs_vs_sham_all{e};
    
    if ~isempty(PEB)
        % 1) Keep only the second column
        PEB.M.X = PEB.M.X(:,2);

        % 2) Keep only the second column in Ep
        PEB.Ep = PEB.Ep(:,2);

        % 3) Prune Cp to keep the lower‐right block
        PEB.Cp = PEB.Cp(89:176, 89:176);

        % 4) pE
        PEB.M.pE = PEB.M.pE(89:176);

        % 5) pC
        PEB.M.pC = PEB.M.pC(89:176, 89:176);

        % 6) Fix Xnames
        if isfield(PEB,'Xnames') && numel(PEB.Xnames) >= 2
            PEB.Xnames = {PEB.Xnames{2}};
        end

        % 7) Store back
        PEB_ctbs_vs_sham_all{e} = PEB;
    end
end

%%
% -----------------------------------------------------------------------
%                           Group level PEB     
%     Second-Level PEB: Modeling Baseline/Transient/Sustained Patterns
% -----------------------------------------------------------------------

% We have PEB_itbs_vs_sham_all and PEB_ctbs_vs_sham_all from the first-level results.
% Each is a {numEpochs x 1} cell array of PEB structures, one per epoch.

% Define the patterns across the 4 epochs
pattern_baseline = [1;  1;  1;  1];   % no change over time
pattern_transient = [1; -1; -1;  1];  % transient change returns to baseline at the end
pattern_sustained = [1; -1; -1; -1];  % sustained change does not return to baseline

M = struct();
M.X = [pattern_baseline, pattern_transient, pattern_sustained];
fields = {'A','H','AN','T','CV'};

%--------------------------------------------------------------------------
% Second-level analysis: itbs vs sham
%--------------------------------------------------------------------------

% PEB_itbs_vs_sham_all is {4 x 1}, one entry per epoch.
% Each entry is a PEB structure from spm_dcm_peb at the first level.

% Perform a "PEB of PEBs"
PEB_itbs_second_level = spm_dcm_peb(PEB_itbs_vs_sham_all, M, fields);

% Perform Bayesian Model Averaging/Comparison at second-level
BMA_itbs_second_level = spm_dcm_peb_bmc(PEB_itbs_second_level);

fprintf('Second-level results (itbs vs sham): \n');
spm_dcm_peb_review(PEB_itbs_second_level);

% Pause before proceeding to the next analysis
input('Press Enter to continue to the ctbs vs sham analysis...');
%--------------------------------------------------------------------------
% Second-level analysis: ctbs vs sham
%--------------------------------------------------------------------------

% Repeat for ctbs vs sham
PEB_ctbs_second_level = spm_dcm_peb(PEB_ctbs_vs_sham_all, M, fields);
BMA_ctbs_second_level = spm_dcm_peb_bmc(PEB_ctbs_second_level);

% Review results 
fprintf('Second-level results (ctbs vs sham):\n');
spm_dcm_peb_review( PEB_ctbs_second_level);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE END :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Full path
full_path = fullfile('D:/signal_data/PEB_of_PEB_results', 'full_PEB_B0_removed.mat');

% Save the workspace variables
save(full_path);
