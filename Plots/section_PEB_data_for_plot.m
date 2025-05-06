% combine_code_results.m
% Combines aggregated connection percentages and thresholded PEB/BMA results
% into a single MAT file.

%% 1. Setup SPM
addpath('C:/Users/Growth fire/Programs/Matlab plugins/spm_25.01.rc3/spm');
spm('defaults','eeg');

%% 2. Parameters
name      = 'cTBS';      % base name for files
threshold = 0.95;        % Pp threshold
T         = 0;           % value at which to evaluate spm_Ncdf

%% 3. Load and threshold PEB/BMA structure
bmaFile     = sprintf('BMA_%s_secondlevel.mat', name);
PEBstruct   = load(bmaFile);
fns         = fieldnames(PEBstruct);
PEB         = PEBstruct.(fns{1});

Ep          = PEB.Ep;           % parameter estimates
Cp_all      = PEB.Cp;           % full covariance matrix
Pnames      = PEB.Pnames;       % parameter names

Cp_marginal = diag(Cp_all);     % marginal variances
Pp          = 1 - spm_Ncdf(T, abs(Ep), Cp_marginal);

% Threshold: zero out Ep and Cp_marginal where Pp < threshold
for i = 1:numel(Pp)
    if Pp(i) < threshold
        Ep(i)           = 0;
        Cp_marginal(i)  = 0;
    end
end

%% 4. Convert numeric arrays to cell arrays
Ep_cell       = num2cell(Ep);
Pp_cell       = num2cell(Pp);
Cp_cell       = num2cell(Cp_marginal);
Pnames_cell   = Pnames;  % already a cell array

PEB_thresholded = struct( ...
    'Ep',     Ep_cell, ...
    'Pp',     Pp_cell, ...
    'Pnames', Pnames_cell, ...
    'Cp',     Cp_cell    ...
);

%% 5. Verify thresholding
idx_pass = find(Pp >= threshold);
idx_Ep   = find(cellfun(@(x) ~isequal(x,0), Ep_cell));
idx_Cp   = find(cellfun(@(x) ~isequal(x,0), Cp_cell));

if isequal(idx_pass, idx_Ep) && isequal(idx_pass, idx_Cp)
    fprintf('Check passed: nonzero Ep/Cp indices match for Pp >= %.2f\n', threshold);
else
    warning('Mismatch in nonzero Ep/Cp indices for Pp >= %.2f\n', threshold);
end

fprintf('Total parameters: %d\n', numel(Pp));
fprintf('Pp >= %.2f: %d parameters\n', threshold, numel(idx_pass));
fprintf('Pp < %.2f: %d parameters set to zero\n', threshold, sum(Pp < threshold));

%% 6. Combine with previously computed percentages and save
% Note: percentageCell must have been created earlier in this script.
Results = struct( ...
    'percentageCell',  percentageCell, ...
    'PEB_thresholded', PEB_thresholded ...
);

outfile = sprintf('%s_final_secondlevel_results.mat', name);
save(outfile, 'Results');
fprintf('Combined results saved to %s\n', outfile);
