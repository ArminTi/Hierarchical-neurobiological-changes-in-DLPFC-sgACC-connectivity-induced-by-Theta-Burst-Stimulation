%% batch_extract_BMA_data_for_plot.m  —  SEGMENT-PRESERVING VERSION
% ------------------------------------------------------------------------------
% Overview
% ------------------------------------------------------------------------------
% This script processes each `BMA_*.mat` file in the current directory and
% writes a companion `plot_<token>.mat` file containing a plotting-friendly
% structure in `Results.PEB_thresholded`.
%
% Key properties of the saved structure (segment-preserving contract):
%   • `Ep`     — length nEp vector of parameter estimates with all segments
%                concatenated in order (no reordering or deduplication).
%   • `Pp`     — length nEp vector of posterior probabilities for `Ep`.
%   • `Cp`     — length nEp vector of thresholded (retained) marginal variances.
%   • `Pnames` — cellstr of length nNames representing a single covariate block
%                (ONE block only; do not duplicate across segments).
%
% Segment logic expected by the Python plotter:
%     nSeg = numel(Ep) / numel(Pnames)
%
% Important design choices (by intent, do NOT change here):
%   • Does NOT expand `Pnames` to match `Ep` length.
%   • Does NOT collapse duplicate names or strip trailing "#k".
%   • Only applies probability thresholding (zeroing), preserving vector length.
%
% Typical workflow:
%   1) Place this script alongside your `BMA_*.mat` (and optionally `Peb_*.mat`).
%   2) Run the script (it will create `plot_<token>.mat` files).
%   3) Use the Python plotting code that splits segments via `nSeg` above.
%
% Dependencies:
%   - MATLAB
%   - SPM12 on the MATLAB path (tested with the 25.01.x family)
%     (Only `spm_Ncdf` is required from SPM in this script.)
%
% Notes:
%   - If a matching `Peb_*.mat` is found, the script also computes an optional
%     `Results.percentageCell` summary (1×nCovariates cell).
%   - Replace any machine-specific paths before committing this script.
% ------------------------------------------------------------------------------

%% 1) Setup SPM (edit local paths as needed; ensure SPM is on the MATLAB path)
addpath('C:/Users/Skibidi/Documents/Programs/spm_25.01.02/spm');  % <- edit if needed
if exist('spmPath','var') && exist(spmPath,'dir')
    addpath(spmPath);
end
if exist('spm','file')
    spm('defaults','eeg');
else
    warning('SPM not found on path; spm_Ncdf must still be resolvable.');
end

%% 2) Parameters
% Posterior probability threshold: parameters with Pp < threshold are set to 0 in Ep/Cp.
% T is the value at which the normal CDF is evaluated (standard approach in SPM).
threshold = 0.95;
T         = 0;

%% 3) Enumerate all BMA files
bmaList = dir('BMA_*.mat');
assert(~isempty(bmaList),'No files matching BMA_*.mat found.');

% Sort newest first (compatible with older MATLAB structures)
if isfield(bmaList,'datenum') && ~isempty([bmaList.datenum])
    [~,ord] = sort([bmaList.datenum], 'descend');
else
    dn = datenum({bmaList.date});
    [~,ord] = sort(dn, 'descend');
end
bmaList = bmaList(ord);

fprintf('Found %d BMA files. Processing each separately...\n', numel(bmaList));

for k = 1:numel(bmaList)
    bmaFile = bmaList(k).name;
    [~, stem] = fileparts(bmaFile);

    % Token extraction: BMA_<token>[_secondlevel] → <token>
    tok = regexp(stem, '^BMA_(.+?)(?:_secondlevel)?$', 'tokens', 'once');
    if ~isempty(tok), baseToken = tok{1}; else, baseToken = stem; end

    fprintf('\n[%d/%d] Processing %s  (token: %s)\n', k, numel(bmaList), bmaFile, baseToken);

    %% 4) Load PEB/BMA and compute Pp (robust to Cp being vector or matrix)
    PEBstruct = load(bmaFile);
    fns = fieldnames(PEBstruct);
    PEB = PEBstruct.(fns{1});

    % ---- Standardize shapes ---------------------------------------------------
    % `Ep` is treated as a column vector (all segments concatenated).
    Ep = PEB.Ep(:);

    % `Cp` may be an NxN covariance matrix or an Nx1 vector of marginal variances.
    Cp_all = PEB.Cp;
    if isvector(Cp_all)
        Cp_marginal = Cp_all(:);
    else
        [r,c] = size(Cp_all);
        if r==c
            Cp_marginal = diag(Cp_all);
        else
            error('Unexpected Cp size %dx%d in %s (neither vector nor square).', r, c, bmaFile);
        end
    end
    if numel(Cp_marginal) ~= numel(Ep)
        error('Cp length (%d) does not match Ep length (%d) in %s.', numel(Cp_marginal), numel(Ep), bmaFile);
    end

    % ---- Normalize Pnames to a single block (DO NOT expand) -------------------
    % Keep Pnames as one covariate block (length nNames). Segment count is inferred
    % from the ratio nEp/nNames; no duplication or "#k" collapsing is performed.
    Pnames = PEB.Pnames;
    if ischar(Pnames)
        Pnames = cellstr(Pnames);
    elseif ~iscell(Pnames)
        try
            Pnames = cellstr(string(Pnames));
        catch
            error('Pnames has unexpected type in %s.', bmaFile);
        end
    end
    Pnames = Pnames(:);            % column cellstr
    nEp    = numel(Ep);
    nNames = numel(Pnames);

    % Segment inference: require that Ep length is a multiple of Pnames length.
    if nNames == 0 || mod(nEp, nNames) ~= 0
        error('Cannot infer segments: numel(Ep)=%d not a multiple of numel(Pnames)=%d in %s.', nEp, nNames, bmaFile);
    end
    nSeg = nEp / nNames;
    fprintf('  Inferred segments: nSeg = %d (nEp=%d, nNames=%d)\n', nSeg, nEp, nNames);

    % Posterior probabilities under the normal model used by SPM
    Pp = 1 - spm_Ncdf(T, abs(Ep), Cp_marginal);

    % Threshold: zero out Ep/Cp where Pp < threshold (lengths are preserved)
    mask  = (Pp >= threshold);
    EpThr = Ep;           EpThr(~mask) = 0;
    CpThr = Cp_marginal;  CpThr(~mask) = 0;

    % Save exactly what the Python plotter expects (segment-preserving layout)
    PEB_thresholded = struct( ...
        'Ep',     EpThr, ...     % length nEp (concatenated segments)
        'Pp',     Pp, ...        % length nEp
        'Pnames', {Pnames}, ...  % length nNames (ONE block)
        'Cp',     CpThr ...      % length nEp
    );

    fprintf('  Thresholded: kept %d / %d params (Pp >= %.2f)\n', nnz(mask), numel(mask), threshold);

    %% 5) (Optional) Pair to Peb_*.mat and compute percentageCell
    % If a matching Peb file is found, compute per-covariate percentages of
    % positive/negative/zero across Pebs; stored in `Results.percentageCell`.
    percentageCell = [];
    pebFile = findPebForToken(baseToken);

    if ~isempty(pebFile)
        fprintf('  Using Peb file: %s\n', pebFile);
        S = load(pebFile);
        fn = fieldnames(S);
        Pebs = S.(fn{1});

        if iscell(Pebs)
            try
                Pebs = cell2mat(Pebs);
            catch
                Pebs = [];
            end
        end

        if ~isempty(Pebs) && isstruct(Pebs) && isfield(Pebs,'Ep')
            [nConnections, nCovariates] = size(Pebs(1).Ep);
            posCount   = zeros(nConnections, nCovariates);
            negCount   = zeros(nConnections, nCovariates);
            zeroCount  = zeros(nConnections, nCovariates);
            totalCount = zeros(nConnections, nCovariates);

            for iPeb = 1:numel(Pebs)
                EpM = Pebs(iPeb).Ep;
                [r,c] = size(EpM);
                assert(r==nConnections && c==nCovariates, ...
                    'Peb %d Ep size %dx%d != %dx%d.', iPeb, r, c, nConnections, nCovariates);
                posCount   = posCount   + double(EpM > 0);
                negCount   = negCount   + double(EpM < 0);
                zeroCount  = zeroCount  + double(EpM == 0);
                totalCount = totalCount + 1;
            end

            percentPositive = 100 * posCount ./ totalCount;
            percentNegative = 100 * negCount ./ totalCount;
            percentZero     = 100 * (totalCount - posCount - negCount) ./ totalCount;

            percentageCell = cell(1, nCovariates);
            for j = 1:nCovariates
                percentageCell{j} = [percentPositive(:,j), percentNegative(:,j), percentZero(:,j)];
            end
            fprintf('  Computed percentageCell from %d Pebs (1 x %d cells).\n', numel(Pebs), nCovariates);
        else
            warning('  Pebs invalid or missing Ep; skipping percentageCell.');
        end
    else
        warning('  No Peb_*.mat found for token "%s"; skipping percentageCell.', baseToken);
    end

    %% 6) Save one output per BMA
    % Primary output is `Results.PEB_thresholded` (segment-preserving).
    Results = struct('PEB_thresholded', PEB_thresholded);
    if ~isempty(percentageCell), Results.percentageCell = percentageCell; end

    outfile = sprintf('plot_%s.mat', baseToken);
    save(outfile, 'Results');
    fprintf('  Saved: %s\n', outfile);
end

fprintf('\nDone. Wrote %d result files.\n', numel(bmaList));

%% ------------ Local helpers ------------
function pebFile = findPebForToken(baseToken)
    % Return the best Peb_*.mat to pair with a given baseToken (case-insensitive).
    pebFile = '';
    allPeb = dir('Peb_*.mat');
    if isempty(allPeb), return; end

    tokSan = sanitize(baseToken);
    core   = sanitize(regexprep(baseToken, '^([^-_\s]+).*', '$1'));

    bestScore = -Inf; bestIdx = 0; newestTime = -Inf;
    for i = 1:numel(allPeb)
        nameSan = sanitize(allPeb(i).name);
        score = 0;
        if contains(nameSan, tokSan), score = score + 2; end
        if ~isempty(core) && contains(nameSan, core), score = score + 1; end

        t = getDatenum(allPeb(i));
        if score > bestScore || (score==bestScore && t>newestTime)
            bestScore = score; bestIdx = i; newestTime = t;
        end
    end
    if bestIdx>0, pebFile = allPeb(bestIdx).name; end
end

function s = sanitize(x)
    x = lower(x);
    s = regexprep(x, '[^a-z0-9]+', '');
end

function d = getDatenum(S)
    if isfield(S,'datenum') && ~isempty(S.datenum), d = S.datenum;
    else, d = datenum(S.date);
    end
end
