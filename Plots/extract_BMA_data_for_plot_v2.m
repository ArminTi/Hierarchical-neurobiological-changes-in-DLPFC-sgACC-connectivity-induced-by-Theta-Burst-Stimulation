%% batch_extract_BMA_data_for_plot.m
% Purpose:
%   Process each BMA_*.mat in the working folder and write a plotting-friendly
%   results file (plot_<token>.mat) containing thresholded parameters.
%
% Steps:
%   - Load PEB/BMA and compute posterior probabilities (Pp) from Ep and Cp.
%   - Normalize/expand Pnames across covariate repetitions without introducing "#k".
%   - Collapse any residual duplicates that differ only by a trailing "#k".
%   - (Optional) Pair to Peb_*.mat and compute percentageCell summaries.
%   - Save Results with the field `PEB_thresholded`.
%
% Requirements:
%   - MATLAB + SPM12 available on the MATLAB path (tested with 25.01.x family).
%
% Notes:
%   - The duplicate-collapse prevents multiple arrows for the same intrinsic
%     connection when BMA exports repeated covariate blocks.
%   - Replace any machine-specific paths before committing to a public repo.

%% 1) Setup SPM (edit local paths as needed; ensure SPM is on the MATLAB path)
addpath('C:/Users/Skibidi/Documents/Programs/spm_25.01.02/spm') % <-- replace with your local SPM path or manage paths externally
if exist(spmPath,'dir') % NOTE: spmPath must exist if used here; otherwise rely on the addpath() above
    addpath(spmPath);
    spm('defaults','eeg');
else
    warning('SPM path not found: %s', spmPath);
end

%% 2) Parameters
threshold = 0.95;   % Posterior probability threshold for retaining parameters
T         = 0;      % Value at which to evaluate spm_Ncdf

%% 3) Enumerate all BMA files
bmaList = dir('BMA_*.mat');
assert(~isempty(bmaList),'No files matching BMA_*.mat found.');

% Sort newest first (compatible with older MATLAB)
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

    % Extract base token: BMA_<token>[_secondlevel]
    tok = regexp(stem, '^BMA_(.+?)(?:_secondlevel)?$', 'tokens', 'once');
    if ~isempty(tok), baseToken = tok{1}; else, baseToken = stem; end

    fprintf('\n[%d/%d] Processing %s  (token: %s)\n', k, numel(bmaList), bmaFile, baseToken);

    %% 4) Load PEB/BMA and compute Pp (robust to Cp being vector or matrix)
    PEBstruct = load(bmaFile);
    fns = fieldnames(PEBstruct);
    PEB = PEBstruct.(fns{1});

    % ---- Standardize shapes ----
    Ep = PEB.Ep(:);                        % column vector
    Cp_all = PEB.Cp;
    % Cp may be NxN covariance OR Nx1 vector of marginal variances
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

    % ---- Normalize Pnames to match Ep length ----
    Pnames = PEB.Pnames;
    if ischar(Pnames)
        Pnames = cellstr(Pnames);
    end
    if iscell(Pnames)
        Pnames = Pnames(:);
    else
        try
            Pnames = cellstr(string(Pnames));
            Pnames = Pnames(:);
        catch
            error('Pnames has unexpected type in %s.', bmaFile);
        end
    end

    nEp    = numel(Ep);
    nNames = numel(Pnames);

    % If Ep contains repeated covariate blocks, expand Pnames by rewriting the
    % "Covariate <j>:" prefix for each repetition (no "#k" suffix is appended).
    if nNames > 0 && mod(nEp, nNames) == 0
        reps = nEp / nNames;
        if reps > 1
            base = Pnames;                 % base labels (e.g., "Covariate 1: H(...)")
            Pnames_rep = cell(nEp,1);
            for j = 1:reps
                idx = (j-1)*nNames + (1:nNames);
                block = base;
                % If there is an explicit "Covariate <num>:" prefix, rewrite it to j
                hasCov = ~cellfun('isempty', regexp(block, '^(?i)\s*covariate\s+\d+\s*:'));
                block(hasCov) = regexprep(block(hasCov), '^(?i)\s*covariate\s+\d+\s*:', sprintf('Covariate %d:', j));
                % Otherwise, prepend a covariate prefix
                block(~hasCov) = cellfun(@(s) sprintf('Covariate %d: %s', j, s), block(~hasCov), 'UniformOutput', false);
                Pnames_rep(idx) = block;
            end
            Pnames = Pnames_rep;
            nNames = numel(Pnames);
            assert(nNames == nEp, 'Internal error: expanded Pnames length mismatch.');
        end
    end

    % If still mismatched, fail fast (do not fabricate labels silently)
    if nNames ~= nEp
        error('Pnames count (%d) != Ep count (%d) in %s after normalization. Refusing to fabricate labels.', nNames, nEp, bmaFile);
    end

    % Posterior probability
    Pp = 1 - spm_Ncdf(T, abs(Ep), Cp_marginal);

    % Threshold (zero-out values below Pp threshold)
    mask  = (Pp >= threshold);
    EpThr = Ep;           EpThr(~mask) = 0;
    CpThr = Cp_marginal;  CpThr(~mask) = 0;

    % Collapse duplicates that only differ by a trailing "#k" (e.g., "...#1", "...#2").
    % This avoids duplicate edges per (covariate, ROI, src→dst) in downstream plots.
    % Adjust `combineEp` if you prefer sum/max/first-nonzero instead of mean.
    baseNames = regexprep(Pnames, '#\s*\d+\s*$', '');   % strip trailing #1/#2 if present
    [uniqNames, ~, grp] = unique(baseNames, 'stable');

    % Choose how to combine duplicates:
    combineEp  = @(x) mean(x);   % alternatives: sum / max / first-nonzero
    combinePp  = @(x) max(x);    % conservative: keep highest probability
    combineCp  = @(x) max(x);    % conservative: keep largest retained variance

    EpThr2 = accumarray(grp, EpThr, [], combineEp);
    Pp2    = accumarray(grp, Pp,    [], combinePp);
    CpThr2 = accumarray(grp, CpThr, [], combineCp);

    % Struct array, one per parameter (all columns same length now)
    PEB_thresholded = struct( ...
        'Ep',     num2cell(EpThr2(:)), ...
        'Pp',     num2cell(Pp2(:)), ...
        'Pnames', uniqNames(:), ...
        'Cp',     num2cell(CpThr2(:)) );

    % Report (raw count vs. unique after duplicate collapse)
    fprintf('  Params (raw): %d | unique after collapse: %d | pass >=%.2f: %d | zeroed: %d\n', ...
        numel(Pp), numel(Pp2), threshold, nnz(mask), nnz(~mask));

    %% 5) Find Peb file for this token (best-effort), compute percentageCell (optional)
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
    Results = struct('PEB_thresholded', PEB_thresholded);
    if ~isempty(percentageCell), Results.percentageCell = percentageCell; end

    outfile = sprintf('plot_%s.mat', baseToken);
    save(outfile, 'Results');
    fprintf('  Saved: %s\n', outfile);
end

fprintf('\nDone. Wrote %d result files.\n', numel(bmaList));

%% ------------ Local helpers ------------
function pebFile = findPebForToken(baseToken)
    % Return the best Peb_*.mat to pair with a given baseToken. Case-insensitive.
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
