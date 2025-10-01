function [all_blinks, eye_movements, muscle, disc, sugg_comps,all_art] = tbs_rseeg_autoica(comps) % comps.topo: channels x IC
%% Parameters
FS          = comps.fsample;
THRES_B     = 2.5; % blink threshold
THRES_L     = 4; % lateral eye movement threshold 
THRES_M     = 0.5; % muscle activity threshold

% Frequencies
FREQ        = [0 150];
MUSCLE_FREQ = [31 100];

% Channels
FP1         = find(strcmp(comps.topolabel,'Fp1'));
FP2         = find(strcmp(comps.topolabel,'Fp2'));
F7          = find(strcmp(comps.topolabel,'F7'));
F8          = find(strcmp(comps.topolabel,'F8'));

%% PSD Calculations
% Needs to be sped up - timeseries is long!
for a = 1:size(comps.trial{1},1)
    [P,F]   = pwelch(comps.trial{1}(a,:),FS,FS/2,FS,FS);
    index   = find(F > FREQ(1) & F < FREQ(2));
    PSD(a,:)= P(index);
end
F = F(index);

%% Z-Scores
% Calculate TOPOLOGY z-scores
ztopo = zscore(comps.topo,0,1);

%% Blink ICs
% Average zscore of Fp1, and Fp2 should be greater than 2.5
blinkmean   = mean(ztopo([FP1,FP2],:));
blinks      = find(abs(blinkmean) > THRES_B);

%% Lateral Eye Movement
lateral     = find(abs(ztopo(F7,:) - ztopo(F8,:)) > THRES_L);

%% Muscle ICs
m_activity  = zeros(1,length(comps.topo));
i           = F >= MUSCLE_FREQ(1) & F <= MUSCLE_FREQ(2);

for m = 1:size(comps.topo,2)
    m_activity(m) = sum(PSD(m,i))/sum(PSD(m,:));
end

muscle      = find(m_activity > THRES_M);

%% FROM EEGLAB'S ADJUST
% comps.topo = channels x components

% Variables memorizing artifacted ICs indexes
blink   = [];
horiz   = [];
vert    = [];
disc    = [];

% Create icaact
cntdata      = horzcat(comps.trial{:});
comps.icaact = reshape(cntdata,[size(comps.topo,2),comps.fsample,length(cntdata)/comps.fsample]);

% Computes IC topographies
topografie = comps.topo'; 

% Normalise topographies
for i = 1:size(comps.topo,2) % number of ICs
    
    ScalingFactor   = norm(topografie(i,:));
    topografie(i,:) = topografie(i,:)/ScalingFactor;
    
    comps.icaact(i,:,:) = ScalingFactor*comps.icaact(i,:,:);
end

% Load EEGLAB-format easycapM1 channel position information
load('chanlocs.mat')

%% Feature extraction
disp(' ')
disp('Features Extraction:')

% GDSF - General Discontinuity Spatial Feature
disp('GDSF - General Discontinuity Spatial Feature...')
GDSF = compute_GD_feat(topografie,chanlocs,size(comps.topo,2));

% SED - Spatial Eye Difference
disp('SED - Spatial Eye Difference...')
[SED,medie_left,medie_right] = computeSED_NOnorm(topografie,...
    chanlocs,size(comps.topo,2));

% SAD - Spatial Average Difference
disp('SAD - Spatial Average Difference...')
[SAD,var_front,var_back,mean_front,mean_back] = computeSAD(topografie,...
    chanlocs,size(comps.topo,2));

% SVD - Spatial Variance Difference between front zone and back zone
diff_var = var_front - var_back;

% Epoch dynamic range, variance and kurtosis
num_epoch = size(comps.icaact,3);
K       = zeros(num_epoch,size(comps.topo,2)); % kurtosis
Vmax    = zeros(num_epoch,size(comps.topo,2)); % variance

disp('Computing variance and kurtosis of all epochs...')
for i = 1:size(comps.topo,2) % number of ICs
    for j = 1:num_epoch              
        Vmax(j,i)   = var(comps.icaact(i,:,j));        
        K(j,i)      = kurt(comps.icaact(i,:,j));
    end  
end

% TK - Temporal Kurtosis
disp('Temporal Kurtosis...')
meanK = zeros(1,size(comps.topo,2));

for i = 1:size(comps.topo,2)
    if num_epoch > 100
        meanK(1,i) = trim_and_mean(K(:,i)); 
    else meanK(1,i) = mean(K(:,i));
    end
end

% MEV - Maximum Epoch Variance
disp('Maximum epoch variance...')
maxvar  = zeros(1,size(comps.topo,2));
meanvar = zeros(1,size(comps.topo,2));

for i=1:size(comps.topo,2)
    if num_epoch > 100
        maxvar(1,i)     = trim_and_max(Vmax(:,i)');
        meanvar(1,i)    = trim_and_mean(Vmax(:,i)');
    else
        maxvar(1,i)     = max(Vmax(:,i));
        meanvar(1,i)    = mean(Vmax(:,i));
    end
end

% MEV in reviewed formulation:
nuovaV = maxvar./meanvar;

%% Thresholds computation
disp('Computing EM thresholds...')

[soglia_K,~,~]      = EM(meanK);
[soglia_SED,~,~]    = EM(SED);
[soglia_SAD,~,~]    = EM(SAD);
[soglia_GDSF,~,~]   = EM(GDSF);
[soglia_V,~,~]      = EM(nuovaV); 

%% Calculate artifacts
% Horizontal Eye Movements (HEM)
horiz = intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
    (find(nuovaV>=soglia_V)));

% Vertical eye movements (VEM)
vert = intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));

% Eye Blink (EB)
blink = intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
    intersect ( find(meanK>=soglia_K),find(diff_var>0) ));

% Generic Discontinuities (GD)
disc = intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));

% Combined artifacts
aic = unique([blink disc horiz vert]);

% Artifact ICs
art         = unique([blink,horiz,vert,disc]);
ic_limit    = round(size(comps.topo,2)*0.15);
if length(art) >= ic_limit
    sugg_comps  = art(1:ic_limit);
elseif length(art) < ic_limit
    sugg_comps  = art;
end

% All combined (within suggested limit)
all_art     = unique([blink,blinks,lateral,horiz,vert,muscle,disc]);

% Outputs
all_blinks      = unique([blinks,blink]);
eye_movements   = unique([lateral,horiz,vert]);

%% Store outputs
comps.hdr.blinks    = all_blinks;
comps.hdr.eyes      = eye_movements;
comps.hdr.muscle    = muscle;
comps.hdr.disc      = disc;
comps.hdr.suggested = sugg_comps;

%% Display output
disp(['Blink components are: ', num2str(all_blinks)]);
disp(['Eye movement components are: ', num2str(eye_movements)]);
disp(['Muscle components are: ', num2str(muscle)]);
disp(['Generic discontinuity components are: ', num2str(disc)]);
disp(['Suggest components for removal are: ', num2str(sugg_comps)]);

