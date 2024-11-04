# Computational analysis of longitudinal electroencephalograms after theta burst stimulation over left DLPFC using hierarchical dynamic causal modelling.

## Steps:

### Data Preprocessing
1. using individual channel locs
2. Bad channel removal
3. ICA with runica
4. RV with DIPFIT
5. Reject components with less than 70% prob of being brain or RV more than 0.15
6.  Epoching data in 2 sec
7.  STFT applied to each epoch (1-50 hz)
8.  Epochs with a Z-score of more than three standard deviations were discarded

### Data Analysis

1. Loading data to SPM using fieldtrip conversion
2. Source localization: using individual MRI and real-time chanel locs
3. Model specification: spectral DCM, Conductance-based Canonical Microcircuit Model (cmm-NMDA)
4. Model estimation
5. Explained Variance: More than 95 is desirable
6. Model Selection: Using model variation and selecting the best model with free energy criteria
7. First Level Parametric empirical bayes
8. Peb of Peb: The influence of TMS protocols on connectivity parameters
9. Cross-Validation
-------------------------------------------------------------------------------------------------
IDS team 23:
Supervisor: Prof. Ali Motie Nasrabadi
Mentor: Armin Toghi
Members: Ghazale ghaffaripour, Hamed moghtaderi, Babak Aliyari

