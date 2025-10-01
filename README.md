# Hierarchical neurobiological changes in DLPFC–sgACC connectivity induced by Theta-Burst Stimulation

## Processing Pipline:
![image](https://github.com/user-attachments/assets/4af15c6e-61ee-47ba-b1f7-6d4d611ba84f)




Reproducible TMS-EEG preprocessing and DCM/PEB analysis.

## Requirements
- MATLAB R2023b
- SPM 25 (25.01.rc3)
- Python 3
- FieldTrip 


## Data
Available: https://doi.org/10.25452/figshare.plus.c.5910329

Citation: Moffa, A. H. et al. Neuromodulatory effects of theta burst stimulation to the prefrontal cortex. Sci Data 9, 717 (2022).

## Repository layout
processing_eeg_pipeline  
|-  tbs_rseeg_pipeline.m          (main pipeline script)  
|-  functions/                     (helper functions used by the pipeline)

dcm_analysis  
|-  I. DCM (DCM model specification and estimation)  
|-  II. PEB (Hierarchical PEB comparison between modalities and time lines)  
     |-Pre_Post (script + output comparison of before and and immidiatly after TMS)   
     |-Over_Time (script + output of sustained, transient, no changes overtime between modalities)  

plots   
|- extract_BMA_data_for_plot (matlab code that extracts the output parameters for plotting)   
|- connectivity_plot (plots the macro meso and micro connections and dynamics)  

## Repository layout

- **processing_eeg_pipeline/**
  - `tbs_rseeg_pipeline.m` — main pipeline script  
  - `functions/` — helper functions used by the pipeline  

- **dcm_analysis/**
  - *I. DCM* — model specification and estimation  
  - *II. PEB* — hierarchical PEB comparison between modalities and timelines  
    - **Pre_Post/** — scripts + output (before vs immediately after TMS)  
    - **Over_Time/** — scripts + output (sustained, transient, no changes over time)  

- **plots/**
  - `extract_BMA_data_for_plot.m` — extracts output parameters for plotting  
  - `connectivity_plot.m` — plots macro, meso, and micro connections/dynamics  

