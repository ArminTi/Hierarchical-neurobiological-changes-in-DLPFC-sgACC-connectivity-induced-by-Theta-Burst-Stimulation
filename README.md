# Hierarchical neurobiological changes in DLPFC–sgACC connectivity induced by Theta-Burst Stimulation

## Processing Pipline:
![image](https://github.com/user-attachments/assets/4af15c6e-61ee-47ba-b1f7-6d4d611ba84f)




This repository contains code that can be used to reproduce a Dynamic Causal Modelling (DCM) analysis and hierarchical Parametric Empirical Bayes (PEB) analysis for computing induced connectivity changes in DLPFC–sgACC brain regions after Theta Burst Stimulation (TMS) in 22 healthy participants.

## Requirements
- MATLAB R2023b
- SPM 25 (25.01.rc3)
- Python 3
- FieldTrip 


## Data
Available: https://doi.org/10.25452/figshare.plus.c.5910329

Citation: Moffa, A. H. et al. Neuromodulatory effects of theta burst stimulation to the prefrontal cortex. Sci Data 9, 717 (2022).


## Repository layout

- **processing_eeg_pipeline/**
  - `tbs_rseeg_pipeline.m` _ main pipeline script  
  - `functions/` _ helper functions used by the pipeline  

- **dcm_analysis/**
  - *I. DCM* _ model specification and estimation  
  - *II. PEB* _ hierarchical PEB comparison between modalities and timelines  
    - **Pre_Post/** _ scripts + output (before vs immediately after TMS)  
    - **Over_Time/** _ scripts + output (sustained, transient, no changes over time)  

- **plots/**
  - `extract_BMA_data_for_plot.m` _ extracts output parameters for plotting  
  - `connectivity_plot.m` _ plots macro, meso, and micro connections/dynamics  

