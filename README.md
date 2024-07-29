# Dynamic-Fluctuations-in-DMN-network-after-TBS

### Steps:

1. semi-automated ICA for artifact rejection with EEGLAB (B.A, F.A)
2. Focusing on 10 second before stimulation and three 10 secondes after stimulation (H.M)
3. t-map for each stimulation condition on theta, alpha and beta range (H.M)
4. Convert .set file to spm object (G.G)
5. add real time channel location and fids (G.G)
6. Source recounstruction with each individual MRI data (G.G)
7. Performing spectral DCM with Sliding Window Approach (Van de steen, 2019, 2021) (G.G)
8. Group level analysis with parametric empirical bayes (PEB) 



-------------------------------------------------------------------------------------------------
Reference

Van de Steen F, Almgren H, Razi A, Friston K, Marinazzo D. Dynamic causal modelling of fluctuating connectivity in resting-state EEG. Neuroimage. 2019 Apr 1;189:476-484. doi: 10.1016/j.neuroimage.2019.01.055. Epub 2019 Jan 26. PMID: 30690158; PMCID: PMC6435216.
