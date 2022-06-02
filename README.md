# 2022WaterDepthAnalysisFinal 

- The code found in this repository was used to run the statistics in Lutek et al. 2022.

## Linear Statistics
### [run in R; found in the folder '1_RLinearStatistics']
*** Lutek_etal_2022_WaterDepth.R runs linear statistics for both kinematic and EMG variables.
*** CleanDataset_EMG.csv contains ELECTROMYOGRAPHY data averaged by trial.
*** CleanDataset_KINE.csv contains KINEMATIC data averaged by trial.


----- Polar Statistics -----
[run in Matlab; found in folder '2_MatlabPolarStatistics']
This code requires functions from the 'Circular Statistics Toolbox (Directional Statistics)' (Berens et al., 2009). It can be downloaded from the Matlab File Exchange here:
https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
All other functions are custom code written by K.L. and E.M.S.
Polar statistics for KINEMATIC variables can be run using the script "KINEVars.m". 
Polar statistics for ELECTROMYOGRAPHY variables can be run using the script "EMGVars.m".


REFERENCES
Berens, P. (2009). CircStat: A Matlab toolbox for circular statistics. Journal of Statistical Software. 31(10). 1-21.
