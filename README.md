# sourdough-manuscript
This repository contains all the code used in the parameter estimation, generalized Lotka-Volterra (gLV) model predictions, and analysis and figure generation for the manuscript "Using Lotka-Volterra models to understand the dynamics of multi-species community assembly in sourdough". 

Species and Strain IDs used in the scripts:
1.	17B2 - Fructilactobacillus sanfranciscensis
2.	0092a - Levilactobacillus brevis
3.	232 - Lactiplantibacillus plantarum
4.	550 - Companilactobacillus paralimentarius
5.	460 - Acetobacter malorum
6.	253 - Saccharomyces cerevisiae
7.	163 - Wickerhamomyces anomalus
8.	228 - Kazakhstania humilis
9.	177 - Kazakhstania servazzii

Raw data information:
Experiment data: All original raw data from the growth curves are in the excel sheets singlegrowthrawdata.xlsx and pairwisegrowthrawdata.xlsx. This includes information about glycerol stocks used, input counts, timepoints and colonies counted. For analysis, the normalized CFU counts were saved in the CFU files GCobs_single_normalized.csv and GCobs_paired_normalized.csv. All results from the community assembly experiment are in the file communityassemblyexperiment.xlsx. For analysis, the data was copied into two CSV files, new_CFU_relative.csv and new_CFU.csv. The python scripts used for parameter estimation and R scripts used for analysis use the CSV files.

Model data: Parameters estimated from single and pairwise growth curves were saved in the files est_single_all_6-3.txt and est_all_pairs_6-3.txt. These were used in the multi-species gLV model predictions.
Growth curve data corresponding to the estimated parameters were also stored in the files allSingSpecTraj_5-3.csv and allPairSpecTraj_6-3.csv. These were used for supplementary figure 2. The file endpoints_new_OG_model_288.csv contains the predictions from the gLV model.

Scripts used:
Part A: python scripts used for parameter estimation
1.	singleModelRejSims.py – We used this script for randomly sampling values of ri and aii. We sampled 5*106 parameter combinations for single-species growth curves. We picked the initial values of Ni and Nj from the observed growth curves. We ran this script on the Tufts HPC Cluster as an array and then combined all the sampled growth curves into a single file using a bash script.
2.	rejSampSing.py – We used this script to obtain the parameters that best fit our observed data by minimizing the log-squared distance between the modeled and observed growth curves at all time points. Once again, we ran this on the Tufts HPC Cluster separately for each species and then concatanated it into a single file to generate the data in est_single_all_6-3.txt. 
3.	pairwiseModelRejSims.py - We used this script for randomly sampling values of aij given previously obtained values of ri and aii. We sampled 107 parameter combinations for pairwise growth curves. We picked the initial values of Ni and Nj from the observed growth curves. We used a similar protocol to generate these as the randomly sampled parameters from single species growth curves.
4.	rejSampPair.py - We used this script to obtain the parameters that best fit our observed data by minimizing the log-squared distance between the modeled and observed growth curves at all time points. Again, used a similar approach as above to generate the data in est_all_pairs_6-3.txt.

Part B: python scripts used for gLV model
1. model10_newPairwise.py - This is the main script we used to write the multi-species gLV model and generate species abundance predictions. Briefly, the script takes as input the growth parameters we estimated in the previous section, and randomly samples them to generate multi-species model predictions. TWe use the function solve_ivp to solve the 9-species differential equation. The other two scripts are small modifications made to this core script.
2. model10_newPairwise-LOO.py - Same as the previous script, except we now added an extra line of code to set each species' initial abundance to 0 and save the species abundances at the end.
3. model11_test_serial_transfer.py - Same as script 1, except this time the script has two models for every set of parameters - with and without bottlenecks. We added the bottlenecks simply by dividing the abundances at the end of a growth cycel by a dilution factor and using them as initial abundances for the new growth cycle. Also, we saved the entire growth curve data (growth over time) instead of just saving the species' abundances at the final timepoint.

Part C: R scripts used for analysis




