# TSSKEA
We propose a Two-Stage Sparse Knowledge-Guided Evolutionary Algorithm (TSSKEA) by incorporating network degree knowledge of molecular interaction network and sparse knowledge from the population’s historical and current non-dominated solutions into multiobjective multimodal optimization process.

This package includes Matlab scripts and several datasets for demo of TSSKEA approach:
‘TSSKEA.m’ is a Matlab function for the routine of experimental analysis. TSSKEA aims to identify personalized biomarkers (multi-modal PDNBs or PDNBs) and personalized edge-network biomarkers (multiple PDENBs contains multi-modal PDENBs and PDENBs) for detecting early warning signal of individual patients in cancer and find the drug target genes which can provide effective information for the early treatment of cancer by analyzing the biomarkers.

The input (case: BRCA) include: 
(1)Path: The path of the user where the 'Main,m' is located.

(2)PGIN_BRCA: BRCA patients' personalized gene interaction network construct by SSN method. i.‘BRCA_i_PGIN.mat’ indicates that the PGIN of the i-th BRCA patient which contains the subnetwork adjacency matrix and the name of the gene in the subnetwork.

(3)PEN_BRCA: BRCA patients' personalized edge network construct by iENA method. i.‘BRCA_i_PEN.mat’ indicates that the PEN of the i-th BRCA patient which contains the edge-network, gene pair, and the name of the gene in the edge-network.

The output results:
BRCA_result: Non-dominated solutions of patient samples obtained by TSSKEA.

i.‘BRCA_sample_i_TSSKEA_boxchart.mat’ stores non-dominated solutions for the i-th patient by running the evolutionary algorithm 30 times each time. 
ii.‘BRCA_sample_i_ TSSKEA PF.mat’ indicates that a group Pareto front solutions of i-th patient obtained by performing non-dominated sort on boxchart. 
iii.‘BRCA_sample_i TSSKEA _PS.mat’ indicates the PDNB\PDENB corresponding to Pareto front solutions. 
iv.‘BRCA_DNB_name.txt: Genes’ name of multi-modal PDNB or PDNB with the biggest score of BRCA patients. i.For example：i-th patient: patient’s ID. PDNB or multi-modal PDNB: the name of genes in PDNB or multi-modal PDNB；
‘BRCA_PDENB_name.txt: Gene pairs’ name of multi-modal PDENB or PDENB with the highest score of BRCA patients. For example：i-th patient: patient’s ID. PDENB or multi-modal PDENB: the name of genes in PDENB or multi-modal PDENB.

Suggestions (1)Hardware suggestions for running this package: Window 10 or above; Matlab 2021 or above; RAM 32G or above.

(2)When users analyzed running this package, please note that: 
i.Users should set the path in the program, firstly. 
ii.Parameter setting of Popnum, Max_CalNum,  Experiment_num, sLower and sUpper will affect the running time. With default parameters, multi-modal EA in TSSKEA takes about 30 minutes to identify multiple PDNBs\PDENBs for a BRCA patient. Users can decrease running time by modifying above parameter. 

