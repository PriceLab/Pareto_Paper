# The geometry of clinical labs and wellness states from deeply phenotyped humans


In this paper we analyzed the clinical lab tests of 3094 individuals that participated in a wellness program (Arivale) using the ParTI package that was published in Hart et al., 2015.
The full Arivale deidentified dataset supporting the findings in this study, and necessary for this analysis, may be accessed for qualified researchers for research purposes upon request. Requests should be sent to Andrew Magis (andrew.magis@isbscience.org).

The article describing the ParTI approach can be found here:
http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3254.html

The documentation of the ParTI package can be found here:
http://www.weizmann.ac.il/mcb/UriAlon/download/ParTI

We downloaded the ParTI package and made the following adjustments, to analyze the data and create the figures:

The ParTI analysis work in two steps: 
1) Finding the archetype positions, calculating the significance of the simplex (P-value of 1000 times randomization process and t-ratio calculation), and calculating the error of the archetypes position. 
2) Enrichment analysis of given features (used to characterize the archetypes).

The first part of finding the archetype position and calculating the t-ratio is computationally intensive, depending on the machine, the size of the data and the two parameters: numIter, maxRuns, (it can take a few hours). 
In the second part, the enrichment analysis, we used 20 different matrixes, with ~13,000 variables. The ParTI analysis transforming the discrete variables into Boolean matrixes. Because of this process, and the huge amount of data, the enrichment analysis also took a lot of time, and many times stopped with errors. Therefore, we performed the analysis in steps: we first run the code to find the archetype position and calculating the P-value, the results were saved into log file, figures, and relevant parameters were saved for further analysis and the enrichment analysis. 

1) Finding the archetype positions, calculating the significance of the simplex (Figure 2)
Each participant in the program had between 1-8 visits. To determine the robustness of the tetrahedron for data selection (Figure 2), we randomly selected one timepoint per participant, and run the first part of the analysis to find the archetypes position and to calculate the significance of the tetrahedron (P-value). because of the randomization process, even running the analysis on the exact data selection - might yield slightly different results. The code that executes this process is in "rand_visit.m". We repeated this process 7 times, and picked the 7th time for further analysis (P-value<0.001, the variables were saved in 'Parti_vars.mat'). The exact position of the archetype and the P-value does not affect the enrichment analysis. 

2) Enrichment analysis (Figure 3, Table 1)
we performed the enrichment analysis also in steps, by calling the "run_part.m" script. This code uploads every matrix separately and the variables that were saved from the previous step, and run the enrichment analysis. 

To run both parts together, there's a need to unite all the matrixes into one, and then run "ParTI.m". 

corr_dist_analyte.m - creates Figure 4, and Figure S6, S7.

corr_pc_analyte.m - generates Table 2.

trajectories_analysis.m - generate the statistics and the results presented in the section titled "Utilizing longitudinal data and the movement on the tetrahedron for early detection of transitions from health to disease state", and Figure 6, and figure S9-S13. 

