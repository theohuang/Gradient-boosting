# Gradient-boosting

## Description

This respository contains code to apply gradient boosting extend MMRpro prediction of carrying mutations of MLH1, MSH2, and MSH6.  All code is written by Theodore Huang except where specified.

## Running the code
In order to reproduce the simulation results, run "GradientBoostingClusterDiff10.R". Note that the code is intended to be used with the Harvard FAS cluster. The BayesMendel R package, which can be found at https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package, is also required to run the code.
    
## The rest of files are all described below:

* Generating Families Functions subfolder: Functions to generate family history data
* MMRpro.gast.R: An extension to MMRpro that allows for adding gastric cancer information
* CheckFamStructure.gast.R, ImputeAge.gast.R, LyteSimple.gast.R: Helper functions for MMRpro.gast that allow for adding gastric cancer information
* Heterozygous Penetrance.R: Helper function that alters the penetrance so that homozygous and heterozygous carriers have the same risk
* Simulation Functions.R: Functions used in the simulations
* Gradient Boosting Functions.R: Functions used for gradient boosting
* pen_gb_gastric10.RData: Object containing the cancer penetrance, where the gastric cancer penetrance is scaled to further separate mutation carriers and non-carriers
* pen_gb_gastric10.RData: Object containing the cancer penetrance, where the gastric cancer penetrance is scaled to further separate mutation carriers and non-carriers, and the colorectal and endometrial cancer penetrances are scaled to have lower penetrance
* Gradient Boosting Simulation Analysis.R: Simulation data analysis
* GB_Sim_Results.RData: Simulation data analysis results, where the data-generating colorectal and endometrial cancer penetrances are not scaled
* GB_Sim_Results_lp.RData: Simulation data analysis results, where the data-generating colorectal and endometrial cancer penetrances are not scaled to be have lower penetrance
* GB_Sim_NGC_Results.RData: Simulation data analysis results for the gradient boosting models, not using gastric cancer information as a feature
* gb_sim.job, gb_sim_lp.job, gb_sim_ngc.job: Files used to run the models on the Harvard FAS Odyssey cluster
