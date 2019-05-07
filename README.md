# Gradient-boosting

## Description

This respository contains code to apply gradient boosting extend MMRpro prediction of carrying mutations of MLH1, MSH2, and MSH6.  All code is written by Theodore Huang except where specified.

## Running the code
In order to reproduce the simulation results, run "GradientBoostingClusterDiff10.R". Note that the code is intended to be used with the Harvard FAS cluster. The BayesMendel R package, which can be found at https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package, is also required to run the code.
    
## The rest of files are all described below:

* Generating Families Functions subfolder: functions to generate family history data
* MMRpro.gast.R: An extension to MMRpro that llows for adding gastric cancer information
* CheckFamStructure.gast.R, ImputeAge.gast.R, LyteSimple.gast.R: helper functions for MMRpro.gast that allow for adding gastric cancer information
* Heterozygous Penetrance.R: helper function that alters the penetrance so that homozygous and heterozygous carriers have the same risk
* Simulation Functions.R: functions used in the simulations
* Gradient Boosting Functions.R: functions used for gradient boosting
* GB_Sim_Diff10_Analysis.R: Simulation data analysis
* Diff10_Sim_Results.RData: Simulation data analysis results
