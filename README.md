# Gradient-boosting

## Description

This respository contains code to apply gradient boosting to extend MMRpro prediction of carrying mutations of MLH1, MSH2, and MSH6.  All code is written by Theodore Huang except where specified.

## Running the code
In order to reproduce the simulation results, run "GradientBoostingSimulaton_lp.R" (for the low penetrance scenario) and "GradientBoostingSimulation" (for the high penetrance scenario). Note that the code is intended to be used with the Harvard FAS Odyssey cluster. The BayesMendel R package, which can be found at https://projects.iq.harvard.edu/bayesmendel/bayesmendel-r-package, is also required to run the code.
    
## The rest of files are all described below:

* Generating Families Functions subfolder: Functions to generate family history data
* MMRpro.gast.R: An extension to MMRpro that allows for adding gastric cancer information
* CheckFamStructure.gast.R, ImputeAge.gast.R, LyteSimple.gast.R: Helper functions for MMRpro.gast that allow for adding gastric cancer information
* Heterozygous Penetrance.R: Helper function that alters the penetrance so that homozygous and heterozygous carriers have the same risk
* Simulation Functions.R: Functions used in the simulations
* Gradient Boosting Functions.R: Functions used for gradient boosting
* pen_gb_gastric10.RData, pen_gb_gastric10_lp.RData, pen_gb_gastric10_lgc.RData: Objects containing the cancer penetrance, where the gastric cancer penetrance is scaled to further separate mutation carriers and non-carriers. "lp" refers to the low penetrance scenario, and "lgc" refers to the low gastric cancer scenario.
* GradientBoostingSimulation.R: Runs the simulations for the high penetrance scenario
* GradientBoostingSimulation_lp.R: Runs the simulations for the low penetrance scenario
* GradientBoostingSimulation_lp_sss.R: Runs the simulations for the low penetrance scenario with small sample size
* GradientBoostingSimulation_Transportability.R: Runs the simulations for the scenario where the training and testing sets are different
* GradientBoostingSimulation_lgc.R: Runs the simulations for the scenario where we have low gastric cancer penetrances and small sample size
* SimFamCharacteristics.R: Obtaining characteristics of the simulated families
* SimFamCharacteristics.RData: Characteristics of the simulated families
* GB_Simulation_Results.R: Combining all the simulation results
* Gradient Boosting Simulation Analysis.R: Summarizing the simulation results
* GB_Simulation_Results.RData: Simulation results
* GB_Comparisons.RData: Simulation results for the comparisons of each model to MMRpro without gastric cancer
