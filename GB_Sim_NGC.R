## Gradient boosting without gastric cancer information
## Last updated: May 2, 2019

rm(list = ls())

library(xgboost)
library(pROC)


source("Gradient Boosting Functions.R")


cancers <- c("ColorC", "EndomC", "GastC")

types <- vector()
n.boot <- 100
covs <- paste("Prop", cancers[1:2], sep = "")
shrink <- 0.1
bag <- 0.5



for(a1 in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Simulation/Low Penetrance/simgb_lp_", a1, ".RData", sep = ""))
  res.gb.50.ngc <- gb.mmr(sim.gb, shrink, bag, M.mmr = 50, M.const = 50, covs, n.boot, types, seed = 1)
  save(res.gb.50.ngc, file = paste(getwd(), "/Gradient Boosting/Simulation/Low Penetrance/No Gastric Cancer/simgb_lp_ngc_", a1, ".RData", sep = ""))
  load(paste(getwd(), "/Gradient Boosting/Simulation/simgb_", a1, ".RData", sep = ""))
  res.gb.50.ngc <- gb.mmr(sim.gb, shrink, bag, M.mmr = 50, M.const = 50, covs, n.boot, types, seed = 1)
  save(res.gb.50.ngc, file = paste(getwd(), "/Gradient Boosting/Simulation/No Gastric Cancer/simgb_ngc_", a1, ".RData", sep = ""))
}