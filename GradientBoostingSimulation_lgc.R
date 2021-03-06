## Gradient Boosting Compared to Mendelian Models
## Running on Odyssey Cluster
## Simulated data
## Using a non-carrier lifetime risk of 0.05 and carrier lifetime risks of 0.5
## Using low penetrance for gastric cancer
## Last updated: August 4, 2020

rm(list = ls())

a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

set.seed(a1)

library(BayesMendel)
library(abind)
library(xgboost)
library(pROC)
library(doParallel)
registerDoParallel(cores = 4)

source("MMRpro.gast.R")
source("ImputeAge.gast.R")
source("CheckFamStructure.gast.R")
source("Simulation Functions.R")
source("Gradient Boosting Functions.R")
source("Heterozygous Penetrance.R")
source("LyteSimple.gast.R")

source(paste(getwd(), "/Generating Families Functions/sim.simFam.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/genCancerPen.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.buildGenoMat.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.linkParents.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.simCurAgeVar.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.simCancerVars.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.buildBranchOfAlleleMats.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/helpers.R", sep = ""))


## Getting the penetrance
cancers <- c("ColorC", "EndomC", "GastC")
genes <- c("MLH1", "MSH2", "MSH6")

## using PanelPRO gastric cancer penetrance
load("pen_gb_gastric10_lgc.RData")



CP <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)

## using allele frequencies of 0.01
af <- setNames(rep(0.01, 3), genes)


## simulating data from the "true" penetrance
start <- Sys.time()
fam.sim <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3, CP, af)
}
print(difftime(Sys.time(), start, units = "secs"))



## Running MMRPRO with and without gastric cancer
## First two are MMRpro with data-generating penetrances (oracle)
## Second two are MMRpro with mispecified CRC and EC penetrances
## Last 4 are MMRpro with misspecified GC penetrances (along with misspecified CRC and EC penetrances)
start <- Sys.time()
mmr.simres <- foreach(i = 1:length(fam.sim), .combine = rbind) %dopar% {
  c(tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, gastric = FALSE, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = c(0.5, 0.5, 0.25)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = c(0.5, 0.5, 0.5)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = c(0.5, 0.5, 2)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = c(0.5, 0.5, 4)),
             error = function(e) rep(NA, 4)))
}
print(difftime(Sys.time(), start, units = "secs"))


mmr.simres <- data.frame(mmr.simres)
names(mmr.simres) <- c("P.MMR.ngc.ocl", "P.MLH1.ngc.ocl", "P.MSH2.ngc.ocl", "P.MSH6.ngc.ocl",
                       "P.MMR.ocl", "P.MLH1.ocl", "P.MSH2.ocl", "P.MSH6.ocl",
                       "P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc",
                       "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6",
                       "P.MMR.025", "P.MLH1.025", "P.MSH2.025", "P.MSH6.025",
                       "P.MMR.05", "P.MLH1.05", "P.MSH2.05", "P.MSH6.05",
                       "P.MMR.2", "P.MLH1.2", "P.MSH2.2", "P.MSH6.2",
                       "P.MMR.4", "P.MLH1.4", "P.MSH2.4", "P.MSH6.4")



sim.gb <- mmr.simres
for(i in 1:nrow(sim.gb)){
  sim.gb$MLH1[i] <- fam.sim[[i]]$MLH1[fam.sim[[i]]$isProband == 1]
  sim.gb$MSH2[i] <- fam.sim[[i]]$MSH2[fam.sim[[i]]$isProband == 1]
  sim.gb$MSH6[i] <- fam.sim[[i]]$MSH6[fam.sim[[i]]$isProband == 1]
  sim.gb$MMR[i] <- as.numeric(sum(c(sim.gb$MLH1[i], sim.gb$MSH2[i], sim.gb$MSH6[i])) > 0)
}


sim.gb$FamID <- 1:nrow(sim.gb)

can.names <- paste("Prop", cancers, sep = "")
sim.gb[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.gb)){
  if(!is.null(fam.sim[[i]])){
    sim.gb[c("PropColorC", "PropGastC")][sim.gb$FamID == i, ] <- colMeans(fam.sim[[i]][paste("isAff", cancers[c(1, 3)], sep = "")])
    sim.gb$PropEndomC[sim.gb$FamID == i] <- mean(fam.sim[[i]]$isAffEndomC[fam.sim[[i]]$Gender == 0])
  }
}
print(difftime(Sys.time(), start, units = "secs"))



## using gradient boosting, incorporating information on gastric cancer
types <- c(".ngc.ocl", ".ocl", ".ngc", "", ".025", ".05", ".2", ".4")
n.cv <- 100
covs <- list(can.names, can.names[1:2])
shrink <- 0.1
bag <- 0.5
M.mmr <- c(25, 50, 100)
M.const <- c(25, 50, 100)

start <- Sys.time()
res.gb <- gb.mmr(sim.gb, shrink, bag, M.mmr, M.const, covs, n.cv, types, seed = a1 + 999)

## small sample size (using first 1000 families)
res.gb.sss <- gb.mmr(sim.gb[1:1000, ], shrink, bag, M.mmr, M.const, covs, n.cv, types, seed = a1 + 999)
difftime(Sys.time(), start, units = "secs")


save(res.gb, res.gb.sss, file = paste0(getwd(), "/Gradient Boosting/Simulation/Low Gastric Cancer/simgb_lgc_", a1, ".RData"))



