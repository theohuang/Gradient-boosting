## Gradient Boosting Compared to Mendelian Models
## Running on Odyssey Cluster
## Using a non-carrier lifetime risk of 0.05 and carrier lifetime risks of 0.5
## Last updated: April 4, 2019

rm(list = ls())

a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

set.seed(a1)

library(BayesMendel)
library(BMmultigene)
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


## Getting the penetrance
cancers <- c("ColorC", "EndomC", "GastC")
genes <- c("MLH1", "MSH2", "MSH6")

load("pen_gb_gastric10.RData")

CP <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)

## using allele frequencies of 0.01
af <- setNames(rep(0.01, 3), genes)


## simulating data from the "true" penetrance
start <- Sys.time()
fam.sim <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3, CP, af)
}
print(difftime(Sys.time(), start, units = "secs"))



## Running MMRPRO with and without gastric cancer, with different scales
mmr.simres <- foreach(i = 1:length(fam.sim), .combine = rbind) %dopar% {
  c(tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = 0.25),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = 0.5),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = 2),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim[[i]], af = af, CP, pwr = 4),
             error = function(e) rep(NA, 4)))
}



mmr.simres <- data.frame(mmr.simres)
names(mmr.simres) <- c("P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc",
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
types <- c(".ngc", "", ".025", ".05", ".2", ".4")
n.boot <- 100
M.mmr <- 50
M.const <- 50
covs <- can.names
shrink <- 0.1
bag <- 0.5


start <- Sys.time()
res.gb.50 <- gb.mmr(sim.gb, shrink, bag, M.mmr = 50, M.const = 50, covs, n.boot, types, seed = a1 + 999)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.gb.25 <- gb.mmr(sim.gb, shrink, bag, M.mmr = 25, M.const = 25, covs, n.boot, types, seed = a1 + 999)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.gb.100 <- gb.mmr(sim.gb, shrink, bag, M.mmr = 100, M.const = 100, covs, n.boot, types, seed = a1 + 999)
difftime(Sys.time(), start, units = "secs")


save(res.gb.50, res.gb.25, res.gb.100, file = paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", a1, ".RData", sep = ""))



