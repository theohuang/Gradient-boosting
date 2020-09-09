## Gradient Boosting Compared to Mendelian Models
## Running on Odyssey Cluster
## Simulated data
## Using a non-carrier lifetime risk of 0.05 and carrier lifetime risks of 0.5
## Training and testing sets are different to evaluate transportability
## Last updated: July 29, 2020

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

## using scaled GC penetrance so that lifetime risk of carriers is 0.5 and
## lifetime risk of non-carrires is 0.05
load("pen_gb_gastric10.RData")

mod.pen <- function(pen, pwr){
  cdf <- 1 - (1 - cumsum(pen))^pwr
  return(c(cdf[1], diff(cdf)))
}

CP.train <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)
CP.test <- genCancerPen(genes, cancers,
                        lapply(penCancersF, function(x) apply(x, 2, mod.pen, pwr = 0.5)),
                        lapply(penCancersM, function(x) apply(x, 2, mod.pen, pwr = 0.5)),
                        maxK = length(genes), age.last = 95)

## using allele frequencies of 0.01
af <- setNames(rep(0.01, 3), genes)


## simulating data from the "true" penetrance
start <- Sys.time()
fam.sim.train <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3 / 2, CP.train, af)
}
fam.sim.test <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3 / 2, CP.test, af)
}
print(difftime(Sys.time(), start, units = "secs"))



## Running MMRPRO with and without gastric cancer
## First two are MMRpro with data-generating penetrances (oracle)
## Second two are MMRpro with mispecified CRC and EC penetrances
## Last 4 are MMRpro with misspecified GC penetrances (along with misspecified CRC and EC penetrances)
start <- Sys.time()
mmr.simres.train <- foreach(i = 1:length(fam.sim.train), .combine = rbind) %dopar% {
  c(tryCatch(mmr.gb.sim(fam.sim.train[[i]], af = af, CP.train, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.train[[i]], af = af, CP.train),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.train[[i]], af = af, CP.train, gastric = FALSE, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.train[[i]], af = af, CP.train, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)))
}
# still using the training penetrances when running MMRpro on the testing data
mmr.simres.test <- foreach(i = 1:length(fam.sim.test), .combine = rbind) %dopar% {
  c(tryCatch(mmr.gb.sim(fam.sim.test[[i]], af = af, CP.train, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.test[[i]], af = af, CP.train),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.test[[i]], af = af, CP.train, gastric = FALSE, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb.sim(fam.sim.test[[i]], af = af, CP.train, pwr = c(0.5, 0.5, 1)),
             error = function(e) rep(NA, 4)))
}
print(difftime(Sys.time(), start, units = "secs"))


mmr.simres.train <- data.frame(mmr.simres.train)
mmr.simres.test <- data.frame(mmr.simres.test)
names(mmr.simres.train) <- names(mmr.simres.test) <- c("P.MMR.ngc.ocl", "P.MLH1.ngc.ocl", "P.MSH2.ngc.ocl", "P.MSH6.ngc.ocl",
                                                       "P.MMR.ocl", "P.MLH1.ocl", "P.MSH2.ocl", "P.MSH6.ocl",
                                                       "P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc",
                                                       "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6")



sim.gb.train <- mmr.simres.train
sim.gb.test <- mmr.simres.test
for(i in 1:nrow(sim.gb.train)){
  sim.gb.train$MLH1[i] <- fam.sim.train[[i]]$MLH1[fam.sim.train[[i]]$isProband == 1]
  sim.gb.train$MSH2[i] <- fam.sim.train[[i]]$MSH2[fam.sim.train[[i]]$isProband == 1]
  sim.gb.train$MSH6[i] <- fam.sim.train[[i]]$MSH6[fam.sim.train[[i]]$isProband == 1]
  sim.gb.train$MMR[i] <- as.numeric(sum(c(sim.gb.train$MLH1[i], sim.gb.train$MSH2[i], sim.gb.train$MSH6[i])) > 0)
}
for(i in 1:nrow(sim.gb.test)){
  sim.gb.test$MLH1[i] <- fam.sim.test[[i]]$MLH1[fam.sim.test[[i]]$isProband == 1]
  sim.gb.test$MSH2[i] <- fam.sim.test[[i]]$MSH2[fam.sim.test[[i]]$isProband == 1]
  sim.gb.test$MSH6[i] <- fam.sim.test[[i]]$MSH6[fam.sim.test[[i]]$isProband == 1]
  sim.gb.test$MMR[i] <- as.numeric(sum(c(sim.gb.test$MLH1[i], sim.gb.test$MSH2[i], sim.gb.test$MSH6[i])) > 0)
}

sim.gb.train$FamID <- 1:nrow(sim.gb.train)
sim.gb.test$FamID <- 1:nrow(sim.gb.test)


can.names <- paste("Prop", cancers, sep = "")
sim.gb.train[can.names] <- NA
sim.gb.test[can.names] <- NA
for(i in 1:nrow(sim.gb.train)){
  if(!is.null(fam.sim.train[[i]])){
    sim.gb.train[c("PropColorC", "PropGastC")][sim.gb.train$FamID == i, ] <- colMeans(fam.sim.train[[i]][paste("isAff", cancers[c(1, 3)], sep = "")])
    sim.gb.train$PropEndomC[sim.gb.train$FamID == i] <- mean(fam.sim.train[[i]]$isAffEndomC[fam.sim.train[[i]]$Gender == 0])
  }
}
for(i in 1:nrow(sim.gb.test)){
  if(!is.null(fam.sim.test[[i]])){
    sim.gb.test[c("PropColorC", "PropGastC")][sim.gb.test$FamID == i, ] <- colMeans(fam.sim.test[[i]][paste("isAff", cancers[c(1, 3)], sep = "")])
    sim.gb.test$PropEndomC[sim.gb.test$FamID == i] <- mean(fam.sim.test[[i]]$isAffEndomC[fam.sim.test[[i]]$Gender == 0])
  }
}



## using gradient boosting, incorporating information on gastric cancer
types <- c(".ngc.ocl", ".ocl", ".ngc", "")
n.cv <- 1
covs <- list(can.names, can.names[1:2])
shrink <- 0.1
bag <- 0.5
M.mmr <- 50
M.const <- 50

start <- Sys.time()
res.gb <- gb.mmr(sim.gb.test, shrink, bag, M.mmr, M.const, covs, n.cv, types, boot = FALSE,
                 dat.train = sim.gb.train, dat.test = sim.gb.test, cv = FALSE, seed = a1 + 999)
difftime(Sys.time(), start, units = "secs")


save(res.gb, file = paste0(getwd(), "/Gradient Boosting/Simulation/Transportability/simgb_trans_", a1, ".RData"))



