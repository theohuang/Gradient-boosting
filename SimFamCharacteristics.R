## Characteristics of simulated families
## Using scaled versions of CRC and EC penetrances
## Low CRC and EC penetrances
## Last updated: September 1, 2020

rm(list = ls())

library(BayesMendel)
library(abind)
library(xgboost)
library(pROC)
library(doParallel)
registerDoParallel(cores = 4)
library(dplyr)
library(data.table)

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


##Low penetrance
## using scaled GC penetrance so that lifetime risk of carriers is 0.5 and
## lifetime risk of non-carrires is 0.05
## using scaled CRC and EC penetrances so that the lifetime risk for
## carriers is 0.2
load("pen_gb_gastric10_lp.RData")
CP.lp <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)
af.lp <- setNames(rep(0.01, 3), genes)

## High penetrance
load("pen_gb_gastric10.RData")
CP.hp <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)
af.hp <- setNames(rep(0.01, 3), genes)


cn.lp <- paste0(c("fs", "male", "pro.crc", "pro.ec", "pro.gc", "pro.crc.age",
                "pro.ec.age", "pro.gc.age", "fam.crc", "fam.ec", "fam.gc",
                "pro.m1", "pro.m2", "pro.m6"), ".lp")
cn.hp <- paste0(c("fs", "male", "pro.crc", "pro.ec", "pro.gc", "pro.crc.age",
                  "pro.ec.age", "pro.gc.age", "fam.crc", "fam.ec", "fam.gc",
                  "pro.m1", "pro.m2", "pro.m6"), ".hp")
res.lp <- setNames(data.frame(matrix(NA, 100, length(cn.lp))), cn.lp)
res.hp <- setNames(data.frame(matrix(NA, 100, length(cn.hp))), cn.hp)


## simulating data from the "true" penetrance
start <- Sys.time()
for(j in 1:100){
  set.seed(j)
  print(j)
  
  ### Low penetrance ###
  fam.sim.lp <- foreach(i = 1:10, .combine = append) %dopar% {
    gen.fam(1e3, CP.lp, af.lp)
  }
  for(i in 1:length(fam.sim.lp)){
    fam.sim.lp[[i]]$FamID <- i
  }
  fam.pro <- rbindlist(lapply(fam.sim.lp, function(x) x[x$isProband == 1, ]))
  fam.rel <- rbindlist(lapply(fam.sim.lp, function(x) x[x$isProband != 1, ]))
  
  ## General
  res.lp$fs.lp[j] <- mean(unlist(lapply(fam.sim.lp, nrow)))
  res.lp$male.lp[j] <- nrow(filter(fam.pro, Gender == 1))
  
  ## Participant cancer history
  res.lp$pro.crc.lp[j] <- nrow(filter(fam.pro, isAffColorC == 1))
  res.lp$pro.crc.age.lp[j] <- mean(filter(fam.pro, isAffColorC == 1)$AgeColorC)
  res.lp$pro.ec.lp[j] <- nrow(filter(fam.pro, isAffEndomC == 1))
  res.lp$pro.ec.age.lp[j] <- mean(filter(fam.pro, isAffEndomC == 1)$AgeEndomC)
  res.lp$pro.gc.lp[j] <- nrow(filter(fam.pro, isAffGastC == 1))
  res.lp$pro.gc.age.lp[j] <- mean(filter(fam.pro, isAffGastC == 1)$AgeGastC)
  
  ## Family cancer history
  res.lp$fam.crc.lp[j] <- nrow(unique(fam.rel %>% filter(isAffColorC == 1) %>% select(FamID)))
  res.lp$fam.ec.lp[j] <- nrow(unique(fam.rel %>% filter(isAffEndomC == 1) %>% select(FamID)))
  res.lp$fam.gc.lp[j] <- nrow(unique(fam.rel %>% filter(isAffGastC == 1) %>% select(FamID)))
  
  ## Participant MMR mutations
  res.lp$pro.m1.lp[j] <- nrow(filter(fam.pro, MLH1 == 1))
  res.lp$pro.m2.lp[j] <- nrow(filter(fam.pro, MSH2 == 1))
  res.lp$pro.m6.lp[j] <- nrow(filter(fam.pro, MSH6 == 1))
  
  
  ### High penetrance ###
  fam.sim.hp <- foreach(i = 1:10, .combine = append) %dopar% {
    gen.fam(1e3, CP.hp, af.hp)
  }
  for(i in 1:length(fam.sim.hp)){
    fam.sim.hp[[i]]$FamID <- i
  }
  fam.pro <- rbindlist(lapply(fam.sim.hp, function(x) x[x$isProband == 1, ]))
  fam.rel <- rbindlist(lapply(fam.sim.hp, function(x) x[x$isProband != 1, ]))
  
  ## General
  res.hp$fs.hp[j] <- mean(unlist(lapply(fam.sim.hp, nrow)))
  res.hp$male.hp[j] <- nrow(filter(fam.pro, Gender == 1))
  
  ## Participant cancer history
  res.hp$pro.crc.hp[j] <- nrow(filter(fam.pro, isAffColorC == 1))
  res.hp$pro.crc.age.hp[j] <- mean(filter(fam.pro, isAffColorC == 1)$AgeColorC)
  res.hp$pro.ec.hp[j] <- nrow(filter(fam.pro, isAffEndomC == 1))
  res.hp$pro.ec.age.hp[j] <- mean(filter(fam.pro, isAffEndomC == 1)$AgeEndomC)
  res.hp$pro.gc.hp[j] <- nrow(filter(fam.pro, isAffGastC == 1))
  res.hp$pro.gc.age.hp[j] <- mean(filter(fam.pro, isAffGastC == 1)$AgeGastC)
  
  ## Family cancer history
  res.hp$fam.crc.hp[j] <- nrow(unique(fam.rel %>% filter(isAffColorC == 1) %>% select(FamID)))
  res.hp$fam.ec.hp[j] <- nrow(unique(fam.rel %>% filter(isAffEndomC == 1) %>% select(FamID)))
  res.hp$fam.gc.hp[j] <- nrow(unique(fam.rel %>% filter(isAffGastC == 1) %>% select(FamID)))
  
  ## Participant MMR mutations
  res.hp$pro.m1.hp[j] <- nrow(filter(fam.pro, MLH1 == 1))
  res.hp$pro.m2.hp[j] <- nrow(filter(fam.pro, MSH2 == 1))
  res.hp$pro.m6.hp[j] <- nrow(filter(fam.pro, MSH6 == 1))
}
print(difftime(Sys.time(), start, units = "secs"))


save(res.lp, res.hp, file = "SimFamCharacteristics.RData")