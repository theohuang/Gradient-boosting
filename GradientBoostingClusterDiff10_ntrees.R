## Gradient Boosting Compared to Mendelian Models
## Running on Odyssey Cluster
## Using a non-carrier lifetime risk of 0.05 and carrier lifetime risks of 0.5
## Finding optimal number of iterations
## Last updated: October 19, 2018

rm(list = ls())

a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

set.seed(a1 + 1e5)

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
  gen.fam(1e3, CP, af, censoring = TRUE)
}
print(difftime(Sys.time(), start, units = "secs"))


## function to run MMRPRO, including options to include gastric cancer
## and to scale the gastric cancer penetrance
mmr.gb <- function(fam, af, CP, gastric = TRUE, scl = 1, pwr = NULL){
  fam$AffectedColon <- fam$isAffColorC
  fam$AffectedEndometrium <- fam$isAffEndomC
  fam$AgeColon <- fam$AgeColorC
  fam$AgeEndometrium <- fam$AgeEndomC
  fam$Twins <- rep(0, nrow(fam))
  fam$Death <- fam$isDead
  fam$AgeDeath <- fam$CurAge
  if(gastric == TRUE){
    fam$AffectedGastric <- fam$isAffGastC
    fam$AgeGastric <- fam$AgeGastC
    return(unlist(MMRpro.gast(fam, counselee.id = fam$ID[fam$isProband == 1],
                              params = MMRparams(allef = list(c(1 - af[1],af[1]),
                                                              c(1 - af[2], af[2]),
                                                              c(1 - af[3], af[3])),
                                                 penetrance.net = hetpen(CP, gastric = TRUE, scl = scl,
                                                                         pwr = pwr)),
                              print = FALSE)@probs))
  } else{
    return(unlist(MMRpro(fam, counselee.id = fam$ID[fam$isProband == 1],
                         params = MMRparams(allef = list(c(1 - af[1],af[1]),
                                                         c(1 - af[2], af[2]),
                                                         c(1 - af[3], af[3])),
                                            penetrance.net = hetpen(CP, gastric = FALSE, scl = scl,
                                                                    pwr = pwr)),
                         print = FALSE)@probs))
  }
}


## Running MMRPRO with and without gastric cancer, with different scales
start <- Sys.time()
mmr.simres <- foreach(i = 1:length(fam.sim), .combine = rbind) %dopar% {
  c(tryCatch(mmr.gb(fam.sim[[i]], af = af, CP, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(fam.sim[[i]], af = af, CP),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(fam.sim[[i]], af = af, CP, pwr = 0.25),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(fam.sim[[i]], af = af, CP, pwr = 0.5),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(fam.sim[[i]], af = af, CP, pwr = 2),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(fam.sim[[i]], af = af, CP, pwr = 4),
             error = function(e) rep(NA, 4)))
}
print(difftime(Sys.time(), start, units = "secs"))



mmr.simres <- data.frame(mmr.simres)
names(mmr.simres) <- c("P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc")



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
for(i in 1:nrow(sim.gb)){
  if(!is.null(fam.sim[[i]])){
    sim.gb[c("PropColorC", "PropGastC")][sim.gb$FamID == i, ] <- colMeans(fam.sim[[i]][paste("isAff", cancers[c(1, 3)], sep = "")])
    sim.gb$PropEndomC[sim.gb$FamID == i] <- mean(fam.sim[[i]]$isAffEndomC[fam.sim[[i]]$Gender == 0])
  }
}

# save(fam.sim, mmr.simres, sim.gb, file = "gb_sim_test.RData")


## using gradient boosting, incorporating information on gastric cancer
nm.list <- c("MMR.ngc", "XGB.mmr", "XGB.const")
n.boot <- 100
M <- 300
covs <- c(can.names)
shrink <- 0.1
bag <- 0.5


######### evaluating using OE

res.gb <- setNames(vector("list", length(nm.list)), nm.list)
# res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
#                                                           c("EO", "AUC", "BS")),
#                                           eval = list(),
#                                           eval_best = setNames(data.frame(matrix(0, n.boot, 2)),
#                                                                c("AUC", "BS"))))
res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                          c("OE", "AUC", "BS")),
                                          eval = list(),
                                          eval_best = setNames(data.frame(matrix(0, n.boot, 1)),
                                                               "OE")))


start <- Sys.time()
for(i in 1:n.boot){
  print(i)
  smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
  train <- sim.gb[smp.train, ]
  test <- sim.gb[-smp.train, ]
  
  ## MMRpro without GC
  res.gb$MMR.ngc$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.ngc)
  # res.gb$MMR.ngc$risk[res.gb$MMR.ngc$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.ngc
  
  ## XGBoost with MMRPRO
  param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                objective = "binary:logistic")
  dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                            base_margin = logit(train$P.MMR.ngc))
  dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                           base_margin = logit(test$P.MMR.ngc))
  watchlist <- list(train = dtrain.mmr, test = dtest.mmr)
  # res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
  #                          early_stopping_rounds = M, eval_metric = "rmse",
  #                          eval_metric = "auc", verbose = 0)
  res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
                           early_stopping_rounds = M, eval_metric = oe.xgb,
                           maximize = FALSE, verbose = 0)
  pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
  res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
  # res.gb$XGB.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
  res.gb$XGB.mmr$eval[[i]] <- res.xgb.mmr$evaluation_log
  # res.gb$XGB.mmr$eval_best[i, ] <- c(which.max(res.xgb.mmr$evaluation_log$test_auc),
  #                                    which.min(res.xgb.mmr$evaluation_log$test_rmse))
  res.gb$XGB.mmr$eval_best[i, ] <- which.min(abs(res.xgb.mmr$evaluation_log$test_OE - 1))
  
  
  ## XGBoost with constant
  dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
  dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
  # watchlist <- list(train = dtrain.const, test = dtest.const)
  # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
  #                            eval_metric = "auc", eval_metric = "rmse")
  # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
  #                            early_stopping_rounds = M, eval_metric = "rmse",
  #                            eval_metric = "auc", verbose = 0)
  res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
                             early_stopping_rounds = M, eval_metric = oe.xgb,
                             maximize = FALSE, verbose = 0)
  pred.xgb.const <- predict(res.xgb.const, dtest.const)
  res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
  # res.gb$XGB.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
  res.gb$XGB.const$eval[[i]] <- res.xgb.const$evaluation_log
  # res.gb$XGB.const$eval_best[i, ] <- c(which.max(res.xgb.const$evaluation_log$test_auc),
  #                                      which.min(res.xgb.const$evaluation_log$test_rmse))
  res.gb$XGB.const$eval_best[i, ] <- which.min(abs(res.xgb.const$evaluation_log$test_OE - 1))
  
}


####### evaluating using AUC and BS

res.gb2 <- setNames(vector("list", length(nm.list)), nm.list)
res.gb2 <- lapply(res.gb2, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                          c("OE", "AUC", "BS")),
                                          eval = list(),
                                          eval_best = setNames(data.frame(matrix(0, n.boot, 2)),
                                                               c("AUC", "BS"))))


for(i in 1:n.boot){
  print(i)
  smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
  train <- sim.gb[smp.train, ]
  test <- sim.gb[-smp.train, ]
  
  ## MMRpro without GC
  res.gb2$MMR.ngc$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.ngc)
  # res.gb2$MMR.ngc$risk[res.gb2$MMR.ngc$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.ngc
  
  ## XGBoost with MMRPRO
  param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                objective = "binary:logistic")
  dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                            base_margin = logit(train$P.MMR.ngc))
  dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                           base_margin = logit(test$P.MMR.ngc))
  watchlist <- list(train = dtrain.mmr, test = dtest.mmr)
  res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
                           early_stopping_rounds = M, eval_metric = "rmse",
                           eval_metric = "auc", verbose = 0)
  # res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
  #                          early_stopping_rounds = M, eval_metric = eo.xgb,
  #                          maximize = FALSE, verbose = 0)
  pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
  res.gb2$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
  # res.gb2$XGB.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
  res.gb2$XGB.mmr$eval[[i]] <- res.xgb.mmr$evaluation_log
  res.gb2$XGB.mmr$eval_best[i, ] <- c(which.max(res.xgb.mmr$evaluation_log$test_auc),
                                     which.min(res.xgb.mmr$evaluation_log$test_rmse))
  # res.gb2$XGB.mmr$eval_best[i, ] <- which.min(res.xgb.mmr$evaluation_log$test_EO)
  
  
  ## XGBoost with constant
  dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
  dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
  # watchlist <- list(train = dtrain.const, test = dtest.const)
  # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
  #                            eval_metric = "auc", eval_metric = "rmse")
  res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
                             early_stopping_rounds = M, eval_metric = "rmse",
                             eval_metric = "auc", verbose = 0)
  # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
  #                            early_stopping_rounds = M, eval_metric = eo.xgb,
  #                            maximize = FALSE, verbose = 0)
  pred.xgb.const <- predict(res.xgb.const, dtest.const)
  res.gb2$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
  # res.gb2$XGB.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
  res.gb2$XGB.const$eval[[i]] <- res.xgb.const$evaluation_log
  res.gb2$XGB.const$eval_best[i, ] <- c(which.max(res.xgb.const$evaluation_log$test_auc),
                                       which.min(res.xgb.const$evaluation_log$test_rmse))
  # res.gb2$XGB.mmr$eval_best[i, ] <- which.min(res.xgb.mmr$evaluation_log$test_EO)
  
}
difftime(Sys.time(), start, units = "secs")

save(res.gb, res.gb2, file = paste(getwd(), "/Gradient Boosting/Diff10/NTrees/simgbpen_", a1, ".RData", sep = ""))



