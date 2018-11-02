## Gradient Boosting Compared to Mendelian Models
## Running on Odyssey Cluster
## Using a non-carrier lifetime risk of 0.05 and carrier lifetime risks of 0.5
## Last updated: October 21, 2018

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
# nm.list <- c("MMR.ngc", "MMR", "MMR.025", "MMR.05", "MMR.2", "MMR.4", "XGB.mmr", "XGB.const")
types <- c(".ngc", "", ".025", ".05", ".2", ".4")
n.boot <- 100
M.mmr <- 35
M.const <- 60
covs <- can.names
shrink <- 0.1
bag <- 0.5


# res.gb <- setNames(vector("list", length(nm.list)), nm.list)
# res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
#                                                           c("EO", "AUC", "BS")), 
#                                           risk = setNames(data.frame(cbind(sim.gb$FamID, matrix(NA, nrow(sim.gb), n.boot))),
#                                                           c("FamID", paste("risk", 1:n.boot, sep = "")))))


start <- Sys.time()
res.gb <- gb.mmr(sim.gb, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")


# for(i in 1:n.boot){
#   print(i)
#   smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
#   train <- sim.gb[smp.train, ]
#   test <- sim.gb[-smp.train, ]
#   
#   ## MMRpro without GC
#   res.gb$MMR.ngc$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.ngc)
#   res.gb$MMR.ngc$risk[res.gb$MMR.ngc$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.ngc
#   
#   ## MMRpro with GC
#   res.gb$MMR$perf[i, ] <- perf.meas(test$MMR, test$P.MMR)
#   res.gb$MMR$risk[res.gb$MMR$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR
#   
#   res.gb$MMR.025$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.025)
#   res.gb$MMR.025$risk[res.gb$MMR.025$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.025
#   
#   res.gb$MMR.05$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.05)
#   res.gb$MMR.05$risk[res.gb$MMR.05$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.05
#   
#   res.gb$MMR.2$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.2)
#   res.gb$MMR.2$risk[res.gb$MMR.2$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.2
#   
#   res.gb$MMR.4$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.4)
#   res.gb$MMR.4$risk[res.gb$MMR.4$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.4
#   
#   ## XGBoost with MMRPRO
#   param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
#                 objective = "binary:logistic")
#   dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
#                             base_margin = logit(train$P.MMR.ngc))
#   dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
#                            base_margin = logit(test$P.MMR.ngc))
#   # watchlist <- list(train = dtrain.mmr, test = dtest.mmr)
#   # res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
#   #                          eval_metric = "auc", eval_metric = "rmse")
#   res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr)
#   # cv <- xgb.cv(params = param, data = dtrain.mmr, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
#   pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
#   res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
#   res.gb$XGB.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
#   
#   ## XGBoost with constant
#   dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
#   dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
#   # watchlist <- list(train = dtrain.const, test = dtest.const)
#   # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
#   #                            eval_metric = "auc", eval_metric = "rmse")
#   res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const)
#   # cv.const <- xgb.cv(params = param, data = dtrain.const, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
#   pred.xgb.const <- predict(res.xgb.const, dtest.const)
#   res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
#   res.gb$XGB.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
# }


save(res.gb, file = paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", a1, ".RData", sep = ""))



