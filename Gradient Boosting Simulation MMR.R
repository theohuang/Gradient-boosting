## Simulation for MMR genes and cancers related to these genes
## Last updated: July 26, 2018


library(BayesMendel)
library(dplyr)
library(data.table)
library(BayesMendelKnowledgeBase)
library(BMmultigene)
library(pROC)
library(xtable)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(rpart)
library(caret)
library(e1071)
library(Metrics)
library(xgboost)
library(purrr)
library(Hmisc)
library(gbm)
library(doParallel)
registerDoParallel(cores = 4)

source('~/Dropbox (Partners HealthCare)/Gradient Boosting/Simulation Functions.R')
source('~/Dropbox (Partners HealthCare)/Gradient Boosting/Gradient Boosting Functions.R')

## Choosing the cancers, genes, and allele frequencies
genes <- c("MLH1", "MSH2", "MSH6")
# finding the cancers for which we have penetrances related to these genes
cancers <- unique(foreach(i = 1:length(genes), .combine = "c") %do% {
  setdiff(unlist(strsplit(ls(name = "package:BayesMendelKnowledgeBase",
                             pattern = paste(genes[i], ".*", sep = "")),
                          ".", fixed = TRUE)), genes[i])
})
# some penetrance objects may be NULL. If there is a cancer where we don't
# really have the penetrance, find it and delete it from the cancers list
can.rem <- vector()
for(i in 1:length(cancers)){
  ct <- 0
  for(j in 1:length(genes)){
    if(!exists(paste(genes[j], ".", cancers[i], sep = ""))){
      ct <- ct + 1
    } else{
      if(is.null(eval(parse(text = paste(genes[j], ".", cancers[i], sep = "")))$risk.table)){
        ct <- ct + 1
      }
    }
  }
  if(ct == length(genes)){
    can.rem <- c(can.rem, cancers[i])
  }
}
cancers <- setdiff(cancers, can.rem)

# af <- c(MMRparams()$allef[[1]][2], MMRparams()$allef[[2]][2], MMRparams()$allef[[3]][3])
# using an allele frequency of 0.01 for each gene to increase the number of positive carriers
af <- rep(0.01, 3)
names(af) <- genes

## Getting the penetrances
penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")
CP <- genCancerPen(mutations, cancers, penCancersF, penCancersM, maxK = length(genes), age.max = 94)

## increasing gastric, pancreatic, ovarian, and brain cancer penetrances to see if gradient boosting does better
# scaling them so that the lifetime risks are all 0.95
penCancersF2 <- gen.pen.list(cancers, genes, "Female")
penCancersM2 <- gen.pen.list(cancers, genes, "Male")
lt.risk <- 0.9
can.scale <- c("GastC", "PancC", "OC", "BrainC")
penCancersF2 <- pen.scale(penCancersF2, can.scale, lt.risk)
penCancersM2 <- pen.scale(penCancersM2, can.scale, lt.risk)
penCancersF2$ColorC <- penet.mmr.net
CP2 <- genCancerPen(mutations, cancers, penCancersF2, penCancersM2, maxK = length(genes), age.max = 94)


## using lifetime risk of 0.5
## using MMRPRO penetrances instead of KnowledgeBase ones for colorectal and endometrial cancers
penCancersF3 <- gen.pen.list(cancers, genes, "Female")
penCancersM3 <- gen.pen.list(cancers, genes, "Male")
lt.risk <- 0.5
can.scale <- c("GastC", "PancC", "OC", "BrainC")
penCancersF3 <- pen.scale(penCancersF3, can.scale, lt.risk)
penCancersM3 <- pen.scale(penCancersM3, can.scale, lt.risk)
penCancersF3$ColorC <- penet.mmr.net$fFX[, c("M000", "M100", "M010", "M001")]
penCancersM3$ColorC <- penet.mmr.net$fMX[, c("M000", "M100", "M010", "M001")]
penCancersF3$EndomC <- penet.mmr.net$fFY[, c("M000", "M100", "M010", "M001")]
CP3 <- genCancerPen(mutations, cancers, penCancersF3, penCancersM3, maxK = length(genes), age.max = 94)




# source('~/Dropbox (Partners HealthCare)/BM_multigene/R/helpers.R')



start <- Sys.time()
fam.sim <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e4, CP, af, seed = i)
}
print(difftime(Sys.time(), start, units = "secs"))



save(fam.sim, file = "gb.famsim.mmr01.RData")


start <- Sys.time()
fam.sim3 <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3, CP3, af, seed = i)
}
print(difftime(Sys.time(), start, units = "secs"))

save(fam.sim3, file = "fam.sim3.RData")


## using a different allele frequency
start <- Sys.time()
fam.sim4 <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3, CP3, af = setNames(rep(0.001, 3), genes), seed = i)
}
print(difftime(Sys.time(), start, units = "secs"))


mmr.gb <- function(fam, af){
  fam$AffectedColon <- fam$isAffColorC
  fam$AffectedEndometrium <- fam$isAffEndomC
  fam$AgeColon <- fam$AgeColorC
  fam$AgeEndometrium <- fam$AgeEndomC
  fam$Twins <- rep(0, nrow(fam))
  fam$Death <- fam$isDead
  fam$AgeDeath <- fam$CurAge
  return(unlist(MMRpro(fam, counselee.id = fam$ID[fam$isProband == 1],
                       params = MMRparams(allef = list(c(1 - af[1],af[1]),
                                                       c(1 - af[2], af[2]),
                                                       c(1 - af[3], af[3]))),
                       print = FALSE)@probs))
}

# 
# start <- Sys.time()
# sim.bm <- foreach(i = 1:5000, .combine = rbind) %dopar% {
#   tryCatch(mmr.gb(fam.sim[[i]], af = af), error = function(e) rep(NA, 4))
# }
# print(difftime(Sys.time(), start, units = "secs"))
# 
# 
# sim.bm <- data.frame(sim.bm)
# sim.bm$FamID <- 1:nrow(sim.bm)
# sim.bm <- sim.bm[, c(ncol(sim.bm), 1:(ncol(sim.bm) - 1))]
# names(sim.bm) <- c("FamID", "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6")
# 
# sim.bm[c("MMR", "MLH1", "MSH2", "MSH6")] <- NA
# start <- Sys.time()
# for(i in 1:nrow(sim.bm)){
#   if(!is.null(fam.sim[[i]])){
#     sim.bm$MLH1[i] <- fam.sim[[i]]$MLH1[fam.sim[[i]]$isProband == 1]
#     sim.bm$MSH2[i] <- fam.sim[[i]]$MSH2[fam.sim[[i]]$isProband == 1]
#     sim.bm$MSH6[i] <- fam.sim[[i]]$MSH6[fam.sim[[i]]$isProband == 1]
#     sim.bm$MMR[i] <- as.numeric(sum(c(sim.bm$MLH1[i], sim.bm$MSH2[i], sim.bm$MSH6[i])) > 0)
#   }
# }
# print(difftime(Sys.time(), start, units = "secs"))
# 
# sim.bm["ColorC"] <- NA
# start <- Sys.time()
# for(i in 1:nrow(sim.bm)){
#   if(!is.null(fam.sim[[i]])){
#     sim.bm$ColorC[i] <- fam.sim[[i]]$isAffColorC[fam.sim[[i]]$isProband == 1]
#   }
# }
# print(difftime(Sys.time(), start, units = "secs"))
# 
# # sim.bm <- filter(sim.bm, !is.na(P.BRCA))
# 
# can.names <- paste("Prop", cancers, sep = "")
# sim.bm[can.names] <- NA
# start <- Sys.time()
# for(i in 1:nrow(sim.bm)){
#   if(!is.null(fam.sim[[i]])){
#     sim.bm[can.names][sim.bm$FamID == i, ] <- colMeans(fam.sim[[i]][paste("isAff", cancers, sep = "")])
#   }
# }
# print(difftime(Sys.time(), start, units = "secs"))
# 
# can.names.num <- paste("Num", cancers, sep = "")
# sim.bm[can.names.num] <- NA
# start <- Sys.time()
# for(i in 1:nrow(sim.bm)){
#   if(!is.null(fam.sim[[i]])){
#     sim.bm[can.names.num][sim.bm$FamID == i, ] <- colSums(fam.sim[[i]][paste("isAff", cancers, sep = "")])
#   }
# }
# print(difftime(Sys.time(), start, units = "secs"))
# 
# 
# save(sim.bm, file = "gb.sim.mmr01.RData")
# 


### increasing penetrance

start <- Sys.time()
sim.bm2 <- foreach(i = 1:length(fam.sim2), .combine = rbind) %dopar% {
  tryCatch(mmr.gb(fam.sim2[[i]], af = af), error = function(e) rep(NA, 4))
}
print(difftime(Sys.time(), start, units = "secs"))


sim.bm2 <- data.frame(sim.bm2)
sim.bm2$FamID <- 1:nrow(sim.bm2)
sim.bm2 <- sim.bm2[, c(ncol(sim.bm2), 1:(ncol(sim.bm2) - 1))]
names(sim.bm2) <- Cs(FamID, P.MMR, P.MLH1, P.MSH2, P.MSH6)

sim.bm2[Cs(MMR, MLH1, MSH2, MSH6)] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2$MLH1[i] <- fam.sim2[[i]]$MLH1[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MSH2[i] <- fam.sim2[[i]]$MSH2[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MSH6[i] <- fam.sim2[[i]]$MSH6[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MMR[i] <- as.numeric(sum(c(sim.bm2$MLH1[i], sim.bm2$MSH2[i], sim.bm2$MSH6[i])) > 0)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

sim.bm2["ColorC"] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2$ColorC[i] <- fam.sim2[[i]]$isAffColorC[fam.sim2[[i]]$isProband == 1]
  }
}
print(difftime(Sys.time(), start, units = "secs"))

# sim.bm2 <- filter(sim.bm2, !is.na(P.BRCA))

can.names <- paste("Prop", cancers, sep = "")
sim.bm2[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2[can.names][sim.bm2$FamID == i, ] <- colMeans(fam.sim2[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

can.names.num <- paste("Num", cancers, sep = "")
sim.bm2[can.names.num] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2[can.names.num][sim.bm2$FamID == i, ] <- colSums(fam.sim2[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))




##### with MMRPRO penetrances

start <- Sys.time()
sim.bm3 <- foreach(i = 1:length(fam.sim3), .combine = rbind) %dopar% {
  tryCatch(mmr.gb(fam.sim3[[i]], af = af), error = function(e) rep(NA, 4))
}
print(difftime(Sys.time(), start, units = "secs"))

sim.bm3 <- data.frame(sim.bm3)
sim.bm3$FamID <- 1:nrow(sim.bm3)
sim.bm3 <- sim.bm3[, c(ncol(sim.bm3), 1:(ncol(sim.bm3) - 1))]
names(sim.bm3) <- Cs(FamID, P.MMR, P.MLH1, P.MSH2, P.MSH6)

sim.bm3[Cs(MMR, MLH1, MSH2, MSH6)] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm3)){
  if(!is.null(fam.sim3[[i]])){
    sim.bm3$MLH1[i] <- fam.sim3[[i]]$MLH1[fam.sim3[[i]]$isProband == 1]
    sim.bm3$MSH2[i] <- fam.sim3[[i]]$MSH2[fam.sim3[[i]]$isProband == 1]
    sim.bm3$MSH6[i] <- fam.sim3[[i]]$MSH6[fam.sim3[[i]]$isProband == 1]
    sim.bm3$MMR[i] <- as.numeric(sum(c(sim.bm3$MLH1[i], sim.bm3$MSH2[i], sim.bm3$MSH6[i])) > 0)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

save(sim.bm3, file = "sim.bm3.RData")

######


##### with new allele freuqency

start <- Sys.time()
sim.bm4 <- foreach(i = 1:length(fam.sim4), .combine = rbind) %dopar% {
  tryCatch(mmr.gb(fam.sim4[[i]], af = af), error = function(e) rep(NA, 4))
}
print(difftime(Sys.time(), start, units = "secs"))

sim.bm4 <- data.frame(sim.bm4)
sim.bm4$FamID <- 1:nrow(sim.bm4)
sim.bm4 <- sim.bm4[, c(ncol(sim.bm4), 1:(ncol(sim.bm4) - 1))]
names(sim.bm4) <- Cs(FamID, P.MMR, P.MLH1, P.MSH2, P.MSH6)

# for(i in 1:length(fam.sim4)){
#   names(fam.sim4[[i]])[which(names(fam.sim4[[i]]) %in% c("a", "b", "c"))] <-
#     genes
# }

sim.bm4[Cs(MMR, MLH1, MSH2, MSH6)] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm4)){
  if(!is.null(fam.sim4[[i]])){
    sim.bm4$MLH1[i] <- fam.sim4[[i]]$MLH1[fam.sim4[[i]]$isProband == 1]
    sim.bm4$MSH2[i] <- fam.sim4[[i]]$MSH2[fam.sim4[[i]]$isProband == 1]
    sim.bm4$MSH6[i] <- fam.sim4[[i]]$MSH6[fam.sim4[[i]]$isProband == 1]
    sim.bm4$MMR[i] <- as.numeric(sum(c(sim.bm4$MLH1[i], sim.bm4$MSH2[i], sim.bm4$MSH6[i])) > 0)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

save(sim.bm4, file = "sim.bm4.RData")

######





start <- Sys.time()
DT <- data.table(fam.sim[1:5000])
DT <- DT[, list(lapply(V1, vert, gene = "mmr"), lapply(V1, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm[Cs(NumVert, NumHorz)] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm$MSI <- 0
sim.bm$MSI[sim.bm$ColorC == 1 & sim.bm$MMR == 1] <-
  rbinom(length(which(sim.bm$ColorC == 1 & sim.bm$MMR == 1)), 1, 0.968)
sim.bm$MSI[sim.bm$ColorC == 1 & sim.bm$MMR == 0] <-
  rbinom(length(which(sim.bm$ColorC == 1 & sim.bm$MMR == 0)), 1, 1 - 0.914)



## increasing penetrance

start <- Sys.time()
DT <- data.table(fam.sim2)
DT <- DT[, list(lapply(fam.sim2, vert, gene = "mmr"), lapply(fam.sim2, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm2[Cs(NumVert, NumHorz)] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm2$MSI <- 0
sim.bm2$MSI[sim.bm2$ColorC == 1 & sim.bm2$MMR == 1] <-
  rbinom(length(which(sim.bm2$ColorC == 1 & sim.bm2$MMR == 1)), 1, 0.968)
sim.bm2$MSI[sim.bm2$ColorC == 1 & sim.bm2$MMR == 0] <-
  rbinom(length(which(sim.bm2$ColorC == 1 & sim.bm2$MMR == 0)), 1, 1 - 0.914)



### using MMRPRO penetrance

start <- Sys.time()
DT <- data.table(fam.sim3)
DT <- DT[, list(lapply(fam.sim3, vert, gene = "mmr"), lapply(fam.sim3, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm3[Cs(NumVert, NumHorz)] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm3$MSI <- 0
sim.bm3$MSI[sim.bm3$ColorC == 1 & sim.bm3$MMR == 1] <-
  rbinom(length(which(sim.bm3$ColorC == 1 & sim.bm3$MMR == 1)), 1, 0.968)
sim.bm3$MSI[sim.bm3$ColorC == 1 & sim.bm3$MMR == 0] <-
  rbinom(length(which(sim.bm3$ColorC == 1 & sim.bm3$MMR == 0)), 1, 1 - 0.914)


# can.names <- paste("Prop", cancers, sep = "")
sim.bm3[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm3)){
  if(!is.null(fam.sim3[[i]])){
    sim.bm3[can.names][sim.bm3$FamID == i, ] <- colMeans(fam.sim3[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

# can.names.num <- paste("Num", cancers, sep = "")
sim.bm3[can.names.num] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm3)){
  if(!is.null(fam.sim3[[i]])){
    sim.bm3[can.names.num][sim.bm3$FamID == i, ] <- colSums(fam.sim3[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))




##### different allele frequency

start <- Sys.time()
DT <- data.table(fam.sim4)
DT <- DT[, list(lapply(fam.sim4, vert, gene = "mmr"), lapply(fam.sim4, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm4[Cs(NumVert, NumHorz)] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm4$MSI <- 0
sim.bm4$MSI[sim.bm4$ColorC == 1 & sim.bm4$MMR == 1] <-
  rbinom(length(which(sim.bm4$ColorC == 1 & sim.bm4$MMR == 1)), 1, 0.968)
sim.bm4$MSI[sim.bm4$ColorC == 1 & sim.bm4$MMR == 0] <-
  rbinom(length(which(sim.bm4$ColorC == 1 & sim.bm4$MMR == 0)), 1, 1 - 0.914)


# can.names <- paste("Prop", cancers, sep = "")
sim.bm4[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm4)){
  if(!is.null(fam.sim4[[i]])){
    sim.bm4[can.names][sim.bm4$FamID == i, ] <- colMeans(fam.sim4[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

# can.names.num <- paste("Num", cancers, sep = "")
sim.bm4[can.names.num] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm4)){
  if(!is.null(fam.sim4[[i]])){
    sim.bm4[can.names.num][sim.bm4$FamID == i, ] <- colSums(fam.sim4[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

save(sim.bm4, file = "sim.bm4.RData")

# sim.gb <- sim.bm4
sim.gb <- sim.bm3
# sim.gb <- sim.bm2
# sim.gb$LO.0 <- log(sim.gb$P.MMR / (1 - sim.gb$P.MMR))


save(sim.gb, file = "sim.gb.incpen.RData")


# covs <- c(can.names[-c(2:3)], "NumVert", "NumHorz")
# # covs <- "."
# train <- sim.gb[1:5000, ]
# test <- sim.gb[-(1:5000), ]
# res <- gbm.r(sim.gb, covs, M = 100, shrink = 0.1, bag = 0.5, train)



# 
# ##### XGBoost ######
# 
# dtrain <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
# dtest <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
# # watchlist <- list(train = dtrain, eval = dtest)
# 
# param <- list(max_depth = 2, eta = 0.01, silent = 1, nthread = 10, subsample = 0.5,
#               objective = "binary:logistic")
# res.xgb <- xgb.train(param, dtrain, nrounds = 500)
# pred.xgb <- predict(res.xgb, dtest)
# 
# sum(pred.xgb) / sum(test$MMR)
# pROC::auc(test$MMR ~ pred.xgb)
# mean((test$MMR - pred.xgb)^2)
# 
# sum(pred.xgb) / sum(train$MMR)
# pROC::auc(train$MMR ~ pred.xgb)
# mean((train$MMR - pred.xgb)^2)
# 
# 
# 
# ########


# df <- data.frame(matrix(0, 100, 12))
# names(df) <- c("EO.GB", "EO.GBM", "EO.XGB", "EO.MMR",
#                "AUC.GB", "AUC.GBM", "AUC.XGB", "AUC.MMR",
#                "BS.GB", "BS.GBM", "BS.XGB", "BS.MMR")
nm.list <- Cs(MMR, linreg.mmr, linreg.const, tree.mmr, tree.const, GBM, XGB.mmr, XGB.const)
n.boot <- 100
M <- 100
covs <- c(can.names, can.names.num, "NumVert", "NumHorz")
shrink <- 0.1
bag <- 0.5


res.gb <- setNames(vector("list", length(nm.list)), nm.list)
res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                          Cs(EO, AUC, BS)), 
                                          risk = setNames(data.frame(cbind(sim.gb$FamID, matrix(NA, nrow(sim.gb), n.boot))),
                                                          c("FamID", paste("risk", 1:n.boot, sep = "")))))
set.seed(1)
start <- Sys.time()
for(i in 1:n.boot){
  print(i)
  smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
  train <- sim.gb[smp.train, ]
  test <- sim.gb[-smp.train, ]
  
  # ## GB using linear regression and initializing with MMRPRO results
  # res.linreg.mmr <- gb.mmr(sim.gb, covs, M, shrink, bag, train, test,
  #                          init = "MMR", method = "linreg", seed = i)
  # res.gb$MMR$perf[i, ] <- res.linreg.mmr$perf.test.mmr
  # res.gb$linreg.mmr$perf[i, ] <- res.linreg.mmr$perf.test
  # res.gb$MMR$risk[res.gb$MMR$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR
  # res.gb$linreg.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- res.linreg.mmr$p.test
  # 
  # ## GB using linear regression and initializing with constant
  # res.linreg.const <- gb.mmr(sim.gb, covs, M, shrink, bag, train, test,
  #                            init = "const", method = "linreg", seed = i)
  # res.gb$linreg.const$perf[i, ] <- res.linreg.const$perf.test
  # res.gb$linreg.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- res.linreg.const$p.test
  # 
  # ## GB using trees and initializing with MMRPRO results
  # res.tree.mmr <- gb.mmr(sim.gb, covs, M, shrink, bag, train, test,
  #                        init = "MMR", method = "tree", seed = i)
  # res.gb$tree.mmr$perf[i, ] <- res.tree.mmr$perf.test
  # res.gb$tree.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- res.tree.mmr$p.test
  # 
  # ## GB using trees and initializing with constant
  # res.tree.const <- gb.mmr(sim.gb, covs, M, shrink, bag, train, test,
  #                          init = "const", method = "tree", seed = i)
  # res.gb$tree.const$perf[i, ] <- res.tree.const$perf.test
  # res.gb$tree.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- res.tree.const$p.test
  # 
  # ## gradient boosting using the gbm R package
  # res.gbm <- gbm(as.formula(paste("MMR ~", paste(covs,
  #                                                collapse = " + "))),
  #                data = train, interaction.depth = 2, shrinkage = shrink,
  #                n.trees = M, distribution = "bernoulli")
  # gbm.pred <- predict(res.gbm, test, M)
  # rp.gbm <- 1 / (1 + exp(-gbm.pred))
  # res.gb$GBM$perf[i, ] <- perf.meas(test$MMR, rp.gbm)
  # res.gb$GBM$risk[sim.gb$FamID %in% test$FamID, i + 1] <- rp.gbm
  # 
  
  ## XGBoost with MMRPRO
  param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                objective = "binary:logistic")
  dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                        base_margin = logit(train$P.MMR))
  dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                       base_margin = logit(test$P.MMR))
  res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M)
  pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
  res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
  res.gb$XGB.mmr$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
  
  ## XGBoost with constant
  dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
  dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
  res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M)
  pred.xgb.const <- predict(res.xgb.const, dtest.const)
  res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
  res.gb$XGB.const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
}
difftime(Sys.time(), start, units = "secs")

save(res.gb, file = "res.gb073118.RData")


res.xgb <- setNames(vector("list", 2), Cs(MMR, const))
res.xgb <- lapply(res.xgb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                          Cs(EO, AUC, BS)), 
                                          risk = setNames(data.frame(cbind(sim.gb$FamID, matrix(NA, nrow(sim.gb), n.boot))),
                                                          c("FamID", paste("risk", 1:n.boot, sep = "")))))


M <- 100
for(i in 1:n.boot){
  param <- list(max_depth = 3, eta = shrink, subsample = bag,
                objective = "binary:logistic", seed = i)
  print(i)
  smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
  train <- sim.gb[smp.train, ]
  test <- sim.gb[-smp.train, ]

  dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                        base_margin = logit(train$P.MMR))
  dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
  dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                           base_margin = logit(test$P.MMR))
  dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
  
  ## constant
  bst <- xgb.train(param, dtrain.const, nrounds = M)
  pred.xgb <- predict(bst, dtest.const)
  res.xgb$const$perf[i, ] <- perf.meas(test$MMR, pred.xgb)
  res.xgb$const$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb
  
  ## MMRPRO
  # setinfo(dtrain, "base_margin", logit(train$P.MMR))
  bst <- xgb.train(param, dtrain.mmr, nrounds = M)
  pred.xgb <- predict(bst, dtest.mmr)
  res.xgb$MMR$perf[i, ] <- perf.meas(test$MMR, pred.xgb)
  res.xgb$MMR$risk[sim.gb$FamID %in% test$FamID, i + 1] <- pred.xgb
}

lapply(res.xgb, function(x) summary(x$perf))
lapply(res.xgb, function(x) apply(x$perf, 2, function(y) setNames(c(sort(y)[5], mean(y), sort(y)[95]),
                                                                 Cs(CI.L, Mean, CI.U))))

mean.xgb <- setNames(vector("list", length(nm.list)), nm.list)
mean.xgb <- lapply(res.xgb, function(x) rowMeans(x$risk[, -1], na.rm = TRUE))




mean.gb <- setNames(vector("list", length(nm.list)), nm.list)
mean.gb <- lapply(res.gb, function(x) rowMeans(x$risk[, -1], na.rm = TRUE))


lapply(res.gb, function(x) summary(x$perf))

lapply(res.gb, function(x) apply(x$perf, 2, function(y) setNames(c(sort(y)[5], mean(y), sort(y)[95]),
                                                                 Cs(CI.L, Mean, CI.U))))

lapply(res.gb000, function(x) apply(x$perf, 2, function(y) setNames(c(sort(y)[5], mean(y), sort(y)[95]),
                                                                 Cs(CI.L, Mean, CI.U))))




save(res.gb, file = "res.gb.pen.mmrpro.RData")


res.gb.perf <- setNames(data.frame(matrix(unlist(lapply(res.gb, function(x) colMeans(x$perf))), nrow = 7, byrow = TRUE)), Cs(EO, AUC, BS))
rownames(res.gb.perf) <- nm.list
xtable(res.gb.perf, digits = 4)

grid.arrange(
  ggplot(data.frame(EO = c(res.tree.mmr$perf.train$EO,
                           res.tree.const$perf.train$EO),
                    Type = rep(1:2, each = nrow(res.tree.mmr$perf.train)),
                    Iteration = rep(1:(M + 1), 2)),
         aes(Iteration, EO)) +
    geom_point(aes(color = as.factor(Type))) +
    scale_color_discrete(name = "Initial Prediction", labels = Cs(MMRPRO, Constant)),
  ggplot(data.frame(AUC = c(res.tree.mmr$perf.train$AUC,
                            res.tree.const$perf.train$AUC),
                    Type = rep(1:2, each = nrow(res.tree.mmr$perf.train)),
                    Iteration = rep(1:(M + 1), 2)),
         aes(Iteration, AUC)) +
    geom_point(aes(color = as.factor(Type))) +
    scale_color_discrete(name = "Initial Prediction", labels = Cs(MMRPRO, Constant)),
  ggplot(data.frame(BS = c(res.tree.mmr$perf.train$BS,
                           res.tree.const$perf.train$BS),
                    Type = rep(1:2, each = nrow(res.tree.mmr$perf.train)),
                    Iteration = rep(1:(M + 1), 2)),
         aes(Iteration, BS)) +
    geom_point(aes(color = as.factor(Type))) +
    scale_color_discrete(name = "Initial Prediction", labels = Cs(MMRPRO, Constant))
)


pl.pred <- vector("list", 20)
ind.pos <- which(train$MMR == 1)[1:10]
ind.neg <- which(train$MMR == 0)[1:10]
for(i in 1:10){
  cl <- Cs(blue, red)[train$MMR[ind.neg[i]] + 1]
  pl.pred[[i]] <- ggplot(data.frame(Prediction = expit(res.tree.mmr$train.lo[ind.neg[i], ]),
                                    Iteration = 0:M),
                         aes(Iteration, Prediction)) +
    geom_point(color = cl) + labs(x = "", y = "")
}
for(i in 1:10){
  cl <- Cs(blue, red)[train$MMR[ind.pos[i]] + 1]
  pl.pred[[i + 10]] <- ggplot(data.frame(Prediction = expit(res.tree.mmr$train.lo[ind.pos[i], ]),
                                    Iteration = 0:M),
                         aes(Iteration, Prediction)) +
    geom_point(color = cl) + labs(x = "", y = "")
}

do.call("grid.arrange", c(pl.pred, ncol= 5))


pl.pred <- vector("list", 20)
ind.pos <- which(train$MMR == 1)[1:10]
ind.neg <- which(train$MMR == 0)[1:10]
for(i in 1:10){
  cl <- Cs(blue, red)[train$MMR[ind.neg[i]] + 1]
  pl.pred[[i]] <- ggplot(data.frame(Prediction = expit(res.tree.const$train.lo[ind.neg[i], ]),
                                    Iteration = 0:M),
                         aes(Iteration, Prediction)) +
    geom_point(color = cl) + labs(x = "", y = "")
}
for(i in 1:10){
  cl <- Cs(blue, red)[train$MMR[ind.pos[i]] + 1]
  pl.pred[[i + 10]] <- ggplot(data.frame(Prediction = expit(res.tree.const$train.lo[ind.pos[i], ]),
                                         Iteration = 0:M),
                              aes(Iteration, Prediction)) +
    geom_point(color = cl) + labs(x = "", y = "")
}

do.call("grid.arrange", c(pl.pred, ncol= 5))


ggplot(data.frame(Prediction = as.vector(t(expit(res.tree.mmr$train.lo[2, ]))),
                  MMR = rep(train$MMR[2], each = M + 1),
                  Iteration = rep(0:M, 2)),
       aes(Iteration, Prediction)) +
  geom_point(aes(color = as.factor(MMR))) +
  scale_color_discrete(name = "MMR", labels = Cs(No, Yes))


# med.gb <- cbind(risk.gb$FamID, apply(risk.gb[, -1], 1, median, na.rm = TRUE))
# med.gbc <- cbind(risk.gbc$FamID, apply(risk.gbc[, -1], 1, median, na.rm = TRUE))
# med.gbm <- cbind(risk.gbm$FamID, apply(risk.gbm[, -1], 1, median, na.rm = TRUE))


fn <- function(x){
  length(which(!is.na(x)))
}

summary(apply(risk.gb[, -1], 1, fn))
summary(apply(risk.gbm[, -1], 1, fn))

data.frame(EO.GB = sum(mean.gb[, 2]) / sum(sim.gb$MMR),
           AUC.GB = pROC::auc(sim.gb$MMR ~ mean.gb[, 2]),
           BS.GB = mean((sim.gb$MMR - mean.gb[, 2])^2),
           EO.GBM = sum(mean.gbm[, 2]) / sum(sim.gb$MMR),
           AUC.GBM = pROC::auc(sim.gb$MMR ~ mean.gbm[, 2]),
           BS.GBM =mean((sim.gb$MMR - mean.gbm[, 2])^2))

data.frame(EO.GB = sum(med.gb[, 2]) / sum(sim.gb$MMR),
           AUC.GB = pROC::auc(sim.gb$MMR ~ med.gb[, 2]),
           BS.GB = mean((sim.gb$MMR - med.gb[, 2])^2),
           EO.GBM = sum(med.gbm[, 2]) / sum(sim.gb$MMR),
           AUC.GBM = pROC::auc(sim.gb$MMR ~ med.gbm[, 2]),
           BS.GBM =mean((sim.gb$MMR - med.gbm[, 2])^2))

# res <- gb.mmr(sim.gb, c(can.names[-c(2:3)], "NumVert", "NumHorz"),
#               M, shrink, bag, seed = 1)



plot.gb <- function(sim.gb, mean.gb, method, init){
  method.type <- c("linreg", "tree")
  init.type <- c("mmr", "const")
  method.lab <- c("LR", "Tree")
  init.lab <- Cs(MMRPRO, Constant)
  if(method == "gbm"){
    yl <- "GBM (Tree), Constant"
  } else if(method == "xgb"){
    yl <- "XGB (Tree), Constant"
  } else{
    yl <- paste(method.lab[method.type == method], ", ",
                init.lab[init.type == init], sep = "")
  }
  ggplot(mutate(sim.gb, GB = mean.gb), aes(P.MMR, GB)) +
    geom_hex(data = filter(mutate(sim.gb, GB = mean.gb), MMR == 0),
             bins = 30) +
    geom_point(data = filter(mutate(sim.gb, GB = mean.gb), MMR == 1),
               aes(color = "red"), size = 0.5) +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "MMRPRO", y = "GB",
         title = yl) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_discrete(labels = "MMR carrier", name = "")
}

ggplot(mutate(sim.gb, Tree.MMR = mean.gb$tree.mmr,
              Tree.Const = mean.gb$tree.const), aes(Tree.Const, Tree.MMR)) +
  geom_hex(data = filter(mutate(sim.gb, Tree.MMR = mean.gb$tree.mmr,
                                Tree.Const = mean.gb$tree.const), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, Tree.MMR = mean.gb$tree.mmr,
                                  Tree.Const = mean.gb$tree.const), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Tree, Constant", y = "Tree, MMRPRO") +
  scale_color_discrete(labels = "MMR carrier", name = "")

grid.arrange(
  plot.gb(sim.gb, mean.gb$linreg.mmr, method = "linreg", init = "mmr"),
  plot.gb(sim.gb, mean.gb$linreg.const, method = "linreg", init = "const"),
  plot.gb(sim.gb, mean.gb$tree.mmr, method = "tree", init = "mmr"),
  plot.gb(sim.gb, mean.gb$tree.const, method = "tree", init = "const"),
  plot.gb(sim.gb, mean.gb$GBM, method = "gbm", init = "const"),
  plot.gb(sim.gb, mean.gb$XGB, method = "xgb", init = "const"),
  layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6))
)

ggplot(data.frame(pred = c(sim.gb$P.MMR, mean.gb$linreg.mmr, mean.gb$linreg.const,
                           mean.gb$tree.mmr, mean.gb$tree.const, mean.gb$GBM,
                           mean.gb$XGB),
                  type = rep(1:7, each = nrow(sim.gb)),
                  MMR = rep(sim.gb$MMR, 7)),
       aes(y = pred, x = as.factor(type), fill = as.factor(MMR))) +
  geom_boxplot() +
  scale_x_discrete(labels = Cs(MMRPRO, LR/MMR, LR/C, Tree/MMR, Tree/C, GBM, XGB)) +
  scale_fill_discrete(labels = Cs(Noncarrier, Carrier), name = "MMR") +
  labs(x = "", y = "Carrier Probability")

ggplot(filter(data.frame(pred = c(sim.gb$P.MMR, mean.gb$linreg.mmr, mean.gb$linreg.const,
                                  mean.gb$tree.mmr, mean.gb$tree.const, mean.gb$GBM),
                         type = rep(1:6, each = nrow(sim.gb)),
                         MMR = rep(sim.gb$MMR, 6)), MMR == 0),
       aes(log(pred))) +
  geom_density(aes(color = as.factor(type)), trim = TRUE) +
  scale_color_discrete(labels = Cs(MMRPRO, LR/MMR, LR/C, Tree/MMR, Tree/C, GBM),
                       name = "Type") +
  labs(x = "", y = "Carrier Probability")


geom_hex(data = filter(mutate(sim.gb, GB = mean.gb), MMR == 0),
         bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GB = mean.gb), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "Gradient Boosting Prediction",
       title = yl) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = "MMR carrier", name = "")




ggplot(mutate(sim.gb, GB = mean.gb$linreg.mmr), aes(P.MMR, GB)) +
  geom_hex(data = filter(mutate(sim.gb, GB = mean.gb$linreg.mmr), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GB = mean.gb$linreg.mmr), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "GB, LinReg, Constant") +
  scale_color_discrete(labels = "MMR carrier", name = "")

ggplot(mutate(sim.gb, GB = mean.gb$linreg.const), aes(P.MMR, GB)) +
  geom_hex(data = filter(mutate(sim.gb, GB = mean.gb$linreg.const), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GB = mean.gb$linreg.const), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "GB, LinReg, Constant") +
  scale_color_discrete(labels = "MMR carrier", name = "")





ggplot(mutate(res.gb.mmr$test, GBM = rp.gbm), aes(P.MMR, GBM)) +
  geom_hex(data = filter(mutate(res.gb.mmr$test, GBM = rp.gbm), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(res.gb.mmr$test, GBM = rp.gbm), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "Gradient Boosting R Package") +
  scale_color_discrete(labels = "MMR carrier", name = "")

summ <-
  rbind(
    summary(res.gb.mmr$p.test[test$MMR == 0]),
    summary(res.gb.mmr$p.test[test$MMR == 1]),
    summary(rp.gbm[test$MMR == 0]),
    summary(rp.gbm[test$MMR == 1])
  )
row.names(summ) <- c("GB, MMR = 0", "GB, MMR = 1",
                     "GBM, MMR = 0", "GBM, MMR = 1")

xtable(summ, digits = 4)





ggplot(mutate(sim.gb, GB = mean.gb[, 2]), aes(P.MMR, GB)) +
  geom_hex(data = filter(mutate(sim.gb, GB = mean.gb[, 2]), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GB = mean.gb[, 2]), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "Gradient Boosting") +
  scale_color_discrete(labels = "MMR carrier", name = "")

ggplot(mutate(sim.gb, GBM = mean.gbm[, 2]), aes(P.MMR, GBM)) +
  geom_hex(data = filter(mutate(sim.gb, GBM = mean.gbm[, 2]), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GBM = mean.gbm[, 2]), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "Gradient Boosting, gbm R package") +
  scale_color_discrete(labels = "MMR carrier", name = "")

ggplot(mutate(sim.gb, GB = mean.gb[, 2], GBM = mean.gbm[, 2]), aes(GB, GBM)) +
  geom_hex(data = filter(mutate(sim.gb, GB = mean.gb[, 2], GBM = mean.gbm[, 2]), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(sim.gb, GB = mean.gb[, 2], GBM = mean.gbm[, 2]), MMR == 1),
             aes(color = "red"), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Gradient Boosting", y = "Gradient Boosting, gbm R package") +
  scale_color_discrete(labels = "MMR carrier", name = "")


summ2 <-
  rbind(
    summary(mean.gb[sim.gb$MMR == 0, 2]),
    summary(mean.gb[sim.gb$MMR == 1, 2]),
    summary(mean.gbm[sim.gb$MMR == 0, 2]),
    summary(mean.gbm[sim.gb$MMR == 1, 2])
  )
row.names(summ2) <- c("GB, MMR = 0", "GB, MMR = 1",
                      "GBM, MMR = 0", "GBM, MMR = 1")
print(summ2, digits = 4)

summ.med <-
  rbind(
    summary(med.gb[sim.gb$MMR == 0, 2]),
    summary(med.gb[sim.gb$MMR == 1, 2]),
    summary(med.gbm[sim.gb$MMR == 0, 2]),
    summary(med.gbm[sim.gb$MMR == 1, 2])
  )
row.names(summ.med) <- c("GB, MMR = 0", "GB, MMR = 1",
                         "GBM, MMR = 0", "GBM, MMR = 1")
print(summ2, digits = 4)



ggplot(mutate(sim.gb, GB = mean.gb[, 2], GBM = mean.gbm[, 2]), aes(GB, GBM)) +
  geom_boxplot()

########



ggplot(melt(mutate(test, GB = p.test, GBM = rp.gbm),
            id.vars = "FamID", measure.vars = c("P.MMR", "GB", "GBM")),
       aes(variable, value)) +
  geom_boxplot()

ggplot(melt(filter(mutate(test, GB = p.test, GBM = rp.gbm), MMR == 0),
            id.vars = "FamID", measure.vars = c("P.MMR", "GB", "GBM")),
       aes(variable, value)) +
  geom_boxplot() +
  labs(title = "Non-carriers")

ggplot(melt(filter(mutate(test, GB = p.test, GBM = rp.gbm), MMR == 1),
            id.vars = "FamID", measure.vars = c("P.MMR", "GB", "GBM")),
       aes(variable, value)) +
  geom_boxplot() +
  labs(title = "Carriers")


ggplot(melt(mutate(test, GB = res.gb.mmr$p.test, GBM = rp.gbm),
            id.vars = "FamID", measure.vars = c("P.MMR", "GB", "GBM")),
       aes(variable, value)) +
  geom_boxplot()







# print(xtable(
#   data.frame(EO.GB = sum(p.test) / sum(test$MMR),
#              EO.MMR = sum(test$P.MMR) / sum(test$MMR),
#              AUC.GB = auc(test$MMR ~ p.test),
#              AUC.MMR = auc(test$MMR ~ test$P.MMR),
#              BS.GB = mean((p.test - test$MMR)^2),
#              BS.MMR = mean((test$P.MMR - test$MMR)^2)),
#   digits = 3), include.rownames = FALSE
# )
# 


