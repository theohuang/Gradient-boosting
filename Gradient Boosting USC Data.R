### Gradient boosting on USC data
### Last updated: September 24, 2018

library(BayesMendel)
library(dplyr)
library(data.table)
library(xgboost)
library(pROC)
library(BayesMendelKnowledgeBase)
library(BMmultigene)
library(doParallel)
registerDoParallel(cores = 4)

source("MMRpro.gast.R")
source("ImputeAge.gast.R")
source("CheckFamStructure.gast.R")
source("Simulation Functions.R")
source("Gradient Boosting Functions.R")
source("Heterozygous Penetrance.R")

## loading USC data
load("~/Dropbox (Partners HealthCare)/USC-Stanford_cohort/results_July12018.Rdata")

usc <- family
usc$FamID <- as.numeric(factor(usc$PedID))
N <- length(unique(usc$FamID)) ## 2000 families

## changing cancer diagnoses which are empty to NA
usc$Cancer.History.Table.Cancer.Diagnosis[which(usc$Cancer.History.Table.Cancer.Diagnosis == "")] <- NA
for(i in 1:9){
  usc[which(usc[, which(names(usc) == paste("Cancer.History.Table.Cancer.Diagnosis_", i, sep = ""))] == ""), 
      which(names(usc) == paste("Cancer.History.Table.Cancer.Diagnosis_", i, sep = ""))] <-
    NA
}

## adding gastric cancer ages
usc$AgeGastric <- NA
usc.ngc <- usc[usc$AffectedGastric == 0, ]
usc.ngc$AgeGastric <- ifelse(usc.ngc$AffectedBreast == 0, usc.ngc$AgeBreast,
                             ifelse(usc.ngc$AffectedOvary == 0, usc.ngc$AgeOvary,
                                    ifelse(usc.ngc$AffectedPancreas == 0, usc.ngc$AgePancreas,
                                           ifelse(usc.ngc$AffectedSkin == 0, usc.ngc$AgeSkin,
                                                  ifelse(usc.ngc$AffectedColon == 0, usc.ngc$AgeColon,
                                                         ifelse(usc.ngc$AffectedEndometrium == 0, usc.ngc$AgeEndometrium,
                                                                99999))))))
usc$AgeGastric[usc$AffectedGastric == 0] <- usc.ngc$AgeGastric
for(i in 1:length(which(usc$AffectedGastric == 1))){
  ind <- which(usc$AffectedGastric == 1)[i]
  if(!is.na(usc$Cancer.History.Table.Cancer.Diagnosis[ind]) &
     usc$Cancer.History.Table.Cancer.Diagnosis[ind] == "Stomach"){
    usc$AgeGastric[ind] <- usc$Cancer.History.Table.Diagnosis.Age[ind]
  } else{
    for(j in 1:9){
      if(!is.na(usc[ind, which(names(usc) == paste("Cancer.History.Table.Cancer.Diagnosis_", j, sep = ""))]) &
         usc[ind, which(names(usc) == paste("Cancer.History.Table.Cancer.Diagnosis_", j, sep = ""))] == "Stomach"){
        usc$AgeGastric[ind] <- as.numeric(usc[ind, which(names(usc) == paste("Cancer.History.Table.Diagnosis.Age_", j, sep = ""))])
      }
    }
  }
}
usc$AgeGastric <- as.numeric(usc$AgeGastric)
usc$AgeGastric[usc$AgeGastric == 0] <- NA


## for family 26, getting rid of some unrelated family members to fix peeling issue
usc <- usc[-which(usc$FamID == 26 & !(usc$ID %in% c(1, 5, 6, 37:41))), ]


###### Penetrances

## MMRpro penetrance without gastric cancer
pen.ngc <- MMRparams()$penetrance.net

## MMRpro penetrance with gastric cancer
cancers <- c("ColorC", "EndomC", "GastC")
genes <- c("MLH1", "MSH2", "MSH6")
penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")
# using the MMRPRO penetrances for colon and endometrial cancer
penCancersF$ColorC <- penet.mmr.net$fFX[, c("M000", "M100", "M010", "M001")]
penCancersM$ColorC <- penet.mmr.net$fMX[, c("M000", "M100", "M010", "M001")]
penCancersF$EndomC <- penet.mmr.net$fFY[, c("M000", "M100", "M010", "M001")]
penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")
CP <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)

penCancersF.crude <- penCancersF
penCancersF.crude$ColorC <- penet.mmr.crude$fFX[, c("M000", "M100", "M010", "M001")]
penCancersF.crude$EndomC <- penet.mmr.crude$fFY[, c("M000", "M100", "M010", "M001")]
penCancersM.crude <- penCancersM
penCancersM.crude$ColorC <- penet.mmr.crude$fMX[, c("M000", "M100", "M010", "M001")]
CP.crude <- genCancerPen(genes, cancers, penCancersF.crude, penCancersM.crude, maxK = length(genes), age.last = 95)

## function to run MMRPRO, including options to include gastric cancer
## and to raise the gastric cancer survival function to a power
mmr.gb <- function(fam, CP, CP.crude, gastric = TRUE, pwr = 1){
  if(gastric == TRUE){
    # fam$AgeGastric <- ifelse(fam$AffectedPancreas != 1, fam$AgePancreas,
    #                          ifelse(fam$AffectedEndometrium != 1, fam$AgeEndometrium,
    #                                 ifelse(fam$AffectedOvary != 1, fam$AgeOvary)))
    return(MMRpro.gast(fam, counselee.id = fam$ID[fam$Proband == 1],
                       params = MMRparams(penetrance.net = hetpen(CP, gastric = TRUE, pwr = pwr),
                                          penetrance.crude = hetpen(CP.crude, gastric = TRUE, pwr = pwr)),
                       print = FALSE)@probs)
  } else{
    return(MMRpro(fam, counselee.id = fam$ID[fam$Proband == 1],
                  params = MMRparams(penetrance.net = hetpen(CP, gastric = FALSE, pwr = pwr),
                                     penetrance.crude = hetpen(CP.crude, gastric = FALSE, pwr = pwr)),
                  print = FALSE)@probs)
  }
}


## Running MMRPRO with and without gastric cancer, with different scales
usc.mmrpro <- foreach(i = 1:N, .combine = rbind) %dopar% {
  unlist(c(tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude, gastric = FALSE),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude, pwr = 0.25),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude, pwr = 0.5),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude, pwr = 2),
             error = function(e) rep(NA, 4)),
    tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, CP.crude, pwr = 4),
             error = function(e) rep(NA, 4))))
}

# for(i in 1:N){
#   print(i)
#   c(tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, gastric = FALSE),
#              error = function(e) rep(NA, 4)),
#     tryCatch(mmr.gb(usc[usc$FamID == i, ], CP),
#              error = function(e) rep(NA, 4)),
#     tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, pwr = 0.25),
#              error = function(e) rep(NA, 4)),
#     tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, pwr = 0.5),
#              error = function(e) rep(NA, 4)),
#     tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, pwr = 2),
#              error = function(e) rep(NA, 4)),
#     tryCatch(mmr.gb(usc[usc$FamID == i, ], CP, pwr = 4),
#              error = function(e) rep(NA, 4)))
# }



usc.mmrpro <- data.frame(usc.mmrpro)
names(usc.mmrpro) <- c("P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc",
                       "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6",
                       "P.MMR.025", "P.MLH1.025", "P.MSH2.025", "P.MSH6.025",
                       "P.MMR.05", "P.MLH1.05", "P.MSH2.05", "P.MSH6.05",
                       "P.MMR.2", "P.MLH1.2", "P.MSH2.2", "P.MSH6.2",
                       "P.MMR.4", "P.MLH1.4", "P.MSH2.4", "P.MSH6.4")


## Seeing if the probands are carriers
usc.gb <- usc.mmrpro
usc.gb$MLH1 <- as.numeric(filter(usc, Proband == 1)$GenTestTable.Gene.1 == "MLH1")
usc.gb$MSH2 <- as.numeric(filter(usc, Proband == 1)$GenTestTable.Gene.1 == "MSH2")
usc.gb$MSH6 <- as.numeric(filter(usc, Proband == 1)$GenTestTable.Gene.1 == "MSH6")
usc.gb$MMR <- as.numeric(rowSums(cbind(usc.gb$MLH1, usc.gb$MSH2, usc.gb$MSH6)) > 0)

## Getting family information for the three cancers
usc.gb$FamID <- 1:N
can.names <- paste("Prop", cancers, sep = "")
usc.gb[can.names] <- NA
for(i in 1:nrow(usc.gb)){
  usc.gb[can.names][usc.gb$FamID == i, ] <- colMeans(usc[usc$FamID == i, ][c("AffectedColon", "AffectedEndometrium", "AffectedGastric")])
}



## Gradient boosting
nm.list <- c("MMR.ngc", "MMR", "MMR.025", "MMR.05", "MMR.2", "MMR.4", "XGB.mmr", "XGB.const")
n.boot <- 100
M <- 100
covs <- c(can.names)
shrink <- 0.1
bag <- 0.5

res.gb <- setNames(vector("list", length(nm.list)), nm.list)
res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                          c("EO", "AUC", "BS")), 
                                          risk = setNames(data.frame(cbind(usc.gb$FamID, matrix(NA, nrow(usc.gb), n.boot))),
                                                          c("FamID", paste("risk", 1:n.boot, sep = "")))))


start <- Sys.time()
for(i in 1:n.boot){
  print(i)
  smp.train <- sample(1:nrow(usc.gb), floor(nrow(usc.gb) / 2))
  train <- usc.gb[smp.train, ]
  test <- usc.gb[-smp.train, ]
  
  ## MMRpro without GC
  res.gb$MMR.ngc$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.ngc)
  res.gb$MMR.ngc$risk[res.gb$MMR.ngc$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.ngc
  
  ## MMRpro with GC
  res.gb$MMR$perf[i, ] <- perf.meas(test$MMR, test$P.MMR)
  res.gb$MMR$risk[res.gb$MMR$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR
  
  res.gb$MMR.025$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.025)
  res.gb$MMR.025$risk[res.gb$MMR.025$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.025
  
  res.gb$MMR.05$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.05)
  res.gb$MMR.05$risk[res.gb$MMR.05$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.05
  
  res.gb$MMR.2$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.2)
  res.gb$MMR.2$risk[res.gb$MMR.2$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.2
  
  res.gb$MMR.4$perf[i, ] <- perf.meas(test$MMR, test$P.MMR.4)
  res.gb$MMR.4$risk[res.gb$MMR.4$risk$FamID %in% test$FamID, i + 1] <- test$P.MMR.4
  
  ## XGBoost with MMRPRO
  param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                objective = "binary:logistic")
  dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                            base_margin = logit(train$P.MMR.ngc))
  dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                           base_margin = logit(test$P.MMR.ngc))
  # watchlist <- list(train = dtrain.mmr, test = dtest.mmr)
  # res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M, watchlist,
  #                          eval_metric = "auc", eval_metric = "rmse")
  res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M)
  # cv <- xgb.cv(params = param, data = dtrain.mmr, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
  pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
  res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
  res.gb$XGB.mmr$risk[usc.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
  
  ## XGBoost with constant
  dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
  dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
  # watchlist <- list(train = dtrain.const, test = dtest.const)
  # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
  #                            eval_metric = "auc", eval_metric = "rmse")
  res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M)
  # cv.const <- xgb.cv(params = param, data = dtrain.const, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
  pred.xgb.const <- predict(res.xgb.const, dtest.const)
  res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
  res.gb$XGB.const$risk[usc.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
}
difftime(Sys.time(), start, units = "secs")

save(res.gb, file = paste(getwd(), "/Gradient Boosting/Diff10/USC/uscgbpen_", a1, ".RData", sep = ""))


