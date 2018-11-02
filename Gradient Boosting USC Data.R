### Gradient boosting on USC data
### Last updated: October 22, 2018

library(BayesMendel)
library(dplyr)
library(data.table)
library(xgboost)
library(pROC)
library(BayesMendelKnowledgeBase)
library(BMmultigene)
library(ggplot2)
library(doParallel)
registerDoParallel(cores = 4)

source("MMRpro.gast.R")
source("ImputeAge.gast.R")
source("CheckFamStructure.gast.R")
source("Simulation Functions.R")
source("Gradient Boosting Functions.R")
source("Heterozygous Penetrance.R")

## loading USC data
load("usc_dat.Rdata")
famid.list <- unique(usc$FamID)
N <- length(famid.list)

###### Penetrances

## MMRpro penetrance without gastric cancer
# net
cancers <- c("ColorC", "EndomC", "GastC")
genes <- c("MLH1", "MSH2", "MSH6")
penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")
# using the MMRPRO penetrances for colon and endometrial cancer
penCancersF$ColorC <- penet.mmr.net$fFX[, c("M000", "M100", "M010", "M001")]
penCancersM$ColorC <- penet.mmr.net$fMX[, c("M000", "M100", "M010", "M001")]
penCancersF$EndomC <- penet.mmr.net$fFY[, c("M000", "M100", "M010", "M001")]
CP <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)
# crude
penCancersF.crude <- penCancersF
penCancersF.crude$ColorC <- penet.mmr.crude$fFX[, c("M000", "M100", "M010", "M001")]
penCancersF.crude$EndomC <- penet.mmr.crude$fFY[, c("M000", "M100", "M010", "M001")]
penCancersM.crude <- penCancersM
penCancersM.crude$ColorC <- penet.mmr.crude$fMX[, c("M000", "M100", "M010", "M001")]
CP.crude <- genCancerPen(genes, cancers, penCancersF.crude, penCancersM.crude, maxK = length(genes), age.last = 95)


## MMRpro penetrance with gastric cancer
load("GCPenetrance.RData")


# save(usc, CP, CP.crude, file = "usc_gastpen.RData")



## function to run MMRPRO, including options to include gastric cancer
## and to raise the gastric cancer survival function to a power
mmr.gb <- function(fam, CP, CP.crude, gt = TRUE, mt = TRUE, gastric = TRUE, pwr = 1){
  
  ## genetic testing (only for family members, not proband)
  # if we turn on genetic testing (gt = TRUE) and at least one family member
  # (not the proband) has genetic testing results, we will use the genetic testing results
  if(gt == TRUE & !all(apply(fam[fam$Proband != 1, c("MLH1", "MSH2", "MSH6")], 1,  function(x) all(x == 0)))){
    germline.testing <- data.frame(MLH1 = fam$MLH1, MSH2 = fam$MSH2, MSH6 = fam$MSH6,
                                   TestOrder = rep(0, nrow(fam)))
    # changing the proband's testing results to unknown
    germline.testing[fam$Proband == 1, ] <- rep(0, 4)
  } else{
    germline.testing <- NULL
  }
  
  ## marker testing
  # if we turn on marker testing (mt = TRUE) and at least one family member
  # (including the proband) has marker testing results, we will use the results
  if(mt == TRUE & !all(fam$MSI.mmr == 0)){
    marker.testing <- data.frame(MSI = fam$MSI.mmr, location = rep(0, nrow(fam)))
  } else{
    marker.testing <- NULL
  }
  
  ## running MMRpro
  if(gastric == TRUE){
    return(MMRpro.gast(fam, counselee.id = fam$ID[fam$Proband == 1], race = fam$race.mmr[1],
                       germline.testing = germline.testing, marker.testing = marker.testing,
                       params = MMRparams(penetrance.net = pen.gc.pwr(penet.mmr.net.gc, pwr),
                                          penetrance.crude = pen.gc.pwr(penet.mmr.crude.gc, pwr)),
                       print = FALSE)@probs)
  } else{
    return(MMRpro(fam, counselee.id = fam$ID[fam$Proband == 1], race = fam$race.mmr[1],
                  germline.testing = germline.testing, marker.testing = marker.testing,
                  params = MMRparams(penetrance.net = hetpen(CP, gastric = FALSE, pwr = pwr),
                                     penetrance.crude = hetpen(CP.crude, gastric = FALSE, pwr = pwr)),
                  print = FALSE)@probs)
  }
}



## Running MMRPRO with and without gastric cancer, with different scales
start <- Sys.time()
usc.mmrpro <- foreach(i = 1:N, .combine = rbind) %dopar% {
  unlist(c(famid.list[i],
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude, gastric = FALSE),
                    error = function(e) rep(NA, 4)),
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude),
                    error = function(e) rep(NA, 4)),
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude, pwr = 0.25),
                    error = function(e) rep(NA, 4)),
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude, pwr = 0.5),
                    error = function(e) rep(NA, 4)),
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude, pwr = 2),
                    error = function(e) rep(NA, 4)),
           tryCatch(mmr.gb(usc[usc$FamID == famid.list[i], ], CP, CP.crude, pwr = 4),
                    error = function(e) rep(NA, 4))))
}
difftime(Sys.time(), start, units = "secs")


usc.mmrpro <- data.frame(usc.mmrpro)
names(usc.mmrpro) <- c("FamID", "P.MMR.ngc", "P.MLH1.ngc", "P.MSH2.ngc", "P.MSH6.ngc",
                       "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6",
                       "P.MMR.025", "P.MLH1.025", "P.MSH2.025", "P.MSH6.025",
                       "P.MMR.05", "P.MLH1.05", "P.MSH2.05", "P.MSH6.05",
                       "P.MMR.2", "P.MLH1.2", "P.MSH2.2", "P.MSH6.2",
                       "P.MMR.4", "P.MLH1.4", "P.MSH2.4", "P.MSH6.4")
# usc.mmrpro$FamID <- famid.list


## Getting family information for the three cancers
usc.gb <- usc.mmrpro
# usc.gb$FamID <- famid.list
usc.gb <- filter(usc.gb, !is.na(P.MMR))
famid.gb <- unique(usc.gb$FamID)
can.names <- paste("Prop", cancers, sep = "")
usc.gb[can.names] <- NA
for(i in 1:nrow(usc.gb)){
  usc.gb[c("PropColorC", "PropGastC")][usc.gb$FamID == famid.gb[i], ] <- colMeans(usc[usc$FamID == famid.gb[i], ][c("AffectedColon", "AffectedGastric")])
  # using proportion of female family members for endometrial cancer
  usc.gb$PropEndomC[usc.gb$FamID == famid.gb[i]] <- mean(usc$AffectedEndometrium[usc$FamID == famid.gb[i] & usc$Gender == 0])
}

## adding the MMR carrier statuses
usc.gb[c("MLH1", "MSH2", "MSH6", "MMR")] <- select(filter(usc, Proband == 1, FamID %in% famid.gb), MLH1, MSH2, MSH6, MMR)

# usc.gb <- filter(usc.gb, !is.na(P.MMR))

## identifying the families who had NA MMRpro results
id.na <- famid.list[which(is.na(usc.mmrpro$P.MMR))]


## removing families with probands who have VUS results for the MMR genes
usc.vus <- usc.gb
id.vus <- c(filter(usc, Proband == 1, GenTestTable.Gene.1 %in% c("MLH1", "MSH2", "MSH6"),
                   GenTestTable.Overall.Result == "VUS")$FamID,
            filter(usc, Proband == 1, GenTestTable.Gene.1 %in% c("MLH1", "MSH2", "MSH6"),
                   GenTestTable.Overall.Result == "VUS - Suspect Del/ Likely Pathogenic")$FamID)
usc.gb <- usc.gb[-which(usc.gb$FamID %in% id.vus), ]

## adding EPCAM and PMS2
usc.ep <- usc.gb
usc.ep$MMR[usc.ep$FamID %in% filter(usc, FamID %in% usc.ep$FamID, Proband == 1,
                                    GenTestTable.Gene.1 %in% c("EPCAM", "PMS2"),
                                    GenTestTable.Overall.Result == "Positive")$FamID] <- 1

## removing families with probands who tested positive for some other gene
id.othergene <- filter(usc, Proband == 1, GenTestTable.Overall.Result == "Positive",
                       FamID %in% usc.gb$FamID, !(GenTestTable.Gene.1 %in% c("MLH1", "MSH2", "MSH6")))$FamID
id.lynchgene <- filter(usc, Proband == 1, GenTestTable.Overall.Result == "Positive",
                       FamID %in% usc.gb$FamID, !(GenTestTable.Gene.1 %in% c("MLH1", "MSH2", "MSH6")),
                       GenTestTable.Gene.1 %in% c("APC", "CDH1", "CHEK2", "EPCAM",
                                                  "MUTYH", "PMS2", "TP53"))$FamID
usc.og <- usc.gb[-which(usc.gb$FamID %in% id.othergene), ]
usc.lg <- usc.gb[-which(usc.gb$FamID %in% id.lynchgene), ]


## Gradient boosting
# nm.list <- c("MMR.ngc", "MMR", "MMR.025", "MMR.05", "MMR.2", "MMR.4", "XGB.mmr", "XGB.const")
# n.boot <- 100
# M. <- 100
# covs <- c(can.names)
# shrink <- 0.1
# bag <- 0.5


types <- c(".ngc", "", ".025", ".05", ".2", ".4")
n.boot <- 1000
M.mmr <- 35
M.const <- 60
covs <- can.names
shrink <- 0.1
bag <- 0.5

start <- Sys.time()
res.gb <- gb.mmr(usc.gb, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.ep <- gb.mmr(usc.ep, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.vus <- gb.mmr(usc.vus, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.og <- gb.mmr(usc.og, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")

start <- Sys.time()
res.lg <- gb.mmr(usc.lg, shrink, bag, M.mmr, M.const, covs, n.boot, types)
difftime(Sys.time(), start, units = "secs")

# start <- Sys.time()
# res.gb <- gb.mmr(usc.vus, shrink, bag, M.mmr, M.const, covs, n.boot, types)
# difftime(Sys.time(), start, units = "secs")



# res.gb <- setNames(vector("list", length(nm.list)), nm.list)
# res.gb <- lapply(res.gb, function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
#                                                           c("EO", "AUC", "BS")), 
#                                           risk = setNames(data.frame(cbind(usc.gb$FamID, matrix(NA, nrow(usc.gb), n.boot))),
#                                                           c("FamID", paste("risk", 1:n.boot, sep = "")))))
# 
# 
# start <- Sys.time()
# for(i in 1:n.boot){
#   print(i)
#   smp.train <- sample(1:nrow(usc.gb), floor(nrow(usc.gb) / 2))
#   train <- usc.gb[smp.train, ]
#   test <- usc.gb[-smp.train, ]
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
#   res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M)
#   # cv <- xgb.cv(params = param, data = dtrain.mmr, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
#   pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
#   res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
#   res.gb$XGB.mmr$risk[usc.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
#   
#   ## XGBoost with constant
#   dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
#   dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
#   # watchlist <- list(train = dtrain.const, test = dtest.const)
#   # res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M, watchlist,
#   #                            eval_metric = "auc", eval_metric = "rmse")
#   res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M)
#   # cv.const <- xgb.cv(params = param, data = dtrain.const, nrounds = M, nfold = 5, metrics = list("auc", "rmse"))
#   pred.xgb.const <- predict(res.xgb.const, dtest.const)
#   res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
#   res.gb$XGB.const$risk[usc.gb$FamID %in% test$FamID, i + 1] <- pred.xgb.const
# }
# difftime(Sys.time(), start, units = "secs")


lapply(res.gb, function(x) colMeans(x$perf))

res.summ <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ[i, c(1, 4, 7)] <- colMeans(res.gb[[i]]$perf)
  res.summ[i, c(2, 5, 8)] <- apply(res.gb[[i]]$perf, 2, function(x) sort(x)[25])
  res.summ[i, c(3, 6, 9)] <- apply(res.gb[[i]]$perf, 2, function(x) sort(x)[975])
}

row.names(res.summ) <- names(res.gb)
colnames(res.summ) <- c("OE mean", "OE 2.5", "OE 97.5",
                        "AUC mean", "AUC 2.5", "AUC 97.5",
                        "BS mean", "BS 2.5", "BS 97.5")
res.summ

## with EPCAM and PMS2
res.summ.ep <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.ep[i, c(1, 4, 7)] <- colMeans(res.ep[[i]]$perf)
  res.summ.ep[i, c(2, 5, 8)] <- apply(res.ep[[i]]$perf, 2, function(x) sort(x)[25])
  res.summ.ep[i, c(3, 6, 9)] <- apply(res.ep[[i]]$perf, 2, function(x) sort(x)[975])
}

row.names(res.summ.ep) <- names(res.ep)
colnames(res.summ.ep) <- c("OE mean", "OE 2.5", "OE 97.5",
                           "AUC mean", "AUC 2.5", "AUC 97.5",
                           "BS mean", "BS 2.5", "BS 97.5")
res.summ.ep


## including VUS and assuming they are negative
res.summ.vus <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.vus[i, c(1, 4, 7)] <- colMeans(res.vus[[i]]$perf)
  res.summ.vus[i, c(2, 5, 8)] <- apply(res.vus[[i]]$perf, 2, function(x) sort(x)[25])
  res.summ.vus[i, c(3, 6, 9)] <- apply(res.vus[[i]]$perf, 2, function(x) sort(x)[975])
}

row.names(res.summ.vus) <- names(res.vus)
colnames(res.summ.vus) <- c("OE mean", "OE 2.5", "OE 97.5",
                            "AUC mean", "AUC 2.5", "AUC 97.5",
                            "BS mean", "BS 2.5", "BS 97.5")
res.summ.vus

## excluding probands who tested positive for other genes
res.summ.og <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.og[i, c(1, 4, 7)] <- colMeans(res.og[[i]]$perf)
  res.summ.og[i, c(2, 5, 8)] <- apply(res.og[[i]]$perf, 2, function(x) sort(x)[25])
  res.summ.og[i, c(3, 6, 9)] <- apply(res.og[[i]]$perf, 2, function(x) sort(x)[975])
}

row.names(res.summ.og) <- names(res.og)
colnames(res.summ.og) <- c("OE mean", "OE 2.5", "OE 97.5",
                           "AUC mean", "AUC 2.5", "AUC 97.5",
                           "BS mean", "BS 2.5", "BS 97.5")
res.summ.og

## excluding probands who tested positive for genes associated with Lynch syndrome
res.summ.lg <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.lg[i, c(1, 4, 7)] <- colMeans(res.lg[[i]]$perf)
  res.summ.lg[i, c(2, 5, 8)] <- apply(res.lg[[i]]$perf, 2, function(x) sort(x)[25])
  res.summ.lg[i, c(3, 6, 9)] <- apply(res.lg[[i]]$perf, 2, function(x) sort(x)[975])
}

row.names(res.summ.lg) <- names(res.lg)
colnames(res.summ.lg) <- c("OE mean", "OE 2.5", "OE 97.5",
                           "AUC mean", "AUC 2.5", "AUC 97.5",
                           "BS mean", "BS 2.5", "BS 97.5")
res.summ.lg




id.ep <- filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                GenTestTable.Gene.1 %in% c("PMS2", "EPCAM", "MUTYH", "STK11", "TP53"))$FamID
id.mmr <- filter(usc, FamID %in% usc.gb$FamID, Proband == 1, MMR == 1)$FamID

plot(usc.gb$P.MMR.ngc, ylab = "MMRpro Prediction without GC")
points(which(usc.gb$FamID %in% id.ep),
       filter(usc.gb, FamID %in% id.ep)$P.MMR.ngc,
       col = "red")
points(which(usc.gb$FamID %in% id.mmr),
       filter(usc.gb, FamID %in% id.mmr)$P.MMR.ngc,
       col = "blue")


risk.mmr <- filter(usc.gb, MMR == 1)$P.MMR.ngc
risk.apc <- filter(usc.gb, FamID %in% filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                                             GenTestTable.Gene.1 == "APC")$FamID)$P.MMR.ngc
risk.epcam <- filter(usc.gb, FamID %in% filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                                               GenTestTable.Gene.1 == "EPCAM")$FamID)$P.MMR.ngc
risk.pms2 <- filter(usc.gb, FamID %in% filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                                              GenTestTable.Gene.1 == "PMS2")$FamID)$P.MMR.ngc
risk.mutyh <- filter(usc.gb, FamID %in% filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                                               GenTestTable.Gene.1 == "MUTYH")$FamID)$P.MMR.ngc
risk.tp53 <- filter(usc.gb, FamID %in% filter(usc, FamID %in% usc.gb$FamID, GenTestTable.Overall.Result == "Positive", Proband == 1,
                                              GenTestTable.Gene.1 == "TP53")$FamID)$P.MMR.ngc

boxplot(risk.mmr, risk.apc, risk.epcam, risk.pms2, risk.mutyh, risk.tp53,
        names = c("MMR", "APC", "EPCAM", "PMS2", "MUTYH", "TP53"),
        main = "MMRpro Predictions (no GC) for Carriers")



risk.xgb.mmr <- rowMeans(res.og$XGB.mmr$risk[, -1], na.rm = TRUE)

plot(usc.gb$P.MMR.ngc, risk.xgb.mmr)
abline(0,1)

df <- dplyr::arrange(data.frame(P.MMR = usc.og$P.MMR.ngc,
                                GB = risk.xgb.mmr,
                                MMR = usc.og$MMR), MMR)
ggplot(df, aes(P.MMR, GB)) +
  geom_point(aes(color = factor(MMR))) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRpro Prediction without Gastric Cancer",
       y = "Mean Gradient Boosting Prediction") +
  scale_color_discrete(name = "MMR Carriers",
                       labels = c("No", "Yes"))


library(pROC)
sum(df$MMR) / sum(df$GB)
pROC::auc(MMR ~ GB, data = df)
mean((df$MMR - df$GB)^2)
sum(df$MMR) / sum(df$P.MMR)
pROC::auc(MMR ~ P.MMR, data = df)
mean((df$MMR - df$P.MMR)^2)


table(filter(usc, Proband == 1, FamID %in% filter(usc.gb, P.MMR.ngc > 0.05)$FamID,
             GenTestTable.Overall.Result == "Positive")$GenTestTable.Gene.1)

table(filter(usc, Proband == 1, FamID %in% usc.gb$FamID[risk.xgb.mmr > 0.05],
             GenTestTable.Overall.Result == "Positive")$GenTestTable.Gene.1)




save(usc, usc.mmrpro, usc.gb, usc.vus, usc.ep, usc.og, usc.lg,
     res.gb, res.vus, res.ep, res.og, res.lg,
     res.summ, res.summ.vus, res.summ.ep, res.summ.og, res.summ.lg,
     file = "GB_USC_Results_101118.RData")


library(xtable)
xtab <- data.frame(round(res.summ.og, 3))
xtab$OE <- xtab$AUC <- xtab$BS <- 0
for(i in 1:8){
  xtab$OE[i] <- paste(xtab$OE.mean[i], " (", xtab[i, 2], ", ", xtab[i, 3], ")", sep = "")
  xtab$AUC[i] <- paste(xtab$AUC.mean[i], " (", xtab[i, 5], ", ", xtab[i, 6], ")", sep = "")
  xtab$BS[i] <- paste(xtab$BS.mean[i], " (", xtab[i, 8], ", ", xtab[i, 9], ")", sep = "")
}
xtab <- xtab[, 12:10]

rownames(xtab) <- c("MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                    "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab)


## plotting the scaled penetrances
ggplot(data.frame(Penetrance = c(penCancersM$GastC[, 1],
                                 homoz.genes(penCancersM$GastC[, 1], 0.25),
                                 homoz.genes(penCancersM$GastC[, 1], 0.5),
                                 homoz.genes(penCancersM$GastC[, 1], 2),
                                 homoz.genes(penCancersM$GastC[, 1], 4)),
                  Type = rep(1:5, each = 94),
                  Age = rep(1:94, 5)),
       aes(Age, Penetrance)) +
  geom_point(aes(color = factor(Type))) +
  scale_color_discrete(name = "Power", labels = c("1", "0.25", "0.5", "2", "4"))

ggplot(data.frame(Penetrance = c(penCancersM$GastC[, 2],
                                 homoz.genes(penCancersM$GastC[, 2], 0.25),
                                 homoz.genes(penCancersM$GastC[, 2], 0.5),
                                 homoz.genes(penCancersM$GastC[, 2], 2),
                                 homoz.genes(penCancersM$GastC[, 2], 4)),
                  Type = rep(1:5, each = 94),
                  Age = rep(1:94, 5)),
       aes(Age, Penetrance)) +
  geom_point(aes(color = factor(Type))) +
  scale_color_discrete(name = "Power", labels = c("1", "0.25", "0.5", "2", "4"))




# save(res.gb, usc.gb, usc.mmrpro, file = "GB_USC_Results.RData")

# save(res.gb, file = paste(getwd(), "/Gradient Boosting/Diff10/USC/uscgbpen_", a1, ".RData", sep = ""))


