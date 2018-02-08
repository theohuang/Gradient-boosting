### Gradient Boosting for Old CGN Data
### Modelling carrier probabilities for BRCA12/2
### Last updated: February 7, 2018

setwd("~/Dropbox (Partners HealthCare)/Gradient Boosting")

## Loading functions
source(paste(getwd(), "/Gradient Boosting Basis Functions.R", sep = ""))

## Loading data
# load(paste(getwd(), "/CGNValidationData/data/PedLevelVars.RData", sep = ""))
# load(paste(getwd(), "/CGNValidationData/data/FamLevelVars.RData", sep = ""))
# load("/Users/Theo/Documents/Harvard/Research with Giovanni Parmigiani/Gradient Boosting/CGNResults012417.RData")
library(CGNValidationData)

## Loading BayesMendel
library(BayesMendel)
library(dplyr)
set.seed(1)

## This runs BRCAPRO, incorporating germline testing
brcapro.fam <- function(fam, seed = 1){
  set.seed(seed)
  ## BRCA12/2 information is coded as a factor in the data, with levels 0 and 1, so
  ## I'm converting it to numeric 0 and 1
  germline.testing <- data.frame(BRCA1 = as.numeric(fam$Germline.BRCA1Result) - 1,
                                 BRCA2 = as.numeric(fam$Germline.BRCA2Result) - 1,
                                 TestOrder = rep(0, nrow(fam)))
  ind.pro <- which(fam$Relation == 1)
  ## not using germline testing results for proband
  germline.testing[ind.pro, ] <- rep(0, 3)
  ## not using germline testing results for the proband's twins
  if(fam$Twins[ind.pro] > 0){
    ind.twin <- which(fam$Twins == fam$Twins[ind.pro] & fam$Relation != 1)
    germline.testing[ind.twin, ] <- rep(0, 3)
  }
  result <- brcapro(family = fam, counselee.id = fam$ID[ind.pro],
                    germline.testing = germline.testing, print = FALSE)
  return(result)
}


## Running BRCAPRO with and without ovarian cancer information
famid <- cbind(unique(PedLevelVars[, c("Center.ID", "File.ID")]), 1:nrow(FamLevelVars))
names(famid)[3] <- "FamID"
cgn <- merge(PedLevelVars, famid, by = c("Center.ID", "File.ID"))
cgn$Death <- NA; cgn$AgeDeath <- NA

FamLevelVars$FamID <- NA
for(i in 1:nrow(FamLevelVars)){
  FamLevelVars$FamID[i] <- cgn$FamID[cgn$Center.ID == FamLevelVars$Center.ID[i] &
                                       cgn$File.ID == FamLevelVars$File.ID[i]][1]
}

# cgn.boost <- data.frame(matrix(0, length(unique(cgn$FamID)), 5))
# names(cgn.boost) <- c("FamID", "P.BRCA", "P.BRCA.NoOC", "P.BRCA.NoAJ", "P.BRCA.NoOCAJ")
# cgn.boost$FamID <- 1:length(unique(cgn$FamID))
# start <- Sys.time()
# for(i in 1:length(unique(cgn$FamID))){
#   cgn.boost$P.BRCA[i] <- brcapro.fam(cgn[cgn$FamID == i, ])@probs[1]
#   dat.nooc <- cgn[cgn$FamID == i, ]
#   dat.nooc$AffectedOvary <- 0
#   dat.nooc$AgeOvary <- dat.nooc$Current.Age
#   cgn.boost$P.BRCA.NoOC[i] <- brcapro.fam(dat.nooc)@probs[1]
#   dat.noaj <- cgn[cgn$FamID == i, ]
#   dat.noaj$ethnic <- NA
#   cgn.boost$P.BRCA.NoAJ[i] <- brcapro.fam(dat.noaj)@probs[1]
#   dat.noocaj <- dat.nooc; dat.noocaj$ethnic <- NA
#   cgn.boost$P.BRCA.NoOCAJ[i] <- brcapro.fam(dat.noocaj)@probs[1]
# }
# print(difftime(Sys.time(), start, units = "secs"))

cgn.boost <- data.frame(FamID = 1:length(unique(cgn$FamID)))
library(data.table)
run.brca <- function(fam){
  dat.nooc <- fam; dat.nooc$AffectedOvary <- 0; dat.nooc$AgeOvary <- dat.nooc$Current.Age
  dat.noaj <- fam; dat.nooc$ethnic <- NA
  dat.noocaj <- dat.nooc; dat.nooc$ethnic <- NA
  return(c(brcapro.fam(fam)@probs[1], brcapro.fam(dat.nooc)@probs[1],
           brcapro.fam(dat.noaj)@probs[1], brcapro.fam(dat.noocaj)@probs[1]))
}
start <- Sys.time()
dt <- data.table(cgn)
dt <- dt[, run.brca(.SD), by = FamID]
print(difftime(Sys.time(), start, units = "secs"))
names(dt)[-1] <- c("P.BRCA", "P.BRCA.NoOC", "P.BRCA.NoAJ", "P.BRCA.NoOCAJ")
cgn.boost <- merge(cgn.boost, dt, by = "FamID")

# famid.twin <- filter(cgn, Relation == 1, Twins > 0)$FamID
# dt.twin <- data.table(filter(cgn, FamID %in% famid.twin))
# dt.twin <- dt.twin[, run.brca(.SD), by = FamID]
# cgn.boost[cgn.boost$FamID %in% famid.twin, ] <- data.frame(dt.twin)


germ <- c(cgn$FamID[cgn$Germline.BRCA1Result == 1],
          cgn$FamID[cgn$Germline.BRCA2Result == 1])

for(i in germ){
  fam <- cgn[cgn$FamID == i, ]
  cid <- fam$ID[fam$Relation == 1]
  cgn.boost$P.BRCA[i] <- brcapro(fam, cid,
                                 print = FALSE)@probs[1]
  dat.nooc <- cgn[cgn$FamID == i, ]
  dat.nooc$AffectedOvary <- 0
  dat.nooc$AgeOvary <- dat.nooc$Current.Age
  cgn.boost$P.BRCA.NoOC[i] <- brcapro(dat.nooc, cid, print = FALSE)@probs[1]
  dat.noaj <- cgn[cgn$FamID == i, ]
  dat.noaj$ethnic <- NA
  cgn.boost$P.BRCA.NoAJ[i] <- brcapro(dat.noaj, cid, print = FALSE)@probs[1]
  dat.noocaj <- dat.nooc; dat.noocaj$ethnic <- NA
  cgn.boost$P.BRCA.NoOCAJ[i] <- brcapro(dat.noocaj, cid, print = FALSE)@probs[1]
}

###################################################
## Now adding in other covariate basis functions ##
###################################################

library(data.table)
DT.ov <- DT.pros <- DT.end <- DT.col <- DT.hp <- DT.nfam <- DT.imm <- DT.ngen <- data.table(cgn)
DT.ov <- DT.ov[, ovary(.SD), by = FamID]
DT.pros <- DT.pros[, prostate(.SD), by = FamID]
DT.end <- DT.end[, endometrium(.SD), by = FamID]
DT.col <- DT.col[, colon(.SD), by = FamID]
DT.hp <- DT.hp[, horz.perc(.SD), by = FamID]
DT.nfam <- DT.nfam[, nrow(.SD), by = FamID]
DT.imm <- DT.imm[, bc.imm(.SD), by = FamID]
DT.ngen <- DT.ngen[, num.gen(.SD), by = FamID]
DT.br <- DT.vert <- DT.horz <- data.table(cgn)
DT.br <- DT.br[, breast(.SD), by = FamID]
DT.vert <- DT.vert[, vert(.SD), by = FamID]
DT.horz <- DT.horz[, horz(.SD), by = FamID]

DT.nov <- DT.npros <- DT.nend <- DT.ncol <- DT.age <- DT.agerel <- data.table(cgn)
DT.nov <- DT.nov[, sum(.SD$AffectedOvary[-which(.SD$Relation == 1)]),
                 by = FamID]
DT.npros <- DT.npros[, sum(.SD$AffectedProstate[-which(.SD$Relation == 1)]),
                     by = FamID]
DT.nend <- DT.nend[, sum(.SD$AffectedEndometrium[-which(.SD$Relation == 1)]),
                   by = FamID]
DT.ncol <- DT.ncol[, sum(.SD$AffectedColon[-which(.SD$Relation == 1)]),
                   by = FamID]
DT.age <- DT.age[, ifelse(.SD$Current.Age[.SD$Relation == 1] == 1, 55,
                          .SD$Current.Age[.SD$Relation == 1]), by = FamID]
DT.agerel <- DT.agerel[, mean(.SD$Current.Age[!(.SD$Current.Age %in% 0:1)],
                              na.rm = TRUE),
                       by = FamID]
DT.hp2 <- data.table(cgn)
DT.hp2 <- DT.hp2[, horz.perc2(.SD), by = FamID]


res.vh <- data.frame(FamID = DT.ov$FamID,
                     AffectedOvary = DT.ov$V1, AffectedProstate = DT.pros$V1,
                     AffectedEndometrium = DT.end$V1, AffectedColon = DT.col$V1,
                     HorzPerc = DT.hp$V1, FamSize = DT.nfam$V1,
                     NumBCImmed = DT.imm$V1, NumGen = DT.ngen$V1,
                     NumOvary = DT.nov$V1, NumProstate = DT.npros$V1,
                     NumEndometrium = DT.nend$V1, NumColon = DT.ncol$V1,
                     AgeProband = DT.age$V1, AvgAgeRel = DT.agerel$V1,
                     AffectedBreast = DT.br$V1, Vert = DT.vert$V1,
                     HorzPerc = DT.hp2$V1)

# cgn.boost2 <- merge(CGN.res.boost, res.vh, by = c("Center.ID", "File.ID"))

cgn.boost2 <- merge(cgn.boost, res.vh, by = "FamID")
cgn.boost2 <- merge(cgn.boost2, FamLevelVars[, c("FamID", "BRCAcarrier")],
                    by = "FamID")


library(gbm)



bern.loss <- function(gamma, h, i, dat){
  ybar <- 2 * dat$BRCA - 1
  pred <- as.numeric(as.matrix(select(dat, one_of(paste("Prob.BRCA.", i - 1, sep = "")))))
  return(log(1 + exp(-2 * ybar * (pred + gamma * h))))
}

logit.loss <- function(y, l, gamma, h){
  l2 <- l + gamma * h
  -sum(y * l2 - log(1 + exp(l2)))
}
lin.loss <- function(y, l, gamma, h){
  l2 <- l + gamma * h
  sum((y - 1 / (1 + exp(-l2)))^2)
}

M <- 100
alpha.int <- alpha.vert <- alpha.horz <- 0
cgn.boost2$P.BRCA.NoOCAJ <- as.numeric(cgn.boost2$P.BRCA.NoOCAJ)
cgn.boost2$P.BRCA.NoOC <- as.numeric(cgn.boost2$P.BRCA.NoOC)
cgn.boost2$P.BRCA.NoAJ <- as.numeric(cgn.boost2$P.BRCA.NoAJ)
cgn.boost2$P.BRCA <- as.numeric(cgn.boost2$P.BRCA)
cgn.boost2$LO.0 <- log(cgn.boost2$P.BRCA.NoOCAJ / (1 - cgn.boost2$P.BRCA.NoOCAJ))
cgn.boost2 <- rename(cgn.boost2, BRCA = BRCAcarrier)
cgn.boost2 <- merge(cgn.boost2, select(filter(cgn, Relation == 1), FamID, ethnic),
                    by = "FamID")
cgn.boost2$ethnic[cgn.boost2$ethnic == "AJ"] <- 1
cgn.boost2$ethnic[cgn.boost2$ethnic == "nonAJ"] <- 0
train <- cgn.boost2[1:1000, ]
test <- cgn.boost2[-(1:1000), ]
shrink <- 0.1
bag <- 0.5
set.seed(1)
for(i in 1:M){
  smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
  l <- as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))))
  # resid <- as.numeric((p - train$BRCA) / (p * (1 - p)))
  resid <- exp(l[smp]) / (1 + exp(l[smp])) - train$BRCA[smp]
  fit.h <- lm(resid ~ factor(AffectedBreast) + AffectedOvary + FamSize +
                NumBCImmed + NumGen + AgeProband + Vert + HorzPerc + ethnic,
              data = train[smp, ])
  h <- predict(fit.h, newdata = train)
  gamma.star <- optimize(f = logit.loss, y = train$BRCA,
                         l = l, h = h, interval = c(-100, 100))$minimum * shrink
  if(i == 1){
    coeff <- coef(fit.h) * gamma.star
  } else{
    coeff <- coeff + coef(fit.h) * gamma.star
  }
  # alpha.int <- alpha.int + gamma.star * coef(fit.h)[1]
  # alpha.vert <- alpha.vert + gamma.star * coef(fit.h)[2]
  # alpha.horz <- alpha.horz + gamma.star * coef(fit.h)[3]
  train$LO <- 0
  train$LO <- as.numeric(as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))) +
                                     gamma.star * h))
  names(train)[ncol(train)] <- paste(names(train)[ncol(train)], ".", i, sep = "")
}

library(pROC)
p.train <- as.numeric(1 / (1 + exp(-as.matrix(select(train, one_of(paste("LO.", M, sep = "")))))))
sum(p.train) / sum(train$BRCA)
sum(train$P.BRCA.NoOCAJ) / sum(train$BRCA)
auc(train$BRCA ~ p.train)
auc(train$BRCA ~ train$P.BRCA.NoOCAJ)
mean((p.train - train$BRCA)^2)
mean((train$P.BRCA.NoOCAJ - train$BRCA)^2)

fit.test <- fit.h
fit.test$coefficients <- coeff
p.test <- 1 / (1 + exp(-test$LO.0 - predict(fit.test, newdata = test)))

# p.test <- 1 / (1 + exp(-test$P.BRCA.NoOCAJ - as.numeric(data.matrix(cbind(rep(1, nrow(test)),
#                                        as.numeric(test$AffectedBreast == 1),
#                                        as.numeric(test$AffectedBreast == 2),
#                                        test[, c("AffectedOvary", "HorzPerc", "FamSize",
#                                                 "NumBCImmed", "NumGen",
#                                                 "AgeProband", "Vert", "Horz", "ethnic")])) %*%
#                                      matrix(coeff, length(coeff), 1))))

sum(p.test) / sum(test$BRCA)
sum(test$P.BRCA.NoOCAJ) / sum(test$BRCA)
auc(test$BRCA ~ p.test)
auc(test$BRCA ~ test$P.BRCA.NoOCAJ)
mean((p.test - test$BRCA)^2)
mean((test$P.BRCA.NoOCAJ - test$BRCA)^2)



fit <- gbm(BRCAcarrier ~ P.BRCA.NoOCAJ + AffectedBreast + AffectedOvary + AffectedProstate +
             AffectedEndometrium + AffectedColon + FamSize + NumGen + AgeProband +
             Vert + Horz + NumBCImmed + NumOvary + NumProstate + NumEndometrium +
             NumColon,
           data = cgn.boost2[1:1000, ],
           distribution = "bernoulli",
           shrinkage = 0.1, n.trees = 1000)

z <- 1 / (1 + exp(-predict(fit, newdata = cgn.boost2[-(1:1000), ], n.trees = 10)))

sum(z) / sum(cgn.boost2$BRCAcarrier[-(1:1000)])
sum(cgn.boost2$P.BRCA.NoOCAJ[-(1:1000)]) / sum(cgn.boost2$BRCAcarrier[-(1:1000)])
auc(cgn.boost2$BRCAcarrier[-(1:1000)] ~ z)
auc(BRCAcarrier ~ P.BRCA.NoOCAJ, data = cgn.boost2[-(1:1000), ])
mean((z - cgn.boost2$BRCAcarrier[-(1:1000)])^2)
mean((cgn.boost2$P.BRCA.NoOCAJ[-(1:1000)] - cgn.boost2$BRCAcarrier[-(1:1000)])^2)


## xgboost
library(xgboost)

train <- xgb.DMatrix(data = data.matrix(cgn.boost2[1:1000, ]),
                     label = cgn.boost2$P.BRCA.NoOCAJ[1:1000])
test <- xgb.DMatrix(data = data.matrix(cgn.boost2[-(1:1000), ]),
                    label = cgn.boost2$P.BRCA.NoOCAJ[-(1:1000)])

train <- sparse.model.matrix(BRCAcarrier ~ AffectedOvary + AffectedProstate +
                               AffectedEndometrium + AffectedColon +
                               HorzPerc + FamSize + NumGen + AgeProband +
                               AffectedBreast + Vert + Horz,
                             data = cgn.boost2[1:1000, ])
test <- sparse.model.matrix(BRCAcarrier ~ AffectedOvary + AffectedProstate +
                              AffectedEndometrium + AffectedColon +
                              HorzPerc + FamSize + NumGen + AgeProband +
                              AffectedBreast + Vert + Horz,
                            data = cgn.boost2[-(1:1000), ])

train <- xgb.DMatrix(data = data.matrix(cgn.boost2[1:1000, 6:22]),
                     label = cgn.boost2$BRCAcarrier[1:1000])
test <- xgb.DMatrix(data = data.matrix(cgn.boost2[-(1:1000), 6:22]),
                    label = cgn.boost2$BRCAcarrier[-(1:1000)])

a <- xgboost(data = train,
             nrounds = 1000,
             params = list(objective = "binary:logistic",
                           subsample = 0.5, eta = 0.01,
                           base_score = cgn.boost2$P.BRCA.NoOCAJ[1:1000],
                           eval_metric = auc))
auc(cgn.boost2$BRCAcarrier[-(1:1000)] ~ predict(a, test))
auc(cgn.boost2$BRCAcarrier[-(1:1000)] ~ cgn.boost2$P.BRCA.NoOCAJ[-(1:1000)])


a <- xgboost(data = train,
             nrounds = 10,
             params = list(objective = "binary:logistic", base_score = cgn.boost2$P.BRCA.NoOCAJ[1:1000]))


