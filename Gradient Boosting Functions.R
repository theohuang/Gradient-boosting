## Functions used for Gradient Boosting
## Last updated: November 28, 2018

## Checking if a family has 2 relatives in a row with BC
## Outputting the number of "vertical" instances (a relative can be a 
## member of multliple instances, for example if two siblings have BC
## and their mother and maternal grandmother also have BC)
vert <- function(fam, gene){
  num.vert <- 0
  if(gene == "brca"){
    id.can <- fam$ID[fam$isAffBC == 1 | fam$isAffOC == 1]
  } else if(gene == "mmr"){
    id.can <- fam$ID[fam$isAffColorC == 1 | fam$isAffEndomC == 1]
  }
  if(length(id.can) <= 1){
    return(num.vert)
  } else{
    if(gene == "brca"){
      for(i in 1:length(id.can)){
        if(nrow(fam[fam$ID == fam$MotherID[fam$ID == id.can[i]] &
                    (fam$isAffBC == 1 | fam$isAffOC == 1), ]) > 0){
          num.vert <- num.vert + 1
        }
        if(nrow(fam[fam$ID == fam$FatherID[fam$ID == id.can[i]] &
                    (fam$isAffBC == 1 | fam$isAffOC == 1), ]) > 0){
          num.vert <- num.vert + 1
        }
      }
    } else if(gene == "mmr"){
      for(i in 1:length(id.can)){
        if(nrow(fam[fam$ID == fam$MotherID[fam$ID == id.can[i]] &
                    (fam$isAffColorC == 1 | fam$isAffEndomC == 1), ]) > 0){
          num.vert <- num.vert + 1
        }
        if(nrow(fam[fam$ID == fam$FatherID[fam$ID == id.can[i]] &
                    (fam$isAffColorC == 1 | fam$isAffEndomC == 1), ]) > 0){
          num.vert <- num.vert + 1
        }
      }
    }
  }
  return(num.vert)
}

## number of people in the family who have either cancer
## with a sibling (could be half-sibling) who also has either cancer
horz <- function(fam, gene){
  num.horz <- 0
  id.used <- vector()
  if(gene == "brca"){
    id.can <- fam$ID[fam$isAffBC == 1 | fam$isAffOC == 1]
  } else if(gene == "mmr"){
    id.can <- fam$ID[fam$isAffColorC == 1 | fam$isAffEndomC == 1]
  }
  if(length(id.can) <= 1){
    return(num.horz)
  } else{
    if(gene == "brca"){
      for(i in 1:length(id.can)){
        if(id.can[i] %in% id.used) break
        # getting the siblings (could be half-siblings) of id.can[i]
        fam.sib <- fam[((fam$MotherID == fam$MotherID[fam$ID == id.can[i]] & fam$MotherID[fam$ID == id.can[i]] != 0) |
                          (fam$FatherID == fam$FatherID[fam$ID == id.can[i]] & fam$FatherID[fam$ID == id.can[i]] != 0)) &
                         (fam$isAffBC == 1 | fam$isAffOC == 1) & fam$ID != id.can[i], ]
        if(nrow(fam.sib) > 0){
          num.horz <- num.horz + nrow(fam.sib) + 1
          id.used <- c(id.used, id.can[i], fam.sib$ID)
        }
      }
    } else if(gene == "mmr"){
      for(i in 1:length(id.can)){
        if(id.can[i] %in% id.used) break
        # getting the siblings (could be half-siblings) of id.can[i]
        fam.sib <- fam[((fam$MotherID == fam$MotherID[fam$ID == id.can[i]] & fam$MotherID[fam$ID == id.can[i]] != 0) |
                          (fam$FatherID == fam$FatherID[fam$ID == id.can[i]] & fam$FatherID[fam$ID == id.can[i]] != 0)) &
                         (fam$isAffColorC == 1 | fam$isAffEndomC == 1) & fam$ID != id.can[i], ]
        if(nrow(fam.sib) > 0){
          num.horz <- num.horz + nrow(fam.sib) + 1
          id.used <- c(id.used, id.can[i], fam.sib$ID)
        }
      }
    }
  }
  return(num.horz)
}

can.prop <- function(fam, can){
  colMeans(fam[paste("isAff", can, sep = "")])
}

can.num <- function(fam, can){
  colSums(fam[paste("isAff", can, sep = "")])
}

homoz.genes <- function(pen, pwr){
  cdf <- 1 - (1 - cumsum(pen))^pwr
  return(c(cdf[1], diff(cdf)))
}

### raising gc penetrance to a power
pen.gc.pwr <- function(pen, pwr){
  pen$fFZ <- apply(pen$fFZ, 2, homoz.genes, pwr = pwr)
  pen$fMZ <- apply(pen$fMZ, 2, homoz.genes, pwr = pwr)
  return(pen)
}


## family size, number of cancers in immediate family,
## age of proband, average age of affected relatives, gender



logit <- function(x){
  log(x / (1 - x))
}
expit <- function(x){
  1 / (1 + exp(-x))
}

logit.loss <- function(y, l, gamma, h){
  # y = outcome
  # l = log-odds (for gradient boosting, previous prediction)
  # gamma = parameter to be optimized for gradient boosting
  # h = current prediction for the residual
  l2 <- l + gamma * h
  return(log(1 + exp(l2)) - y * l2)
}

logit.loss2 <- function(y, rho, h){
  l2 <- h + rho
  return(log(1 + exp(l2)) - y * l2)
}

sum.logit.loss <- function(y, l, gamma, h){
  sum(logit.loss(y, l, gamma, h))
}

sum.logit.loss2 <- function(y, rho, h){
  sum(logit.loss2(y, rho, h))
}

resid.loss <- function(y, l){
  return(y - 1 / (1 + exp(-l)))
}


perf.meas <- function(outcome, prediction){
  c(sum(outcome) / sum(prediction),
    pROC::auc(outcome ~ prediction),
    sqrt(mean((outcome - prediction)^2)))
}

oe.xgb <- function(preds, dtrain){
  labels <- getinfo(dtrain, "label")
  err <- sum(labels) / sum(preds) 
  return(list(metric = "OE", value = err))
}

auc.xgb <- function(preds, dtrain){
  labels <- getinfo(dtrain, "label")
  err <- pROC::auc(labels ~ preds)
  return(list(metric = "AUC", value = err))
}

bs.xgb <- function(preds, dtrain){
  labels <- getinfo(dtrain, "label")
  err <- sqrt(mean((labels - preds)^2))
  return(list(metric = "BS", value = err))
}

perf.xgb <- function(preds, dtrain){
  labels <- getinfo(dtrain, "label")
  err <- list(sum(labels) / sum(preds), pROC::auc(labels ~ preds), sqrt(mean((labels - preds)^2)))
  return(list(metric = list("OE", "AUC", "BS"), value = err))
}

## function to run MMRPRO, including options to include gastric cancer
## and to scale the gastric cancer penetrance
mmr.gb.sim <- function(fam, af, CP, gastric = TRUE, scl = 1, pwr = NULL){
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


gb.mmr <- function(dat, shrink, bag, M.mmr, M.const, covs, n.boot, types, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  res.gb <- lapply(setNames(vector("list", length(types) + 2), c(paste("MMR", types, sep = ""), "XGB.mmr", "XGB.const")),
                   function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                    c("OE", "AUC", "BS")), 
                                    risk = setNames(data.frame(cbind(dat$FamID, matrix(NA, nrow(dat), n.boot))),
                                                    c("FamID", paste("risk", 1:n.boot, sep = "")))))
  for(i in 1:n.boot){
    smp.train <- sample(1:nrow(dat), floor(nrow(dat) / 2))
    train <- dat[smp.train, ]
    test <- dat[-smp.train, ]
    
    print(i)
    for(j in 1:length(types)){
      res.gb[[j]]$perf[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
      res.gb[[j]]$risk[dat$FamID %in% test$FamID, i + 1] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
    }
    
    ## XGBoost with MMRPRO
    param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                  objective = "binary:logistic")
    dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                              base_margin = logit(train$P.MMR.ngc))
    dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                             base_margin = logit(test$P.MMR.ngc))
    res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr)
    pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
    res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
    res.gb$XGB.mmr$risk[dat$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
    
    ## XGBoost with constant
    dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
    dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
    res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const)
    pred.xgb.const <- predict(res.xgb.const, dtest.const)
    res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
    res.gb$XGB.const$risk[dat$FamID %in% test$FamID, i + 1] <- pred.xgb.const
  }
  
  return(res.gb)
}


gb.mmr.cv <- function(dat, shrink, bag, M.mmr, M.const, covs, n.boot, types, metric = "auc", nfolds = 5, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  res.gb <- lapply(setNames(vector("list", length(types) + 2), c(paste("MMR", types, sep = ""), "XGB.mmr", "XGB.const")),
                   function(x) list(perf = setNames(data.frame(matrix(0, n.boot, 3)),
                                                    c("OE", "AUC", "BS")), 
                                    risk = setNames(data.frame(cbind(dat$FamID, matrix(NA, nrow(dat), n.boot))),
                                                    c("FamID", paste("risk", 1:n.boot, sep = ""))),
                                    optM = rep(0, n.boot)))
  for(i in 1:n.boot){
    repeat{
      ## we want at least 5 cases in the training set
      smp.train <- sample(1:nrow(dat), floor(nrow(dat) / 2))
      train <- dat[smp.train, ]
      test <- dat[-smp.train, ]
      if(sum(train$MMR) >= 5){
        break
      }
    }
    
    print(i)
    for(j in 1:length(types)){
      res.gb[[j]]$perf[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
      res.gb[[j]]$risk[dat$FamID %in% test$FamID, i + 1] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
    }
    
    ## XGBoost with MMRPRO
    param <- list(max_depth = 2, eta = shrink, silent = 1, subsample = bag,
                  objective = "binary:logistic")
    dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR,
                              base_margin = logit(train$P.MMR.ngc))
    dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR,
                             base_margin = logit(test$P.MMR.ngc))
    # finding the optimal M using 5-fold cross validation and minimizing the metric
    repeat{
      ## we want cases in each "training" set of the cross validation
      optM.mmr <- tryCatch(xgb.cv(param, dtrain.mmr, nrounds = M.mmr, nfold = nfolds,
                                  metrics = metric, verbose = FALSE, early_stopping_rounds = 10)$best_iteration,
                           error = function(e) NA)
      if(!is.na(optM.mmr)){
        break
      }
    }
    res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = optM.mmr)
    pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
    res.gb$XGB.mmr$perf[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
    res.gb$XGB.mmr$risk[dat$FamID %in% test$FamID, i + 1] <- pred.xgb.mmr
    res.gb$XGB.mmr$optM[i] <- optM.mmr
    
    ## XGBoost with constant
    dtrain.const <- xgb.DMatrix(as.matrix(train[, covs]), label = train$MMR)
    dtest.const <- xgb.DMatrix(as.matrix(test[, covs]), label = test$MMR)
    repeat{
      ## we want cases in each "training" set of the cross validation
      optM.const <- tryCatch(xgb.cv(param, dtrain.const, nrounds = M.const, nfold = nfolds,
                                    metrics = metric, verbose = FALSE, early_stopping_rounds = 10)$best_iteration,
                             error = function(e) NA)
      if(!is.na(optM.const)){
        break
      }
    }
    res.xgb.const <- xgb.train(param, dtrain.const, nrounds = optM.const)
    pred.xgb.const <- predict(res.xgb.const, dtest.const)
    res.gb$XGB.const$perf[i, ] <- perf.meas(test$MMR, pred.xgb.const)
    res.gb$XGB.const$risk[dat$FamID %in% test$FamID, i + 1] <- pred.xgb.const
    res.gb$XGB.const$optM[i] <- optM.const
  }
  
  return(res.gb)
}




# gb.mmr <- function(sim.gb, covs, M, shrink, bag, train, test, init = "MMR", seed = NULL){
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#   # smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
#   # train <- sim.gb[smp.train, ]
#   # test <- sim.gb[-smp.train, ]
#   # train$LO.0 <- test$LO.0 <- optimize(f = logit.loss, y = train$MMR,
#   #                                     gamma = 0, h = 0,
#   #                                     interval = c(-100, 100))$minimum
#   if(init == "MMR"){
#     ## using MMRPRO predictions as initial values
#     train$LO.0 <- logit(train$P.MMR)
#     test$LO.0 <- logit(test$P.MMR)
#   } else if(init == "const"){
#     ## using the sample log-odds of the training set as the initial values
#     train$LO.0 <- logit(mean(train$MMR))
#     test$LO.0 <- logit(mean(train$MMR))
#   }
#   
#   auc.gb <- rep(0, M)
#   for(i in 1:M){
#     smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
#     l <- as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))))
#     resid <- resid.loss(train$MMR[smp], l[smp])
#     # resid <- train$MMR[smp] - exp(l[smp]) / (1 + exp(l[smp]))
#     fit.h <- lm(as.formula(paste("resid ~",
#                                  paste(covs, collapse = " + "))),
#                 data = train[smp,])
#     h <- predict(fit.h, newdata = train)
#     gamma.star <- optimize(f = sum.logit.loss, y = train$MMR,
#                            l = l, h = h, interval = c(-100, 100))$minimum * shrink
#     
#     if(i == 1){
#       coeff <- coef(fit.h) * gamma.star
#     } else{
#       coeff <- coeff + coef(fit.h) * gamma.star
#     }
#     
#     train$LO <- 0
#     train$LO <- as.numeric(as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))) +
#                                        gamma.star * h))
#     
#     ## computing AUC for each iteration
#     pred <- 1 / (1 + exp(-train$LO))
#     auc.gb[i] <- pROC::auc(train$MMR ~ pred)
#     
#     names(train)[which(names(train) == "LO")] <- paste("LO.", i, sep = "")
#   }
#   
#   fit.test <- fit.h
#   fit.test$coefficients <- coeff
#   p.test <- 1 / (1 + exp(-test$LO.0 - predict(fit.test, newdata = test)))
#   
#   df <- c(rbind(perf.meas(test$MMR, p.test),
#                 perf.meas(test$MMR, test$P.MMR)))
#   df <- data.frame(matrix(df, 1, 6))
#   names(df) <- c(rbind(paste(c("EO", "AUC", "BS"), ".GB", sep = ""),
#                        paste(c("EO", "AUC", "BS"), ".MMR", sep = "")))
#   return(list(df = df, p.test = p.test, test = test, auc = auc.gb))
# }

# gb.mmr.i <- function(covs, shrink = 0.1, bag = 0.5,
#                      train, test, method = "tree", maxdepth = 2){
#   
#   ## bagging -- randomly selecting a portion of the training set to fit the model
#   smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
#   
#   ## obtaining the residuals and fitting the model
#   resid <- resid.loss(train$MMR[smp], train.lo[smp, i])
#   if(method == "tree"){
#     fit.h[[i]] <- rpart(as.formula(paste("resid ~",
#                                          paste(covs, collapse = " + "))),
#                         data = train[smp,], control = rpart.control(maxdepth = maxdepth))
#   } else if(method == "linreg"){
#     fit.h[[i]] <- lm(as.formula(paste("resid ~",
#                                       paste(covs, collapse = " + "))),
#                      data = train[smp,])
#   }
#   
#   ## getting the fitted values for the whole training set
#   h <- predict(fit.h[[i]], newdata = train)
#   
#   ## getting the optimal gamma parameter
#   gamma[i] <- optimize(f = sum.logit.loss, y = train$MMR,
#                        l = train.lo[, i],
#                        h = h, interval = c(-100, 100))$minimum
#   
#   # ## storing the coefficients of the linear regression
#   # if(method == "linreg"){
#   #   if(i == 1){
#   #     coeff <- coef(fit.h[[i]]) * shrink * gamma[i]
#   #   } else{
#   #     coeff <- coeff + coef(fit.h[[i]]) * shrink * gamma[i]
#   #   }
#   # }
#   
#   ## updating the log-odds of the training set
#   train.lo[, i + 1] <- train.lo[, i] + shrink * gamma[i] * h
#   
#   ## computing performance measures for each iteration
#   pred <- 1 / (1 + exp(-train.lo[, i + 1]))
#   perf.train[i + 1, ] <- perf.meas(train$MMR, pred)
#   
#   ## updating log-odds for testing set
#   resid.test[, i] <- predict(fit.h[[i]], newdata = test)
#   test.lo[, i + 1] <- test.lo[, i] + shrink * gamma[i] * resid.test[, i]
#   return(fit.h, train.lo, test.lo, resid.test, fit.h, gamma)
# }

# gb.mmr <- function(sim.gb, covs, M = 100, shrink = 0.1, bag = 0.5,
#                    train, test, init = "MMR", method = "tree", maxdepth = 2, seed = NULL){
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#   
#   ## keeping track of the inputs of the function
#   inputs <- list(sim.gb = sim.gb, covs = covs, M = M, shrink = shrink,
#                  bag = bag, train = train, test = test, init = init,
#                  method = method, seed = seed)
#   
#   ## initializing the log-odds of the training and testing sets
#   train.lo <- test.lo <- matrix(0, nrow(train), M + 1)
#   if(init == "MMR"){
#     ## using MMRPRO predictions as initial values
#     train.lo[, 1] <- logit(train$P.MMR)
#     test.lo[, 1] <- logit(test$P.MMR)
#   } else if(init == "const"){
#     ## using the sample log-odds of the training set as the initial values
#     train.lo[, 1] <- logit(mean(train$MMR))
#     test.lo[, 1] <- logit(mean(train$MMR))
#   }
#   
#   gamma <- vector()
#   fit.h <- vector("list", M)
#   resid.test <- matrix(0, nrow(test), M)
#   perf.train <- data.frame(matrix(0, M + 1, 3))
#   names(perf.train) <- Cs(EO, AUC, BS)
#   
#   perf.train[1, ] <- perf.meas(train$MMR, 1 / (1 + exp(-train.lo[, 1])))
#   
#   for(i in 1:M){
#     
#     ## bagging -- randomly selecting a portion of the training set to fit the model
#     smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
#     
#     ## obtaining the residuals and fitting the model
#     resid <- resid.loss(train$MMR[smp], train.lo[smp, i])
#     if(method == "tree"){
#       fit.h[[i]] <- rpart(as.formula(paste("resid ~",
#                                            paste(covs, collapse = " + "))),
#                           data = train[smp,], control = rpart.control(maxdepth = maxdepth))
#     } else if(method == "linreg"){
#       fit.h[[i]] <- lm(as.formula(paste("resid ~",
#                                         paste(covs, collapse = " + "))),
#                        data = train[smp,])
#     }
#     
#     ## getting the fitted values for the whole training set
#     h <- predict(fit.h[[i]], newdata = train)
#     
#     ## getting the optimal gamma parameter
#     gamma[i] <- optimize(f = sum.logit.loss, y = train$MMR,
#                          l = train.lo[, i],
#                          h = h, interval = c(-100, 100))$minimum
#     
#     # ## storing the coefficients of the linear regression
#     # if(method == "linreg"){
#     #   if(i == 1){
#     #     coeff <- coef(fit.h[[i]]) * shrink * gamma[i]
#     #   } else{
#     #     coeff <- coeff + coef(fit.h[[i]]) * shrink * gamma[i]
#     #   }
#     # }
#     
#     ## updating the log-odds of the training set
#     train.lo[, i + 1] <- train.lo[, i] + shrink * gamma[i] * h
#     
#     ## computing performance measures for each iteration
#     pred <- 1 / (1 + exp(-train.lo[, i + 1]))
#     perf.train[i + 1, ] <- perf.meas(train$MMR, pred)
#     
#     ## updating log-odds for testing set
#     resid.test[, i] <- predict(fit.h[[i]], newdata = test)
#     test.lo[, i + 1] <- test.lo[, i] + shrink * gamma[i] * resid.test[, i]
#   }
#   
#   ## getting the predictions for the testing set
#   p.test <- 1 / (1 + exp(-test.lo[, M + 1]))
#   
#   # fit.test <- fit.h[[1]]
#   # fit.test$coefficients <- coeff
#   # p.test <- 1 / (1 + exp(-test.lo[, 1] - predict(fit.test, newdata = test)))
#   
#   perf.test <- data.frame(t(setNames(perf.meas(test$MMR, p.test),
#                                      Cs(EO, AUC, BS))))
#   
#   perf.test.mmr <- data.frame(t(setNames(perf.meas(test$MMR, test$P.MMR),
#                                          Cs(EO, AUC, BS))))
#   
#   return(list(inputs = inputs, perf.test = perf.test, perf.train = perf.train,
#               perf.test.mmr = perf.test.mmr, p.test = p.test,
#               train.lo = train.lo, test.lo = test.lo,
#               gamma = gamma, resid.test = resid.test))
# }
# 
# 
# 
# gbm.treeboost <- function(sim.gb, covs, M, shrink, bag, train, test, init = "MMR", seed = NULL){
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#   
#   train.lo <- test.lo <- matrix(0, nrow(train), M + 1)
#   
#   if(init == "MMR"){
#     ## using MMRPRO predictions as initial values
#     train.lo[, 1] <- logit(train$P.MMR)
#     test.lo[, 1] <- logit(test$P.MMR)
#     # train$LO.0 <- logit(train$P.MMR)
#     # test$LO.0 <- logit(test$P.MMR)
#   } else if(init == "const"){
#     ## using the sample log-odds of the training set as the initial values
#     train.lo[, 1] <- logit(mean(train$MMR))
#     test.lo[, 1] <- logit(mean(train$MMR))
#     # train$LO.0 <- logit(mean(train$MMR))
#     # test$LO.0 <- logit(mean(train$MMR))
#   }
#   
#   auc.gb <- rep(0, M)
#   gamma <- fit.h <- nodes <- vector("list", M)
#   for(i in 1:M){
#     smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
#     l <- train.lo[, i]
#     # l <- as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))))
#     resid <- resid.loss(train$MMR[smp], train.lo[smp, i])
#     fit.h[[i]] <- rpart(as.formula(paste("resid ~",
#                                          paste(covs, collapse = " + "))),
#                         data = train[smp,], control = rpart.control(maxdepth = 3))
#     h <- predict(fit.h[[i]], newdata = train)
#     nodes[[i]] <- as.numeric(names(table(fit.h[[i]]$where)))
#     gamma[[i]] <- rep(0, length(nodes[[i]]))
#     for(k in seq_along(nodes[[i]])){
#       gamma[[i]][k] <- optimize(f = sum.logit.loss, y = train$MMR[fit.h[[i]]$where == nodes[[i]][k]],
#                                 l = train.lo[fit.h[[i]]$where == nodes[[i]][k], i],
#                                 h = 1, interval = c(-100, 100))$minimum
#     }
#     
#     gamma.train <- rep(0, nrow(train))
#     gamma.train[smp] <- gamma[[i]][sapply(fit.h[[i]]$where, function(x) which(nodes[[i]] == x))]
#     gamma.train[-smp] <- gamma[[i]][sapply(unlist(lapply(split(train[-smp, ], 1:(nrow(train) - length(smp))),
#                                                          function(x) rpart:::pred.rpart(fit.h[[i]], x))),
#                                            function(x) which(nodes[[i]] == x))]
#     
#     # train$LO <- 0
#     # train$LO <- as.numeric(as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))) +
#     #                                    shrink * gamma.train))
#     
#     train.lo[, i + 1] <- train.lo[, i] + shrink * gamma.train
#     
#     ## computing AUC for each iteration
#     pred <- 1 / (1 + exp(-train.lo[, i + 1]))
#     auc.gb[i] <- pROC::auc(train$MMR ~ pred)
#     
#     # names(train)[which(names(train) == "LO")] <- paste("LO.", i, sep = "")
#   }
#   
#   # lo.test <- matrix(0, nrow(test), M + 1)
#   # lo.test[, 1] <- test$LO.0
#   resid.test <- matrix(0, nrow(test), M)
#   for(i in 1:M){
#     resid.test[, i] <- predict(fit.h[[i]], newdata = test)
#     test.lo[, i + 1] <- test.lo[, i] + shrink *
#       gamma[[i]][sapply(unlist(lapply(split(test, 1:nrow(test)),
#                                       function(x) rpart:::pred.rpart(fit.h[[i]], x))),
#                         function(x) which(nodes[[i]] == x))]
#   }
#   p.test <- 1 / (1 + exp(-test.lo[, M + 1]))
#   
#   perf <- c(rbind(perf.meas(test$MMR, p.test),
#                   perf.meas(test$MMR, test$P.MMR)))
#   perf <- data.frame(matrix(perf, 1, 6))
#   names(perf) <- c(rbind(paste(Cs(EO, AUC, BS), ".GB", sep = ""),
#                          paste(Cs(EO, AUC, BS), ".MMR", sep = "")))
#   return(list(perf = perf, p.test = p.test, train.lo = train.lo,
#               test.lo = test.lo, auc = auc.gb,
#               nodes = nodes, gamma = gamma, resid.test = resid.test))
# }
# 
# 
# 
# ## Trying to replicate the gbm R package
# gbm.r <- function(sim.gb, covs, M, shrink, bag, train,
#                   depth = 1, seed = NULL){
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#   # l <- rep(optimize(f = sum.logit.loss2, y = train$MMR, h = 0,
#   #                   interval = c(-100, 100))$minimum, nrow(train))
#   l <- log(train$P.MMR / (1 - train$P.MMR))
#   for(i in 1:M){
#     smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
#     resid <- resid.loss(train$MMR[smp], l[smp])
#     fit.h <- rpart(as.formula(paste("resid ~",
#                                     paste(covs, collapse = " + "))),
#                    data = train[smp, ], method = "anova",
#                    control = rpart.control(maxdepth = depth))
#     h <- predict(fit.h, newdata = train, type = "vector")
#     ind1 <- which(eval(parse(text = paste("train$", row.names(fit.h$splits)[1], sep = ""))) <
#                     fit.h$splits[1, "index"])
#     rho1 <- optimize(f = sum.logit.loss2, y = train$MMR[ind1],
#                      h = h[ind1], interval = c(-100, 100))$minimum * shrink
#     l[ind1] <- l[ind1] + rho1 * shrink
#     ind2 <- which(eval(parse(text = paste("train$", row.names(fit.h$splits)[1], sep = ""))) >=
#                     fit.h$splits[1, "index"])
#     rho2 <- optimize(f = sum.logit.loss2, y = train$MMR[ind2],
#                      h = h[ind2], interval = c(-100, 100))$minimum * shrink
#     l[ind2] <- l[ind2] + rho2 * shrink
#   }
#   return(1 / (1 + exp(-l)))
# }


