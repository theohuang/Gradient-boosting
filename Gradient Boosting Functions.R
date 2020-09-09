## Functions used for Gradient Boosting
## Last updated: September 9, 2020



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
  logitp <- logit(prediction)
  ind.fin <- !is.infinite(logitp)
  f <- tryCatch(glm(outcome[ind.fin] ~ logitp[ind.fin], family = "binomial"),
                error = function(e) NA)
  f2<-	tryCatch(glm(outcome[ind.fin] ~ offset(logitp[ind.fin]), family = "binomial"),
                error = function(e) NA)
  return(c(tryCatch(sum(outcome) / sum(prediction), error = function(e) NA),
           tryCatch(pROC::auc(outcome ~ prediction), error = function(e) NA),
           sqrt(mean((outcome - prediction)^2)),
           tryCatch(coef(f2)[1], error = function(e) NA),
           tryCatch(coef(f)[2]), error = function(e) NA))
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
mmr.gb.sim <- function(fam, af, CP, gastric = TRUE, scl = c(1, 1, 1), pwr = c(1, 1, 1)){
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

## isotonic regression (code modified from https://www.analyticsvidhya.com/blog/2016/07/platt-scaling-isotonic-regression-minimize-logloss-error/)
fit.isoreg <- function(iso, x0){
  o = iso$o
  if (is.null(o))
    o = 1:length(x)
  x = iso$x[o]
  y = iso$yf
  ind = cut(x0, breaks = unique(x), labels = FALSE, include.lowest = TRUE)
  min.x <- min(x)
  max.x <- max(x)
  adjusted.knots <- iso$iKnots[c(1, which(iso$yf[iso$iKnots] > 0))]
  fits = sapply(seq(along = x0), function(i) {
    j = ind[i]
    
    # Handles the case where unseen data is outside range of the training data
    if (is.na(j)) {
      if (x0[i] > max.x) j <- length(x)
      else if (x0[i] < min.x) j <- 1
    }
    
    # Find the upper and lower parts of the step
    upper.step.n <- min(which(adjusted.knots > j))
    upper.step <- adjusted.knots[upper.step.n]
    lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )
    
    # Perform a liner interpolation between the start and end of the step
    denom <- x[upper.step] - x[lower.step]
    denom <- ifelse(denom == 0, 1, denom)
    val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
    
    # Ensure we bound the probabilities to [0, 1]
    val <- ifelse(val > 1, max.x, val)
    val <- ifelse(val < 0, min.x, val)
    val <- ifelse(is.na(val), max.x, val) # Bit of a hack, NA when at right extreme of distribution
    val
  })
  fits
}


gb.mmr <- function(dat, shrink, bag, M.mmr, M.const, covs, n.cv, types,
                   p.train = 0.5, boot = TRUE, n.boot = 100,
                   dat.train = NULL, dat.test = NULL, cv = TRUE,
                   reportRisk = TRUE, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(!is.null(dat.train) & !is.null(dat.test)){
    n.cv <- 1
    warning("Using n.cv = 1 because user provided dat.train and dat.test")
  }
  covs.vec <- unlist(lapply(covs, function(x) paste(substr(x, 5, 5), collapse = "")))
  perf.metrics <- c("OE", "AUC", "BS", "Intercept", "Slope")
  if(length(types) > 0){
    if(boot == FALSE){
      res.gb <- lapply(setNames(vector("list", length(types) + length(M.mmr) * length(covs) +
                                         length(M.const) * length(covs) + 2),
                                c(paste0("MMR", types), "platt", "isoreg",
                                  paste0("XGB.mmr.", as.vector(outer(M.mmr, covs.vec, paste, sep = "."))),
                                  paste0("XGB.const.", as.vector(outer(M.const, covs.vec, paste, sep = "."))))),
                       function(x) list(perf = list(cv = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics)), 
                                        risk = list(cv = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                  c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))))))
    } else{
      res.gb <- lapply(setNames(vector("list", length(types) + length(M.mmr) * length(covs) +
                                         length(M.const) * length(covs) + 2),
                                c(paste0("MMR", types), "platt", "isoreg",
                                  paste0("XGB.mmr.", as.vector(outer(M.mmr, covs.vec, paste, sep = "."))),
                                  paste0("XGB.const.", as.vector(outer(M.const, covs.vec, paste, sep = "."))))),
                       function(x) list(perf = list(apparent = setNames(rep(0, length(perf.metrics)), perf.metrics),
                                                    boot = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                    boot.orig = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                    boot.est = setNames(rep(0, length(perf.metrics)), perf.metrics),
                                                    cv = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics)), 
                                        risk = list(apparent = setNames(data.frame(cbind(dat$FamID, dat$MMR, NA)), c("FamID", "MMR", "risk")),
                                                    boot = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                    c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))),
                                                    boot.orig = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                         c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))),
                                                    cv = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                  c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))))))
    }
  } else{
    if(boot == FALSE){
      res.gb <- lapply(setNames(vector("list", length(types) + length(M.mmr) * length(covs) +
                                         length(M.const) * length(covs) + 2),
                                c(paste0("MMR", types), "platt", "isoreg",
                                  paste0("XGB.mmr.", as.vector(outer(M.mmr, covs.vec, paste, sep = "."))),
                                  paste0("XGB.const.", as.vector(outer(M.const, covs.vec, paste, sep = "."))))),
                       function(x) list(perf = list(cv = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics)), 
                                        risk = list(cv = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                  c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))))))
    } else{
      res.gb <- lapply(setNames(vector("list", length(types) + length(M.mmr) * length(covs) +
                                         length(M.const) * length(covs) + 2),
                                c(paste0("MMR", types), "platt", "isoreg",
                                  paste0("XGB.mmr.", as.vector(outer(M.mmr, covs.vec, paste, sep = "."))),
                                  paste0("XGB.const.", as.vector(outer(M.const, covs.vec, paste, sep = "."))))),
                       function(x) list(perf = list(apparent = setNames(rep(0, length(perf.metrics)), perf.metrics),
                                                    boot = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                    boot.orig = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                    boot.est = setNames(rep(0, length(perf.metrics)), perf.metrics),
                                                    cv = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics)), 
                                        risk = list(apparent = setNames(data.frame(cbind(dat$FamID, dat$MMR, NA)), c("FamID", "MMR", "risk")),
                                                    boot = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                    c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))),
                                                    boot.orig = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                         c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))),
                                                    cv = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                                  c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))))))
    }
  }
  
  ## XGBoost parameters
  param <- list(max_depth = 2, eta = shrink, verbosity = 0, subsample = bag,
                objective = "binary:logistic")
  
  ####  if training and testing data are separate ####
  if(!is.null(dat.train) & !is.null(dat.test)){
    train <- dat.train
    test <- dat.test
    i <- 1
    
    print(i)
    if(length(types) > 0){
      for(j in 1:length(types)){
        res.gb[[j]]$perf$cv[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
      }
    }
    
    ## MMRpro with Platt scaling
    mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
    pred.platt <- predict(mod.platt, test, type = "response")
    res.gb$platt$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.platt
    res.gb$platt$perf$cv[i, ] <- perf.meas(test$MMR, pred.platt)
    
    ## MMRpro with isotonic regression
    mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
    pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
    res.gb$isoreg$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.ir
    res.gb$isoreg$perf$cv[i, ] <- perf.meas(test$MMR, pred.ir)
    
    ## XGBoost with MMRPRO
    for(j in 1:length(M.mmr)){
      for(k in 1:length(covs)){
        dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                  base_margin = logit(train$P.MMR.ngc))
        dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                 base_margin = logit(test$P.MMR.ngc))
        res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
        pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.mmr
      }
    }
    
    ## XGBoost with constant
    for(j in 1:length(M.const)){
      for(k in 1:length(covs)){
        dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
        dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
        res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
        pred.xgb.const <- predict(res.xgb.const, dtest.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.const
      }
    }
  }
  
  if(boot == TRUE){ ## using bootstrap method from Steyerberg et al. 2001
    
    ## apparent perfromance (training on all data, testing on all data) ##
    train <- test <- dat
    
    if(length(types) > 0){
      for(j in 1:length(types)){
        res.gb[[j]]$perf$apparent <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$risk$apparent$risk <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
      }
    }
    
    ## MMRpro with Platt scaling
    mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
    pred.platt <- predict(mod.platt, test, type = "response")
    res.gb$platt$risk$apparent$risk <- pred.platt
    res.gb$platt$perf$apparent <- perf.meas(test$MMR, pred.platt)
    
    ## MMRpro with isotonic regression
    mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
    pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
    res.gb$isoreg$risk$apparent$risk <- pred.ir
    res.gb$isoreg$perf$apparent <- perf.meas(test$MMR, pred.ir)
    
    ## XGBoost with MMRPRO
    for(j in 1:length(M.mmr)){
      for(k in 1:length(covs)){
        dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                  base_margin = logit(train$P.MMR.ngc))
        dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                 base_margin = logit(test$P.MMR.ngc))
        res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
        pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$apparent <- perf.meas(test$MMR, pred.xgb.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$apparent$risk <- pred.xgb.mmr
      }
    }
    
    ## XGBoost with constant
    for(j in 1:length(M.const)){
      for(k in 1:length(covs)){
        dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
        dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
        res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
        pred.xgb.const <- predict(res.xgb.const, dtest.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$apparent <- perf.meas(test$MMR, pred.xgb.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$apparent$risk <- pred.xgb.const
      }
    }
    
    
    ## bootstrap performance ## 
    test.orig <- dat
    ind.famid <- which(names(dat) == "FamID")
    for(i in 1:n.boot){
      smp.boot <- sample(nrow(dat), replace = TRUE)
      train <- test <- dat[smp.boot, ]
      ind.fam <- !duplicated(test$FamID) ## avoiding duplicate families when outputting the risk (since some families are sampled multiple times durimg the bootstrap)
      
      print(i)
      if(length(types) > 0){
        for(j in 1:length(types)){
          res.gb[[j]]$perf$boot[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
          res.gb[[j]]$perf$boot.orig[i, ] <- perf.meas(test.orig$MMR, test.orig[, which(names(test.orig) == paste("P.MMR", types[j], sep = ""))])
          
          res.gb[[j]]$risk$boot[dat$FamID %in% test$FamID, i + 2] <- test[ind.fam, which(names(test) == paste0("P.MMR", types[j]))]
          res.gb[[j]]$risk$boot.orig[dat$FamID %in% test.orig$FamID, i + 2] <- test.orig[, which(names(test.orig) == paste0("P.MMR", types[j]))]
        }
      }
      
      ## MMRpro with Platt scaling
      mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
      pred.platt <- predict(mod.platt, test, type = "response")
      pred.platt.orig <- predict(mod.platt, test.orig, type = "response")
      
      res.gb$platt$risk$boot[dat$FamID %in% test$FamID, i + 2] <- pred.platt[ind.fam]
      res.gb$platt$perf$boot[i, ] <- perf.meas(test$MMR, pred.platt)
      res.gb$platt$risk$boot.orig[, i + 2] <- pred.platt.orig
      res.gb$platt$perf$boot.orig[i, ] <- perf.meas(test.orig$MMR, pred.platt.orig)
      
      ## MMRpro with isotonic regression
      mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
      pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
      pred.ir.orig <- fit.isoreg(mod.ir, test.orig$P.MMR.ngc)
      res.gb$isoreg$risk$boot[dat$FamID %in% test$FamID, i + 2] <- pred.ir[ind.fam]
      res.gb$isoreg$perf$boot[i, ] <- perf.meas(test$MMR, pred.ir)
      res.gb$isoreg$risk$boot.orig[, i + 2] <- pred.ir.orig
      res.gb$isoreg$perf$boot.orig[i, ] <- perf.meas(test.orig$MMR, pred.ir.orig)
      
      ## XGBoost with MMRPRO
      for(j in 1:length(M.mmr)){
        for(k in 1:length(covs)){
          dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                    base_margin = logit(train$P.MMR.ngc))
          dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                   base_margin = logit(test$P.MMR.ngc))
          res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
          pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$boot[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$boot[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.mmr[ind.fam]
          
          dtest.mmr.orig <- xgb.DMatrix(as.matrix(test.orig[, covs[[k]]]), label = test.orig$MMR,
                                        base_margin = logit(test.orig$P.MMR.ngc))
          pred.xgb.mmr.orig <- predict(res.xgb.mmr, dtest.mmr.orig)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$boot.orig[i, ] <- perf.meas(test.orig$MMR, pred.xgb.mmr.orig)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$boot.orig[, i + 2] <- pred.xgb.mmr.orig
        }
      }
      
      ## XGBoost with constant
      for(j in 1:length(M.const)){
        for(k in 1:length(covs)){
          dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
          dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
          res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
          pred.xgb.const <- predict(res.xgb.const, dtest.const)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$boot[i, ] <- perf.meas(test$MMR, pred.xgb.const)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$boot[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.const[ind.fam]
          
          dtest.const.orig <- xgb.DMatrix(as.matrix(test.orig[, covs[[k]]]), label = test.orig$MMR)
          pred.xgb.const.orig <- predict(res.xgb.const, dtest.const.orig)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$boot.orig[i, ] <- perf.meas(test.orig$MMR, pred.xgb.const.orig)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$boot.orig[, i + 2] <- pred.xgb.const.orig
        }
      }
    }
    
    ## apparent - average(bootstrap - test) ##
    for(j in 1:length(res.gb)){
      res.gb[[j]]$perf$boot.est <- res.gb[[j]]$perf$apparent -
        colMeans(res.gb[[j]]$perf$boot - res.gb[[j]]$perf$boot.orig)
    }
  }
  
  
  ### Monte Carlo cross-validation ###
  if(cv == TRUE){
    for(i in 1:n.cv){
      ## sampling p.train proportion of data for training, and the rest is used for testing
      smp.train <- sample(1:nrow(dat), floor(nrow(dat) * p.train))
      train <- dat[smp.train, ]
      test <- dat[-smp.train, ]
      
      print(i)
      if(length(types) > 0){
        for(j in 1:length(types)){
          res.gb[[j]]$perf$cv[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
          res.gb[[j]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
        }
      }
      
      ## MMRpro with Platt scaling
      mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
      pred.platt <- predict(mod.platt, test, type = "response")
      res.gb$platt$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.platt
      res.gb$platt$perf$cv[i, ] <- perf.meas(test$MMR, pred.platt)
      
      ## MMRpro with isotonic regression
      mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
      pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
      res.gb$isoreg$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.ir
      res.gb$isoreg$perf$cv[i, ] <- perf.meas(test$MMR, pred.ir)
      
      ## XGBoost with MMRPRO
      for(j in 1:length(M.mmr)){
        for(k in 1:length(covs)){
          dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                    base_margin = logit(train$P.MMR.ngc))
          dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                   base_margin = logit(test$P.MMR.ngc))
          res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
          pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
          res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.mmr
        }
      }
      
      ## XGBoost with constant
      for(j in 1:length(M.const)){
        for(k in 1:length(covs)){
          dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
          dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
          res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
          pred.xgb.const <- predict(res.xgb.const, dtest.const)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.const)
          res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.const
        }
      }
    }
  }
  
  ## only reporting the average of the CV risks to reduce the size of res.gb
  if(reportRisk == FALSE){
    for(j in 1:length(res.gb)){
      if(cv == TRUE){
        res.gb[[j]]$risk <- list(cv = cbind(MMR = dat$MMR, risk = rowMeans(res.gb[[j]]$risk$cv[, -1], na.rm = TRUE)))
      } else{
        res.gb[[j]]$risk <- list(cv = cbind(MMR = dat$MMR, risk = res.gb[[j]]$risk$cv[, -1]))
      }
    }
  }
  
  return(res.gb)
}


gb.mmr.gc <- function(dat, shrink, bag, M.mmr, M.const, covs, n.cv, types,
                      p.train = 0.5, dat.train = NULL, dat.test = NULL,
                      reportRisk = TRUE, seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(!is.null(dat.train) & !is.null(dat.test)){
    n.cv <- 1
    warning("Using n.cv = 1 because user provided dat.train and dat.test")
  }
  covs.vec <- unlist(lapply(covs, function(x) paste(substr(x, 5, 5), collapse = "")))
  perf.metrics <- c("OE", "AUC", "BS", "Intercept", "Slope")
  res.gb <- lapply(setNames(vector("list", length(types) + length(M.mmr) * length(covs) +
                                     length(M.const) * length(covs) + 2),
                            c(paste0("MMR", types), "platt", "isoreg",
                              paste0("XGB.mmr.", as.vector(outer(M.mmr, covs.vec, paste, sep = "."))),
                              paste0("XGB.const.", as.vector(outer(M.const, covs.vec, paste, sep = "."))))),
                   function(x) list(perf = list(cv = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                cv.gc = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics),
                                                cv.ngc = setNames(data.frame(matrix(0, n.cv, length(perf.metrics))), perf.metrics)), 
                                    risk = list(cv = setNames(data.frame(cbind(dat$FamID, dat$MMR, matrix(NA, nrow(dat), n.cv))),
                                                              c("FamID", "MMR", paste("risk", 1:n.cv, sep = ""))))))
  
  
  ## XGBoost parameters
  param <- list(max_depth = 2, eta = shrink, verbosity = 0, subsample = bag,
                objective = "binary:logistic")
  
  ####  if training and testing data are separate ####
  if(!is.null(dat.train) & !is.null(dat.test)){
    train <- dat.train
    test <- dat.test
    i <- 1
    
    print(i)
    if(length(types) > 0){
      for(j in 1:length(types)){
        res.gb[[j]]$perf$cv[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, filter(test, PropGastC > 0)[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, filter(test, PropGastC == 0)[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
      }
    }
    
    ## MMRpro with Platt scaling
    mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
    pred.platt <- predict(mod.platt, test, type = "response")
    res.gb$platt$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.platt
    res.gb$platt$perf$cv[i, ] <- perf.meas(test$MMR, pred.platt)
    res.gb$platt$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.platt[which(test$PropGastC > 0)])
    res.gb$platt$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.platt[which(test$PropGastC == 0)])
    
    ## MMRpro with isotonic regression
    mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
    pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
    res.gb$isoreg$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.ir
    res.gb$isoreg$perf$cv[i, ] <- perf.meas(test$MMR, pred.ir)
    res.gb$isoreg$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.ir[which(test$PropGastC > 0)])
    res.gb$isoreg$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.ir[which(test$PropGastC == 0)])
    
    ## XGBoost with MMRPRO
    for(j in 1:length(M.mmr)){
      for(k in 1:length(covs)){
        dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                  base_margin = logit(train$P.MMR.ngc))
        dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                 base_margin = logit(test$P.MMR.ngc))
        res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
        pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.xgb.mmr[which(test$PropGastC > 0)])
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.xgb.mmr[which(test$PropGastC == 0)])
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.mmr
      }
    }
    
    ## XGBoost with constant
    for(j in 1:length(M.const)){
      for(k in 1:length(covs)){
        dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
        dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
        res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
        pred.xgb.const <- predict(res.xgb.const, dtest.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.xgb.const[which(test$PropGastC > 0)])
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.xgb.const[which(test$PropGastC == 0)])
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.const
      }
    }
  }
  
  
  ### Monte Carlo cross-validation ###
  for(i in 1:n.cv){
    ## sampling p.train proportion of data for training, and the rest is used for testing
    smp.train <- sample(1:nrow(dat), floor(nrow(dat) * p.train))
    train <- dat[smp.train, ]
    test <- dat[-smp.train, ]
    
    print(i)
    if(length(types) > 0){
      for(j in 1:length(types)){
        res.gb[[j]]$perf$cv[i, ] <- perf.meas(test$MMR, test[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, filter(test, PropGastC > 0)[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, filter(test, PropGastC == 0)[, which(names(test) == paste("P.MMR", types[j], sep = ""))])
        res.gb[[j]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- test[, which(names(test) == paste("P.MMR", types[j], sep = ""))]
      }
    }
    
    ## MMRpro with Platt scaling
    mod.platt <- glm(MMR ~ P.MMR.ngc, family = "binomial", data = train)
    pred.platt <- predict(mod.platt, test, type = "response")
    res.gb$platt$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.platt
    res.gb$platt$perf$cv[i, ] <- perf.meas(test$MMR, pred.platt)
    res.gb$platt$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.platt[which(test$PropGastC > 0)])
    res.gb$platt$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.platt[which(test$PropGastC == 0)])
    
    ## MMRpro with isotonic regression
    mod.ir <- isoreg(train$P.MMR.ngc, train$MMR)
    pred.ir <- fit.isoreg(mod.ir, test$P.MMR.ngc)
    res.gb$isoreg$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.ir
    res.gb$isoreg$perf$cv[i, ] <- perf.meas(test$MMR, pred.ir)
    res.gb$isoreg$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.ir[which(test$PropGastC > 0)])
    res.gb$isoreg$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.ir[which(test$PropGastC == 0)])
    
    ## XGBoost with MMRPRO
    for(j in 1:length(M.mmr)){
      for(k in 1:length(covs)){
        dtrain.mmr <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR,
                                  base_margin = logit(train$P.MMR.ngc))
        dtest.mmr <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR,
                                 base_margin = logit(test$P.MMR.ngc))
        res.xgb.mmr <- xgb.train(param, dtrain.mmr, nrounds = M.mmr[j])
        pred.xgb.mmr <- predict(res.xgb.mmr, dtest.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.mmr)
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.xgb.mmr[which(test$PropGastC > 0)])
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.xgb.mmr[which(test$PropGastC == 0)])
        res.gb[[paste0("XGB.mmr.", M.mmr[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.mmr
      }
    }
    
    ## XGBoost with constant
    for(j in 1:length(M.const)){
      for(k in 1:length(covs)){
        dtrain.const <- xgb.DMatrix(as.matrix(train[, covs[[k]]]), label = train$MMR)
        dtest.const <- xgb.DMatrix(as.matrix(test[, covs[[k]]]), label = test$MMR)
        res.xgb.const <- xgb.train(param, dtrain.const, nrounds = M.const[j])
        pred.xgb.const <- predict(res.xgb.const, dtest.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv[i, ] <- perf.meas(test$MMR, pred.xgb.const)
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv.gc[i, ] <- perf.meas(filter(test, PropGastC > 0)$MMR, pred.xgb.const[which(test$PropGastC > 0)])
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$perf$cv.ngc[i, ] <- perf.meas(filter(test, PropGastC == 0)$MMR, pred.xgb.const[which(test$PropGastC == 0)])
        res.gb[[paste0("XGB.const.", M.const[j], ".", covs.vec[k])]]$risk$cv[dat$FamID %in% test$FamID, i + 2] <- pred.xgb.const
      }
    }
  }
  
  ## only reporting the average of the CV risks to reduce the size of res.gb
  if(reportRisk == FALSE){
    for(j in 1:length(res.gb)){
      if(cv == TRUE){
        res.gb[[j]]$risk <- list(cv = cbind(MMR = dat$MMR, risk = rowMeans(res.gb[[j]]$risk$cv[, -1], na.rm = TRUE)))
      } else{
        res.gb[[j]]$risk <- list(cv = cbind(MMR = dat$MMR, risk = res.gb[[j]]$risk$cv[, -1]))
      }
    }
  }
  
  return(res.gb)
}