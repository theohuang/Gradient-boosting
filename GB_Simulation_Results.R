## Getting simulation results
## Last updated: September 9, 2020

rm(list=ls())

perf.metrics <- c("OE", "AUC", "BS", "Intercept", "Slope")

## Low penetrance
res.lp <- res.lp.sss <- lapply(vector("list", 22), function(x) list(boot.est = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                    cv = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                    risk.cv = vector("list", 100)))
for(i in 1:100){
  print(i)
  load(paste(getwd(), "/Gradient Boosting/Simulation/Low Penetrance/simgb_lp_", i, ".RData", sep = ""))
  for(j in 1:length(res.gb)){
    res.lp[[j]]$boot.est[i, ] <- res.gb[[j]]$perf$boot.est
    res.lp[[j]]$cv[i, ] <- colMeans(res.gb[[j]]$perf$cv)
    res.lp[[j]]$risk.cv[[i]] <- res.gb[[j]]$risk$cv
    res.lp.sss[[j]]$boot.est[i, ] <- res.gb.sss[[j]]$perf$boot.est
    res.lp.sss[[j]]$cv[i, ] <- colMeans(res.gb.sss[[j]]$perf$cv)
    res.lp.sss[[j]]$risk.cv[[i]] <- res.gb.sss[[j]]$risk$cv
  }
}
names(res.lp) <- names(res.lp.sss) <- names(res.gb)

## Improvement frequency
res.lp.comp <- vector("list", 2)
res.lp.comp[[1]] <- res.lp.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.lp.comp[[1]]) <- rownames(res.lp.comp[[2]]) <- names(res.lp)
names(res.lp.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.lp.comp$boot.est$OE[i] <- mean(abs(res.lp[[i]]$boot.est$OE - 1) < abs(res.lp$MMR.ngc$boot.est$OE - 1))
  res.lp.comp$boot.est$AUC[i] <- mean(res.lp[[i]]$boot.est$AUC > res.lp$MMR.ngc$boot.est$AUC)
  res.lp.comp$boot.est$BS[i] <- mean(res.lp[[i]]$boot.est$BS < res.lp$MMR.ngc$boot.est$BS)
  res.lp.comp$boot.est$Intercept[i] <- mean(abs(res.lp[[i]]$boot.est$Intercept) < abs(res.lp$MMR.ngc$boot.est$Intercept))
  res.lp.comp$boot.est$Slope[i] <- mean(abs(res.lp[[i]]$boot.est$Slope - 1) < abs(res.lp$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.lp.comp$cv$OE[i] <- mean(abs(res.lp[[i]]$cv$OE - 1) < abs(res.lp$MMR.ngc$cv$OE - 1))
  res.lp.comp$cv$AUC[i] <- mean(res.lp[[i]]$cv$AUC > res.lp$MMR.ngc$cv$AUC)
  res.lp.comp$cv$BS[i] <- mean(res.lp[[i]]$cv$BS < res.lp$MMR.ngc$cv$BS)
  res.lp.comp$cv$Intercept[i] <- mean(abs(res.lp[[i]]$cv$Intercept) < abs(res.lp$MMR.ngc$cv$Intercept))
  res.lp.comp$cv$Slope[i] <- mean(abs(res.lp[[i]]$cv$Slope - 1) < abs(res.lp$MMR.ngc$cv$Slope - 1))
}

res.lp.sss.comp <- vector("list", 2)
res.lp.sss.comp[[1]] <- res.lp.sss.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.lp.sss.comp[[1]]) <- rownames(res.lp.sss.comp[[2]]) <- names(res.lp.sss)
names(res.lp.sss.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.lp.sss.comp$boot.est$OE[i] <- mean(abs(res.lp.sss[[i]]$boot.est$OE - 1) < abs(res.lp.sss$MMR.ngc$boot.est$OE - 1))
  res.lp.sss.comp$boot.est$AUC[i] <- mean(res.lp.sss[[i]]$boot.est$AUC > res.lp.sss$MMR.ngc$boot.est$AUC)
  res.lp.sss.comp$boot.est$BS[i] <- mean(res.lp.sss[[i]]$boot.est$BS < res.lp.sss$MMR.ngc$boot.est$BS)
  res.lp.sss.comp$boot.est$Intercept[i] <- mean(abs(res.lp.sss[[i]]$boot.est$Intercept) < abs(res.lp.sss$MMR.ngc$boot.est$Intercept))
  res.lp.sss.comp$boot.est$Slope[i] <- mean(abs(res.lp.sss[[i]]$boot.est$Slope - 1) < abs(res.lp.sss$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.lp.sss.comp$cv$OE[i] <- mean(abs(res.lp.sss[[i]]$cv$OE - 1) < abs(res.lp.sss$MMR.ngc$cv$OE - 1))
  res.lp.sss.comp$cv$AUC[i] <- mean(res.lp.sss[[i]]$cv$AUC > res.lp.sss$MMR.ngc$cv$AUC)
  res.lp.sss.comp$cv$BS[i] <- mean(res.lp.sss[[i]]$cv$BS < res.lp.sss$MMR.ngc$cv$BS)
  res.lp.sss.comp$cv$Intercept[i] <- mean(abs(res.lp.sss[[i]]$cv$Intercept) < abs(res.lp.sss$MMR.ngc$cv$Intercept))
  res.lp.sss.comp$cv$Slope[i] <- mean(abs(res.lp.sss[[i]]$cv$Slope - 1) < abs(res.lp.sss$MMR.ngc$cv$Slope - 1))
}

## High penetrance
res.hp <- res.hp.sss <- lapply(vector("list", 22), function(x) list(boot.est = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                    cv = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                    risk.cv = vector("list", 100)))
for(i in 1:100){
  print(i)
  load(paste(getwd(), "/Gradient Boosting/Simulation/simgb_", i, ".RData", sep = ""))
  for(j in 1:length(res.gb)){
    res.hp[[j]]$boot.est[i, ] <- res.gb[[j]]$perf$boot.est
    res.hp[[j]]$cv[i, ] <- colMeans(res.gb[[j]]$perf$cv)
    # res.hp[[j]]$risk.cv[[i]] <- res.gb[[j]]$risk$cv
    res.hp.sss[[j]]$boot.est[i, ] <- res.gb.sss[[j]]$perf$boot.est
    res.hp.sss[[j]]$cv[i, ] <- colMeans(res.gb.sss[[j]]$perf$cv)
    # res.hp.sss[[j]]$risk.cv[[i]] <- res.gb.sss[[j]]$risk$cv
  }
}
names(res.hp) <- names(res.hp.sss) <- names(res.gb)

## Improvement frequency
res.hp.comp <- vector("list", 2)
res.hp.comp[[1]] <- res.hp.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.hp.comp[[1]]) <- rownames(res.hp.comp[[2]]) <- names(res.hp)
names(res.hp.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.hp.comp$boot.est$OE[i] <- mean(abs(res.hp[[i]]$boot.est$OE - 1) < abs(res.hp$MMR.ngc$boot.est$OE - 1))
  res.hp.comp$boot.est$AUC[i] <- mean(res.hp[[i]]$boot.est$AUC > res.hp$MMR.ngc$boot.est$AUC)
  res.hp.comp$boot.est$BS[i] <- mean(res.hp[[i]]$boot.est$BS < res.hp$MMR.ngc$boot.est$BS)
  res.hp.comp$boot.est$Intercept[i] <- mean(abs(res.hp[[i]]$boot.est$Intercept) < abs(res.hp$MMR.ngc$boot.est$Intercept))
  res.hp.comp$boot.est$Slope[i] <- mean(abs(res.hp[[i]]$boot.est$Slope - 1) < abs(res.hp$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.hp.comp$cv$OE[i] <- mean(abs(res.hp[[i]]$cv$OE - 1) < abs(res.hp$MMR.ngc$cv$OE - 1))
  res.hp.comp$cv$AUC[i] <- mean(res.hp[[i]]$cv$AUC > res.hp$MMR.ngc$cv$AUC)
  res.hp.comp$cv$BS[i] <- mean(res.hp[[i]]$cv$BS < res.hp$MMR.ngc$cv$BS)
  res.hp.comp$cv$Intercept[i] <- mean(abs(res.hp[[i]]$cv$Intercept) < abs(res.hp$MMR.ngc$cv$Intercept))
  res.hp.comp$cv$Slope[i] <- mean(abs(res.hp[[i]]$cv$Slope - 1) < abs(res.hp$MMR.ngc$cv$Slope - 1))
}

res.hp.sss.comp <- vector("list", 2)
res.hp.sss.comp[[1]] <- res.hp.sss.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.hp.sss.comp[[1]]) <- rownames(res.hp.sss.comp[[2]]) <- names(res.hp.sss)
names(res.hp.sss.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.hp.sss.comp$boot.est$OE[i] <- mean(abs(res.hp.sss[[i]]$boot.est$OE - 1) < abs(res.hp.sss$MMR.ngc$boot.est$OE - 1))
  res.hp.sss.comp$boot.est$AUC[i] <- mean(res.hp.sss[[i]]$boot.est$AUC > res.hp.sss$MMR.ngc$boot.est$AUC)
  res.hp.sss.comp$boot.est$BS[i] <- mean(res.hp.sss[[i]]$boot.est$BS < res.hp.sss$MMR.ngc$boot.est$BS)
  res.hp.sss.comp$boot.est$Intercept[i] <- mean(abs(res.hp.sss[[i]]$boot.est$Intercept) < abs(res.hp.sss$MMR.ngc$boot.est$Intercept))
  res.hp.sss.comp$boot.est$Slope[i] <- mean(abs(res.hp.sss[[i]]$boot.est$Slope - 1) < abs(res.hp.sss$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.hp.sss.comp$cv$OE[i] <- mean(abs(res.hp.sss[[i]]$cv$OE - 1) < abs(res.hp.sss$MMR.ngc$cv$OE - 1))
  res.hp.sss.comp$cv$AUC[i] <- mean(res.hp.sss[[i]]$cv$AUC > res.hp.sss$MMR.ngc$cv$AUC)
  res.hp.sss.comp$cv$BS[i] <- mean(res.hp.sss[[i]]$cv$BS < res.hp.sss$MMR.ngc$cv$BS)
  res.hp.sss.comp$cv$Intercept[i] <- mean(abs(res.hp.sss[[i]]$cv$Intercept) < abs(res.hp.sss$MMR.ngc$cv$Intercept))
  res.hp.sss.comp$cv$Slope[i] <- mean(abs(res.hp.sss[[i]]$cv$Slope - 1) < abs(res.hp.sss$MMR.ngc$cv$Slope - 1))
}

## Assessing transportability (training and testing sets come from different distributions)
res.trans <- lapply(vector("list", 10), function(x) list(cv = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                         risk.cv = vector("list", 100)))
for(i in 1:100){
  print(i)
  load(paste(getwd(), "/Gradient Boosting/Simulation/Transportability/simgb_trans_", i, ".RData", sep = ""))
  for(j in 1:length(res.gb)){
    res.trans[[j]]$cv[i, ] <- colMeans(res.gb[[j]]$perf$cv)
    res.trans[[j]]$risk.cv[[i]] <- res.gb[[j]]$risk$cv
  }
}
names(res.trans) <- names(res.gb)

## Improvement frequency
res.trans.comp <- setNames(data.frame(matrix(0, 10, length(perf.metrics))), perf.metrics)
rownames(res.trans.comp) <- names(res.trans)
for(i in 1:10){
  ## cv
  res.trans.comp$OE[i] <- mean(abs(res.trans[[i]]$cv$OE - 1) < abs(res.trans$MMR.ngc$cv$OE - 1))
  res.trans.comp$AUC[i] <- mean(res.trans[[i]]$cv$AUC > res.trans$MMR.ngc$cv$AUC)
  res.trans.comp$BS[i] <- mean(res.trans[[i]]$cv$BS < res.trans$MMR.ngc$cv$BS)
  res.trans.comp$Intercept[i] <- mean(abs(res.trans[[i]]$cv$Intercept) < abs(res.trans$MMR.ngc$cv$Intercept))
  res.trans.comp$Slope[i] <- mean(abs(res.trans[[i]]$cv$Slope - 1) < abs(res.trans$MMR.ngc$cv$Slope - 1))
}

## Low gastric cancer (not scaled to separate carriers and noncarriers)
res.lgc <- res.lgc.sss <- lapply(vector("list", 22), function(x) list(boot.est = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                      cv = setNames(data.frame(matrix(0, 100, length(perf.metrics))), perf.metrics),
                                                                      risk.cv = vector("list", 100)))
for(i in 1:100){
  print(i)
  load(paste(getwd(), "/Gradient Boosting/Simulation/Low Gastric Cancer/simgb_lgc_", i, ".RData", sep = ""))
  for(j in 1:length(res.gb)){
    res.lgc[[j]]$boot.est[i, ] <- res.gb[[j]]$perf$boot.est
    res.lgc[[j]]$cv[i, ] <- colMeans(res.gb[[j]]$perf$cv)
    res.lgc[[j]]$risk.cv[[i]] <- res.gb[[j]]$risk$cv
    res.lgc.sss[[j]]$boot.est[i, ] <- res.gb.sss[[j]]$perf$boot.est
    res.lgc.sss[[j]]$cv[i, ] <- colMeans(res.gb.sss[[j]]$perf$cv)
    res.lgc.sss[[j]]$risk.cv[[i]] <- res.gb.sss[[j]]$risk$cv
  }
}
names(res.lgc) <- names(res.lgc.sss) <- names(res.gb)

## Improvement frequency
res.lgc.comp <- vector("list", 2)
res.lgc.comp[[1]] <- res.lgc.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.lgc.comp[[1]]) <- rownames(res.lgc.comp[[2]]) <- names(res.lgc)
names(res.lgc.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.lgc.comp$boot.est$OE[i] <- mean(abs(res.lgc[[i]]$boot.est$OE - 1) < abs(res.lgc$MMR.ngc$boot.est$OE - 1))
  res.lgc.comp$boot.est$AUC[i] <- mean(res.lgc[[i]]$boot.est$AUC > res.lgc$MMR.ngc$boot.est$AUC)
  res.lgc.comp$boot.est$BS[i] <- mean(res.lgc[[i]]$boot.est$BS < res.lgc$MMR.ngc$boot.est$BS)
  res.lgc.comp$boot.est$Intercept[i] <- mean(abs(res.lgc[[i]]$boot.est$Intercept) < abs(res.lgc$MMR.ngc$boot.est$Intercept))
  res.lgc.comp$boot.est$Slope[i] <- mean(abs(res.lgc[[i]]$boot.est$Slope - 1) < abs(res.lgc$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.lgc.comp$cv$OE[i] <- mean(abs(res.lgc[[i]]$cv$OE - 1) < abs(res.lgc$MMR.ngc$cv$OE - 1))
  res.lgc.comp$cv$AUC[i] <- mean(res.lgc[[i]]$cv$AUC > res.lgc$MMR.ngc$cv$AUC)
  res.lgc.comp$cv$BS[i] <- mean(res.lgc[[i]]$cv$BS < res.lgc$MMR.ngc$cv$BS)
  res.lgc.comp$cv$Intercept[i] <- mean(abs(res.lgc[[i]]$cv$Intercept) < abs(res.lgc$MMR.ngc$cv$Intercept))
  res.lgc.comp$cv$Slope[i] <- mean(abs(res.lgc[[i]]$cv$Slope - 1) < abs(res.lgc$MMR.ngc$cv$Slope - 1))
}

res.lgc.sss.comp <- vector("list", 2)
res.lgc.sss.comp[[1]] <- res.lgc.sss.comp[[2]] <- setNames(data.frame(matrix(0, 22, length(perf.metrics))), perf.metrics)
rownames(res.lgc.sss.comp[[1]]) <- rownames(res.lgc.sss.comp[[2]]) <- names(res.lgc.sss)
names(res.lgc.sss.comp) <- c("boot.est", "cv")
for(i in 1:22){
  ## bootstrap
  res.lgc.sss.comp$boot.est$OE[i] <- mean(abs(res.lgc.sss[[i]]$boot.est$OE - 1) < abs(res.lgc.sss$MMR.ngc$boot.est$OE - 1))
  res.lgc.sss.comp$boot.est$AUC[i] <- mean(res.lgc.sss[[i]]$boot.est$AUC > res.lgc.sss$MMR.ngc$boot.est$AUC)
  res.lgc.sss.comp$boot.est$BS[i] <- mean(res.lgc.sss[[i]]$boot.est$BS < res.lgc.sss$MMR.ngc$boot.est$BS)
  res.lgc.sss.comp$boot.est$Intercept[i] <- mean(abs(res.lgc.sss[[i]]$boot.est$Intercept) < abs(res.lgc.sss$MMR.ngc$boot.est$Intercept))
  res.lgc.sss.comp$boot.est$Slope[i] <- mean(abs(res.lgc.sss[[i]]$boot.est$Slope - 1) < abs(res.lgc.sss$MMR.ngc$boot.est$Slope - 1))
  
  ## cv
  res.lgc.sss.comp$cv$OE[i] <- mean(abs(res.lgc.sss[[i]]$cv$OE - 1) < abs(res.lgc.sss$MMR.ngc$cv$OE - 1))
  res.lgc.sss.comp$cv$AUC[i] <- mean(res.lgc.sss[[i]]$cv$AUC > res.lgc.sss$MMR.ngc$cv$AUC)
  res.lgc.sss.comp$cv$BS[i] <- mean(res.lgc.sss[[i]]$cv$BS < res.lgc.sss$MMR.ngc$cv$BS)
  res.lgc.sss.comp$cv$Intercept[i] <- mean(abs(res.lgc.sss[[i]]$cv$Intercept) < abs(res.lgc.sss$MMR.ngc$cv$Intercept))
  res.lgc.sss.comp$cv$Slope[i] <- mean(abs(res.lgc.sss[[i]]$cv$Slope - 1) < abs(res.lgc.sss$MMR.ngc$cv$Slope - 1))
}

save(res.lp.comp, res.lp.sss.comp, res.hp.comp, res.hp.sss.comp,
     res.trans.comp, res.lgc.comp, res.lgc.sss.comp,
     file = "GB_Comparisons.RData")

#### Summary tables ####

## Low penetrance ##
res.summ.lp.boot <- matrix(0, length(res.lp), length(perf.metrics) * 3)
for(i in 1:length(res.lp)){
  res.summ.lp.boot[i, seq(1, ncol(res.summ.lp.boot), 3)] <- colMeans(res.lp[[i]]$boot.est, na.rm = TRUE)
  res.summ.lp.boot[i, seq(2, ncol(res.summ.lp.boot), 3)] <- apply(res.lp[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lp.boot[i, seq(3, ncol(res.summ.lp.boot), 3)] <- apply(res.lp[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.lp.cv <- matrix(0, length(res.lp), length(perf.metrics) * 3)
for(i in 1:length(res.lp)){
  res.summ.lp.cv[i, seq(1, ncol(res.summ.lp.boot), 3)] <- colMeans(res.lp[[i]]$cv, na.rm = TRUE)
  res.summ.lp.cv[i, seq(2, ncol(res.summ.lp.boot), 3)] <- apply(res.lp[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lp.cv[i, seq(3, ncol(res.summ.lp.boot), 3)] <- apply(res.lp[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.lp.boot <- mutate(data.frame(res.summ.lp.boot), OE.IF = res.lp.comp$boot.est$OE,
                           AUC.IF = res.lp.comp$boot.est$AUC, BS.IF = res.lp.comp$boot.est$BS,
                           Intercept.IF = res.lp.comp$boot.est$Intercept, Slope.IF = res.lp.comp$boot.est$Slope)
res.summ.lp.cv <- mutate(data.frame(res.summ.lp.cv), OE.IF = res.lp.comp$cv$OE,
                         AUC.IF = res.lp.comp$cv$AUC, BS.IF = res.lp.comp$cv$BS,
                         Intercept.IF = res.lp.comp$cv$Intercept, Slope.IF = res.lp.comp$cv$Slope)
res.summ.lp.boot <- res.summ.lp.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.lp.cv <- res.summ.lp.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]

row.names(res.summ.lp.boot) <- row.names(res.summ.lp.cv) <- names(res.lp)
colnames(res.summ.lp.boot) <- colnames(res.summ.lp.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975", "IF"), paste, sep = ".")))

# small sample size
res.summ.lp.sss.boot <- matrix(0, length(res.lp.sss), length(perf.metrics) * 3)
for(i in 1:length(res.lp.sss)){
  res.summ.lp.sss.boot[i, seq(1, ncol(res.summ.lp.boot), 3)] <- colMeans(res.lp.sss[[i]]$boot.est, na.rm = TRUE)
  res.summ.lp.sss.boot[i, seq(2, ncol(res.summ.lp.boot), 3)] <- apply(res.lp.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lp.sss.boot[i, seq(3, ncol(res.summ.lp.boot), 3)] <- apply(res.lp.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.lp.sss.cv <- matrix(0, length(res.lp.sss), length(perf.metrics) * 3)
for(i in 1:length(res.lp.sss)){
  res.summ.lp.sss.cv[i, seq(1, ncol(res.summ.lp.boot), 3)] <- colMeans(res.lp.sss[[i]]$cv, na.rm = TRUE)
  res.summ.lp.sss.cv[i, seq(2, ncol(res.summ.lp.boot), 3)] <- apply(res.lp.sss[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lp.sss.cv[i, seq(3, ncol(res.summ.lp.boot), 3)] <- apply(res.lp.sss[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.lp.sss.boot <- mutate(data.frame(res.summ.lp.sss.boot), OE.IF = res.lp.sss.comp$boot.est$OE,
                               AUC.IF = res.lp.sss.comp$boot.est$AUC, BS.IF = res.lp.sss.comp$boot.est$BS,
                               Intercept.IF = res.lp.sss.comp$boot.est$Intercept, Slope.IF = res.lp.sss.comp$boot.est$Slope)
res.summ.lp.sss.cv <- mutate(data.frame(res.summ.lp.sss.cv), OE.IF = res.lp.sss.comp$cv$OE,
                             AUC.IF = res.lp.sss.comp$cv$AUC, BS.IF = res.lp.sss.comp$cv$BS,
                             Intercept.IF = res.lp.sss.comp$cv$Intercept, Slope.IF = res.lp.sss.comp$cv$Slope)
res.summ.lp.sss.boot <- res.summ.lp.sss.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.lp.sss.cv <- res.summ.lp.sss.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]


row.names(res.summ.lp.sss.boot) <- row.names(res.summ.lp.sss.cv) <- names(res.lp.sss)
colnames(res.summ.lp.sss.boot) <- colnames(res.summ.lp.sss.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975"), paste, sep = ".")))


## High penetrance ##
res.summ.hp.boot <- matrix(0, length(res.hp), length(perf.metrics) * 3)
for(i in 1:length(res.hp)){
  res.summ.hp.boot[i, seq(1, ncol(res.summ.hp.boot), 3)] <- colMeans(res.hp[[i]]$boot.est, na.rm = TRUE)
  res.summ.hp.boot[i, seq(2, ncol(res.summ.hp.boot), 3)] <- apply(res.hp[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.hp.boot[i, seq(3, ncol(res.summ.hp.boot), 3)] <- apply(res.hp[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.hp.cv <- matrix(0, length(res.hp), length(perf.metrics) * 3)
for(i in 1:length(res.hp)){
  res.summ.hp.cv[i, seq(1, ncol(res.summ.hp.boot), 3)] <- colMeans(res.hp[[i]]$cv, na.rm = TRUE)
  res.summ.hp.cv[i, seq(2, ncol(res.summ.hp.boot), 3)] <- apply(res.hp[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.hp.cv[i, seq(3, ncol(res.summ.hp.boot), 3)] <- apply(res.hp[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.hp.boot <- mutate(data.frame(res.summ.hp.boot), OE.IF = res.hp.comp$boot.est$OE,
                           AUC.IF = res.hp.comp$boot.est$AUC, BS.IF = res.hp.comp$boot.est$BS,
                           Intercept.IF = res.hp.comp$boot.est$Intercept, Slope.IF = res.hp.comp$boot.est$Slope)
res.summ.hp.cv <- mutate(data.frame(res.summ.hp.cv), OE.IF = res.hp.comp$cv$OE,
                         AUC.IF = res.hp.comp$cv$AUC, BS.IF = res.hp.comp$cv$BS,
                         Intercept.IF = res.hp.comp$cv$Intercept, Slope.IF = res.hp.comp$cv$Slope)
res.summ.hp.boot <- res.summ.hp.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.hp.cv <- res.summ.hp.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]


row.names(res.summ.hp.boot) <- row.names(res.summ.hp.cv) <- names(res.hp)
colnames(res.summ.hp.boot) <- colnames(res.summ.hp.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975"), paste, sep = ".")))

# small sample size
res.summ.hp.sss.boot <- matrix(0, length(res.hp.sss), length(perf.metrics) * 3)
for(i in 1:length(res.hp.sss)){
  res.summ.hp.sss.boot[i, seq(1, ncol(res.summ.hp.boot), 3)] <- colMeans(res.hp.sss[[i]]$boot.est, na.rm = TRUE)
  res.summ.hp.sss.boot[i, seq(2, ncol(res.summ.hp.boot), 3)] <- apply(res.hp.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.hp.sss.boot[i, seq(3, ncol(res.summ.hp.boot), 3)] <- apply(res.hp.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.hp.sss.cv <- matrix(0, length(res.hp.sss), length(perf.metrics) * 3)
for(i in 1:length(res.hp.sss)){
  res.summ.hp.sss.cv[i, seq(1, ncol(res.summ.hp.boot), 3)] <- colMeans(res.hp.sss[[i]]$cv, na.rm = TRUE)
  res.summ.hp.sss.cv[i, seq(2, ncol(res.summ.hp.boot), 3)] <- apply(res.hp.sss[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.hp.sss.cv[i, seq(3, ncol(res.summ.hp.boot), 3)] <- apply(res.hp.sss[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.hp.sss.boot <- mutate(data.frame(res.summ.hp.sss.boot), OE.IF = res.hp.sss.comp$boot.est$OE,
                               AUC.IF = res.hp.sss.comp$boot.est$AUC, BS.IF = res.hp.sss.comp$boot.est$BS,
                               Intercept.IF = res.hp.sss.comp$boot.est$Intercept, Slope.IF = res.hp.sss.comp$boot.est$Slope)
res.summ.hp.sss.cv <- mutate(data.frame(res.summ.hp.sss.cv), OE.IF = res.hp.sss.comp$cv$OE,
                             AUC.IF = res.hp.sss.comp$cv$AUC, BS.IF = res.hp.sss.comp$cv$BS,
                             Intercept.IF = res.hp.sss.comp$cv$Intercept, Slope.IF = res.hp.sss.comp$cv$Slope)
res.summ.hp.sss.boot <- res.summ.hp.sss.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.hp.sss.cv <- res.summ.hp.sss.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]


row.names(res.summ.hp.sss.boot) <- row.names(res.summ.hp.sss.cv) <- names(res.hp.sss)
colnames(res.summ.hp.sss.boot) <- colnames(res.summ.hp.sss.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975"), paste, sep = ".")))




## Assessing transportabiliity

res.summ.trans.cv <- matrix(0, length(res.trans), length(perf.metrics) * 3)
for(i in 1:length(res.trans)){
  res.summ.trans.cv[i, seq(1, ncol(res.summ.trans.cv), 3)] <- colMeans(res.trans[[i]]$cv, na.rm = TRUE)
  res.summ.trans.cv[i, seq(2, ncol(res.summ.trans.cv), 3)] <- apply(res.trans[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.trans.cv[i, seq(3, ncol(res.summ.trans.cv), 3)] <- apply(res.trans[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.trans.cv <- mutate(data.frame(res.summ.trans.cv), OE.IF = res.trans.comp$OE,
                         AUC.IF = res.trans.comp$AUC, BS.IF = res.trans.comp$BS,
                         Intercept.IF = res.trans.comp$Intercept, Slope.IF = res.trans.comp$Slope)
res.summ.trans.cv <- res.summ.trans.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]


row.names(res.summ.trans.cv) <- names(res.trans)
colnames(res.summ.trans.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975", "IF"), paste, sep = ".")))

## Low gastric cancer
res.summ.lgc.boot <- matrix(0, length(res.lgc), length(perf.metrics) * 3)
for(i in 1:length(res.lgc)){
  res.summ.lgc.boot[i, seq(1, ncol(res.summ.lgc.boot), 3)] <- colMeans(res.lgc[[i]]$boot.est, na.rm = TRUE)
  res.summ.lgc.boot[i, seq(2, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lgc.boot[i, seq(3, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.lgc.cv <- matrix(0, length(res.lgc), length(perf.metrics) * 3)
for(i in 1:length(res.lgc)){
  res.summ.lgc.cv[i, seq(1, ncol(res.summ.lgc.boot), 3)] <- colMeans(res.lgc[[i]]$cv, na.rm = TRUE)
  res.summ.lgc.cv[i, seq(2, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lgc.cv[i, seq(3, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.lgc.boot <- mutate(data.frame(res.summ.lgc.boot), OE.IF = res.lgc.comp$boot.est$OE,
                           AUC.IF = res.lgc.comp$boot.est$AUC, BS.IF = res.lgc.comp$boot.est$BS,
                           Intercept.IF = res.lgc.comp$boot.est$Intercept, Slope.IF = res.lgc.comp$boot.est$Slope)
res.summ.lgc.cv <- mutate(data.frame(res.summ.lgc.cv), OE.IF = res.lgc.comp$cv$OE,
                         AUC.IF = res.lgc.comp$cv$AUC, BS.IF = res.lgc.comp$cv$BS,
                         Intercept.IF = res.lgc.comp$cv$Intercept, Slope.IF = res.lgc.comp$cv$Slope)
res.summ.lgc.boot <- res.summ.lgc.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.lgc.cv <- res.summ.lgc.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]



row.names(res.summ.lgc.boot) <- row.names(res.summ.lgc.cv) <- names(res.lgc)
colnames(res.summ.lgc.boot) <- colnames(res.summ.lgc.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975"), paste, sep = ".")))

# small sample size
res.summ.lgc.sss.boot <- matrix(0, length(res.lgc.sss), length(perf.metrics) * 3)
for(i in 1:length(res.lgc.sss)){
  res.summ.lgc.sss.boot[i, seq(1, ncol(res.summ.lgc.boot), 3)] <- colMeans(res.lgc.sss[[i]]$boot.est, na.rm = TRUE)
  res.summ.lgc.sss.boot[i, seq(2, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lgc.sss.boot[i, seq(3, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc.sss[[i]]$boot.est, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}
res.summ.lgc.sss.cv <- matrix(0, length(res.lgc.sss), length(perf.metrics) * 3)
for(i in 1:length(res.lgc.sss)){
  res.summ.lgc.sss.cv[i, seq(1, ncol(res.summ.lgc.boot), 3)] <- colMeans(res.lgc.sss[[i]]$cv, na.rm = TRUE)
  res.summ.lgc.sss.cv[i, seq(2, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc.sss[[i]]$cv, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  res.summ.lgc.sss.cv[i, seq(3, ncol(res.summ.lgc.boot), 3)] <- apply(res.lgc.sss[[i]]$cv, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
}

res.summ.lgc.sss.boot <- mutate(data.frame(res.summ.lgc.sss.boot), OE.IF = res.lgc.sss.comp$boot.est$OE,
                           AUC.IF = res.lgc.sss.comp$boot.est$AUC, BS.IF = res.lgc.sss.comp$boot.est$BS,
                           Intercept.IF = res.lgc.sss.comp$boot.est$Intercept, Slope.IF = res.lgc.sss.comp$boot.est$Slope)
res.summ.lgc.sss.cv <- mutate(data.frame(res.summ.lgc.sss.cv), OE.IF = res.lgc.sss.comp$cv$OE,
                         AUC.IF = res.lgc.sss.comp$cv$AUC, BS.IF = res.lgc.sss.comp$cv$BS,
                         Intercept.IF = res.lgc.sss.comp$cv$Intercept, Slope.IF = res.lgc.sss.comp$cv$Slope)
res.summ.lgc.sss.boot <- res.summ.lgc.sss.boot[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]
res.summ.lgc.sss.cv <- res.summ.lgc.sss.cv[, c(1:3, 16, 4:6, 17, 7:9, 18, 10:12, 19, 13:15, 20)]



row.names(res.summ.lgc.sss.boot) <- row.names(res.summ.lgc.sss.cv) <- names(res.lgc.sss)
colnames(res.summ.lgc.sss.boot) <- colnames(res.summ.lgc.sss.cv) <-  as.vector(t(outer(perf.metrics, c("mean", "025", "975"), paste, sep = ".")))



library(xtable)
library(dplyr)

### Organizing in tables ###

### Low penetrance

# cv
xtab.lp.cv <- data.frame(matrix(NA, nrow(res.summ.lp.cv), length(perf.metrics) * 3))
names(xtab.lp.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lp.cv) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.lp.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.lp.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lp.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lp.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lp.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lp.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}


# boot
xtab.lp.boot <- data.frame(matrix(NA, nrow(res.summ.lp.boot), length(perf.metrics) * 3))
names(xtab.lp.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lp.boot) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.lp.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.lp.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lp.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lp.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lp.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lp.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

## small sample size
# cv
xtab.lp.sss.cv <- data.frame(matrix(NA, nrow(res.summ.lp.sss.cv), length(perf.metrics) * 3))
names(xtab.lp.sss.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lp.sss.cv) <- rownames(res.lp.sss.comp$cv)
for(i in 1:nrow(xtab.lp.sss.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.lp.sss.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lp.sss.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lp.sss.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lp.sss.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lp.sss.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}
# boot
xtab.lp.sss.boot <- data.frame(matrix(NA, nrow(res.summ.lp.sss.boot), length(perf.metrics) * 3))
names(xtab.lp.sss.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lp.sss.boot) <- rownames(res.hp.comp$cv)
for(i in 1:nrow(xtab.lp.sss.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.lp.sss.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lp.sss.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lp.sss.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lp.sss.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lp.sss.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

### High penetrance
# cv
xtab.hp.cv <- data.frame(matrix(NA, nrow(res.summ.hp.cv), length(perf.metrics) * 3))
names(xtab.hp.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.hp.cv) <- rownames(res.hp.comp$cv)
for(i in 1:nrow(xtab.hp.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.hp.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.hp.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.hp.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.hp.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.hp.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}
# boot
xtab.hp.boot <- data.frame(matrix(NA, nrow(res.summ.hp.boot), length(perf.metrics) * 3))
names(xtab.hp.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.hp.boot) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.hp.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.hp.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.hp.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.hp.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.hp.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.hp.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

## small sample size
# cv
xtab.hp.sss.cv <- data.frame(matrix(NA, nrow(res.summ.hp.sss.cv), length(perf.metrics) * 3))
names(xtab.hp.sss.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.hp.sss.cv) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.hp.sss.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.hp.sss.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.hp.sss.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.hp.sss.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.hp.sss.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.hp.sss.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}
# boot
xtab.hp.sss.boot <- data.frame(matrix(NA, nrow(res.summ.hp.sss.boot), length(perf.metrics) * 3))
names(xtab.hp.sss.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.hp.sss.boot) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.hp.sss.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.hp.sss.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.hp.sss.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.hp.sss.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.hp.sss.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.hp.sss.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

### Assessing transportability
xtab.trans.cv <- data.frame(matrix(NA, nrow(res.summ.trans.cv), length(perf.metrics) * 3))
names(xtab.trans.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.trans.cv) <- rownames(res.trans.comp)
for(i in 1:nrow(xtab.trans.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.trans.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.trans.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.trans.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.trans.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.trans.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

### Low gastric cancer
# cv
xtab.lgc.cv <- data.frame(matrix(NA, nrow(res.summ.lgc.cv), length(perf.metrics) * 3))
names(xtab.lgc.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lgc.cv) <- rownames(res.lgc.comp$cv)
for(i in 1:nrow(xtab.lgc.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.lgc.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lgc.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lgc.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lgc.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lgc.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}
# boot
xtab.lgc.boot <- data.frame(matrix(NA, nrow(res.summ.lgc.boot), length(perf.metrics) * 3))
names(xtab.lgc.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lgc.boot) <- rownames(res.lp.comp$cv)
for(i in 1:nrow(xtab.lgc.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.lgc.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lgc.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lgc.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lgc.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lgc.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}

## small sample size
# cv
xtab.lgc.sss.cv <- data.frame(matrix(NA, nrow(res.summ.lgc.sss.cv), length(perf.metrics) * 3))
names(xtab.lgc.sss.cv) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lgc.sss.cv) <- rownames(res.lgc.sss.comp$cv)
for(i in 1:nrow(xtab.lgc.sss.cv)){
  for(j in 1:length(perf.metrics)){
    xtab.lgc.sss.cv[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lgc.sss.cv[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lgc.sss.cv[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lgc.sss.cv[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lgc.sss.cv[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}
# boot
xtab.lgc.sss.boot <- data.frame(matrix(NA, nrow(res.summ.lgc.sss.boot), length(perf.metrics) * 3))
names(xtab.lgc.sss.boot) <- apply(expand.grid(c("mean", "CI", "IF"), perf.metrics)[, 2:1], 1, paste0, collapse = ".")
rownames(xtab.lgc.sss.boot) <- rownames(res.lgc.comp$cv)
for(i in 1:nrow(xtab.lgc.sss.boot)){
  for(j in 1:length(perf.metrics)){
    xtab.lgc.sss.boot[i, paste0(perf.metrics[j], c(".mean", ".CI", ".IF"))] <-
      c(format(round(as.numeric(as.character(res.summ.lgc.sss.boot[i, paste0(perf.metrics[j], ".mean")])), 3), nsmall = 3),
        paste0("(", format(round(as.numeric(as.character(res.summ.lgc.sss.boot[i, paste0(perf.metrics[j], ".025")])), 3), nsmall = 3),
               ", ", format(round(as.numeric(as.character(res.summ.lgc.sss.boot[i, paste0(perf.metrics[j], ".975")])), 3), nsmall = 3), ")"),
        round(as.numeric(as.character(res.summ.lgc.sss.boot[i, paste0(perf.metrics[j], ".IF")]))*100))
  }
}



save(res.summ.lp.boot, res.summ.lp.cv, res.summ.lp.sss.boot, res.summ.lp.sss.cv,
     res.summ.hp.boot, res.summ.hp.cv, res.summ.hp.sss.boot, res.summ.hp.sss.cv,
     res.summ.trans.cv,
     res.summ.lgc.boot, res.summ.lgc.cv, res.summ.lgc.sss.boot, res.summ.lgc.sss.cv,
     xtab.lp.boot, xtab.lp.cv, xtab.lp.sss.boot, xtab.lp.sss.cv,
     xtab.hp.boot, xtab.hp.cv, xtab.hp.sss.boot, xtab.hp.sss.cv,
     xtab.trans.cv,
     xtab.lgc.boot, xtab.lgc.cv, xtab.lgc.sss.boot, xtab.lgc.sss.cv,
     file = "GB_Simulation_Results.RData")
