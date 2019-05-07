## Analyzing Gradient Boosting Simulation
## Last updated: May 7, 2019

rm(list = ls())

res.25 <- res.50 <- res.100 <- vector("list", 10)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Simulation/simgb_", i, ".RData", sep = ""))
  for(j in 1:10){
    res.25[[j]] <- rbind(res.25[[j]], colMeans(res.gb.25[[j]]$perf))
    res.50[[j]] <- rbind(res.50[[j]], colMeans(res.gb.50[[j]]$perf))
    res.100[[j]] <- rbind(res.100[[j]], colMeans(res.gb.100[[j]]$perf))
  }
}

## Using low-penetrant CRC and EC
res.25.lp <- res.50.lp <- res.100.lp <- vector("list", 10)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Simulation/Low Penetrance/simgb_lp_", i, ".RData", sep = ""))
  for(j in 1:10){
    res.25.lp[[j]] <- rbind(res.25.lp[[j]], colMeans(res.gb.25[[j]]$perf))
    res.50.lp[[j]] <- rbind(res.50.lp[[j]], colMeans(res.gb.50[[j]]$perf))
    res.100.lp[[j]] <- rbind(res.100.lp[[j]], colMeans(res.gb.100[[j]]$perf))
  }
}

## Not using gastric cancer information
res.50.lp.ngc <- res.50.ngc <- vector("list", 2)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Simulation/Low Penetrance/No Gastric Cancer/simgb_lp_ngc_", i, ".RData", sep = ""))
  for(j in 1:2){
    res.50.lp.ngc[[j]] <- rbind(res.50.lp.ngc[[j]], colMeans(res.gb.50.ngc[[j]]$perf))
  }
  load(paste(getwd(), "/Gradient Boosting/Simulation/No Gastric Cancer/simgb_ngc_", i, ".RData", sep = ""))
  for(j in 1:2){
    res.50.ngc[[j]] <- rbind(res.50.ngc[[j]], colMeans(res.gb.50.ngc[[j]]$perf))
  }
}



res.summ.25 <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.25[i, c(1, 4, 7)] <- colMeans(res.25[[i]])
  res.summ.25[i, c(2, 5, 8)] <- apply(res.25[[i]], 2, function(x) sort(x)[3])
  res.summ.25[i, c(3, 6, 9)] <- apply(res.25[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.25) <- names(res.gb.25)
colnames(res.summ.25) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")

res.summ.50 <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.50[i, c(1, 4, 7)] <- colMeans(res.50[[i]])
  res.summ.50[i, c(2, 5, 8)] <- apply(res.50[[i]], 2, function(x) sort(x)[3])
  res.summ.50[i, c(3, 6, 9)] <- apply(res.50[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.50) <- names(res.gb.50)
colnames(res.summ.50) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")


res.summ.100 <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.100[i, c(1, 4, 7)] <- colMeans(res.100[[i]])
  res.summ.100[i, c(2, 5, 8)] <- apply(res.100[[i]], 2, function(x) sort(x)[3])
  res.summ.100[i, c(3, 6, 9)] <- apply(res.100[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.100) <- names(res.gb.100)
colnames(res.summ.100) <- c("OE mean", "OE 3", "OE 97",
                            "AUC mean", "AUC 3", "AUC 97",
                            "BS mean", "BS 3", "BS 97")

res.summ.25.lp <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.25.lp[i, c(1, 4, 7)] <- colMeans(res.25.lp[[i]])
  res.summ.25.lp[i, c(2, 5, 8)] <- apply(res.25.lp[[i]], 2, function(x) sort(x)[3])
  res.summ.25.lp[i, c(3, 6, 9)] <- apply(res.25.lp[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.25.lp) <- names(res.gb.25)
colnames(res.summ.25.lp) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")

res.summ.50.lp <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.50.lp[i, c(1, 4, 7)] <- colMeans(res.50.lp[[i]])
  res.summ.50.lp[i, c(2, 5, 8)] <- apply(res.50.lp[[i]], 2, function(x) sort(x)[3])
  res.summ.50.lp[i, c(3, 6, 9)] <- apply(res.50.lp[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.50.lp) <- names(res.gb.50)
colnames(res.summ.50.lp) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")


res.summ.100.lp <- matrix(0, 10, 9)
for(i in 1:10){
  res.summ.100.lp[i, c(1, 4, 7)] <- colMeans(res.100.lp[[i]])
  res.summ.100.lp[i, c(2, 5, 8)] <- apply(res.100.lp[[i]], 2, function(x) sort(x)[3])
  res.summ.100.lp[i, c(3, 6, 9)] <- apply(res.100.lp[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.100.lp) <- names(res.gb.100)
colnames(res.summ.100.lp) <- c("OE mean", "OE 3", "OE 97",
                            "AUC mean", "AUC 3", "AUC 97",
                            "BS mean", "BS 3", "BS 97")

save(res.50, res.25, res.100,
     res.summ.50, res.summ.25, res.summ.100,
     file = "Diff10_Sim_Results.RData")

save(res.50.lp, res.25.lp, res.100.lp,
     res.summ.50.lp, res.summ.25.lp, res.summ.100.lp,
     file = "Diff10_Sim_Results_lp.RData")


res.summ.50.lp.ngc <- matrix(0, 2, 9)
for(i in 1:2){
  res.summ.50.lp.ngc[i, c(1, 4, 7)] <- colMeans(res.50.lp.ngc[[i]])
  res.summ.50.lp.ngc[i, c(2, 5, 8)] <- apply(res.50.lp.ngc[[i]], 2, function(x) sort(x)[3])
  res.summ.50.lp.ngc[i, c(3, 6, 9)] <- apply(res.50.lp.ngc[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.50.lp.ngc) <- names(res.gb.50.ngc)
colnames(res.summ.50.lp.ngc) <- c("OE mean", "OE 3", "OE 97",
                                 "AUC mean", "AUC 3", "AUC 97",
                                 "BS mean", "BS 3", "BS 97")

res.summ.50.ngc <- matrix(0, 2, 9)
for(i in 1:2){
  res.summ.50.ngc[i, c(1, 4, 7)] <- colMeans(res.50.ngc[[i]])
  res.summ.50.ngc[i, c(2, 5, 8)] <- apply(res.50.ngc[[i]], 2, function(x) sort(x)[3])
  res.summ.50.ngc[i, c(3, 6, 9)] <- apply(res.50.ngc[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.50.ngc) <- names(res.gb.50.ngc)
colnames(res.summ.50.ngc) <- c("OE mean", "OE 3", "OE 97",
                               "AUC mean", "AUC 3", "AUC 97",
                               "BS mean", "BS 3", "BS 97")

save(res.50.lp.ngc, res.50.ngc, res.summ.50.lp.ngc, res.summ.50.ngc,
     file = "GB_Sim_NGC_Results.RData")



load("Diff10_Sim_Results.RData")
load("Diff10_Sim_Results_lp.RData")
load("GB_Sim_NGC_Results.RData")



### Summary tables

library(xtable)
library(dplyr)

xtab.25 <- data.frame(round(res.summ.25, 3))
xtab.25$OE <- xtab.25$AUC <- xtab.25$BS <- 0
for(i in 1:10){
  xtab.25$OE[i] <- paste(xtab.25$OE.mean[i], " (", xtab.25[i, 2], ", ", xtab.25[i, 3], ")", sep = "")
  xtab.25$AUC[i] <- paste(xtab.25$AUC.mean[i], " (", xtab.25[i, 5], ", ", xtab.25[i, 6], ")", sep = "")
  xtab.25$BS[i] <- paste(xtab.25$BS.mean[i], " (", xtab.25[i, 8], ", ", xtab.25[i, 9], ")", sep = "")
}
xtab.25 <- xtab.25[, 12:10]

rownames(xtab.25) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                       "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.25)

xtab.50 <- data.frame(round(rbind(res.summ.50, res.summ.50.ngc), 3))
xtab.50$OE <- xtab.50$AUC <- xtab.50$BS <- 0
for(i in 1:nrow(xtab.50)){
  xtab.50$OE[i] <- paste(xtab.50$OE.mean[i], " (", xtab.50[i, 2], ", ", xtab.50[i, 3], ")", sep = "")
  xtab.50$AUC[i] <- paste(xtab.50$AUC.mean[i], " (", xtab.50[i, 5], ", ", xtab.50[i, 6], ")", sep = "")
  xtab.50$BS[i] <- paste(xtab.50$BS.mean[i], " (", xtab.50[i, 8], ", ", xtab.50[i, 9], ")", sep = "")
}
xtab.50 <- xtab.50[, 12:10]

rownames(xtab.50) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                       "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant",
                       "GB, MMR, no GC", "GB, Constant, no GC")
xtable(xtab.50)


xtab.100 <- data.frame(round(res.summ.100, 3))
xtab.100$OE <- xtab.100$AUC <- xtab.100$BS <- 0
for(i in 1:10){
  xtab.100$OE[i] <- paste(xtab.100$OE.mean[i], " (", xtab.100[i, 2], ", ", xtab.100[i, 3], ")", sep = "")
  xtab.100$AUC[i] <- paste(xtab.100$AUC.mean[i], " (", xtab.100[i, 5], ", ", xtab.100[i, 6], ")", sep = "")
  xtab.100$BS[i] <- paste(xtab.100$BS.mean[i], " (", xtab.100[i, 8], ", ", xtab.100[i, 9], ")", sep = "")
}
xtab.100 <- xtab.100[, 12:10]

rownames(xtab.100) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                        "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.100)


## gradient boosting with 50 iterations (with and without MMRpro) vs
## MMRpro (with and without gastric cancer)
xtable(xtab.50[c(3, 12, 11, 4, 10, 9, 1, 2), ], digits = 3)

## gradient boosting with 25, 50, and 100 iterations (with and without MMRpro)
xtable(rbind(xtab.25[9:10, ], xtab.50[9:10, ], xtab.100[9:10, ])[c(1, 3, 5, 2, 4, 6), ])

## MMRpro with misspecified penetrances
xtable(xtab.50[5:8, ], digits = 3)



### Low penetrance CRC and EC penetrances

xtab.25.lp <- data.frame(round(res.summ.25.lp, 3))
xtab.25.lp$OE <- xtab.25.lp$AUC <- xtab.25.lp$BS <- 0
for(i in 1:10){
  xtab.25.lp$OE[i] <- paste(xtab.25.lp$OE.mean[i], " (", xtab.25.lp[i, 2], ", ", xtab.25.lp[i, 3], ")", sep = "")
  xtab.25.lp$AUC[i] <- paste(xtab.25.lp$AUC.mean[i], " (", xtab.25.lp[i, 5], ", ", xtab.25.lp[i, 6], ")", sep = "")
  xtab.25.lp$BS[i] <- paste(xtab.25.lp$BS.mean[i], " (", xtab.25.lp[i, 8], ", ", xtab.25.lp[i, 9], ")", sep = "")
}
xtab.25.lp <- xtab.25.lp[, 12:10]

rownames(xtab.25.lp) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25.lp", "MMR, 0.5",
                       "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.25.lp)

xtab.50.lp <- data.frame(round(rbind(res.summ.50.lp, res.summ.50.lp.ngc), 3))
xtab.50.lp$OE <- xtab.50.lp$AUC <- xtab.50.lp$BS <- 0
for(i in 1:nrow(xtab.50.lp)){
  xtab.50.lp$OE[i] <- paste(xtab.50.lp$OE.mean[i], " (", xtab.50.lp[i, 2], ", ", xtab.50.lp[i, 3], ")", sep = "")
  xtab.50.lp$AUC[i] <- paste(xtab.50.lp$AUC.mean[i], " (", xtab.50.lp[i, 5], ", ", xtab.50.lp[i, 6], ")", sep = "")
  xtab.50.lp$BS[i] <- paste(xtab.50.lp$BS.mean[i], " (", xtab.50.lp[i, 8], ", ", xtab.50.lp[i, 9], ")", sep = "")
}
xtab.50.lp <- xtab.50.lp[, 12:10]

rownames(xtab.50.lp) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25.lp", "MMR, 0.5",
                       "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant",
                       "GB, MMR, no GC", "GB, Constant, no GC")
xtable(xtab.50.lp)


xtab.100.lp <- data.frame(round(res.summ.100.lp, 3))
xtab.100.lp$OE <- xtab.100.lp$AUC <- xtab.100.lp$BS <- 0
for(i in 1:10){
  xtab.100.lp$OE[i] <- paste(xtab.100.lp$OE.mean[i], " (", xtab.100.lp[i, 2], ", ", xtab.100.lp[i, 3], ")", sep = "")
  xtab.100.lp$AUC[i] <- paste(xtab.100.lp$AUC.mean[i], " (", xtab.100.lp[i, 5], ", ", xtab.100.lp[i, 6], ")", sep = "")
  xtab.100.lp$BS[i] <- paste(xtab.100.lp$BS.mean[i], " (", xtab.100.lp[i, 8], ", ", xtab.100.lp[i, 9], ")", sep = "")
}
xtab.100.lp <- xtab.100.lp[, 12:10]

rownames(xtab.100.lp) <- c("MMR, no GC, oracle", "MMR, with GC, oracle", "MMR, no GC", "MMR, with GC", "MMR, 0.25.lp", "MMR, 0.5",
                        "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.100.lp)


## gradient boosting with 50 iterations (with and without MMRpro) vs
## MMRpro (with and without gastric cancer)
xtable(xtab.50.lp[c(3, 12, 11, 4, 10, 9, 1, 2), ], digits = 3)

## gradient boosting with 25, 50, and 100 iterations (with and without MMRpro)
xtable(rbind(xtab.25.lp[9:10, ], xtab.50.lp[9:10, ], xtab.100.lp[9:10, ])[c(1, 3, 5, 2, 4, 6), ])

## MMRpro with misspecified penetrances
xtable(xtab.50.lp[5:8, ], digits = 3)



xtab.50.lp.ngc <- data.frame(round(res.summ.50.lp.ngc, 3))
xtab.50.lp.ngc$OE <- xtab.50.lp.ngc$AUC <- xtab.50.lp.ngc$BS <- 0
for(i in 1:2){
  xtab.50.lp.ngc$OE[i] <- paste(xtab.50.lp.ngc$OE.mean[i], " (", xtab.50.lp.ngc[i, 2], ", ", xtab.50.lp.ngc[i, 3], ")", sep = "")
  xtab.50.lp.ngc$AUC[i] <- paste(xtab.50.lp.ngc$AUC.mean[i], " (", xtab.50.lp.ngc[i, 5], ", ", xtab.50.lp.ngc[i, 6], ")", sep = "")
  xtab.50.lp.ngc$BS[i] <- paste(xtab.50.lp.ngc$BS.mean[i], " (", xtab.50.lp.ngc[i, 8], ", ", xtab.50.lp.ngc[i, 9], ")", sep = "")
}
xtab.50.lp.ngc <- xtab.50.lp.ngc[, 12:10]

rownames(xtab.50.lp.ngc) <- c("GB, MMR", "GB, Constant")
xtable(xtab.50.lp.ngc)



## plotting the difference between the scaled GC penetrances
# library(BMmultigene)
library(ggplot2)
load("pen_gb_gastric10.RData")
# CP <- genCancerPen(genes, cancers, penCancersF, penCancersM, maxK = length(genes), age.last = 95)
homoz.genes <- function(pen, pwr){
  cdf <- 1 - (1 - cumsum(pen))^pwr
  return(c(cdf[1], diff(cdf)))
}

ggplot(data.frame(Penetrance = c(penCancersM$GastC[, 1],
                                 homoz.genes(penCancersM$GastC[, 1], 0.25),
                                 homoz.genes(penCancersM$GastC[, 1], 0.5),
                                 homoz.genes(penCancersM$GastC[, 1], 2),
                                 homoz.genes(penCancersM$GastC[, 1], 4)),
                  Type = rep(1:5, each = 94),
                  Age = rep(1:94, 5)),
       aes(Age, Penetrance)) +
  geom_point(aes(color = factor(Type))) +
  scale_color_discrete(name = "Power", labels = c("1", "0.25", "0.5", "2", "4")) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 97)) +
  ggtitle("Non-carriers of the MMR mutations") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data.frame(Penetrance = c(penCancersM$GastC[, 2],
                                 homoz.genes(penCancersM$GastC[, 2], 0.25),
                                 homoz.genes(penCancersM$GastC[, 2], 0.5),
                                 homoz.genes(penCancersM$GastC[, 2], 2),
                                 homoz.genes(penCancersM$GastC[, 2], 4)),
                  Type = rep(1:5, each = 94),
                  Age = rep(1:94, 5)),
       aes(Age, Penetrance)) +
  geom_point(aes(color = factor(Type))) +
  scale_color_discrete(name = "Power", labels = c("1", "0.25", "0.5", "2", "4")) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 97)) +
  ggtitle(expression(paste("Carriers of ", italic("MLH1"), "Only", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5))




