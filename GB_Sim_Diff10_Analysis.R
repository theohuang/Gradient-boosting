## Analyzing Gradient Boosting Simulation
## Last updated: April 4, 2019

rm(list = ls())

res.25 <- vector("list", 8)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", i, ".RData", sep = ""))
  for(j in 1:8){
    res.25[[j]] <- rbind(res.25[[j]], colMeans(res.gb.25[[j]]$perf))
  }
}

res.50 <- vector("list", 8)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", i, ".RData", sep = ""))
  for(j in 1:8){
    res.50[[j]] <- rbind(res.50[[j]], colMeans(res.gb.50[[j]]$perf))
  }
}

res.100 <- vector("list", 8)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", i, ".RData", sep = ""))
  for(j in 1:8){
    res.100[[j]] <- rbind(res.100[[j]], colMeans(res.gb.100[[j]]$perf))
  }
}



res.summ.25 <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.25[i, c(1, 4, 7)] <- colMeans(res.25[[i]])
  res.summ.25[i, c(2, 5, 8)] <- apply(res.25[[i]], 2, function(x) sort(x)[3])
  res.summ.25[i, c(3, 6, 9)] <- apply(res.25[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.25) <- names(res.gb.25)
colnames(res.summ.25) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")

res.summ.50 <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.50[i, c(1, 4, 7)] <- colMeans(res.50[[i]])
  res.summ.50[i, c(2, 5, 8)] <- apply(res.50[[i]], 2, function(x) sort(x)[3])
  res.summ.50[i, c(3, 6, 9)] <- apply(res.50[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.50) <- names(res.gb.50)
colnames(res.summ.50) <- c("OE mean", "OE 3", "OE 97",
                        "AUC mean", "AUC 3", "AUC 97",
                        "BS mean", "BS 3", "BS 97")


res.summ.100 <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ.100[i, c(1, 4, 7)] <- colMeans(res.100[[i]])
  res.summ.100[i, c(2, 5, 8)] <- apply(res.100[[i]], 2, function(x) sort(x)[3])
  res.summ.100[i, c(3, 6, 9)] <- apply(res.100[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ.100) <- names(res.gb.100)
colnames(res.summ.100) <- c("OE mean", "OE 3", "OE 97",
                           "AUC mean", "AUC 3", "AUC 97",
                           "BS mean", "BS 3", "BS 97")

save(res.50, res.25, res.100,
     res.summ.50, res.summ.25, res.summ.100,
     file = "Diff10_Sim_Results.RData")

library(xtable)
library(dplyr)
xtable(res.summ.50, digits = 3)

xtab.25 <- data.frame(round(res.summ.25, 3))
xtab.25$OE <- xtab.25$AUC <- xtab.25$BS <- 0
for(i in 1:8){
  xtab.25$OE[i] <- paste(xtab.25$OE.mean[i], " (", xtab.25[i, 2], ", ", xtab.25[i, 3], ")", sep = "")
  xtab.25$AUC[i] <- paste(xtab.25$AUC.mean[i], " (", xtab.25[i, 5], ", ", xtab.25[i, 6], ")", sep = "")
  xtab.25$BS[i] <- paste(xtab.25$BS.mean[i], " (", xtab.25[i, 8], ", ", xtab.25[i, 9], ")", sep = "")
}
xtab.25 <- xtab.25[, 12:10]

rownames(xtab.25) <- c("MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                       "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.25)

xtab.50 <- data.frame(round(res.summ.50, 3))
xtab.50$OE <- xtab.50$AUC <- xtab.50$BS <- 0
for(i in 1:8){
  xtab.50$OE[i] <- paste(xtab.50$OE.mean[i], " (", xtab.50[i, 2], ", ", xtab.50[i, 3], ")", sep = "")
  xtab.50$AUC[i] <- paste(xtab.50$AUC.mean[i], " (", xtab.50[i, 5], ", ", xtab.50[i, 6], ")", sep = "")
  xtab.50$BS[i] <- paste(xtab.50$BS.mean[i], " (", xtab.50[i, 8], ", ", xtab.50[i, 9], ")", sep = "")
}
xtab.50 <- xtab.50[, 12:10]

rownames(xtab.50) <- c("MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                    "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.50)


xtab.100 <- data.frame(round(res.summ.100, 3))
xtab.100$OE <- xtab.100$AUC <- xtab.100$BS <- 0
for(i in 1:8){
  xtab.100$OE[i] <- paste(xtab.100$OE.mean[i], " (", xtab.100[i, 2], ", ", xtab.100[i, 3], ")", sep = "")
  xtab.100$AUC[i] <- paste(xtab.100$AUC.mean[i], " (", xtab.100[i, 5], ", ", xtab.100[i, 6], ")", sep = "")
  xtab.100$BS[i] <- paste(xtab.100$BS.mean[i], " (", xtab.100[i, 8], ", ", xtab.100[i, 9], ")", sep = "")
}
xtab.100 <- xtab.100[, 12:10]

rownames(xtab.100) <- c("MMR, no GC", "MMR, with GC", "MMR, 0.25", "MMR, 0.5",
                    "MMR, 2", "MMR, 4", "GB, MMR", "GB, Constant")
xtable(xtab.100)


xtab.all <- rbind(xtab.25, xtab.50[7:8, ], xtab.100[7:8, ])
rownames(xtab.all)[7:12] <- c("GB, MMR, 25", "GB, Constant, 25",
                              "GB, MMR, 50", "GB, Constant, 50",
                              "GB, MMR, 100", "GB, Constant, 100")
xtab.all
xtable(xtab.all)

## gradient boosting with 50 iterations (with and without MMRpro) vs
## MMRpro (with and without gastric cancer)
xtable(xtab.50[c(7:8, 1:2), ], digits = 3)

## gradient boosting with 25, 50, and 100 iterations (with and without MMRpro)
xtable(rbind(xtab.25[7:8, ], xtab.50[7:8, ], xtab.100[7:8, ])[c(1, 3, 5, 2, 4, 6), ])

## MMRpro with misspecified penetrances
xtable(xtab.50[3:6, ], digits = 3)

## plotting the difference between the scaled GC penetrances
library(BMmultigene)
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




