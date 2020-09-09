## Analyzing Gradient Boosting Simulation
## Last updated: September 9, 2020

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(xtable)
library(dplyr)
library(ggplot2)
library(gridExtra)

load("GB_Simulation_Results")


#### Tables in paper #####

### Table in main text ###
mods <- c("MMR.ngc", "XGB.const.50.CE", "XGB.mmr.50.CE",
          "MMR", "XGB.const.50.CEG", "XGB.mmr.50.CEG",
          "MMR.ngc.ocl", "MMR.ocl",
          "platt", "isoreg")
xtable(rbind(xtab.lp.cv[mods, 1:9], xtab.hp.cv[mods, 1:9]))

### Supplementary tables ###

## small sample size
xtable(rbind(xtab.lp.sss.cv[mods, 1:9], xtab.hp.sss.cv[mods, 1:9]))

## number of iterations
xtable(rbind(xtab.lp.cv[paste0("XGB.", as.vector(t(outer(c("mmr", "const"), c("25", "50", "100"), paste, sep = "."))),
                               ".CEG"), 1:9],
             xtab.hp.cv[paste0("XGB.", as.vector(t(outer(c("mmr", "const"), c("25", "50", "100"), paste, sep = "."))),
                               ".CEG"), 1:9]))

## exponent on GC survival function (extent of penetrance mispecification)
xtable(rbind(xtab.lp.cv[paste0("MMR.", c("025", "05", "2", "4")), 1:9],
             xtab.hp.cv[paste0("MMR.", c("025", "05", "2", "4")), 1:9]))

## bootstrap
xtable(rbind(xtab.lp.boot[mods, 1:9], xtab.hp.boot[mods, 1:9]))

## transportability
xtable(xtab.trans.cv[mods, 1:9])

## calibration intercept/slope
xtable(rbind(xtab.lp.cv[mods, 10:15], xtab.hp.cv[mods, 10:15]))

## low gastric cancer
xtable(xtab.lgc.sss.cv[mods, 1:9])


## calibration plot (by deciles)
load("~/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Gradient Boosting/Simulation/Low Penetrance/LSS/simgb_lp_lss_1.RData")
plots.cal <- vector("list", length(mods))
mods.plot <- c("(1) Mend. w/ GC", "(2) GB w/o GC", "(3) GB w/o GC, w/ Mend.",
               "(4) Mend. w/ GC", "(5) GB w/ GC", "(6) GB w/ GC, w/ Mend.",
               "(7) Ocl. Mend. w/o GC", "(8) Ocl. Mend. w/o GC",
               "(9) Platt scaling", "(10) Isotonic reg.")
for(i in 1:length(plots.cal)){
  for(j in 1:10){
    if(j == 1){
      rp <- res.gb[[mods[i]]]$risk$cv[!is.na(res.gb[[mods[i]]]$risk$cv$risk1), c("MMR", "risk1")]
      names(rp) <- c("MMR", "risk")
    } else{
      rp <- rbind(rp, setNames(res.gb[[mods[i]]]$risk$cv[!is.na(res.gb[[mods[i]]]$risk$cv[, paste0("risk", j)]), c("MMR", paste0("risk", j))],
                               c("MMR", "risk")))
    }
  }
  rp <- mutate(rp, quantile = ntile(risk, 10))
  rp.dec <- setNames(data.frame(matrix(0, 10, 4)), c("RiskAvg", "Prop", "CI.lo", "CI.hi"))
  for(k in 1:10){
    rp.dec.k <- filter(rp, quantile == k)
    p <- mean(rp.dec.k$MMR)
    rp.dec$RiskAvg[k] <- mean(rp.dec.k$risk)
    rp.dec$Prop[k] <- p
    rp.dec$CI.lo[k] <- p - 1.96 * sqrt(p * (1 - p) / nrow(rp.dec.k))
    rp.dec$CI.hi[k] <- p + 1.96 * sqrt(p * (1 - p) / nrow(rp.dec.k))
  }
  plots.cal[[i]] <- ggplot(rp.dec, aes(RiskAvg, Prop)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI.lo, ymax = CI.hi)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(title = mods.plot[i], x = "", y = "") +
    theme(plot.title = element_text(hjust = 0.5))
}

grid.arrange(plots.cal[[1]], plots.cal[[2]], plots.cal[[3]],
             plots.cal[[4]], plots.cal[[5]], plots.cal[[6]],
             plots.cal[[7]], plots.cal[[8]], plots.cal[[9]], plots.cal[[10]],
             nrow = 3, ncol = 4, bottom = "Average Risk in Decile",
             left = "Proportion of MMR Carriers")





## plotting the difference between the scaled GC penetrances
library(ggplot2)
load("pen_gb_gastric10.RData")
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
