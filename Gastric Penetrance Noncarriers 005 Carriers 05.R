## Generating the penetrances for gradient boosting simulations
## Using the MMRPRO penetrances for colorectal and endometrial cancer for
## non-carriers and carriers of one mutation. Using the BayesMendel Knowledge Base
## penetrances for gastric cancer, scaled so that the lifetime risk for non-carriers
## is 0.05 and the lifetime risk for carriers os 0.5.
## Last updated: September 18, 2018

library(BayesMendel)
library(BayesMendelKnowledgeBase)
source('Simulation Functions.R')

cancers <- c("ColorC", "EndomC", "GastC")
genes <- c("MLH1", "MSH2", "MSH6")

penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")

## changing the lifetime risk for gastric cancer to 0.5
lt.risk.nc <- 0.05; lt.risk.c <- 0.5
can.scale <- "GastC"
penCancersF <- pen.scale(penCancersF, can.scale, lt.risk.nc, lt.risk.c)
penCancersM <- pen.scale(penCancersM, can.scale, lt.risk.nc, lt.risk.c)
# using the MMRPRO penetrances for colon and endometrial cancer
penCancersF$ColorC <- penet.mmr.net$fFX[, c("M000", "M100", "M010", "M001")]
penCancersM$ColorC <- penet.mmr.net$fMX[, c("M000", "M100", "M010", "M001")]
penCancersF$EndomC <- penet.mmr.net$fFY[, c("M000", "M100", "M010", "M001")]


save(penCancersF, penCancersM, file = "pen_gb_gastric10.RData")
