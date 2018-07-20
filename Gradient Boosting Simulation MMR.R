## Simulation for MMR genes and cancers related to these genes
## Last updated: July 19, 2018


library(BayesMendel)
library(dplyr)
library(data.table)
library(BayesMendelKnowledgeBase)
library(BMmultigene)
library(pROC)
library(xtable)
library(ggplot2)
library(doParallel)
registerDoParallel(cores = 4)

## generating the penetrance matrix for a given cancer and gender
## creates a matrix where each column is the penetrance of cancer and gender
## for a gene
gen.pen <- function(cancer, genes, gender){
  
  if(gender == "Female"){
    cancer.oppgen <- "ProstC"
  } else{
    cancer.oppgen <- c("OC", "Endom")
  }
  if(cancer %in% cancer.oppgen){
    pen.nc <- rep(0, 94)
  } else{
    for(i in 1:length(genes)){
      if(exists(paste(genes[i], ".", cancer, sep = ""))){
        if(!is.null(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))))){
          if(length(which(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = "")))$gender == gender)) > 0){
            pen.nc <- filter(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))),
                             carrier == "No", (mastectomy == "No" | is.na(mastectomy)),
                             (oophorectomy == "No" | is.na(oophorectomy)),
                             gender == gender)$risk[1:94]
            if(sum(pen.nc) >= 1){
              pen.nc <- pen.nc * 0.99999 / sum(pen.nc)
            }
            break
          }
        }
      }
    }
  }
  if(!exists("pen.nc")){
    pen.nc <- rep(0, 94)
  }
  
  pen <- replicate(length(genes) + 1, pen.nc, cbind())
  
  for(i in 1:length(genes)){
    if(exists(paste(genes[i], ".", cancer, sep = ""))){
      if(!is.null(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))))){
        if(length(which(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = "")))$gender == gender)) > 0){
          pen[, i + 1] <- filter(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))),
                                 carrier == "Yes", (mastectomy == "No" | is.na(mastectomy)),
                                 (oophorectomy == "No" | is.na(oophorectomy)),
                                 gender == gender)$risk[1:94]
          if(sum(pen[, i + 1]) >= 1){
            pen[, i + 1] <- pen[, i + 1] * 0.99999 / sum(pen[, i + 1])
          }
        }
      }
    }
  }
  return(pen)
}

## generating a list of the penetrances
gen.pen.list <- function(cancers, genes, gender){
  pen <- list()
  for(i in 1:length(cancers)){
    pen[[i]] <- gen.pen(cancers[i], genes, gender)
  }
  names(pen) <- cancers
  return(pen)
}

## generating the families
gen.fam <- function(n.sim, CP, af, seed = 1, age.max = 94){
  set.seed(seed)
  fam.sim <- list()
  for(i in 1:n.sim){
    nSibsPatern <- sample(0:3, 2, replace = TRUE)
    nSibsMatern <- sample(0:3, 2, replace = TRUE)
    nSibs <- sample(0:3, 2, replace = TRUE)
    nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs) + 1), replace = TRUE), sum(nSibs) + 1, 2)
    fam.sim[[i]] <- tryCatch(sim.simFam(nSibsPatern, nSibsMatern, nSibs,
                                        nGrandchild, af, CP, includeGeno = TRUE, age.max = age.max),
                             error = function(e) NULL)
  }
  return(fam.sim)
}

## Choosing the cancers, genes, and allele frequencies
genes <- c("MLH1", "MSH2", "MSH6")
# finding the cancers for which we have penetrances related to these genes
cancers <- unique(foreach(i = 1:length(genes), .combine = "c") %do% {
  setdiff(unlist(strsplit(ls(pos = 8, pattern = paste(genes[i], ".*", sep = "")),
                          ".", fixed = TRUE)), genes[i])
})
# some penetrance objects may be NULL. If there is a cancer where we don't
# really have the penetrance, find it and delete it from the cancers list
can.rem <- vector()
for(i in 1:length(cancers)){
  ct <- 0
  for(j in 1:length(genes)){
    if(!exists(paste(genes[j], ".", cancers[i], sep = ""))){
      ct <- ct + 1
    } else{
      if(is.null(eval(parse(text = paste(genes[j], ".", cancers[i], sep = "")))$risk.table)){
        ct <- ct + 1
      }
    }
  }
  if(ct == length(genes)){
    can.rem <- c(can.rem, cancers[i])
  }
}
cancers <- setdiff(cancers, can.rem)

# af <- c(MMRparams()$allef[[1]][2], MMRparams()$allef[[2]][2], MMRparams()$allef[[3]][3])
# using an allele frequency of 0.01 for each gene to increase the number of positive carriers
af <- rep(0.01, 3)
names(af) <- genes

## Getting the penetrances
penCancersF <- gen.pen.list(cancers, genes, "Female")
penCancersM <- gen.pen.list(cancers, genes, "Male")
CP <- genCancerPen(mutations, cancers, penCancersF, penCancersM, maxK = length(genes), age.max = 94)

## increasing gastric, pancreatic, ovarian, and brain cancer penetrances to see if gradient boosting does better
# scaling them so that the lifetime risks are all 0.95
penCancersF2 <- gen.pen.list(cancers, genes, "Female")
penCancersM2 <- gen.pen.list(cancers, genes, "Male")
penCancersF2$GastC[, -1] <- penCancersF2$GastC[, -1] %*% diag(0.95 / colSums(penCancersF2$GastC[, -1]))
penCancersM2$GastC[, -1] <- penCancersM2$GastC[, -1] %*% diag(0.95 / colSums(penCancersM2$GastC[, -1]))
penCancersF2$PancC[, -1] <- penCancersF2$PancC[, -1] %*% diag(0.95 / colSums(penCancersF2$PancC[, -1]))
penCancersM2$PancC[, -1] <- penCancersM2$PancC[, -1] %*% diag(0.95 / colSums(penCancersM2$PancC[, -1]))
penCancersF2$OC[, -1] <- penCancersF2$OC[, -1] %*% diag(0.95 / colSums(penCancersF2$OC[, -1]))
penCancersF2$BrainC[, -1] <- penCancersF2$BrainC[, -1] %*% diag(0.95 / colSums(penCancersF2$BrainC[, -1]))
penCancersM2$BrainC[, -1] <- penCancersM2$BrainC[, -1] %*% diag(0.95 / colSums(penCancersM2$BrainC[, -1]))
CP2 <- genCancerPen(mutations, cancers, penCancersF2, penCancersM2, maxK = length(genes), age.max = 94)


## calculate average family size
famsize <- 0
for(i in 1:length(fam.sim)){
  famsize <- famsize + nrow(fam.sim[[i]])
}
famsize / length(fam.sim)



# source('~/Dropbox (Partners HealthCare)/BM_multigene/R/helpers.R')



start <- Sys.time()
fam.sim <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e4, CP, af, seed = i)
}
print(difftime(Sys.time(), start, units = "secs"))



save(fam.sim, file = "gb.famsim.mmr01.RData")


start <- Sys.time()
fam.sim2 <- foreach(i = 1:10, .combine = append) %dopar% {
  gen.fam(1e3, CP2, af, seed = i)
}
print(difftime(Sys.time(), start, units = "secs"))


# brca.gb <- function(fam){
#   fam$AffectedBreast <- fam$isAffBC
#   fam$AffectedOvary <- fam$isAffOC
#   fam$AgeBreast <- fam$AgeBC
#   fam$AgeOvary <- fam$AgeOC
#   fam$AgeBreastContralateral <- fam$CurAge
#   fam$Twins <- rep(0, nrow(fam))
#   fam$ethnic <- NA
#   fam$Death <- fam$isDead
#   fam$AgeDeath <- fam$CurAge
#   return(unlist(brcapro(fam, counselee.id = fam$ID[fam$isProband == 1],
#                         print = FALSE)@probs)[1:3])
# }

mmr.gb <- function(fam){
  fam$AffectedColon <- fam$isAffColorC
  fam$AffectedEndometrium <- fam$isAffEndomC
  fam$AgeColon <- fam$AgeColorC
  fam$AgeEndometrium <- fam$AgeEndomC
  fam$Twins <- rep(0, nrow(fam))
  fam$Death <- fam$isDead
  fam$AgeDeath <- fam$CurAge
  return(unlist(MMRpro(fam, counselee.id = fam$ID[fam$isProband == 1],
                       params = MMRparams(allef = list(c(0.99, 0.01), c(0.99, 0.01), c(0.99, 0.01))),
                       print = FALSE)@probs))
}

# mmr.gb2 <- function(fam){
#   fam$AffectedColon <- fam$isAffColorC
#   fam$AffectedEndometrium <- fam$isAffEndomC
#   fam$AgeColon <- fam$AgeColorC
#   fam$AgeEndometrium <- fam$AgeEndomC
#   fam$Twins <- rep(0, nrow(fam))
#   fam$Death <- fam$isDead
#   fam$AgeDeath <- fam$CurAge
#   return(unlist(MMRpro(fam, counselee.id = fam$ID[fam$isProband == 1],
#                        print = FALSE)@probs))
# }

# brcammr.gb <- function(fam){
#   c(tryCatch(brca.gb(fam), error = function(e) rep(NA, 3)),
#     tryCatch(mmr.gb(fam), error = function(e) rep(NA, 4)))
# }

start <- Sys.time()
sim.bm <- foreach(i = 1:5000, .combine = rbind) %dopar% {
  tryCatch(mmr.gb(fam.sim[[i]]), error = function(e) rep(NA, 4))
}
print(difftime(Sys.time(), start, units = "secs"))


sim.bm <- data.frame(sim.bm)
sim.bm$FamID <- 1:nrow(sim.bm)
sim.bm <- sim.bm[, c(ncol(sim.bm), 1:(ncol(sim.bm) - 1))]
names(sim.bm) <- c("FamID", "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6")

sim.bm[c("MMR", "MLH1", "MSH2", "MSH6")] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm)){
  if(!is.null(fam.sim[[i]])){
    sim.bm$MLH1[i] <- fam.sim[[i]]$MLH1[fam.sim[[i]]$isProband == 1]
    sim.bm$MSH2[i] <- fam.sim[[i]]$MSH2[fam.sim[[i]]$isProband == 1]
    sim.bm$MSH6[i] <- fam.sim[[i]]$MSH6[fam.sim[[i]]$isProband == 1]
    sim.bm$MMR[i] <- as.numeric(sum(c(sim.bm$MLH1[i], sim.bm$MSH2[i], sim.bm$MSH6[i])) > 0)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

sim.bm["ColorC"] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm)){
  if(!is.null(fam.sim[[i]])){
    sim.bm$ColorC[i] <- fam.sim[[i]]$isAffColorC[fam.sim[[i]]$isProband == 1]
  }
}
print(difftime(Sys.time(), start, units = "secs"))

# sim.bm <- filter(sim.bm, !is.na(P.BRCA))

can.names <- paste("Prop", cancers, sep = "")
sim.bm[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm)){
  if(!is.null(fam.sim[[i]])){
    sim.bm[can.names][sim.bm$FamID == i, ] <- colMeans(fam.sim[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

can.names.num <- paste("Num", cancers, sep = "")
sim.bm[can.names.num] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm)){
  if(!is.null(fam.sim[[i]])){
    sim.bm[can.names.num][sim.bm$FamID == i, ] <- colSums(fam.sim[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))


save(sim.bm, file = "gb.sim.mmr01.RData")



### increasing penetrance

start <- Sys.time()
sim.bm2 <- foreach(i = 1:length(fam.sim2), .combine = rbind) %dopar% {
  tryCatch(mmr.gb(fam.sim2[[i]]), error = function(e) rep(NA, 4))
}
print(difftime(Sys.time(), start, units = "secs"))


sim.bm2 <- data.frame(sim.bm2)
sim.bm2$FamID <- 1:nrow(sim.bm2)
sim.bm2 <- sim.bm2[, c(ncol(sim.bm2), 1:(ncol(sim.bm2) - 1))]
names(sim.bm2) <- c("FamID", "P.MMR", "P.MLH1", "P.MSH2", "P.MSH6")

sim.bm2[c("MMR", "MLH1", "MSH2", "MSH6")] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2$MLH1[i] <- fam.sim2[[i]]$MLH1[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MSH2[i] <- fam.sim2[[i]]$MSH2[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MSH6[i] <- fam.sim2[[i]]$MSH6[fam.sim2[[i]]$isProband == 1]
    sim.bm2$MMR[i] <- as.numeric(sum(c(sim.bm2$MLH1[i], sim.bm2$MSH2[i], sim.bm2$MSH6[i])) > 0)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

sim.bm2["ColorC"] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2$ColorC[i] <- fam.sim2[[i]]$isAffColorC[fam.sim2[[i]]$isProband == 1]
  }
}
print(difftime(Sys.time(), start, units = "secs"))

# sim.bm2 <- filter(sim.bm2, !is.na(P.BRCA))

can.names <- paste("Prop", cancers, sep = "")
sim.bm2[can.names] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2[can.names][sim.bm2$FamID == i, ] <- colMeans(fam.sim2[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))

can.names.num <- paste("Num", cancers, sep = "")
sim.bm2[can.names.num] <- NA
start <- Sys.time()
for(i in 1:nrow(sim.bm2)){
  if(!is.null(fam.sim2[[i]])){
    sim.bm2[can.names.num][sim.bm2$FamID == i, ] <- colSums(fam.sim2[[i]][paste("isAff", cancers, sep = "")])
  }
}
print(difftime(Sys.time(), start, units = "secs"))


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


start <- Sys.time()
DT <- data.table(fam.sim[1:5000])
DT <- DT[, list(lapply(V1, vert, gene = "mmr"), lapply(V1, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm[c("NumVert", "NumHorz")] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm$MSI <- 0
sim.bm$MSI[sim.bm$ColorC == 1 & sim.bm$MMR == 1] <-
  rbinom(length(which(sim.bm$ColorC == 1 & sim.bm$MMR == 1)), 1, 0.968)
sim.bm$MSI[sim.bm$ColorC == 1 & sim.bm$MMR == 0] <-
  rbinom(length(which(sim.bm$ColorC == 1 & sim.bm$MMR == 0)), 1, 1 - 0.914)



## increasing penetrance

start <- Sys.time()
DT <- data.table(fam.sim2)
DT <- DT[, list(lapply(fam.sim2, vert, gene = "mmr"), lapply(fam.sim2, horz, gene = "mmr"))]
print(difftime(Sys.time(), start, units = "secs"))
sim.bm2[c("NumVert", "NumHorz")] <- cbind(unlist(DT$V1), unlist(DT$V2))

sim.bm2$MSI <- 0
sim.bm2$MSI[sim.bm2$ColorC == 1 & sim.bm2$MMR == 1] <-
  rbinom(length(which(sim.bm2$ColorC == 1 & sim.bm2$MMR == 1)), 1, 0.968)
sim.bm2$MSI[sim.bm2$ColorC == 1 & sim.bm2$MMR == 0] <-
  rbinom(length(which(sim.bm2$ColorC == 1 & sim.bm2$MMR == 0)), 1, 1 - 0.914)



logit.loss <- function(y, l, gamma, h){
  l2 <- l + gamma * h
  -sum(y * l2 - log(1 + exp(l2)))
}

M <- 100
shrink <- 0.1
bag <- 0.5

sim.gb <- sim.bm2
sim.gb$LO.0 <- log(sim.gb$P.MMR / (1 - sim.gb$P.MMR))


gb.mmr <- function(sim.gb, covs, M, shrink, bag, seed = 1){
  set.seed(seed)
  smp.train <- sample(1:nrow(sim.gb), floor(nrow(sim.gb) / 2))
  train <- sim.gb[smp.train, ]
  test <- sim.gb[-smp.train, ]
  for(i in 1:M){
    smp <- sort(sample(1:nrow(train), ceiling(bag * nrow(train))))
    l <- as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))))
    resid <- exp(l[smp]) / (1 + exp(l[smp])) - train$MMR[smp]
    fit.h <- lm(as.formula(paste("resid ~",
                                 paste(covs, collapse = " + "))),
                data = train[smp,])
    # fit.h.m <- lm(resid.m ~ NumVert.M + NumHorz.M,
    #               data = train[smp,])
    h <- predict(fit.h, newdata = train)
    gamma.star <- optimize(f = logit.loss, y = train$MMR,
                           l = l, h = h, interval = c(-100, 100))$minimum * shrink
    
    if(i == 1){
      coeff <- coef(fit.h) * gamma.star
    } else{
      coeff <- coeff + coef(fit.h) * gamma.star
    }
    
    train$LO <- 0
    train$LO <- as.numeric(as.matrix(select(train, one_of(paste("LO.", i - 1, sep = ""))) +
                                       gamma.star * h))
    names(train)[c(ncol(train) - 1, ncol(train))] <- paste(names(train)[c(ncol(train) - 1, ncol(train))], ".", i, sep = "")
  }
  
  
  p.train <- as.numeric(1 / (1 + exp(-as.matrix(select(train, one_of(paste("LO.", M, sep = "")))))))
  
  sum(p.train) / sum(train$MMR)
  sum(train$P.MMR) / sum(train$MMR)
  
  auc(train$MMR ~ p.train)
  auc(train$MMR ~ train$P.MMR)
  
  mean((p.train - train$MMR)^2)
  mean((train$P.MMR - train$MMR)^2)
  
  
  fit.test <- fit.h
  fit.test$coefficients <- coeff
  p.test <- 1 / (1 + exp(-test$LO.0 - predict(fit.test, newdata = test)))
  
  
  return(list(df = data.frame(EO.GB = sum(p.test) / sum(test$MMR),
                              EO.MMR = sum(test$P.MMR) / sum(test$MMR),
                              AUC.GB = auc(test$MMR ~ p.test),
                              AUC.MMR = auc(test$MMR ~ test$P.MMR),
                              BS.GB = mean((p.test - test$MMR)^2),
                              BS.MMR = mean((test$P.MMR - test$MMR)^2)),
              p.test = p.test, test = test))
}


start <- Sys.time()
df <- data.frame(matrix(0, 100, 6))
names(df) <- c("EO.GB", "EO.MMR", "AUC.GB", "AUC.MMR", "BS.GB", "BS.MMR")
for(i in 1:100){
  df[i, ] <- gb.mmr(sim.gb, c(can.names[-c(2:3)], "NumVert", "NumHorz"),
                    M, shrink, bag, seed = i)$df
}
colMeans(df)
difftime(Sys.time(), start, units = "secs")



res <- gb.mmr(sim.gb, c(can.names[-c(2:3)], "NumVert", "NumHorz"),
              M, shrink, bag, seed = 1)

ggplot(mutate(res$test, GB = res$p.test), aes(P.MMR, GB)) +
  geom_hex(data = filter(mutate(res$test, GB = res$p.test), MMR == 0),
           bins = 30) +
  geom_point(data = filter(mutate(res$test, GB = res$p.test), MMR == 1),
             aes(color = "red"), size = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "MMRPRO", y = "Gradient Boosting") +
  scale_color_discrete(labels = "MMR carrier", name = "")



# print(xtable(
#   data.frame(EO.GB = sum(p.test) / sum(test$MMR),
#              EO.MMR = sum(test$P.MMR) / sum(test$MMR),
#              AUC.GB = auc(test$MMR ~ p.test),
#              AUC.MMR = auc(test$MMR ~ test$P.MMR),
#              BS.GB = mean((p.test - test$MMR)^2),
#              BS.MMR = mean((test$P.MMR - test$MMR)^2)),
#   digits = 3), include.rownames = FALSE
# )
# 


