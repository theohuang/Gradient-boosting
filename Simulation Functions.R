## Functions used for simulating families
## Last updated: September 10, 2018

## generating the penetrance matrix for a given cancer and gender
## creates a matrix where each column is the penetrance of cancer and gender
## for a gene
gen.pen <- function(cancer, genes, gender, age.max = 94){
  
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
                             gender == gender)$risk[1:age.max]
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
    pen.nc <- rep(0, age.max)
  }
  
  pen <- replicate(length(genes) + 1, pen.nc, cbind())
  
  for(i in 1:length(genes)){
    if(exists(paste(genes[i], ".", cancer, sep = ""))){
      if(!is.null(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))))){
        if(length(which(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = "")))$gender == gender)) > 0){
          pen[, i + 1] <- filter(eval(parse(text = paste(genes[i], ".", cancer, "$risk.table", sep = ""))),
                                 carrier == "Yes", (mastectomy == "No" | is.na(mastectomy)),
                                 (oophorectomy == "No" | is.na(oophorectomy)),
                                 gender == gender)$risk[1:age.max]
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
gen.fam <- function(n.sim, CP, af, seed = NULL, age.max = 94, age.min = 2,
                    censoring = TRUE){
  if(!is.null(seed)) set.seed(seed)
  fam.sim <- list()
  for(i in 1:n.sim){
    nSibsPatern <- sample(0:3, 2, replace = TRUE)
    nSibsMatern <- sample(0:3, 2, replace = TRUE)
    nSibs <- sample(0:3, 2, replace = TRUE)
    nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs) + 1), replace = TRUE), sum(nSibs) + 1, 2)
    fam.sim[[i]] <- tryCatch(sim.simFam(nSibsPatern, nSibsMatern, nSibs,
                                         nGrandchild, af, CP, includeGeno = TRUE, age.max = age.max,
                                        age.min = age.min, censoring = censoring),
                             error = function(e) NULL)
  }
  return(fam.sim)
}


pen.scale <- function(pen.obj, can, lt.risk.nc, lt.risk.c){
  ## This scales the penetrance so that the new lifetime risk for non-carriers is lt.risk.nc,
  ## and the new lifetime risk for carriers if lt.risk.c
  lt.scale <- function(pen, lt.risk){
    pen %*% diag(lt.risk / colSums(pen))
  }
  for(i in 1:length(can)){
    if(sum(pen.obj[can[i]][[1]][, 1]) != 0){ 
      pen.obj[can[i]][[1]][, 1] <- pen.obj[can[i]][[1]][, 1] * lt.risk.nc / sum(pen.obj[can[i]][[1]][, 1])
    }
  }
  for(i in 1:length(can)){
    if(all(colSums(pen.obj[can[i]][[1]][, -1]) != 0)){ 
      pen.obj[can[i]][[1]][, -1] <- lt.scale(pen.obj[can[i]][[1]][, -1], lt.risk.c)
    }
  }
  return(pen.obj)
}


