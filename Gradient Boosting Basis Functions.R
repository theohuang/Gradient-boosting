### Functions to check if a family is "vertical" or "horoizontal"
### Last updated: January 30, 2018


## Checking if a family has 3 relatives in a row with BC
## Outputting the number of "vertical" instances (a relative can be a 
## member of multliple instances, for example if two siblings have BC
## and their mother and maternal grandmother also have BC)
vert <- function(fam){
  num.vert <- 0
  id.bc <- fam$ID[fam$AffectedBreast == 1]
  if(length(id.bc) < 3){
    return(num.vert)
  }
  for(i in 1:length(id.bc)){
    id.f <- fam$FatherID[fam$ID == id.bc[i]]
    if(id.f %in% fam$ID){
      if(fam$AffectedBreast[fam$ID == id.f] == 1){
        id.pgf <- fam$FatherID[fam$ID == id.f]
        if(id.pgf %in% fam$ID){
          if(fam$AffectedBreast[fam$ID == id.pgf] == 1) num.vert <- num.vert + 1
        }
        id.pgm <- fam$MotherID[fam$ID == id.f]
        if(id.pgm %in% fam$ID){
          if(fam$AffectedBreast[fam$ID == id.pgm] == 1) num.vert <- num.vert + 1
        }
      }
    }
    id.m <- fam$MotherID[fam$ID == id.bc[i]]
    if(id.m %in% fam$ID){
      if(fam$AffectedBreast[fam$ID == id.m] == 1){
        id.mgf <- fam$FatherID[fam$ID == id.m]
        if(id.mgf %in% fam$ID){
          if(fam$AffectedBreast[fam$ID == id.mgf] == 1) num.vert <- num.vert + 1
        }
        id.mgm <- fam$MotherID[fam$ID == id.m]
        if(id.mgm %in% fam$ID){
          if(fam$AffectedBreast[fam$ID == id.mgm] == 1) num.vert <- num.vert + 1
        }
      }
    }
  }
  return(num.vert)
}

## Checking if a family has 2 relatives in the same generation with BC
## Outputting the number of "horizontal" instances
horz <- function(fam){
  num.horz <- 0
  id.used <- vector()
  id.bc <- fam$ID[fam$AffectedBreast == 1]
  if(length(id.bc) < 2){
    return(num.horz)
  }
  for(i in 1:length(id.bc)){
    if(id.bc[i] %in% id.used) next ## avoiding double-counting relatives
    id.rel <- vector()
    id.f <- fam$FatherID[fam$ID == id.bc[i]]
    if(id.f %in% fam$ID){
      id.rel <- fam$ID[fam$FatherID == id.f]
    }
    id.m <- fam$MotherID[fam$ID == id.bc[i]]
    if(id.m %in% fam$ID){
      id.rel <- unique(c(id.rel, fam$ID[fam$MotherID == id.m]))
    }
    id.used <- c(id.used, id.rel)
    num.rel.bc <- sum(fam$AffectedBreast[fam$ID %in% id.rel])
    if(num.rel.bc >= 2){
      num.horz <- num.horz + 1
    }
  }
  return(num.horz)
}

breast <- function(fam){
  br <- as.numeric(fam$AffectedBreast[which(fam$Relation == 1)])
  if(is.na(br)){
    return(0)
  } else{
    return(br)
  }
}

ovary <- function(fam){
  ov <- as.numeric(fam$AffectedOvary[which(fam$Relation == 1)])
  if(is.na(ov)){
    return(0)
  } else{
    return(ov)
  }
}

prostate <- function(fam){
  pros <- as.numeric(fam$AffectedProstate[which(fam$Relation == 1)])
  if(is.na(pros)){
    return(0)
  } else{
    return(pros)
  }
}

endometrium <- function(fam){
  endo <- as.numeric(fam$AffectedEndometrium[which(fam$Relation == 1)])
  if(is.na(endo)){
    return(0)
  } else{
    return(endo)
  }
}

colon <- function(fam){
  colo <- as.numeric(fam$AffectedColon[which(fam$Relation == 1)])
  if(is.na(colo)){
    return(0)
  } else{
    return(colo)
  }
}

## Number of horiztonal siblings divided by number of siblings
## If there are more than one horizontal case, take the average
horz.perc <- function(fam){
  per.horz <- 0
  id.used <- vector()
  id.bc <- fam$ID[fam$AffectedBreast == 1]
  if(length(id.bc) < 2){
    return(per.horz)
  }
  for(i in 1:length(id.bc)){
    if(id.bc[i] %in% id.used) next ## avoiding double-counting relatives
    id.rel <- vector()
    id.f <- fam$FatherID[fam$ID == id.bc[i]]
    if(id.f %in% fam$ID){
      id.rel <- fam$ID[fam$FatherID == id.f]
    }
    id.m <- fam$MotherID[fam$ID == id.bc[i]]
    if(id.m %in% fam$ID){
      id.rel <- unique(c(id.rel, fam$ID[fam$MotherID == id.m]))
    }
    id.used <- c(id.used, id.rel)
    num.rel.bc <- sum(fam$AffectedBreast[fam$ID %in% id.rel])
    if(num.rel.bc >= 2){
      per.horz <- c(per.horz, num.rel.bc / length(id.rel))
    }
  }
  return(mean(per.horz))
}

## Number of horiztonal siblings divided by number of siblings
## If there are more than one horizontal case, take the average
horz.perc2 <- function(fam){
  ind.pro <- which(fam$Relative == 1)
  return(mean(fam$AffectedBreast[fam$MotherID == fam$MotherID[ind.pro] |
                            fam$FatherID == fam$FatherID[ind.pro]]))
}

## number of immediate family members with breast cancer
bc.imm <- function(fam){
  id.bc <- fam$ID[fam$AffectedBreast == 1]
  num.imm <- length(which(fam$Relation[fam$ID %in% id.bc] %in% 2:4))
  return(num.imm)
}

  
num.gen <- function(fam){
  up.1 <- as.numeric(length(which(fam$Relation %in% c(4, 6, 8))) >= 1)
  up.2 <- as.numeric(length(which(fam$Relation %in% c(5, 7))) >= 1)
  down.1 <- as.numeric(length(which(fam$Relation %in% c(3, 13))) >= 1)
  return(up.1 + up.2 + down.1 + 1)
}  
  
