### Changing the penetrance so that homozygous carriers are the same as
## heterozygous carriers
## Last updated: September 7, 2018

hetpen <- function(CP, gastric = TRUE, scl = 1, pwr = NULL){
  
  # net2crude <- function(pen.c.net, pen.d.crude){
  #   pen.d.net <- c(pen.d.crude[1], pen.d.crude[-1] / (1 - cumsum(pen.c.net)[-1]))
  #   return(pen.c.net * (1 - cumsum(pen.d.net)))
  # }
  # 
  # mult.genes <- function(pen){
  #   cdf <- 1 - apply(apply(pen, 2, function(x) 1 - cumsum(x)), 1, prod)
  #   return(c(cdf[1], diff(cdf)))
  # }
  # 
  homoz.genes <- function(pen, pwr){
    cdf <- 1 - (1 - cumsum(pen))^pwr
    return(c(cdf[1], diff(cdf)))
  }
  
  homo2hetero <- function(pen, geno){
    for(i in 2:4){
      substr(geno, i, i) <- ifelse(substr(geno, i, i) == "2", "1", substr(geno, i, i))
    }
    return(pen[, which(colnames(pen) == geno)])
  }
  
  
  
  penet.mmr.net.g <- penet.mmr.net
  # penet.mmr.crude.g <- penet.mmr.crude
  
  pen.crc.f <- CP$cancerFDens[1:94, , "ColorC"]
  pen.crc.m <- CP$cancerMDens[1:94, , "ColorC"]
  pen.ec.f <- CP$cancerFDens[1:94, , "EndomC"]
  pen.ec.m <- CP$cancerMDens[1:94, , "EndomC"]
  if(gastric == TRUE){
    pen.gc.f <- CP$cancerFDens[1:94, , "GastC"]
    pen.gc.m <- CP$cancerMDens[1:94, , "GastC"]
    
    ## scaling the gastric cancer penetrance or raising the survival function to a power
    pen.gc.m <- pen.gc.m * scl
    pen.gc.f <- pen.gc.f * scl
    if(!is.null(pwr)){
      pen.gc.m <- apply(pen.gc.m, 2, homoz.genes, pwr = pwr)
      pen.gc.f <- apply(pen.gc.f, 2, homoz.genes, pwr = pwr)
    }
  }
  

  

  
  ## adding gastric cancer penetrance to MMRpro penetrance
  
  # female CRC
  penet.mmr.net.g$fFX[, "M000"] <- pen.crc.f[, "none"]
  penet.mmr.net.g$fFX[, "M100"] <- pen.crc.f[, "MLH1"]
  penet.mmr.net.g$fFX[, "M010"] <- pen.crc.f[, "MSH2"]
  penet.mmr.net.g$fFX[, "M001"] <- pen.crc.f[, "MSH6"]
  penet.mmr.net.g$fFX[, "M110"] <- pen.crc.f[, "MLH1.MSH2"]
  penet.mmr.net.g$fFX[, "M101"] <- pen.crc.f[, "MLH1.MSH6"]
  penet.mmr.net.g$fFX[, "M011"] <- pen.crc.f[, "MSH2.MSH6"]
  penet.mmr.net.g$fFX[, "M111"] <- pen.crc.f[, "MLH1.MSH2.MSH6"]
  for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
    penet.mmr.net.g$fFX[, i] <- homo2hetero(penet.mmr.net.g$fFX, colnames(penet.mmr.net.g$fFX)[i])
  }
  
  # male CRC
  penet.mmr.net.g$fMX[, "M000"] <- pen.crc.m[, "none"]
  penet.mmr.net.g$fMX[, "M100"] <- pen.crc.m[, "MLH1"]
  penet.mmr.net.g$fMX[, "M010"] <- pen.crc.m[, "MSH2"]
  penet.mmr.net.g$fMX[, "M001"] <- pen.crc.m[, "MSH6"]
  penet.mmr.net.g$fMX[, "M110"] <- pen.crc.m[, "MLH1.MSH2"]
  penet.mmr.net.g$fMX[, "M101"] <- pen.crc.m[, "MLH1.MSH6"]
  penet.mmr.net.g$fMX[, "M011"] <- pen.crc.m[, "MSH2.MSH6"]
  penet.mmr.net.g$fMX[, "M111"] <- pen.crc.m[, "MLH1.MSH2.MSH6"]
  for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
    penet.mmr.net.g$fMX[, i] <- homo2hetero(penet.mmr.net.g$fMX, colnames(penet.mmr.net.g$fMX)[i])
  }
  
  # female EC
  penet.mmr.net.g$fFY[, "M000"] <- pen.ec.f[, "none"]
  penet.mmr.net.g$fFY[, "M100"] <- pen.ec.f[, "MLH1"]
  penet.mmr.net.g$fFY[, "M010"] <- pen.ec.f[, "MSH2"]
  penet.mmr.net.g$fFY[, "M001"] <- pen.ec.f[, "MSH6"]
  penet.mmr.net.g$fFY[, "M110"] <- pen.ec.f[, "MLH1.MSH2"]
  penet.mmr.net.g$fFY[, "M101"] <- pen.ec.f[, "MLH1.MSH6"]
  penet.mmr.net.g$fFY[, "M011"] <- pen.ec.f[, "MSH2.MSH6"]
  penet.mmr.net.g$fFY[, "M111"] <- pen.ec.f[, "MLH1.MSH2.MSH6"]
  for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
    penet.mmr.net.g$fFY[, i] <- homo2hetero(penet.mmr.net.g$fFY, colnames(penet.mmr.net.g$fFY)[i])
  }
  
  # male EC
  penet.mmr.net.g$fMY[, "M000"] <- pen.ec.m[, "none"]
  penet.mmr.net.g$fMY[, "M100"] <- pen.ec.m[, "MLH1"]
  penet.mmr.net.g$fMY[, "M010"] <- pen.ec.m[, "MSH2"]
  penet.mmr.net.g$fMY[, "M001"] <- pen.ec.m[, "MSH6"]
  penet.mmr.net.g$fMY[, "M110"] <- pen.ec.m[, "MLH1.MSH2"]
  penet.mmr.net.g$fMY[, "M101"] <- pen.ec.m[, "MLH1.MSH6"]
  penet.mmr.net.g$fMY[, "M011"] <- pen.ec.m[, "MSH2.MSH6"]
  penet.mmr.net.g$fMY[, "M111"] <- pen.ec.m[, "MLH1.MSH2.MSH6"]
  for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
    penet.mmr.net.g$fMY[, i] <- homo2hetero(penet.mmr.net.g$fMY, colnames(penet.mmr.net.g$fMY)[i])
  }
  
  if(gastric == TRUE){
    # female GC
    penet.mmr.net.g$fFZ <- penet.mmr.net$fFX
    penet.mmr.net.g$fFZ[, "M000"] <- pen.gc.f[, "none"]
    penet.mmr.net.g$fFZ[, "M100"] <- pen.gc.f[, "MLH1"]
    penet.mmr.net.g$fFZ[, "M010"] <- pen.gc.f[, "MSH2"]
    penet.mmr.net.g$fFZ[, "M001"] <- pen.gc.f[, "MSH6"]
    penet.mmr.net.g$fFZ[, "M110"] <- pen.gc.f[, "MLH1.MSH2"]
    penet.mmr.net.g$fFZ[, "M101"] <- pen.gc.f[, "MLH1.MSH6"]
    penet.mmr.net.g$fFZ[, "M011"] <- pen.gc.f[, "MSH2.MSH6"]
    penet.mmr.net.g$fFZ[, "M111"] <- pen.gc.f[, "MLH1.MSH2.MSH6"]
    for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
      penet.mmr.net.g$fFZ[, i] <- homo2hetero(penet.mmr.net.g$fFZ, colnames(penet.mmr.net.g$fFZ)[i])
    }
    
    # male GC
    penet.mmr.net.g$fMZ <- penet.mmr.net$fMX
    penet.mmr.net.g$fMZ[, "M000"] <- pen.gc.m[, "none"]
    penet.mmr.net.g$fMZ[, "M100"] <- pen.gc.m[, "MLH1"]
    penet.mmr.net.g$fMZ[, "M010"] <- pen.gc.m[, "MSH2"]
    penet.mmr.net.g$fMZ[, "M001"] <- pen.gc.m[, "MSH6"]
    penet.mmr.net.g$fMZ[, "M110"] <- pen.gc.m[, "MLH1.MSH2"]
    penet.mmr.net.g$fMZ[, "M101"] <- pen.gc.m[, "MLH1.MSH6"]
    penet.mmr.net.g$fMZ[, "M011"] <- pen.gc.m[, "MSH2.MSH6"]
    penet.mmr.net.g$fMZ[, "M111"] <- pen.gc.m[, "MLH1.MSH2.MSH6"]
    for(i in (1:27)[-c(1, 2, 4, 5, 10, 11, 13, 14)]){
      penet.mmr.net.g$fMZ[, i] <- homo2hetero(penet.mmr.net.g$fMZ, colnames(penet.mmr.net.g$fMZ)[i])
    }
  }
  
  # # crude
  # penet.mmr.crude.g$fFZ <- penet.mmr.crude$fFX
  # penet.mmr.crude.g$fMZ <- penet.mmr.crude$fFX
  # for(i in 1:ncol(penet.mmr.crude.g$fFZ)){
  #   penet.mmr.crude.g$fFZ[, i] <- net2crude(penet.mmr.net.g$fFZ[, i],
  #                                           pen.d.gcf.crude)
  #   penet.mmr.crude.g$fMZ[, i] <- net2crude(penet.mmr.net.g$fMZ[, i],
  #                                           pen.d.gcm.crude)
  # }
  
  return(penet.mmr.net.g)
}


