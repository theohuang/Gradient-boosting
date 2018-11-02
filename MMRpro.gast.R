# Last update:  March 2008
# Use peeling algorithm, v2.0

## version 2.1-5


## MMRPRO that takes gastric cancer information
## Last updated: August 23, 2018

MMRpro.gast <- function(family, counselee.id=1, race = "Unknown", germline.testing=NULL, 
                        marker.testing=NULL, params=MMRparams(), print=TRUE, 
                        imputeAges=TRUE, imputeRelatives=TRUE, net = TRUE) {
  
  age.max <- 94
  
  pro.agemax <- 0
  ## if proband is older than the maximum age - 1, then we only return carrier probabilities
  if(max(family[which(family$ID == counselee.id), c("AgeColon", "AgeEndometrium", "AgeGastric")], na.rm = TRUE) > age.max - 1){
    pro.agemax <- 1
  }
  
  # probability of inherited genetic susceptibility
  #
  # family: fam hist in MMRpro format. If 12 column, add 0 to msh6 and msi
  #         if 13 column, the last col is assumed to be msh6.
  #         To enter MMRpro files use readfamfilter.
  
  # output: joint probabiilty of genotypes
  
  warnmsg = ""
  if (imputeRelatives==TRUE & imputeAges==FALSE) {
    imputeAges=TRUE
    warnmsg = paste(warnmsg, "Warning: imputeAges was input as false, but has been set to TRUE.  If imputeRelatives=TRUE then by default imputeAges must also be TRUE.")
  }
  if (nchar(warnmsg) > 0) {print(warnmsg)}
  
  proband.current.age <- max(family[which(family$ID == counselee.id), c("AgeColon", "AgeEndometrium")])
  
  
  # Set up allele frequencies to work with new peelingR
  allef <- list()
  allef["Lynch"] <- list(params$allef)
  family$ethnic <- rep("Lynch", nrow(family))
  
  
  family <- CheckFamStructure.gast("MMRpro",fff=family, counselee.id=counselee.id, germline.testing=germline.testing, marker.testing=marker.testing, 
                                   params=params, imputeAges=imputeAges, imputeRelatives=imputeRelatives)
  if(is.character(family)) {return(family)}
  
  ## truncating ages past age.max to age.max, since we don't have data past age age.max
  if("AgeColon" %in% colnames(family)){
    family$AgeColon[family$AgeColon > age.max] <- age.max
  }
  if("AgeEndometrium" %in% colnames(family)){
    family$AgeEndometrium[family$AgeEndometrium > age.max] <- age.max
  }
  if("AgeGastric" %in% colnames(family)){
    family$AgeGastric[family$AgeGastric > age.max] <- age.max
  }
  
  # if the counselee has an identical twin, make sure the ID used in the calculation is that of the counselee and not his/her twin so peeling does not error out
  originalcounselee.id <- counselee.id
  if (family$Twins[family$ID==counselee.id]>0) {
    twins <- family[family$Twins==family$Twins[family$ID==counselee.id],]
    counselee.id <- twins$ID[1]
  }
  
  # Get the race-specific baseline estimates #
  
  
  if (race=="Unknown") {
    params$penetrance.net <- params$penetrance.net
    params$penetrance.crude <- params$penetrance.crude
    haz.d <- haz.death.crude
  }else{
    if (race=="Asian") {
      params$penetrance.net$fFX[,"M000"] <- baseline.race.mmr.net.gc$fFX[,"Asian"]
      params$penetrance.net$fFY[,"M000"] <- baseline.race.mmr.net.gc$fFY[,"Asian"]
      params$penetrance.net$fFZ[,"M000"] <- baseline.race.mmr.net.gc$fFZ[,"Asian"]
      params$penetrance.net$fMX[,"M000"] <- baseline.race.mmr.net.gc$fMX[,"Asian"]
      params$penetrance.net$fMY[,"M000"] <- baseline.race.mmr.net.gc$fMY[,"Asian"]
      params$penetrance.net$fMZ[,"M000"] <- baseline.race.mmr.net.gc$fMZ[,"Asian"]
      
      params$penetrance.crude$fFX[,"M000"] <- baseline.race.mmr.crude.gc$fFX[,"Asian"]
      params$penetrance.crude$fFY[,"M000"] <- baseline.race.mmr.crude.gc$fFY[,"Asian"]
      params$penetrance.crude$fFZ[,"M000"] <- baseline.race.mmr.crude.gc$fFZ[,"Asian"]
      params$penetrance.crude$fMX[,"M000"] <- baseline.race.mmr.crude.gc$fMX[,"Asian"]
      params$penetrance.crude$fMY[,"M000"] <- baseline.race.mmr.crude.gc$fMY[,"Asian"]
      params$penetrance.crude$fMZ[,"M000"] <- baseline.race.mmr.crude.gc$fMZ[,"Asian"]
      
      haz.d <- haz.death.race.crude.gc$Asian
    }else{
      if (race=="Black") {
        params$penetrance.net$fFX[,"M000"] <- baseline.race.mmr.net.gc$fFX[,"Black"]
        params$penetrance.net$fFY[,"M000"] <- baseline.race.mmr.net.gc$fFY[,"Black"]
        params$penetrance.net$fFZ[,"M000"] <- baseline.race.mmr.net.gc$fFZ[,"Black"]
        params$penetrance.net$fMX[,"M000"] <- baseline.race.mmr.net.gc$fMX[,"Black"]
        params$penetrance.net$fMY[,"M000"] <- baseline.race.mmr.net.gc$fMY[,"Black"]
        params$penetrance.net$fMZ[,"M000"] <- baseline.race.mmr.net.gc$fMZ[,"Black"]
        
        params$penetrance.crude$fFX[,"M000"] <- baseline.race.mmr.crude.gc$fFX[,"Black"]
        params$penetrance.crude$fFY[,"M000"] <- baseline.race.mmr.crude.gc$fFY[,"Black"]
        params$penetrance.crude$fFZ[,"M000"] <- baseline.race.mmr.crude.gc$fFZ[,"Black"]
        params$penetrance.crude$fMX[,"M000"] <- baseline.race.mmr.crude.gc$fMX[,"Black"]
        params$penetrance.crude$fMY[,"M000"] <- baseline.race.mmr.crude.gc$fMY[,"Black"]
        params$penetrance.crude$fMZ[,"M000"] <- baseline.race.mmr.crude.gc$fMZ[,"Black"]
        
        haz.d <- haz.death.race.crude.gc$Black
      }else{
        if (race=="Hispanic") {
          params$penetrance.net$fFX[,"M000"] <- baseline.race.mmr.net.gc$fFX[,"Hispanic"]
          params$penetrance.net$fFY[,"M000"] <- baseline.race.mmr.net.gc$fFY[,"Hispanic"]
          params$penetrance.net$fFZ[,"M000"] <- baseline.race.mmr.net.gc$fFZ[,"Hispanic"]
          params$penetrance.net$fMX[,"M000"] <- baseline.race.mmr.net.gc$fMX[,"Hispanic"]
          params$penetrance.net$fMY[,"M000"] <- baseline.race.mmr.net.gc$fMY[,"Hispanic"]
          params$penetrance.net$fMZ[,"M000"] <- baseline.race.mmr.net.gc$fMZ[,"Hispanic"]
          
          params$penetrance.crude$fFX[,"M000"] <- baseline.race.mmr.crude.gc$fFX[,"Hispanic"]
          params$penetrance.crude$fFY[,"M000"] <- baseline.race.mmr.crude.gc$fFY[,"Hispanic"]
          params$penetrance.crude$fFZ[,"M000"] <- baseline.race.mmr.crude.gc$fFZ[,"Hispanic"]
          params$penetrance.crude$fMX[,"M000"] <- baseline.race.mmr.crude.gc$fMX[,"Hispanic"]
          params$penetrance.crude$fMY[,"M000"] <- baseline.race.mmr.crude.gc$fMY[,"Hispanic"]
          params$penetrance.crude$fMZ[,"M000"] <- baseline.race.mmr.crude.gc$fMZ[,"Hispanic"]
          
          haz.d <- haz.death.race.crude.gc$Hispanic
        }else{
          if (race=="NativeAmerican") {
            params$penetrance.net$fFX[,"M000"] <- baseline.race.mmr.net.gc$fFX[,"NativeAmerican"]
            params$penetrance.net$fFY[,"M000"] <- baseline.race.mmr.net.gc$fFY[,"NativeAmerican"]
            params$penetrance.net$fFZ[,"M000"] <- baseline.race.mmr.net.gc$fFZ[,"NativeAmerican"]
            params$penetrance.net$fMX[,"M000"] <- baseline.race.mmr.net.gc$fMX[,"NativeAmerican"]
            params$penetrance.net$fMY[,"M000"] <- baseline.race.mmr.net.gc$fMY[,"NativeAmerican"]
            params$penetrance.net$fMZ[,"M000"] <- baseline.race.mmr.net.gc$fMZ[,"NativeAmerican"]
            
            params$penetrance.crude$fFX[,"M000"] <- baseline.race.mmr.crude.gc$fFX[,"NativeAmerican"]
            params$penetrance.crude$fFY[,"M000"] <- baseline.race.mmr.crude.gc$fFY[,"NativeAmerican"]
            params$penetrance.crude$fFZ[,"M000"] <- baseline.race.mmr.crude.gc$fFZ[,"NativeAmerican"]
            params$penetrance.crude$fMX[,"M000"] <- baseline.race.mmr.crude.gc$fMX[,"NativeAmerican"]
            params$penetrance.crude$fMY[,"M000"] <- baseline.race.mmr.crude.gc$fMY[,"NativeAmerican"]
            params$penetrance.crude$fMZ[,"M000"] <- baseline.race.mmr.crude.gc$fMZ[,"NativeAmerican"]
            
            haz.d <- haz.death.race.crude.gc$NativeAmerican
          }else{
            if (race=="White") {
              params$penetrance.net$fFX[,"M000"] <- baseline.race.mmr.net.gc$fFX[,"White"]
              params$penetrance.net$fFY[,"M000"] <- baseline.race.mmr.net.gc$fFY[,"White"]
              params$penetrance.net$fFZ[,"M000"] <- baseline.race.mmr.net.gc$fFZ[,"White"]
              params$penetrance.net$fMX[,"M000"] <- baseline.race.mmr.net.gc$fMX[,"White"]
              params$penetrance.net$fMY[,"M000"] <- baseline.race.mmr.net.gc$fMY[,"White"]
              params$penetrance.net$fMZ[,"M000"] <- baseline.race.mmr.net.gc$fMZ[,"White"]
              
              params$penetrance.crude$fFX[,"M000"] <- baseline.race.mmr.crude.gc$fFX[,"White"]
              params$penetrance.crude$fFY[,"M000"] <- baseline.race.mmr.crude.gc$fFY[,"White"]
              params$penetrance.crude$fFZ[,"M000"] <- baseline.race.mmr.crude.gc$fFZ[,"White"]
              params$penetrance.crude$fMX[,"M000"] <- baseline.race.mmr.crude.gc$fMX[,"White"]
              params$penetrance.crude$fMY[,"M000"] <- baseline.race.mmr.crude.gc$fMY[,"White"]
              params$penetrance.crude$fMZ[,"M000"] <- baseline.race.mmr.crude.gc$fMZ[,"White"]
              
              haz.d <- haz.death.race.crude.gc$White
            }
          }}}}}
  
  
  psize <- dim(family)[1] #size of the family
  ndis <- 3 # colon, endo, gastric
  nloci <- 3
  ngen <- 27 # 3^nloci
  horizon <- length(params$penetrance.net$fFX[,1])
  
  
  colon <- 1
  endo <- 2
  gastric <- 3
  DIS <- matrix(0,psize,ndis)   # disease status
  DIS[,colon][family[,"AffectedColon"]>0] <- 1    # colon
  DIS[,endo][family[,"AffectedEndometrium"]>0] <- 1   # endo
  DIS[,gastric][family[,"AffectedGastric"]>0] <- 1   # gastric
  
  
  ages.cancer <- data.frame(family[,"AgeColon",drop=FALSE],
                            family[,"AgeEndometrium",drop=FALSE],
                            family[,"AgeGastric",drop=FALSE])
  maxage <- apply(ages.cancer,1,max)
  # AGE <- matrix(0,psize,ndis)   
  # AGE[,colon] <- family[,"AgeColon"]   
  # AGE[,endo] <- family[,"AgeEndometrium"]
  # AGE[,gastric] <- family[,"AgeGastric"]
  
  
  males <- (1:psize)[family[,"Gender"]==1]
  females <- (1:psize)[family[,"Gender"]==0]
  
  PPf <- FFf <- array(0,c(horizon,ndis,ngen))  # penetrances & cumu, females
  PPm <- FFm <- array(0,c(horizon,ndis,ngen))  # penetrances & cumu, males
  
  PPf[,colon,] <- as.matrix(params$penetrance.net$fFX)
  PPf[,endo,] <- as.matrix(params$penetrance.net$fFY)
  PPf[,gastric,] <- as.matrix(params$penetrance.net$fFZ)
  
  PPm[,colon,] <- as.matrix(params$penetrance.net$fMX)
  PPm[,endo,] <- as.matrix(params$penetrance.net$fMY)
  PPm[,gastric,] <- as.matrix(params$penetrance.net$fMZ)
  
  
  FFf[,colon,] <- apply(PPf[,colon,],2,cumsum)
  FFf[,endo, ] <- apply(PPf[,endo, ],2,cumsum)
  FFf[,gastric, ] <- apply(PPf[,gastric, ],2,cumsum)
  FFm[,colon,] <- apply(PPm[,colon,],2,cumsum)
  FFm[,endo, ] <- apply(PPm[,endo, ],2,cumsum)
  FFm[,gastric, ] <- apply(PPm[,gastric, ],2,cumsum)
  
  
  femalescrc <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==1]
  malescrc <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==1]
  femalesnocrc <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==0]
  malesnocrc <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==0]
  
  ## Make a within-*pro function that computes needed components
  ## for capturing the LIK and ped objects required by peeling
  
  runPeelingMMRpro = function(family, params, ngen,
                              germline.testing, 
                              marker.testing, psize, counselee.id, allef)
  {
    
    ## first getting the ages
    AGE <- matrix(0,psize,ndis)   
    AGE[,colon] <- family[,"AgeColon"]   
    AGE[,endo] <- family[,"AgeEndometrium"]
    AGE[,gastric] <- family[,"AgeGastric"]
    
    
    ########################################################
    # No tumor location info #
    
    
    if (is.null(marker.testing$location)) {
      
      AFF <- array(0,c(psize,ndis,ngen)) # likelihood if affected
      for ( gg in 1:ngen ) { for ( dd in 1:ndis ) {
        AFF[females,dd,gg] <- PPf[AGE[females,dd],dd,gg]
        AFF[males,dd,gg] <- PPm[AGE[males,dd],dd,gg]
      }}
      
      UNA <- array(0,c(psize,ndis,ngen)) # likelihood if UNAffeected
      for ( gg in 1:ngen ) { for ( dd in 1:ndis ){
        UNA[females,dd,gg] <- (1-FFf[AGE[females,dd],dd,gg])
        UNA[males,dd,gg] <- (1-FFm[AGE[males,dd],dd,gg])
      }}
      
    }
    #########################################################
    # Tumor location info #
    
    if (!is.null(marker.testing$location)) {
      
      femalescrc.prox <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==1 & family[,"location"]==1]
      femalescrc.dist <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==1 & family[,"location"]==2]
      femalescrc.unk <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==1 & family[,"location"]==0]
      femalesunaff <- (1:psize)[family[,"Gender"]==0 & DIS[,colon]==0]
      
      malescrc.prox <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==1 & family[,"location"]==1]
      malescrc.dist <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==1 & family[,"location"]==2]
      malescrc.unk <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==1 & family[,"location"]==0]
      malesunaff <- (1:psize)[family[,"Gender"]==1 & DIS[,colon]==0]
      
      p.prox.carriers <- 0.873
      p.dist.carriers <- 0.126
      p.prox.noncarriers <- 0.375
      p.dist.noncarriers <- 0.625
      
      AFF <- array(0,c(psize,ndis,ngen)) # likelihood if affected
      
      for (dd in 1:ndis) { # genotype profile of [0,0,0]: non-carriers
        AFF[femalescrc.prox,dd,1] <- PPf[AGE[femalescrc.prox,dd],dd,1]*p.prox.noncarriers 
        AFF[malescrc.prox,dd,1] <- PPm[AGE[malescrc.prox,dd],dd,1]*p.prox.noncarriers 
        AFF[femalescrc.dist,dd,1] <- PPf[AGE[femalescrc.dist,dd],dd,1]*p.dist.noncarriers 
        AFF[malescrc.dist,dd,1] <- PPm[AGE[malescrc.dist,dd],dd,1]*p.dist.noncarriers 
        AFF[femalescrc.unk,dd,1] <- PPf[AGE[femalescrc.unk,dd],dd,1] 
        AFF[malescrc.unk,dd,1] <- PPm[AGE[malescrc.unk,dd],dd,1] 
        AFF[femalesunaff,dd,1] <- PPf[AGE[femalesunaff,dd],dd,1] 
        AFF[malesunaff,dd,1] <- PPm[AGE[malesunaff,dd],dd,1] 
      }
      for ( gg in 2:ngen ) { for ( dd in 1:ndis ) {
        AFF[femalescrc.prox,dd,gg] <- PPf[AGE[femalescrc.prox,dd],dd,gg]*p.prox.carriers 
        AFF[malescrc.prox,dd,gg] <- PPm[AGE[malescrc.prox,dd],dd,gg]*p.prox.carriers 
        AFF[femalescrc.dist,dd,gg] <- PPf[AGE[femalescrc.dist,dd],dd,gg]*p.dist.carriers 
        AFF[malescrc.dist,dd,gg] <- PPm[AGE[malescrc.dist,dd],dd,gg]*p.dist.carriers 
        AFF[femalescrc.unk,dd,gg] <- PPf[AGE[femalescrc.unk,dd],dd,gg] 
        AFF[malescrc.unk,dd,gg] <- PPm[AGE[malescrc.unk,dd],dd,gg] 
        AFF[femalesunaff,dd,gg] <- PPf[AGE[femalesunaff,dd],dd,gg] 
        AFF[malesunaff,dd,gg] <- PPm[AGE[malesunaff,dd],dd,gg] 
      }}
      
      UNA <- array(0,c(psize,ndis,ngen)) # likelihood if UNAffeected
      for ( gg in 1:ngen ) { for ( dd in 1:ndis ){
        UNA[females,dd,gg] <- (1-FFf[AGE[females,dd],dd,gg])
        UNA[males,dd,gg] <- (1-FFm[AGE[males,dd],dd,gg])
      }}
      
    }
    
    ###########################################################	 
    # Contribution of germline testing results. Pr(Test|gene) #
    
    TES1 <- TES2 <- TES3 <- matrix(1,psize,ngen)
    if(!is.null(germline.testing) & (sum(family[,"MLH1"])>0 | sum(family[,"MSH2"])>0 |
                                     sum(family[,"MSH6"])>0)){
      
      # Create a vector indicating which columns of the penetrance are informative of that gene's carrier status
      mlh1cols <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27)
      msh2cols <- c(4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27)
      msh6cols <- 10:27
      
      if (sum(germline.testing$TestOrder)==0) { # Unknown test order
        whogottested <- which(family[,"MLH1"]>0 | family[,"MSH2"]>0 | family[,"MSH6"]>0)
        ntested <- length(whogottested)
        if (is.element(counselee.id, whogottested)) {
          germline.testing$TestOrder[which(family$ID==counselee.id)] <- ntested # make the proband as the last person who got tested
          germline.testing$TestOrder[setdiff(whogottested, counselee.id)] <- sample(1:(ntested-1),size=ntested-1)
        }
        if (!is.element(counselee.id, whogottested)) {
          germline.testing$TestOrder[whogottested] <- sample(1:ntested, size=ntested)
        }
        
        print("Warning: Testing order has been imputed (see help(MMRpro) for details).")
      }
      
      
      # Need to group members with test results as follows:
      # 'pre.pos.neg': those who test negative before anyone tests positive
      # 'first.pos': first person to test positive
      # 'after.pos' and 'after.neg': anyone tested after the first positive
      
      # MLH1 #
      positives.mlh1 <- which(family[,"MLH1"]==1)
      negatives.mlh1 <- which(family[,"MLH1"]==2)
      
      if (length(positives.mlh1)>0) {
        first.pos.mlh1 <- which(family[,"MLH1"]==1 & germline.testing$TestOrder==min(germline.testing$TestOrder[positives.mlh1]))
        pre.pos.neg.mlh1 <- which(family[,"MLH1"]==2 & germline.testing$TestOrder<germline.testing$TestOrder[first.pos.mlh1])
        after.pos.mlh1 <- setdiff(positives.mlh1, first.pos.mlh1)
        after.neg.mlh1 <- setdiff(negatives.mlh1, pre.pos.neg.mlh1)
      }
      if (length(positives.mlh1)==0) {
        pre.pos.neg.mlh1 <- negatives.mlh1
        first.pos.mlh1 <- after.pos.mlh1 <- after.neg.mlh1 <- NULL
      }
      
      # MSH2 #
      positives.msh2 <- which(family[,"MSH2"]==1)
      negatives.msh2 <- which(family[,"MSH2"]==2)
      
      if (length(positives.msh2)>0) {
        first.pos.msh2 <- which(family[,"MSH2"]==1 & germline.testing$TestOrder==min(germline.testing$TestOrder[positives.msh2]))
        pre.pos.neg.msh2 <- which(family[,"MSH2"]==2 & germline.testing$TestOrder<germline.testing$TestOrder[first.pos.msh2])
        after.pos.msh2 <- setdiff(positives.msh2, first.pos.msh2)
        after.neg.msh2 <- setdiff(negatives.msh2, pre.pos.neg.msh2)
      }
      if (length(positives.msh2)==0) {
        pre.pos.neg.msh2 <- negatives.msh2
        first.pos.msh2 <- after.pos.msh2 <- after.neg.msh2 <- NULL
      }
      
      # MSH6 #
      positives.msh6 <- which(family[,"MSH6"]==1)
      negatives.msh6 <- which(family[,"MSH6"]==2)
      
      if (length(positives.msh6)>0) {
        first.pos.msh6 <- which(family[,"MSH6"]==1 & germline.testing$TestOrder==min(germline.testing$TestOrder[positives.msh6]))
        pre.pos.neg.msh6 <- which(family[,"MSH6"]==2 & germline.testing$TestOrder<germline.testing$TestOrder[first.pos.msh6])
        after.pos.msh6 <- setdiff(positives.msh6, first.pos.msh6)
        after.neg.msh6 <- setdiff(negatives.msh6, pre.pos.neg.msh6)
      }
      if (length(positives.msh6)==0) {
        pre.pos.neg.msh6 <- negatives.msh6
        first.pos.msh6 <- after.pos.msh6 <- after.neg.msh6 <- NULL
      }
      
      mlh1test <- family[,"MLH1"]
      msh2test <- family[,"MSH2"]
      msh6test <- family[,"MSH6"]
      
      if (length(first.pos.mlh1)>0) TES1[first.pos.mlh1,mlh1cols] <- params$sensitivity1
      if (length(pre.pos.neg.mlh1)>0) TES1[pre.pos.neg.mlh1,mlh1cols] <- 1-params$sensitivity1  
      if (length(pre.pos.neg.mlh1)>0) TES1[pre.pos.neg.mlh1,-mlh1cols] <- params$specificity1
      if (length(first.pos.mlh1)>0) TES1[first.pos.mlh1,-mlh1cols] <- 1-params$specificity1
      
      if (length(after.pos.mlh1)>0) TES1[after.pos.mlh1,mlh1cols] <- 1 - 2*allef[[family$ethnic[after.pos.mlh1]]][[1]][2]
      if (length(after.neg.mlh1)>0) TES1[after.neg.mlh1,mlh1cols] <- 2*allef[[family$ethnic[after.neg.mlh1]]][[1]][2] 
      if (length(after.neg.mlh1)>0) TES1[after.neg.mlh1,-mlh1cols] <- params$specificity1
      if (length(after.pos.mlh1)>0) TES1[after.pos.mlh1,-mlh1cols] <- 1-params$specificity1      
      
      if (length(first.pos.msh2)>0) TES2[first.pos.msh2,msh2cols] <- params$sensitivity2
      if (length(pre.pos.neg.msh2)>0) TES2[pre.pos.neg.msh2,msh2cols] <- 1-params$sensitivity2  
      if (length(pre.pos.neg.msh2)>0) TES2[pre.pos.neg.msh2,-msh2cols] <- params$specificity2
      if (length(first.pos.msh2)>0) TES2[first.pos.msh2,-msh2cols] <- 1-params$specificity2
      
      if (length(after.pos.msh2)>0) TES2[after.pos.msh2,msh2cols] <- 1 - 2*allef[[family$ethnic[after.pos.msh2]]][[2]][2]
      if (length(after.neg.msh2)>0) TES2[after.neg.msh2,msh2cols] <- 2*allef[[family$ethnic[after.neg.msh2]]][[2]][2] 
      if (length(after.neg.msh2)>0) TES2[after.neg.msh2,-msh2cols] <- params$specificity2
      if (length(after.pos.msh2)>0) TES2[after.pos.msh2,-msh2cols] <- 1-params$specificity2      
      
      if (length(first.pos.msh6)>0) TES3[first.pos.msh6,msh6cols] <- params$sensitivity3
      if (length(pre.pos.neg.msh6)>0) TES3[pre.pos.neg.msh6,msh6cols] <- 1-params$sensitivity3  
      if (length(pre.pos.neg.msh6)>0) TES3[pre.pos.neg.msh6,-msh6cols] <- params$specificity3
      if (length(first.pos.msh6)>0) TES3[first.pos.msh6,-msh6cols] <- 1-params$specificity3
      
      if (length(after.pos.msh6)>0) TES3[after.pos.msh6,msh6cols] <- 1 - 2*allef[[family$ethnic[after.pos.msh6]]][[3]][2]
      if (length(after.neg.msh6)>0) TES3[after.neg.msh6,msh6cols] <- 2*allef[[family$ethnic[after.neg.msh6]]][[3]][2] 
      if (length(after.neg.msh6)>0) TES3[after.neg.msh6,-msh6cols] <- params$specificity3
      if (length(after.pos.msh6)>0) TES3[after.pos.msh6,-msh6cols] <- 1-params$specificity3      
      
    }
    
    TES <- TES1*TES2*TES3
    
    ######################################################	 
    # Contribution of MSI testing results. Pr(Test|gene) #
    
    msi1 <- msi2 <- msi6 <- matrix(1,psize,ngen)
    if(!is.null(marker.testing)){
      msi <- family[,"MSI"]
      
      # Create a vector indicating which columns of the penetrance are informative of that gene's carrier status
      mlh1cols <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27)
      msh2cols <- c(4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27)
      msh6cols <- 10:27
      
      #we assume that MSI status is informative of MLH1 and MSH2 carrier status, then independently of MSH6 carrier status
      msi1[msi==1,mlh1cols] <- params$msi1.sens
      msi1[msi==2,mlh1cols] <- 1-params$msi1.sens
      msi1[msi==1,-mlh1cols] <- 1-params$msi1.spec
      msi1[msi==2,-mlh1cols] <- params$msi1.spec
      
      msi2[msi==1,msh2cols] <- params$msi2.sens
      msi2[msi==2,msh2cols] <- 1-params$msi2.sens
      msi2[msi==1,-msh2cols] <- 1-params$msi2.spec
      msi2[msi==2,-msh2cols] <- params$msi2.spec
      
      msi6[msi==1,msh6cols] <- params$msi6.sens
      msi6[msi==2,msh6cols] <- 1-params$msi6.sens
      msi6[msi==1,-msh6cols] <- 1-params$msi6.spec
      msi6[msi==2,-msh6cols] <- params$msi6.spec
    }
    MSI <- msi1*msi2 * msi6
    
    #-----------------------------------------------#
    # Compute each person's likelihood contribution #
    
    LIK <- matrix(1,psize,ngen)
    for (gg in 1:ngen) {
      LIK[,gg] <- apply( DIS * AFF[,,gg] + (1-DIS) * UNA[,,gg], 1, prod ) *
        TES[,gg] * MSI[,gg]
    }
    
    # If there are twins in the family, collapse them here #
    if (sum(family$Twins)>0) {
      if (length(unique(family$Twins))>=1 & !any(family$Twins==0)) {loopID <- length(unique(family$Twins))}
      if (length(unique(family$Twins))>1 & any(family$Twins==0)) {loopID <- length(unique(family$Twins))-1}
      family$NewMotherIndex = family$NewFatherIndex = 0
      
      # Because we are collapsing twins, need to replace any offspring's mid and fid
      for (iiiii in 1:loopID) {
        temp <- family[family$Twins==unique(family$Twins[family$Twins>0])[iiiii],]
        if (temp$Gender[1]==0) {
          family[family$Twins==unique(family$Twins[family$Twins>0])[iiiii],"NewMotherIndex"] = temp$ID[1]}
        if (temp$Gender[1]==1) {
          family[family$Twins==unique(family$Twins[family$Twins>0])[iiiii],"NewFatherIndex"] = temp$ID[1]}
      }
      
      twins = family$ID[family$Twins>0]
      offspring = family$ID[is.element(family$MotherID, twins) | is.element(family$FatherID, twins)]
      if (length(offspring)>0) {
        for (ooo in 1:length(offspring)){
          motherIsTwin = family$Twins[family$ID==family$MotherID[family$ID==offspring[ooo]]]
          fatherIsTwin = family$Twins[family$ID==family$FatherID[family$ID==offspring[ooo]]]
          if (length(motherIsTwin)!=0) {
            if (motherIsTwin!=0) {
              family$MotherID[family$ID==offspring[ooo]] = family$NewMotherIndex[family$ID==family$MotherID[family$ID==offspring[ooo]]]}
          }
          if (length(fatherIsTwin)!=0) {
            if (fatherIsTwin!=0) {
              family$FatherID[family$ID==offspring[ooo]] = family$NewFatherIndex[family$ID==family$FatherID[family$ID==offspring[ooo]]]}
          }
        }
      }
      # get the first twin's ped info for peeling
      firsttwin = which(family$Twins>0)
      firsttwin = firsttwin[match(unique(family$Twins[family$Twins!=0]), family$Twins[family$Twins!=0])]
      nottwins = which(family$Twins==0)
      collapsedtwins = sort(c(firsttwin,nottwins)) 
      ped = family[collapsedtwins,c("ID", "Gender", "FatherID", "MotherID", "ethnic")]
      # set up new LIK
      newLIK = LIK
      twinIDs = unique(family$Twins[family$Twins!=0])
      for (ttt in 1:length(twinIDs)) {
        twinrows = which(family$Twins==twinIDs[ttt])
        newLIK[twinrows[1],] = LIK[twinrows[1],] * LIK[twinrows[2],]
        newLIK[twinrows[2],] = rep(NA,ngen)
      }
      # Now take out empty rows
      LIK = newLIK[apply(newLIK,1,function(x){return(!any(is.na(x)))}),]
    }
    
    if (sum(family$Twins)==0) {
      ped = family[,c("ID", "Gender", "FatherID", "MotherID", "ethnic")]
    }
    
    #-------------------------------------------------# 
    # Compute posterior probability of carrier status #
    
    
    # peelingR  (peeling in R)
    #  outputs: numerator of BayesThereom calculation
    #           weight used for likelihood
    #  logp <- peelingR(allef,LIK,ped=family[,c("ID", "Gender", "FatherID", "MotherID","ethnic")],counselee.id=counselee.id,nloci=nloci)
    #  post <- rescale(logp$num)
    #  ll <- log(sum(logp$num))+logp$weight
    #
    # peelingRC (peeling in C)
    #  outputs: posterior probability
    #
    # When moving to peelingRC, decided to leave ll slot but set to NULL
    # dyn.load('peelingC')
    post <- peelingRC(allef=allef,LIK,ped=ped,counselee.id=counselee.id,nloci=nloci)
    
    return(post)
    
  }
  
  ## --------------
  ## Calculation for LIK start here.  This is where we start the multiple age imputation process ---------------###
  
  missingage = any(family$AgeColon==1 | family$AgeEndometrium==1 |
                     family$AgeGastric == 1)
  
  if (missingage & imputeAges==FALSE) {
    agewarning = "Warning: Age imputation has been turned off, but there are 
    unknown ages of some unaffected family members.  
    You may want to get more information about family member ages and re-run 
    the calculation, or set imputeAges=TRUE."
    post = runPeelingMMRpro(family=family, 
                            params=params, ngen=ngen,
                            germline.testing=germline.testing,
                            marker.testing=marker.testing, psize=psize, 
                            counselee.id=counselee.id, allef=allef)
    
  }
  
  if (missingage & imputeAges==TRUE){
    agewarning = "Warning: Unknown ages of some unaffected and affected family members have been imputed.  You may want to get more information about family member ages and re-run the calculation."
    ifamily = ImputeAge.gast(fff=family, params=params, model="MMRpro") #ImputeAge is run only once and returns a matrix of nIter possible affection ages
    postList = NULL
    nIter = params$nIter
    
    # Run each imputed age through runPeelingMMRpro if there are affected relatives
    # who had age imputed
    if (length(ifamily$mem_bc)>0 | length(ifamily$mem_oc)>0 | length(ifamily$mem_gc) > 0) {
      for (iIter in 1:nIter){
        # creating imputed fff
        nf = ifamily$fff
        
        if (length(ifamily$mem_bc)>0) {
          nf$AgeColon[ifamily$mem_bc] = ifamily$age_bc[iIter,ifamily$mem_bc]  
        }
        
        if (length(ifamily$mem_oc)>0) {
          nf$AgeEndometrium[ifamily$mem_oc]<-ifamily$age_oc[iIter,ifamily$mem_oc]  
        }
        
        if (length(ifamily$mem_gc)>0) {
          nf$AgeGastric[ifamily$mem_gc]<-ifamily$age_gc[iIter,ifamily$mem_gc]  
        }
        
        postList = rbind(postList,
                         runPeelingMMRpro(family=nf, 
                                          params=params, ngen=ngen,
                                          germline.testing=germline.testing,
                                          marker.testing=marker.testing, psize=psize, 
                                          counselee.id=counselee.id, allef=allef)
        )
      }
      
      post = apply(postList, 2, mean)
    }
    
    if (length(ifamily$mem_bc)==0 & length(ifamily$mem_oc)==0 & length(ifamily$mem_gc) == 0) {
      post = runPeelingMMRpro(family=ifamily$fff, 
                              params=params, ngen=ngen,
                              germline.testing=germline.testing,
                              marker.testing=marker.testing, psize=psize, 
                              counselee.id=counselee.id, allef=allef)
    }
    
  }
  
  if (!missingage | imputeAges==FALSE) { # No missing ages or impute flag is turned off, run peeling once with observed input data.
    agewarning=""
    post = runPeelingMMRpro(family=family, 
                            params=params, ngen=ngen,
                            germline.testing=germline.testing,
                            marker.testing=marker.testing, psize=psize, 
                            counselee.id=counselee.id, allef=allef)
    
  }
  
  
  #### Age imputation loop would stop here, output results
  #### ---------------------------------------------------####
  
  
  # save family object to be output
  if (!missingage | imputeAges==FALSE) {family= family}
  if (missingage & imputeAges==TRUE) { family = ifamily$fff}
  
  
  ll   <- NULL
  
  post <- array(post, rep(3, 3))
  dn <- list(0)
  dn[[1]] <- c("MLH10", "MLH11", "MLH12")
  dn[[2]] <- c("MSH20", "MSH21", "MSH22")
  dn[[3]] <- c("MSH60", "MSH61", "MSH62")
  dimnames(post) <- dn
  
  marg <- rep(0,3)
  marg[1] <- 1-sum(post[1,,])
  marg[2] <- 1-sum(post[,1,])
  marg[3] <- 1-sum(post[,,1])
  # 
  # if (print) {
  #   cat("The probability of being a carrier is",1-post[1,1,1], "\n an MLH1 carrier", marg[1], "\n an MSH2 carrier", marg[2], "\n an MSH6 carrier", marg[3],"\n")
  # }
  
  # current.age <- maxage[which(family$ID==originalcounselee.id)]
  # if(proband.current.age + params$age.by > age.max){
  #   ages <- age.max
  # } else{
  #   ages <- seq(proband.current.age+params$age.by, ifelse(params$age.to>=proband.current.age+params$age.by, params$age.to, age.max), by=params$age.by)
  # }
  
  # if(pro.agemax == 1){
  #   predictions <- future.risk <- data.frame()
  #   warning(paste("Proband is older than ", age.max - 1, ". Future risk prediction is not calculated.", sep = ""))
  # } else{
  #   
  #   predictions <- NULL
  #   ages.by1 <- (current.age+1):max(ages)
  #   
  #   
  #   if(family[which(family$ID==originalcounselee.id),"Gender"]==0){#female
  #     hFX0 <- hFX1 <- hFX2 <- hFX6 <- hFY0 <- hFY1 <- hFY2 <- hFY6 <- rep(0,max(ages))
  #     future.risk <- data.frame(hFX0,hFX1,hFX2,hFX6,hFY0,hFY1,hFY2,hFY6)
  #     
  #     if(net == FALSE){
  #       ## crude risk
  #       future.risk$hFX0[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFX[,"M000"],haz.d$femaleCRC)[ages.by1],haz.d$femaleCRC[ages.by1])
  #       future.risk$hFY0[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFY[,"M000"],haz.d$uterus)[ages.by1],haz.d$uterus[ages.by1])
  #       future.risk$hFZ0[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFZ[,"M000"],haz.d$femaleGC)[ages.by1],haz.d$femaleGC[ages.by1])
  #       future.risk$hFX1[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFX[,"M100"],haz.d$femaleCRC)[ages.by1],haz.d$femaleCRC[ages.by1])
  #       future.risk$hFY1[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFY[,"M100"],haz.d$uterus)[ages.by1],haz.d$uterus[ages.by1])
  #       future.risk$hFZ1[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFZ[,"M100"],haz.d$femaleGC)[ages.by1],haz.d$femaleGC[ages.by1])
  #       future.risk$hFX2[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFX[,"M010"],haz.d$femaleCRC)[ages.by1],haz.d$femaleCRC[ages.by1])
  #       future.risk$hFY2[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFY[,"M010"],haz.d$uterus)[ages.by1],haz.d$uterus[ages.by1])
  #       future.risk$hFZ2[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFZ[,"M010"],haz.d$femaleGC)[ages.by1],haz.d$femaleGC[ages.by1])
  #       future.risk$hFX6[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFX[,"M001"],haz.d$femaleCRC)[ages.by1],haz.d$femaleCRC[ages.by1])
  #       future.risk$hFY6[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFY[,"M001"],haz.d$uterus)[ages.by1],haz.d$uterus[ages.by1])
  #       future.risk$hFZ6[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance.crude$fFZ[,"M001"],haz.d$femaleGC)[ages.by1],haz.d$femaleGC[ages.by1])
  #     } else{
  #       ## net risk
  #       future.risk$hFX0[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFX[,"M000"], proband.current.age, ages.by1)
  #       future.risk$hFY0[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFY[,"M000"], proband.current.age, ages.by1)
  #       future.risk$hFZ0[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFZ[,"M000"], proband.current.age, ages.by1)
  #       future.risk$hFX1[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFX[,"M100"], proband.current.age, ages.by1)
  #       future.risk$hFY1[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFY[,"M100"], proband.current.age, ages.by1)
  #       future.risk$hFZ1[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFZ[,"M100"], proband.current.age, ages.by1)
  #       future.risk$hFX2[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFX[,"M010"], proband.current.age, ages.by1)
  #       future.risk$hFY2[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFY[,"M010"], proband.current.age, ages.by1)
  #       future.risk$hFZ2[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFZ[,"M010"], proband.current.age, ages.by1)
  #       future.risk$hFX6[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFX[,"M001"], proband.current.age, ages.by1)
  #       future.risk$hFY6[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFY[,"M001"], proband.current.age, ages.by1)
  #       future.risk$hFZ6[ages.by1] <- calc.future.risk.nd(params$penetrance.net$fFZ[,"M001"], proband.current.age, ages.by1)
  #     }
  #     
  #     
  #     X.risk <- future.risk$hFX1[ages]*marg[1] + future.risk$hFX2[ages]*marg[2] + future.risk$hFX6[ages]*marg[3] + future.risk$hFX0[ages]*post[1,1,1]
  #     Y.risk <- future.risk$hFY1[ages]*marg[1] + future.risk$hFY2[ages]*marg[2] + future.risk$hFY6[ages]*marg[3] + future.risk$hFY0[ages]*post[1,1,1]
  #     Z.risk <- future.risk$hFZ1[ages]*marg[1] + future.risk$hFZ2[ages]*marg[2] + future.risk$hFZ6[ages]*marg[3] + future.risk$hFZ0[ages]*post[1,1,1]
  #     
  #   }else{#male
  #     hMX0 <- hMX1 <- hMX2 <- hMX6 <- hMY0 <- hMY1 <- hMY2 <- hMY6 <- rep(0,max(ages))
  #     future.risk <- data.frame(hMX0,hMX1,hMX2,hMX6,hMY0,hMY1,hMY2,hMY6)
  #     
  #     if(net == FALSE){
  #       ## crude risk
  #       future.risk$hMX0[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMX[,"M000"],haz.d$maleCRC)[ages.by1],haz.d$maleCRC[ages.by1])
  #       future.risk$hMX1[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMX[,"M100"],haz.d$maleCRC)[ages.by1],haz.d$maleCRC[ages.by1])
  #       future.risk$hMX2[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMX[,"M010"],haz.d$maleCRC)[ages.by1],haz.d$maleCRC[ages.by1])
  #       future.risk$hMX6[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMX[,"M001"],haz.d$maleCRC)[ages.by1],haz.d$maleCRC[ages.by1])
  #       future.risk$hMZ0[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMZ[,"M000"],haz.d$maleGC)[ages.by1],haz.d$maleGC[ages.by1])
  #       future.risk$hMZ1[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMZ[,"M100"],haz.d$maleGC)[ages.by1],haz.d$maleGC[ages.by1])
  #       future.risk$hMZ2[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMZ[,"M010"],haz.d$maleGC)[ages.by1],haz.d$maleGC[ages.by1])
  #       future.risk$hMZ6[ages.by1] <- calc.future.risk(dens2haz.cs(params$penetrance$fMZ[,"M001"],haz.d$maleGC)[ages.by1],haz.d$maleGC[ages.by1])
  #     } else{
  #       ## net risk
  #       future.risk$hMX0[ages.by1] <- calc.future.risk.nd(params$penetrance$fMX[,"M000"], proband.current.age, ages.by1)
  #       future.risk$hMX1[ages.by1] <- calc.future.risk.nd(params$penetrance$fMX[,"M100"], proband.current.age, ages.by1)
  #       future.risk$hMX2[ages.by1] <- calc.future.risk.nd(params$penetrance$fMX[,"M010"], proband.current.age, ages.by1)
  #       future.risk$hMX6[ages.by1] <- calc.future.risk.nd(params$penetrance$fMX[,"M001"], proband.current.age, ages.by1)
  #       future.risk$hMZ0[ages.by1] <- calc.future.risk.nd(params$penetrance$fMZ[,"M000"], proband.current.age, ages.by1)
  #       future.risk$hMZ1[ages.by1] <- calc.future.risk.nd(params$penetrance$fMZ[,"M100"], proband.current.age, ages.by1)
  #       future.risk$hMZ2[ages.by1] <- calc.future.risk.nd(params$penetrance$fMZ[,"M010"], proband.current.age, ages.by1)
  #       future.risk$hMZ6[ages.by1] <- calc.future.risk.nd(params$penetrance$fMZ[,"M001"], proband.current.age, ages.by1)
  #     }
  #     
  #     X.risk <- future.risk$hMX1[ages]*marg[1] + future.risk$hMX2[ages]*marg[2] + future.risk$hMX6[ages]*marg[3] + future.risk$hMX0[ages]*post[1,1,1] 
  #     Y.risk <- rep(0, length(ages))
  #     Z.risk <- future.risk$hMZ1[ages]*marg[1] + future.risk$hMZ2[ages]*marg[2] + future.risk$hMZ6[ages]*marg[3] + future.risk$hMZ0[ages]*post[1,1,1] 
  #     
  #   }
  #   
  #   # If proband already has colorectal or endometrial cancer, set future risk estimates to NA
  #   if (family$AffectedColon[which(family$ID==originalcounselee.id)]==1) {X.risk <- rep(NA, length(ages))}
  #   if (family$AffectedEndometrium[which(family$ID==originalcounselee.id)]==1) {Y.risk <- rep(NA, length(ages))}
  #   if (family$AffectedGastric[which(family$ID==originalcounselee.id)]==1) {Z.risk <- rep(NA, length(ages))}
  #   
  #   
  #   predictions <- data.frame(ages, colorectal.risk=X.risk, endometrial.risk=Y.risk, gastric.risk = Z.risk)
  #   colnames(predictions) <- c("By age", "Colorectal Ca Risk", "Endometrial Ca Risk", "Gastric Ca Risk")
  #   
  #   if (print) {
  #     cat("The risks of developing cancers are", "\n")
  #     print(predictions)
  #   }
  # }
  predictions <- data.frame()
  future.risk <- data.frame()
  
  probs <- data.frame(non=1-post[1,1,1],MLH1=marg[1],MSH2=marg[2],MSH6=marg[3])
  names(probs) <- c("Pr(Being a carrier)", "Pr(MLH1 mutation)", "Pr(MSH2 mutation)", "Pr(MSH6)")
  
  
  return(new(Class="BayesMendel", family=family, posterior=post, probs=probs, predictions=predictions, counselee.id=originalcounselee.id, loglik=ll, future.risk=future.risk))
}

