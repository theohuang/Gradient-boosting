CheckFamStructure.gast  <-  function(model,fff,counselee.id,germline.testing=NULL,marker.testing=NULL,
                                     oophorectomy=NULL, mastectomy=NULL, 
                                     imputeAges, imputeRelatives, params){
  
  
  #### THIS IS VERSION 2.1-2, adding Lyte simple and correcting age replacements
  
  #### THIS IS VERSION 2.1, removing twins collapse ####
  
  #### THIS IS VERSION 2.0-812, new input counselee.id ###
  
  # counselee.id is required for the brcapro model.
  # If only the proband has ethnicity of 'AJ', then all NAs are set to 'AJ'.
  # Otherwise, NAs are set to 'nonAJ'
  
  #### THIS IS VERSION 2.0-4, addition of MELAPRO ###
  
  # Check for obvious errors in family structure and disease reporting
  # turn vector into matrix for pedigrees with one person
  #
  # Columns needed in the basic family structure
  # BRCAPRO: ID, Gender, FatherID, MotherID, AffectedBreast, AffectedOvary, AgeBreast, AgeOvary,
  # AgeBreastContralateral, Twins
  #
  # MMRPRO: ID, Gender, FatherID, MotherID, AffectedColon, AffectedEndometrium, AgeColon, AgeEndometrium
  #
  # PANCPRO: ID, Gender, FatherID, MotherID, AffectedPancreas, AgePancreas
  #
  
  # Additional columns:
  # germline.testing -- BRCA1, BRCA2 for BRCAPRO
  #                     MLH1, MSH2 for MMRpro
  # marker.testing -- ER, CK14, CK5.6, PR for BRCAPRO
  #                -- MSI, location for MMRpro
  # oophorectomy -- Oophorectomy, AgeOophorectomy
  #
  # mastectomy -- Mastectomy, AgeMastectomy
  
  ## version 2.1-5 : now we only use the other cancer age to impute missing ages if
  ## the individual is unaffected for both cancers (i.e., the ages are current ages and not cancer ages)
  ## also for 2.1-5, if the user inputs germline.testing or marker.testing information
  ## and the family already has that information in the columns, we now delete
  ## the columns and use the inputted information (in case there are inconsistencies between the two)
  
  errormsg <- ""
  warnmsg <- ""
  psize <- nrow(fff)
  
  age.max <- 94
  
  #--------------------------------------------------------------------#
  # New to version 2.1: allowing "NA" for AffectedX and AgeX variables #
  
  if (model=="brcapro"){
    UnkBrAffAge1 = which(fff$AffectedBreast==1 & (is.na(fff$AgeBreast) | fff$AgeBreast==1))
    if (length(UnkBrAffAge1)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkBrAffAge1], " has unknown breast cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnkBrAffAge2 = which(fff$AffectedBreast==2 & (is.na(fff$AgeBreast) | is.na(fff$AgeBreastContralateral) | fff$AgeBreast==1 | fff$AgeBreastContralateral==1 | fff$AgeBreastContralateral==0))
    if (length(UnkBrAffAge2)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkBrAffAge2], " has unknown breast cancer diagnosis age(s).  Please input a diagnosis age(s) or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnkOvAffAge = which(fff$AffectedOvary==1 & (is.na(fff$AgeOvary) | fff$AgeOvary==1))
    if (length(UnkOvAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkOvAffAge], " has unknown ovarian cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnknownStatus = which(is.na(fff$AffectedBreast) & is.na(fff$AffectedOvary) & (!is.na(fff$AgeBreast) | !is.na(fff$AgeOvary) | !is.na(fff$AgeBreastContralateral)))
    fff$AffectedBreast[is.na(fff$AffectedBreast)] = 0
    fff$AffectedOvary[is.na(fff$AffectedOvary)] = 0
    fff$AgeBreast[is.na(fff$AgeBreast)] = 1
    fff$AgeOvary[is.na(fff$AgeOvary)] = 1
    if (!is.null(mastectomy)) {mastectomy$AgeMastectomy[is.na(mastectomy$AgeMastectomy)] = 1 }
    if (!is.null(oophorectomy)) {oophorectomy$AgeOophorectomy[is.na(oophorectomy$AgeOophorectomy)] = 1}
    fff$AgeBreastContralateral[is.na(fff$AgeBreastContralateral) & (fff$AffectedBreast==0 | fff$AffectedBreast==1)] = 0
    fff$AgeBreastContralateral[is.na(fff$AgeBreastContralateral) & fff$AffectedBreast==2] = 1
    #fff$AgeBreast[UnknownStatus] = fff$AgeOvary[UnknownStatus] = 1
    #fff$AgeBreastContralateral[UnknownStatus] = 0
    
    # Go through some simpler age checks here first, and output a warning message about changes that were made
    # If one age parameter is missing but the other is known, replace with this age.
    if (any(fff$AgeBreast==1 & fff$AffectedBreast == 0 & fff$AgeOvary>1 & fff$AffectedOvary == 0)) {
      paste(warnmsg, "Warning: AgeBreast=1 has been replaced with AgeOvary, which was input as greater than 1, for at least one family member.")
    }
    fff$AgeBreast[fff$AgeBreast==1 & fff$AffectedBreast == 0 & fff$AgeOvary>1 & fff$AffectedOvary == 0] = fff$AgeOvary[fff$AgeBreast==1 & fff$AffectedBreast == 0 & fff$AgeOvary>1 & fff$AffectedOvary == 0]
    
    if (any(fff$AgeBreast>1 & fff$AffectedBreast == 0 & fff$AgeOvary==1 & fff$AffectedOvary == 0)) {
      paste(warnmsg, "Warning: AgeOvary=1 has been replaced with AgeBreast, which was input as greater than 1, for at least one family member.")
    }
    fff$AgeOvary[fff$AgeBreast>1 & fff$AffectedBreast == 0 & fff$AgeOvary==1 & fff$AffectedOvary == 0] = fff$AgeBreast[fff$AgeBreast>1 & fff$AffectedBreast == 0 & fff$AgeOvary==1 & fff$AffectedOvary == 0]
    
    if (any(fff$AffectedBreast>0 & fff$AgeBreast > fff$AgeOvary & fff$AffectedOvary==0)) {
      paste(warnmsg, "Warning: at least one family member with breast cancer has a current age that is less than the breast cancer diagnosis age.  Current age (AgeOvary) has been replaced with breast cancer age (AgeBreast).")
    }
    fff$AgeOvary[fff$AffectedBreast>0 & fff$AgeBreast > fff$AgeOvary & fff$AffectedOvary==0] = fff$AgeBreast[fff$AffectedBreast>0 & fff$AgeBreast > fff$AgeOvary & fff$AffectedOvary==0]
    
    if (any(fff$AffectedOvary==1 & fff$AgeBreast < fff$AgeOvary & fff$AffectedBreast==0)) {paste(warnmsg, "Warning: at least one family member with ovarian cancer has a current age that is less than the ovarian cancer diagnosis age.  Current age (AgeBreast) has been replaced with ovarian cancer age (AgeOvary).")}
    fff$AgeBreast[fff$AffectedOvary==1 & fff$AgeBreast < fff$AgeOvary & fff$AffectedBreast==0] = fff$AgeOvary[fff$AffectedOvary==1 & fff$AgeBreast < fff$AgeOvary & fff$AffectedBreast==0]
    
    # If mastectomy age is missing, replace with age breast or age contralateral
    # If both are missing, impute 
    if (!is.null(mastectomy)) {
      if (any((mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2)) {paste(warnmsg, "Warning: Unknown mastectomy age for those with mastectomy has been replaced with current age, via AgeBreast, which also may have been imputed, see other warnings.")}
      mastectomy$AgeMastectomy[(mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2] = fff$AgeBreast[(mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2]
      
      if (any((mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1)) {paste(warnmsg, "Warning: Unknown mastectomy age for those with mastectomy has been replaced with age of contralateral breast cancer, via AgeBreastContralateral, which also may have been imputed, see other warnings.")}
      mastectomy$AgeMastectomy[(mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1] = fff$AgeBreastContralateral[(mastectomy$AgeMastectomy==1|mastectomy$AgeMastectomy==0) & mastectomy$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1]
    }
    
    # If oophorectomy age is missing, replace with age ovary
    # If both are missing, impute 
    if (!is.null(oophorectomy)) {
      if (any((oophorectomy$AgeOophorectomy==1 | oophorectomy$AgeOophorectomy==0) & oophorectomy$Oophorectomy==1 & fff$AgeOvary!=1)) {paste(warnmsg, "Warning: Unknown oophorectomy age for those with oophorectomy has been replaced with current age, via AgeOvary, which also may have been imputed, see other warnings.")}
      oophorectomy$AgeOophorectomy[(oophorectomy$AgeOophorectomy==1 | oophorectomy$AgeOophorectomy==0) & oophorectomy$Oophorectomy==1 & fff$AgeOvary!=1] = fff$AgeOvary[(oophorectomy$AgeOophorectomy==1 | oophorectomy$AgeOophorectomy==0) & oophorectomy$Oophorectomy==1 & fff$AgeOvary!=1]
    }    
    
    if (length(UnknownStatus)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[UnknownStatus], " have a known age input but unknown cancer status.  The calculation would be improved with an updated affection status for these individuals.", sep="")
    }
  }
  
  if (model=="MMRpro"){
    UnkCRCAffAge = which(fff$AffectedColon==1 & (is.na(fff$AgeColon) | fff$AgeColon==1))
    if (length(UnkCRCAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkCRCAffAge], " has unknown colorectal cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnkECAffAge = which(fff$AffectedEndometrium==1 & (is.na(fff$AgeEndometrium) | fff$AgeEndometrium==1))
    if (length(UnkECAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkECAffAge], " has unknown endometrial cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnkGCAffAge = which(fff$AffectedGastric==1 & (is.na(fff$AgeGastric) | fff$AgeGastric==1))
    if (length(UnkGCAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkGCAffAge], " has unknown gastric cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnknownStatus = which(is.na(fff$AffectedColon) & is.na(fff$AffectedEndometrium) & is.na(fff$AffectedGastric) &
                            (!is.na(fff$AgeColon) | !is.na(fff$AgeEndometrium) | !is.na(fff$AgeGastric)))
    AffectedAge1 = which((fff$AffectedColon==1 & (is.na(fff$AgeColon) | fff$AgeColon == 1)) |
                           (fff$AffectedEndometrium==1 & (is.na(fff$AgeEndometrium) | fff$AgeEndometrium == 1)) |
                           (fff$AffectedGastric == 1 & (is.na(fff$AgeGastric) | fff$AgeGastric == 1)))
    fff$AffectedColon[is.na(fff$AffectedColon)] = 0
    fff$AffectedEndometrium[is.na(fff$AffectedEndometrium)] = 0
    fff$AffectedGastric[is.na(fff$AffectedGastric)] = 0
    fff$AgeColon[is.na(fff$AgeColon)] = 1
    fff$AgeEndometrium[is.na(fff$AgeEndometrium)] = 1
    fff$AgeGastric[is.na(fff$AgeGastric)] = 1
    fff$AgeColon[UnknownStatus] = fff$AgeEndometrium[UnknownStatus] = fff$AgeGastric[UnknownStatus] = 1
    if (length(UnknownStatus)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[UnknownStatus], " have a known age input but unknown cancer status.  The calculation would be improved with an updated affection status for these individuals.", sep="")
    }
    if (length(AffectedAge1)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[AffectedAge1], " have an unknown age of diagnosis.  The calculation would be improved with an estimate of the affection age included.", sep="")
    }
    
    if (any(fff$AgeColon==1 & fff$AffectedColon == 0 & fff$AgeEndometrium>1 & fff$AffectedEndometrium == 0)) {paste(warnmsg, "Warning: AgeColon=1 has been replaced with AgeEndometrium, which was input as greater than 1, for at least one family member.")}
    fff$AgeColon[fff$AgeColon==1 & fff$AffectedColon == 0 & fff$AgeEndometrium>1 & fff$AffectedEndometrium == 0] = fff$AgeEndometrium[fff$AgeColon==1 & fff$AffectedColon == 0 & fff$AgeEndometrium>1 & fff$AffectedEndometrium == 0]
    if (any(fff$AgeColon==1 & fff$AffectedColon == 0 & !(fff$AgeEndometrium & fff$AffectedEndometrium == 0) & fff$AgeGastric > 1 & fff$AffectedGastric == 0)) {paste(warnmsg, "Warning: AgeColon=1 has been replaced with AgeEndometrium, which was input as greater than 1, for at least one family member.")}
    fff$AgeColon[fff$AgeColon==1 & fff$AffectedColon == 0 & !(fff$AgeEndometrium & fff$AffectedEndometrium == 0) & fff$AgeGastric > 1 & fff$AffectedGastric == 0] = fff$AgeGastric[fff$AgeColon==1 & fff$AffectedColon == 0 & !(fff$AgeEndometrium & fff$AffectedEndometrium == 0) & fff$AgeGastric > 1 & fff$AffectedGastric == 0]
    
    
    
    if (any(fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0)) {paste(warnmsg, "Warning: AgeEndometrium=1 has been replaced with AgeColon, which was input as greater than 1, for at least one family member.")}
    fff$AgeEndometrium[fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0] = fff$AgeColon[fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0]
    if (any(!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0 &
            fff$AgeGastric > 1 & fff$AgeGastric == 0)) {paste(warnmsg, "Warning: AgeEndometrium=1 has been replaced with AgeColon, which was input as greater than 1, for at least one family member.")}
    fff$AgeEndometrium[!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0 &
                         fff$AgeGastric > 1 & fff$AgeGastric == 0] = fff$AgeGastric[!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeEndometrium==1 & fff$AffectedEndometrium == 0 &
                                                                                      fff$AgeGastric > 1 & fff$AgeGastric == 0]
    
    if (any(fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeGastric==1 & fff$AffectedGastric == 0)) {paste(warnmsg, "Warning: AgeGastric=1 has been replaced with AgeColon, which was input as greater than 1, for at least one family member.")}
    fff$AgeGastric[fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeGastric==1 & fff$AffectedGastric == 0] = fff$AgeColon[fff$AgeColon>1 & fff$AffectedColon == 0 & fff$AgeGastric==1 & fff$AffectedGastric == 0]
    if (any(!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeGastric==1 & fff$AffectedEGastric == 0 &
            fff$AgeEndometrium > 1 & fff$AgeEndometrium == 0)) {paste(warnmsg, "Warning: AgeGastric=1 has been replaced with AgeEndometrium, which was input as greater than 1, for at least one family member.")}
    fff$AgeGastric[!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeGastric==1 & fff$AffectedEGastric == 0 &
                     fff$AgeEndometrium > 1 & fff$AgeEndometrium == 0] = fff$AgeEndometrium[!(fff$AgeColon>1 & fff$AffectedColon == 0) & fff$AgeGastric==1 & fff$AffectedEGastric == 0 &
                                                                                              fff$AgeEndometrium > 1 & fff$AgeEndometrium == 0]
    
    
    
  }
  if (model=="pancpro"){
    UnkPancAffAge = which(fff$AffectedPancreas==1 & (is.na(fff$AgePancreas) | fff$AgePancreas==1))
    if (length(UnkPancAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkPancAffAge], " has unknown pancreas cancer diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnknownStatus = which(is.na(fff$AffectedPancreas) & !is.na(fff$AgePancreas))
    AffectedAge1 = which(fff$AffectedPancreas==1 & is.na(fff$AgePancreas))
    fff$AffectedPancreas[is.na(fff$AffectedPancreas)] = 0
    fff$AgePancreas[is.na(fff$AgePancreas)] = 1
    fff$AgePancreas[UnknownStatus] = 1
    if (length(UnknownStatus)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[UnknownStatus], " have a known age input but unknown cancer status.  The calculation would be improved with an updated affection status for these individuals.", sep="")
    }
    if (length(AffectedAge1)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[AffectedAge1], " have an unknown age of diagnosis.  The calculation would be improved with an estimate of the affection age included.", sep="")
    }
  }
  if (model=="melapro"){
    UnkSkinAffAge = which(fff$AffectedSkin==1 & (is.na(fff$AgeSkin) | fff$AgeSkin==1))
    if (length(UnkSkinAffAge)>0 & imputeAges==FALSE){
      errormsg <- paste("Error: ID ", fff$ID[UnkSkinAffAge], " has unknown melanoma diagnosis age.  Please input a diagnosis age or turn on age imputation with imputeAges=TRUE.")
      return(errormsg)
    }
    UnknownStatus = which(is.na(fff$AffectedSkin) & !is.na(fff$AgeSkin))
    AffectedAge1 = which(fff$AffectedSkin>0 & is.na(fff$AgeSkin))
    fff$AffectedSkin[is.na(fff$AffectedSkin)] = 0
    fff$AgeSkin[is.na(fff$AgeSkin)] = 1
    fff$AgeSkin[UnknownStatus] = 1
    if (length(UnknownStatus)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[UnknownStatus], " have a known age input but unknown cancer status.  The calculation would be improved with an updated affection status for these individuals.", sep="")
    }
    if (length(AffectedAge1)>0) {
      warnmsg = paste(warnmsg, "Warning: individuals", fff$id[AffectedAge1], " have an unknown age of diagnosis.  The calculation would be improved with an estimate of the affection age included.", sep="")
    }
  }
  
  
  #--------------------------------------------------------------#
  # Check if the inputs are data frames and have the right names #
  
  if(model=="brcapro"){
    family.names <- c("ID", "Gender", "FatherID", "MotherID", "AffectedBreast", 
                      "AffectedOvary", "AgeBreast", "AgeOvary", 
                      "AgeBreastContralateral", "Twins", "ethnic")
    gene.names <- c("BRCA1","BRCA2","TestOrder")
    marker.names <- c("ER","CK14","CK5.6","PR","HER2")
    oophorectomy.names <- c("Oophorectomy","AgeOophorectomy")
    mastectomy.names <- c("Mastectomy", "AgeMastectomy")
    
  }
  
  if(model=="MMRpro"){
    family.names <- c("ID", "Gender", "FatherID", "MotherID", "AffectedColon", 
                      "AffectedEndometrium", "AffectedGastric",
                      "AgeColon", "AgeEndometrium", "AgeGastric",
                      "Twins", "ethnic")
    gene.names <- c("MLH1","MSH2","MSH6","TestOrder")
    marker.names <- c("MSI", "location")
  }
  
  if(model=="pancpro"){
    family.names <- c("ID", "Gender", "FatherID", "MotherID", "AffectedPancreas", 
                      "AgePancreas", "Twins", "ethnic")
    gene.names <- "PANC"
  }
  
  if(model=="melapro"){
    family.names <- c("ID", "Gender", "FatherID", "MotherID", "AffectedSkin", 
                      "AgeSkin", "Twins", "ethnic")
    gene.names <- c("P16","TestOrder")
  }
  
  
  #--------------#
  # Family input #
  
  if(!is.data.frame(fff)){
    errormsg <- "Error: family history is not a data frame"
    return(errormsg)
  }
  names.input <- names(fff)
  temp <- family.names%in%names.input
  if(length(family.names)!=sum(temp)){
    missing.columns <- family.names[!temp]
    errormsg <- paste("Error: column ", missing.columns," is missing.",sep="")
    return(errormsg)
  }
  
  psize <- dim(fff)[1]  # get the number of individuals in the family, need to check this later for other inputs
  
  #--------------------------------#
  # Ethnicity input (brcapro only) #
  if (model=="brcapro") {
    ethnicity.names <- c("AJ", "nonAJ", "Other", "Italian")
    
    # If at least one member is AJ, set all missing fam members to AJ.
    if(class(fff$ethnic) == 'factor') fff$ethnic <- as.character(fff$ethnic)
    if(sum(fff$ethnic=="AJ",na.rm=TRUE)>0) {fff$ethnic[is.na(fff$ethnic)] <- "AJ"} else {fff$ethnic[is.na(fff$ethnic)] <- "nonAJ"}
    
    if (length(setdiff(fff$ethnic, ethnicity.names))>0) {
      errormsg <- "Error: invalid input for ethnicity.  See help(brcapro) for options allowed."
      return(errormsg)
    }
  }
  
  #------------------------# 
  # Germline.testing input #
  
  if(!is.null(germline.testing)){
    if(!is.data.frame(germline.testing)){
      errormsg <- "Error: germline.testing is not a data frame."
      return(errormsg)
    }
    if(dim(germline.testing)[1]!=psize){
      errormsg <- "Error: the number of rows for germline.testing is different from family input."
      return(errormsg)
    }
    names.input <- names(germline.testing)
    temp <- gene.names%in%names.input
    if(length(gene.names)!=sum(temp)){
      missing.columns <- gene.names[!temp]
      errormsg <- paste("Error: column ",missing.columns," is missing.",sep="")
      return(errormsg)
    }
    
    # if fff already has columns with the genetic testing information, replace it with
    # the germline.testing data frame information
    fff[, gene.names] <- NULL
    fff <- data.frame(fff,germline.testing)
  }
  
  #----------------------# 
  # marker.testing input #
  
  if(!is.null(marker.testing)){
    if(!is.data.frame(marker.testing)){
      errormsg <- "Error: marker.testing is not a data frame."
      return(errormsg)
    }
    if(dim(marker.testing)[1]!=psize){
      errormsg <- "Error: the number of rows for marker.testing is different from family input."
      return(errormsg)
    }
    names.input <- names(marker.testing)
    temp <- marker.names%in%names.input
    if(length(marker.names)!=sum(temp)){
      missing.columns <- marker.names[!temp]
      errormsg <- paste("Error: column ",missing.columns," is missing.",sep="")
      return(errormsg)
    }
    
    # if fff already has columns with the marker testing information, replace it with
    # the marker.testing data frame information
    fff[, marker.names] <- NULL
    fff <- data.frame(fff,marker.testing)
  }
  
  #--------------------#
  # oophorectomy input #
  
  if(!is.null(oophorectomy)){
    if(!is.data.frame(oophorectomy)){
      errormsg <- "Error: oophorectomy is not a data frame."
      return(errormsg)
    }
    if(dim(oophorectomy)[1]!=psize){
      errormsg <- "Error: the number of rows for oophorectomy is different from family input."
      return(errormsg)
    }
    names.input <- names(oophorectomy)
    temp <- oophorectomy.names%in%names.input
    if(length(oophorectomy.names)!=sum(temp)){
      missing.columns <- oophorectomy.names[!temp]
      errormsg <- paste("Error: column ",missing.columns," is missing.",sep="")
      return(errormsg)
    }
    # cannot have age=0 for interventions.
    oophorectomy$AgeOophorectomy[oophorectomy$AgeOophorectomy==0] = 1
    # cannot have unknown intervention ages for unaffected members with
    # unknown age
    # Removing this requirement as of v2.1-2 now that we impute age.
    #if (any(oophorectomy$Oophorectomy==1 & oophorectomy$AgeOophorectomy==1 & fff$AgeBreast==1 & fff$AgeOvary==1 & fff$AgeBreastContralateral==1)) {
    #  errormsg <- "Error: oophorectomy may not be included for individuals with no known age."
    #  return(errormsg)
    #}
    
    # if fff already has columns with the oophorectomy information, replace it with
    # the oophorectomy data frame information
    fff[, oophorectomy.names] <- NULL
    fff <- data.frame(fff,oophorectomy)
  }
  
  #--------------------#
  # Mastectomy input #
  
  if(!is.null(mastectomy)){
    if(!is.data.frame(mastectomy)){
      errormsg <- "Error: mastectomy is not a data frame."
      return(errormsg)
    }
    if(dim(mastectomy)[1]!=psize){
      errormsg <- "Error: the number of rows for mastectomy is different from family input."
      return(errormsg)
    }
    names.input <- names(mastectomy)
    temp <- mastectomy.names%in%names.input
    if(length(mastectomy.names)!=sum(temp)){
      missing.columns <- mastectomy.names[!temp]
      errormsg <- paste("Error: column ",missing.columns," is missing.",sep="")
      return(errormsg)
    }
    mastectomy$AgeMastectomy[mastectomy$AgeMastectomy==0] = 1    
    # cannot have unknown intervention ages for unaffected members with
    # unknown age
    #if (any(mastectomy$Mastectomy==1 & mastectomy$AgeMastectomy==1 & fff$AgeBreast==1 & fff$AgeOvary==1 & fff$AgeBreastContralateral==1)) {
    #  errormsg <- "Error: mastectomy may not be included for individuals with no known age."
    #  return(errormsg)
    #}
    
    # if fff already has columns with the mastectomy information, replace it with
    # the mastectomy data frame information
    fff[, mastectomy.names] <- NULL
    fff <- data.frame(fff,mastectomy)
  }
  
  #------------------------#
  # Biomarker checks #
  
  if((model=="brcapro" | model=="brcapancpro") & !is.null(marker.testing)){
    indx <- apply(marker.testing[,c("ER","CK14","CK5.6","PR")],1,sum)
    if(sum(indx >0 & fff[,"AffectedBreast"]==0)>0) {
      errormsg <-
        paste("At least 1 member in the family without breast cancer has breast cancer markers.")
      return(errormsg)
    }
    
    ## if ER is there, don't use PR
    for(i in 1:nrow(marker.testing)) 
      if(marker.testing[i,"ER"]!=0) marker.testing[i,"PR"] <- 0
      
  } 
  
  
  #----------------------------------------------#
  # Include checks specific to oophorectomy here #
  
  # Column names: Oophorectomy, OophorectomyAge, 15-Breast cancer
  # after oophorectomy, 16-Ovarian after oophorectomy, 17-BCage1 after oophorectomy,
  # 18-OCage after oophorectomy, 19-BCage2 after oophorectomy
  
  if((model=="brcapro") &!is.null(oophorectomy)){
    # Check valid oophorectomy status numbers
    if( any(fff[,"Oophorectomy"]!=0 & fff[,"Oophorectomy"]!=1) )
      errormsg <- paste(errormsg,"Error: oophorectomy status out of range.  ")
    
    # These people got oophorectomies
    oophors <- which(fff[,"Oophorectomy"]==1)
    
    # Only do this if there are any relatives with oophorectomy
    if ( length(oophors)>0 ) {
      
      # Subset out those with oophorectomy, note the age of oophorectomy
      # drop=F means do not squeeze out any singleton indices, i.e. don't coerce row
      # matrices to vectors
      oofff <- fff[oophors,,drop=FALSE]
      
      # Flag males with oophorectomy
      if(sum(oofff[,"Gender"]==1)>0) { #male members present
        if(sum(oofff[oofff[,"Gender"]==1, "Oophorectomy"])>0)
          errormsg <- paste(errormsg, "Error: male with oophorectomy.  ")
      }
      
      # Check for age consistency:
      # max(pre-oophorectomy ages) <= oophorectomy age <= all post-oophorectomy ages
      # Use drop=F in case only one person with oophorectomy: keep row matrix a matrix
      ages.cancer <- data.frame(oofff[,"AgeBreast",drop=FALSE],oofff[,"AgeOvary",drop=FALSE],oofff[,"AgeBreastContralateral",drop=FALSE])
      maxage <- apply(ages.cancer,1,max)
      ooage <- oofff[,"AgeOophorectomy",drop=FALSE]
      
      if ( all(maxage < ooage) )
        warnmsg <- paste(warnmsg,
                         "Warning: There is no event or censoring included after oophorectomy age.  This could make sense if the person has bilateral breast and ovarian cancer before oophorectomy.  If so, just drop the oophorectomy from this person's row because oophorectomy won't affect the calculations very much.")
      
      # You can't have oophorectomy after ovarian cancer, that's meaningless
      # Leaving this in --> if this happens, it's because the user put it in.
      if ( any(oofff[,"AffectedOvary"]==1 & oofff[,"AgeOvary"]<oofff[,"AgeOophorectomy"]) )
        errormsg <-
        paste(errormsg,"BRCAPRO only handles oophorectomy occurring before ovarian cancer, not after")
    }
  } # End oophorectomy checks
  
  #----------------------------------------------#
  # Include checks specific to mastectomy here #
  
  if((model=="brcapro") &!is.null(mastectomy)){
    # Check valid mastectomy status numbers
    if( any(fff[,"Mastectomy"]!=0 & fff[,"Mastectomy"]!=1) )
      errormsg <- paste(errormsg,"Error: mastectomy status out of range.  ")
    
    # These people got mastectomies
    mast <- which(fff[,"Mastectomy"]==1)
    
    # Only do this if there are any relatives with mastectomy
    if ( length(mast)>0 ) {
      
      # Subset out those with mastectomy, note the age of mastectomy
      # drop=F means do not squeeze out any singleton indices, i.e. don't coerce row
      # matrices to vectors
      mmm <- fff[mast,,drop=FALSE]
      
      # Check that those with mastectomy did not have prior bilateral cancer
      prior.cancer <- mmm[,"AffectedBreast",drop=FALSE]
      ages.cancer <- mmm[,"AgeBreastContralateral",drop=FALSE]
      mmage <- mmm[,"AgeMastectomy",drop=FALSE]
      
      if (any(ages.cancer < mmage & prior.cancer==2))
        errormsg <- paste(errormsg,
                          "Error: Mastectomy occurred after bilateral breast cancer. If so, just drop the mastectomy from this person's row because mastectomy won't affect the calculations very much.")
    }
    
  } # End mastectomy checks
  
  
  if(model=="brcapro"){
    # Winsorize ages at age.max
    ages.cancer <- data.frame(fff[,"AgeBreast",drop=FALSE],fff[,"AgeOvary",drop=FALSE],fff[,"AgeBreastContralateral",drop=FALSE])
    if(sum(ages.cancer>age.max)>0) warnmsg <- paste(warnmsg, "Warning: age(s) older than age.max has been converted to age.max!")
    fff[fff[,"AgeBreast"]>age.max,"AgeBreast"] <- age.max
    fff[fff[,"AgeOvary"]>age.max,"AgeOvary"] <- age.max
    fff[fff[,"AgeBreastContralateral"]>age.max,"AgeBreastContralateral"] <- age.max
    
    # Change censoring ages of zero to one
    # Check that no ages of onset less than one
    fff[fff[,"AffectedBreast"]==0 & fff[,"AgeBreast"]==0,"AgeBreast"] <- 1
    fff[fff[,"AffectedOvary"]==0 & fff[,"AgeOvary"]==0,"AgeOvary"] <- 1
    if(sum(fff[,"AffectedBreast"]==1 & fff[,"AgeBreast"]==0)>0 | sum(fff[,"AffectedOvary"]==1 & fff[,"AgeOvary"]==0)>0 | sum(fff[,"AffectedBreast"]==2 & fff[,"AgeBreast"]==0)>0) errormsg <- paste(errormsg, "Error: cancer age at onset is 0!  ") 
    
    # Valid cancer status and genotype status numbers
    if(sum((fff[,"AffectedBreast"]!=0 & fff[,"AffectedBreast"]!=1 & fff[,"AffectedBreast"]!=2) | (fff[,"AffectedOvary"]!=0 & fff[,"AffectedOvary"]!=1))>0)
      errormsg <- paste(errormsg,"Error: cancer status out of range.  ")
    
    if(!is.null(germline.testing)){
      if(sum((fff[,"BRCA1"]!=0&fff[,"BRCA1"]!=1&fff[,"BRCA1"]!=2)|
             (fff[,"BRCA2"]!=0&fff[,"BRCA2"]!=1&fff[,"BRCA2"]!=2))>0)
        errormsg <- paste(errormsg,"Error: test result out of range. ")
    }
    
    # If bilateral breast cancer, make sure the first age is less than the second age
    if ( any(fff[,"AffectedBreast"]==2) ) {
      bilaterals = which(fff[,"AffectedBreast"]==2)
      if ( any(fff[bilaterals,"AgeBreast"] > fff[bilaterals,"AgeBreastContralateral"] & fff[bilaterals,"AgeBreastContralateral"]!=1)) {
        errormsg <- paste(errormsg,"Error: Age of 1st breast cancer > Age of 2nd breast cancer.")}
      if ( any(fff[bilaterals,"AgeBreast"] > fff[bilaterals,"AgeBreastContralateral"] & 
               fff[bilaterals,"AgeBreastContralateral"][fff[bilaterals,"AgeBreast"] > fff[bilaterals,"AgeBreastContralateral"]] >1)) {
        errormsg <- paste(errormsg,"Error: Age of 1st breast cancer > Age of 2nd breast cancer.")
      }
    }
    
    if(sum(fff[,"Gender"]==1)>0) {#male members present
      if(sum(fff[fff[,"Gender"]==1, "AffectedOvary"], na.rm=T)>0) {
        errormsg <- paste(errormsg, "Error: male ovarian cancer.  ")
      }}
  }
  
  if(model=="MMRpro"){
    # Winsorize ages at age.max
    ages.cancer <- data.frame(fff[,"AgeColon",drop=FALSE],fff[,"AgeEndometrium",drop=FALSE], fff[,"AgeGastric",drop=FALSE])
    if(sum(ages.cancer>age.max)>0) warnmsg <- paste(warnmsg, "Warning: age(s) older than age.max has been converted to age.max!")
    fff[fff[,"AgeColon"]>age.max,"AgeColon"] <- age.max
    fff[fff[,"AgeEndometrium"]>age.max,"AgeEndometrium"] <- age.max
    fff[fff[,"AgeGastric"]>age.max,"AgeGastric"] <- age.max
    
    
    # Change censoring ages of zero to one
    # Check that no ages of onset less than one
    fff[fff[,"AffectedColon"]==0 & fff[,"AgeColon"]==0,"AgeColon"] <- 1
    fff[fff[,"AffectedEndometrium"]==0 & fff[,"AgeEndometrium"]==0,"AgeEndometrium"] <- 1
    fff[fff[,"AffectedGastric"]==0 & fff[,"AgeGastric"]==0,"AgeGastric"] <- 1
    
    if(sum(fff[,"AffectedColon"]==1 & fff[,"AgeColon"]==0)>0 |
       sum(fff[,"AffectedEndometrium"]==1 & fff[,"AgeEndometrium"]==0)>0 |
       sum(fff[,"AffectedGastric"]==1 & fff[,"AgeGastric"]==0)>0) errormsg <- paste(errormsg, "Error: cancer age at onset is 0!  ") 
    
    # Valid cancer status and genotype status numbers
    if(any(fff[,"AffectedEndometrium"]!=0 & fff[,"AffectedEndometrium"]!=1)) {
      errormsg <- paste(errormsg,"Error: cancer status out of range.  ")}
    
    if(!is.null(germline.testing)){
      if(sum((fff[,"MLH1"]!=0&fff[,"MLH1"]!=1&fff[,"MLH1"]!=2)|
             (fff[,"MSH2"]!=0&fff[,"MSH2"]!=1&fff[,"MSH2"]!=2))>0)
        errormsg <- paste(errormsg,"Error: test result out of range. ")
    }
    
    if(sum(fff[,"Gender"]==1)>0) #male members present
      if(sum(fff[fff[,"Gender"]==1, "AffectedEndometrium"])>0)
        errormsg <- paste(errormsg, "Error: male endometrial cancer.  ")
    
    # PREMM allows AffectedColon to be any integer, round down to 1 here.
    fff$AffectedColon[fff$AffectedColon>=1] = 1
  }
  
  if(model=="pancpro" ){
    # Winsorize ages at age.max
    ages.cancer <- fff[,"AgePancreas",drop=FALSE]
    if(sum(ages.cancer>age.max)>0) warnmsg <- paste(warnmsg, "Warning: age(s) older than age.max has been converted to age.max!")
    fff[fff[,"AgePancreas"]>age.max,"AgePancreas"] <- age.max
    
    # Change censoring ages of zero to one
    # Check that no ages of onset less than one
    fff[fff[,"AffectedPancreas"]==0 & fff[,"AgePancreas"]==0,"AgePancreas"] <- 1
    if(sum(fff[,"AffectedPancreas"]==1 & fff[,"AgePancreas"]==0)>0 ) errormsg <- paste(errormsg, "Error: cancer age at onset is 0!  ") 
    
    # Valid cancer status and genotype status numbers
    if(sum(fff[,"AffectedPancreas"]!=0 & fff[,"AffectedPancreas"]!=1 )>0)
      errormsg <- paste(errormsg,"Error: cancer status out of range.  ")
  }
  
  if(model=="melapro" ) {
    # Winsorize ages at age.max
    ages.cancer <- fff[,"AgeSkin",drop=FALSE]
    if(sum(ages.cancer>age.max)>0) warnmsg <- paste(warnmsg, "Warning: age(s) older than age.max has been converted to age.max!")
    fff[fff[,"AgeSkin"]>age.max,"AgeSkin"] <- age.max
    
    # Change censoring ages of zero to one
    # Check that no ages of onset less than one
    fff[fff[,"AffectedSkin"]==0 & fff[,"AgeSkin"]==0,"AgeSkin"] <- 1
    if(sum(fff[,"AffectedSkin"]==1 & fff[,"AgeSkin"]==0)>0 ) errormsg <- paste(errormsg, "Error: cancer age at onset is 0!  ") 
    
  }
  
  #------------------------------------------------------------#
  # Check that all mothers are female and all fathers are male #
  
  mother <- father <- rep(NA, psize)
  for (i in 1:psize) {
    mother[i] <- ifelse(length(intersect(fff$ID[i], fff$MotherID)>0), 1, 0)
    father[i] <- ifelse(length(intersect(fff$ID[i], fff$FatherID)>0), 1, 0)
  }
  
  if (any(fff$Gender[which(mother==1)]==1)) {errormsg <- paste(errormsg, "Error: mother is not female.  ")}
  if (any(fff$Gender[which(father==1)]==0)) {errormsg <- paste(errormsg, "Error: father is not male.  ")}
  
  #-------------------------------------------------------------------------------------------#
  # Check that all mother and father IDs listed in the pedigree exist as their own ID as well #
  
  # If someone has a mother or father who is not listed as an ID, set that mother or father ID value to zero.
  # Cannot refer to someone in mother/father ID who is not entered in the pedigree.
  
  mothers <- unique(fff$MotherID)
  fathers <- unique(fff$FatherID)
  ids <- unique(fff$ID)
  
  if (sum(setdiff(fathers, ids)) > 0) {
    zeroes <- which(setdiff(fathers, ids)==0)
    if (length(zeroes)>0) {who <- setdiff(fathers, ids)[-zeroes]} else{
      if (length(zeroes)==0) {who <- setdiff(fathers, ids)}}
    for (ww in 1:length(who)) {
      fff$FatherID[fff$FatherID==who[ww]] <- 0
    }
  }
  
  if (sum(setdiff(mothers, ids)) > 0) {
    zeroes <- which(setdiff(mothers, ids)==0)
    if (length(zeroes)>0) {who <- setdiff(mothers, ids)[-zeroes]} else{
      if (length(zeroes)==0) {who <- setdiff(mothers, ids)}}
    for (ww in 1:length(who)) {
      fff$MotherID[fff$MotherID==who[ww]] <- 0
    }
  }
  
  #------------------------------#
  # Check specific relationships #
  
  # First need to know who someone's first cousins are
  
  getotherrels <- function(ped, indid) {
    cid <- msibs <- psibs <- NULL
    mid <- ped[ped$ID==indid,"MotherID"]
    fid <- ped[ped$ID==indid,"FatherID"]
    
    if (mid!=0) {
      mgmid <- ped[ped$ID==mid,"MotherID"]
      mgfid <- ped[ped$ID==mid,"FatherID"]
      if (mgmid!=0 & mgfid!=0) 
      {msibs <- ped[(ped$FatherID==mgfid | ped$MotherID==mgmid) & ped$ID!=mid,"ID"]} #full siblings of the mother
      if (mgmid!=0 & mgfid==0) 
      {msibs <- ped[ped$MotherID==mgmid & ped$ID!=mid,"ID"]} #full and/or half siblings of the mother
      if (mgmid==0 & mgfid!=0) 
      {msibs <- ped[ped$FatherID==mgfid & ped$ID!=mid,"ID"]} #full and/or half siblings of the mother
      
      if (length(msibs)!=0) {
        for (jjj in 1:length(msibs)) {
          if(ped[ped$ID==msibs[jjj],"Gender"]==1){
            temp <- ped$FatherID==msibs[jjj]
            cid <- c(cid, ped[temp,"ID"])
          }else{
            temp <- ped$MotherID==msibs[jjj]
            cid <- c(cid, ped[temp,"ID"])
          }
        } 
      }
      
    }
    
    if (fid!=0) {
      pgmid <- ped[ped$ID==fid,"MotherID"]
      pgfid <- ped[ped$ID==fid,"FatherID"]
      if (pgmid!=0 | pgfid!=0) {psibs <- ped[(ped$FatherID==pgfid | ped$MotherID==pgmid) & ped$ID!=fid,"ID"]} #full and half siblings of the father
      if (pgmid!=0 & pgfid==0) 
      {psibs <- ped[ped$MotherID==pgmid & ped$ID!=fid,"ID"]} #full and/or half siblings of the mother
      if (pgmid==0 & pgfid!=0) 
      {psibs <- ped[ped$FatherID==pgfid & ped$ID!=fid,"ID"]} #full and/or half siblings of the mother
      
      if (length(psibs)!=0) {
        for (jjj in 1:length(psibs)) {
          if(ped[ped$ID==psibs[jjj],"Gender"]==1){
            temp <- ped$FatherID==psibs[jjj]
            cid <- c(cid, ped[temp,"ID"])
          }else{
            temp <- ped$MotherID==psibs[jjj]
            cid <- c(cid, ped[temp,"ID"])
          }
        } 
        
      } 
    }
    
    if (length(cid)==0) {cid <- 0}
    if (length(msibs)==0) {msibs <- 0}
    if (length(psibs)==0) {psibs <- 0}
    out <- list(cid, msibs, psibs)
    names(out) <- c("cousid", "msibs", "psibs")
    return(out)
    
  }
  
  ped <- fff[,c("ID", "Gender", "FatherID", "MotherID")]
  
  for (pp in 1:psize) {
    tempid <- findids(ped, ped$ID[pp])
    tempid.others <- getotherrels(ped, ped$ID[pp])
    tempid$cousid <- tempid.others$cousid
    tempid$msibs <- tempid.others$msibs
    tempid$psibs <- tempid.others$psibs
    
    # Check to make sure there are no married siblings
    if (length(intersect(tempid$sid, tempid$fsibid))>0) {
      if (intersect(tempid$sid, tempid$fsibid)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her sibling.  ", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
    # Check to make sure there are no married first cousins
    if (length(intersect(tempid$sid, tempid$cousid))>0) {
      if (intersect(tempid$sid, tempid$cousid)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her cousin.  ", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
    # Check to make sure paternal aunt/uncle is not married to maternal aunt/uncle
    temp <- table(tempid$cousid)
    if (sum(temp)/length(temp) > 1) {
      temperror <- paste("Error: ID ", ped$ID[pp], "'s paternal aunt/uncle is married to maternal aunt/uncle.  ", sep="")
      errormsg <- paste(errormsg, temperror)
    }
    
    # Check to make sure parent and child are not married
    if (length(intersect(tempid$sid, tempid$mid))>0) {
      if (intersect(tempid$sid, tempid$mid)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her mother.  ", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
    if (length(intersect(tempid$sid, tempid$fid))>0) {
      if (intersect(tempid$sid, tempid$fid)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her father.  ", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
    if (length(intersect(tempid$sid, tempid$msibs))>0) {
      if (intersect(tempid$sid, tempid$msibs)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her aunt/uncle.", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
    if (length(intersect(tempid$sid, tempid$psibs))>0) {
      if (intersect(tempid$sid, tempid$psibs)!=0) {
        temperror <- paste("Error: ID ", ped$ID[pp], " is married to his/her aunt/uncle.", sep="")
        errormsg <- paste(errormsg, temperror)
      }}
    
  }
  
  # Removed the collapsing of twins data, to be done later in each *pro.R
  # Do some twins checks here
  # If the twins column is not there, put it there as all zeroes
  if (is.null(fff$Twins)) {fff$Twins = 0}
  # Make sure twins are entered in multiples of 2
  twincheck <- sum(fff$Twins>0)
  if (twincheck %% 2 != 0) {errormsg <- paste(errormsg, "Error: pedigree contains an odd number of twins specified.  Because only identical twins are allowed in the calculation, they must be identified in pairs.")}
  
  if (sum(fff$Twins)>0 & (twincheck %% 2 == 0 | !any(table(fff$Twins[fff$Twins>0])==1))) {
    
    notwins <- fff[fff$Twins==0,]
    twins <- fff[fff$Twins>0,]
    newtwins <- NULL
    if (length(unique(fff$Twins))==1) {loopID <- 1}
    if (length(unique(fff$Twins))>1 & !any(fff$Twins==0)) {loopID <- length(unique(fff$Twins))}
    if (length(unique(fff$Twins))>1 & any(fff$Twins==0)) {loopID <- length(unique(fff$Twins))-1}
    for (iii in 1:loopID) {
      temp <- fff[fff$Twins==unique(fff$Twins[fff$Twins>0])[iii],]
      if (temp$Gender[1]!=temp$Gender[2]) {errormsg <- paste(errormsg, "Error: identical twins must be the same gender.")}
      if (temp$MotherID[1]!=temp$MotherID[2]) {errormsg <- paste(errormsg, "Error: identical twins must have the same mother.")}
      if (temp$FatherID[1]!=temp$FatherID[2]) {errormsg <- paste(errormsg, "Error: identical twins must have the same father.")}
      # If testing was done,
      if (!is.null(germline.testing) & model=='brcapro') {
        if (temp$BRCA1[1]==1 & temp$BRCA1[2]==2 | 
            temp$BRCA1[1]==2 & temp$BRCA1[2]==1 | 
            temp$BRCA2[1]==1 & temp$BRCA2[2]==2 | 
            temp$BRCA2[1]==2 & temp$BRCA2[2]==1) 
        {errormsg <- c(errormsg, 
                       "Error: identical twins should not have different germline testing results.")
        }
      }
      if (!is.null(germline.testing) & model=='MMRpro') {
        if (temp$MLH1[1]==1 & temp$MLH1[2]==2 | 
            temp$MLH1[1]==2 & temp$MLH1[2]==1 | 
            temp$MSH2[1]==1 & temp$MSH2[2]==2 | 
            temp$MSH2[1]==2 & temp$MSH2[2]==1 |
            temp$MSH6[1]==1 & temp$MSH6[2]==2 | 
            temp$MSH6[1]==2 & temp$MSH6[2]==1) 
        {errormsg <- c(errormsg, 
                       "Error: identical twins should not have different germline testing results.")
        }
      }
    }
  }
  
  
  
  #--------------------------------------------------------#
  # Use MakeRelationship to make Relation codes of members #
  
  fff$Relation = MakeRelationship(ped=fff[,c("ID", "Gender", "MotherID", "FatherID")], counselee.id=counselee.id)
  
  #----------------------------------------------------------------------------------------------------#
  # Check if pedigree has 1) only affected relatives, or 2) only affected relatives plus the 6 nuclear
  # relatives (parents, grandparents) who are unaffected, unknown age.
  # If so, and if impute=TRUE then run Lyte Simple # (new to version 2.1-2)
  
  if (model=="brcapro") {
    nunaff = sum(1*(fff$AffectedBreast==0 & fff$AffectedOvary==0 & fff$ID!=counselee.id))
    nunaff.extended = sum(1*(fff$AffectedBreast==0 & fff$AffectedOvary==0 & fff$ID!=counselee.id & 
                               fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
    nunaff.nuclear = sum(1*(fff$AffectedBreast==0 & fff$AffectedOvary==0 & fff$ID!=counselee.id & 
                              fff$AgeBreast==1 & fff$AgeOvary==1 & 
                              (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    n.aff.nuclear = sum(1*((fff$AffectedBreast>0 | fff$AffectedOvary==1) & fff$ID!=counselee.id & 
                             (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    if (nunaff==0) {runLyteSimple=TRUE} else 
      if (nunaff>0 & nunaff.extended==0 & (nunaff.nuclear + n.aff.nuclear == 6)) {runLyteSimple=TRUE} else{
        runLyteSimple=FALSE
      }
  }
  
  if (model=="MMRpro") {
    nunaff = sum(1*(fff$AffectedColon==0 & fff$AffectedEndometrium==0 & fff$AffectedGastric==0 & fff$ID!=counselee.id))
    nunaff.extended = sum(1*(fff$AffectedColon==0 & fff$AffectedEndometrium==0 & fff$AffectedGastric==0 & fff$ID!=counselee.id & 
                               fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
    nunaff.nuclear = sum(1*(fff$AffectedColon==0 & fff$AffectedEndometrium==0 & fff$AffectedGastric==0 & fff$ID!=counselee.id & 
                              fff$AgeColon==1 & fff$AgeEndometrium==1 & fff$AgeGastric==1 &
                              (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    n.aff.nuclear = sum(1*((fff$AffectedColon>0 | fff$AffectedEndometrium==1 & fff$AffectedGastric==1) & fff$ID!=counselee.id & 
                             (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    if (nunaff==0) {runLyteSimple=TRUE} else 
      if (nunaff>0 & nunaff.extended==0 & (nunaff.nuclear + n.aff.nuclear == 6)) {runLyteSimple=TRUE} else{
        runLyteSimple=FALSE
      }
  }
  if (model=="melapro") {
    nunaff = sum(1*(fff$AffectedSkin==0 & fff$ID!=counselee.id))
    nunaff.extended = sum(1*(fff$AffectedSkin==0 & fff$ID!=counselee.id & 
                               fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
    nunaff.nuclear = sum(1*(fff$AffectedSkin==0 & fff$ID!=counselee.id & fff$AgeSkin==1 & 
                              (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    n.aff.nuclear = sum(1*(fff$AffectedSkin>0 & fff$ID!=counselee.id & 
                             (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    if (nunaff==0) {runLyteSimple=TRUE} else 
      if (nunaff>0 & nunaff.extended==0 & (nunaff.nuclear + n.aff.nuclear == 6)) {runLyteSimple=TRUE} else{
        runLyteSimple=FALSE
      } 
  }
  if (model=="pancpro") {
    nunaff = sum(1*(fff$AffectedPancreas==0 & fff$ID!=counselee.id))
    nunaff.extended = sum(1*(fff$AffectedPancreas==0 & fff$ID!=counselee.id & 
                               fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
    nunaff.nuclear = sum(1*(fff$AffectedPancreas==0 & fff$ID!=counselee.id & fff$AgePancreas==1 & 
                              (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    n.aff.nuclear = sum(1*(fff$AffectedPancreas>0 & fff$ID!=counselee.id & 
                             (fff$Relation==4 | fff$Relation==5 | fff$Relation==7)))
    if (nunaff==0) {runLyteSimple=TRUE} else 
      if (nunaff>0 & nunaff.extended==0 & (nunaff.nuclear + n.aff.nuclear == 6)) {runLyteSimple=TRUE} else{
        runLyteSimple=FALSE
      }
  }
  
  
  
  if (model=="brcapro" & is.null(mastectomy)) {
    fff$Mastectomy = 0
    fff$AgeMastectomy=1
  }
  if (model=="brcapro" & is.null(oophorectomy)) {
    fff$Oophorectomy = 0
    fff$AgeOophorectomy = 1
  }
  if (model=="brcapro" & is.null(germline.testing)) {
    fff$BRCA1 = 0
    fff$BRCA2 = 0
    fff$TestOrder = 0
  }
  if (model=="brcapro" & is.null(marker.testing)) {
    fff$ER = 0
    fff$PR = 0
    fff$CK14 = 0
    fff$CK5.6 = 0
    fff$HER2 = 0
  }
  if (model=="MMRpro" & is.null(germline.testing)) {
    fff$MLH1 = 0
    fff$MSH2 = 0
    fff$TestOrder = 0
    fff$MSH6 = 0
  }
  if (model=="MMRpro" & is.null(marker.testing)) {
    fff$MSI = 0
    fff$location = 0
  }
  if (model=="melapro" & is.null(germline.testing)) {
    fff$P16 = 0
    fff$TestOrder = 0
  }
  
  if(model == "MMRpro"){
    if (runLyteSimple & imputeRelatives==TRUE) { #run Lyte simple
      fff = LyteSimple.gast(fff, params=params, model=model)  
    }
  } else{
    if (runLyteSimple & imputeRelatives==TRUE) { #run Lyte simple
      fff = LyteSimple(fff, params=params, model=model)  
    }
  }
  
  # Check for errors in family, only useful when lyte simple is not run.
  if (!runLyteSimple | imputeRelatives==FALSE) {
    famcheck <- kinship2::familycheck(famid=rep(1,nrow(fff)), id=fff$ID, father.id=fff$FatherID, mother.id=fff$MotherID)
    if (famcheck$split>1) {warnmsg <- "Warning: This family appears to be a set of disjoint trees.  Some member(s) are not linked directly to the proband through FatherID and/or MotherID and have been dropped from the calculation.  Please re-check your pedigree; Are you missing some of the people?"}
  }
  
  
  #--------------#
  # Final output #
  
  if (nchar(errormsg) > 0) {return(errormsg)}
  if (nchar(warnmsg) > 0) {print(warnmsg)}
  
  return(fff)
  
}
