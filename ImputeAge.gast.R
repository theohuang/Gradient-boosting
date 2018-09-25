### This is for version 2.1 ##
### Takes in a single fff in BM pedigree format and imputes ages
### of fff members with age=1
### Returns the fff in the same format
### A family with any age set to 1 will come through this function

## Main code written by Philamer Atienza, structured for BM by Amanda B.

## version 2.1-5
## using crude penetrance by default

ImputeAge.gast = function(fff, params, model, net = FALSE) {
  
  warnmsg = NULL
  
  age.max <- 94
  
  if(net == TRUE){
    params$penetrance <- params$penetrance.net
    params$CBCpenetrance <- cbc.net
  } else{
    params$penetrance <- params$penetrance.crude
    params$CBCpenetrance <- cbc.crude
  }
  
  # Get constants for ages of relatives
  col.med.1.sis=params$col.med.1.sis
  col.med.2.sis=params$col.med.2.sis
  col.med.5.sis=params$col.med.5.sis
  col.med.1.bro=params$col.med.1.bro
  col.med.2.bro=params$col.med.2.bro
  col.med.5.bro=params$col.med.5.bro
  col.med.1.paunt=params$col.med.1.paunt
  col.med.2.paunt=params$col.med.2.paunt
  col.med.5.paunt=params$col.med.5.paunt
  col.med.1.puncl=params$col.med.1.puncl
  col.med.2.puncl=params$col.med.2.puncl
  col.med.5.puncl=params$col.med.5.puncl
  col.med.1.maunt=params$col.med.1.maunt
  col.med.2.maunt=params$col.med.2.maunt
  col.med.5.maunt=params$col.med.5.maunt
  col.med.1.muncl=params$col.med.1.muncl
  col.med.2.muncl=params$col.med.2.muncl
  col.med.5.muncl=params$col.med.5.muncl
  col.med.1.mom=params$col.med.1.mom
  col.med.1.dad=params$col.med.1.dad
  col.med.1.pgrma=params$col.med.1.pgrma
  col.med.1.pgrpa=params$col.med.1.pgrpa
  col.med.1.mgrma=params$col.med.1.mgrma
  col.med.1.mgrpa=params$col.med.1.mgrpa
  col.med.1.daught=params$col.med.1.daught
  col.med.2.daught=params$col.med.2.daught
  col.med.5.daught=params$col.med.5.daught
  col.med.1.son=params$col.med.1.son
  col.med.2.son=params$col.med.2.son
  col.med.5.son=params$col.med.5.son
  col.med.1.niecenephew=params$col.med.1.niecenephew
  col.med.2.niecenephew=params$col.med.2.niecenephew
  col.med.5.niecenephew=params$col.med.5.niecenephew
  
  # Get needed constants  
  if (model=="brcapro"){cnslage <- max(fff[which(fff$Relation==1),c("AgeBreast", "AgeOvary", "AgeBreastContralateral")])}
  if (model=="MMRpro"){cnslage <- max(fff[which(fff$Relation==1),c("AgeColon", "AgeEndometrium", "AgeGastric")])}
  if (model=="melapro"){cnslage <- fff[which(fff$Relation==1),c("AgeSkin")]}
  if (model=="pancpro"){cnslage <- fff[which(fff$Relation==1),c("AgePancreas")]}
  
  if (model=="brcapro") {ca1 = "AffectedBreast"; ca2 = "AffectedOvary"; age1 = "AgeBreast"; age2 = "AgeOvary"; nc="B00"}
  if (model=="MMRpro") {ca1 = "AffectedColon"; ca2 = "AffectedEndometrium"; ca3 = "AffectedGastric"; age1 = "AgeColon"; age2 = "AgeEndometrium"; age3 = "AgeGastric"; nc="M000"}
  if (model=="melapro") {ca1 = "AffectedSkin"; age1 = "AgeSkin"; nc = "P160"}
  if (model=="pancpro") {ca1 = "AffectedPancreas"; age1 = "AgePancreas"; nc = "P0"}
  
  cnsl    = fff$ID[fff$Relation == 1]
  cnslbc  = fff[which(fff$ID == cnsl), ca1] > 0
  if (model=="brcapro") {cnsloc  = fff[which(fff$ID == cnsl),ca2]}
  if (model=="melapro" | model=="pancpro") {cnsloc = FALSE}
  if(model == "MMRpro"){
    cnslec <- fff[which(fff$ID == cnsl), ca2]
    cnslgc <- fff[which(fff$ID == cnsl), ca3]
  }
  nhus    = sum(as.logical(fff$Relation == 14)) 			# Husband
  nsis    = sum(as.logical((fff$Relation == 2|fff$Relation == 16) & fff$Gender == 0))	# Sister
  nbro    = sum(as.logical((fff$Relation == 2|fff$Relation == 16) & fff$Gender == 1))	# Brother
  ndaught = sum(as.logical(fff$Relation == 3 & fff$Gender == 0))	# Daughter
  nson    = sum(as.logical(fff$Relation == 3 & fff$Gender == 1))	# Son
  nmom    = sum(as.logical(fff$Relation == 4 & fff$Gender == 0))	# Mother
  ndad    = sum(as.logical(fff$Relation == 4 & fff$Gender == 1))	# Father
  npgrma  = sum(as.logical(fff$Relation == 5 & fff$Gender == 0))	# Paternal Grandmother
  npgrpa  = sum(as.logical(fff$Relation == 5 & fff$Gender == 1))	# Paternal Grandfather
  npaunt  = sum(as.logical(fff$Relation == 6 & fff$Gender == 0))	# Paternal Aunt
  npuncl  = sum(as.logical(fff$Relation == 6 & fff$Gender == 1))	# Paternal Uncle
  nmgrma  = sum(as.logical(fff$Relation == 7 & fff$Gender == 0))	# Maternal Grandmother
  nmgrpa  = sum(as.logical(fff$Relation == 7 & fff$Gender == 1))	# Maternal Grandfather
  nmaunt  = sum(as.logical(fff$Relation == 8 & fff$Gender == 0))	# Maternal Aunt
  nmuncl  = sum(as.logical(fff$Relation == 8 & fff$Gender == 1))	# Maternal Uncle
  nniecenephew = sum(as.logical(fff$Relation == 13)) # Niece or nephew
  
  if(model == "MMRpro"){
    summaryfff = cbind(cnsl,cnslage,cnslbc,cnslec, cnslgc, nhus,nsis,nbro,ndaught,nson,nmom,ndad,npgrma,npgrpa,npaunt,npuncl,nmgrma,nmgrpa,nmaunt,nmuncl,nniecenephew)
  } else{
    summaryfff = cbind(cnsl,cnslage,cnslbc,cnsloc,nhus,nsis,nbro,ndaught,nson,nmom,ndad,npgrma,npgrpa,npaunt,npuncl,nmgrma,nmgrpa,nmaunt,nmuncl,nniecenephew)
  }
  
  combi = merge(fff,summaryfff)
  
  ### Impute Age for Unaffected members from Median in Colorectal Data. Proband age untouched. 
  ### ___________________
  
  # Get the indeces of the family members with unknown age and unaffected
  # Will apply Lyte+ estimates to these individuals only
  # uua=unaffected,unknown age
  
  if (model=='brcapro') {combi$uua = combi[,ca1]==0 & combi[,ca2]==0 & fff[,age1]==1 & fff[,age2]==1}
  if (model=="melapro" | model=="pancpro") {combi$uua = combi[,ca1]==0 & combi[,age1]==1}
  if(model == "MMRpro"){
    combi$uua <- combi[, ca1] == 0 & combi[, ca2] == 0 & combi[, ca3] == 0 &
      fff[, age1] == 1 & fff[, age2] == 1 & fff[, age3] == 1
  }
  
  combi$nsis <- as.numeric(as.character(combi$nsis))
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis==1]=col.med.1.sis
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=2 & combi$nsis<=4]=col.med.2.sis
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=5]=col.med.5.sis
  
  combi$nbro <- as.numeric(as.character(combi$nbro))
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro==1]=col.med.1.bro
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=2 & combi$nbro<=4]=col.med.2.bro
  combi[,age1][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=5]=col.med.5.bro
  
  combi$npaunt <- as.numeric(as.character(combi$npaunt))
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt==1]=col.med.1.paunt
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=2 & combi$npaunt<=4]=col.med.2.paunt
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=5]=col.med.5.paunt
  
  combi$npuncl <- as.numeric(as.character(combi$npuncl))
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl==1]=col.med.1.puncl
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=2 & combi$npuncl<=4]=col.med.2.puncl
  combi[,age1][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=5]=col.med.5.puncl
  
  combi$nmaunt <- as.numeric(as.character(combi$nmaunt))
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt==1]=col.med.1.maunt
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=2 & combi$nmaunt<=4]=col.med.2.maunt
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=5]=col.med.5.maunt
  
  combi$nmuncl <- as.numeric(as.character(combi$nmuncl))
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl==1]=col.med.1.muncl
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=2 & combi$nmuncl<=4]=col.med.2.muncl
  combi[,age1][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=5]=col.med.5.muncl
  
  combi$nmom <- as.numeric(as.character(combi$nmom))
  combi[,age1][combi$uua==1 & combi$Relation==4 & combi$Gender==0 & combi$nmom==1]=col.med.1.mom
  
  combi$ndad <- as.numeric(as.character(combi$ndad))
  combi[,age1][combi$uua==1 & combi$Relation==4 & combi$Gender==1 & combi$ndad==1]=col.med.1.dad
  
  combi$npgrma <- as.numeric(as.character(combi$npgrma))
  combi[,age1][combi$uua==1 & combi$Relation==5 & combi$Gender==0 & combi$npgrma==1]=col.med.1.pgrma
  
  combi$npgrpa <- as.numeric(as.character(combi$npgrpa))
  combi[,age1][combi$uua==1 & combi$Relation==5 & combi$Gender==1 & combi$npgrpa==1]=col.med.1.pgrpa
  
  combi$nmgrma <- as.numeric(as.character(combi$nmgrma))
  combi[,age1][combi$uua==1 & combi$Relation==7 & combi$Gender==0 & combi$nmgrma==1]=col.med.1.mgrma
  
  combi$nmgrpa <- as.numeric(as.character(combi$nmgrpa))
  combi[,age1][combi$uua==1 & combi$Relation==7 & combi$Gender==1 & combi$nmgrpa==1]=col.med.1.mgrpa
  
  combi$ndaught <- as.numeric(as.character(combi$ndaught))
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught==1]=col.med.1.daught
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=2 & combi$ndaught<=4]=col.med.2.daught
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=5]=col.med.5.daught
  
  combi$nson <- as.numeric(as.character(combi$nson))
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson==1]=col.med.1.son
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=2 & combi$nson<=4]=col.med.2.son
  combi[,age1][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=5]=col.med.5.son
  
  combi$nniecenephew <- as.numeric(as.character(combi$nniecenephew))
  combi[,age1][combi$uua==1 & combi$Relation==13 & combi$nniecenephew==1]=col.med.1.niecenephew
  combi[,age1][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=2 & combi$nniecenephew<=4]=col.med.2.niecenephew
  combi[,age1][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=5]=col.med.5.niecenephew
  
  # If relatives of 3rd degree or higher or brother/sister in-laws and spouses are 
  # unaffected unknown age, use age of proband
  combi[,age1][combi$uua==1 & (combi$Relation==0 | combi$Relation==14 | combi$Relation==15)] = cnslage
  
  if (model=='brcapro' | model=='MMRpro') {
    
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis==1]=col.med.1.sis
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=2 & combi$nsis<=4]=col.med.2.sis
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=5]=col.med.5.sis
    
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro==1]=col.med.1.bro
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=2 & combi$nbro<=4]=col.med.2.bro
    combi[,age2][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=5]=col.med.5.bro
    
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt==1]=col.med.1.paunt
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=2 & combi$npaunt<=4]=col.med.2.paunt
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=5]=col.med.5.paunt
    
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl==1]=col.med.1.puncl
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=2 & combi$npuncl<=4]=col.med.2.puncl
    combi[,age2][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=5]=col.med.5.puncl
    
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt==1]=col.med.1.maunt
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=2 & combi$nmaunt<=4]=col.med.2.maunt
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=5]=col.med.5.maunt
    
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl==1]=col.med.1.muncl
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=2 & combi$nmuncl<=4]=col.med.2.muncl
    combi[,age2][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=5]=col.med.5.muncl
    
    combi[,age2][combi$uua==1 & combi$Relation==4 & combi$Gender==0 & combi$nmom==1]=col.med.1.mom
    
    combi[,age2][combi$uua==1 & combi$Relation==4 & combi$Gender==1 & combi$ndad==1]=col.med.1.dad
    combi[,age2][combi$uua==1 & combi$Relation==5 & combi$Gender==0 & combi$npgrma==1]=col.med.1.pgrma
    combi[,age2][combi$uua==1 & combi$Relation==5 & combi$Gender==1 & combi$npgrpa==1]=col.med.1.pgrpa
    combi[,age2][combi$uua==1 & combi$Relation==7 & combi$Gender==0 & combi$nmgrma==1]=col.med.1.mgrma
    combi[,age2][combi$uua==1 & combi$Relation==7 & combi$Gender==1 & combi$nmgrpa==1]=col.med.1.mgrpa
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught==1]=col.med.1.daught
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=2 & combi$ndaught<=4]=col.med.2.daught
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=5]=col.med.5.daught
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson==1]=col.med.1.son
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=2 & combi$nson<=4]=col.med.2.son
    combi[,age2][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=5]=col.med.5.son
    combi[,age2][combi$uua==1 & combi$Relation==13 & combi$nniecenephew==1]=col.med.1.niecenephew
    combi[,age2][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=2 & combi$nniecenephew<=4]=col.med.2.niecenephew
    combi[,age2][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=5]=col.med.5.niecenephew
    combi[,age2][combi$uua==1 & (combi$Relation==0 | combi$Relation==14 | combi$Relation==15)] = cnslage
    
  }
  
  if(model == "MMRpro"){
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis==1]=col.med.1.sis
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=2 & combi$nsis<=4]=col.med.2.sis
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==0 & combi$nsis>=5]=col.med.5.sis
    
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro==1]=col.med.1.bro
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=2 & combi$nbro<=4]=col.med.2.bro
    combi[,age3][combi$uua==1 & (combi$Relation==2|combi$Relation==16) & combi$Gender==1 & combi$nbro>=5]=col.med.5.bro
    
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt==1]=col.med.1.paunt
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=2 & combi$npaunt<=4]=col.med.2.paunt
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==0 & combi$npaunt>=5]=col.med.5.paunt
    
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl==1]=col.med.1.puncl
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=2 & combi$npuncl<=4]=col.med.2.puncl
    combi[,age3][combi$uua==1 & combi$Relation==6 & combi$Gender==1 & combi$npuncl>=5]=col.med.5.puncl
    
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt==1]=col.med.1.maunt
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=2 & combi$nmaunt<=4]=col.med.2.maunt
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==0 & combi$nmaunt>=5]=col.med.5.maunt
    
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl==1]=col.med.1.muncl
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=2 & combi$nmuncl<=4]=col.med.2.muncl
    combi[,age3][combi$uua==1 & combi$Relation==8 & combi$Gender==1 & combi$nmuncl>=5]=col.med.5.muncl
    
    combi[,age3][combi$uua==1 & combi$Relation==4 & combi$Gender==0 & combi$nmom==1]=col.med.1.mom
    
    combi[,age3][combi$uua==1 & combi$Relation==4 & combi$Gender==1 & combi$ndad==1]=col.med.1.dad
    combi[,age3][combi$uua==1 & combi$Relation==5 & combi$Gender==0 & combi$npgrma==1]=col.med.1.pgrma
    combi[,age3][combi$uua==1 & combi$Relation==5 & combi$Gender==1 & combi$npgrpa==1]=col.med.1.pgrpa
    combi[,age3][combi$uua==1 & combi$Relation==7 & combi$Gender==0 & combi$nmgrma==1]=col.med.1.mgrma
    combi[,age3][combi$uua==1 & combi$Relation==7 & combi$Gender==1 & combi$nmgrpa==1]=col.med.1.mgrpa
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught==1]=col.med.1.daught
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=2 & combi$ndaught<=4]=col.med.2.daught
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==0 & combi$ndaught>=5]=col.med.5.daught
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson==1]=col.med.1.son
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=2 & combi$nson<=4]=col.med.2.son
    combi[,age3][combi$uua==1 & combi$Relation==3 & combi$Gender==1 & combi$nson>=5]=col.med.5.son
    combi[,age3][combi$uua==1 & combi$Relation==13 & combi$nniecenephew==1]=col.med.1.niecenephew
    combi[,age3][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=2 & combi$nniecenephew<=4]=col.med.2.niecenephew
    combi[,age3][combi$uua==1 & combi$Relation==13 & combi$nniecenephew>=5]=col.med.5.niecenephew
    combi[,age3][combi$uua==1 & (combi$Relation==0 | combi$Relation==14 | combi$Relation==15)] = cnslage
  }
  
  ### Get only needed columns
  ### _____________________________________________________________________________________
  
  if(model == "MMRpro"){
    nonfamilyvars = c("cnsl","cnslage","cnslbc","cnslec", "cnslgc","nhus","nsis","nbro","ndaught","nson","nmom","ndad","npgrma","npgrpa","npaunt","npuncl","nmgrma","nmgrpa","nmaunt","nmuncl","nniecenephew", "uua")
  } else{
    nonfamilyvars = c("cnsl","cnslage","cnslbc","cnsloc","nhus","nsis","nbro","ndaught","nson","nmom","ndad","npgrma","npgrpa","npaunt","npuncl","nmgrma","nmgrpa","nmaunt","nmuncl","nniecenephew", "uua")
  }
  fff=combi[,-which(is.element(colnames(combi),nonfamilyvars))]
  
  
  if (model=='brcapro') {
    # If mastectomy age is missing, replace with age breast or age contralateral
    # If both are missing, impute below
    if (any((fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2)) {warnmsg = c(warnmsg, "Warning: Unknown mastectomy age for those with mastectomy has been replaced with current age, via AgeBreast, which also may have been imputed, see other warnings.")}
    fff$AgeMastectomy[(fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2] = fff$AgeBreast[(fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AgeBreast!=1 & fff$AffectedBreast!=2]
    
    if (any((fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1)) {warnmsg = c(warnmsg, "Warning: Unknown mastectomy age for those with mastectomy has been replaced with age of contralateral breast cancer, via AgeBreastContralateral, which also may have been imputed, see other warnings.")}
    fff$AgeMastectomy[(fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1] = fff$AgeBreastContralateral[(fff$AgeMastectomy==1|fff$AgeMastectomy==0) & fff$Mastectomy==1 & fff$AffectedBreast==2 & fff$AgeBreastContralateral!=1]
    
    
    # If oophorectomy age is missing, replace with age ovary
    # If both are missing, impute below
    if (any((fff$AgeOophorectomy==1 | fff$AgeOophorectomy==0) & fff$Oophorectomy==1 & fff$AgeOvary!=1)) {warnmsg = c(warnmsg, "Warning: Unknown oophorectomy age for those with oophorectomy has been replaced with current age, via AgeOvary, which also may have been imputed, see other warnings.")}
    fff$AgeOophorectomy[(fff$AgeOophorectomy==1 | fff$AgeOophorectomy==0) & fff$Oophorectomy==1 & fff$AgeOvary!=1] = fff$AgeOvary[(fff$AgeOophorectomy==1 | fff$AgeOophorectomy==0) & fff$Oophorectomy==1 & fff$AgeOvary!=1]
  }
  
  ### --------------------------------------------------------------
  ### Now apply Danielle's multiple imputation approach for imputing
  ### ages of affected relatives when affection age is input as 1 (unknown)
  
  #number of samples for multiple imputation
  nIter=params$nIter
  
  # if a user specifies a set seed, put it here.
  if (!is.null(params$myseed)) {set.seed(params$myseed)}
  
  #Relatives 
  fm<-dim(fff)[1]
  
  #creating variables to store imputed ages, use bc=crc=pc=skinca; oc=ec
  age_bc<-matrix(0,nIter,fm)
  age_bcc<-matrix(0,nIter,fm)
  age_oc<-matrix(0,nIter,fm)
  age_gc<-matrix(0,nIter,fm)
  age_mast = matrix(0,nIter,fm)
  age_ooph = matrix(0,nIter,fm)
  
  #probabilities of cancer by age, should be a vector from 1:age.max. 
  seer_bc_females<-params$penetrance$fFX[,nc]
  seer_bc_males<-params$penetrance$fMX[,nc]
  seer_oc_females<-params$penetrance$fFY[,nc]
  seer_gc_females<-params$penetrance$fFZ[,nc]
  seer_gc_males<-params$penetrance$fMZ[,nc]
  
  if (model=='brcapro') {
    seer_bcc_under40_females<-params$CBCpenetrance$fFX.Under40[,"B00"]
    seer_bcc_under40_males<-params$CBCpenetrance$fMX.Under40[,"B00"]
    seer_bcc_over40_females<-params$CBCpenetrance$fFX.Over40[,"B00"]
    seer_bcc_over40_males<-params$CBCpenetrance$fMX.Over40[,"B00"]
  }
  
  #standardizing these probabilities (in case they aren't)
  pb_bc_females<-seer_bc_females/sum(seer_bc_females, na.rm = TRUE)
  pb_oc_females<-seer_oc_females/sum(seer_oc_females, na.rm = TRUE)
  pb_bc_males<-seer_bc_males/sum(seer_bc_males, na.rm = TRUE)
  pb_gc_females<-seer_gc_females/sum(seer_gc_females, na.rm = TRUE)
  pb_gc_males<-seer_gc_males/sum(seer_gc_males, na.rm = TRUE)
  
  
  if (model=='brcapro') {
    pb_bcc_under40_females<-seer_bcc_under40_females/sum(seer_bcc_under40_females, na.rm = TRUE)
    pb_bcc_over40_females<-seer_bcc_over40_females/sum(seer_bcc_over40_females, na.rm = TRUE)
    pb_bcc_under40_males<-seer_bcc_under40_males/sum(seer_bcc_under40_males, na.rm = TRUE)
    pb_bcc_over40_males<-seer_bcc_over40_males/sum(seer_bcc_over40_males, na.rm = TRUE)
  }
  
  ####################################################################################################
  # Scenario 1: BC (brcapro), CRC (mmrpro), PC (pancpro) or SkinCa (melapro)
  
  #determining which fff members you would sample, assume age=1 means unknown age. 
  if (model=='brcapro' | model=='MMRpro') {
    mem_bc = which(fff[,ca1]>=1 & fff[,age1]==1 & !(fff[, ca2] == 1 & fff[, age2] == 1))
    mem_bc_females<-which(fff[,ca1]>=1 & fff[,age1]==1 & !(fff[, ca2] == 1 & fff[, age2] == 1) & fff$Gender==0)
    mem_bc_males<-which(fff[,ca1]>=1 & fff[,age1]==1 & !(fff[, ca2] == 1 & fff[, age2] == 1) & fff$Gender==1)
  }
  
  if (model=='pancpro' | model=='melapro') {
    mem_bc = which(fff[,ca1]>=1 & fff[,age1]==1)
    mem_bc_females<-which(fff[,ca1]>=1 & fff[,age1]==1 & fff$Gender==0)
    mem_bc_males<-which(fff[,ca1]>=1 & fff[,age1]==1  & fff$Gender==1)
  }
  
  #sampling n potential BC ages for each individual with missing BC age. 
  if (length(mem_bc_females)>0){
    age_bc[,mem_bc_females] = sample(age.max,nIter*length(mem_bc_females),replace=TRUE,prob=pb_bc_females)
    # age_oc[,mem_bc_females] = age_bc[,mem_bc_females]
  }
  if (length(mem_bc_males)>0){
    age_bc[,mem_bc_males] = sample(age.max,nIter*length(mem_bc_males),replace=TRUE,prob=pb_bc_males)
    # age_oc[,mem_bc_males] = age_bc[,mem_bc_males]
  }
  
  #################################################################
  # Scenario 2: OC (brcapro) or EC (mmrpro) #
  
  if (model=='brcapro' | model=='MMRpro') {
    #determining which fff members you would sample, assume age=1 means unknown age. 
    mem_oc<-which(fff[,ca2]==1 & fff[,age2]==1 & !(fff[, ca1] >= 1 & fff[, age1] == 1))
    #sampling n potential OC ages for each individual with missing OC age. 
    if (length(mem_oc)>0){
      age_oc[,mem_oc]<-sample(age.max,nIter*length(mem_oc),replace=TRUE, prob=pb_oc_females)
    }
  }
  
  #################################################################
  # Scenario 3: GC (mmrpro) #
  
  if (model=='MMRpro') {
    mem_gc = which(fff[,ca3]>=1 & fff[,age3]==1)
    mem_gc_females<-which(fff[,ca3]>=1 & fff[,age3]==1 & fff$Gender==0)
    mem_gc_males<-which(fff[,ca3]>=1 & fff[,age3]==1 & fff$Gender==1)
    
    #sampling n potential GC ages for each individual with missing GC age. 
    if (length(mem_gc_females)>0){
      age_gc[,mem_gc_females] = sample(age.max,nIter*length(mem_gc_females),replace=TRUE,prob=pb_gc_females)
      # age_oc[,mem_gc_females] = age_gc[,mem_gc_females]
    }
    if (length(mem_gc_males)>0){
      age_gc[,mem_gc_males] = sample(age.max,nIter*length(mem_gc_males),replace=TRUE,prob=pb_gc_males)
      # age_oc[,mem_gc_males] = age_gc[,mem_gc_males]
    }
  }
  
  ##########################################################
  # Scenario 4: BC and OC (brcapro) or CRC and EC (mmrpro) #  
  if (model=='brcapro') {jointdsn = params$BrOvJointDsn}
  if (model=='MMRpro') {jointdsn = params$CrcEcJointDsn}
  #sampling n potential OC and BC ages simultaneously for each individual with missing BC and OC ages.
  if (model=='brcapro' | model=='MMRpro') {
    mem_bc_oc = which(fff[,ca2]==1 & fff[,age2]==1 & fff[,ca1]>0 & fff[,age1]==1)
    if (length(mem_bc_oc)>0) {
      mysample = sample(10000,nIter*length(mem_bc_oc), replace=TRUE, prob=jointdsn$Probability)
      age_bc[,mem_bc_oc] = jointdsn[,age1][mysample]
      age_oc[,mem_bc_oc] = jointdsn[,age2][mysample]
    }
  }
  
  if (model=='brcapro') {
    # If age of mastectomy is unknown, replace with imputed BC age. 
    mem_mast = which(fff$AgeMastectomy==1 & fff$Mastectomy==1 & fff[,age1]==1 & fff[,ca1]==1)
    if (length(mem_mast)>0) {
      age_mast[,mem_mast] = age_bc[,mem_mast]
    }
    
    # If age of oophorectomy is unknown, replace with imputed OC age
    mem_ooph = which(fff$Oophorectomy==1 & fff$AgeOophorectomy==1 & fff[,age2]==1 & fff[,ca2]==1)
    if (length(mem_ooph)>0) {
      age_ooph[,mem_ooph] = age_oc[,mem_ooph]
    }
  }
  
  ############################
  # Scenario 5: Bilateral BC #
  
  if (model=='brcapro') {
    #determining which fff members you would sample, assume age=1 means unknown age.
    mem_bcc = which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1)
    mem_bcc_females<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$Gender==0)
    mem_bcc_males<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$Gender==1)
    
    #for those who also have known age for first diangosis age.  
    mem_bcc1_under40_females<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast!=1 & fff$AgeBreast<40 & fff$Gender==0)
    mem_bcc1_over40_females<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast!=1 & fff$AgeBreast>=40 & fff$Gender==0)
    mem_bcc1_under40_males<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast!=1 & fff$AgeBreast<40 & fff$Gender==1)
    mem_bcc1_over40_males<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast!=1 & fff$AgeBreast>=40 & fff$Gender==1)
    
    #sampling n potential BCC ages for each individual with missing BCC age. 
    if (length(mem_bcc1_under40_females)>0){
      for (mmm in mem_bcc1_under40_females) {
        agefirst<-fff$AgeBreast[mmm]
        pb_bcc<-seer_bcc_under40_females[1:(age.max - agefirst + 1)]/sum(seer_bcc_under40_females[1:(age.max - agefirst + 1)], na.rm = TRUE)
        numyears<-agefirst:age.max
        if (agefirst==age.max) {
          age_bcc[,mmm] = age.max
        }
        if (agefirst<age.max) {
          age_bcc[,mmm]<-sample(numyears,nIter,replace=TRUE, prob=pb_bcc)
        }
      }
    }
    if (length(mem_bcc1_over40_females)>0){
      for (mmm in mem_bcc1_over40_females) {
        agefirst<-fff$AgeBreast[mmm]
        pb_bcc<-seer_bcc_over40_females[1:(age.max - agefirst + 1)]/sum(seer_bcc_over40_females[1:(age.max - agefirst + 1)])
        numyears<-agefirst:age.max
        if (agefirst==age.max) {
          age_bcc[,mmm] = age.max
        }
        if (agefirst<age.max) {
          age_bcc[,mmm]<-sample(numyears,nIter,replace=TRUE, prob=pb_bcc)
        }
      }
    }
    if (length(mem_bcc1_under40_males)>0){
      for (mmm in mem_bcc1_under40_males) {
        agefirst<-fff$AgeBreast[mmm]
        pb_bcc<-seer_bcc_under40_males[1:(age.max - agefirst + 1)]/sum(seer_bcc_under40_males[1:(age.max - agefirst + 1)])
        numyears<-agefirst:age.max
        if (agefirst==age.max) {
          age_bcc[,mmm] = age.max
        }
        if (agefirst<age.max) {
          age_bcc[,mmm]<-sample(numyears,nIter,replace=TRUE, prob=pb_bcc)
        }
      }
    }
    if (length(mem_bcc1_over40_males)>0){
      for (mmm in mem_bcc1_over40_males) {
        agefirst<-fff$AgeBreast[mmm]
        pb_bcc<-seer_bcc_over40_males[1:(age.max - agefirst + 1)]/sum(seer_bcc_over40_males[1:(age.max - agefirst + 1)])
        numyears<-agefirst:age.max
        if (agefirst==age.max) {
          age_bcc[,mmm] = age.max
        }
        if (agefirst<age.max) {
          age_bcc[,mmm]<-sample(numyears,nIter,replace=TRUE, prob=pb_bcc)
        }
      }
    }
    
    #for those who also have missing first diangosis age.  
    # Use their imputed age of 1st BC to decide which dsn to draw from.
    mem_bcc2<-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast==1)
    mem_bcc2_females <-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast==1 & fff$Gender==0)
    mem_bcc2_males <-which(fff$AffectedBreast==2 & fff$AgeBreastContralateral==1 & fff$AgeBreast==1 & fff$Gender==1)
    
    #sampling n potential BCC ages for each individual with missing BCC age. 
    if (length(mem_bcc2_females)>0){
      for (i in mem_bcc2_females) {
        for (j in 1:nIter){
          agefirst<-age_bc[j,i]
          if (agefirst <= 40){
            pb_bcc<-seer_bcc_under40_females[1:(age.max - agefirst + 1)]/sum(seer_bcc_under40_females[1:(age.max - agefirst + 1)])
          }
          if (agefirst > 40){
            pb_bcc<-seer_bcc_over40_females[1:(age.max - agefirst + 1)]/sum(seer_bcc_over40_females[1:(age.max - agefirst + 1)])
          }
          numyears<-agefirst:age.max
          if (length(numyears)==1) {
            age_bcc[j,i] = age.max
          }
          if (length(numyears)>1) {
            age_bcc[j,i]<-sample(numyears,1,prob=pb_bcc)
          }
        }
      }
    }
    
    if (length(mem_bcc2_males)>0){
      for (i in mem_bcc2_males) {
        for (j in 1:nIter){
          agefirst<-age_bc[j,i]
          if (agefirst <= 40){
            pb_bcc<-seer_bcc_under40_males[1:(age.max - agefirst + 1)]/sum(seer_bcc_under40_males[1:(age.max - agefirst + 1)])
          }
          if (agefirst > 40){
            pb_bcc<-seer_bcc_over40_males[1:(age.max - agefirst + 1)]/sum(seer_bcc_over40_males[1:(age.max - agefirst + 1)])
          }
          numyears<-agefirst:age.max
          if (length(numyears)==1) {
            age_bcc[j,i] = age.max
          }
          if (length(numyears)>1) {
            age_bcc[j,i]<-sample(numyears,1,prob=pb_bcc)
          }
        }
      }
    }
    
    # If age of mastectomy is unknown and bilateral BC, replace with imputed BCC age. 
    mem_mast_bcc = which(fff$AgeMastectomy==1 & fff$Mastectomy==1 & fff$AgeBreast==1 & fff$AffectedBreast==2)
    if (length(mem_mast_bcc)>0) {
      age_mast[,mem_mast_bcc] = age_bcc[,mem_mast_bcc]
    }
    
  } #end of if model=brcapro
  
  # this part of the function would stop now and return
  # fff with imputed ages on unaffected
  # age_bc, age_oc and age_bcc to be used in a loop in *pro.R
  
  if (any(nchar(warnmsg)>0)) {print(warnmsg)}
  
  if (model=='brcapro'){
    ifamily = list(fff, mem_bc, mem_oc, mem_bcc, mem_mast, mem_ooph, age_bc, age_oc, age_bcc, age_mast, age_ooph)
    names(ifamily) = c('fff', 'mem_bc', 'mem_oc', 'mem_bcc', 'mem_mast', 'mem_ooph', 'age_bc', 'age_oc', 'age_bcc', 'age_mast', 'age_ooph')
  }
  if (model=='MMRpro'){
    ifamily = list(fff, mem_bc, mem_oc, mem_gc, age_bc, age_oc, age_gc)
    names(ifamily) = c('fff', 'mem_bc', 'mem_oc', "mem_gc",  'age_bc', 'age_oc', "age_gc")
  }
  if (model=='pancpro' | model=='melapro') {
    ifamily = list(fff, mem_bc, age_bc)
    names(ifamily) = c('fff', 'mem_bc', 'age_bc')    
  }
  
  
  
  return(ifamily)
  
}
















