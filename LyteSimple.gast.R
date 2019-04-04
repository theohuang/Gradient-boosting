## Returns the fff in the same format, with relatives added
## All families will come through this function

## Main code written by Philamer Atienza, structured for BayesMendel by Amanda

LyteSimple.gast = function(fff, params, model) {

### Load Colorectal Data
### _____________________________________________________________________________________
col.med.sis=1
col.med.bro=1
col.med.paunt=1
col.med.puncl=1
col.med.maunt=1
col.med.muncl=1
col.med.mom=1
col.med.dad=1
col.med.pgrma=1
col.med.pgrpa=1
col.med.mgrma=1
col.med.mgrpa=1
col.med.daught=1
col.med.son=1
col.med.niecenephew=2

warnmsg = NULL

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

if (model=="brcapro") {
  ca1 = "AffectedBreast"
  ca2 = "AffectedOvary"
  age1 = "AgeBreast"
  age2 = "AgeOvary"
  age3 = "AgeBreastContralateral"
  ncancer = 2
}

if (model=="MMRpro") {
  ca1 = "AffectedColon"
  ca2 = "AffectedEndometrium"
  ca3 = "AffectedGastric"
  age1 = "AgeColon"
  age2 = "AgeEndometrium"
  age3 = "AgeGastric"
  ncancer = 3
}

if (model=="pancpro"){
  ca1 = "AffectedPancreas"
  age1 = "AgePancreas"
  ncancer = 1
}

if (model=="melapro") {
  ca1 = "AffectedSkin"
  age1 = "AgeSkin"
  ncancer = 1
}

if (ncancer==1){
  family=subset(fff,!(fff[,ca1]==0 & fff$Relation!=1 & fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
}
if (ncancer==2){
  family=subset(fff,!(fff[,ca1]==0 & fff[,ca2]==0 & fff$Relation!=1 & fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
}
if (ncancer==3){
  family=subset(fff,!(fff[,ca1]==0 & fff[,ca2]==0 & fff[, ca3] == 0 & fff$Relation!=1 & fff$Relation!=4 & fff$Relation!=5 & fff$Relation!=7))
}


  cnslees = family[family$Relation==1,]
  cnslID  = family$ID[family$Relation==1] #cnslees$FamilyID

### Get info on number of relatives per family
### _____________________________________________________________________________________

#NA -> family[,ca1][family[,ca1] == 1]
#NA -> family[,ca2][family[,ca2] == 1]
#famunique = unique(family$FamilyID) 

### Get only needed columns
### _____________________________________________________________________________________

if (model=="brcapro") {
 family_vars = c("ID", "Relation", "Gender", "FatherID", "MotherID", "AffectedBreast", "AffectedOvary",
	       "AgeBreast", "AgeOvary", "AgeBreastContralateral", "Twins", "ethnic",
         "BRCA1", "BRCA2", "TestOrder", "ER", "PR", "CK14", "CK5.6", "HER2")
}
if (model=="MMRpro"){
  family_vars = c("ID", "Relation", "Gender", "FatherID", "MotherID", "AffectedColon", "AffectedEndometrium",
                  "AffectedGastric", "AgeColon", "AgeEndometrium", "AgeGastric", "Twins", "ethnic",
                  "MLH1", "MSH2", "MSH6", "TestOrder", "MSI", "location")
}
if (model=='pancpro'){
  family_vars = c("ID", "Relation", "Gender", "FatherID", "MotherID", "AffectedPancreas",
                  "AgePancreas", "Twins", "ethnic")
}
if (model=='melapro'){
  family_vars = c("ID", "Relation", "Gender", "FatherID", "MotherID", "AffectedSkin",
                  "AgeSkin", "Twins", "ethnic",
                  "P16", "TestOrder")
}

 fam=family[,family_vars]

### Add family members according to Median number in Colorectal Data
### _____________________________________________________________________________________

cnsl  = fam$ID[fam$Relation == 1]
cnslsex = fam$Gender[fam$Relation == 1]

nhus    = sum(as.logical(fam$Relation == 14)) 			# Husband
nsis    = sum(as.logical((fam$Relation == 2|fam$Relation == 16) & fam$Gender == 0)) # Sister
nbro    = sum(as.logical((fam$Relation == 2|fam$Relation == 16) & fam$Gender == 1)) # Brother
ndaught = sum(as.logical(fam$Relation == 3 & fam$Gender == 0))	# Daughter
nson    = sum(as.logical(fam$Relation == 3 & fam$Gender == 1))	# Son
nmom    = sum(as.logical(fam$Relation == 4 & fam$Gender == 0))	# Mother
ndad    = sum(as.logical(fam$Relation == 4 & fam$Gender == 1))	# Father
npgrma  = sum(as.logical(fam$Relation == 5 & fam$Gender == 0))	# Paternal Grandmother
npgrpa  = sum(as.logical(fam$Relation == 5 & fam$Gender == 1))	# Paternal Grandfather
npaunt  = sum(as.logical(fam$Relation == 6 & fam$Gender == 0))	# Paternal Aunt
npuncl  = sum(as.logical(fam$Relation == 6 & fam$Gender == 1))	# Paternal Uncle
nmgrma  = sum(as.logical(fam$Relation == 7 & fam$Gender == 0))	# Maternal Grandmother
nmgrpa  = sum(as.logical(fam$Relation == 7 & fam$Gender == 1))	# Maternal Grandfather
nmaunt  = sum(as.logical(fam$Relation == 8 & fam$Gender == 0))	# Maternal Aunt
nmuncl  = sum(as.logical(fam$Relation == 8 & fam$Gender == 1))	# Maternal Uncle

nfam    = dim(fam)[1]
nhsd    = sum(as.logical(fam$Relation == 1|fam$Relation == 14|fam$Relation == 3))
ntwin   = sum(as.logical(fam$Relation == 16))
ntwin1  = sum(as.logical(fam$Twins == 1)) 			# Twins
ntwin2  = sum(as.logical(fam$Twins == 2)) 			# Twins
ntwin3  = sum(as.logical(fam$Twins == 3)) # Twins
ntwin4  = sum(as.logical(fam$Twins == 4)) # Twins
ntwin5  = sum(as.logical(fam$Twins == 5)) # Twins
ntwin6  = sum(as.logical(fam$Twins == 6)) # Twins
ntwin7  = sum(as.logical(fam$Twins == 7)) # Twins
ntwin8  = sum(as.logical(fam$Twins == 8)) # Twins
ntwin9  = sum(as.logical(fam$Twins == 9)) # Twins
ntwin10 = sum(as.logical(fam$Twins == 10)) # Twins

dad     = fam$FatherID[fam$Relation == 1]
mom     = fam$MotherID[fam$Relation == 1]
pgrpa   = fam$FatherID[fam$Relation == 4 & fam$Gender == 1]
pgrma   = fam$MotherID[fam$Relation == 4 & fam$Gender == 1]
mgrpa   = fam$FatherID[fam$Relation == 4 & fam$Gender == 0]
mgrma   = fam$MotherID[fam$Relation == 4 & fam$Gender == 0]

if (nhus>0) {
  spouse = fam$ID[fam$Relation == 14]
  spousesex = fam$Gender[fam$Relation == 14]
}
if (nhus==0) {
  spouse=0 
  if (cnslsex==1) {spousesex = 0}
  if (cnslsex==0) {spousesex = 1}
}
if (spouse==0 & (ndaught>0 | nson>0)) {
  if (cnslsex==1) {fam$MotherID[fam$Relation == 3]=0}
  if (cnslsex==0) {fam$FatherID[fam$Relation == 3]=0}
}

max.ID  = max(fam$ID,fam$FatherID,fam$MotherID)
if (ntwin1==1){
  fam$Twins[fam$Twins==1]=0
}

if (ntwin2==1){
  fam$Twins[fam$Twins==2]=0
}
if (ntwin3==1){
  fam$Twins[fam$Twins==3]=0
}
if (ntwin4==1){
  fam$Twins[fam$Twins==4]=0
}
if (ntwin5==1){
  fam$Twins[fam$Twins==5]=0
}
if (ntwin6==1){
  fam$Twins[fam$Twins==6]=0
}
if (ntwin7==1){
  fam$Twins[fam$Twins==7]=0
}
if (ntwin8==1){
  fam$Twins[fam$Twins==8]=0
}
if (ntwin9==1){
  fam$Twins[fam$Twins==9]=0
}
if (ntwin10==1){
  fam$Twins[fam$Twins==10]=0
}

# need to make model-specific from here down

tempfam = NULL

    if (nfam==1|nfam==nhsd|mom==0|nmom<col.med.mom){
        diff=col.med.mom-nmom
         for (x in 1:diff){
          mom=ID=max.ID=max.ID+1
          fam$MotherID[fam$Relation==1]=mom
          mgrpa=FatherID=max.ID=mom+1
          mgrma=MotherID=max.ID=mom+2
          Relation=4
          Gender=0
          FatherID=MotherID=0
          add=data.frame(ID,Relation,Gender,FatherID,MotherID)
          tempfam=rbind(tempfam,add)
        }
      }

        if (nfam==1|nfam==nhsd|dad==0|ndad<col.med.dad){
           diff=col.med.dad-ndad
           for (x in 1:diff){
             dad=ID=max.ID=max.ID+1
              fam$FatherID[fam$Relation==1]=dad
              pgrpa=FatherID=max.ID=dad+1
              pgrma=MotherID=max.ID=dad+2
              Relation=4
              Gender=1
              FatherID=MotherID=0
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nfam==1|nfam==nhsd|nmgrma<col.med.mgrma){
           diff=col.med.mgrma-nmgrma
           for (x in 1:diff){
              mgrma=ID=max.ID=max.ID+1
              Relation=7
              Gender=0
              FatherID=MotherID=0
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nfam==1|nfam==nhsd|nmgrpa<col.med.mgrpa){
           diff=col.med.mgrpa-nmgrpa
           for (x in 1:diff){
              mgrpa=ID=max.ID=max.ID+1
              Relation=7
              Gender=1
              FatherID=MotherID=0
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nfam==1|nfam==nhsd|npgrma<col.med.pgrma){
           diff=col.med.pgrma-npgrma
           for (x in 1:diff){
              pgrma=ID=max.ID=max.ID+1
              Relation=5
              Gender=0
              FatherID=MotherID=0
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nfam==1|nfam==nhsd|npgrpa<col.med.pgrpa){
           diff=col.med.pgrpa-npgrpa
           for (x in 1:diff){
              pgrpa=ID=max.ID=max.ID+1
              Relation=5
              Gender=1
              FatherID=MotherID=0
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nsis<col.med.sis){
           diff=col.med.sis-nsis
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=2
              Gender=0
              FatherID=fam$FatherID[fam$Relation==1] 
              MotherID=fam$MotherID[fam$Relation==1]
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nbro<col.med.bro){
           diff=col.med.bro-nbro
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=2
              Gender=1
              FatherID=fam$FatherID[fam$Relation==1]
              MotherID=fam$MotherID[fam$Relation==1]
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (npaunt<col.med.paunt){
           diff=col.med.paunt-npaunt
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=6
              Gender=0
              if (any(fam$ID==dad)) {
              FatherID=fam$FatherID[fam$ID==dad]
              MotherID=fam$MotherID[fam$ID==dad]
              }
              if (!any(fam$ID==dad)) {
                FatherID = tempfam$FatherID[tempfam$ID==dad]
                MotherID=tempfam$MotherID[tempfam$ID==dad]
              }
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (npuncl<col.med.puncl){
           diff=col.med.puncl-npuncl
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=6
              Gender=1
              if (any(fam$ID==dad)) {
                FatherID=fam$FatherID[fam$ID==dad]
                MotherID=fam$MotherID[fam$ID==dad]
              }
              if (!any(fam$ID==dad)) {
                FatherID = tempfam$FatherID[tempfam$ID==dad]
                MotherID=tempfam$MotherID[tempfam$ID==dad]
              }
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nmaunt<col.med.maunt){
           diff=col.med.maunt-nmaunt
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=8
              Gender=0
              if (any(fam$ID==mom)) {
                FatherID=fam$FatherID[fam$ID==mom]
                MotherID=fam$MotherID[fam$ID==mom]
              }
              if (!any(fam$ID==mom)) {
                FatherID = tempfam$FatherID[tempfam$ID==mom]
                MotherID=tempfam$MotherID[tempfam$ID==mom]
              }
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nmuncl<col.med.muncl){
           diff=col.med.muncl-nmuncl
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=8
              Gender=1
              if (any(fam$ID==mom)) {
                FatherID=fam$FatherID[fam$ID==mom]
                MotherID=fam$MotherID[fam$ID==mom]
              }
              if (!any(fam$ID==mom)) {
                FatherID = tempfam$FatherID[tempfam$ID==mom]
                MotherID=tempfam$MotherID[tempfam$ID==mom]
              }
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (ndaught<col.med.daught){
           diff=col.med.daught-ndaught
           for (x in 1:diff){
              ID=max.ID=max.ID+1
              Relation=3
              Gender=0
	      if (fam$Gender[fam$Relation == 1] == 0){MotherID=cnsl
                 FatherID=spouse} 
              else{FatherID=cnsl
                 MotherID=spouse}
             add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

        if (nson<col.med.son){
           diff=col.med.son-nson
           for (x in 1:diff){
               ID=max.ID=max.ID+1
              Relation=3
              Gender=1
	      if (fam$Gender[fam$Relation == 1] == 0){MotherID=cnsl
                 FatherID=spouse} 
              else{FatherID=cnsl
                 MotherID=spouse}
              add=data.frame(ID,Relation,Gender,FatherID,MotherID)
              tempfam=rbind(tempfam,add)
           }
        }

tempfam$ethnic=fam$ethnic[fam$Relation==1]

if (model=="brcapro") {
  tempfam$AffectedBreast = tempfam$AffectedOvary = tempfam$AgeBreastContralateral = tempfam$Twins = 0
  tempfam$AgeBreast = tempfam$AgeOvary = 1
  tempfam$BRCA1 = tempfam$BRCA2 = tempfam$TestOrder = 0
  tempfam$ER = tempfam$PR = tempfam$CK14 = tempfam$CK5.6 = tempfam$HER2 = 0
  tempfam = tempfam[,family_vars]
}

if (model=="MMRpro") {
  tempfam$AffectedColon = tempfam$AffectedEndometrium = tempfam$AffectedGastric = tempfam$Twins = 0
  tempfam$AgeColon = tempfam$AgeEndometrium = tempfam$AgeGastric = 1
  tempfam$MSI = tempfam$location = 0
  tempfam$MLH1 = tempfam$MSH2 = tempfam$MSH6 = tempfam$TestOrder = 0
  tempfam = tempfam[,family_vars]
}

if (model=="pancpro") {
  tempfam$AffectedPancreas = tempfam$Twins = 0
  tempfam$AgePancreas = 1
  tempfam = tempfam[,family_vars]
}

if (model=="melapro") {
  tempfam$AffectedSkin = tempfam$Twins = 0
  tempfam$AgeSkin = 1
  tempfam$P16 = tempfam$TestOrder = 0
  tempfam = tempfam[,family_vars]
}

fam = rbind(fam, tempfam)

if (fam$FatherID[fam$Relation == 1]==0) {
  fam$FatherID[fam$Relation == 1]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
  fam$FatherID[fam$Relation == 2]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
  fam$FatherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
}

if (fam$MotherID[fam$Relation == 1]==0) {
  fam$MotherID[fam$Relation == 1]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
  fam$MotherID[fam$Relation == 2]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
  fam$MotherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
}

if (any(fam$FatherID[fam$Relation == 2]==0)) {
  fam$FatherID[fam$Relation == 2]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
  fam$FatherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
}

if (any(fam$MotherID[fam$Relation == 2]==0)) {
  fam$MotherID[fam$Relation == 2]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
  fam$MotherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
}

if (ntwin>0) {if (fam$FatherID[fam$Relation == 16]==0) {
  fam$FatherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 1]
}}

if (ntwin>0) {if (fam$MotherID[fam$Relation == 16]==0) {
  fam$MotherID[fam$Relation == 16]=fam$ID[fam$Relation == 4 & fam$Gender == 0]
}}

if (fam$FatherID[fam$Relation == 4 & fam$Gender == 1]==0) {
  fam$FatherID[fam$Relation == 4 & fam$Gender == 1]=fam$ID[fam$Relation == 5 & fam$Gender == 1]
  fam$FatherID[fam$Relation == 6]=fam$ID[fam$Relation == 5 & fam$Gender == 1]
}

if (fam$MotherID[fam$Relation == 4 & fam$Gender == 1]==0) {
  fam$MotherID[fam$Relation == 4 & fam$Gender == 1]=fam$ID[fam$Relation == 5 & fam$Gender == 0]
  fam$MotherID[fam$Relation == 6]=fam$ID[fam$Relation == 5 & fam$Gender == 0]
}

if (fam$FatherID[fam$Relation == 4 & fam$Gender == 0]==0) {
  fam$FatherID[fam$Relation == 4 & fam$Gender == 0]=fam$ID[fam$Relation == 7 & fam$Gender == 1]
  fam$FatherID[fam$Relation == 8]=fam$ID[fam$Relation == 7 & fam$Gender == 1]
}

if (fam$MotherID[fam$Relation == 4 & fam$Gender == 0]==0) {
  fam$MotherID[fam$Relation == 4 & fam$Gender == 0]=fam$ID[fam$Relation == 7 & fam$Gender == 0]
  fam$MotherID[fam$Relation == 8]=fam$ID[fam$Relation == 7 & fam$Gender == 0]
}

if (cnslsex==1) {fam$FatherID[fam$Relation == 3]=fam$ID[fam$Relation == 1]}
if (cnslsex==0) {fam$MotherID[fam$Relation == 3]=fam$ID[fam$Relation == 1]}

if (cnslsex==1 & spousesex==1) { 
       fam$MotherID[fam$Relation == 3]=0}
if (cnslsex==0 & spousesex==0) {
       fam$FatherID[fam$Relation == 3]=0}

    ninlaw = sum(as.logical(fam$Relation == 15)) 		# Brother or sister in law
    nn     = sum(as.logical(fam$Relation == 13))		# Niece/Nephew
    nnlist = as.numeric(fam$ID[fam$Relation==13]) 		# Niece/Nephew IDs

    if (nn>0) {
      holdtwin=0
      nntwin=0
    for (x in 1:nn) {

       nndad = fam$FatherID[fam$ID==nnlist[[x]]]
       nnmom = fam$MotherID[fam$ID==nnlist[[x]]]
       nntwin = fam$Twins[fam$ID==nnlist[[x]]]
       bro   = sum(as.numeric(fam$ID==nndad & (fam$Relation==2|fam$Relation==16) & fam$Gender==1))
       sis   = sum(as.numeric(fam$ID==nnmom & (fam$Relation==2|fam$Relation==16) & fam$Gender==0))

       dist.sib  = table(fam$ID[fam$Relation==2|fam$Relation==16]) 		# Brother/Sister
       pdist.sib = dist.sib / sum(dist.sib)
       xdist.sib = as.numeric(names(dist.sib))

       if (length(xdist.sib) > 1) sib = sample(xdist.sib,1,prob=pdist.sib) 
       else                       sib = xdist.sib
       if (holdtwin!=0 & holdtwin==nntwin) sib = holdsib
       parentSex=fam$Gender[fam$ID==sib]

       if (bro==0 & sis==0) {
          if (parentSex==0) {fam$MotherID[fam$ID==nnlist[[x]]]=sib
                             fam$FatherID[fam$ID==nnlist[[x]]]=0}
          else              {fam$MotherID[fam$ID==nnlist[[x]]]=0
                             fam$FatherID[fam$ID==nnlist[[x]]]=sib}
       }

       else if (ninlaw==0) {if (bro==0) {
                                 fam$FatherID[fam$ID==nnlist[[x]]]=0}
                            if (sis==0) {
                                 fam$MotherID[fam$ID==nnlist[[x]]]=0}
       }

       else if (ninlaw!=0) {
          if (parentSex==0) {
             if (sis==0) fam$MotherID[fam$ID==nnlist[[x]]]=sib
             inlaw=fam$FatherID[fam$ID==nnlist[[x]]]
             broInlaw=sum(as.numeric(fam$ID==inlaw & fam$Relation==15 & fam$Gender==1))
             if (broInlaw==0) fam$FatherID[fam$ID==nnlist[[x]]]=0
             }
          else {
             if (bro==0) fam$FatherID[fam$ID==nnlist[[x]]]=sib
             inlaw=fam$MotherID[fam$ID==nnlist[[x]]]
             sisInlaw=sum(as.numeric(fam$ID==inlaw & fam$Relation==15 & fam$Gender==0))
             if (sisInlaw==0) fam$MotherID[fam$ID==nnlist[[x]]]=0
             }
       }

       if (nntwin > 0) {holdtwin = nntwin
          holdsib = sib}
      }
    }

  return(fam)
} #end of LyteSimple


