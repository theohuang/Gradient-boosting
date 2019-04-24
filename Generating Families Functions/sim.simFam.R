#' Simulate a Family/Pedigree Matrix
#' 
#' Returns a completed family matrix with the following columns for each member: 
#' \itemize{
#'   \item \code{ID} = Member identifier
#'   \item \code{MotherID} = Mother's identifier number
#'   \item \code{FatherID} = Father's identifier number
#'   \item \code{Gender} = 1 for males, 0 for females
#'   \item \code{isProband} = Indicator for whether or not person is proband
#'   \item \code{CurAge} = Family member's current age, or age at death if dead
#'   \item \code{isAffCancer} variables = Indicators for whether or not person 
#'   had cancer type \code{Cancer}, e.g. \code{isAffBC}, \code{isAffO}, etc, plus
#'   an indicator \code{isAffAny} for having any cancer
#'   \item \code{AgeCancer} variables = Age of cancer diagnosis for cancer type 
#'   \code{Cancer}, e.g. \code{isAffBC}, \code{isAffO}, etc, plus the age of 
#'   first cancer diagnosis \code{AgeAny}. 0 if person did not have that cancer 
#'   type
#'   \item \code{isDead} = Indicator for whether or not person is dead
#' }
#' @param nSibsPatern vector of length \code{2}, indicating the number of sisters 
#' and brothers the father has (does not include father) 
#' @param nSibsMatern vector of length \code{2}, indicating the number of sisters 
#' and brothers the mother has (does not include mother) 
#' @param nSibs vector of length \code{2}, indicating the number of sisters and 
#' brothers the proband has (does not include proband) 
#' @param nGrandchild \code{sum(nSibs)+1} by \code{2} matrix, indicating the number 
#' of daughters and number of sons that the proband and each of her siblings have. 
#' Each row corresponds to the children of one of the members of \code{nChild}, with 
#' the proband as the first row. Alternatively, \code{nGrandchild} can be passed in 
#' as a vector of length \code{2}, indicating the number of daughters and number of 
#' sons that the proband and each of her siblings have, when they each have the same 
#' number.
#' @param q allele frequencies for each gene of interest (named vector)
#' @param CP list of cancer penetrance matrices, separated by each gender
#' @param includeGeno boolean flag, indicating whether or not to include genotype 
#' matrix in output. Defaults to FALSE, and was mostly used for troubleshooting. 
#' @param age.max maximum age to consider. Defaults to 110.
#' @param age.min minimum age to consider. Defaults to 1.  
#' @param includeGrandparents boolean flag, indicating whether to drop the 
#' grandparents from the final pedigree. Defaults to TRUE
#' @param censoring if FALSE, then will assume everyone without any cancers lives
#' to the maximum age. Defaults to TRUE.
#' @param genderPro Can be set to "Female" or "Male" to designate the gender of the proband.
#' If left as NULL, the proband's gender will be randomly generated.
#' @param genoMat genotype matrix for the family. This is a matrix where the number of rows is
#' the number of family members and the number of columns is the number of mutations. If
#' left as NULL (the default), the genotypes will be generated.
#' @param CurAge vector of current ages for the family. If left as NULL (the default),
#' the current ages will be generated.
#' @param affTime if TRUE, the output will include the cancer times (which may or may not be
#' observed due to censoring). Defaults to FALSE
#' @details Assumes naming conventions for cancers are consistent between arguments.
#' @family simulations export
#' @export
sim.simFam = function(nSibsPatern, nSibsMatern, nSibs, nGrandchild, q, CP, 
                      includeGeno=FALSE, age.max = 110, age.min = 1, includeGrandparents=TRUE,
                      censoring = TRUE, genderPro = NULL, genoMat = NULL,
                      CurAge = NULL, affTime = FALSE) {
  # Add the father to the number of male children in his branch
  if (length(nSibsPatern) != 2) {
    stop("nSibsPatern needs to be a numeric vector of length 2")
  }
  nChildPatern = nSibsPatern
  nChildPatern[2] = nChildPatern[2] + 1
  
  # Add the mother to the number of female children in her branch
  if (length(nSibsMatern) != 2) {
    stop("nSibsMatern needs to be a numeric vector of length 2")
  }
  nChildMatern = nSibsMatern
  nChildMatern[1] = nChildMatern[1] + 1
  
  # Add the proband to the number of female or male children in her branch/generation
  if (length(nSibs) != 2) {
    stop("nSibs needs to be a numeric vector of length 2")
  }
  nChild = nSibs
  if(is.null(genderPro)){
    genderPro <- ifelse(rbinom(1, 1, 0.5) == 1, "Female", "Male")
  }
  if(!is.null(genderPro) & !(genderPro %in% c("Female", "Male"))){
    genderPro <- ifelse(rbinom(1, 1, 0.5) == 1, "Female", "Male")
  }
  if(genderPro == "Female"){
    nChild[1] <- nChild[1] + 1
  } else{
    nChild[2] <- nChild[2] + 1
  }
  
  # If nGrandchild is passed in as a vector, use those numbers for each person in nChild
  if (is.null(nrow(nGrandchild))) {
    if (length(nGrandchild) != 2) {
      stop("If nGrandchild is passed in as a vector, it needs to have length 2")
    }
    nGrandchildInBranches = do.call(rbind, rep(list(nGrandchild), length=sum(nChild)))
    # Error if nGrandchild is a matrix with number of rows not equal to the number of 
    # people in nChild
  } else if (nrow(nGrandchild) != sum(nChild)) {
    stop("Number of rows in nGrandchild do not equal number of children in nChild.")
  } else {
    if (ncol(nGrandchild) != 2) {
      stop("If nGrandchild is passed in as a matrix, it needs to have 2 columns")
    }
    nGrandchildInBranches = nGrandchild
  }
  
  # Print warning message if includeGrandparents is FALSE but parents still have siblings
  if (includeGrandparents==FALSE && (nSibsPatern > 0 || nSibsMatern > 0)) {
    warning("You specified includeGrandparents==FALSE, but the mother and/or father 
            still has siblings. These will be included in the family pedigree as 
            unlinked members.")
  }
  
  
  # Generate genotype matrix for family
  if(is.null(genoMat)){
    genoMat = sim.buildGenoMat(q, nChildPatern, nChildMatern, 
                               nChild, nGrandchildInBranches)
  }
  # Total number of people in family
  N = nrow(genoMat)
  
  # Link Father ID, Mother ID, and genders for each branch of family
  # Paternal grandparents and their children 
  # (including father, who is the first male child)
  GenderParentIDsPatern = sim.linkParents(nChildPatern[1], nChildPatern[2], 
                                          includeParents=includeGrandparents)
  lastID = nrow(GenderParentIDsPatern) - 2*(!includeGrandparents)
  
  # Maternal grandparents and their children 
  # (including mother, who is the first female child)
  GenderParentIDsMatern = sim.linkParents(nChildMatern[1], nChildMatern[2], 
                                          lastID=lastID, 
                                          includeParents=includeGrandparents)
  lastID = lastID + nrow(GenderParentIDsMatern) - 2*(!includeGrandparents)
  
  # Proband and siblings (proband is the first female or male child listed)
  if(genderPro == "Female"){
    probandID <- lastID + 1
  } else{
    probandID <- lastID + nChild[1] + 1
  }
  GenderParentIDsChild = sim.linkParents(nChild[1], nChild[2], 
                                         mothID=GenderParentIDsMatern$ID[3], 
                                         fathID=GenderParentIDsPatern$ID[3+nChildPatern[1]], 
                                         lastID=lastID)
  lastID = lastID + nrow(GenderParentIDsChild)
  
  # Proband and siblings' spouses and children
  GenderParentIDsGrandchild = data.frame()
  for (i in 1:nrow(nGrandchildInBranches)) {
    # Mother is blood relative in family, father is spouse
    if (GenderParentIDsChild$Gender[i] == 0) {
      GenderParentIDsGrandchild = rbind(GenderParentIDsGrandchild, 
                                        sim.linkParents(nGrandchildInBranches[i,1], nGrandchildInBranches[i,2], 
                                                        mothID=GenderParentIDsChild$ID[i], 
                                                        lastID=lastID+nrow(GenderParentIDsGrandchild)))
      # Father is blood relative in family, mother is spouse
    } else {
      GenderParentIDsGrandchild = rbind(GenderParentIDsGrandchild, 
                                        sim.linkParents(nGrandchildInBranches[i,1], nGrandchildInBranches[i,2], 
                                                        fathID=GenderParentIDsChild$ID[i], 
                                                        lastID=lastID+nrow(GenderParentIDsGrandchild)))
    }
  }
  
  # Combine genders and parent IDs for all family members
  GenderParentIDs = rbind(GenderParentIDsPatern, GenderParentIDsMatern,
                          GenderParentIDsChild, GenderParentIDsGrandchild)
  
  
  if(is.null(CurAge)){
    if(censoring == TRUE){
      # Generate current ages for each branch of family
      # Paternal grandparents and their children 
      # (including father, who is the first male child)
      CurAgePatern = sim.simCurAgeVar(sum(nChildPatern))
      
      # Maternal grandparents and their children 
      # (including mother, who is the first female child)
      CurAgeMatern = sim.simCurAgeVar(sum(nChildMatern))
      
      # Proband and siblings
      CurAgeChild = sim.simCurAgeVar(sum(nChild), CurAgePatern[3], CurAgeMatern[3])
      
      # Proband and siblings' spouses and children
      CurAgeGrandchild = unlist(lapply(1:nrow(nGrandchildInBranches), function(i){
        # Mother is blood relative in family, father is spouse
        if (GenderParentIDsChild$Gender[i] == 0) {
          sim.simCurAgeVar(sum(nGrandchildInBranches[i,]), mothAge=CurAgeChild[i])
          # Father is blood relative in family, mother is spouse
        } else {
          sim.simCurAgeVar(sum(nGrandchildInBranches[i,]), fathAge=CurAgeChild[i])
        }
      }))
      
      # Combine current ages for all family members
      CurAge = c(CurAgePatern, CurAgeMatern, CurAgeChild, CurAgeGrandchild)
      
      # Ensure that all current ages fall between age.min and age.max
      CurAge[CurAge < age.min] = age.min
      CurAge[CurAge > age.max] = age.max
    } else{
      CurAge = rep(age.max, N)
    }
  }
  
  # Get Affected and Age columns for each cancer, plus the isDead column
  Cancers = sim.simCancerVars(genoMat, GenderParentIDs$Gender, CurAge, CP,
                              age.max = age.max, censoring = censoring,
                              affTime = affTime)
  
  # Designate the proband
  isProband = rep(0, length=N)
  isProband[probandID] = 1
  
  # Put the family pedigree together, including the genotype matrix for all 
  # members, if necessary
  if (includeGeno==TRUE) {
    fam = data.frame(GenderParentIDs, isProband, Cancers, genoMat)
  } else {
    fam = data.frame(GenderParentIDs, isProband, Cancers)
  }
  
  # Drop the grandparents, if necessary
  if (includeGrandparents==FALSE) {
    fam = fam[fam$ID!=0,]
  }
  
  return(fam)
}
