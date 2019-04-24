#' Generate Unique IDs and Genders in a Branch and Link Parents' IDs
#' 
#' Returns a data frame with 4 columns: \code{ID}, \code{Gender}, \code{MotherID}, 
#' and \code{FatherID}. There will be a row for each of the progeny, plus a row for  
#' any parent that were simulated/not blood a relation of the main family, i.e. 
#' spouses or the grandparents. \code{ID} stores the individuals' ID, \code{Gender} 
#' stores their genders, and \code{MotherID} and \code{FatherID} stores their 
#' parent's IDs. A parental ID of 0 corresponds to a parent who is not part of the 
#' family matrix. 
#' @param nDaughters number of daughters in branch
#' @param nSons number of sons in branch
#' @param mothID mother's ID, if known (optional)
#' @param fathID father's ID, if known (optional)
#' @param lastID the largest ID number that is already in use. 
#' @param includeParents boolean flag, indicating whether to assign IDs to the 
#' parents. Defaults to TRUE
#' Defaults to 0 (i.e. new IDs will start from 1)
#' @details If the \code{mothID} and/or \code{fathID} is not known, the mother 
#' and/or father is assumed to have been simulated. 
#' @family simulations
sim.linkParents = function(nDaughters, nSons, mothID, fathID, 
                           lastID=0, includeParents=TRUE) {
  # Generate IDs for all children, as well as parents if necessary
  ID = (lastID + 1):(lastID + missing(mothID)*includeParents + 
                       missing(fathID)*includeParents + 
                       nDaughters + nSons)
  
  # Initialize vectors
  Gender = numeric(0) # gender, 0 for females and 1 for males
  MotherID = numeric(0) # Mother IDs
  FatherID = numeric(0) # Father IDs
  
  # If no mother ID is supplied, add her info
  if (missing(mothID)) {
    # Assign mother an ID, if necessary
    if (includeParents==TRUE) { 
      mothID = lastID+1
      lastID = mothID
    } else { 
      mothID = 0
      ID = c(0, ID)
    }
    
    # Add mother to the Gender vector as 0
    Gender = c(Gender, 0) 
    # Add 0 for her parent IDs
    FatherID = c(FatherID, 0)
    MotherID = c(MotherID, 0)
  } else if (includeParents==FALSE) {
    mothID = 0
  }
  
  # If no father ID is supplied, add his info
  if (missing(fathID)) {
    # Assign father an ID, if necessary
    if (includeParents==TRUE) { 
      fathID = lastID+1
      lastID = fathID
    } else { 
      fathID = 0
      ID = c(0, ID)
    }
    
    # Add father to the Gender vector as 1
    Gender = c(Gender, 1)
    # Add 0 for his parent IDs
    FatherID = c(FatherID, 0)
    MotherID = c(MotherID, 0)
  } else if (includeParents==FALSE) {
    fathID = 0
  }
  
  # Add daughters' genders (0) and sons' genders (1)
  Gender = c(Gender, rep(0, nDaughters), rep(1, nSons))
  # Add children's parent IDs
  MotherID = c(MotherID, rep(mothID, nDaughters + nSons))
  FatherID = c(FatherID, rep(fathID, nDaughters + nSons))
  
  return(data.frame(ID, Gender, MotherID, FatherID))
}