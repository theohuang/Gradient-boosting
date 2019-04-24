#' Helper Function to Generate Possible Combinations of k Simultaneous Mutations
#' 
#' Returns a vector of genotypes for possible combinations of k simultaneous 
#' mutations based on an input vector of mutation names.
#' @param k number of simultaneous mutations
#' @param mutations character vector of mutation names
#' @param indepGenes factor vector of the same length as \code{mutations}, 
#' with each level representing a different independent gene.
#' @family multigene
getKCombs = function(k, mutations, indepGenes){
  if (length(levels(indepGenes))==length(mutations)) {
    # Case if there are no regions (i.e. each element in the vector of mutations 
    # corresponds to an independent genes)
    mutKCombs = unlist(combn(mutations, k, FUN=paste, collapse="."))
  } else { 
    # Case if there are regions
    indepCombs = combn(levels(indepGenes), k)
    mutKCombs = unlist(apply(indepCombs, 2, function(x){
      do.call(paste, c(expand.grid(lapply(1:k, function(i){
        mutations[indepGenes==x[i]]})), sep="."))
    }))
  }
  
  return(mutKCombs)
}


#' Helper Function to Generate Possible Genotypes
#' 
#' Returns a vector of possible genotypes based on an input vector of mutation 
#' names, assuming no more than \code{maxK} simultaneous mutations. Possible  
#' genotypes include the baseline, single mutations, and all mutation 
#' combinations up to \code{maxK} (in that order) for a total of 
#' \code{1 + length(mutations) + [length(mutations) choose 2] + 
#' [length(mutations) choose 3] + ...[length(mutations) choose maxK]} 
#' possibilities. 
#' @param mutations character vector of mutation names
#' @param indepGenes factor vector of the same length as \code{mutations}, 
#' with each level representing a different independent gene (optional). 
#' Defaults to the number of levels being equal to the length of the vector, i.e.
#' every mutation corresponds to an independent gene. 
#' @param maxK maximum number of simulataneous mutations to consider. 
#' Defaults to 2. 
#' @details The optional \code{indepGenes} factor variable allows for the extension 
#' to region-specific penetrances.  
#' @family multigene
genPossibleGeno = function(mutations, indepGenes=as.factor(1:length(mutations)), 
                           maxK=2) {
  # Drop unused factor levels in indepGenes
  indepGenes = droplevels(as.factor(indepGenes))
  
  # Ensure that mutations and indepGenes have the same length. 
  if (length(mutations) != length(indepGenes)) {
    stop("Length of mutations needs to be same sam as the length of indepGenes, 
         the factor variable indicating independent genes.")
  }
  
  # Generate mutation combinations if more than one mutation is provided
  if (length(mutations) > 1 && length(levels(indepGenes)) > 1 && maxK>1) {
    mutCombs = unlist(sapply(2:maxK, getKCombs, mutations, indepGenes))
  } else { # If there is only one mutation/independent gene, there are no combinations
    mutCombs = NULL
  }
  
  # Possible genotypes are baseline, one mutation, or 2 up to 
  # maxK combinations of mutations
  out = c("none", mutations, mutCombs)
  return(out)
}


#' Generate a List of Cancer Penetrance Matrices
#'
#' Returns a list of penetrance matrices (survival functions and density
#' functions) for cancer-related deaths, separated by each sex. If there are 
#' multiple mutations, use the minimum of waiting time to cancer onset (age of 
#' cancer diagnosis). The output list also includes the input vector of cancer 
#' names and the vector of possible mutations.
#' @param mutations names of mutations to use
#' @param cancers names of cancers to use
#' @param penCancersF list of cancer penetrance matrices for females
#' @param penCancersM list of cancer penetrance matrices for males
#' @param indepGenes factor vector of the same length as \code{mutations}, 
#' with each level representing a different independent gene (optional). 
#' Defaults to the number of levels being equal to the length of the vector, i.e.
#' every mutation corresponds to an independent gene. 
#' @param maxK maximum number of simulataneous mutations to consider. 
#' @param age.last maximum age to consider plus one (the "age" where you do not develop the disease)
#' Defaults to 2. 
#' @details \code{penCancersF} and \code{penCancersM} are each assumed to be of 
#' length \code{length(mutations)+1} and consist of literature estimates for the 
#' penetrances for non-carriers followed by the penetrances for each of the 
#' mutations in `mutations` (in the same order). The optional \code{indepGenes} 
#' factor variable allows for the extension to region-specific penetrances.  
#' @family multigene exported 
#' @export
genCancerPen = function(mutations, cancers, penCancersF, penCancersM, 
                        indepGenes=as.factor(1:length(mutations)), maxK=2, age.last = 111){
  # Possible genotypes
  PG = genPossibleGeno(mutations, indepGenes, maxK)
  
  # Survival and density matrices
  cancerFSur = array(NA, dim=c(age.last, length(PG), length(cancers)), 
                     dimnames=list(1:age.last, PG, cancers))
  cancerMSur = array(NA, dim=c(age.last, length(PG), length(cancers)), 
                     dimnames=list(1:age.last, PG, cancers))
  cancerFDens = array(NA, dim=c(age.last, length(PG), length(cancers)), 
                      dimnames=list(1:age.last, PG, cancers))
  cancerMDens = array(NA, dim=c(age.last, length(PG), length(cancers)), 
                      dimnames=list(1:age.last, PG, cancers))
  
  # Loop through the cancers
  # Non-carrier and single mutation penetrance should be taken from literature
  # This assumes very nice, idealized data that comes in the same format and order
  # Even if that's not the case, it should be straightforward matching and assigning. 
  for (k in 1:length(cancers)){
    cancerFDens[1:(age.last-1),1:(1+length(mutations)),k] = penCancersF[[cancers[k]]]
    cancerMDens[1:(age.last-1),1:(1+length(mutations)),k] = penCancersM[[cancers[k]]]
  }
  
  # Get the penetrances for the last age in the density matrices
  # Then use the density matrices to get the survival matrices. 
  for(j in 1:(1+length(mutations))){
    cancerFDens[age.last,j,] = 1 - safeApply(cancerFDens[1:(age.last - 1),j,], 2, sum)
    cancerMDens[age.last,j,] = 1 - safeApply(cancerMDens[1:(age.last - 1),j,], 2, sum)
    cancerFSur[,j,] = 1 - safeApply(cancerFDens[,j,], 2, cumsum)
    cancerMSur[,j,] = 1 - safeApply(cancerMDens[,j,], 2, cumsum)
  }
  
  # Loop through the cancers to get the estimates for the double mutations
  if (length(mutations) > 1 && length(levels(indepGenes)) > 1) {
    for (k in 1:length(cancers)){
      # Break up all the mutation pairs into the two mutations
      mutPairs = strsplit(PG[(2+length(mutations)):length(PG)], "\\.")
      # Loop through the double mutations
      for(j in (2+length(mutations)):length(PG)){
        mutPair = mutPairs[[j-(1+length(mutations))]]
        # Survival = product of individual mutations
        cancerFSur[,j,k] = apply(cancerFSur[,mutPair,k], 1, prod)
        cancerMSur[,j,k] = apply(cancerMSur[,mutPair,k], 1, prod)
        # Go back to density
        cancerFDens[,j,k] = -diff(c(1, cancerFSur[,j,k]))
        cancerMDens[,j,k] = -diff(c(1, cancerMSur[,j,k]))
      }
    }
  }
  
  # Output as a list
  CP = list()
  CP$PG = PG # Possible genotypes
  CP$cancers = cancers # Names of cancers, supplied as input
  CP$cancerFSur = cancerFSur # Female survival matrix
  CP$cancerMSur = cancerMSur # Male survival matrix
  CP$cancerFDens = cancerFDens # Female density matrix
  CP$cancerMDens = cancerMDens # Male density matrix
  
  return(CP)
}
