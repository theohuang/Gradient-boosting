#' Apply Functions Over Array Margins or Vector
#'
#' Identical to the \code{apply} function, except if X is a vector, it simply applies 
#' \code{FUN} to X. 
#' @param X an array that may be a vector
#' @param MARGIN a vector giving the subscripts which the function will be applied over.
#' @param FUN the function to be applied
#' @param ... optional arguments to \code{FUN}
#' @keywords apply
safeApply = function(X, MARGIN, FUN, ...){
  if (is.vector(X)) {
    return(FUN(X))
  } else {
    return(apply(X, MARGIN, FUN, ...))
  }
}