#' interp_weightedMean
#'
#' Interpolation by consideration
#' @param X, and n by m numeric matrix to with missing values to be completed
#' @param sim, a list of dist objects indicating the dissimilarity of each element
#' @param index, a numeric vector of the index of matrix X to which each dist object corresponds. max(index) must be less than or equal to the max index of the matrix X
#' @param alpha, a numeric vector of values, 0<=alpha<=1, indicating the relative weight of each dist object.

interp_weightedMean <- function(X,sim,index,alpha){
  if(!all(unlist(lapply(sim,function(x) attr(x,'Size')))==dim(X)[index])){stop('all dist objects in sim must be dist objects')}
  if(!all(unlist(lapply(sim,is.dist)))){stop('all members of sim must be dist objects')}
  if(!is.list(sim)){stop('sim must be a list')}
  for(i in sim){}
}
