#' sim_linear_func
#' 
#' sim_linear_func linear conversion of dist object to a similarity matrix
#' @param D dist object
#' @param c numeric value to linearly scale similarity
#' @return simiarity matrix 
#' @export 
sim_linear_func <- function(D,c=1){
  1-(c*as.matrix(D))
}

#' sim_exp_func
#' 
#' sim_exp_func exponential conversion of dist object to a similarity matrix
#' @param D dist object
#' @param e numeric value to exponentially scale similarity
#' @param c numeric value to linearly scale similarity
#' @return simiarity matrix 
#' @export 
sim_exp_func <- function(D,e=2,c=1){
  e^(1-(c*as.matrix(D)))
}