#' dist_nonparametric
#' 
#' dist_nonparametric is a non-parametric distance function
#' @param X numeric matrix
#' @return dist object 
#' @export 
dist_nonparametric <- function(X){
  if(!(is.matrix(X) | is.numeric(X))){stop('X must be a numeric matrix')}
  dist(X,method='minkowski') 
}

#' dist_binary
#' 
#' dist_binary is a jaccard distance function
#' @param X numeric matrix
#' @return dist object 
#' @export 
dist_binary <- function(X){
  if(!(is.matrix(X) | is.numeric(X))){stop('X must be a numeric matrix')}
  dist(X,method='binary') 
}

#' dist_integer
#' 
#' dist_integer is a manhattan distance function
#' @param X numeric matrix
#' @return dist object 
#' @export 
dist_integer <- function(X){
  if(!(is.matrix(X) | is.numeric(X))){stop('X must be a numeric matrix')}
  dist(X,method='manhattan') 
}