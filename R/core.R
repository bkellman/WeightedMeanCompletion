#' interp_weightedMean_n
#'
#' Interpolation by consideration
#' @param X, and n by m numeric matrix to with missing values to be completed
#' @param sim, a list of dist objects indicating the dissimilarity of each element
#' @param index, a numeric vector of the index of matrix X to which each dist object corresponds. max(index) must be less than or equal to the max index of the matrix X
#' @param alpha, a numeric vector of values, 0<=alpha<=1, indicating the relative weight of each dist object.

interp_weightedMean_n <- function(X,sim,index,alpha){
  if(!all(unlist(lapply(sim,function(x) attr(x,'Size')))==dim(X)[index])){stop('all dist objects in sim must be dist objects')}
  if(!all(unlist(lapply(sim,is.dist)))){stop('all members of sim must be dist objects')}
  if(!is.list(sim)){stop('sim must be a list')}
  for(i in sim){}
}

#' interp_weightedMean
#'
#' Interpolation by consideration of row and column similarity
#' @param X, n by m numeric matrix to with missing values (NA) to be completed
#' @param t, integer indicating the number of smoothing iteractions to perform. t=2 by default
#' @param alpha, (optional) 0<=alpha<=1. Proportion row weights should contribute relative to col weights. aplha=.6 by default
#' @param D_i, (optional) a distance matrix indicating the distance between rows. Uses euclidean distance by default. NA to skip row similarity consideration.
#' @param D_j, (optional) a distance matrix indicating the distance between columns. Uses euclidean distance by default. NA to skip row similarity consideration.
#' @param sim_func_i, (optional) a similarity function describing the conversion of D_i to a similarity metric. Uses $1-10*dist$ (sim_linear_func(x,10)) by default.
#' @param sim_func_j, (optional) a similarity function describing the conversion of D_j to a similarity metric. Uses $100^(1-dist)$ (sim_exp_func(x,100)) by default.
#' @export
#' @author Ben Kellman
#' @example 
#' X = data.matrix(mtcars)
#' X[sample(prod(dim(X)),20)] = NA
#' X_out1 = interp_weightedMean(X,t=5)
#' X_out2 = interp_weightedMean(X,t=10,D_i=dist_integer(X),sim_func_i=sim_exp_func)
interp_weightedMean <- function(X,t=2,alpha=.5,D_i=NULL,D_j=NULL,sim_func_i=NULL,sim_func_j=NULL){
  # Set defaults
  if(is.null(D_i)){ D_i = dist(X) }
  if(class(D_i)!='dist'){ if(is.matrix(D_i) { D_i = as.dist(D_i) }else{stop('D_i must be a matrix or a distance object')} }
  if(is.null(D_j)){ D_i = dist(t(X)) }
  if(class(D_j)!='dist'){ if(is.matrix(D_j) { D_j = as.dist(D_j) }else{stop('D_j must be a matrix or a distance object')} }
  if(is.null(sim_func_i)){ sim_func_i = function(x) sim_lin_func(x,e=10) }
  if(is.null(sim_func_j)){ sim_func_j = function(x) sim_exp_func(x,e=100) }
  
  #X_orig = X
  if(!(is.matrix(X) | is.numeric(X))){stop('X must be a numeric matrix')}
  if(class(D_i)!='dist' | class(D_j)!='dist'){stop('D_i and D_j must be dist objects')}
  if(alpha<0 | alpha>1){stop('alpha is not: 0<=alpha<=1')}
  # normailize dist object
  if(any(D_i>1)){D_i = D_i/max(D_i)}
  if(any(D_j>1)){D_j = D_j/max(D_j)}
  # get similarity from distance objects
  W_i = sim_func_i(D_i)
  W_j = sim_func_j(D_j)
  # get normalization factors
  Omega_i = diag(unique(dim(W_i)))
  diag(Omega_i) = 1/apply(W_i,1,function(k) sum(k,na.rm=TRUE))
  W_i_star = Omega_i %*% W_i
  Omega_j = diag(unique(dim(W_j)))
  diag(Omega_j) = 1/apply(W_j,2,function(k) sum(k,na.rm=TRUE))
  W_j_star = Omega_j %*% W_j
  
  # smooth
  X[is.na(X)]=0
  while(t>0){
    t = t-1
    X =  alpha * W_i_star %*% X  +  (1-alpha) * X %*% W_j_star
  }
  X
}




