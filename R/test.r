gen_test <- function(n=10,m=10){
  mat = matrix(runif(1:(n*m)),n,m,dimnames = list(1:n,letters[1:m]))
  d_i = dist(mat)
  d_j = dist(t(mat))
  
}

test<-function(dat=mtcars){
  X_orig=data.matrix(dat)
  X = X_orig
  X[indx<-sample(prod(dim(X)),40)] = NA
  X_out1 = interp_weightedMean(X,t=2,alpha = .5,
                               D_i = dist(X),D_j=dist_integer(t(X)),
                               sim_func_i=sim_exp_func,sim_func_j=sim_exp_func)
  cor( data.matrix(mtcars)[indx] , X_out1[indx] )
}