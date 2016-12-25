gen_test <- function(n=10,m=10){
  mat = matrix(runif(1:(n*m)),n,m,dimnames = list(1:n,letters[1:m]))
  d_i = dist(mat)
  d_j = dist(t(mat))
  
}







