kunstack<-function (k, mat) 
{
  retval <- matrix(NA, nrow = k, ncol = length(mat)/k)
  for (i in 1:(length(mat)/k)) {
    retval[1:k, i] <- mat[((i - 1) * k + 1):(i * k), 1]
  }
  return(retval)
}