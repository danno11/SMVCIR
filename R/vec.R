vec<-function (mat) 
{
  retval <- matrix(NA, nrow = ncol(mat) * nrow(mat), ncol = 1)
  for (i in 1:ncol(mat)) {
    retval[((i - 1) * nrow(mat) + 1):((i - 1) * nrow(mat) + 
                                        nrow(mat)), 1] <- mat[, i]
  }
  return(retval)
}