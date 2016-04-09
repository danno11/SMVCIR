eststackvar<-function (delta5act) 
{
  gprop <- delta5act$gprop
  g <- delta5act$g
  k <- delta5act$k
  ef <- matrix(0, nrow = g * k + g * k * k + g * k, ncol = 1)
  ef[(1:(g * k + g * k * k)), 1] <- delta5act$ef
  IDP <- diag(rep(1, g * k + g * k * k + g * k))
  Perm <- IDP
  for (i in 1:g) {
    for (z in 1:k) {
      Perm[g * k + k * k * (i - 1) + 1 + (z - 1) * (k + 
                                                      1), ] <- IDP[g * k + g * k * k + z + (i - 1) * 
                                                                     k, ]
      Perm[g * k + g * k * k + z + (i - 1) * k, ] <- IDP[g * 
                                                           k + k * k * (i - 1) + 1 + (z - 1) * (k + 1), 
                                                         ]
    }
  }
  ef <- Perm %*% ef
  return(list(ef = ef, delta5act = delta5act, Perm = Perm))
}