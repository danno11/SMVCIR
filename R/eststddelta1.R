eststddelta1<-function (clt1act) 
{
  gprop <- clt1act$gprop
  g <- clt1act$g
  k <- clt1act$k
  ef <- clt1act$eu
  for (i in 1:g) {
    ef[(g + (i - 1) * k + 1):(g + i * k), 1] <- ef[(g + (i - 
                                                           1) * k + 1):(g + i * k), 1]/gprop[i]
    ef[(g + g * k + (i - 1) * k * k + 1):(g + g * k + i * 
                                            k * k), 1] <- ef[(g + g * k + (i - 1) * k * k + 1):(g + 
                                                                                                  g * k + i * k * k), 1]/gprop[i]
  }
  d1 <- diag(rep(1, g + g * k + g * k * k + k + k * k))
  duvdp <- matrix(0, nrow = g * k + g * k * k, ncol = g)
  duvdmv <- diag(rep(1, times = g * k + g * k * k))
  NROWduvdp <- g * k + g * k * k
  NROWduvdmv <- g * k + g * k * k
  NCOLduvdp <- g
  NCOLduvdmv <- g * k + g * k * k
  clt1eu = as.double(clt1act$eu)
  gprop = as.double(gprop)
  duvdp = as.double(duvdp)
  duvdmv = as.double(duvdmv)
  gi = as.integer(g)
  ki = as.integer(k)
  g = gi
  k = ki
  NROWduvdp = g * k + g * k * k
  NROWduvdmv = g * k + g * k * k
  NCOLduvdp = g
  NCOLduvdmv = g * k + g * k * k
  for (l in 1:g) {
    for (i in 1:k) {
      for (j in 1:g) {
        if (j == l) {
          duvdp[(((l - 1) * k + i) - 1) + (j - 1) * NROWduvdp + 
                  1] = -clt1eu[g + (l - 1) * k + i]/((gprop[j]) * 
                                                       (gprop[j]))
        }
      }
    }
  }
  for (l in 1:g) {
    for (i in 1:k) {
      for (f in 1:k) {
        for (j in 1:g) {
          if (j == l) {
            duvdp[((g * k + (l - 1) * k * k + (i - 1) * 
                      k + f) - 1) + (j - 1) * NROWduvdp + 1] = -clt1eu[g + 
                                                                         g * k + (l - 1) * k * k + (i - 1) * k + 
                                                                         f]/((gprop[j]) * (gprop[j]))
          }
        }
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      duvdmv[(((i - 1) * k + j) - 1) + ((i - 1) * k + j - 
                                          1) * NROWduvdmv + 1] = 1/gprop[i]
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        duvdmv[((g * k + (i - 1) * k * k + (j - 1) * 
                   k + l) - 1) + (g * k + (i - 1) * k * k + (j - 
                                                               1) * k + l - 1) * NROWduvdmv + 1] = 1/gprop[i]
      }
    }
  }
  dim(duvdp) <- c(g * k + g * k * k, g)
  dim(duvdmv) <- c(g * k + g * k * k, g * k + g * k * k)
  d1[(g + 1):(g + g * k + g * k * k), 1:g] <- duvdp
  d1[(g + 1):(g + g * k + g * k * k), (g + 1):(g + g * k + 
                                                 g * k * k)] <- duvdmv
  return(list(ef = ef, clt1act = clt1act, d1 = d1))
}