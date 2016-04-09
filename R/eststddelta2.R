eststddelta2<-function (delta1act) 
{
  k <- delta1act$k
  g <- delta1act$g
  ef <- delta1act$ef
  d2 <- diag(rep(1, g + g * k + g * k * k + k + k * k))
  dsigdmu <- matrix(0, nrow = g * k * k, ncol = g * k)
  dsigmardmumar <- matrix(0, nrow = k * k, ncol = k)
  ef = as.double(ef)
  dsigdmu = as.double(dsigdmu)
  dsigmardmumar = as.double(dsigmardmumar)
  d1ef = as.double(delta1act$ef)
  ki = as.integer(k)
  gi = as.integer(g)
  g = gi
  k = ki
  NCOLdsigdmu = g * k
  NROWdsigdmu = g * k * k
  NROWdsigmardmumar = k * k
  NCOLdsigmardmumar = k
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        ef[g + g * k + (i - 1) * k * k + (j - 1) * k + 
             l] = ef[g + g * k + (i - 1) * k * k + (j - 
                                                      1) * k + l] - ef[g + (i - 1) * k + j] * ef[g + 
                                                                                                   (i - 1) * k + l]
      }
    }
  }
  for (j in 1:k) {
    for (l in 1:k) {
      ef[g + g * k + g * k * k + k + (j - 1) * k + l] = ef[g + 
                                                             g * k + g * k * k + k + (j - 1) * k + l] - ef[g + 
                                                                                                             g * k + g * k * k + j] * ef[g + g * k + g * k * 
                                                                                                                                           k + l]
    }
  }
  for (l in 1:g) {
    for (i in 1:k) {
      for (f in 1:k) {
        for (j in 1:k) {
          if ((j == i) && (i != f)) {
            dsigdmu[(((l - 1) * k * k + (i - 1) * k + 
                        f) - 1) + ((l - 1) * k + j - 1) * NROWdsigdmu + 
                      1] = -d1ef[g + (l - 1) * k + f]
          }
          if ((j == i) && (i == f)) {
            dsigdmu[(((l - 1) * k * k + (i - 1) * k + 
                        f) - 1) + ((l - 1) * k + j - 1) * NROWdsigdmu + 
                      1] = -2 * d1ef[g + (l - 1) * k + i]
          }
          if ((j == f) && (i != f)) {
            dsigdmu[(((l - 1) * k * k + (i - 1) * k + 
                        f) - 1) + ((l - 1) * k + j - 1) * NROWdsigdmu + 
                      1] = -d1ef[g + (l - 1) * k + i]
          }
        }
      }
    }
  }
  for (i in 1:k) {
    for (f in 1:k) {
      for (j in 1:k) {
        if ((j == i) && (i != f)) {
          dsigmardmumar[(((i - 1) * k + f) - 1) + (j - 
                                                     1) * NROWdsigmardmumar + 1] = -d1ef[g + g * 
                                                                                           k + g * k * k + f]
        }
        if ((j == i) && (i == f)) {
          dsigmardmumar[(((i - 1) * k + f) - 1) + (j - 
                                                     1) * NROWdsigmardmumar + 1] = -2 * d1ef[g + 
                                                                                               g * k + g * k * k + i]
        }
        if ((j == f) && (i != f)) {
          dsigmardmumar[(((i - 1) * k + f) - 1) + (j - 
                                                     1) * NROWdsigmardmumar + 1] = -d1ef[g + g * 
                                                                                           k + g * k * k + i]
        }
      }
    }
  }
  dim(dsigdmu) <- c(g * k * k, g * k)
  dim(dsigmardmumar) <- c(k * k, k)
  dim(ef) <- dim(delta1act$ef)
  d2[(g + g * k + 1):(g + g * k + g * k * k), (g + 1):(g + 
                                                         g * k)] <- dsigdmu
  d2[(g + g * k + g * k * k + k + 1):(g + g * k + g * k * k + 
                                        k + k * k), (g + g * k + g * k * k + 1):(g + g * k + 
                                                                                   g * k * k + k)] <- dsigmardmumar
  return(list(ef = ef, delta1act = delta1act, d2 = d2))
}