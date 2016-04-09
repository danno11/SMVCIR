eststddelta4<-function (delta3act) 
{
  k <- delta3act$k
  g <- delta3act$g
  gprop <- delta3act$gprop
  ef <- delta3act$ef
  wvef <- ef
  d4 <- diag(rep(1, g + g * k + g * k * k))
  dscdp <- matrix(0, nrow = g * k * k, ncol = g)
  dscds <- matrix(0, nrow = g * k * k, ncol = g * k * k)
  ki = as.integer(k)
  gi = as.integer(g)
  gprop = as.double(gprop)
  d3ef = as.double(delta3act$ef)
  ef = as.double(ef)
  wvef = as.double(wvef)
  dscdp = as.double(dscdp)
  dscds = as.double(dscds)
  k = ki
  g = gi
  NROWdscdp = g * k * k
  NROWdscds = g * k * k
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        weightvar = 0
        for (m in 1:g) {
          weightvar = weightvar + gprop[m] * ef[g + g * 
                                                  k + (m - 1) * k * k + (j - 1) * k + l]
        }
        wvef[g + g * k + (i - 1) * k * k + (j - 1) * 
               k + l] = weightvar
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        ef[g + g * k + (i - 1) * k * k + (j - 1) * k + 
             l] = ef[g + g * k + (i - 1) * k * k + (j - 
                                                      1) * k + l] - wvef[g + g * k + (i - 1) * k * 
                                                                           k + (j - 1) * k + l]
      }
    }
  }
  for (l in 1:g) {
    for (q in 1:k) {
      for (f in 1:k) {
        for (j in 1:g) {
          dscdp[(((l - 1) * k * k + (q - 1) * k + f) - 
                   1) + (j - 1) * NROWdscdp + 1] = -d3ef[g + 
                                                           g * k + (j - 1) * k * k + (q - 1) * k + f]
        }
      }
    }
  }
  for (l in 1:g) {
    for (q in 1:k) {
      for (h in 1:g) {
        for (m in 1:k) {
          for (f in 1:k) {
            for (j in 1:k) {
              if (((h == l) && (m == q)) & (j == f)) {
                dscds[(((l - 1) * k * k + (q - 1) * k + 
                          f) - 1) + ((h - 1) * k * k + (m - 1) * 
                                       k + j - 1) * NROWdscds + 1] = 1 - gprop[l]
              }
              if (((h != l) & (m == q)) & (j == f)) {
                dscds[(((l - 1) * k * k + (q - 1) * k + 
                          f) - 1) + ((h - 1) * k * k + (m - 1) * 
                                       k + j - 1) * NROWdscds + 1] = -gprop[h]
              }
            }
          }
        }
      }
    }
  }
  dim(dscdp) <- c(g * k * k, g)
  dim(dscds) <- c(g * k * k, g * k * k)
  dim(ef) <- dim(delta3act$ef)
  d4[(g + g * k + 1):(g + g * k + g * k * k), 1:g] <- dscdp
  d4[(g + g * k + 1):(g + g * k + g * k * k), (g + g * k + 
                                                 1):(g + g * k + g * k * k)] <- dscds
  return(list(ef = ef, delta3act = delta3act, d4 = d4))
}