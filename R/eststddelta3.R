eststddelta3<-function (delta2act) 
{
  k <- delta2act$k
  g <- delta2act$g
  ef <- t(t(delta2act$ef[1:(g + g * k + g * k * k), 1]))
  d3 <- matrix(0, nrow = g + g * k + g * k * k, ncol = g + 
                 g * k + g * k * k + k + k * k)
  dmuzdmu <- diag(rep(1, g * k))
  dsigzdsig <- diag(rep(1, g * k * k))
  dmuzdmumar <- matrix(0, nrow = g * k, ncol = k)
  dmuzdsigmar <- matrix(0, nrow = g * k, ncol = k * k)
  dsigzdsigmar <- matrix(0, nrow = g * k * k, ncol = k * k)
  ki = as.integer(k)
  gi = as.integer(g)
  ef = as.double(ef)
  d2ef = as.double(delta2act$ef)
  dmuzdmu = as.double(dmuzdmu)
  dsigzdsig = as.double(dsigzdsig)
  dmuzdmumar = as.double(dmuzdmumar)
  dmuzdsigmar = as.double(dmuzdsigmar)
  dsigzdsigmar = as.double(dsigzdsigmar)
  k = ki
  g = gi
  NROWdmuzdmu = g * k
  NROWdsigzdsig = g * k * k
  NROWdmuzdmumar = g * k
  NROWdmuzdsigmar = g * k
  NROWdsigzdsigmar = g * k * k
  for (i in 1:g) {
    for (j in 1:k) {
      ef[g + (i - 1) * k + j] = (ef[g + (i - 1) * k + j] - 
                                   d2ef[g + g * k + g * k * k + j])
    }
  }
  for (i in 1:k) {
    for (j in 1:g) {
      ef[g + (j - 1) * k + i] = ef[g + (j - 1) * k + i]/sqrt(d2ef[g + 
                                                                    g * k + g * k * k + k + 1 + (i - 1) * (k + 1)])
    }
  }
  for (lg in 1:g) {
    for (i in 1:k) {
      for (j in 1:k) {
        (ef[g + g * k + (lg - 1) * k * k + (i - 1) * 
              k + j] = ef[g + g * k + (lg - 1) * k * k + 
                            (i - 1) * k + j]/(sqrt(d2ef[g + g * k + g * 
                                                          k * k + k + 1 + (i - 1) * (k + 1)]) * sqrt(d2ef[g + 
                                                                                                            g * k + g * k * k + k + 1 + (j - 1) * (k + 
                                                                                                                                                     1)])))
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      dmuzdmu[(((i - 1) * k + j) - 1) + ((i - 1) * k + 
                                           j - 1) * NROWdmuzdmu + 1] = 1/sqrt(d2ef[g + g * 
                                                                                     k + g * k * k + k + 1 + (j - 1) * (k + 1)])
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        dsigzdsig[(((i - 1) * k * k + (j - 1) * k + l) - 
                     1) + ((i - 1) * k * k + (j - 1) * k + l - 1) * 
                    NROWdsigzdsig + 1] = 1/(sqrt(d2ef[g + g * k + 
                                                        g * k * k + k + 1 + (j - 1) * (k + 1)]) * sqrt(d2ef[g + 
                                                                                                              g * k + g * k * k + k + 1 + (l - 1) * (k + 
                                                                                                                                                       1)]))
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      dmuzdmumar[(((i - 1) * k + j) - 1) + (j - 1) * NROWdmuzdmumar + 
                   1] = -1/sqrt(d2ef[g + g * k + g * k * k + k + 
                                       1 + (j - 1) * (k + 1)])
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        if (j == l) {
          dmuzdsigmar[(((i - 1) * g + l) - 1) + ((j - 
                                                    1) * k + l - 1) * NROWdmuzdsigmar + 1] = (-0.5 * 
                                                                                                (d2ef[g + (i - 1) * k + l] - d2ef[g + g * 
                                                                                                                                    k + g * k * k + l])/(d2ef[g + g * k + g * 
                                                                                                                                                                k * k + k + 1 + (l - 1) * (k + 1)] * sqrt(d2ef[g + 
                                                                                                                                                                                                                 g * k + g * k * k + k + 1 + (l - 1) * (k + 
                                                                                                                                                                                                                                                          1)])))
        }
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (f in 1:k) {
        for (l in 1:k) {
          for (m in 1:k) {
            if (f == m) {
              if ((f == j) && (j != l)) {
                dsigzdsigmar[(((i - 1) * k * k + (j - 
                                                    1) * k + l) - 1) + ((f - 1) * k + m - 
                                                                          1) * NROWdsigzdsigmar + 1] = (-0.5 * 
                                                                                                          d2ef[g + g * k + (i - 1) * k * k + 
                                                                                                                 (j - 1) * k + l]/((d2ef[g + g * k + 
                                                                                                                                           g * k * k + k + 1 + (j - 1) * (k + 
                                                                                                                                                                            1)] * sqrt(d2ef[g + g * k + g * k * 
                                                                                                                                                                                              k + k + 1 + (j - 1) * (k + 1)])) * 
                                                                                                                                     sqrt(d2ef[g + g * k + g * k * k + k + 
                                                                                                                                                 1 + (l - 1) * (k + 1)])))
              }
              if ((f == l) && (l != j)) {
                dsigzdsigmar[(((i - 1) * k * k + (j - 
                                                    1) * k + l) - 1) + ((f - 1) * k + m - 
                                                                          1) * NROWdsigzdsigmar + 1] = (-0.5 * 
                                                                                                          d2ef[g + g * k + (i - 1) * k * k + 
                                                                                                                 (j - 1) * k + l]/((d2ef[g + g * k + 
                                                                                                                                           g * k * k + k + 1 + (l - 1) * (k + 
                                                                                                                                                                            1)] * sqrt(d2ef[g + g * k + g * k * 
                                                                                                                                                                                              k + k + 1 + (l - 1) * (k + 1)])) * 
                                                                                                                                     sqrt(d2ef[g + g * k + g * k * k + k + 
                                                                                                                                                 1 + (j - 1) * (k + 1)])))
              }
              if ((f == j) && (j == l)) {
                dsigzdsigmar[(((i - 1) * k * k + (j - 
                                                    1) * k + l) - 1) + ((f - 1) * k + m - 
                                                                          1) * NROWdsigzdsigmar + 1] = (-d2ef[g + 
                                                                                                                g * k + (i - 1) * k * k + (j - 1) * 
                                                                                                                k + l]/(d2ef[g + g * k + g * k * k + 
                                                                                                                               k + 1 + (l - 1) * (k + 1)] * d2ef[g + 
                                                                                                                                                                   g * k + g * k * k + k + 1 + (l - 1) * 
                                                                                                                                                                   (k + 1)]))
              }
            }
          }
        }
      }
    }
  }
  dim(dmuzdmu) <- c(g * k, g * k)
  dim(dsigzdsig) <- c(g * k * k, g * k * k)
  dim(dmuzdmumar) <- c(g * k, k)
  dim(dmuzdsigmar) <- c(g * k, k * k)
  dim(dsigzdsigmar) <- c(g * k * k, k * k)
  dim(ef) <- c(g + g * k + g * k * k, 1)
  d3[1:g, 1:g] <- diag(rep(1, g))
  d3[(g + 1):(g + g * k), (g + 1):(g + g * k)] <- dmuzdmu
  d3[(g + g * k + 1):(g + g * k + g * k * k), (g + g * k + 
                                                 1):(g + g * k + g * k * k)] <- dsigzdsig
  d3[(g + 1):(g + g * k), (g + g * k + g * k * k + 1):(g + 
                                                         g * k + g * k * k + k)] <- dmuzdmumar
  d3[(g + 1):(g + g * k), (g + g * k + g * k * k + k + 1):(g + 
                                                             g * k + g * k * k + k + k * k)] <- dmuzdsigmar
  d3[(g + g * k + 1):(g + g * k + g * k * k), (g + g * k + 
                                                 g * k * k + k + 1):(g + g * k + g * k * k + k + k * k)] <- dsigzdsigmar
  return(list(ef = ef, delta2act = delta2act, d3 = d3))
}