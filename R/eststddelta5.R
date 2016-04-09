eststddelta5<-function (delta4act) 
{
  gprop <- delta4act$gprop
  k <- delta4act$k
  g <- delta4act$g
  ef <- as.matrix(delta4act$ef[(g + 1):(g + g * k + g * k * 
                                          k), 1])
  for (i in 1:g) {
    ef[((i - 1) * k + 1):(i * k), 1] <- sqrt(gprop[i]) * 
      ef[((i - 1) * k + 1):(i * k), 1]
    ef[(g * k + (i - 1) * k * k + 1):(g * k + i * k * k), 
       1] <- sqrt(gprop[i]) * ef[(g * k + (i - 1) * k * 
                                    k + 1):(g * k + i * k * k), 1]
  }
  d5 <- matrix(0, nrow = g * k + g * k * k, ncol = g + g * 
                 k + g * k * k)
  dvdp <- matrix(0, nrow = g * k, ncol = g)
  dDeldp <- matrix(0, nrow = g * k * k, ncol = g)
  dvdmu <- diag(rep(1, g * k))
  dDeldsc <- diag(rep(1, g * k * k))
  gprop = as.double(gprop)
  ki = as.integer(k)
  gi = as.integer(g)
  ef = as.double(ef)
  d4ef = as.double(delta4act$ef)
  dvdp = as.double(dvdp)
  dDeldp = as.double(dDeldp)
  dvdmu = as.double(dvdmu)
  dDeldsc = as.double(dDeldsc)
  g = gi
  k = ki
  NROWdvdp = g * k
  NROWdDeldp = g * k * k
  NROWdvdmu = g * k
  NROWdDeldsc = g * k * k
  for (j in 1:g) {
    for (m in 1:k) {
      dvdp[(((j - 1) * k + m) - 1) + (j - 1) * NROWdvdp + 
             1] = d4ef[g + (j - 1) * k + m]/(2 * sqrt(gprop[j]))
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (m in 1:k) {
        dDeldp[(((i - 1) * k * k + (j - 1) * k + m) - 
                  1) + (i - 1) * NROWdDeldp + 1] = d4ef[g + g * 
                                                          k + (i - 1) * k * k + (j - 1) * k + m]/(2 * 
                                                                                                    sqrt(gprop[i]))
      }
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      dvdmu[(((i - 1) * k + j) - 1) + ((i - 1) * k + j - 
                                         1) * NROWdvdmu + 1] = sqrt(gprop[i])
    }
  }
  for (i in 1:g) {
    for (j in 1:k) {
      for (l in 1:k) {
        dDeldsc[(((i - 1) * k * k + (j - 1) * k + l) - 
                   1) + ((i - 1) * k * k + (j - 1) * k + l - 1) * 
                  NROWdDeldsc + 1] = sqrt(gprop[i])
      }
    }
  }
  dim(ef) <- c(g * k + g * k * k, 1)
  dim(dvdp) <- c(g * k, g)
  dim(dDeldp) <- c(g * k * k, g)
  dim(dvdmu) <- c(g * k, g * k)
  dim(dDeldsc) <- c(g * k * k, g * k * k)
  d5[1:(g * k), 1:g] <- dvdp
  d5[1:(g * k), (g + 1):(g + g * k)] <- dvdmu
  d5[(g * k + 1):(g * k + g * k * k), 1:g] <- dDeldp
  d5[(g * k + 1):(g * k + g * k * k), (g + g * k + 1):(g + 
                                                         g * k + g * k * k)] <- dDeldsc
  return(list(ef = ef, delta4act = delta4act, d5 = d5))
}