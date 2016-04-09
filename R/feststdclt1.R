feststdclt1<-function (groups, dat) 
{
  g <- groups
  k <- ncol(dat) - 1
  realizes = cbind(matrix(NA, nrow = nrow(dat), ncol = g), 
                   dat[, 1:k], matrix(NA, nrow = nrow(dat), ncol = g * k + 
                                        k * (k + 1)/2 + g * k * (k + 1)/2))
  for (i in 1:g) {
    realizes[, i] = 1 * (dat[, k + 1] == i)
  }
  for (i in 1:g) {
    realizes[, (g + k + (i - 1) * k + 1):(g + k + i * k)] = (1 * 
                                                               (dat[, k + 1] == i)) * realizes[, (g + 1):(g + k)]
  }
  d = 1
  for (i in 1:k) {
    for (j in 1:i) {
      realizes[, (g + k + g * k + d)] = realizes[, g + 
                                                   i] * realizes[, g + j]
      d = d + 1
    }
  }
  for (git in 1:g) {
    for (i in 1:k) {
      for (j in 1:i) {
        realizes[, (g + k + g * k + d)] = (1 * (dat[, 
                                                    k + 1] == git)) * realizes[, g + i] * realizes[, 
                                                                                                   g + j]
        d = d + 1
      }
    }
  }
  meanit <- t(t(colMeans(realizes)))
  dmeanit <- dim(meanit)
  covit <- var(realizes)
  dcovit <- dim(covit)
  NCOL = g + g * k + g * k * k + k + k * k
  NCOLc = g + k + g * k + k * (k + 1)/2 + g * k * (k + 1)/2
  eu <- matrix(0, nrow = g + g * k + g * k * k + k + k * k, 
               ncol = 1)
  vu <- matrix(0, nrow = g + g * k + g * k * k + k + k * k, 
               ncol = g + g * k + g * k * k + k + k * k)
  meanit = as.double(meanit)
  covit = as.double(covit)
  eu = as.double(eu)
  vu = as.double(vu)
  d = 1
  for (i in 1:k) {
    for (j in 1:i) {
      eu[g + g * k + g * k * k + k + (i - 1) * k + j] = meanit[g + 
                                                                 k + g * k + d]
      eu[g + g * k + g * k * k + k + (j - 1) * k + i] = meanit[g + 
                                                                 k + g * k + d]
      d = d + 1
    }
  }
  for (git in 1:g) {
    for (i in 1:k) {
      for (j in 1:i) {
        eu[g + g * k + (git - 1) * k * k + (i - 1) * 
             k + j] = meanit[(g + k + g * k + d)]
        eu[g + g * k + (git - 1) * k * k + (j - 1) * 
             k + i] = meanit[(g + k + g * k + d)]
        d = d + 1
      }
    }
  }
  dijg = 0
  for (git in 0:g) {
    for (i in 1:k) {
      for (j in 1:i) {
        dijg = dijg + 1
        for (im in 1:k) {
          for (gim in 0:g) {
            if (gim == 0) {
              if (git == 0) {
                vu[(g + g * k + g * k * k + im - 1) * 
                     NCOL + g + g * k + g * k * k + k + 
                     (i - 1) * k + j] = covit[(g + k + g * 
                                                 k + dijg - 1) * NCOLc + g + gim * k + 
                                                im]
                vu[(g + g * k + g * k * k + k + (i - 
                                                   1) * k + j - 1) * NCOL + g + g * k + 
                     g * k * k + im] = covit[(g + k + g * 
                                                k + dijg - 1) * NCOLc + g + gim * k + 
                                               im]
                vu[(g + g * k + g * k * k + im - 1) * 
                     NCOL + g + g * k + g * k * k + k + 
                     (j - 1) * k + i] = covit[(g + k + g * 
                                                 k + dijg - 1) * NCOLc + g + gim * k + 
                                                im]
                vu[(g + g * k + g * k * k + k + (j - 
                                                   1) * k + i - 1) * NCOL + g + g * k + 
                     g * k * k + im] = covit[(g + k + g * 
                                                k + dijg - 1) * NCOLc + g + gim * k + 
                                               im]
              }
              else {
                vu[(g + g * k + g * k * k + im - 1) * 
                     NCOL + g + g * k + (git - 1) * k * 
                     k + (i - 1) * k + j] = covit[(g + k + 
                                                     g * k + dijg - 1) * NCOLc + g + gim * 
                                                    k + im]
                vu[(g + g * k + (git - 1) * k * k + (i - 
                                                       1) * k + j - 1) * NCOL + g + g * k + 
                     g * k * k + im] = covit[(g + k + g * 
                                                k + dijg - 1) * NCOLc + g + gim * k + 
                                               im]
                vu[(g + g * k + g * k * k + im - 1) * 
                     NCOL + g + g * k + (git - 1) * k * 
                     k + (j - 1) * k + i] = covit[(g + k + 
                                                     g * k + dijg - 1) * NCOLc + g + gim * 
                                                    k + im]
                vu[(g + g * k + (git - 1) * k * k + (j - 
                                                       1) * k + i - 1) * NCOL + g + g * k + 
                     g * k * k + im] = covit[(g + k + g * 
                                                k + dijg - 1) * NCOLc + g + gim * k + 
                                               im]
              }
            }
            else {
              if (git == 0) {
                vu[(g + (gim - 1) * k + im - 1) * NCOL + 
                     g + g * k + g * k * k + k + (i - 1) * 
                     k + j] = covit[(g + k + g * k + dijg - 
                                       1) * NCOLc + g + gim * k + im]
                vu[(g + g * k + g * k * k + k + (i - 
                                                   1) * k + j - 1) * NCOL + g + (gim - 
                                                                                   1) * k + im] = covit[(g + k + g * k + 
                                                                                                           dijg - 1) * NCOLc + g + gim * k + im]
                vu[(g + (gim - 1) * k + im - 1) * NCOL + 
                     g + g * k + g * k * k + k + (j - 1) * 
                     k + i] = covit[(g + k + g * k + dijg - 
                                       1) * NCOLc + g + gim * k + im]
                vu[(g + g * k + g * k * k + k + (j - 
                                                   1) * k + i - 1) * NCOL + g + (gim - 
                                                                                   1) * k + im] = covit[(g + k + g * k + 
                                                                                                           dijg - 1) * NCOLc + g + gim * k + im]
              }
              else {
                vu[(g + (gim - 1) * k + im - 1) * NCOL + 
                     g + g * k + (git - 1) * k * k + (i - 
                                                        1) * k + j] = covit[(g + k + g * k + 
                                                                               dijg - 1) * NCOLc + g + gim * k + im]
                vu[(g + g * k + (git - 1) * k * k + (i - 
                                                       1) * k + j - 1) * NCOL + g + (gim - 
                                                                                       1) * k + im] = covit[(g + k + g * k + 
                                                                                                               dijg - 1) * NCOLc + g + gim * k + im]
                vu[(g + (gim - 1) * k + im - 1) * NCOL + 
                     g + g * k + (git - 1) * k * k + (j - 
                                                        1) * k + i] = covit[(g + k + g * k + 
                                                                               dijg - 1) * NCOLc + g + gim * k + im]
                vu[(g + g * k + (git - 1) * k * k + (j - 
                                                       1) * k + i - 1) * NCOL + g + (gim - 
                                                                                       1) * k + im] = covit[(g + k + g * k + 
                                                                                                               dijg - 1) * NCOLc + g + gim * k + im]
              }
            }
          }
        }
        for (gip in 1:g) {
          if (git == 0) {
            vu[(g + g * k + g * k * k + k + (i - 1) * 
                  k + j - 1) * NCOL + gip] = covit[(g + k + 
                                                      g * k + dijg - 1) * NCOLc + gip]
            vu[(gip - 1) * NCOL + g + g * k + g * k * 
                 k + k + (j - 1) * k + i] = covit[(g + k + 
                                                     g * k + dijg - 1) * NCOLc + gip]
            vu[(gip - 1) * NCOL + g + g * k + g * k * 
                 k + k + (i - 1) * k + j] = covit[(g + k + 
                                                     g * k + dijg - 1) * NCOLc + gip]
            vu[(g + g * k + g * k * k + k + (j - 1) * 
                  k + i - 1) * NCOL + gip] = covit[(g + k + 
                                                      g * k + dijg - 1) * NCOLc + gip]
          }
          else {
            vu[(g + g * k + (git - 1) * k * k + (i - 
                                                   1) * k + j - 1) * NCOL + gip] = covit[(g + 
                                                                                            k + g * k + dijg - 1) * NCOLc + gip]
            vu[(gip - 1) * NCOL + g + g * k + (git - 
                                                 1) * k * k + (j - 1) * k + i] = covit[(g + 
                                                                                          k + g * k + dijg - 1) * NCOLc + gip]
            vu[(gip - 1) * NCOL + g + g * k + (git - 
                                                 1) * k * k + (i - 1) * k + j] = covit[(g + 
                                                                                          k + g * k + dijg - 1) * NCOLc + gip]
            vu[(g + g * k + (git - 1) * k * k + (j - 
                                                   1) * k + i - 1) * NCOL + gip] = covit[(g + 
                                                                                            k + g * k + dijg - 1) * NCOLc + gip]
          }
        }
        odijg = 0
        for (ogit in 0:g) {
          for (oi in 1:k) {
            for (oj in 1:oi) {
              odijg = odijg + 1
              if (git == 0) {
                if (ogit == 0) {
                  vu[(g + g * k + g * k * k + k + (i - 
                                                     1) * k + j - 1) * NCOL + g + g * 
                       k + g * k * k + k + (oi - 1) * k + 
                       oj] = covit[(g + k + g * k + dijg - 
                                      1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (j - 
                                                     1) * k + i - 1) * NCOL + g + g * 
                       k + g * k * k + k + (oi - 1) * k + 
                       oj] = covit[(g + k + g * k + dijg - 
                                      1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oi - 
                                                     1) * k + oj - 1) * NCOL + g + g * 
                       k + g * k * k + k + (i - 1) * k + 
                       j] = covit[(g + k + g * k + dijg - 
                                     1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oi - 
                                                     1) * k + oj - 1) * NCOL + g + g * 
                       k + g * k * k + k + (j - 1) * k + 
                       i] = covit[(g + k + g * k + dijg - 
                                     1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (i - 
                                                     1) * k + j - 1) * NCOL + g + g * 
                       k + g * k * k + k + (oj - 1) * k + 
                       oi] = covit[(g + k + g * k + dijg - 
                                      1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (j - 
                                                     1) * k + i - 1) * NCOL + g + g * 
                       k + g * k * k + k + (oj - 1) * k + 
                       oi] = covit[(g + k + g * k + dijg - 
                                      1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oj - 
                                                     1) * k + oi - 1) * NCOL + g + g * 
                       k + g * k * k + k + (i - 1) * k + 
                       j] = covit[(g + k + g * k + dijg - 
                                     1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oj - 
                                                     1) * k + oi - 1) * NCOL + g + g * 
                       k + g * k * k + k + (j - 1) * k + 
                       i] = covit[(g + k + g * k + dijg - 
                                     1) * NCOLc + g + k + g * k + odijg]
                }
                else {
                  vu[(g + g * k + g * k * k + k + (i - 
                                                     1) * k + j - 1) * NCOL + g + g * 
                       k + (ogit - 1) * k * k + (oi - 1) * 
                       k + oj] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + g * k * k + k + (j - 
                                                     1) * k + i - 1) * NCOL + g + g * 
                       k + (ogit - 1) * k * k + (oi - 1) * 
                       k + oj] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oi - 1) * k + oj - 1) * NCOL + g + 
                       g * k + g * k * k + k + (i - 1) * 
                       k + j] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oi - 1) * k + oj - 1) * NCOL + g + 
                       g * k + g * k * k + k + (j - 1) * 
                       k + i] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (i - 
                                                     1) * k + j - 1) * NCOL + g + g * 
                       k + (ogit - 1) * k * k + (oj - 1) * 
                       k + oi] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + g * k * k + k + (j - 
                                                     1) * k + i - 1) * NCOL + g + g * 
                       k + (ogit - 1) * k * k + (oj - 1) * 
                       k + oi] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oj - 1) * k + oi - 1) * NCOL + g + 
                       g * k + g * k * k + k + (i - 1) * 
                       k + j] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oj - 1) * k + oi - 1) * NCOL + g + 
                       g * k + g * k * k + k + (j - 1) * 
                       k + i] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                }
              }
              else {
                if (ogit == 0) {
                  vu[(g + g * k + (git - 1) * k * k + 
                        (i - 1) * k + j - 1) * NCOL + g + 
                       g * k + g * k * k + k + (oi - 1) * 
                       k + oj] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (j - 1) * k + i - 1) * NCOL + g + 
                       g * k + g * k * k + k + (oi - 1) * 
                       k + oj] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + g * k * k + k + (oi - 
                                                     1) * k + oj - 1) * NCOL + g + g * 
                       k + (git - 1) * k * k + (i - 1) * 
                       k + j] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oi - 
                                                     1) * k + oj - 1) * NCOL + g + g * 
                       k + (git - 1) * k * k + (j - 1) * 
                       k + i] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (i - 1) * k + j - 1) * NCOL + g + 
                       g * k + g * k * k + k + (oj - 1) * 
                       k + oi] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (j - 1) * k + i - 1) * NCOL + g + 
                       g * k + g * k * k + k + (oj - 1) * 
                       k + oi] = covit[(g + k + g * k + 
                                          dijg - 1) * NCOLc + g + k + g * k + 
                                         odijg]
                  vu[(g + g * k + g * k * k + k + (oj - 
                                                     1) * k + oi - 1) * NCOL + g + g * 
                       k + (git - 1) * k * k + (i - 1) * 
                       k + j] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                  vu[(g + g * k + g * k * k + k + (oj - 
                                                     1) * k + oi - 1) * NCOL + g + g * 
                       k + (git - 1) * k * k + (j - 1) * 
                       k + i] = covit[(g + k + g * k + dijg - 
                                         1) * NCOLc + g + k + g * k + odijg]
                }
                else {
                  vu[(g + g * k + (git - 1) * k * k + 
                        (i - 1) * k + j - 1) * NCOL + g + 
                       g * k + (ogit - 1) * k * k + (oi - 
                                                       1) * k + oj] = covit[(g + k + g * 
                                                                               k + dijg - 1) * NCOLc + g + k + g * 
                                                                              k + odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (j - 1) * k + i - 1) * NCOL + g + 
                       g * k + (ogit - 1) * k * k + (oi - 
                                                       1) * k + oj] = covit[(g + k + g * 
                                                                               k + dijg - 1) * NCOLc + g + k + g * 
                                                                              k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oi - 1) * k + oj - 1) * NCOL + g + 
                       g * k + (git - 1) * k * k + (i - 
                                                      1) * k + j] = covit[(g + k + g * 
                                                                             k + dijg - 1) * NCOLc + g + k + g * 
                                                                            k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oi - 1) * k + oj - 1) * NCOL + g + 
                       g * k + (git - 1) * k * k + (j - 
                                                      1) * k + i] = covit[(g + k + g * 
                                                                             k + dijg - 1) * NCOLc + g + k + g * 
                                                                            k + odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (i - 1) * k + j - 1) * NCOL + g + 
                       g * k + (ogit - 1) * k * k + (oj - 
                                                       1) * k + oi] = covit[(g + k + g * 
                                                                               k + dijg - 1) * NCOLc + g + k + g * 
                                                                              k + odijg]
                  vu[(g + g * k + (git - 1) * k * k + 
                        (j - 1) * k + i - 1) * NCOL + g + 
                       g * k + (ogit - 1) * k * k + (oj - 
                                                       1) * k + oi] = covit[(g + k + g * 
                                                                               k + dijg - 1) * NCOLc + g + k + g * 
                                                                              k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oj - 1) * k + oi - 1) * NCOL + g + 
                       g * k + (git - 1) * k * k + (i - 
                                                      1) * k + j] = covit[(g + k + g * 
                                                                             k + dijg - 1) * NCOLc + g + k + g * 
                                                                            k + odijg]
                  vu[(g + g * k + (ogit - 1) * k * k + 
                        (oj - 1) * k + oi - 1) * NCOL + g + 
                       g * k + (git - 1) * k * k + (j - 
                                                      1) * k + i] = covit[(g + k + g * 
                                                                             k + dijg - 1) * NCOLc + g + k + g * 
                                                                            k + odijg]
                }
              }
            }
          }
        }
      }
    }
  }
  eu <- as.matrix(eu, nrow = g + g * k + g * k * k + k + k * 
                    k, ncol = 1)
  vu <- as.matrix(vu, nrow = g + g * k + g * k * k + k + k * 
                    k, ncol = g + g * k + g * k * k + k + k * k)
  dim(vu) <- c(g + g * k + g * k * k + k + k * k, g + g * k + 
                 g * k * k + k + k * k)
  dim(meanit) <- dmeanit
  dim(covit) <- dcovit
  eu[1:g, 1] = meanit[1:g, 1]
  eu[(g + g * k + g * k * k + 1):(g + g * k + g * k * k + k), 
     1] = meanit[(g + 1):(g + k), 1]
  eu[(g + 1):(g + g * k), 1] = meanit[(g + k + 1):(g + k + 
                                                     g * k), 1]
  vu[1:(g + g * k), 1:(g + g * k)] = covit[c((1:g), ((g + k + 
                                                        1):(g + k + g * k))), c((1:g), ((g + k + 1):(g + k + 
                                                                                                       g * k)))]
  vu[(g + g * k + g * k * k + 1):(g + g * k + g * k * k + k), 
     (g + g * k + g * k * k + 1):(g + g * k + g * k * k + 
                                    k)] = covit[(g + 1):(g + k), (g + 1):(g + k)]
  vu[(g + g * k + g * k * k + 1):(g + g * k + g * k * k + k), 
     1:(g + g * k)] = covit[(g + 1):(g + k), c((1:g), ((g + 
                                                          k + 1):(g + k + g * k)))]
  vu[1:(g + g * k), (g + g * k + g * k * k + 1):(g + g * k + 
                                                   g * k * k + k)] = covit[c(1:g, (g + k + 1):(g + k + g * 
                                                                                                 k)), (g + 1):(g + k)]
  return(list(eu = eu, vu = vu, gprop = eu[1:g, 1], g = g, 
              k = k))
}