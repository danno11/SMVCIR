#'Main SMVCIR function
#'
#'Build a Sliced Mean Variance Covariance Inverse Regression model.
#'
#'@param group A character string specifying the name of the class variable of interest in your dataset.
#'@param data  A data frame (including your group variable).
#'@param pdt Percentage of SMVCIR discrimination desired from dimension
#'@param level Level of dimensionality test
#'@param test Types of tests, can use either, neither, or both.
#'@param empsimss Number of draws to use in performing each dimensionality test.
#'@param scree_plot If TRUE, a scree plot of cumulative percentage of variation explained by discriminant dimensions is produced.
#'@export


smvcir<-function (group, data, pdt = 100, level = 0.05, test = FALSE, scree_plot=FALSE)
{
  empTest = FALSE
  apempTest = FALSE
  if (test == TRUE) {
    apempTest = TRUE
  }
  compcases <- data.frame(complete.cases(data))
  compcases <- cbind(compcases, t(t(1:nrow(data))))
  data_nu <- na.omit(data)
  namelist <- attr(data, "names")
  groupindex <- 0
  for (i in 1:length(namelist)) {
    if (namelist[i] == group) {
      groupindex <- i
      break
    }
  }
  if (groupindex == 0) {
    return(NAN)
  }
  if (groupindex == dim(data_nu)[2]) {
    stage1readydat <- data.matrix(data_nu)
  } else if (groupindex == 1) {
    stage1readydat <- data.matrix(data_nu[, c(2:dim(data_nu)[2],
                                              1)])
  } else {
    stage1readydat <- data.matrix(data_nu[, c(1:(groupindex -
                                                   1), (groupindex + 1):dim(data_nu)[2], groupindex)])
  }
  rm(data_nu)
  n <- dim(stage1readydat)[1]
  k <- dim(stage1readydat)[2] - 1
  g <- max(stage1readydat[, k + 1])
  class_labels<-levels(data[,groupindex]) #####class_labels
  standat <- cbind(scale(stage1readydat[, 1:k], center = TRUE,
                         scale = TRUE), stage1readydat[, k + 1])
  xbar <- t(t(colMeans(stage1readydat[, 1:k])))
  sigmaminusoh <- matrix(0, nrow = k, ncol = k)
  diag(sigmaminusoh) <- t(apply(stage1readydat[, 1:k], 2, sd))^(-1)
  ind <- stage1readydat[, (k + 1)]
  nupart = matrix(rep(0, k), ncol = 1)
  deltapart = matrix(rep(0, k), ncol = 1)
  storesigzi = matrix(rep(0, k), ncol = 1)
  storecapdel = matrix(rep(0, k), ncol = 1)
  storecapdel0 = matrix(rep(0, k), ncol = 1)
  storevecd = matrix(rep(0, k), ncol = 1)
  ni = rep(0, g)
  for (i in 1:g) {
    xg = stage1readydat[ind == i, 1:k]
    ni[i] = length(xg[, 1])
    xbari = t(t(colMeans(xg)))
    zbari = sigmaminusoh %*% (xbari - xbar)
    sigzi <- sigmaminusoh %*% var(xg) %*% sigmaminusoh
    storesigzi = cbind(storesigzi, sigzi)
    vecd = diag(sigzi)
    storevecd = cbind(storevecd, vecd)
    nui = sqrt(ni[i]/n) * zbari
    nupart = cbind(nupart, nui)
  }
  nupart = nupart[, 2:(g + 1)]
  storevecd = storevecd[, 2:(g + 1)]
  storesigzi = storesigzi[, 2:(k * g + 1)]
  sigzbar = matrix(rep(0, k^2), ncol = k)
  ib = 1
  ie = k
  for (i in 1:g) {
    sigzbar = sigzbar + (ni[i]/n) * storesigzi[, ib:ie]
    ib = ib + k
    ie = ie + k
  }
  sigzbard = diag(sigzbar)
  ib = 1
  ie = k
  for (i in 1:g) {
    deli = sqrt(ni[i]/n) * (storevecd[, i] - sigzbard)
    deltapart = cbind(deltapart, deli)
    capdeli = sqrt(ni[i]/n) * (storesigzi[, ib:ie] - sigzbar)
    capdeli0 = capdeli - diag(diag(capdeli))
    storecapdel <- cbind(storecapdel, capdeli)
    storecapdel0 <- cbind(storecapdel0, capdeli0)
    ib = ib + k
    ie = ie + k
  }
  deltapart <- deltapart[, 2:(g + 1)]
  storecapdel <- storecapdel[, 2:(k * g + 1)]
  storecapdel0 <- storecapdel0[, 2:(k * g + 1)]
  spansetF = cbind(storecapdel0, deltapart, nupart)
  kernelF = spansetF %*% t(spansetF)
  if (empTest == TRUE | apempTest == TRUE) {
    estclt <- feststdclt1(g, stage1readydat)
    estgrpmom <- eststddelta1(list(eu = estclt$eu, gprop = estclt$gprop,
                                   g = estclt$g, k = estclt$k))
    gprop = estclt$gprop
    estvar <- eststddelta2(list(ef = estgrpmom$ef, g = g,
                                k = k, gprop = gprop))
    eststd <- eststddelta3(list(ef = estvar$ef, g = g, k = k,
                                gprop = gprop))
    estcentvar <- eststddelta4(list(ef = eststd$ef, g = g,
                                    k = k, gprop = gprop))
    estpropwt <- eststddelta5(list(ef = estcentvar$ef, g = g,
                                   k = k, gprop = gprop))
    eststackvarint <- eststackvar(list(ef = estpropwt$ef,
                                       g = g, k = k, gprop = gprop))
    vals <- eigen(kernelF, symmetric = TRUE)$values
    or <- rev(order(abs(vals)))
    evalues <- vals[or]
    esteigensmvcir1 <- evalues
  }
  tmp = svd(spansetF)
  sumsing = rep(0, k)
  for (j in 1:k) {
    sumsing[j] = sum(tmp$d[1:j])
  }
  sumsing = (100/sumsing[k]) * sumsing
  scree = cbind(1:k, sumsing)
  dpDIM = -1
  for (j in 1:k) {
    if (sumsing[j] >= pdt & dpDIM == -1) {
      dpDIM = j
    }
  }
  if (dpDIM == -1) {
    dpDIM = k
  }
  if (scree_plot == TRUE) {
    x11()
    plot(c(0, scree[, 1]), c(0, scree[, 2]), xlab = "Singular value",
         ylab = "Percentage", ylim = c(0, 100), type = "l",
         main = "Scree Plot SV %")
    abline(h = pdt)
  }
  emppvals <- matrix(NA, nrow = k, ncol = 1)
  apemppvals <- matrix(NA, nrow = k, ncol = 1)
  if (empTest == TRUE | apempTest == TRUE) {
    vf <- matrix(0, nrow = g * k + g * k * k + g * k, ncol = g *
                   k + g * k * k + g * k)
    vf[(1:(g * k + g * k * k)), (1:(g * k + g * k * k))] <- (estpropwt$d5 %*%
                                                               estcentvar$d4 %*% eststd$d3 %*% estvar$d2 %*% estgrpmom$d1 %*%
                                                               estclt$vu %*% t(estgrpmom$d1) %*% t(estvar$d2) %*%
                                                               t(eststd$d3) %*% t(estcentvar$d4) %*% t(estpropwt$d5))
    for (i in 1:k) {
      tstat <- n * sum(esteigensmvcir1[(i):length(esteigensmvcir1)])
      tsmvcird <- i - 1
      bigK <- kunstack(k, as.matrix(eststackvarint$ef))
      sv <- svd(bigK, nu = nrow(bigK), nv = ncol(bigK))
      gammanot <- sv$u[, (tsmvcird + 1):(ncol(sv$u))]
      psinot <- sv$v[, (tsmvcird + 1):ncol(sv$v)]
      dcmat <- ((t(psinot) %x% t(gammanot)) %*% eststackvarint$Perm %*%
                  vf %*% t(eststackvarint$Perm) %*% (psinot %x%
                                                       gammanot))
      if (empTest == TRUE) {
        devalues <- eigen(dcmat, symmetric = TRUE)$values
        exceed = 0
        for (j in 1:empsimss) {
          realz <- t(devalues) %*% t(t(rchisq(length(devalues),
                                              1)))
          if (realz >= tstat) {
            exceed <- exceed + 1
          }
        }
        pvalue <- exceed/empsimss
        emppvals[i, 1] <- pvalue
      }
      if (apempTest == TRUE) {
        trdcmat <- sum(eigen(dcmat, symmetric = TRUE)$values)
        trdcmat2 <- sum(eigen(dcmat %*% t(dcmat), symmetric = TRUE)$values)
        d_num <- round((trdcmat^2)/trdcmat2)
        scalecorrectstat <- tstat * ((trdcmat/d_num)^(-1))
        pvalue <- 1 - pchisq(scalecorrectstat, d_num)
        apemppvals[i, 1] <- pvalue
      }
    }
  }
  tmp.eig = eigen(kernelF)
  wmati = standat[, 1:k] %*% tmp.eig$vectors
  wmati = cbind(wmati, standat[, k + 1])
  stdcoeffmat <- matrix(NA, nrow = k, ncol = k)
  x <- cbind(scale(stage1readydat[, 1:k], center = TRUE, scale = TRUE),
             rep(1, nrow(stage1readydat[, 1:k])))
  for (i in 1:k) {
    stdcoeffmat[1:k, i] <- t(t(lm.fit(y = wmati[, i], x = x)$coefficients[1:k]))[1:k,
                                                                                 1]
    norm <- sum(stdcoeffmat[1:k, i] * stdcoeffmat[1:k, i])
    stdcoeffmat[1:k, i] <- stdcoeffmat[1:k, i]/sqrt(norm)
  }
  if (empTest == TRUE | apempTest == TRUE) {
    if (empTest == TRUE) {
      empDIM = -1
      for (i in 1:k) {
        if (emppvals[i] >= level & empDIM == -1) {
          empDIM = i - 1
        }
      }
      if (empDIM == -1) {
        empDIM = k
      }
    }
    else {
      empDIM = -1
    }
    if (apempTest == TRUE) {
      apempDIM = -1
      for (i in 1:k) {
        if (apemppvals[i] >= level & apempDIM == -1) {
          apempDIM = i - 1
        }
      }
      if (apempDIM == -1) {
        apempDIM = k
      }
    }
    else {
      apempDIM = -1
    }
    chosenDIM = max(empDIM, apempDIM)
  }
  printit1 <- matrix("", 7, 1)
  printit1[1, 1] <- "SMVCIR"
  printit1[2, 1] <- paste("# Groups:                ", g, sep = "")
  printit1[3, 1] <- paste("# Predictors:            ", k, sep = "")
  printit1[4, 1] <- paste("Observations used:     ", nrow(standat),
                         sep = "")
  printit1[5, 1] <- paste("Total Observations:    ", nrow(data),
                         sep = "")
  if (empTest == TRUE | apempTest == TRUE) {
    printit1[6, 1] <- paste("Dimension: ", chosenDIM, ", at level: ",
                           level, sep = "")
  }
  printit1[7, 1] <- paste("Dimension: ", dpDIM, ", provides ",
                         round(sumsing[dpDIM], 0), "% Discrimination", sep = "")
  rownames(printit1) <- rep("", 7)
  colnames(printit1) <- rep("", 1)
  print(printit1, quote = FALSE)  ##First output for summary function
  if (empTest == TRUE | apempTest == TRUE) {
    printit2 <- matrix("", 1, 1)
    printit2[1, 1] <- "Dimensionality Test P-Values"
    rownames(printit2) <- c("")
    colnames(printit2) <- c("")
    print(printit2, quote = FALSE)
    if (apempTest == TRUE & empTest == TRUE) {
      pvalmat <- round(cbind(emppvals, apemppvals), 3)
      colnames(pvalmat) <- c("Empirical", "Approximate Empirical")
    }
    else if (apempTest == TRUE) {
      pvalmat <- round(t(t(apemppvals)), 3)
      colnames(pvalmat) <- c("Approximate Empirical")
    }
    else {
      pvalmat <- round(t(t(emppvals)), 3)
      colnames(pvalmat) <- c("Empirical")
    }
    rownamesIt <- 0:(k - 1)
    if (empDIM > -1) {
      if (empDIM < k) {
        rownamesIt[empDIM + 1] <- paste("E ", rownamesIt[empDIM +
                                                           1], sep = "")
      }
    }
    if (apempDIM > -1) {
      if (apempDIM < k) {
        rownamesIt[apempDIM + 1] <- paste("AE ", rownamesIt[apempDIM +
                                                              1], sep = "")
      }
    }
    rownames(pvalmat) <- rownamesIt
    print(pvalmat, quote = FALSE)
  }
  prednames <- rep("", k)
  i = 1
  for (j in 1:(k + 1)) {
    if (names(data)[j] != group) {
      prednames[i] <- names(data)[j]
      i = i + 1
    }
  }
  rownames(stdcoeffmat) <- prednames
  colnames(stdcoeffmat) <- paste("D", 1:k, sep = "")
  #printit <- matrix("", 1, 1)
  #printit[1, 1] <- "Standardized Coefficients"
  #rownames(printit) <- c("")
  #colnames(printit) <- c("")
  #print(printit, quote = FALSE)
  # if (empTest == TRUE | apempTest == TRUE) {
  #   print(round(stdcoeffmat[, 1:max(min(chosenDIM, dpDIM),
  #                                  1)], 3), quote = FALSE)
  #}
  #else {
  #  print(round(stdcoeffmat[, 1:(max(dpDIM, 1))], 3), quote = FALSE)
  #}
  colnames(stage1readydat) <- c(prednames, group)
  colnames(standat) <- c(prednames, group)
  transdat <- data.frame(wmati)
  names(transdat) <- c(paste("D", rep(1:k), sep = ""), group)
  spanFnames <- character(length = ncol(spansetF))
  nb <- ncol(spansetF)
  for (p in 1:nb) {
    if (p <= g * k) {
      if (floor(p/k) < p/k) {
        tempgroup <- floor(p/k) + 1
        covcolm <- p - floor(p/k) * k
      }
      else {
        tempgroup <- floor(p/k)
        covcolm = k
      }
      spanFnames[p] <- paste("C", tempgroup, ".", covcolm,
                             sep = "")
    }
    else if (p <= g * k + g) {
      tempgroup <- p - g * k
      spanFnames[p] <- paste("V", tempgroup, sep = "")
    }
    else {
      tempgroup <- p - g * k - g
      spanFnames[p] <- paste("M", tempgroup, sep = "")
    }
  }
  colnames(spansetF) <- spanFnames
  if (empTest == TRUE | apempTest == TRUE) {
    if (min(chosenDIM, dpDIM) > 1) {
      TRANScm <- cor(transdat[, 1:(max(min(chosenDIM, dpDIM),
                                       1))])
      maxit <- max(TRANScm - diag(diag(TRANScm)))
      minit <- min(TRANScm - diag(diag(TRANScm)))
      if (abs(minit) > maxit) {
        maxTRANScm <- minit
      }
      else {
        maxTRANScm <- maxit
      }
    }
    else {
      TRANScm <- 0
      maxTRANScm <- 0
    }
  }
  else {
    if (dpDIM > 1) {
      TRANScm <- cor(transdat[, 1:dpDIM])
      maxit <- max(TRANScm - diag(diag(TRANScm)))
      minit <- min(TRANScm - diag(diag(TRANScm)))
      if (abs(minit) > maxit) {
        maxTRANScm <- minit
      }
      else {
        maxTRANScm <- maxit
      }
    }
    else {
      TRANScm <- 0
      maxTRANScm <- 0
    }
  }
  printit <- matrix("SMVCIR dimensions should have low correlations.",
                    1, 1)
  rownames(printit) <- c("")
  colnames(printit) <- c("")
  print(printit, quote = FALSE)
  printit <- matrix(paste("Maximum SMVCIR dimension correlation: ",
                          maxTRANScm, sep = ""), 1, 1)
  rownames(printit) <- c("")
  colnames(printit) <- c("")
  print(printit, quote = FALSE)
  printit <- matrix("SMVCIR correlations.", 1, 1)
  rownames(printit) <- c("")
  colnames(printit) <- c("")
  #print(printit, quote = FALSE)
  #print(TRANScm, quote = FALSE) remove printing of discrim coords (too large)
  if (empTest == FALSE & apempTest == FALSE) {
    chosenDIM = NA
  }
  c_mat<-tmp.eig$vectors
  muhat_ls<-aggregate(scale(stage1readydat[,1:k]), list(stage1readydat[,k+1]), mean)
  muhat_ls<-matrix(NA, nrow = g, ncol = k)
  for(i in 1:g){
    muhat_ls[i,1:k] <- colMeans(standat[ind==i,1:k])
  }
  rownames(muhat_ls)<-levels(data[,which(names(data)%in%group)])
  colnames(muhat_ls)<-prednames
  muhat_z<-muhat_ls%*%c_mat
  rownames(muhat_z)<-rownames(muhat_ls)
  colnames(muhat_z)<-paste("D", 1:k, sep = "")
  sighatx<-cov(standat[,1:k])
  sighatz<-t(c_mat)%*%sighatx%*%c_mat
  class.props<-matrix(NA, ncol = g)
  for(i in 1:g){
    class.props[1,i]<-mean(standat[,k+1]==i)
  }
 colnames(class.props)<-class_labels
 colnames(sighatz)<-rownames(sighatz)<-c(paste("D", rep(1:k), sep = ""))

  #if (plot == TRUE) {
  #  if (GL) {
  #   plot(smv, dimensions = 1:3, GL = GL)
  # }
  # else {
  #   if (empTest == TRUE | apempTest == TRUE) {
  #     plot(smv, dimensions = 1:max(min(chosenDIM, dpDIM),
  #                                  1), GL = GL)
  #   }
  #   else {
  #     plot(smv, dimensions = 1:max(dpDIM, 1), GL = GL)
  #   }
  # }
  #}#
  transdat[,k+1]<-factor(class_labels[transdat[,k+1]])

  smv <- list(groups = g, predictors = k, statdim = chosenDIM,
              sighatz = sighatz, muhat_z = muhat_z, groupindex = groupindex,
              class.props=class.props,
              muhat_ls = muhat_ls, xbar = xbar, sighatx = sighatx,
              dimCorr = TRANScm, maxTRANScm = maxTRANScm,
              direct = transdat, compcases = compcases, spansetF = spansetF,
              call = match.call(), coefficients = stdcoeffmat, originalx = stage1readydat[,1:k],
              kernelF = kernel)
  attr(smv, "class") <- "smvcir"
  return(smv)
}
