plotit2D<-function (wmat, dimensions, groups, method, build_svm, kernel)
{
  if(build_svm==TRUE){
    if(!require("e1071")){
      stop("The package e1071 must be installed.")
    }
    if(!require("gplots")){
      stop("The package gplots must be installed.")
    }
    if(is.null(kernel)){
      kernel<-"radial"
    }

     k = length(wmat[1, ]) - 1
     mod_obj<-svm(as.formula(paste(names(wmat)[k+1],"~ .")), kernel = kernel, data = wmat)
     slice<-lapply(wmat[,1:k], median)

  }
  def.par <- par(no.readonly = TRUE)
  par(mai = c(0, 0, 0, 0))
  totcol <- length(dimensions)
  mat <- matrix(nrow = totcol, ncol = totcol)
  ordering <- 1
  for (i in 1:totcol) {
    for (j in 1:totcol) {
      mat[i, j] <- ordering
      ordering <- ordering + 1
    }
  }
  layout(mat, widths = rep(0.5, ncol(mat)), heights = rep(0.5,
                                                          ncol(mat)), respect = TRUE)
  colors <- c("black", "red", "blue", "orange", "purple", "brown",
              "green", "pink", "yellow", "aquamarine")
  gsymbols <- c(0, 15, 2, 17, 3, 4, 8, 1, 16, 6)
  for (i in 1:totcol) {
    for (j in 1:totcol) {
      if (i == j) {
        if (i != totcol) {
          plot(wmat[, dimensions[j]], wmat[, dimensions[j]],
               type = "n", yaxt = "n", xlab = "", ylab = "")
          legend("center", legend = c(dimensions[j]),
                 bty = "n", title = method)
          axis(side = 4)
        }
        else {
          plot(wmat[, dimensions[j]], wmat[, dimensions[j]],
               type = "n", yaxt = "n", xaxt = "n", xlab = "",
               ylab = "")
          legend("center", legend = c(j), bty = "n",
                 title = method)
          axis(side = 3)
        }
      }
      else if (i > j) {
          x <- c(0, 1)
          plot(x, type = "n", axes = FALSE, xlab = "",
               ylab = "")
          legend("center", legend = c(" "), bty = "n",
                 title = "")

      }
      else {
        coords = dimensions
        k = length(wmat[1, ]) - 1
        n = length(wmat[, 1])
        zmat = wmat[, 1:k]
        indi = wmat[, (k + 1)]
        numg = max(as.numeric(indi)) ###add as.numeric()
        ssgroup = rep(0, numg)
        for (p in 1:numg) {
          see = wmat[indi == levels(indi)[p], ]
          ssgroup[p] = length(see[, 1])
        }
        xl = c(min(zmat), max(zmat))
        yl = c(min(zmat), max(zmat))
        newindi <- order(indi)
        z1 <- zmat[, coords[i]]
        z2 <- zmat[, coords[j]]
        zz1 <- z1[newindi]
        zz2 <- z2[newindi]
        plot(zz2, zz1, xlab = "", ylab = "", pch = " ",
             xlim = xl, ylim = yl, axes = FALSE)
        if(build_svm==TRUE){
          xr <- seq(xl[1], xl[2], length = 50)
          yr <- seq(yl[1], yl[2], length = 50)
          formala<-as.formula(paste(names(zmat)[i], "~", names(zmat)[j]))
          lis <- c(list(yr), list(xr), slice)
          names(lis)[1:2] <- colnames(wmat)[c(i,j)]
          new <- expand.grid(lis)[, labels(terms(mod_obj))]
          preds <- predict(mod_obj, new)
          .filled.contour(xr,
                          yr,
                          matrix(as.numeric(preds),nrow = length(xr), byrow = TRUE),
                          levels = 1:length(levels(preds))+1,
                          col=colorpanel(3, low = "gray70", high = "gray45"))
        }



        ic <- 0
        coli <- 0
        for (p in 1:numg) {
          there <- 0
          for (m in 1:length(groups)) {
            if (groups[m] == p) {
              there <- 1
              break
            }
          }
          start <- ic + 1
          stop <- ic + ssgroup[p]
          coli <- coli + 1
          if (there == 1) {
            points(zz2[start:stop], zz1[start:stop],
                   col = colors[p], pch = gsymbols[p])
            abline(lsfit(zz2[start:stop], zz1[start:stop]),
                   col = colors[p])
            there <- 0
          }
          ic <- ic + ssgroup[p]
        }
      }
    }
  }
  par(def.par)
}
