plotit3D<-function (wmat, coords, groups, ID, compcases, build_svm, groupcol, kernel, svmModel)
{
  if(!require("rgl")){
    stop("The package rgl must be installed.")
  }

  k = length(wmat[1, ]) - 1
  n = length(wmat[, 1])
  zmat = wmat[, 1:k]
  indi = wmat[, as.numeric(k + 1)]     ####add as.numeric model matrix??
  a = 1 * (compcases[, 1])
  oorder <- compcases[a == 1, 2]
  numg = length(unique(as.numeric(indi)))
  ssgroup = rep(0, numg)
  for (i in 1:numg) {
    see = wmat[indi==levels(indi)[i], ]###levels around ind
    ssgroup[i] = length(see[, 1])
  }
  xl = c(min(zmat), max(zmat))
  yl = c(min(zmat), max(zmat))
  newindi <- order(indi)
  z1 <- zmat[, coords[1]]
  z2 <- zmat[, coords[2]]
  z3 <- zmat[, coords[3]]
  zz1 <- z1[newindi]
  zz2 <- z2[newindi]
  zz3 <- z3[newindi]
  zoorder <- oorder[newindi]
  xl = c(min(zz1), max(zz1)) ###changes zz2->zz1 here
  yl = c(min(zz2), max(zz2))
  zl = c(min(zz3), max(zz3))
  colors <- c("black", "red", "blue", "orange", "purple", "brown",
              "green", "pink", "yellow", "aquamarine")
  open3d()
  plot3d(x = zz1, y = zz2, z = zz3, xlab = "", ylab = "", zlab = "",   ###and here
         type = "n")
  decorate3d(xl, yl, zl, xlab = paste("D", coords[1], sep = ""),
             ylab = paste("D", coords[2], sep = ""), zlab = paste("D",
                                                                  coords[3], sep = ""))
  ic <- 0
  coli <- 0
  for (i in 1:numg) {
    there <- 0
    for (j in 1:length(groups)) {
      if (groups[j] == i) {
        there <- 1
        break
      }
    }
    start <- ic + 1
    stop <- ic + ssgroup[i]
    if (there == 1) {
      points3d(zz1[start:stop], zz2[start:stop], zz3[start:stop],
               color = colors[i], size = 7)
      there <- 0
    }
    ic <- ic + ssgroup[i]
  }
  x <- zz1
  y <- zz2
  z <- zz3
  kept <- rep(0, n)###change length(x) to length(n)
  dispdat <- data.frame(zoorder, x, y, z, kept)
  names(dispdat) <- c("Obs #", paste("D", coords[1], sep = ""),
                      paste("D", coords[2], sep = ""), paste("D", coords[3],
                                                             sep = ""), "keep")
  while (ID) {
    f <- select3d(button = c("right"))
    dispdat[, 5] <- f(dispdat[, 2], dispdat[, 3], dispdat[,
                                                          4])
    rgl.clear()
    decorate3d(xl, yl, zl, xlab = paste("D", coords[1], sep = ""),
               ylab = paste("D", coords[2], sep = ""), zlab = paste("D",
                                                                    coords[3], sep = ""))
    ic <- 0
    coli <- 0
    for (i in 1:numg) {
      there <- 0
      for (j in 1:length(groups)) {
        if (groups[j] == i) {
          there <- 1
          break
        }
      }
      start <- ic + 1
      stop <- ic + ssgroup[i]
      if (there == 1) {
        dispdatred <- dispdat[start:stop, ]
        dispdatredKeep <- subset(dispdatred, keep ==
                                   1)
        if (nrow(dispdatredKeep) > 0) {
          pdKeep <- as.matrix(dispdatredKeep[, c(1, 2,
                                                 3, 4)])
          rownames(pdKeep) <- rep("", nrow(pdKeep))
          print(pdKeep)
        }
        dispdatredNKeep <- subset(dispdatred, keep ==
                                    0)
        points3d(dispdatredNKeep[, 2], dispdatredNKeep[,
                                                       3], dispdatredNKeep[, 4], color = colors[i],
                 size = 7)
        points3d(dispdatredKeep[, 2], dispdatredKeep[,
                                                     3], dispdatredKeep[, 4], color = "gray", size = 7)
        there <- 0
      }
      ic <- ic + ssgroup[i]
    }
  }

  if(build_svm==TRUE){
    if(!require("e1071")){
      stop("The package e1071 needs to be installed.")
    }
    if(!require("misc3d")){
      stop("The package misc3d must be installed.")
    }
    if(numg > 2){
      stop("svm can only be built here for binary classification.")
    }
    if(is.null(kernel)){
      kernel<-"radial"
    }

    svmx<-wmat[,dimensions]
    svmy<-wmat[,groupcol]

    if(!is.null(svmModel)){
      svm_mod<-svmModel
    } else{
      svm_mod<-svm(x = svmx, y = svmy, kernel=kernel, type = "C-classification", scale = FALSE)
    }
    nnew = 50
    newdat.list = lapply(svmx, function(svmx) seq(min(svmx), max(svmx), len=nnew))
    newdat      = expand.grid(newdat.list)
    newdat.pred = predict(svm_mod, newdata=newdat, decision.values=T)
    newdat.dv   = attr(newdat.pred, 'decision.values')
    newdat.dv   = array(newdat.dv, dim=rep(nnew, 3))
    # Fit/plot an isosurface to the decision boundary
    contour3d(newdat.dv, level=0, x=newdat.list[[coords[1]]], y=newdat.list[[coords[2]]], z=newdat.list[[coords[3]]], add=T)
  }
}
