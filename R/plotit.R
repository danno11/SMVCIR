plotit<-function (x, dimensions = c(1, 2), groups = 1:10, GL = FALSE,
          ID = FALSE, build_svm = FALSE, kernel = NULL, svmModel, method = "SMVCIR")
{
  if(class(x)=="smvcir"){
    groupnum<-x$groups
    groupcol<-ncol(x$direct)
    compcases<-x$compcases
    datmat<-x<-x$direct      ####add x double assignment

  } else{
    groupcol<-NA
    for (i in 1:ncol(x)){
      if (class(x[,i])=="factor"){
        groupcol[i]<-i
      } else{
        groupcol[i]<-NA
      }
    }
    groupcol<-groupcol[!is.na(groupcol)]
    if(length(groupcol)>1)
      stop("Only numerical predictors allowed for smvcir, appears you have more than one variable of class=='factor'")
        if(length(groupcol)==0)
          stop("Missing a dependent group variable of class=='factor' in your coordinates")

    groupnum<-length(unique(x[,groupcol]))
    datmat<-x
    #datmat[,groupcol]<-as.numeric(datmat[,groupcol])
    compcases <- data.frame(complete.cases(datmat)) ###look ahead check line 72
    compcases <- cbind(compcases, t(t(1:nrow(datmat))))
  }



  if (length(groups) >= groupnum) {
    bx <- matrix("Rendering all Groups", 1, 1)
    rownames(bx) <- ""
    colnames(bx) <- ""
    print(noquote(bx))
    groups <- 1:(groupnum)
  }
  colors <- c("black", "red", "blue", "orange", "purple", "brown",
              "green", "pink", "yellow", "aquamarine")
  gsymbols <- c(0, 15, 2, 17, 3, 4, 8, 1, 16, 6)
  gsymbolsExpl <- c("Empty Square", "Fill Square", "Empty Triangle",
                    "Fill Triangle", "Plus", "X", "Asterisk", "Empty Diamond",
                    "Fill Diamond", "Upside Down Triangle")
  dispmat = c()
  for (i in 1:length(groups)) {
    if ((groups[i] > groupnum) || (groups[i] < 1)) {
      print(noquote("User entered illegal group"))
      return(NaN)
    }
    else {
      dispmat <- rbind(dispmat, c(levels(x[,groupcol])[i], colors[groups[i]],  #####change first item of c()
                                  gsymbolsExpl[groups[i]]))
    }
  }
  rownames(dispmat) <- rep("", nrow(dispmat))
  colnames(dispmat) <- c("Group", "Color", "Symbol")
  print(noquote(dispmat))
  if (GL == TRUE) {
    if (length(dimensions) < 3) {
      dimensions = c(1, 2, 3)
    }
    else if (length(dimensions) != 3) {
      print(noquote("Warning: More than 3 dimensions specified"))
      dimensions = c(dimensions[1], dimensions[2], dimensions[3])
    }
  }
  for (i in 1:length(dimensions)) {
    if ((dimensions[i] < 1)) {
      print(noquote("User entered illegal dimension"))
      print(i)
      return(NaN)
    }
  }
  if (length(dimensions) == 1) {
    plotit1D(wmat = datmat, dimension = dimensions, groups = groups)
  }
  else {
    if (GL == TRUE) {
      plotit3D(wmat = datmat, coords = dimensions, groups = groups,
               ID = ID, compcases = compcases, build_svm = build_svm, groupcol = groupcol, kernel=kernel, svmModel = svmModel)
    }
    else {
      plotit2D(wmat = datmat, dimensions = dimensions,
               groups = groups, method = method, build_svm = build_svm, kernel = kernel, svmModel = svmModel)
    }
  }
}
