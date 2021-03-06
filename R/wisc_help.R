#'Wisconsin Breast Cancer dataset
#'
#'@docType data
#'
#'@usage data(wisc)
#'
#'@keywords datasets
#'
#'@source \href{http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Original)}{UCI Machine Learning Repository}
#'
#'
#'Class Variable Tumor Type "tumor": "malignant" or "benign"
#'Attributes used:
#'
#'1. Clump Thickness
#'
#'2. Uniformity of cell size (1-10)
#'
#'3. Uniformity of Cell Shape (1-10)
#'
#'4. Marginal adhesion (1-10)
#
#'5. Single Epithelial Cell Size (1-10)
#'
#'6. Bare Nuclei        (1 - 10)
#'
#'7. Bland Chromatin (1-10)
#'
#'8. Normal Nucleoli (1-10)
#'
#'9. Mitoses(1-10)
#'
#'10. Type of Tumor (malignant, benign)
#'
#' @examples
#' data(wisc)
#' library(caret)
#' train<-createDataPartition(wisc$tumor, p = .8, list = F)###Create a training set using 80% of dataset
#' wisc.smv<-smvcir("tumor", data = wisc[train,], test = T) ###Build smvcir model on training set
#' bcpreds<-predict(wisc.smv, newdata = wisc[-train,], type = "prob")
#' head(bcpreds)  ###probability estimates
#'
#'
#' ###Get Coordinates
#' coords<-predict(wisc.smv, newdata = wisc, coordinates_only = TRUE)
#' coords$tumor<-wisc$tumor
#'
#'
#' plotSVM3d<-function(x, y, kernel = "radial", ...){
#'  open3d()
#'  plot3d(x, col = as.numeric(y)+1)
#'  svm_mod<-svm(x = x, y = y, kernel =paste(kernel), type = "C-classification", ...)
#'  n=100
#'  nnew = 50
#'  newdat.list = lapply(x, function(x) seq(min(x), max(x), len=nnew))
#'  newdat      = expand.grid(newdat.list)
#'  newdat.pred = predict(svm_mod, newdata=newdat, decision.values=T)
#'  newdat.dv   = attr(newdat.pred, 'decision.values')
#'  newdat.dv   = array(newdat.dv, dim=rep(nnew, 3))
#'  # Fit/plot an isosurface to the decision boundary
#'  contour3d(newdat.dv, level=0, x=newdat.list[[1]], y=newdat.list[[2]], z=newdat.list[[3]], add=T)
#'  return(list(svm_mod = svm_mod))
#'  }
#'
#'
#' plotSVM3d(x = coords[,1:3], coords[,10], kernel = "linear") ###Visualize a support vector machine model with our coordinates
#' plotSVM3d(x = coords[,1:3], coords[,10], kernel = "radial")
"wisc"
