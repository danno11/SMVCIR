#'Plot smvcir discriminant coordinates
#'
#'Plots specified dimensions of smvcir space.
#'
#'@param x Either an smvcir model object or a data frame of discriminant coordinates.
#'@param dimensions Dimensions to be plotted. Must have three if a 3d plot is desired.  Otherwise pairwise plots are produced.
#'@param GL If GL = TRUE and three dimensions above are specified, a three-dimensional plot using the rgl package will be produced.
#'@param build_svm If build_svm=TRUE, a support vector machine model is built using the dimensions specificied as covariates,
#'and can be visualized with respect to the plotted dimensions.
#'@param svmModel If desired, this option may be used to plot new coordinates along with an svm model built from other data coordinates
#'(possibly desirable if training/testing splits are used).  The svmModel specified must also use the same covariates as those to be plotted.
#'See example below.
#'
#'
#'@examples
#'
#'library(caret)
#'train<-createDataPartition(pima$diabetes, p = .8, list = FALSE)
#'pim.smv<-smvcir("diabetes", data = pima[train,], test = T) ###Build smvcir model on training set
#'pred_coords<-predict(pim.smv, newdata = pima[-train,], coordinates_only = TRUE)
#'pred_coords$diabetes<-pima$diabetes[-train]
#'
#'plotcoords(pred_coords, dimensions = 1:8)
#'plotcoords(pred_coords, dimensions = 1:3, GL = TRUE, build_svm = TRUE)
#'plotcoords(pred_coords, dimensions = 1:3, build_svm = TRUE)
#'
#'data(banknote)
#'###Create training rows for Banknote data
#'train<-sample(1:nrow(banknote), nrow(banknote)*.8, replace = FALSE)
#'####Build smvcir model on training set
#'banksmv<-smvcir("Y", data = banknote[train,], test = TRUE)
#'####Build svm model on first three dimensions of smvcir discriminant coordinates
#'svmModel<-svm(Y~., data = banksmv$direct[,c(1:3,7) ], probability = TRUE)
#'
#'###Get discriminant coordinates from test set
#'bcoords<-predict(banksmv, newdata = banknote[-train,], coordinates_only = TRUE, maxdim = 3)
#'bcoords$Y<-banknote[-train,]$Y
#'
#'####Plot test set coordinates with respect to svm model built from training set
#'plotcoords(bcoords, dimensions = 1:3, GL = TRUE, build_svm = TRUE, svmModel = svmModel)
#'
#'@export
plotcoords<-function (x, dimensions = c(1, 2), groups = 1:10, GL = FALSE,
                       ID = FALSE, build_svm = FALSE, kernel = NULL, svmModel = NULL)
{

  plotit(x, dimensions, groups, GL, ID, build_svm, kernel, svmModel, method = "SMVCIR")
}
