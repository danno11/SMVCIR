#'Plot smvcir discriminant coordinates
#'
#'Plots specified dimensions of smvcir space.
#'
#'@param model smvcir model to plot coordinates
#'@param dimensions Dimensions to be plotted. Must be exactly three if a 3d plot is desired.  Otherwise pairwise plots are produced.
#'@param GL If GL = TRUE and three dimensions above are specified, a three-dimensional plot using the rgl package will be produced.
#'@param build_svm If build_svm=TRUE, a support vector machine model is built and can be visualized with respect to the
#'plotted dimensions.
#'
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
#'
#'
#'@export
plotcoords<-function (x, dimensions = c(1, 2), groups = 1:10, GL = FALSE,
                       ID = FALSE, build_svm = FALSE, kernel = NULL)
{

  plotit(x, dimensions, groups, GL, ID, build_svm, kernel, method = "SMVCIR")
}
