#'Pima Indian diabetes dataset
#'
#'@docType data
#'
#'@usage data(pima)
#'
#'@keywords datasets
#'
#'@references UCI Machine Learning Repository
#'
#'@source \href{http://archive.ics.uci.edu/ml/datasets/Pima+Indians+Diabetes}{UCI Machine Learning Repository}
#'
#'
#'Class Variable: "diabetes" 0 = no diabetes, 1 = diabetes
#'
#'Attributes used:
#'
#'1. Number of times pregnant
#'
#'2. Plasma glucose concentration a 2 hours in an oral glucose tolerance test
#'
#'3. Diastolic blood pressure (mm Hg)
#'
#'4. Triceps skin fold thickness (mm)
#'
#'5. 2-Hour serum insulin (mu U/ml)
#'
#'6. Body mass index (weight in kg/(height in m)^2)
#'
#'7. Diabetes pedigree function
#'
#'8. Age (years)
#'
#'
#' @examples
#' data(pima)
#' library(caret)
#' train<-createDataPartition(pima$diabetes, p = .8, list = F)###Create a training set using 80% of dataset
#' pim.smv<-smvcir("diabetes", data = pima[train,], test = T) ###Build smvcir model on training set
#' preds<-predict(pim.smv, newdata = pima[-train,], type = "class")
#' table(preds, pima$diabetes[-train])  ###Check accuracy
#'
#'
#' ###Get Coordinates
#' pred_coords<-predict(pim.smv, newdata = pima, coordinates_only = TRUE)
#' pred_coords$diabetes<-pima$diabetes
#'
#' library(e1071)
#' svm_mod<-svm(diabetes~., data = pred_coords[train,], kernel = "radial")  ###Build an SVM model and check accuracy
#' svmp<-predict(svm_mod, newdata = pred_coords[-train,])
#' confusionMatrix(svmp, pred_coords$diabetes[-train], positive = "1")
"pima"




