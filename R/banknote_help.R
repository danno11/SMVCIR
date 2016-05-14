#'Banknote dataset
#'
#'@docType data
#'
#'@usage data(banknote)
#'
#'@keywords datasets
#'
#'
#'@references A Modern Approach to Regression with R by Dr. Charles Scheather
#'
#'@source \href{http://www.stat.tamu.edu/~sheather/book/data_sets.php}{Dataset Index Book Webpage}
#'
#'
#'Class Variable  "Y": "0" or "1"
#'Attributes used:
#'
#'1. Length
#'
#'2. Left
#'
#'3. Right
#'
#'4. Bottom
#
#'5. Top
#'
#'6. Diagonal
#'
#'7. Y
#'
#' @examples
#' data(banknote)
#'
#' banksmv<-smvcir(group = "Y", data = banknote, test = TRUE)
#'
#' plotcoords(banksmv, dimensions = 1:3, build_svm = TRUE)
#'
"banknote"
