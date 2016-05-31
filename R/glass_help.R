#'Glass dataset
#'
#'@docType data
#'
#'@usage data(glass)
#'
#'@keywords datasets
#'
#'@source \href{http://archive.ics.uci.edu/ml/datasets/Glass+Identification}{UCI Machine Learning Repository}
#'
#'
#'Type of glass: (class attribute)
#'
#'-- 1 building_windows_float_processed
#'
#'-- 2 building_windows_non_float_processed
#'
#'-- 3 vehicle_windows_float_processed
#'
#'-- 4 vehicle_windows_non_float_processed (none in this database)
#'
#'-- 5 containers
#'
#'-- 6 tableware
#'
#'-- 7 headlamps
#'
#' Attributes used:
#'
#' 1. RI: refractive index
#'
#' 2. Na: Sodium (unit measurement: weight percent in corresponding oxide, as are attributes 4-10)
#'
#' 3. Mg: Magnesium
#'
#' 4. Al: Aluminum
#'
#' 5. Si: Silicon
#'
#' 6. K: Potassium
#'
#' 7. Ca: Calcium
#'
#' 8. Ba: Barium
#'
#' 9. Fe: Iron
#'
#' @examples
#'
#'data(glass)
#'table(glass$Type)
#'glass.smv<-smvcir(group = "Type", data = glass, test = TRUE)
#'summary(glass.smv)
#'
#'plotcoords(glass.smv, dimensions = 1:3)
#'
"glass"
