#'@export
summary.smvcir<-function(model,...){
  cat("Call: \n")
  print(model$call)
  print(model$summary1)
  cat("\nDimensionality Test P-Values\n")
  print(model$pvalmat)
  cat("\nCorrelations between Dimensions \n")
  print(model$dimCorr)
  cat("\nSMVCIR dimensions should have low correlations.
      \nMaximum SMVCIR dimension correlation: ",model$maxTRANScm)
  cat("\nCoefficients:\n")
  print(model$coefficients)
}
