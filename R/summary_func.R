summary.smvcir<-function(model,...){
  cat("Call: \n")
  print(model$call)
  print(model$summary1)
  print(model$summary2)
  print(model$pvalmat)
  print(model$dimCorr)
  cat("\nCoefficients:\n")
  print(model$coefficients)

}
