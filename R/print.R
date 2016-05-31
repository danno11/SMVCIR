#'@export
print.smvcir<-function(model,...){
  cat("Call: \n")
  print(model$call)
  print(model$summary1)
  cat("\nDimensionality Test P-Values\n")
  print(model$pvalmat)
}
