#'Predict new data using SMVCIR
#'
#'Takes a set of observations and projects them into the SMVCIR discriminant space based on a previously built SMVCIR model object.
#'Probabilities are calculated by measuring the mahalonobis distance of an observation from each of the class centroids based on the
#'original SMVCIR model object.
#'@param type Either "prob" for class probabilities, or "class" for class predictions.
#'@param coordinates_only  Default is FALSE.  If TRUE, SMVCIR discriminant coordinates are returned and
#'can be used towards other classification techniques, such as SVM's or logistic regression.
#'@param maxdim Desired number of smvcir dimensions.  A dimensionality test from previously built model
#'may indicate a sufficient number of dimensions, by default, all possible dimensions are included.
#'@export

predict.smvcir<-function(model, newdata, type = "prob", coordinates_only = FALSE, method = "centroid",
                         maxdim = model$predictors){


  muhat_z<-model$muhat_z[,1:maxdim]
  sighatz<-model$sighatz[1:maxdim,1:maxdim]
  c_mat<-model$coefficients
  k<-model$predictors
  g<-model$groups


  compcases <- data.frame(complete.cases(newdata))
  compcases <- cbind(compcases, t(t(1:nrow(newdata))))
  data_nu <- na.omit(newdata)
  namelist <- attr(newdata, "names")
  groupindex <- model$groupindex


  if (groupindex == dim(data_nu)[2]) {
    stage1readydat <- data.matrix(data_nu)
  }  else if (groupindex == 1) {
    stage1readydat <- data.matrix(data_nu[, c(2:dim(data_nu)[2],
                                              1)])
  }  else {
    stage1readydat <- data.matrix(data_nu[, c(1:(groupindex -
                                                   1), (groupindex + 1):dim(data_nu)[2], groupindex)])
  }



  rm(data_nu)
  n <- dim(stage1readydat)[1]

  standat <- scale(stage1readydat[,1:k], center = model$xbar, scale = sqrt(diag(cov(model$originalx))))

  xx<-standat%*%c_mat
  colnames(xx)<-c(paste("D", rep(1:k), sep = ""))
  xx<-xx[,1:maxdim]
  ###Mahalanobis Distance from centroid

  centroid.distem<-matrix(NA, nrow(xx), g)
  for(k in 1:nrow(xx)){
    for(i in 1:g){
      centroid.distem[k,i]<-mahalanobis(xx[k,], center = muhat_z[i,], cov = sighatz, tol = 1e-20)
    }
  }
  if(type=="class"){
    colnames(centroid.distem)<-levels(newdata[,groupindex])
    classpreds<-apply(centroid.distem, 1, which.min)
    classpreds<-colnames(centroid.distem)[classpreds]
    return(classpreds)
  }

  ###Probabilities based on chi-square distance
  prob_mat<-matrix(NA, nrow(centroid.distem), ncol(centroid.distem))
  for(k in 1:nrow(centroid.distem)){
    for(i in 1:ncol(centroid.distem)){
      prob_mat[k,i]<-pchisq(centroid.distem[k,i], df = nrow(sighatz), lower.tail = F)
    }
  }

  ind.prob<-matrix(NA, nrow(centroid.distem), ncol(centroid.distem))
  for(i in 1:ncol(centroid.distem)){
    ind.prob[,i] <- prob_mat[,i]/apply(prob_mat, 1, sum)
  }

  rownames(ind.prob)<-rownames(xx)
  colnames(ind.prob)<-levels(newdata[,groupindex])

  if(coordinates_only == FALSE){
    return(prob = ind.prob)
  } else {

    return(data.frame(xx))
  }
}
