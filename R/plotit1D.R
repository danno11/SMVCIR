plotit1D<-function (wmat, dimension, groups = groups)
{
  graphdat <- wmat[, c(dimension, ncol(wmat))]
  zecolor <- rep(0, nrow(wmat))
  graphdat <- data.frame(graphdat, zecolor)
  names(graphdat)[1] <- "sortit"
  graphdat <- graphdat[order(graphdat$sortit), ]
  names(graphdat)[1] <- paste("D", dimension, sep = "")
  names(graphdat)[2] <- names(wmat)[ncol(wmat)]
  #zecolor <- rep(0, nrow(graphdat))
  #graphdat <- data.frame(graphdat, zecolor)
  names(graphdat)[2] <- "group"
  graphdat$levs<-as.numeric(graphdat$group)

  graphdat$zecolor <- colors()[(graphdat$levs == 1) * 24 +
                                 (graphdat$levs == 2) * 552 + (graphdat$levs == 3) *
                                 26 + (graphdat$levs == 4) * 498 + (graphdat$levs == 5) * 547 + (graphdat$levs == 6) * 32 +
                                 (graphdat$levs == 7) * 254 + (graphdat$levs == 8) * 536 + (graphdat$levs == 9) * 652 +
                                 (graphdat$levs == 10) * 8]


  dotchart(graphdat[, 1], groups = graphdat$group, xlab = paste("D",
                                                                dimension, sep = ""), lcolor = "white", gcolor = "black",
           color = graphdat$zecolor)
}
