estfinalmatD<-function (tsmvcird, stackvaract) 
{
  k <- stackvaract$delta5act$delta4act$delta3act$delta2act$delta1act$clt1act$k
  bigK <- kunstack(k, as.matrix(stackvaract$ef))
  sv <- svd(bigK, nu = nrow(bigK), nv = ncol(bigK))
  gammanot <- sv$u[, (tsmvcird + 1):(ncol(sv$u))]
  psinot <- sv$v[, (tsmvcird + 1):ncol(sv$v)]
  varfinalmatD <- ((t(psinot) %x% t(gammanot)) %*% stackvaract$vf %*% 
                     (psinot %x% gammanot))
  return(list(varfinalmatD = varfinalmatD, stackvaract = stackvaract, 
              tsmvcird = tsmvcird))
}