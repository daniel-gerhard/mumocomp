mumotest <-
function(mlist, K, B){
  n <- nrow(model.matrix(mlist[[1]]))
  psilist <- lapply(mlist, function(fm) t(bread(fm) %*% t(estfun(fm))))
  psi <- t(matrix(unlist(psilist), nrow=n))
  est <- as.vector(K %*% unlist(lapply(mlist, function(fm) coef(summary(fm))[,1])))
  covorig <- ((psi %*% t(psi)) / n)
  covar <- K %*% covorig %*% t(K)
  stat <- est / (sqrt(diag(covar)) * (1/sqrt(n)))
  corr <- cov2cor(covorig)
  bs <- cmr(psi, corr, n, K, B)
  bmax <- apply(bs, 1, function(x) max(abs(x)))
  alter <- sapply(stat, function(x) bmax > abs(x))
  pv <- apply(alter, 2, function(x) (sum(x) + 1)/(length(x) + 1))
  return(pv)
}
