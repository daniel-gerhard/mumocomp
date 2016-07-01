mumocomp <-
function(mlist, K, B){
  # number of obs (same for every model)
  n <- nrow(model.matrix(mlist[[1]]))
  # evaluate score functions using sandwich
  psilist <- lapply(mlist, function(fm) t(bread(fm) %*% t(estfun(fm))))
  # contrast coefficient estimates
  est <- as.vector(K %*% unlist(lapply(mlist, function(fm) coef(summary(fm))[,1])))
  # sandwich covariance of parameters within and between models
  psi <- t(matrix(unlist(psilist), nrow=n))
  covorig <- ((psi %*% t(psi)) / n)
  corr <- cov2cor(covorig)
  # covariance of parameter linear combinations
  covar <- K %*% covorig %*% t(K)
  # test statistic
  stat <- est / (sqrt(diag(covar)) * (1/sqrt(n)))
  # resampling of test statistic
  bs <- cmr(psi, corr, n, K, B)
  # maximum test statistics
  bmax <- apply(bs, 1, function(x) max(abs(x)))
  # adjusted p-values
  alter <- sapply(stat, function(x) bmax > abs(x))
  pv <- apply(alter, 2, function(x) (sum(x) + 1)/(length(x) + 1))
  return(pv)
}
