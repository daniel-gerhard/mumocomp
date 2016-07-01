mumotest <-
function(mlist, K, B){
  namek <- rownames(K)
  if (is.null(namek)) namek <- paste("C", 1:nrow(K), sep="") 
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

  out <- list(coefficients=est,
              vcov=covar,
              statistics=stat,
              resamp=bs,
              names=namek,
              N=n)
  class(out) <- "mumo"
  return(out)
}
