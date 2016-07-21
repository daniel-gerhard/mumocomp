mumotest <- function(mlist, K, B, margin=0){
  namek <- rownames(K)
  if (is.null(namek)) namek <- paste("C", 1:nrow(K), sep="") 
  if (length(margin) == 1 & length(margin) < nrow(K)) margin <- rep(margin, each=nrow(K))
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
  stat <- (est - margin) / (sqrt(diag(covar)) * (1/sqrt(n)))
  # resampling of test statistic
  bs <- cmr(psi, corr, n, K, B, margin)

  out <- list(coefficients=est,
              vcov=covar,
              statistics=stat,
              resamp=bs,
              names=namek,
              N=n,
              margin=margin)
  class(out) <- "mumo"
  return(out)
}
