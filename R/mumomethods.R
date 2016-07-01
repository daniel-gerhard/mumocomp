print.mumo <- function(x, digits= max(3, getOption("digits") - 3), ...){
  cat("\t Multiple Model Contrast Tests\n\n") 
  dd <- data.frame(Estimates=x$coefficients,
                   Std.Err=sqrt(diag(x$vcov)/x$N),
                   row.names=x$names)
  printCoefmat(dd, digits=digits, cs.ind=1:2, tst.ind=NULL)
}

summary.mumo <- function(object, digits= max(3, getOption("digits") - 3), ...){
  # maximum test statistics
  bmax <- apply(object$resamp, 1, function(x) max(abs(x)))
  # adjusted p-values
  alter <- sapply(object$statistics, function(x) bmax > abs(x))
  pv <- apply(alter, 2, function(x) (sum(x) + 1)/(length(x) + 1))
  
  cat("\t Multiple Model Contrast Tests\n\n") 
  dd <- data.frame(Estimate=object$coefficients,
                   Std.Err=sqrt(diag(object$vcov)/object$N),
                   Statistic=object$statistics,
                   pvalue=pv,
                   row.names=object$names)
  printCoefmat(dd, digits=digits, cs.ind=1:2, tst.ind=3, has.Pvalue=TRUE) 
}


confint.mumo <- function(object, parm, level = 0.95, ...){
  bmax <- apply(object$resamp, 1, function(x) max(abs(x)))
  quant <- quantile(bmax, level)
  upp <- object$coefficients + quant * sqrt(diag(object$vcov)/object$N) 
  low <- object$coefficients - quant * sqrt(diag(object$vcov)/object$N) 
  
  cat(paste("\t Simultaneous", level, "Confidence Intervals\n\n")) 
  cid <- data.frame(Estimate=object$coefficients,
                    lower=low,
                    upper=upp,
                   row.names=object$names)
  return(cid)
}
