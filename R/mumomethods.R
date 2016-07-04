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
  
  dd <- data.frame(Estimate=object$coefficients,
                   Std.Err=sqrt(diag(object$vcov)/object$N),
                   Statistic=object$statistics,
                   pvalue=pv)
  rownames(dd) <- object$names
  out <- list(mumo=object, results=dd)
  class(out) <- "summary.mumo"
  return(out)
}

print.summary.mumo <- function(x, digits= max(3, getOption("digits") - 3), ...){
  cat("\t Multiple Model Contrast Tests\n\n") 
  printCoefmat(x$results, digits=digits, cs.ind=1:2, tst.ind=3, has.Pvalue=TRUE) 
}


confint.mumo <- function(object, parm, level = 0.95, ...){
  bmax <- apply(object$resamp, 1, function(x) max(abs(x)))
  quant <- quantile(bmax, level)
  upp <- object$coefficients + quant * sqrt(diag(object$vcov)/object$N) 
  low <- object$coefficients - quant * sqrt(diag(object$vcov)/object$N) 
  
  
  cid <- data.frame(Estimate=object$coefficients,
                    lower=low,
                    upper=upp)
  rownames(cid) <- object$names
  out <- list(mumo=object, 
              level=level,
              confint=cid)
  class(out) <- "confint.mumo"
  return(out)
}

print.confint.mumo <- function(x, digits=max(3, getOption("digits") - 3), ...){
  cat(paste("\t Simultaneous", x$level, "Confidence Intervals\n\n")) 
  print(x$confint, digits=digits, ...) 
}
