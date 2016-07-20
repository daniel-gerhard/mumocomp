\name{mumo methods}
\alias{confint.mumo}
\alias{summary.mumo}
\alias{print.mumo}
\alias{print.confint.mumo}
\alias{print.summary.mumo}
\title{Methods for mumo objects}
\description{Resampling-based testing of general linear hypotheses for parameters of multiple models}
\usage{
\method{confint}{mumo}(object, parm, level = 0.95, ...)
\method{summary}{mumo}(object, digits=max(3, getOption("digits") - 3), ...)
\method{print}{mumo}(x, digits=max(3, getOption("digits") - 3), ...)
\method{print}{confint.mumo}(x, digits=max(3, getOption("digits") - 3), ...)
\method{print}{summary.mumo}(x, digits=max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{mumo object}
  \item{x}{mumo or related object}
  \item{parm}{ignore this}
  \item{level}{confidence level}
  \item{digits}{Number of digits to show}
  \item{...}{dots}
}

\value{An object of class confint.mumo or summary.mumo}
\author{Daniel Gerhard}

\keyword{multivariate}
\keyword{htest }% use one of  RShowDoc("KEYWORDS")
