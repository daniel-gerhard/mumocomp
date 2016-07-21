\name{cmr}
\alias{cmr}
\title{Conditional Multiplier Resampling}
\description{Multimodel resampling based on score functions}
\usage{
cmr(psi, X, n, K, reps, margin)
}
\arguments{
  \item{psi}{Matrix with evaluated scores}
  \item{X}{Correlation Matrix}
  \item{n}{Number of observations}
  \item{K}{Contrast matrix}
  \item{reps}{Number of resampling iterations}
  \item{margin}{Vector of test margins}
}

\value{A matrix with test statistics}
\author{Daniel Gerhard}


\keyword{multivariate}
\keyword{htest }% use one of  RShowDoc("KEYWORDS")
