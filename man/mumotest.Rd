\name{mumotest}
\alias{mumotest}
\title{Simultaneous Inference for Parameters in Multiple Models}
\description{Resampling-based testing of general linear hypotheses for parameters of multiple models}
\usage{
mumotest(mlist, K, B)
}
\arguments{
  \item{mlist}{List of model objects}
  \item{K}{A contrast matrix (Number of columns = Number of coefficients in all models)}
  \item{B}{Number of iterations}
}

\value{An object of class mumo}
\author{Daniel Gerhard}

\examples{
### data example from Hasler & Hothorn (2011), 
### A Dunnett-Type Procedure for Multiple Endpoints,
### The International Journal of Biostatistics: Vol. 7: Iss. 1, Article 3.
### DOI: 10.2202/1557-4679.1258
# but using two-sided inference

### see ?coagulation
data("coagulation", package = "SimComp")

### Marginal Models for each endpoint
(m1 <- lm(Thromb.count ~ Group, data = coagulation))
(m2 <- lm(ADP ~ Group, data = coagulation))
(m3 <- lm(TRAP ~ Group, data = coagulation))

### Dunnett contrast for comparisons to a control for each endpoint
K <- rbind("Thromb: B - S"=c(0,1,0,0,0,0,0,0,0),
           "Thromb: H - S"=c(0,0,1,0,0,0,0,0,0),
           "ADP: B - S"   =c(0,0,0,0,1,0,0,0,0),
           "ADP: H - S"   =c(0,0,0,0,0,1,0,0,0),
           "TRAP: B - S"  =c(0,0,0,0,0,0,0,1,0),
           "TRAP: H - S"  =c(0,0,0,0,0,0,0,0,1))

### Resampling of contrast test statistics
mm <- mumotest(list(m1, m2, m3), K, B=10000)

### Adjusted p-values
summary(mm)

### Simultaneous confidence intervals
confint(mm)
}

\keyword{multivariate}
\keyword{htest }% use one of  RShowDoc("KEYWORDS")
