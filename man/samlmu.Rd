\name{L-moments}
\alias{samlmu}
\title{Compute Sample L-moments}
\description{Compute the sample L-moments - unbiased version.}

\usage{samlmu(x, nmom = 4, sort.data = TRUE)}
\arguments{
  \item{x}{a vector of data}
  \item{nmom}{a numeric value giving the number of sample L-moments to
    be computed}
  \item{sort.data}{a locical which specifies if the vector of data x
    should be sorted or not.}
}
\value{
  This function returns a vector of length 'nmom' corresponding to the
  sample L-moments. Note that for orders greater or equal than 3 it is the L-moments
  ratio that is sample L-CV, sample L-Sknewness, \ldots 
}

\references{
  Hosking, J. R. M. (1990) L-moment analysis and estimation of
  order statistics. Journal of the Royal Statistical Society Series
  B, \bold{52}: 105--124.  
}
\examples{
x <- runif(50)
samlmu(x, nmom = 5)
}
\keyword{htest}
