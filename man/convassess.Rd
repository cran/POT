\name{convassess}
\alias{convassess}
\alias{convassess.uvpot}
\alias{convassess.bvpot}
\title{Convergence Assessment for Fitted Objects}
\description{
  This function uses several starting values to assess if the fitted
  object is sensitive to the starting values choice.
}
\usage{
\method{convassess}{uvpot}(fitted, n = 50)
\method{convassess}{bvpot}(fitted, n = 50)
}
\arguments{
  \item{fitted}{A fitted object of class \code{'uvpot'} or
    \code{'bvpot'}. Generally, an object return by \code{fitgpd} or
    \code{fitbvgpd}.}
  \item{n}{The number of starting values to be tested.}
}
\value{
  Graphics: the considered starting values, the objective values derived
  from numerical optimizations and traceplots for all estimated
  parameters. In addition, it returns invisibly all these informations.
}
\details{
  The starting values are defined using the unbiased probability
  weighted moments fitted on \code{n} bootstrap samples. 
}
\seealso{
  \code{\link{fitgpd}}, \code{\link{convassess.bvpot}}
}
\examples{
##Univariate Case
x <- rgpd(30, 0, 1, 0.2)
med <- fitgpd(x, 0, "med")
convassess(med)
##Bivariate Case
x <- rbvgpd(50, model = "log", alpha = 0.5, mar1 = c(0, 1, 0.2))
log <- fitbvgpd(x, c(0,0))
convassess(log)
}
\author{Mathieu Ribatet}
\keyword{htest}