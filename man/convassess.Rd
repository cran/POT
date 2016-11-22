\name{convassess}
\alias{convassess}
\alias{convassess.uvpot}
\alias{convassess.bvpot}
\alias{convassess.mcpot}
\title{Convergence Assessment for Fitted Objects}
\description{

\code{convassess} is a generic function used to assess the convergence of 
the estimation procedure to the global maximum. The function invokes particular \code{methods} 
which depend on the \code{\link{class}} of the first argument.
This function uses several starting values to assess the sensitiveness of the
fitted object with respect to starting values.
}
\usage{
convassess(fitted, n = 50)
  
\method{convassess}{uvpot}(fitted, n = 50)
\method{convassess}{mcpot}(fitted, n = 50)
\method{convassess}{bvpot}(fitted, n = 50)
}
\arguments{
  \item{fitted}{A fitted object. When using the POT package, an object 
  of class \code{'uvpot'}, \code{'mcpot'}
    or \code{'bvpot'}. Generally, an object return by \code{fitgpd},
    \code{fitmcgpd} or \code{fitbvgpd}.}
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
  \code{\link{fitgpd}}, \code{\link{fitmcgpd}}, \code{\link{fitbvgpd}}
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