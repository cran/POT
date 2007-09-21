\name{anova.uvpot}
\alias{anova.uvpot}
\title{Anova Tables: Univariate Case}
\description{
  Computes analysis of deviance for ``uvpot'' object 
}
\usage{
\method{anova}{uvpot}(object, object2, \dots)
}
\arguments{
  \item{object, object2}{Two objects of class ``uvpot'', most often return of the
    \code{\link{fitgpd}} function.}
  \item{\dots}{Other options to be passed to the \code{\link{anova}}
    function.}
}
\value{
  This function returns an object of class anova. These objects
  represent analysis-of-deviance tables.  
}
\seealso{
  \code{\link{anova}}, \code{\link{anova.bvpot}}
}
\examples{
x <- rgpd(1000, 0, 1, -0.15)
M0 <- fitgpd(x, 0, shape = -0.15)
M1 <- fitgpd(x, 0)
anova(M0, M1) 
}
\author{Mathieu Ribatet}
\keyword{models}
