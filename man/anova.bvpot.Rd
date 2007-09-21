\name{anova.bvpot}
\alias{anova.bvpot}
\title{Anova Tables: Bivariate Case}
\description{
  Computes analysis of deviance for ``bvpot'' object 
}
\usage{
\method{anova}{bvpot}(object, object2, \dots, half = FALSE)
}
\arguments{
  \item{object, object2}{Two objects of class ``bvpot'', most often return of the
    \code{\link{fitbvgpd}} function.}
  \item{\dots}{Other options to be passed to the \code{\link{anova}}
    function.}
  \item{half}{Logical. For some non-regular testing problems the deviance
    difference is known to be one half of a chi-squared random
    variable. Set half to \code{TRUE} in these cases.}
}
\value{
  This function returns an object of class anova. These objects
  represent analysis-of-deviance tables.  
}
\section{Warning}{
  Circumstances may arise such that the asymptotic distribution of the
  test statistic is not chi-squared. In particular, this occurs when the
  smaller model is constrained at the edge of the parameter space. It is
  up to the user recognize this, and to interpret the output correctly. 
  
  In some cases the asymptotic distribution is known to be one half of a
  chi-squared; you can set \code{half = TRUE} in these cases.
}
\seealso{
  \code{\link{anova}}, \code{\link{anova.uvpot}}
}
\examples{
x <- rgpd(1000, 0, 1, -0.25)
y <- rgpd(1000, 2, 0.5, 0)
M0 <- fitbvgpd(cbind(x,y), c(0, 2))
M1 <- fitbvgpd(cbind(x,y), c(0,2), model = "alog")
anova(M0, M1)

##Non regular case
M0 <- fitbvgpd(cbind(x,y), c(0, 2))
M1 <- fitbvgpd(cbind(x,y), c(0, 2), alpha = 1)
anova(M0, M1, half = TRUE)
}
\author{Mathieu Ribatet (Alec Stephenson for the ``Warning'' case)}
\keyword{models}
