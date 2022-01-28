\name{confint.uvpot}
\alias{confint.uvpot}
\title{Generic Function to Compute (Profile) Confidence Intervals}
\description{
  Compute (profile) confidence intervals for the scale, shape GPD
  parameters and also for GPD quantiles.
}
\usage{
\method{confint}{uvpot}(object, parm, level = 0.95, \dots, range, prob,
 prof = TRUE)
}
\arguments{
  \item{object}{\code{R} object given by function \code{\link{fitgpd}}.}
  \item{parm}{Charater string specifies for which variable confidence
    intervals are computed. One of \code{"quant"}, \code{"scale"} or
    \code{"shape"}.}
  \item{level}{Numeric. The confidence level.}
  \item{\dots}{Optional parameters. See details.}  
  \item{range}{Vector of dimension two. It gives the lower and upper
    bound on which the profile likelihood is performed. Only required
    when "prof = TRUE".}
  \item{prob}{The probability of non exceedance.}
  \item{prof}{Logical. If \code{TRUE} (the default), profile confidence
    intervals are computed. Otherwise, it is Fisher ones.}  
}
\value{
  Returns a vector of the lower and upper bound for the (profile) confidence
  interval.  Moreover, a graphic of the profile likelihood function is
  displayed when \code{prof = TRUE}.
}
\details{
  Additional options can be passed using "\dots" in the function
  call. Possibilites are related to the specific functions:
  \code{link{gpd.fiscale}}, \code{link{gpd.fishape}},
  \code{link{gpd.firl}}, \code{link{gpd.pfscale}},
  \code{link{gpd.pfshape}}, \code{link{gpd.pfrl}}.
}
\seealso{
  \code{link{gpd.fiscale}}, \code{link{gpd.fishape}},
  \code{link{gpd.firl}}, \code{link{gpd.pfscale}},
  \code{link{gpd.pfshape}} and \code{link{gpd.pfrl}}
}
\examples{
x <- rgpd(100, 0, 1, 0.25)
fmle <- fitgpd(x, 0)
confint(fmle, prob = 0.2)
confint(fmle, parm = "shape")
}
\author{Mathieu Ribatet}
\keyword{htest}