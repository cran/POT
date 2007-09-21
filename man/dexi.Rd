\name{dexi}
\alias{dexi}
\title{Compute the Density of the Extremal Index}
\description{
  Compute the density of the extremal index using simulations from a
  fitted markov chain model.
}
\usage{
dexi(x, n.sim = 1000, n.mc = length(x$data), plot = TRUE, \dots)
}
\arguments{
  \item{x}{A object of class \code{"mcpot"} - most often the returned
    object of the \code{\link{fitmcgpd}} function.}
  \item{n.sim}{The number of simulation of Markov chains.}
  \item{n.mc}{The length of the simulated Markov chains.}
  \item{plot}{Logical. If \code{TRUE} (default), the density of the
    extremal index is plotted.}
  \item{\dots}{Optional parameters to be passed to the
    \code{\link{plot}} function.}
}
\value{
  The function returns a optionally plot of the kernel density estimate
  of the extremal index. In addition, the vector of extremal index
  estimations is returned invisibly.
}
\details{
  The Markov chains are simulated using the \code{\link{simmc}}
  function to obtained dependent realisations \eqn{u_i}{u_i} of standard
  uniform realizations. Then, they are transformed to correspond to the
  parameter of the fitted markov chain model. Thus, if \eqn{u, \sigma,
    \xi}{u, sigma, xi} is the location, scale and shape parameters ; and
  \eqn{\lambda}{lambda} is the probability of exceedance of \eqn{u}{u},
  then by defining :
  
  \deqn{\sigma_* = \xi * \frac{u}{\lambda^{-\xi} - 1}}{sigma_* = xi * u
    / (lambda^(-xi) - 1)}

  the realizations \eqn{y_i = qgpd(u_i, 0, \sigma_*, \xi)}{y_i =
    qgpd(u_i, 0, sigma_*, xi)} are distributed such as the probability
  of exceedance of \eqn{u}{u} is equal to \eqn{\lambda}{lambda}.

  At last, the extremal index for each generated Markov chain is
  estimated using the estimator of Ferro and Segers (2003) (and thus
  avoid any declusterization).
}
\references{
  Fawcett L., and Walshaw D. (2006) Markov chain models for extreme wind
  speed. Environmetrics, \bold{17}:(8) 795--809.
  

  Ferro, C. and Segers, J. (2003) Inference for clusters of extreme
  values. Journal of the Royal Statistical Society. Series B
  \bold{65}:(2)  545--556.

  
  Ledford A., and Tawn, J. (1996) Statistics for near Independence in
  Multivariate Extreme Values. Biometrika, \bold{83} 169--187.
  

  Smith, R., and Tawn, J., and Coles, S. (1997) Markov chain models for
  threshold exceedances. Biometrika, \bold{84} 249--268.
}
\seealso{
  \code{\link{simmc}}, \code{\link{fitmcgpd}}, \code{\link{fitexi}}
}
\examples{
mc <- simmc(1000, alpha = 0.25)
mc <- qgpd(mc, 0, 1, 0.25)
log <- fitmcgpd(mc, 2, shape = 0.25, scale = 1)
dexi(log, n.sim = 200)
}
\author{Mathieu Ribatet}
\keyword{models}
