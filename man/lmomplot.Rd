\name{Threshold Selection}
\alias{lmomplot}
\title{Threshold Selection: The L-moments Plot}
\description{
  Plots of sample L-Skewness ans L-Kurtosis estimates at various
  thresholds for peaks over threshold modelling, using the Generalized
  Pareto parametrization.}
\usage{
lmomplot(data, u.range, identify = TRUE, \dots)
}
\arguments{
  \item{data}{A numeric vector.}
  \item{u.range}{A numeric vector of length two, giving the limits for
    the thresholds at which the model is fitted.}
  \item{identify}{Logical. If \code{TRUE}, points on the plot are
    identify using \code{\link{identify}} function.}
  \item{\dots}{Other arguments to be passed to the model fit
    function \code{fitgpd}.}
}
\details{
  For each thresholds, sample L-skewness and L-kurtosis are computed. If
  data are GP distributed, one have :
  \deqn{ \tau_4 = \frac{\tau_3 \left( 1 + 5 \tau_3 \right)}{5 + \tau_3}
  }{%
    tau_4 = tau_3 ( 1 + 5 tau_3 ) / ( 5 + tau_3 )}
  So, a threshold is acceptable if sample \eqn{\left(\tau_3,
    \tau_4\right)} are near the theoretical curve.
}
\section{Warnings}{
  L-moments plot are really difficult to interpret.  It can help us to
  say if the GP distribution is suited to model data.
}
\seealso{\code{\link{fitgpd}}, \code{\link{mrlplot}},
  \code{\link{tcplot}}}

\references{Hosking, J. R. M. and Wallis, J. R. (1997) \emph{Regional
    Frequency Analysis}. Cambridge University Press.

  Begueria, S. (2005) Uncertainties in partial duration series
    modelling of extremes related to the choice of the threshold value.
  Journal of Hydrology, 303\bold{(1-4)}: 215--230.
}
\examples{
data(ardieres)
ardieres <- clust(ardieres, 4, 10 / 365, clust.max = TRUE)
flows <- ardieres[, "obs"]
lmomplot(flows, identify = FALSE)
}
\author{Mathieu Ribatet}
\keyword{hplot}
