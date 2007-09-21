\name{clust}
\alias{clust}
\title{Identify Extreme Clusters within a Time Series}
\description{
  A function to identify clusters of exceedances of a time series.
}
\usage{
clust(data, u, tim.cond = 1, clust.max = FALSE, plot = FALSE,
only.excess = TRUE, \dots)
}
\arguments{
  \item{data}{A matrix/data.frame with two columns. Columns names \bold{must
      be} ``obs'' for observations and ``time'' for the associated date of
    each observation.}
  \item{u}{Numeric. A value giving the threshold.}
  \item{tim.cond}{A time condition to ensure independence between
    events. Should be in the same unit than \code{data[,"time"]}.}
  \item{clust.max}{Logical. If \code{FALSE} (the default), a list
    containing the clusters of exceedances is returned. Else, a matrix
    containing the cluster maxima and related dates is returned.}
  \item{plot}{Logical. If \code{TRUE}, identified
    clusters are displayed. Else (the default), no plot is produced.}
  \item{only.excess}{Logical. If \code{TRUE} (the default), only
    exceedances are plotted. Else, all observations are displayed.}
  \item{\dots}{Optional parameters to be passed in \code{\link{plot}}
    function.}
}
\value{
  If \code{clust.max} is \code{FALSE}, a list containing the clusters of
  exceedances is returned. Else, a matrix containing the cluster maxima,
  related dates and indices are returned.

  In any case, the returned object has an attribute \code{exi} giving
  an estimation of the Extremal Index, that is the inverse of the
  average cluster size.
}
\details{
  The clusters of exceedances are defined as follows:
  \item The first exceedance initiates the first cluster;
  \item The first observation under the threshold \code{u} ``ends'' the
  current cluster unless \code{tim.cond} does not hold;
  \item The next exceedance initiates a new cluster;
  \item The process is iterated as needed.

  
  This function differs from the function \code{clusters} of \code{evd}
  Package as independence condition i.e. \code{tim.cond} could be a
  ``temporal'' condition. That is, two events are considered independent
  if the inter-arrival time is greater than a fixed duration.

  However, it is also possible to used the ``index'' independence as in
  \code{\link{clust}} by setting \code{data[,"time"] =
  1:length(data[,"obs"])}.
}
\seealso{
  \code{clusters} of package \code{evd}.
}
\examples{
data(ardieres)
par(mfrow=c(1,2))
clust(ardieres, 4, 10 / 365)
clust(ardieres, 4, 10 / 365, clust.max = TRUE)
clust(ardieres, 4, 10 / 365, clust.max = TRUE, plot = TRUE)
##The same but with optional arguments passed to function ``plot''
clust(ardieres, 4, 10 / 365, clust.max = TRUE, plot = TRUE,
xlab = "Time (Years)", ylab = "Flood discharges",
xlim = c(1972, 1980))
}
\author{Mathieu Ribatet}
\keyword{manip}
