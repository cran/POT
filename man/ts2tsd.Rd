\name{ts2tsd}
\alias{ts2tsd}
\title{Mobile Window on a Time Series}
\description{
  This function performs a mobile average windows on the whole time
  series. Thus, if the time series represents flood discharges, it
  returns the averaged discharges over a specific duration.
}
\usage{
ts2tsd(ts, d, vol = FALSE, method = "linear")
}
\arguments{
  \item{ts}{The time series. It consists of two columns: one named
    \code{"time"} and the second \code{"obs"}.}
  \item{d}{Numeric which corresponds of the duration for the mobile
    window.}
  \item{vol}{Logical. If \code{FALSE} -the default, average values are
    computed, else volumes.}
  \item{method}{Specifies the interpolation method to be used.  Choices
    are \code{"linear"} or \code{"constant"}.}
}
\value{
  Returns a time series like object \code{ts}. In particular
  \code{ts[,"time"]} and \code{tsd[,"time"]} are identical.
}
\details{
  A mobile windows of length \code{d} is performed on the whole time
  sire. The ``discrete'' time series in first transformed in a
  function; interpolation are obtained using the \code{\link{approx}}
  function. Thus, if f(t) is the function representing the time series,
  volume over duration \code{d} is defined by:
  \deqn{vol(t) = \int_{t-d/2}^{t+d/2} f(u)du}{vol(t) =
    \int_{t-d/2}^{t+d/2} f(u)du}
  while average values are:
  \deqn{ave(t) = \frac{1}{d}\int_{t-d/2}^{t+d/2} f(u)du}{ave(t) = 1/d
    \int_{t-d/2}^{t+d/2} f(u)du}
}
\section{Warnings}{
  Please note that as the time series is interpolated, caution should be
  taken if the method to interpolate is not efficient. 
  
  Note that object \code{d} should have the same unit than
  \code{ts[,"time"]}.
}
\seealso{
  \code{\link{approx}}
}
\examples{
data(ardieres)
tsd <- ts2tsd(ardieres, 3 / 365)
plot(ardieres, type = "l", col = "blue")
lines(tsd, col = "green")
}
\author{Mathieu Ribatet}
\keyword{ts}
