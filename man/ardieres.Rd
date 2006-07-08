\name{Flood Flows}
\alias{ardieres}
\title{High Flood Flows of river Ardieres at Beaujeu}
\usage{data(ardieres)}
\description{
  A data frame containing flood discharges, in units of cubic meters per
  second, of the Ardieres River at Beaujeu (FRANCE), over a period of 33
  years and the related date of those events.}
\format{A data frame with two columns: "time" and "obs".}
\examples{
data(ardieres)
plot(ardieres, xlab = "Time (Years)", ylab = expression(paste("Flood
discharges ", m^2/s, sep="")), type = "l")
}
\author{Mathieu Ribatet}
\keyword{datasets}
