\name{Return Periods Tools}
\alias{rp2prob}
\alias{prob2rp}
\title{Tools to convert return periods to probability and vice versa}
\description{
  Compute return period from probability of non exceedance and vice versa.
}
\usage{
rp2prob(retper, npy)
prob2rp(prob, npy)
}
\arguments{
  \item{retper}{The return period.}
  \item{prob}{the probability of non exceedance.}
  \item{npy}{The mean \bold{N}umber of events \bold{p}er \bold{y}ear
    (block).}
}
\value{
  Returns a table with mean numbers of events per year, return periods
  and probabilities of non exceedance associated.
}
\examples{
rp2prob(50, 1.8)
prob2rp(0.6, 2.2)
}
\author{Mathieu Ribatet}
\keyword{htest}