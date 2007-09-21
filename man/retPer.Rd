\name{Return Periods Tools}
\alias{rp2prob}
\alias{prob2rp}
\title{Converts Return Periods to Probability and Vice Versa}
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
\details{
  The return period is defined by:
  \deqn{T = \frac{1}{npy (1-p)}}{%
    T = 1 / (npy (1-p))}
  where \eqn{npy} is the mean number of events per year (block), \eqn{p}
  is the probability of non exceedance.
}
\examples{
rp2prob(50, 1.8)
prob2rp(0.6, 2.2)
}
\author{Mathieu Ribatet}
\keyword{htest}
