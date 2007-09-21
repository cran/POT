\name{Clusters}
\alias{exiplot}
\title{Extremal Index Plot}
\description{
  Plot estimates of the Extremal Index
}
\usage{
exiplot(data, u.range, tim.cond = 1, n.u = 50, xlab, ylab, \dots)
}
\arguments{
  \item{data}{A matrix/data.frame with two columns. Columns names \bold{must
      be} ``obs'' for observations and ``time'' for the associated date of
    each observation.}
  \item{u.range}{A numeric vector of length 2. Specify the range of
    threshold for which the Extremal Index is estimated.}
  \item{tim.cond}{A time condition to ensure independence between
    events. Should be in the same unit that \code{data[,"time"]}.}
  \item{n.u}{Numeric. The number of thresholds at which the Extremal
    Index is estimated.}
  \item{xlab, ylab}{Optional character strings to label the x and y
    axis.}
  \item{\dots}{Optional options to be passed to the \code{\link{plot}}
    function.} 
}
\value{
  Returns invisibly a matrix with two columns. The first one ``thresh''
  giving the threshold and the second one ``exi'' the related Extremal
  Index estimate.
}

\seealso{
  \code{\link{clust}}
}
\author{Mathieu Ribatet}
\keyword{htest}
