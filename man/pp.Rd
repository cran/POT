\name{pp}
\alias{pp}
\alias{pp.uvpot}
\title{Probability Probability Plot}
\description{
\code{pp} is a generic function used to show probability-probability plot. 
The function invokes particular \code{methods} 
which depend on the \code{\link{class}} of the first argument.
So the function makes a probability probability plot for univariate POT models. 
}
\usage{
pp(object, \dots)

\method{pp}{uvpot}(object, main, xlab, ylab, ci = TRUE, \dots)
}
\arguments{
  \item{object}{A fitted object. When using the POT package, an object 
  of class \code{'uvpot'}. Most often, the
    return of the \code{\link{fitgpd}} function.}
  \item{main}{The title of the graphic. If missing, the title is set to
\code{"Probability plot"}.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
set to \code{"Empirical"} and \code{"Model"} respectively.}
  \item{ci}{Logical. If \code{TRUE} (the default), 95\% intervals are
    plotted.}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
}
\value{
  A graphical window.
}
\details{
  The probability probability plot consists of plotting the theoretical
  probabilities in function of the empirical model ones. The theoretical
  probabilities are computed from the fitted GPD, while the empirical
  ones are computing from a particular plotting position estimator. This
  plotting position estimator is suited for the GPD case (Hosking,
  1995) and are defined by:

  \deqn{p_{j:n} = \frac{j - 0.35}{n}}{p_{j:n} = (j - 0.35) / n}
  where \eqn{n}{n} is the total number of observations.

  If the theoretical model is correct, then points should be ``near''
  the line \eqn{y=x}{y=x}.
}
\seealso{
  \code{\link{qq}}, \code{\link{qq.uvpot}}
}
\examples{
x <- rgpd(75, 1, 2, 0.1)
pwmb <- fitgpd(x, 1, "pwmb")
pp(pwmb)
}
\references{
  Hosking, J. R. M. and Wallis, J. R. (1995). A comparison of unbiased
  and plotting-position estimators of L moments. Water Resources
  Research. \bold{31}(8): 2019--2025.
}
\author{Mathieu Ribatet}
\keyword{hplot}
