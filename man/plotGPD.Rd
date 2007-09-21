\name{plot.uvpot}
\alias{plot.uvpot}

\title{Graphical Diagnostic: the Univariate GPD Model}
\description{Produces QQ-plot, Probability Plot and a Density Plot of
  the fitted model versus the empirical one.  Another function computes
  the Return Level Plot of the fitted model.}

\usage{
\method{plot}{uvpot}(x, npy, main, which = 1:4, ask = nb.fig <
length(which) && dev.interactive(),ci = TRUE, \dots)
}
\arguments{
  \item{x}{A fitted object of class \code{'uvpot'}. Generally, an object
    return by \code{fitgpd}}
  \item{npy}{The mean \bold{N}umber of events \bold{P}er \bold{Y}ear - or
    more generally a block.}
  \item{main}{optional. A string vector corresponding to the title of
    each plot.}
  \item{which}{a numeric vector which specifies which plot must be drawn
    : \code{'1'} for Probability Plot, \code{'2'} for QQ-Plot,\code{'3'}
    for Density Plot and \code{'4'} for a Return Level Plot.}
  \item{ask}{Logical. If \code{TRUE}, user is asked before each plot.}
  \item{ci}{Logical. If \code{TRUE}, the simulated 95\% confidence
    interval is plotted.}
  \item{\dots}{Other parameters to pass to the \code{\link{plot}} function.}
}
\examples{
data(ardieres)
ardieres <- clust(ardieres, 4, 10 / 365, clust.max = TRUE)
fitted <- fitgpd(ardieres[, "obs"], 6, 'mle')
npy <- fitted$nat / 33.4 ##33.4 is total record length (in year)
par(mfrow=c(2,2))
plot(fitted, npy = npy)
}
\author{Mathieu Ribatet}
\keyword{hplot}
