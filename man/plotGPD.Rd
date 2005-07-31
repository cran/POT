\name{Plot GPD fit}
\alias{plotgpd}
\alias{qq.gpd}
\alias{pp.gpd}
\alias{dens.gpd}
\alias{retlev.gpd}

\title{Graphical Functions to compare the GP Fitted Model to Empirical One}
\description{Produces QQ-plot, Probability Plot and a Denstiy Plot of
  the fitted model versus the empirical one. Another function computes
  the Return Level Plot of the fitted model.}

\usage{plotgpd(fitted, mu, main, which = 1:4, ask = nb.fig < length(which) && 
                     dev.interactive(),ci = TRUE, ...)
}
\arguments{
  \item{fitted}{A fitted object of class \code{'gpd'}. Generally, an object
    return by \code{fitgpd}}
  \item{mu}{The mean number that occur in a block -generally a block
    consits of a year.}
  \item{main}{optional. A string vector corresponding to the title of
    each plot.}
  \item{which}{a numeric vector which specifies which plot must be drawn
    : \code{'1'} for Probability Plot, \code{'2'} for QQ-Plot,\code{'3'}
    for Density Plot and \code{'4'} for a Return Level Plot.}
  \item{ask}{Logical. If \code{TRUE}, user is asked before each plot.}
  \item{ci}{Logical. If \code{TRUE}, the simulated 95\% confidence
    interval is plotted.}
  \item{...}{Other parameters to pass to the \code{\link{plot}} function.}
}
\examples{
data(ardieres)
fitted <- fitgpd(ardieres[,"flows"], 6, 'mle')
mu <- fitted$nhigh / 33.4
plotgpd(fitted, mu = mu)
}		   
\keyword{hplot}