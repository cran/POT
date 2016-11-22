\name{simmcpot}
\alias{simmcpot}
\title{Simulate an Markov Chain with a Fixed Extreme Value Dependence
  from a Fitted mcpot Object}
\description{
Simulate a synthetic Markov chain from a fitted \code{'mcpot'} object.
}
\usage{
simmcpot(fitted, plot = TRUE, \dots)
}
\arguments{
\item{fitted}{An object of class \code{'mcpot'}; most often the returned
    object of the \code{\link{fitmcgpd}} function.}
  \item{plot}{Logical. If \code{TRUE} (the default), the simulated
    Markov chain is plotted.}
  \item{\dots}{Other optional arguments to be passed to the
    \code{\link{plot}} function.}
}
\details{
  The simulated Markov chain is computed as follows:
  \enumerate{
    \item Simulate a Markov chain \code{prob} with uniform margins on
    (0,1) and with the fixed extreme value dependence given by
    \code{fitted};
    \item For all \code{prob} such as \eqn{prob \leq 1 - pat}{prob <= 1
      - pat}, set \eqn{mc = NA}{mc = NA} (where \code{pat} is given by
    \code{fitted$pat});
    \item For all \code{prob} such as \eqn{prob \geq 1 - pat}{prob >= 1
      - pat}, set \eqn{prob2 = \frac{prob - 1 + pat}{pat}}{prob2 = (prob
      - 1 + pat) / pat}. Thus, \code{prob2} are uniformly distributed on
    (0,1);
    \item For all \code{prob2}, set \code{mc = qgpd(prob2, thresh,
      scale, shape)}, where \code{thresh, scale, shape} are given by the
    \code{fitted$threshold, fitted$param["scale"]} and
    \code{fitted$param["shape"]} respectively.
  }
}
\value{
 A Markov chain which has the same features as the fitted object. If
 \code{plot = TRUE}, the Markov chain is plotted.
}
\author{Mathieu Ribatet}
\seealso{\code{\link{fitmcgpd}}, \code{\link{simmc}}}
\examples{
data(ardieres)
flows <- ardieres[,"obs"]

Mclog <- fitmcgpd(flows, 5)
par(mfrow = c(1,2))
idx <- which(flows <= 5)
flows[idx] <- NA
plot(flows, main = "Ardieres Data")
flowsSynth <- simmcpot(Mclog, main = "Simulated Data")
}
\keyword{models}

