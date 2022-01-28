\name{pickdep}
\alias{pickdep}
\title{The Pickands' Dependence Function}
\description{
  Return and optionally plot the Pickands' dependence function.
}
\usage{
pickdep(object, main, bound = TRUE, plot = TRUE, \dots)
}
\arguments{
  \item{object}{A object of class \code{bvpot}. Usually, \code{object}
    is the return of function \code{\link{fitbvgpd}}.}
  \item{main}{May be missing. If present, the plot title.}
  \item{bound}{Logical. Should the perfect dependent and independent
    case bounds be plotted?}
  \item{plot}{Logical. Should the dependence function be plotted?}
  \item{\dots}{Optional parameters to be passed to the
    \code{\link{lines}} function.}
}
\value{
  The function returns an invisible function: the Pickands' dependence
  function. Moreover, the returned object has an attribute which
  specifies the model for the bivariate extreme value distribution.

  If \code{plot = TRUE}, then the dependence function is plotted.
}
\details{
  It is common to parametrize a bivariate extreme value distribution
  according to the Pickands' representation (Pickands, 1981). That is,
  if \eqn{G}{G} is any bivariate extreme value distribution, then it has
  the following parametrization:
  \deqn{G\left(y_1,y_2\right) = \exp\left[- \left(\frac{1}{z_1} +
    \frac{1}{z_2} \right) A\left( \frac{z_2}{z_1+z_2} \right)
    \right]}{G(y_1, y_2) = exp[ -(1/z_1 + 1/z_2) A(z_2 / (z_1 + z_2))]}
  where \eqn{z_i}{z_i} are unit Frechet.

  \eqn{A}{A} is the Pickands' dependence function. It has the following
  properties:
  \itemize{
    \item \eqn{A}{A} is defined on [0,1];
    \item \eqn{A(0)=A(1)=1}{A(0) = A(1) = 0};
    \item \eqn{\max \left(w, 1-w \right) \leq A(w) \leq 1, \quad \forall
      w}{max(w, 1-w) <= A(w) <= 1, for all w};
    \item \eqn{A}{A} is a convex function;
    \item For two independent (unit Frechet) random variables, \eqn{A(w)
      = 1, \quad \forall w}{A(w) = 1, for all w};
    \item For two perfectly dependent (unit Frechet) random variables,
    \eqn{A(w) = \max (w, 1-w)}{A(w) = max(w, 1-w)}.
  }
}
\references{
  Pickands, J. (1981) Multivariate Extreme Value Distributions
  \emph{Proceedings 43rd Session International Statistical Institute}
}
\examples{
x <- rbvgpd(1000, alpha = 0.9, model = "mix", mar1 = c(0,1,0.25),
 mar2 = c(2,0.5,0.1))
Mmix <- fitbvgpd(x, c(0,2), "mix")
pickdep(Mmix)
}
\author{Mathieu Ribatet}
\keyword{models}
\keyword{hplot}
