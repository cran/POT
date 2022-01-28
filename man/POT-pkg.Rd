\name{POT-package}
\alias{POT}
\alias{POT-package}
\title{Overview of the \pkg{POT} package}

\description{

The \pkg{POT} package aims to provide operational tools to analyze peak over threshold. 
This package relies on the extreme value theory (EVT) to model the tail of any continuous distribution.
Tail modelling, in particular POT modelling, is of great importance for many financial 
and environmental applications.

The \pkg{POT} package was first committed to the CRAN in April 2005 and is 
still in development. 
The package is hosted in \href{https://pot.r-forge.r-project.org/}{R-forge}.
Since November 2016, the package has a new maintainer.

The main motivation was to provide practical tools for 
probabilistic modelling of high flood flows. However, the strength of the EVT is 
that results do not depend on the process to be modelled. 
Thus, one can use the \pkg{POT} package to analyze precipitations, 
floods, financial times series, earthquakes and so on... 


The \pkg{POT} package can perform univariate and bivariate extreme value analysis; 
first order Markov chains can also be considered. 
For instance, the (univariate) GPD is currently fitted using 18 estimators. 
These estimators rely on three different techniques:
\itemize{
\item Likelihood maximization: MLE, LME, MPLE
\item Moment Approaches: MOM, PWM, MED
\item Distance Minimization: MDPD and MGF estimators.
}

Contrary to the univariate case, there is no finite parametrisation to model bivariate 
exceedances over thresholds. The POT packages allows 6 parametrisation for the 
bivariate GPD: the logistic, negative logistic and mixed models - 
with their respective asymmetric counterparts.

Lastly, first order Markov chains can be fitted using the bivariate GPD for the joint distribution of two consecutive observations. 

We have written a package \href{../doc/POT.pdf}{vignette} to help new users. 
This users guide is a part of the package - just run \code{vignette("POT")} 
once the package is loaded.


}

\author{ 
Mathieu Ribatet and Christophe Dutang.
}

