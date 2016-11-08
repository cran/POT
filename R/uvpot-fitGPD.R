#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################

## In this function, several methods to estimates the GPD parameters
## are available:
##   1) Moments Estimator
##   2) Unbiased Probability Weighted Moment (PWMU) Estimator
##   3) Biased Probability Weighted Moment (PWMB) Estimator
##   4) Maximum Likelihood Estimator
##   5) Pickands' Estimator
##   6) Minimum Density Power Divergence Estimator
##   7) Method of Medians Estimator
##   8) Likelihood Moment Estimator
##   9) Maximum Goodness-of-Fit Estimator
##   10) Maximum penalized likelihood estimator

## A generic function for estimate the GPD parameters
fitgpd <- function(data, threshold, est = "mle", ...){
  threshold.call <- deparse(threshold)
  if (!(est %in% c("moments", "pwmb", "pwmu", "mle", "pickands",
                   "mdpd", "med", "lme", "mgf", "mple")))
    stop("Unknown estimator. Please check the ``est'' argument.")
  
  fitted <- switch(est, 'moments' = gpdmoments(data, threshold),
                   'pwmb' = gpdpwmb(data, threshold, ...),
                   'pwmu' = gpdpwmu(data, threshold),
                   'mle' = gpdmle(data, threshold, ...),
                   'pickands' = gpdpickands(data, threshold),
                   'mdpd' = gpdmdpd(data, threshold, ...),
                   'med' = gpdmed(data, threshold, ...),
                   'lme' = gpdlme(data, threshold, ...),
                   'mgf' = gpdmgf(data, threshold, ...),
                   'mple' = gpdmple(data, threshold, ...)
                   )
  fitted$threshold.call <- threshold.call
  class(fitted) <- c("uvpot","pot")
  return(fitted)
}




