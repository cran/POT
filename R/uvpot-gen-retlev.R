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

## This file contains several functions to plot Peaks Over a Threshold.

retlev.uvpot <- function(fitted, npy, main, xlab,
                         ylab, xlimsup, ci = TRUE, points = TRUE,
                         ...){
  ## Plot the return level plot of a POT model fitted
  ## Input : ``fitted'' is a POT fitted model, result of function
  ##         ``fitgpd'' or ``fitpp''
  ##         npy is the mean number of events per block -generally
  ##         per year- or equivalently the mean -intensity- of the
  ##         Poisson processus.
  if (!inherits(fitted, "uvpot"))
    stop("Use only with 'uvpot' objects")
  #if(!inherits(fitted, "uvpot"))
  #  stop("retlev.uvpot is only valid for uvpot objects")
  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceed
  loc <- fitted$threshold[1]
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  
  n <- fitted$nat
  
  pot.fun <- function(T){
    p <- rp2prob(T, npy)[,"prob"]
    return(qgpd(p,loc,scale,shape))
  }

  eps <- 10^(-3)

  if (!is.null(fitted$noy))
    npy <- n / fitted$noy

  else
    if (missing(npy)){
      warning("Argument ``npy'' is missing. Setting it to 1.")
      npy <- 1
    }
  if (missing(main)) main <- 'Return Level Plot'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- prob2rp((n - .35)/n, npy)[,"retper"] 
  
  plot(pot.fun, from= 1 / npy + eps, to = xlimsup, log='x',
       xlab = xlab, ylab = ylab, main = main, ...)

  if (points){
    p_emp <- (1:n -.35) / n
    points(1 / ( npy * (1 - p_emp) ), sort( data ), pch = 1)
  }

  if (ci){
    p_emp <- (1:n - .35 ) / n
    samp <- rgpd(1000*n, loc, scale, shape)
    samp <- matrix(samp, n, 1000)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- samp[25,]
    ci_sup <- samp[975,]
    lines( 1 / ( npy * (1 - p_emp) ), ci_inf, lty = 2)
    lines( 1 / ( npy * (1 - p_emp) ), ci_sup, lty = 2)
  }

  invisible(pot.fun)
}


