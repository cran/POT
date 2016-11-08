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

## This file contains several functions to simulate pseudo-random
## numbers GPD distributed, evaluate the Cumulative Distribution
## Function, evaluate the Quantile Function and the density.

rgpd <- function(n, loc = 0, scale = 1, shape = 0){
  
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (shape == 0) 
        return(loc + scale * rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1)/shape)

  }

qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
    if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
        1) 
      stop("`p' must contain probabilities in (0,1)")
    if (min(scale) < 0) 
      stop("invalid scale")
    if (length(shape) != 1) 
      stop("invalid shape")
    if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
      stop("invalid lambda")
    if (any(p < lambda))
      stop("``p'' must satisfy ``p >= lambda''")
    if (lower.tail) 
        p <- 1 - p
    p <- p / (1 - lambda)
    if (shape == 0) 
      return(loc - scale * log(p))
    else return(loc + scale * (p^(-shape) - 1)/shape)

  }

dgpd <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE){
  
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if (shape == 0) {
        d[index] <- log(1/scale[index]) - d[index]
        d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) * log(1 + 
            shape * d[index])
        d[!index] <- -Inf
    }
    if (!log) 
        d <- exp(d)

    return(d)

  }

pgpd <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
      stop("invalid lambda")
    q <- pmax(q - loc, 0)/scale
    if (shape == 0) 
        p <- 1 - (1 - lambda) * exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - (1 - lambda) * p^(-1/shape)
    }
    if (!lower.tail) 
        p <- 1 - p

    return(p)

  }





