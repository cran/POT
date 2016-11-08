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


##   5) Pickands' Estimator



##Pickand's Estimator
gpdpickands <- function(data, threshold){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'pickands'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  excess <- sort(exceed - threshold)
  
  n <- length(excess)
  xn.2 <- excess[ceiling(n/2)]
  x3n.4 <- excess[ceiling(.75*n)]
  d <- xn.2^2 / (2 * xn.2 - x3n.4)
  
  shape <- -log(xn.2 / (x3n.4 - xn.2) ) / log(2)
  scale <- -shape * d
  
  if ( (max(excess) * shape) > -scale)
    message <- "Estimates are valid"
  
  else
    message <- "Estimates are not valid"
  
  estim <- param <- c(scale = scale, shape = shape)
  std.err <- var.cov <- corr <- NULL
  convergence <- counts <- NA
  var.thresh <- FALSE
  
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              nat = nat, pat = pat, convergence = convergence,
              corr = corr, counts = counts, exceed = exceed,
              scale = scale, var.thresh = var.thresh, est = "pickands"))
}


