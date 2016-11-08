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


##   1) Moments Estimator


## Moments Estimator

gpdmoments <- function(data, threshold){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'moments'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( nat == 0 )
    stop("None observation above the specified threshold !!!")
  
  exceed <- sort(exceed)
  
  loc <- threshold
  
  ## Evaluate the excess above the threshold 
  exces <- exceed - loc
  
  m <- mean(exces)
  v <- var(exces)
  
  scale <- m / 2 * ( m^2 / v +1 )
  shape <- - ( m^2 / v -1 ) / 2
  
  estim <- param <- c(scale  = scale, shape = shape)
  convergence <- counts <- NA
  
  a11 <- 2*scale^2 * ( 1 - 6*shape + 12*shape^2)
  a12 <- - scale * (1-2*shape) * (1-4*shape+12*shape^2)
  a21 <- a12
  a22 <- (1-2*shape)^2 * (1-shape+6*shape^2)
  
  var.cov <- (1 - shape)^2 / ( (1-2*shape)*(1-3*shape)*(1-4*shape)*nat ) *
    matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)                    
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.25 ) message <- 'Assymptotic theory assumptions
for standard error may not be fullfilled !'
  else message <- NULL
  
  var.thresh <- FALSE
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              nat = nat, pat = pat, convergence = convergence,
              corr= corr, counts = counts, exceed = exceed,
              scale=scale, var.thresh = var.thresh, est = "moments"))
}


