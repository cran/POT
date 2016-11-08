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



##   3) Biased Probability Weighted Moment (PWMB) Estimator



##PWMB Estimator

gpdpwmb <- function(data, threshold, a=0.35, b=0, hybrid = FALSE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'pwmb'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( nat == 0 )
    stop("None observation above the specified threshold !!!")
  
  exceed <- sort(exceed)
  
  loc <- threshold
  
  excess <- exceed - loc
  
  m <- mean(excess)
  n <- length(excess)
  p <- (1:n - a) / (n + b)
  
  t <- sum((1-p)*excess)/n
  
  shape <- - m / (m- 2*t ) + 2
  scale <- 2 * m * t / (m - 2*t )
  est <- 'PWMB'

  if (hybrid)
    if ( (max(excess) >= (-scale / shape)) & (shape < 0) ){
      shape <- -scale / max(excess)
      est <- 'PWMB Hybrid'
    }
  
  estim <- c(scale  = scale, shape = shape)
  param <-  c(scale = scale, shape = shape)
  convergence <- NA
  counts <- NA
  
  a11 <- scale^2 * (7-18*shape+11*shape^2-2*shape^3)
  a12 <- - scale * (2-shape) * (2-6*shape+7*shape^2-2*shape^3)
  a21 <- a12
  a22 <- (1-shape) * (2 -shape)^2 * (1-shape+2*shape^2)
  
  var.cov <- 1 / ( (1-2*shape) * (3-2*shape)*nat ) *
    matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.5 )
    message <- "Assymptotic theory assumptions for standard error may not be fullfilled !"
  else message <- NULL
  
  var.thresh <- FALSE
  
  return(list(fitted.values = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              corr = corr, convergence = convergence, counts = counts,
              nat = nat, pat = pat, exceed = exceed,
              scale=scale, var.thresh = var.thresh, est = est))
}



