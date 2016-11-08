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



##   8) Likelihood Moment Estimator





##Likelihood moment estimation
gpdlme <- function(x, threshold, r = -.5, start, ...,
                   method = "BFGS"){

  nn <- length(x)
  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])

  nat <- length(exceed)
  if (!nat) 
    stop("no data above threshold")

  pat <- nat/nn
  
  excess <- exceed - threshold
  fun <- function(x){
    if (x >= 1/max(excess))
      return(1e6)
    p <- r / mean(log(1 - x * excess))
    abs(mean((1 - x * excess)^p) - 1 / (1 - r))
  }

   if (missing(start))
    start <- list(x = -1)
   
  opt <- optim(start, fun, hessian = TRUE, ..., method = method)

  if (opt$convergence != 0){
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  counts <- opt$counts
  b <- opt$par
  zero <- opt$value
  
  shape <- mean(log(1 - b*excess))
  scale <- - shape / b

  param <- c(scale, shape)
  names(param) <- c("scale", "shape")

  a11 <- scale^2 * (2 + ((r - shape)^2 + 2 * shape) /
                    (1 - 2 * r))
  a12 <- scale * (1 + (r^2 + shape^2 + shape) /
                  (1 - 2 * r))
  a22 <- (1 - r) * (1 + (2*shape^2 - 2 * shape + r) /
                    (1 - 2 * r))

  var.cov <- matrix(c(a11, a12, a12, a22), 2) / nat
  colnames(var.cov) <- c("scale", "shape")
  rownames(var.cov) <- c("scale", "shape")

  std.err <- sqrt(diag(var.cov))

  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c("scale", "shape")
  rownames(corr) <- c("scale", "shape")

  if (shape < -0.5) 
        message <- "Assymptotic theory assumptions\nfor standard error may not be fullfilled !"
    else message <- NULL
  
  var.thresh <- FALSE
  return(list(fitted.values = param, std.err = std.err, std.err.type = "expected",
              var.cov = var.cov, param = param, message = message, data = x,
              threshold = threshold, corr = corr, convergence = opt$convergence,
              counts = counts, nat = nat, pat = pat, exceed = exceed, scale = scale,
              var.thresh = var.thresh, est = "LME", opt.value = opt$value))
}



