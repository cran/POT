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


##   7) Method of Medians Estimator




##Medians estimation for the GPD ( Peng, L. and Welsh, A. (2002) )
gpdmed <- function(x, threshold, start, tol = 10^-3, maxit = 500,
                   show.trace = FALSE){
  
  if ( length(unique(threshold)) != 1){
    warning("Threshold must be a single numeric value for est = 'med'. Taking only the first value !!!")
    threshold <- threshold[1]
  }
  
  nn <- length(x)
  
  threshold <- rep(threshold, length.out = nn)
  
  high <- (x > threshold) & !is.na(x)
  threshold <- as.double(threshold[high])
  exceed <- as.double(x[high])
  nat <- length(exceed)
  
  excess <- exceed - threshold
  
  if(!nat) stop("no data above threshold")
  
  pat <- nat/nn
  
  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0.1)
    start["scale"] <- mean(exceed) - min(threshold)
    
  }

  start <- c(scale = start$scale, shape = start$shape)
  iter <- 1
  
  trace <- round(start, 3)
  
  ##Definition of a function to solve
  f <- function(x, y){
    -log(x)/y - (1+y)/y^2 * (1 - x^y) + log(x + .5)/y +
      (1+y)/y^2 * (1 - (x+.5)^y)
  }
  
  
  while (iter < maxit){
    ##If we have a non feasible point, we move back to feasible region
    if ( (start[2] < 0) & (max(excess) >= (-start[1] / start[2])))
      start[2] <- -start[1] / max(excess) + .1
    
    r1 <- start[2] * median(excess) / (2^start[2] - 1) - start[1]
    
    a <- log( 1 + start[2] * excess / start[1] ) / start[2]^2
    b <- (1 + start[2]) * excess / (start[1]*start[2] +
                                    start[2]^2 * excess)
    
    if (start[2] <= -1)
      y1 <- .5
    
    else{
      opt <- uniroot(f, c(10^-12, .5), y = start[2])
      y1 <- opt$root
    }
    
    r2 <- median(a - b) + log(y1)/start[2] + (1 + start[2]) /
      start[2]^2 * (1 - y1^start[2])
    
    next.point <- c(r1, r2) + start
    
    if (sqrt(sum( (next.point - start)^2) ) < tol)
      break
    
    trace <- rbind(trace, next.point)
    iter <- iter + 1
    start <- next.point
    
  }
  
  if(iter == maxit) opt$convergence <- "iteration limit reached"
  
  else opt$convergence <- "successful"
  
  opt$counts <- iter - 1
  names(opt$counts) <- "function"
  
  shape <- start[2]
  scale <- start[1]
  
  param <- c(scale = scale, shape = shape)
  names(param) <- c("scale", "shape")
  
  std.err <- std.err.type <- var.cov <- corr <- NULL
  
  var.thresh <- FALSE
  
  
  if (show.trace){
    if (iter >= 2)
      rownames(trace) <- c("Init. Val.", 1:(iter-1))

    print(round(trace, 3))
  }
  
  list(fitted.values = param, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = NULL, param = param,
       deviance = NULL, corr = corr, convergence = opt$convergence,
       counts = opt$counts, message = opt$message, threshold = threshold,
       nat = nat, pat = pat, data = x, exceed = exceed,
       scale = scale, var.thresh = var.thresh, est = "MEDIANS",
       opt.value = opt$f.root)
  
}

