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



##   9) Maximum Goodness-of-Fit Estimator



## Maximum goodness-of-fit estimator
gpdmgf <- function(x, threshold, start, stat, ...,
                   method = "BFGS", warn.inf = TRUE){

  nn <- length(x)
  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])

  nat <- length(exceed)
  if (!nat) 
    stop("no data above threshold")

  if (!(stat %in% c("KS","CM","AD","ADR","ADL","AD2R","AD2L",
                   "AD2")))
    stop("`stat' must be one of 'KS','CM','AD','ADR','ADL','AD2R','AD2L', 'AD2'.")

  pat <- nat/nn
  param <- c("scale", "shape")

  excess <- exceed - threshold

  excess <- sort(excess)

  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0)
    start$scale <- mean(exceed) - min(threshold)
    
    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  if (stat == "KS")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else
        1 / 2 / nat + max(abs(pgpd(excess, 0, scale, shape) -
                              ppoints(nat)))
    }
    
  if (stat == "CM")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else
        1/12/nat + sum((pgpd(excess, 0, scale, shape) -
                        ppoints(nat))^2)
    }

  if (stat == "AD")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if ((pgpd(max(excess),  0, scale, shape) >= 1) ||
            (pgpd(min(excess),  0, scale, shape) <= 0))
         1e6

        else
          -nat - mean((2*1:nat - 1) *
                      (log(pgpd(excess, 0, scale, shape)) +
                       log(1 - pgpd(rev(excess), 0, scale, shape))))
      }
    }

  if (stat == "ADR")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(max(excess),  0, scale, shape) >= 1)
         1e6

        else
          nat / 2 - 2 * sum(pgpd(excess, 0, scale, shape)) -
            mean((2 * 1:nat -1)*log(1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }
  
  if (stat == "ADL")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(min(excess),  0, scale, shape) <= 0)
         1e6

        else
          - 3 * nat / 2 + 2 * sum(pgpd(excess, 0, scale, shape)) -
            mean((2 * 1:nat -1)*log(pgpd(excess, 0, scale, shape)))
      }
    }

  if (stat == "AD2R")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(max(excess),  0, scale, shape) >= 1)
         1e6

        else
          2 * sum(log(1 - pgpd(excess, 0, scale, shape))) +
            mean((2* 1:nat - 1) / (1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }

  if (stat == "AD2L")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if (pgpd(min(excess),  0, scale, shape) <= 0)
         1e6

        else
          2 * sum(log(pgpd(excess, 0, scale, shape))) +
            mean((2 * 1:nat - 1) / pgpd(excess, 0, scale, shape))
      }
    }

  if (stat == "AD2")
    fun <- function(scale, shape){
      if (scale <= 0)
        1e6

      else{

        if ((pgpd(max(excess),  0, scale, shape) >= 1) ||
            (pgpd(min(excess),  0, scale, shape) <= 0))
         1e6

        else
          2 * sum(log(pgpd(excess, 0, scale, shape)) +
                  log(1 - pgpd(excess, 0, scale, shape))) +
                    mean((2 * 1:nat - 1) / pgpd(excess, 0, scale, shape) +
                         (2 * 1:nat - 1) /
                         (1 - pgpd(rev(excess), 0, scale, shape)))
      }
    }

  nm <- names(start)
  l <- length(nm)
  f <- formals(fun)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(fun) <- c(f[m], f[-m])
  mgf <- function(p, ...) fun(p, ...)
  
  if(l > 1)
    body(mgf) <- parse(text = paste("fun(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if( warn.inf && do.call("mgf", start.arg) == 1e6 )
    warning("Maximum goodness-of-fit function is infinite at starting values")
  
  opt <- optim(start, mgf, hessian = TRUE, ..., method = method)
  
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  tol <- .Machine$double.eps^0.5
    
  param <- c(opt$par, unlist(fixed.param))
  scale <- param["scale"]
  
  var.thresh <- !all(threshold == threshold[1])

  if (!var.thresh)
    threshold <- threshold[1]

  std.err <- std.err.type <- corr.mat <- NULL
  var.cov <- NULL

  list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
       var.cov = var.cov, fixed = unlist(fixed.param), param = param,
       corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
       message = opt$message, threshold = threshold, nat = nat, pat = pat,
       data = x, exceed = exceed, scale = scale, var.thresh = var.thresh,
       est = "MGF", opt.value = opt$value, stat = stat)
}


