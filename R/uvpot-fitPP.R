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

fitpp <- function(data, threshold, noy = length(data) / 365.25, start, ...,
                  std.err.type = "observed", corr = FALSE,
                  method = "BFGS", warn.inf = TRUE){

  if (all(c("observed", "none") != std.err.type))
    stop("``std.err.type'' must be one of 'observed' or 'none'")
  
  std.err.type <- match.arg(std.err.type, c("observed", "none"))
  nlpp <- function(loc, scale, shape) { 
    -.C(POT_do_pplik, exceed, nat, loc, scale, shape,
        threshold, noy, dns = double(1))$dns
  }

  noy <- as.double(noy)
  nn <- length(data)
  
  high <- (data > threshold) & !is.na(data)
  threshold <- as.double(threshold)
  exceed <- as.double(data[high])
  nat <- as.integer(length(exceed))

  if(!nat)
    stop("no data above threshold")
  
  pat <- nat/nn
  param <- c("loc", "scale", "shape")
  
  if(missing(start)) {
    
    start <- list(loc = 0, scale = 0, shape = 0)
    start$scale <- sqrt(6 * var(exceed))/pi
    start$loc <- mean(exceed) + (log(noy) - 0.58) * start$scale    

    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nlpp)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlpp) <- c(f[m], f[-m])
  nllh <- function(p, ...) nlpp(p, ...)
  
  if(l > 1)
    body(nllh) <- parse(text = paste("nlpp(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if( warn.inf && do.call("nllh", start.arg) == 1e6 )
    warning("negative log-likelihood is infinite at starting values")
  
  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
  
  if ((opt$convergence != 0) || (opt$value == 1e6)) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  tol <- .Machine$double.eps^0.5
  
  if(std.err.type == "observed") {
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular.")
      std.err.type <- "none"
      std.err <- corr.mat <- var.cov <- NULL
    }else
    {
      var.cov <- structure(solve(var.cov, tol = tol), dimnames = list(nm,nm))
      std.err <- diag(var.cov)
      names(std.err) <- nm
      if(any(std.err <= 0)){
        warning("observed information matrix has non positive diagonal terms.")
        std.err.type <- "none"
        std.err <- corr.mat <- var.cov <- NULL
      }else
      {
        std.err <- sqrt(std.err)
        if(corr) 
        {
          .mat <- diag(1/std.err, nrow = length(std.err))
          corr.mat <- structure(.mat %*% var.cov %*% .mat,
                                dimnames = list(nm,nm))
          diag(corr.mat) <- rep(1, length(std.err))
        }else 
        {
          corr.mat <- NULL
        }
      }
    }
  }else # if(std.err.type == "none")
    std.err <- corr.mat <- var.cov <- NULL
  
  param <- c(opt$par, unlist(fixed.param))

  ##Transform the point process parameter to the GPD ones
  scale <- param["scale"]  + param["shape"] *
    (threshold - param["loc"])
  
  var.thresh <- !all(threshold == threshold[1])

  if (!var.thresh)
    threshold <- threshold[1]
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, threshold = threshold,
                 nat = nat, pat = pat, data = data, exceed = exceed, scale = scale,
                 var.thresh = var.thresh, est = "MLE", logLik = -opt$value,
                 opt.value = opt$value)

  fitted$threshold.call <- deparse(threshold)
  class(fitted) <- c("uvpot","pot")
  return(fitted)
}
