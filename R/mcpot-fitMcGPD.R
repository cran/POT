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

fitmcgpd <- function (data, threshold, model = "log", start, ...,
                      std.err.type = "observed", corr = FALSE,
                      warn.inf = TRUE, method = "BFGS"){

  if (all(c("observed", "none") != std.err.type))
    stop("``std.err.type'' must be one of ``observed'' or ``none''")
  
  std.err.type <- match.arg(std.err.type, c("observed", "none"))
  model <- match.arg(model, c("log", "alog", "nlog", "anlog", "mix", "amix"))
  data <- as.double(data)
  threshold <- as.double(threshold)
  
  if (all(data<=threshold))
    stop("No data above threshold.")

  if (any(is.na(data))){
    warning("NAs present in data. Replacing them by the threshold.")
    data[is.na(data)] <- threshold
  }
 
  n <- as.integer(length(data))

  data1 <- data[-1]
  data2 <- data[-n]
  call <- match.call()
  
  idx1 <- data1>threshold
  exceed1 <- data1[idx1]
  nat1 <- as.integer(sum(idx1))
  pat1 <- nat1 / (n-1)
  
  idx2 <- data2>threshold
  exceed2 <- data2[idx2]
  nat2 <- as.integer(sum(idx2))
  pat2 <- nat2 / (n-1)

  nat <- as.integer(sum(idx1 & idx2))
  nat <- c(nat1, nat2, nat)
  pat <- c(pat1, pat2, nat[3]/(n-1))
  
  data3 <- data[2:(n-1)]
  idx3 <- data3 > threshold
  exceed3 <- data3[idx3]
  nat3 <- as.integer(sum(idx3))
  pat3 <- nat3 / (n - 2)
  
  ##Now reformat data to keep only realizations which at least one
  ##margin observation exceed the threshold.
  idx <- idx1 | idx2
  data1 <- data1[idx]
  data2 <- data2[idx]
  nn <- as.integer(sum(idx))
  nat <- c(nat, nn)
  
  param <- c("scale", "shape", "alpha", "asCoef1", "asCoef2", "asCoef")
  
  ##Creating suited negative log-likelihood according to the
  ##specified model
  if (model == "log"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmclog, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1))$dns
  }
  
  if (model == "alog"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmcalog, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef1,
        asCoef2, dns = double(1))$dns
    #param <- c(param, "asCoef1", "asCoef2")
  }
  
  if (model == "nlog"){
   nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmcnlog, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1))$dns
  }
  if (model == "anlog"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmcanlog, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef1,
        asCoef2, dns = double(1))$dns
    #param <- c(param, "asCoef1", "asCoef2")
  }
  if (model == "mix"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmcmix, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1))$dns
  }
  if (model == "amix"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
    -.C(POT_do_gpdmcamix, data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef,
        dns = double(1))$dns
    #param <- c(param, "asCoef")
  }    
  # never defined by M. Ribatet in C code
  # if (model == "amixtest"){
  #   nlpot <- function(scale, shape, alpha, asCoef1, asCoef2, asCoef)
  #   -.C(POT_do_gpdmcamixtest, data1, data2, exceed3, as.integer(n-1),
  #       as.integer(nn), as.integer(n-2), as.integer(nat3),
  #       pat3, threshold, scale, shape, alpha, asCoef,
  #       dns = double(1))$dns
  #   #param <- c(param, "asCoef")
  #   model <- "amix"
  # }    

  ##Creating suited starting values according to the chosen
  ##model (if needed) that is MLE estimates on marginal data
  if (missing(start)){
    start <- list(scale = 0, shape = 0)
    temp <- fitgpd(data, threshold, est = "pwmu")$param
    names(temp) <- NULL
    start$scale <- temp[1]
    start$shape <- temp[2]
            
    if (model == "log")
      start <- c(start, list(alpha = 0.75))
    if (model == "nlog")
      start <- c(start, list(alpha = 0.6))
    if (model == "alog")
      start <- c(start, list(alpha = 0.65, asCoef1 = 0.75,
                             asCoef2 = 0.75))
    if (model == "anlog")
      start <- c(start, list(alpha = 0.8, asCoef1 = 0.75,
                             asCoef2 = 0.75))
    if (model == "mix")
      start <- c(start, list(alpha = 0.25))
    if (model == "amix")
      start <- c(start, list(alpha = 0.75, asCoef = 0))
  }
  if (!is.list(start)) 
    stop("`start' must be a named list")
  
  #removed fixed parameters
  start <- start[!names(start) %in% names(list(...))]
  #get fixed parameters
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (any(!names(start) %in% param)) 
    stop("unspecified parameters in starting values")
  if (any(!names(fixed.param) %in% param)) 
    stop("unspecified parameters in fixed parameters values")
  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  
  nstart <- names(start)
  lstart <- length(nstart)
  f <- formals(nlpot)
  #names(f) <- param
  m <- match(nstart, param)
  m <- m[!is.na(m)]
  
  if (any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  #reorder parameters
  formals(nlpot) <- c(f[m], f[-m])
  
  nllh <- function(p, ...)
   nlpot(p, ...)
  
  if (lstart > 1) 
    body(nllh) <- parse(text = paste("nlpot(", paste("p[", 
                          1:lstart, "]", collapse = ", "),
                          ", ...)"))                                                    
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  
  if (warn.inf && do.call("nllh", start.arg) == 1e+06) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)

  if (opt$convergence != 0) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
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
      var.cov <- structure(solve(var.cov, tol = tol), dimnames = list(nstart,nstart))
      std.err <- diag(var.cov)
      names(std.err) <- nstart
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
                                dimnames = list(nstart,nstart))
          diag(corr.mat) <- rep(1, length(std.err))
        }else 
        {
          corr.mat <- NULL
        }
      }
    }
  }else # if(std.err.type == "none")
    std.err <- corr.mat <- var.cov <- NULL

  names.fixed <- names(fixed.param)
  fixed.param <- unlist(fixed.param)
  names(fixed.param) <- names.fixed
  
  param <- c(opt$par, unlist(fixed.param))
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, var.cov = var.cov,
                 fixed = unlist(fixed.param), param = param, deviance = 2*opt$value,
                 corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
                 message = opt$message, threshold = threshold, nat = nat3, pat = pat3,
                 data = data, exceed = exceed3, call = call, est = "MLE",
                 model = model, logLik = -opt$value, var.thresh = FALSE,
                 opt.value = -opt$value)

  chi <- 2 * (1 - pickdep(fitted, plot = FALSE)(0.5))
  fitted <- c(fitted, list(chi = chi))
  class(fitted) <- c("mcpot", "uvpot", "pot")
  return(fitted)
}

