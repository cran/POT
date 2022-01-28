#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet   
#   Copyright (c) 2022 Christophe Dutang => add names
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

fitbvgpd <- function (data, threshold, model = "log", start, ...,
                      cscale = FALSE, cshape = FALSE,
                      std.err.type = "observed", corr = FALSE,
                      warn.inf = TRUE, method = "BFGS"){

  if (!(std.err.type %in% c("observed", "none")))
    stop("``std.err.type'' must be one of ``observed'' or ``none''")

  std.err.type <- match.arg(std.err.type, c("observed", "none"))
  threshold <- as.double(threshold)
  data1 <- as.double(data[,1])
  data2 <- as.double(data[,2])
  call <- match.call()
  
  n1 <- length(data1)
  idx1 <- (data1>threshold[1]) & !is.na(data1)
  exceed1 <- data1[idx1]
  nat1 <- sum(idx1)
  pat1 <- nat1 / (n1)
  
  if (!nat1)
    stop("No data above threshold for margin 1")

  n2 <- length(data2)
  idx2 <- (data2>threshold[2]) & !is.na(data2)
  exceed2 <- data2[idx2]
  nat2 <- sum(idx2)
  pat2 <- nat2 / (n2)

  if (!nat2)
    stop("No data above threshold for margin 2")
  
  n <- n1

  nat <- sum(idx1 & idx2)
  nat <- c(nat1, nat2, nat)
  pat <- c(pat1, pat2, nat[3]/n)
  
  
  if (any(is.na(data1))){
    warning("NAs present in data1. Replacing them by the threshold.")
    data1[is.na(data1)] <- threshold[1]
  }
  if (any(is.na(data2))){
    warning("NAs present in data2. Replacing them by the threshold.")
    data2[is.na(data2)] <- threshold[2]
  }

  ##Now reformat data to keep only realizations which at least one
  ##margin observation exceed the threshold.
  idx <- idx1 | idx2
  data1 <- data1[idx]
  data2 <- data2[idx]
  nn <- sum(idx)
  nat <- c(nat, nn)
  names(nat) <- c("Exceedance nb marg 1", "Exceedance nb marg 2",
                  "Exceedance nb both marg", "Exceedance nb any marg")
  names(pat) <- c("Exceedance prop marg 1", "Exceedance prop marg 2",
                  "Exceedance prop both marg")
  
  param <- c("scale1", "shape1")

  if(!cscale)
    param <- c(param, "scale2")

  if (!cshape)
    param <- c(param, "shape2")
  
  param <- c(param, "alpha")
  
  ##Creating suited starting values according to the chosen
  ##model (if needed) that is MLE estimates on marginal data
  if (missing(start)){
    start <- list(scale1 = 0, shape1 = 0)
    temp <- fitgpd(data1, threshold[1], est = "pwmu")$param
    start$scale1 <- temp[1]
    start$shape1 <- temp[2]

    temp <- fitgpd(data2, threshold[2], est = "pwmu")$param
    if (!cscale)
      start$scale2 <- temp[1]

    if (!cshape)
      start$shape2 <- temp[2]
        
    if (model == "log")
      start <- c(start, list(alpha = 0.75))
    if (model == "nlog")
      start <- c(start, list(alpha = 0.6))
    if (model == "alog"){
      start <- c(start, list(alpha = 0.65, asCoef1 = 0.75,
                             asCoef2 = 0.75))
      param <- c(param, "asCoef1", "asCoef2")
    }
    if (model == "anlog"){
      start <- c(start, list(alpha = 0.8, asCoef1 = 0.75,
                             asCoef2 = 0.75))
      param <- c(param, "asCoef1", "asCoef2")
    }
    if (model == "mix")
      start <- c(start, list(alpha = 0.25))
    if (model == "amix"){
      start <- c(start, list(alpha = 0.75, asCoef = 0))
      param <- c(param, "asCoef")
    }
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

  ##Creating suited negative log-likelihood according to the
  ##specified model
  if (model == "log")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C(POT_do_gpdbvlog, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, dns = double(1))$dns
    }
  if (model == "nlog")      
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C(POT_do_gpdbvnlog, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, dns = double(1))$dns
    }
  
  if (model == "alog")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C(POT_do_gpdbvalog, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef1,
          asCoef2, dns = double(1))$dns
    }
  
  if (model == "anlog")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
    
      -.C(POT_do_gpdbvanlog, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef1,
          asCoef2, dns = double(1))$dns
    }
   
  if (model == "mix")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C(POT_do_gpdbvmix, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha,
          dns = double(1))$dns
    }
  
  if (model == "amix")   
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2, asCoef){
       if (cscale) scale2 <- scale1
       if (cshape) shape2 <- shape1
    
      -.C(POT_do_gpdbvamix, data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef,
          dns = double(1))$dns
     }
  
  nstart <- names(start)
  lstart <- length(nstart)
  f <- formals(nlbvpot)
  f <- f[c(TRUE, TRUE, !cscale, !cshape, TRUE)]
  names(f) <- param
  m <- match(nstart, param)
  
  if (any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  #reorder parameters
  formals(nlbvpot) <- c(f[m], f[-m])
  nllh <- function(p, ...) 
    nlbvpot(p, ...)

  if (lstart > 1) 
    body(nllh) <- parse(text = paste("nlbvpot(", paste("p[", 
                          1:lstart, "]", collapse = ", "), ", ...)"))
  
  #fixed.param <- list(...)[names(list(...)) %in% param]

  #if (any(!(param %in% c(nstart, names(fixed.param))))) 
  #  stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  
  if (warn.inf && do.call("nllh", start.arg) == 1e+06) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)

  if ((opt$convergence != 0) || (opt$value == 1e+06)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"

  tol <- sqrt(.Machine$double.eps)

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

  
  param <- c(opt$par, unlist(fixed.param))
  
  fittedres <- list(fitted.values = opt$par, std.err = std.err, var.cov = var.cov,
                 fixed = unlist(fixed.param), param = param, deviance = 2*opt$value,
                 corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
                 message = opt$message, threshold = threshold, nat = nat, pat = pat,
                 data = data, exceed1 = exceed1, exceed2 = exceed2, call = call,
                 est = "MLE", model = model, logLik = -opt$value, opt.value = -opt$value)

  chi <- 2 * (1 - pickdep(fittedres, plot = FALSE)(0.5))
  fittedres <- c(fittedres, list(chi = chi))
  class(fittedres) <- c("bvpot", "pot")
  return(fittedres)
}
