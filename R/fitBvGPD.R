fitbvgpd <- function (data, threshold, model = "log", start, ...,
                      cscale = FALSE, cshape = FALSE,
                      std.err.type = "observed", corr = FALSE,
                      warn.inf = TRUE, method = "BFGS"){

  if (all(c("observed", "none") != std.err.type))
    stop("``std.err.type'' must be one of ``observed'' or ``none''")

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

  start <- start[!(param %in% names(list(...)))]

  if (!is.list(start)) 
    stop("`start' must be a named list")
  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  ##Creating suited negative log-likelihood according to the
  ##specified model
  if (model == "log")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C("gpdbvlog", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, dns = double(1),
          PACKAGE = "POT")$dns
    }
  if (model == "nlog")      
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C("gpdbvnlog", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, dns = double(1),
          PACKAGE = "POT")$dns
    }
  
  if (model == "alog")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C("gpdbvalog", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef1,
          asCoef2, dns = double(1), PACKAGE = "POT")$dns
    }
  
  if (model == "anlog")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef1, asCoef2){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
    
      -.C("gpdbvanlog", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef1,
          asCoef2, dns = double(1), PACKAGE = "POT")$dns
    }
   
  if (model == "mix")
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha){
      if (cscale) scale2 <- scale1
      if (cshape) shape2 <- shape1
      
      -.C("gpdbvmix", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha,
          dns = double(1), PACKAGE = "POT")$dns
    }
  
  if (model == "amix")   
    nlbvpot <- function(scale1, shape1, scale2, shape2, alpha,
                        asCoef){
       if (cscale) scale2 <- scale1
       if (cshape) shape2 <- shape1
    
      -.C("gpdbvamix", data1, data2, n, nn, pat1, pat2, threshold,
          scale1, shape1, scale2, shape2, alpha, asCoef,
          dns = double(1), PACKAGE = "POT")$dns
     }
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nlbvpot)
  f <- f[c(TRUE, TRUE, !cscale, !cshape, TRUE)]
  names(f) <- param
  m <- match(nm, param)
  
  if (any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlbvpot) <- c(f[m], f[-m])
  nllh <- function(p, ...) nlbvpot(p, ...)

  if (l > 1) 
    body(nllh) <- parse(text = paste("nlbvpot(", paste("p[", 
                          1:l, "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]

  if (any(!(param %in% c(nm, names(fixed.param))))) 
    stop("unspecified parameters")
  
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
      return
    }
    
    if (std.err.type == "observed"){
      var.cov <- solve(var.cov, tol = tol)
      
      std.err <- diag(var.cov)
      if(any(std.err <= 0)){
        warning("observed information matrix is singular.")
        std.err.type <- "none"
      }

      else{
        std.err <- sqrt(std.err)
        
        if(corr) {
          .mat <- diag(1/std.err, nrow = length(std.err))
          corr.mat <- structure(.mat %*% var.cov %*% .mat,
                                dimnames = list(nm,nm))
          diag(corr.mat) <- rep(1, length(std.err))
        }
        else {
          corr.mat <- NULL
        }
        
        colnames(var.cov) <- nm
        rownames(var.cov) <- nm
        names(std.err) <- nm
      }
    }

    if(std.err.type == "none")
      std.err <- corr.mat <- var.cov <- NULL
  }
  
  param <- c(opt$par, unlist(fixed.param))
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, var.cov = var.cov,
                 fixed = unlist(fixed.param), param = param, deviance = 2*opt$value,
                 corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
                 message = opt$message, threshold = threshold, nat = nat, pat = pat,
                 data = data, exceed1 = exceed1, exceed2 = exceed2, call = call,
                 est = "MLE", model = model, logLik = -opt$value)

  chi <- 2 * (1 - pickdep(fitted, plot = FALSE)(0.5))
  fitted <- c(fitted, list(chi = chi))
  class(fitted) <- c("bvpot", "pot")
  return(fitted)
}
