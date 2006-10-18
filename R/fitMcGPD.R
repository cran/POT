fitmcgpd <- function (data, threshold, model = "log", start, ...,
                      obs.fish = TRUE, corr = FALSE,
                      warn.inf = TRUE, method = "BFGS"){

  data <- as.double(data)
  
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
  data3 <- data3[idx3]
  nat3 <- as.integer(sum(idx3))
  pat3 <- nat3 / (n - 2)
  
  ##Now reformat data to keep only realizations which at least one
  ##margin observation exceed the threshold.
  idx <- idx1 | idx2
  data1 <- data1[idx]
  data2 <- data2[idx]
  nn <- as.integer(sum(idx))
  nat <- c(nat, nn)
  
  param <- c("scale", "shape", "alpha")
  
  ##Creating suited starting values according to the chosen
  ##model (if needed) that is MLE estimates on marginal data
  if (missing(start)){
    start <- list(scale = 0, shape = 0)
    temp <- gpdmle(data, threshold)$param
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

  start <- start[!(param %in% names(list(...)))]

  if (!is.list(start)) 
    stop("`start' must be a named list")
  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  ##The marginal censored likelihood
  nlcpot <- function(scale, shape)
    -.C("cgpdlik", data3, as.integer(n - 2), as.integer(nat3),
        pat3, threshold, scale, shape, dns = double(1),
        PACKAGE = "POT")$dns
  
  ##Creating suited negative log-likelihood according to the
  ##specified model
  if (model == "log"){
    nlbvpot <- function(scale, shape, alpha)
      -.C("gpdbvlog", data1, data2, as.integer(n-1),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, dns = double(1),
          PACKAGE = "POT")$dns
    nllh.temp <- function(scale, shape, alpha){
      jllk <- nlbvpot(scale, shape, alpha)
      mjllk <- nlcpot(scale, shape)
      
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }
  
  if (model == "nlog"){
    nlbvpot <- function(scale, shape, alpha)
      -.C("gpdbvnlog", data1, data2, as.integer(n-1),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, dns = double(1),
          PACKAGE = "POT")$dns
    nllh.temp <- function(scale, shape, alpha){
      jllk <- nlbvpot(scale, shape, alpha)
      mjllk <- nlcpot(scale, shape)
      
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }
  
  if (model == "alog"){
    nlbvpot <- function(scale, shape, alpha, asCoef1, asCoef2)
      -.C("gpdbvalog", data1, data2, as.integer(n-1),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, asCoef1, asCoef2,
          dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef1", "asCoef2")
    nllh.temp <- function(scale, shape, alpha, asCoef1, asCoef2){
      jllk <- nlbvpot(scale, shape, alpha, asCoef1, asCoef2)
      mjllk <- nlcpot(scale, shape)
       
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }
  if (model == "anlog"){
    nlbvpot <- function(scale, shape, alpha, asCoef1, asCoef2)
      -.C("gpdbvanlog", data1, data2, as.integer(n),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, asCoef1, asCoef2,
          dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef1", "asCoef2")
    nllh.temp <- function(scale, shape, alpha, asCoef1, asCoef2){
      jllk <- nlbvpot(scale, shape, alpha, asCoef1, asCoef2)
      mjllk <- nlcpot(scale, shape)
      
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }
  if (model == "mix"){
    nlbvpot <- function(scale, shape, alpha)
      -.C("gpdbvmix", data1, data2, as.integer(n-1),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, dns = double(1),
          PACKAGE = "POT")$dns
    nllh.temp <- function(scale, shape, alpha){
      jllk <- nlbvpot(scale, shape, alpha)
      mjllk <- nlcpot(scale, shape)
      
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }  
  if (model == "amix"){
    nlbvpot <- function(scale, shape, alpha, asCoef)
      -.C("gpdbvamix", data1, data2, as.integer(n-1),
          as.integer(nn), pat1, pat2, c(threshold, threshold),
          scale, shape, scale, shape, alpha, asCoef,
          dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef")
    nllh.temp <- function(scale, shape, alpha, asCoef){
      jllk <- nlbvpot(scale, shape, alpha, asCoef)
      mjllk <- nlcpot(scale, shape)
      
      if ( (jllk == 1e6) || (mjllk == 1e6))
        return(1e6)
      
      return(jllk - mjllk)
    }
  }    

  nm <- names(start)
  l <- length(nm)
  f <- formals(nlbvpot)
  names(f) <- param
  m <- match(nm, param)
  
  if (any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlbvpot) <- c(f[m], f[-m])

  nllh <- function(p, ...)
   nllh.temp(p, ...)
  
  if (l > 1) 
    body(nllh) <- parse(text = paste("nllh.temp(", paste("p[", 
                          1:l, "]", collapse = ", "),
                          ", ...)"))                                                    
  
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

  if(obs.fish) {
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular.")
      obs.fish <- FALSE
      return
    }
    
    if (obs.fish){
      var.cov <- solve(var.cov, tol = tol)
      
      std.err <- diag(var.cov)
      if(any(std.err <= 0)){
        warning("observed information matrix is singular.")
        obs.fish <- FALSE
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

    if(!obs.fish)
      std.err <- corr.mat <- var.cov <- NULL
  }
  
  param <- c(opt$par, unlist(fixed.param))
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, var.cov = var.cov,
                 fixed = unlist(fixed.param), param = param, deviance = 2*opt$value,
                 corr = corr.mat, convergence = opt$convergence, counts = opt$counts,
                 message = opt$message, threshold = threshold, nat = nat3, pat = pat3,
                 data = data, exceed = data3, call = call,
                 type = "MLE", model = model, logLik = -opt$value, var.thresh = FALSE)

  chi <- 2 * (1 - pickdep(fitted, plot = FALSE)(0.5))
  fitted <- c(fitted, list(chi = chi))
  class(fitted) <- c("mcpot", "uvpot", "pot")
  return(fitted)
}
