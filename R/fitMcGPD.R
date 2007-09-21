fitmcgpd <- function (data, threshold, model = "log", start, ...,
                      std.err.type = "observed", corr = FALSE,
                      warn.inf = TRUE, method = "BFGS"){

  if (all(c("observed", "none") != std.err.type))
    stop("``std.err.type'' must be one of ``observed'' or ``none''")
  
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
  
  param <- c("scale", "shape", "alpha")
  
  ##Creating suited negative log-likelihood according to the
  ##specified model
  if (model == "log"){
    nlpot <- function(scale, shape, alpha)
    -.C("gpdmclog", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1), PACKAGE = "POT")$dns
  }
  
  if (model == "alog"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2)
    -.C("gpdmcalog", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef1,
        asCoef2, dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef1", "asCoef2")
  }
  
  if (model == "nlog"){
   nlpot <- function(scale, shape, alpha)
    -.C("gpdmcnlog", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1), PACKAGE = "POT")$dns
  }
  if (model == "anlog"){
    nlpot <- function(scale, shape, alpha, asCoef1, asCoef2)
    -.C("gpdmcanlog", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef1,
        asCoef2, dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef1", "asCoef2")
  }
  if (model == "mix"){
    nlpot <- function(scale, shape, alpha)
    -.C("gpdmcmix", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha,
        dns = double(1), PACKAGE = "POT")$dns
  }
  if (model == "amix"){
    nlpot <- function(scale, shape, alpha, asCoef)
    -.C("gpdmcamix", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef,
        dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef")
  }    
  if (model == "amixtest"){
    nlpot <- function(scale, shape, alpha, asCoef)
    -.C("gpdmcamixtest", data1, data2, exceed3, as.integer(n-1),
        as.integer(nn), as.integer(n-2), as.integer(nat3),
        pat3, threshold, scale, shape, alpha, asCoef,
        dns = double(1), PACKAGE = "POT")$dns
    param <- c(param, "asCoef")
    model <- "amix"
  }    

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

  start <- start[!(param %in% names(list(...)))]

  if (!is.list(start)) 
    stop("`start' must be a named list")
  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  nm <- names(start)
  l <- length(nm)
  f <- formals(nlpot)
  names(f) <- param
  m <- match(nm, param)
  
  if (any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlpot) <- c(f[m], f[-m])

  nllh <- function(p, ...)
   nlpot(p, ...)
  
  if (l > 1) 
    body(nllh) <- parse(text = paste("nlpot(", paste("p[", 
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
  }

  if(std.err.type == "none")
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
                 model = model, logLik = -opt$value, var.thresh = FALSE)

  chi <- 2 * (1 - pickdep(fitted, plot = FALSE)(0.5))
  fitted <- c(fitted, list(chi = chi))
  class(fitted) <- c("mcpot", "uvpot", "pot")
  return(fitted)
}

dexi <- function(x, n.sim = 1000, n.mc = length(x$data),
                 plot = TRUE, ...){
  
  thresh <- x$threshold
  scale <- x$param["scale"]
  shape <- x$param["shape"]
  alpha <- x$param["alpha"]
  pat <- x$pat
  model <- x$model
  
  scale.new <- shape * thresh / (pat^(-shape) - 1)

  if (model %in% c("log", "nlog"))
    param <- list(alpha = alpha)

  if (model %in% c("alog", "anlog"))
    param <- list(alpha = alpha, asCoef1 = x$param["asCoef1"],
                  asCoef2 = x$param["asCoef2"])

  if (model == "mix")
    param <- list(alpha = alpha, asCoef = 0)
  
  if (model == "amix")
    param <- list(alpha = alpha, asCoef = x$param["asCoef"])

  param <- c(param, list(n = n.mc, model = model))

  exi <- rep(0, n.sim)
  mc <- rep(0, n.mc)
  
  for (i in 1:n.sim){
    mc <- do.call("simmc", param)
    mc <- qgpd(mc, 0, scale.new, shape)

    while(sum(mc > thresh) < 2){
      mc <- do.call("simmc", param)
      mc <- qgpd(mc, 0, scale.new, shape)
    }
    
    exi[i] <- fitexi(mc, thresh)$exi
  }

  if (plot)
    plot(density(exi, bw = sd(exi) /2), ...)
  
  invisible(exi)
}
