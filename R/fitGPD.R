## In this file, several functions to estimates the GPD parameters
## are available:
##   1) Moments Estimator
##   2) Unbiased Probability Weighted Moment (PWMU) Estimator
##   3) Biased Probability Weighted Moment (PWMB) Estimator
##   4) Maximum Likelihood Estimator


## A generic function for estimate the GPD parameters
fitgpd <- function(data,threshold,method, ...){
  fitted <- switch(method, 'moments' = gpdmoments(data,threshold, ...),
                 'pwmb' = gpdpwmb(data, threshold, ...),
                 'pwmu' = gpdpwmu(data, threshold, ...),
                 'mle' = gpdmle(data, threshold, ...)
                 )
  printpot(fitted)
}

## Moments Estimator

gpdmoments <- function(data,threshold){

  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( length(exceed) == 0 )
    stop("None observation above the specified threshold !!!")

  exceed <- sort(exceed)
  
  loc <- threshold

  ## Evaluate the excess above the threshold 
  exces <- exceed - loc

  m <- mean(exces)
  v <- var(exces)

  scale <- m / 2 * ( m^2 / v +1 )
  shape <- - ( m^2 / v -1 ) / 2

  estim <- c(scale  = scale, shape = shape)
  param <-  c(scale = scale, shape =shape)
  convergence <- NA
  counts <- NA

  a11 <- 2*scale^2 * ( 1 - 6*shape + 12*shape^2)
  a12 <- - scale * (1-2*shape) * (1-4*shape+12*shape^2)
  a21 <- a12
  a22 <- (1-2*shape)^2 * (1-shape+6*shape^2)

  var.cov <- (1 - shape)^2 / ( (1-2*shape)*(1-3*shape)*(1-4*shape)*nat ) * matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )

  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)                    
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
  
  if ( shape > 0.25 ) message <- 'Assymptotic theory assumptions
for standard error may not be fullfilled !'
  else message <- NULL
  
  return(list(estimate = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              nhigh = nat, nat = nat, pat = pat, convergence = convergence,
              corr= corr, counts = counts, exceedances = exceed, scale=scale))
}

##PWMB Estimator

gpdpwmb <- function(data,threshold,a=0.35,b=0){

  exceed <- data[data>threshold]
  nat <- length( exceed )
  pat <- nat / length( data )
  
  if ( length(exceed) == 0 )
    stop("None observation above the specified threshold !!!")

  exceed <- sort(exceed)
  
  loc <- threshold
  
  exces <- exceed - loc

  m <- mean(exces)
  n <- length(exces)
  p <- 0

  for (i in 1:n){
    p[i] <- (i-a)/(n+b)
  }

  t <- sum((1-p)*exces)/n

  shape <- - m / (m- 2*t ) + 2
  scale <- 2 * m * t / (m - 2*t )
  type <- 'PWM'

  if (exces[n]<scale/shape){

    message <- 'Hybride'
    shape <- - scale / exces[n]

  }
  else message <- NULL
  
  estim <- c(scale  = scale, shape = shape)
  param <-  c(scale = scale, shape =shape)
  convergence <- NA
  counts <- NA

  a11 <- scale^2 * (7-18*shape+11*shape^2-2*shape^3)
  a12 <- - scale * (2-shape) * (2-6*shape+7*shape^2-2*shape^3)
  a21 <- a12
  a22 <- (1-shape) * (2 -shape)^2 * (1-shape+2*shape^2)

  var.cov <- 1 / ( (1-2*shape) * (3-2*shape)*nat ) * matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
      
  if ( shape > 0.5 ) message <- paste(message,"\n Assymptotic theory assumptions
for standard error may not be fullfilled !", sep='')
  
  return(list(estimate = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = type, threshold = threshold,
              corr = corr, convergence = convergence, counts = counts,
              nhigh = nat, nat = nat, pat = pat,
              exceedances = exceed, scale=scale))
}


## PWMU Estimator
## First, we need a function which computes the samples L-moments

samlmu <- function (x, nmom = 4, sort.data = TRUE)
{
    xok <- x[!is.na(x)]
    n <- length(xok)
    if (nmom <= 0) return(numeric(0))
    if (nmom <= 2) rnames <- paste("l", 1:nmom, sep = "_")
    else rnames <- c("l_1", "l_2", paste("t", 3:nmom, sep = "_"))
    lmom <- rep(NA, nmom)
    names(lmom) <- rnames
    if (n == 0) return(lmom)
    if (sort.data == TRUE) xok <- sort(xok)
    nmom.actual <- min(nmom, n)
    lmom[1] <- mean(xok)
    if (nmom.actual == 1) return(lmom)
    temp <- seq(1-n, n-1, by = 2)
    p1 <- rep(1, n)
    p <- temp/(n-1)
    lmom[2] <- mean(xok * p)
    if (nmom.actual == 2) return(lmom)
    if (xok[1] == xok[n]) {
        warning("all data values equal")
        return(lmom)
    }
    for (j in 3:nmom.actual) {
        p2 <- p1
        p1 <- p
        p <- ((2*j-3)*temp*p1 - (j-2)*(n+j-2)*p2) / ((j-1)*(n-j+1))
        lmom[j] <- mean(xok * p)/lmom[2]
    }
    return(lmom)
}

gpdpwmu <- function(data,threshold){

  exceed <- data[data>threshold]

  if ( length(exceed) == 0 )
    stop("None observation above the specified threshold !!!")

  exceed <- sort(unique(exceed))
  nat <- length( exceed )
  pat <- nat / length( data )
  
  loc <- threshold
  
  lmoments <- samlmu(exceed, nmom=2, sort.data = FALSE)
  shape <- - (lmoments[1] - loc)/lmoments[2] + 2
  scale <- (1 - shape)*(lmoments[1] - loc)
  names(shape) <- NULL
  names(scale) <- NULL

  estim <- c(scale  = scale, shape = shape)
  param <-  c(scale = scale, shape =shape)
  convergence <- NA
  counts <- NA
  a11 <- scale^2 * (7-18*shape+11*shape^2-2*shape^3)
  a12 <- - scale * (2-shape) * (2-6*shape+7*shape^2-2*shape^3)
  a21 <- a12
  a22 <- (1-shape) * (2 -shape)^2 * (1-shape+2*shape^2)

  var.cov <- 1 / ( (1-2*shape) * (3-2*shape)*nat ) * matrix(c(a11,a21,a12,a22),2)
  colnames(var.cov) <- c('scale','shape')
  rownames(var.cov) <- c('scale','shape')
  std.err <- sqrt( diag(var.cov) )
  
  .mat <- diag(1/std.err, nrow = length(std.err))
  corr <- structure(.mat %*% var.cov %*% .mat)                    
  diag(corr) <- rep(1, length(std.err))
  colnames(corr) <- c('scale','shape')
  rownames(corr) <- c('scale','shape')
      
  if ( shape > 0.5 ) message <- "Assymptotic theory assumptions
for standard error may not be fullfilled !"
  else message <- NULL
  
  return(list(estimate = estim, std.err = std.err, var.cov = var.cov,
              param = param, message = message, threshold = threshold,
              corr = corr, convergence = convergence, counts = counts,
              nhigh = nat, nat = nat, pat = pat,
              exceedances = exceed, scale=scale))
}

## The last two fucntions came from the evd package. The gpd.mle function
## corresponds to the fpot function. Nevertheless, it was sligthly modified
## to simplify it. So, this function is a ligther version of fpot.
## So, I'm very gratefull to Alec Stephenson.

gpdmle <- function(x, threshold, start,...,
                    std.err = TRUE, corr = FALSE,
                    method = "BFGS", warn.inf = TRUE){
  
  nlpot <- function(scale, shape) { 
    .C("nlgpd",
       exceed, nhigh, threshold, scale, shape, dns = double(1),
       PACKAGE = "POT")$dns
  }

  nn <- length(x)
    
  extind <- r <- NULL
  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])
  nhigh <- nat <- length(exceed)
    
  if(!nhigh) stop("no data above threshold")
  
  pat <- nat/nn
  param <- c("scale", "shape")
  
  if(missing(start)) {
    
    start <- list(scale = 0, shape = 0)
    start$scale <- mean(exceed) - threshold
   
    start <- start[!(param %in% names(list(...)))]
    
  }
  
  if(!is.list(start)) 
    stop("`start' must be a named list")
  
  if(!length(start))
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nlpot)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nlpot) <- c(f[m], f[-m])
  nllh <- function(p, ...) nlpot(p, ...)
  
  if(l > 1)
    body(nllh) <- parse(text = paste("nlpot(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)
  if(warn.inf && do.call("nllh", start.arg) == 1e6)
    warning("negative log-likelihood is infinite at starting values")
  
  opt <- optim(start, nllh, hessian = TRUE, ..., method = method)
  
  if (opt$convergence != 0) {
    warning("optimization may not have succeeded")
    if(opt$convergence == 1) opt$convergence <- "iteration limit reached"
  }
  
  else opt$convergence <- "successful"
  
  if(std.err) {
    tol <- .Machine$double.eps^0.5
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)) 
      stop("observed information matrix is singular; use std.err = FALSE")
    var.cov <- solve(var.cov, tol = tol)
    colnames(var.cov) <- nm
    std.err <- diag(var.cov)
    if(any(std.err <= 0))
      stop("observed information matrix is singular; use std.err = FALSE")
    std.err <- sqrt(std.err)
    names(std.err) <- nm
    if(corr) {
      .mat <- diag(1/std.err, nrow = length(std.err))
      corr <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
      diag(corr) <- rep(1, length(std.err))
    }
    else {
      corr <- NULL
      var.cov <- NULL
    }
  }
  
  else std.err <- corr <- var.cov <- NULL
  param <- c(opt$par, unlist(fixed.param))
  scale <- param["scale"]
  
  list(estimate = opt$par, std.err = std.err, var.cov = var.cov, fixed =
       unlist(fixed.param), param = param, deviance = 2*opt$value,
       corr = corr, convergence = opt$convergence, counts =
       opt$counts, message = opt$message, threshold = threshold, nhigh = nhigh, nat = nat, pat = pat, data = x, exceedances
       = exceed, scale = scale)
}

"printpot" <-  function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:", deparse(x$call), "\n")
    cat("Deviance:", x$deviance, "\n")

    cat("\nThreshold:", round(x$threshold, digits), "\n")
    cat("Number Above:", x$nat, "\n")
    cat("Proportion Above:", round(x$pat, digits), "\n")
    if(!is.null(x$extind)) {
      cat("\nClustering Interval:", x$r, "\n")
      if(is.finite(x$ulow)) {
        cat("Lower Threshold:", round(x$ulow, digits), "\n")
        cat("Lower Clustering Interval:", x$rlow, "\n")
      }
      cat("Number of Clusters:", x$nhigh, "\n")
      cat("Extremal Index:", round(x$extind, digits), "\n")
    }
    
    cat("\nEstimates\n") 
    print.default(format(x$estimate, digits = digits), print.gap = 2, 
        quote = FALSE)
    if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
    if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
    if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
    cat("\nOptimization Information\n")
    cat("  Convergence:", x$convergence, "\n")
    cat("  Function Evaluations:", x$counts["function"], "\n")
    if(!is.na(x$counts["gradient"]))
        cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
    if(!is.null(x$message)) cat("  Message:", x$message, "\n")
    cat("\n")
    invisible(x)
}
