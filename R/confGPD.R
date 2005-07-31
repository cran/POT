## This file contains several function to make (profile)
## confidence intervals on parameter of the GP distribution
## or on return level.

## Compute the profile confidence interval for the shape parameter
gpd.pfshape <- function(fitted, range, xlab, ylab,
                         conf = 0.95, nrang = 100,
                         vert.lines =TRUE, ...){

  cat('If there is some troubles try to put vert.lines = FALSE or change
 the range...\n')
  
  exceed<- fitted$exceedances
  threshold <- fitted$threshold
  nhigh <- fitted$nhigh

  ## First define a function who compute the profile log-likelihood
  ## for the shape parameter.
  gpd.plikshape <- function(scale){
    
    .C("nlgpd",
       exceed, nhigh, threshold, scale, shape, dns = double(1),
       PACKAGE = "POT")$dns
  }

  llik <- NULL
  int.shape <- seq(range[1], range[2], length = nrang)
  init <- fitted$scale
  
  for (shape in int.shape){
    opt <- optim(init, gpd.plikshape, method ="BFGS")
    param <- opt$par
    llik <- c(llik, -opt$value)
  }

  if ( missing(xlab) ) xlab <- 'Shape Parameter'
  if ( missing(ylab) ) ylab <- 'Profile Log-likelihood'
  
  plot(int.shape, llik, type='l', xlab = xlab, ylab = ylab, ...)
  
  llikmax <- - fitted$deviance / 2
  b.conf <- llikmax - 0.5 * qchisq(conf, 1)
  
  abline( h = llikmax)
  abline( h = b.conf)

  ## A special part to compute the bound of the profile likelihood
  ## confidence interval

  shape.mle <- fitted$param[2]
  index.neg1 <- which(llik <= b.conf & int.shape < shape.mle )
  index.neg1 <- max(index.neg1)

  index.neg2 <- which(llik <= b.conf & int.shape > shape.mle )
  index.neg2 <- min(index.neg2)

  index.pos1 <- which(llik >= b.conf & int.shape < shape.mle )
  index.pos1 <- min(index.pos1)

  index.pos2 <- which(llik >= b.conf & int.shape > shape.mle )
  index.pos2 <- max(index.pos2)

  conf.inf <- mean(int.shape[c(index.neg1,index.pos1)])
  conf.sup <- mean(int.shape[c(index.neg2,index.pos2)])
  
  if (vert.lines) abline(v = c(conf.inf, conf.sup) )
  
  return(c(conf.inf = conf.inf, conf.sup = conf.sup))
}

## Compute the profile confidence interval for the scale parameter
gpd.pfscale <- function(fitted, range, xlab, ylab,
                         conf = 0.95, nrang = 100,
                         vert.lines = TRUE, ...){

  cat('If there is some troubles try to put vert.lines = FALSE or change
 the range...\n')
  
  exceed<- fitted$exceedances
  threshold <- fitted$threshold
  nhigh <- fitted$nhigh

  ## First define a function who compute the profile log-likelihood
  ## for the scale parameter.
  gpd.plikscale <- function(shape){
    
    .C("nlgpd",
       exceed, nhigh, threshold, scale, shape, dns = double(1),
       PACKAGE = "POT")$dns
  }

  llik <- NULL
  int.scale <- seq(range[1], range[2], length = nrang)
  init <- fitted$param[2]
  
  for (scale in int.scale){
    opt <- optim(init, gpd.plikscale, method ="BFGS")
    param <- opt$par
    llik <- c(llik, -opt$value)
  }

  if ( missing(xlab) ) xlab <- 'Scale Parameter'
  if ( missing(ylab) ) ylab <- 'Profile Log-likelihood'
  
  plot(int.scale, llik, type='l', xlab = xlab, ylab = ylab, ...)
  
  llikmax <- - fitted$deviance / 2
  b.conf <- llikmax - 0.5 * qchisq(conf, 1)
  
  abline( h = llikmax)
  abline( h = b.conf)

  ## A special part to compute the bound of the profile likelihood
  ## confidence interval

  scale.mle <- fitted$scale
  index.neg1 <- which(llik <= b.conf & int.scale < scale.mle )
  index.neg1 <- max(index.neg1)

  index.neg2 <- which(llik <= b.conf & int.scale > scale.mle )
  index.neg2 <- min(index.neg2)

  index.pos1 <- which(llik >= b.conf & int.scale < scale.mle )
  index.pos1 <- min(index.pos1)

  index.pos2 <- which(llik >= b.conf & int.scale > scale.mle )
  index.pos2 <- max(index.pos2)

  conf.inf <- mean(int.scale[c(index.neg1,index.pos1)])
  conf.sup <- mean(int.scale[c(index.neg2,index.pos2)])

  if (vert.lines) abline(v = c(conf.inf, conf.sup) )
  
  return(c(conf.inf = conf.inf, conf.sup = conf.sup))
}

## Compute the profile confidence interval for the selected return level
gpd.pfrl <- function(fitted, retper, mu, range, xlab, ylab,
                         conf = 0.95, nrang = 100,
                         vert.lines = TRUE, ...){
  
  cat('If there is some troubles try to put vert.lines = FALSE or change
 the range...\n')

  ## First define a function who compute the profile log-likelihood
  ## for the return level. We need a reparametrization of the log-likelihood
  ## function.
  
  eps <- .Machine$double.eps^0.5
  gpd.plikrl <- function(shape){

    if ( abs(shape) < eps )
      scale <- (retlev - threshold) / log(1 - prob)
    else
      scale <- (retlev - threshold) * shape / ( (1 - prob)^(-shape) - 1 )
    .C("nlgpd",
       exceed, nhigh, threshold, scale, shape,
       dns = double(1), PACKAGE = "POT")$dns
  }

  exceed<- fitted$exceedances
  threshold <- fitted$threshold
  nhigh <- fitted$nhigh
  scale.fit <- fitted$scale
  shape.fit <- fitted$param[2]

  if ( range[1] <= threshold)
    stop("The lower bound on Return Level range is incompatible !")
  
  llik <- NULL
  int.retlev <- seq(range[1], range[2], length = nrang)
  retlev.fit <- qgpd(1 - 1 / (mu*retper), threshold, scale.fit, shape.fit)

  prob <- 1 - 1 / (mu * retper)
  for (retlev in int.retlev){
    opt <- optim(retlev.fit, gpd.plikrl, method ="BFGS")
    param <- opt$par
    llik <- c(llik, -opt$value)
  }

  if ( missing(xlab) ) xlab <- 'Return Level'
  if ( missing(ylab) ) ylab <- 'Profile Log-likelihood'
  
  plot(int.retlev, llik, type='l', xlab = xlab, ylab = ylab, ...)
  
  llikmax <- - fitted$deviance / 2
  b.conf <- llikmax - 0.5 * qchisq(conf, 1)
  
  abline( h = llikmax)
  abline( h = b.conf)

  ## A special part to compute the bound of the profile likelihood
  ## confidence interval

  index.neg1 <- which(llik <= b.conf & int.retlev < retlev.fit )
  index.neg1 <- max(index.neg1)

  index.neg2 <- which(llik <= b.conf & int.retlev > retlev.fit )
  index.neg2 <- min(index.neg2)

  index.pos1 <- which(llik >= b.conf & int.retlev < retlev.fit )
  index.pos1 <- min(index.pos1)

  index.pos2 <- which(llik >= b.conf & int.retlev > retlev.fit )
  index.pos2 <- max(index.pos2)

  conf.inf <- mean(int.retlev[c(index.neg1,index.pos1)])
  conf.sup <- mean(int.retlev[c(index.neg2,index.pos2)])

  if (vert.lines) abline(v = c(conf.inf, conf.sup) )
  
  return(c(conf.inf = conf.inf, conf.sup = conf.sup))
}

## Compute the confidence interval given by asymptotic theory
## i.e. the expected information matrix of Fisher (observed in the MLE case)
## for the shape parameter
gpd.fishape <- function(fitted, conf = 0.95){

  nhigh <- fitted$nhigh
  se.shape <- fitted$std.err[2]
  shape.mle <- fitted$param[2]
  
  conf.inf <- shape.mle + qnorm( (1-conf) / 2 ) * se.shape
  conf.sup <- shape.mle - qnorm( (1-conf) / 2 ) * se.shape

  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )

  return(int.conf)
}

## Compute the confidence interval given by asymptotic theory
## i.e. the expected information matrix of Fisher (observed in the MLE case)
## for the scale parameter 
gpd.fiscale <- function(fitted, conf = 0.95){

  nhigh <- fitted$nhigh
  se.scale <- fitted$std.err[1]
  scale.mle <- fitted$scale
  
  conf.inf <- scale.mle + qnorm( (1-conf) / 2 ) * se.scale
  conf.sup <- scale.mle - qnorm( (1-conf) / 2 ) * se.scale

  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )

  return(int.conf)
}

## Compute the confidence interval given by asymptotic theory and
## the Delta-Method for a specified return level.
gpd.firl <- function(fitted, retper, mu, conf = 0.95){

  prob <- 1 - 1 / (mu * retper)
  scale.fit <- fitted$scale
  shape.fit <- fitted$param[2]
  threshold <- fitted$threshold
  rl.fit <- qgpd(prob, threshold, scale.fit, shape.fit)

  varcov <- fitted$corr
  if (is.null(varcov))
    stop("The correlation matrix should be present in object `fitted'!\n
Use `corr = TRUE' in `fitgpd' function.")
  diag(varcov) <- fitted$std.err^2

  eps <- .Machine$double.eps^0.5
  if ( abs(shape.fit) <= eps)
    grad.rl <- c(-log( prob ), 0)
  else
    grad.rl <- c((prob^shape.fit - 1) / shape.fit,
                 - log(prob) * scale.fit / ( prob^shape.fit * shape.fit)
                 - prob^( ( 1 / shape.fit - 1 ) * scale.fit ) / shape.fit^2 )

  var.rl <- t(grad.rl) %*% varcov %*% grad.rl

  conf.inf <- rl.fit + qnorm( (1-conf) / 2 ) * sqrt( var.rl )
  conf.sup <- rl.fit - qnorm( (1-conf) / 2 ) * sqrt( var.rl )

  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )

  return(int.conf)
}
  
