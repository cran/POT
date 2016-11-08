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



## Compute the profile confidence interval for the shape parameter
gpd.pfshape <- function(fitted, range, xlab, ylab,
                        conf = 0.95, nrang = 100,
                        vert.lines = TRUE, ...){
  
  cat('If there is some troubles try to put vert.lines = FALSE or change
 the range...\n')
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")

  if (fitted$var.thresh)
    warning("Becarefull, you specify a varying threshold...\n")
  
  exceed<- fitted$exceed
  threshold <- fitted$threshold
  nat <- fitted$nat

  if (!fitted$var.thresh)
    threshold <- rep(threshold, nat)
  
  ## First define a function who compute the profile log-likelihood
  ## for the shape parameter.
  gpd.plikshape <- function(scale){
    -.C("gpdlik", exceed, nat, threshold, scale, shape, dns = double(1),
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
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")

  if (fitted$var.thresh)
    warning("Becarefull, you specify a varying threshold...\n")
  
  exceed<- fitted$exceed
  threshold <- fitted$threshold
  nat <- fitted$nat

  if (!fitted$var.thresh)
    threshold <- rep(threshold, nat)
  
  ## First define a function who compute the profile log-likelihood
  ## for the scale parameter.
  gpd.plikscale <- function(shape){
    -.C("gpdlik", exceed, nat, threshold, scale, shape, dns = double(1),
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
gpd.pfrl <- function(fitted, prob, range, thresh, xlab, ylab,
                     conf = 0.95, nrang = 100,
                     vert.lines = TRUE, ...){
  
  cat('If there is some troubles try to put vert.lines = FALSE or change
 the range...\n')
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")

  if (fitted$var.thresh)
    warning("Becarefull, you specify a varying threshold...\n")
  
  exceed <- fitted$exceed
  threshold <- fitted$threshold
  nat <- fitted$nat
  scale.fit <- fitted$scale
  shape.fit <- fitted$param[2]

  if (!fitted$var.thresh)
    threshold <- rep(threshold, nat)
  
  if (fitted$var.thresh & missing(thresh))
    stop("You must specify a particular threshold ``thresh'' when ``fitted'' has a varying threshold")
  
  if (missing(thresh))
    thresh <- threshold[1]
  
  if (!any(threshold == thresh))
    warning("``thresh'' isn't a particular threshold value in ``fitted''...")
  
  if ( range[1] <= thresh)
    stop("The lower bound on Return Level range is incompatible !")
    
  ## First define a function who compute the profile log-likelihood
  ## for the return level. We need a reparametrization of the log-likelihood
  ## function.
  
  eps <- .Machine$double.eps^0.5
  
  
  gpd.plikrl <- function(shape){
    if ( abs(shape) < eps )
      scale <- (retlev - thresh) / log(1 - prob)
    else
      scale <- (retlev - thresh) * shape / ( (1 - prob)^(-shape) - 1 )
    -.C("gpdlik", exceed, nat, threshold, scale, shape,
        dns = double(1), PACKAGE = "POT")$dns
  }
   
  llik <- NULL
  int.retlev <- seq(range[1], range[2], length = nrang)
  retlev.fit <- qgpd(prob, thresh, scale.fit, shape.fit)
  
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

