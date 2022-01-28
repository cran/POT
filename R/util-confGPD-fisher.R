#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet                  
#   Copyright (c) 2022 Christophe Dutang => replace fitted to object
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




## Compute the confidence interval given by asymptotic theory
## i.e. the expected information matrix of Fisher (observed in the MLE case)
## for the shape parameter
gpd.fishape <- function(object, conf = 0.95){
  if (!inherits(object, "pot"))
    stop("Use only with 'pot' objects")

  nat <- object$nat
  se.shape <- object$std.err[2]
  shape.mle <- object$param[2]
  
  conf.inf <- shape.mle + qnorm( (1-conf) / 2 ) * se.shape
  conf.sup <- shape.mle - qnorm( (1-conf) / 2 ) * se.shape
  
  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )

  return(int.conf)
}

## Compute the confidence interval given by asymptotic theory
## i.e. the expected information matrix of Fisher (observed in the MLE case)
## for the scale parameter 
gpd.fiscale <- function(object, conf = 0.95){
  if (!inherits(object, "pot"))
    stop("Use only with 'pot' objects")
  
  nat <- object$nat
  se.scale <- object$std.err[1]
  scale.mle <- object$scale
  
  conf.inf <- scale.mle + qnorm( (1-conf) / 2 ) * se.scale
  conf.sup <- scale.mle - qnorm( (1-conf) / 2 ) * se.scale

  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )
  
  return(int.conf)
}

## Compute the confidence interval given by asymptotic theory and
## the Delta-Method for a specified return level.
gpd.firl <- function(object, prob, conf = 0.95){
  if (!inherits(object, "pot"))
    stop("Use only with 'pot' objects")
  
  scale.fit <- object$scale
  shape.fit <- object$param[2]
  threshold <- object$threshold
  rl.fit <- qgpd(prob, threshold, scale.fit, shape.fit)
  
  varcov <- object$var.cov
  if (is.null(varcov))
    stop("The correlation matrix should be present in object `object'!\n
Use `corr = TRUE' in `fitgpd' function.")
  diag(varcov) <- object$std.err^2
  
  eps <- .Machine$double.eps^0.5
  if ( abs(shape.fit) <= eps)
    grad.rl <- c(-log(1 - prob), 0)
  else
    grad.rl <- c(((1-prob)^(-shape.fit) - 1) / shape.fit,
                 -log(1 - prob) * scale.fit * (1 - prob)^(-shape.fit) / shape.fit -
                 ((1-prob)^(-shape.fit) - 1) * scale.fit / shape.fit^2 )
  
  var.rl <- t(grad.rl) %*% varcov %*% grad.rl
  
  conf.inf <- rl.fit + qnorm( (1-conf) / 2 ) * sqrt( var.rl )
  conf.sup <- rl.fit - qnorm( (1-conf) / 2 ) * sqrt( var.rl )
  
  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )
  
  return(int.conf)
}

confint.uvpot <- function(object, parm, level = 0.95, ...,
                          range, prob, prof = TRUE){

  if(!inherits(object, "uvpot"))
    stop("Use only with 'uvpot' objects")
  if (missing(parm)) 
    parm <- "quant"
  
  if (!(parm %in% c("quant","scale","shape")))
    stop("``parm'' must specify one of ``quant'', ``scale'' or ``shape''.")

  if ((parm == "quant") && missing(prob))
    stop("``prob'' must be specified when ``parm = 'quant'''.")
  
  if (prof){
    if (missing(range)){
      tmp <- confint(object, prof = FALSE, level = level,
                     parm = parm, prob = prob)
      range <- c(tmp[1] * 0.9, tmp[2] * 1.1)
    }
      
    ci <- switch(parm, "scale" = gpd.pfscale(object, range,
                        conf = level, ...), "shape" =
                 gpd.pfshape(object, range, conf = level, ...),
                 "quant" = gpd.pfrl(object, prob, range,
                   conf = level, ...))
  }

  else
    ci <- switch(parm, "scale" = gpd.fiscale(object, conf = level),
                 "shape" = gpd.fishape(object, conf = level),
                 "quant" = gpd.firl(object, prob, conf = level))

  return(ci)
}
