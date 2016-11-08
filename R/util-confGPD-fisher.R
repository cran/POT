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




## Compute the confidence interval given by asymptotic theory
## i.e. the expected information matrix of Fisher (observed in the MLE case)
## for the shape parameter
gpd.fishape <- function(fitted, conf = 0.95){
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")

  nat <- fitted$nat
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
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")
  
  nat <- fitted$nat
  se.scale <- fitted$std.err[1]
  scale.mle <- fitted$scale
  
  conf.inf <- scale.mle + qnorm( (1-conf) / 2 ) * se.scale
  conf.sup <- scale.mle - qnorm( (1-conf) / 2 ) * se.scale

  int.conf <- c( conf.inf = conf.inf, conf.sup = conf.sup )
  
  return(int.conf)
}

## Compute the confidence interval given by asymptotic theory and
## the Delta-Method for a specified return level.
gpd.firl <- function(fitted, prob, conf = 0.95){
  if (!inherits(fitted, "pot"))
    stop("Use only with 'pot' objects")
  
  scale.fit <- fitted$scale
  shape.fit <- fitted$param[2]
  threshold <- fitted$threshold
  rl.fit <- qgpd(prob, threshold, scale.fit, shape.fit)
  
  varcov <- fitted$var.cov
  if (is.null(varcov))
    stop("The correlation matrix should be present in object `fitted'!\n
Use `corr = TRUE' in `fitgpd' function.")
  diag(varcov) <- fitted$std.err^2
  
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
