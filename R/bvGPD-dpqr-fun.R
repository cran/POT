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

##A generic function to simulate from a bivariate extreme value
##distribution
rbvgpd <- function(n, alpha, model = "log", asCoef, asCoef1,
                   asCoef2, mar1 = c(0,1,0), mar2 = mar1){
  if (missing(alpha))
    stop("``alpha must be present.")

  if ((model %in% c("alog", "anlog")) && (missing(asCoef1) ||
                                          missing(asCoef2)))
    stop("``asCoef1'' and ``asCoef2'' must be present.")

  if ((model == "amix") && missing(asCoef))
    stop("``asCoef'' must be present.")

  if ((model %in% c("log", "alog")) && ((alpha <= 0) ||
                                        (alpha > 1)))
    stop("``alpha'' must be in ]0, 1] with this model.")

  if ((model %in% c("nlog", "anlog")) && (alpha <= 0))
    stop("``alpha'' must be positive with this model.")

  if ((model == "mix") && ((alpha < 0) || (alpha > 1)))
    stop("``alpha'' must be in [0,1] with this model.")

  if ((model == "amix") &&
      ((alpha < 0) || (alpha + 2 * asCoef >1) ||
       (alpha + 3 * asCoef < 0)))
    stop("``alpha'' and ``asCoef'' are not valid. See the doc.")

  if ((model == "amix")){
    alpha <- alpha + 3 * asCoef
    asCoef <- - asCoef
  }

  alpha <- as.double(alpha)

  switch(model,
         log = rbvlog(n = n, alpha = alpha, mar1 = mar1,
           mar2 = mar2),
         alog = rbvalog(n = n, alpha = alpha, asCoef1 = asCoef1,
           asCoef2 = asCoef2, mar1 = mar1, mar2 = mar2),
         nlog = rbvnlog(n = n, alpha = alpha, mar1 = mar1,
           mar2 = mar2),
         anlog = rbvanlog(n = n, alpha = alpha, asCoef1 = asCoef1,
           asCoef2 = asCoef2, mar1 = mar1, mar2 = mar2),
         mix = rbvmix(n = n, alpha = alpha, mar1 = mar1,
           mar2 = mar2),
         amix = rbvamix(n = n, alpha = alpha, asCoef = asCoef,
           mar1 = mar1, mar2 = mar2)
         )
}

##Uses Algorithm 1.1 in Stephenson(2003)
rbvlog <- function(n, alpha, mar1 = c(0,1,0), mar2 = mar1){
  
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0 ||
     alpha > 1) stop("invalid argument for `alpha'")
  
  sim <- .C("rbvlog_shi", as.integer(n), as.double(alpha),
            sim = double(2*n), PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has unit Frechet margins, thus
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])

  return(sim)
}

## Uses Algorithm 1.2 in Stephenson(2003)
rbvalog <- function(n, alpha, asCoef1, asCoef2, mar1 = c(0,1,0),
                    mar2 = mar1){
  asy <- c(asCoef1, asCoef2)
  
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0 ||
     alpha > 1) stop("invalid argument for `alpha'")
  
  if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
     max(asy) > 1) stop("invalid argument for `asy'")
  
  if(alpha == 1 || any(asy == 0)) {
    asy <- c(0,0)
    alpha <- 1
  }
  
  sim <- .C("rbvalog_shi", as.integer(n), as.double(alpha),
            as.double(asy), sim = double(2*n),
            PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has unit Frechet margins, thus
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])

  return(sim)
}

rbvnlog <- function(n, alpha, mar1 = c(0,1,0), mar2 = mar1){
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0)
    stop("invalid argument for `alpha'")
  
  sim <- .C("rbvnlog", as.integer(n), as.double(alpha),
            sim = runif(2*n), PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has uniform(0,1) margins, thus convert them to unit
  ##Frechet ones
  sim <- - 1 / log(sim)
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])

  return(sim)
}

rbvanlog <- function(n, alpha, asCoef1, asCoef2, mar1 = c(0,1,0),
                     mar2 = mar1){
  asy <- c(asCoef1, asCoef2)
  
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0)
    stop("invalid argument for `alpha'")
  
  sim <- .C("rbvanlog", as.integer(n), as.double(alpha),
            as.double(asy), sim = runif(2*n), PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has uniform(0,1) margins, thus convert them to unit
  ##Frechet ones
  sim <- - 1 / log(sim)
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])
  return(sim)
}

rbvmix <- function(n, alpha, mar1 = c(0,1,0), mar2 = mar1){
  if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for ``alpha''")

  if (alpha < 0)
    stop("``alpha'' must be non-negative")

  sim <- .C("rbvmix", as.integer(n), as.double(alpha),
            sim = runif(2*n), PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has uniform(0,1) margins, thus convert them to unit
  ##Frechet ones
  sim <- - 1 / log(sim)
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])

  return(sim)
}

rbvamix <- function(n, alpha, asCoef, mar1 = c(0,1,0),
                    mar2 = mar1){
  if(length(alpha) != 1 || mode(alpha) != "numeric")
        stop("invalid argument for ``alpha''")

  if (alpha < 0)
    stop("``alpha'' must be non-negative")

  if((alpha + asCoef) > 1)
    stop("``alpha'' + ``asCoef'' cannot be greater than one")
  
  if((alpha + 2*asCoef) > 1)
    stop("`alpha' + `2*asCoef' cannot be greater than one")
  
  if((alpha + 3*asCoef) < 0)
    stop("`alpha' + `3*asCoef' must be non-negative")
  
  sim <- .C("rbvamix", as.integer(n), as.double(alpha),
            as.double(asCoef), sim = runif(2*n),
            PACKAGE = "POT")$sim
  sim <- matrix(sim, nrow = n, ncol = 2, byrow = TRUE)

  ##sim has uniform(0,1) margins, thus convert them to unit
  ##Frechet ones
  sim <- - 1 / log(sim)
  sim[,1] <- frech2gpd(sim[,1], mar1[1], mar1[2], mar1[3])
  sim[,2] <- frech2gpd(sim[,2], mar2[1], mar2[2], mar2[3])

  return(sim)
}

##A generic function to compute (non) exceedande probabilities for
##bivariate POT models
pbvgpd <- function(q, alpha, model = "log", asCoef, asCoef1,
                   asCoef2, mar1 = c(0,1,0), mar2 = mar1,
                   lower.tail = TRUE){
  if (missing(alpha))
    stop("``alpha must be present.")

  if ((model %in% c("alog", "anlog")) && (missing(asCoef1) ||
                                          missing(asCoef2)))
    stop("``asCoef1'' and ``asCoef2'' must be present.")

  if ((model == "amix") && missing(asCoef))
    stop("``asCoef'' must be present.")

  if ((model %in% c("log", "alog")) && ((alpha <= 0) ||
                                        (alpha > 1)))
    stop("``alpha'' must be in ]0, 1] with this model.")

  if ((model %in% c("nlog", "anlog")) && (alpha <= 0))
    stop("``alpha'' must be positive with this model.")

  if ((model == "mix") && ((alpha < 0) || (alpha > 1)))
    stop("``alpha'' must be in [0,1] with this model.")

  if ((model == "amix") &&
      ((alpha < 0) || (alpha + 2 * asCoef >1) ||
       (alpha + 3 * asCoef < 0)))
    stop("``alpha'' and ``asCoef'' are not valid. See the doc.")

  if ((model == "amix")){
    alpha <- alpha + 3 * asCoef
    asCoef <- - asCoef
  }

  alpha <- as.double(alpha)

  switch(model,
         log = pbvlog(q = q, alpha = alpha, mar1 = mar1,
           mar2 = mar2, lower.tail = lower.tail),
         alog = pbvalog(q = q, alpha = alpha, asCoef1 = asCoef1,
           asCoef2 = asCoef2, mar1 = mar1, mar2 = mar2,
           lower.tail = lower.tail),
         nlog = pbvnlog(q = q, alpha = alpha, mar1 = mar1,
           mar2 = mar2, lower.tail = lower.tail),
         anlog = pbvanlog(q = q, alpha = alpha, asCoef1 = asCoef1,
           asCoef2 = asCoef2, mar1 = mar1, mar2 = mar2,
           lower.tail = lower.tail),
         mix = pbvmix(q = q, alpha = alpha, mar1 = mar1,
           mar2 = mar2, lower.tail = lower.tail),
         amix = pbvamix(q = q, alpha = alpha, asCoef = asCoef,
           mar1 = mar1, mar2 = mar2, lower.tail = lower.tail)
         )
}

pbvlog <- function(q, alpha, mar1 = c(0,1,0), mar2 = mar1,
                   lower.tail = TRUE){
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0 ||
     alpha > 1) stop("invalid argument for `alpha'")

  if(is.null(dim(q))) dim(q) <- c(1,2)

  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])
  
  v <- apply(q^(-1/alpha), 1, sum)^alpha
  pp <- exp(-v)
  
  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}

pbvalog <- function(q, alpha, asCoef1, asCoef2, mar1 = c(0,1,0),
                    mar2 = mar1, lower.tail = TRUE){

  asy <- c(asCoef1, asCoef2)
  
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0 ||
     alpha > 1) stop("invalid argument for `alpha'")
  
  if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
     max(asy) > 1) stop("invalid argument for `asy'")
  
  if(alpha == 1 || any(asy == 0)) {
    asy <- c(0,0)
    alpha <- 1
  }
  
  if(is.null(dim(q))) dim(q) <- c(1,2)

  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])

  asy <- rep(asy,rep(nrow(q),2))
  v <- apply((q/asy)^(-1/alpha), 1, sum)^alpha +
    apply((1-asy) / q, 1, sum)
  pp <- exp(-v)
  
  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}


pbvnlog <- function(q, alpha, mar1 = c(0,1,0), mar2 = mar1,
                    lower.tail = TRUE){
  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0)
    stop("invalid argument for `alpha'")
  
  if(is.null(dim(q))) dim(q) <- c(1,2)
  
  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])
  
  v <- apply(1/q, 1, sum) - apply(q^alpha, 1, sum)^(-1/alpha)
  pp <- exp(-v)
  
  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}

pbvanlog <- function(q, alpha, asCoef1, asCoef2, mar1 = c(0,1,0),
                    mar2 = mar1, lower.tail = TRUE){

  asy <- c(asCoef1, asCoef2)

  if(length(alpha) != 1 || mode(alpha) != "numeric" || alpha <= 0)
    stop("invalid argument for `alpha'")
  
  if(length(asy) != 2 || mode(asy) != "numeric" || min(asy) < 0 ||
     max(asy) > 1)
    stop("invalid argument for `asy'")
  
  if(is.null(dim(q))) dim(q) <- c(1,2)

  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])

  asy <- rep(asy,rep(nrow(q),2))
  v <- apply( 1/q, 1, sum) - apply((q/asy)^alpha, 1, sum)^(-1/alpha)
  pp <- exp(-v)

  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}

pbvmix <- function(q, alpha, mar1 = c(0,1,0), mar2 = mar1,
                   lower.tail = TRUE){

  if(length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")

  if((alpha < 0) || (alpha > 1))
    stop("`alpha' must be in [0,1].")

  if(is.null(dim(q))) dim(q) <- c(1,2)

  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])

  v <- apply( 1/q, 1, sum) - alpha / apply(q, 1, sum)
  pp <- exp(-v)

  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}

pbvamix <- function(q, alpha, asCoef, mar1 = c(0,1,0),
                    mar2 = mar1, lower.tail = TRUE){

  if(length(alpha) != 1 || mode(alpha) != "numeric")
    stop("invalid argument for `alpha'")

  if(length(asCoef) != 1 || mode(asCoef) != "numeric")
    stop("invalid argument for `asCoef'")
  
  if(alpha < 0)
    stop("`alpha' must be non-negative.")

  if((alpha + asCoef) > 1)
    stop("`alpha' + `asCoef' cannot be greater than one")
  
  if((alpha + 2*asCoef) > 1)
    stop("`alpha' + `2*asCoef' cannot be greater than one")
  
  if((alpha + 3*asCoef) < 0)
    stop("`alpha' + `3*asCoef' must be non-negative")
  
  if(is.null(dim(q))) dim(q) <- c(1,2)

  q[,1] <- gpd2frech(q[,1], mar1[1], mar1[2], mar1[3])
  q[,2] <- gpd2frech(q[,2], mar2[1], mar2[2], mar2[3])

  coef <- c(alpha + asCoef, alpha + 2 * asCoef)
  dim(coef) <- c(1,2)
  
  v <- apply( 1/q, 1, sum) - apply(coef*q, 1, sum) /
    apply(q, 1, sum)^2
  pp <- exp(-v)

  if(!lower.tail)
    pp <- 1 - pgpd(1/q[,1], lower.tail = FALSE) -
      pgpd(1/q[,2], lower.tail = FALSE) + pp
  
  return(pp)
}

  
  
