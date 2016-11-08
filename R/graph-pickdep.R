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



##A function to plot the Pickands dependance function
pickdep <- function(fitted, main, bound = TRUE, plot = TRUE,
                    ...){

  if(!is.list(fitted))
    if (!inherits(fitted, "bvpot"))
      stop("Use only with 'bvpot' objects")
  else if(!all(c("model", "param") %in% names(fitted))) #function is used in fitbvgpd()
    stop("Use only with 'bvpot' objects")
  
  model <- fitted$model
  alpha <- fitted$param["alpha"]

  if (model == "log"){
    ##Logistic case :
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- ((1-w)^(1/alpha) + w^(1/alpha))^alpha
      }

      else
        ans <- ((1-w)^(1/alpha) + w^(1/alpha))^alpha
      
      return(ans)
    }
  }
  
  if (model == "nlog"){
    ##Negative logistic case:
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- 1 - ((1-w)^(-alpha) + w^(-alpha))^(-1/alpha)
      }

      else
        ans <- 1 - ((1-w)^(-alpha) + w^(-alpha))^(-1/alpha)
      
      return(ans)
    }
  }

  if (model == "alog"){
    ##Asymetric logistic case:
    asCoef1 <- fitted$param["asCoef1"]
    asCoef2 <- fitted$param["asCoef2"]
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- (1 - asCoef1)*(1-w) + (1 - asCoef2) * w +
          ( (asCoef1 * (1-w))^(1/alpha) + (asCoef2 * w)^(1/alpha) )^alpha
      }

      else
        ans <- (1 - asCoef1)*(1-w) + (1 - asCoef2) * w +
          ( (asCoef1 * (1-w))^(1/alpha) + (asCoef2 * w)^(1/alpha) )^alpha
    
      return(ans)
    }
  }

  if (model == "anlog"){
    ##Asymetric negatif logistic case:
    asCoef1 <- fitted$param["asCoef1"]
    asCoef2 <- fitted$param["asCoef2"]
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- 1 - ( ((1-w)*asCoef1)^(-alpha) +
                          (w*asCoef2)^(-alpha) )^(-1/alpha)
      }

      else
        ans <- 1 - ( ((1-w)*asCoef1)^(-alpha) +
                    (w*asCoef2)^(-alpha) )^(-1/alpha)
      
      return(ans)
    }
  }

  if (model == "mix"){
    ##Mixed model:
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- 1 - w * (1-w) * alpha
      }

      else
        ans <- 1 - w * (1-w) * alpha
      
      return(ans)
    }
  }

   if (model == "amix"){
    ##Asymetric Mixed model:
     asCoef <- fitted$param["asCoef"]
    A <- function(w){
      ans <- rep(NA, length(w))
      idx <- which((w <= 0) | (w > 1))

      if (length(idx) > 0){
        w <- w[-idx]
        ans[-idx] <- 1 - (alpha + 2 * asCoef) * w +
          (alpha + 3 * asCoef)* w^2 - asCoef * w^3
      }

      else
        ans <-  1 - (alpha + 2 * asCoef) * w +
          (alpha + 3 * asCoef)* w^2 - asCoef * w^3
      
      return(ans)
    }
  }

  if (plot){
    if (missing(main))
        main <- "Pickands' Dependence Function"
    
    plot(A, ylim = c(0.5, 1), xlim = c(0,1), main = main,
         type = "n")
    
    if (bound){
      lines(x= c(0,1), y = c(1,1), col = "grey", ...)
      lines(x = c(0,0.5,1), y = c(1, 0.5, 1), col = "grey", ...)
    }
    plot(A, add = TRUE)
  }
  
  attributes(A) <- list(model = model)
  invisible(A)
}

