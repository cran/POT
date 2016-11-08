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




convassess.mcpot <- function(fitted, n = 50){

  if (!inherits(fitted, "mcpot"))
    stop("Use only with 'mcpot' objects")
  
  #if (!("mcpot" %in% class(fitted)))
  #  stop("``fitted'' must be of class ``mcpot''")

  if (fitted$model != "log")
    return(cat("This function is only implemented for model ``log''.\nStill work to do...\n"))

  nat <- fitted$nat
  thresh <- fitted$threshold
  param <- fitted$param
  data <- fitted$data
  nobs <- length(data)
  model <- fitted$model
  
  fun <- function()
    fitmcgpd(data, thresh, start = start, model = model,
             std.err.type = "none")

  est <- startValues <- matrix(NA, ncol = 3, nrow = n)
  colnames(est) <- colnames(startValues) <- c("scale", "shape",
                                              "alpha")

  optValues <- rep(NA, n)
  
  for (i in 1:n){

    idx <- sample(1:nobs, size = nat, replace = TRUE)
    x <- data[idx]
    startValues[i,1:2] <- fitgpd(x, thresh[1], "pwmu",
                                 hybrid = TRUE)$param
    startValues[i,3] <- runif(1, 0.2, 0.8)
    start <- list(scale = startValues[i,"scale"],
                  shape = startValues[i,"shape"],
                  alpha = startValues[i,"alpha"])
    
    fit <- fun()
    optValues[i] <- fit$opt.value
    
    est[i,] <- fit$param
  }

  idx <- which(optValues == -1e6)
  optValues[idx] <- NA
  
  par(mfrow=c(3,3))
  
  plot(startValues[,c("scale","shape")], xlab = "Scale",
       ylab = "Shape", main = "Starting Values",
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,c("scale","shape")])
    points(startValues[idx,c("scale","shape")], col = "red",
           pch = 15)
  }
  
  else
    points(startValues[,c("scale","shape")])

  ##Starting Values For Alpha
  plot(startValues[,"alpha"], xlab = "Index",
       ylab = expression(alpha),
       main = expression(paste("Starting Values for ", alpha)),
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,"alpha"])
    points(startValues[idx,"alpha"], col = "red", pch = 15)
  }
  
  else
    points(startValues[,"alpha"])

   ##Scale Marge 1 Trace Plot
  plot(1:n, est[,"scale"], xlab = "Index", ylab = "Scale", type = "n",
       main = "Scale Trace Plot")
  abline(h=fitted$param["scale"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"scale"])
    points(idx, est[idx,"scale"], col = "red", pch = 15)
  }

  else
    points(est[,"scale"])

  ##Shape Marge 1 Trace Plot
  plot(1:n, est[,"shape"], xlab = "Index", ylab = "Shape", type = "n",
       main = "Shape Trace Plot")
  abline(h=fitted$param["shape"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"shape"])
    points(idx, est[idx,"shape"], col = "red", pch = 15)
  }

  else
    points(est[,"shape"])

  ##Alpha Trace Plot
  plot(1:n, est[,"alpha"], xlab = "Index", ylab = expression(alpha), type = "n",
       main = expression(paste(alpha, " Trace Plot")))
  abline(h=fitted$param["alpha"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"alpha"])
    points(idx, est[idx,"alpha"], col = "red", pch = 15)
  }

  else
    points(est[,"alpha"])

  ##Deviance Trace Plot
  plot(1:n, optValues, xlab = "Index", ylab = "Deviance", type = "n",
       main = "Deviance Trace Plot")
  abline(h=fitted$deviance, col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], optValues[-idx])
    points(idx, rep(fitted$deviance, length(idx)), col = "red",
           pch = 15)
  }

  else
    points(optValues)

  
  invisible(list(start.values = startValues, opt.values = optValues,
                 param = est))
    
}
