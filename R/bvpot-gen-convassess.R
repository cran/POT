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




convassess.bvpot <- function(fitted, n = 50){
  if (!inherits(fitted, "bvpot"))
    stop("Use only with 'bvpot' objects")

  #if (!("bvpot" %in% class(fitted)))
  #  stop("``fitted'' must be of class ``bvpot''")

  if (fitted$model != "log")
    return(cat("This function is only implemented for model ``log''.\nStill work to do...\n"))

  nat <- fitted$nat
  thresh <- fitted$threshold
  param <- fitted$param
  data <- fitted$data
  nobs <- nrow(data)
  model <- fitted$model
  
  fun <- function()
    fitbvgpd(data, thresh, start = start, model = model,
             std.err.type = "none")

  est <- startValues <- matrix(NA, ncol = 5, nrow = n)
  colnames(est) <- colnames(startValues) <- c("scale1", "shape1",
                                              "scale2", "shape2",
                                              "alpha")

  optValues <- rep(NA, n)
  
  for (i in 1:n){

    idx <- sample(1:nobs, size = nat, replace = TRUE)
    x <- data[idx,]
    startValues[i,1:2] <- fitgpd(x[,1], thresh[1], "pwmu",
                                 hybrid = TRUE)$param
    startValues[i,3:4] <- fitgpd(x[,2], thresh[1], "pwmu",
                                 hybrid = TRUE)$param
    startValues[i,5] <- runif(1, 0.2, 0.8)
    start <- list(scale1 = startValues[i,"scale1"],
                  shape1 = startValues[i,"shape1"],
                  scale2 = startValues[i,"scale2"],
                  shape2 = startValues[i,"shape2"],
                  alpha = startValues[i,"alpha"])
    
    fit <- fun()
    optValues[i] <- fit$opt.value
    
    est[i,] <- fit$param
  }

  idx <- which(optValues == 2e6)
  optValues[idx] <- NA
  
  par(mfrow=c(3,3))
  
  ##Starting Values Marge 1
  plot(startValues[,c("scale1","shape1")], xlab = "Scale Marge 1",
       ylab = "Shape Marge 1", main = "Starting Values Marge 1",
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,c("scale1","shape1")])
    points(startValues[idx,c("scale1","shape1")], col = "red",
           pch = 15)
  }
  
  else
    points(startValues[,c("scale1","shape1")])

  ##Starting Values Marge 2
  plot(startValues[,c("scale2","shape2")], xlab = "Scale Marge 2",
       ylab = "Shape Marge 2", main = "Starting Values Marge 2",
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,c("scale2","shape2")])
    points(startValues[idx,c("scale2","shape2")], col = "red",
           pch = 15)
  }
  
  else
    points(startValues[,c("scale2","shape2")])

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
  plot(1:n, est[,"scale1"], xlab = "Index", ylab = "Scale Marge 1", type = "n",
       main = "Scale Marge 1 Trace Plot")
  abline(h=fitted$param["scale1"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"scale1"])
    points(idx, est[idx,"scale1"], col = "red", pch = 15)
  }

  else
    points(est[,"scale1"])

  ##Shape Marge 1 Trace Plot
  plot(1:n, est[,"shape1"], xlab = "Index", ylab = "Shape Marge 1", type = "n",
       main = "Shape Marge 1 Trace Plot")
  abline(h=fitted$param["shape1"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"shape1"])
    points(idx, est[idx,"shape1"], col = "red", pch = 15)
  }

  else
    points(est[,"shape1"])

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

  ##Scale Marge 2 Trace Plot
  plot(1:n, est[,"scale2"], xlab = "Index", ylab = "Scale Marge 2", type = "n",
       main = "Scale Marge 2 Trace Plot")
  abline(h=fitted$param["scale2"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"scale2"])
    points(idx, est[idx,"scale2"], col = "red", pch = 15)
  }

  else
    points(est[,"scale2"])

  ##Shape Marge 2 Trace Plot
  plot(1:n, est[,"shape2"], xlab = "Index", ylab = "Shape Marge 2", type = "n",
       main = "Shape Marge 2 Trace Plot")
  abline(h=fitted$param["shape2"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"shape2"])
    points(idx, est[idx,"shape2"], col = "red", pch = 15)
  }

  else
    points(est[,"shape2"])

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

