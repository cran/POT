#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet            
#   Copyright (c) 2022 Christophe Dutang => replace fitted to object and check sample()
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



convassess.uvpot <- function(object, n = 50)
{

  if (!inherits(object, "uvpot"))
    stop("Use only with 'uvpot' objects")

  if (!(object$est %in% c("MLE","LME","MPLE","MEDIANS","MDPD",
                          "MGF")))
    stop("\n    There's no objective function to optimize!")

  thresh <- object$threshold
  scale <- object$param["scale"]
  shape <- object$param["shape"]
  data <- object$data
  nobs <- length(data)

  if (object$est == "MLE")
    fun <- function()
      fitgpd(data, thresh, "mle", start = start,
             std.err.type = "none")

  if (object$est == "LME"){
    fun <- function()
      fitgpd(data, thresh, "lme", start = start)
    startsLME <- seq(-2, 0, length = n)
  }

  if (object$est == "MEDIANS")
    fun <- function()
      fitgpd(data, thresh, "med", start = start)
      
  if (object$est == "MPLE")
    fun <- function()
      fitgpd(data, thresh, "mple", start = start,
             std.err.type = "none")

  if (object$est == "MDPD")
    fun <- function() 
      fitgpd(data, thresh, "mdpd", start = start)

  if (object$est == "MGF")
    fun <- function()
      fitgpd(data, thresh, "mgf", stat = object$stat,
             start = start)

  est <- startValues <- matrix(NA, ncol = 2, nrow = n)
  colnames(est) <- colnames(startValues) <- c("scale", "shape")
  optValues <- rep(NA, n)
  
  size.boot <- max(nobs, 5)
  
  for (i in 1:n)
  {
    if (object$est == "LME")
    {
      start <- list(x = startsLME[i])
    }else
    {
      x <- sample(data, size = size.boot, replace = TRUE)
      startValues[i,] <- fitgpd(x, thresh, "pwmu", hybrid = TRUE)$param
      start <- list(scale = startValues[i,"scale"],
                    shape = startValues[i,"shape"])
    }

    fit <- fun()
    optValues[i] <- fit$opt.value

    est[i,] <- fit$param
  }

  idx <- which(optValues == 1e6)
  optValues[idx] <- NA
  par(mfrow=c(2,2))

  if(object$est == "LME")
    plot(startsLME, xlab = "Index", ylab = expression(b[0]),
         main = "Starting Values")
  else{
    plot(startValues, xlab = "Scale", ylab = "Shape", main = "Starting Values",
         type = "n")
    
    
    if (length(idx) > 0){
      points(startValues[-idx,])
      points(startValues[idx,], col = "red", pch = 15)
    }
    
    else
      points(startValues)
  }
  
  plot(1:n, optValues, xlab = "Index", ylab = "Objective", type = "n",
       main = "Objective Value Trace Plot")
  abline(h=object$opt.value, col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], optValues[-idx])
    points(idx, rep(object$opt.value, length(idx)), col = "red",
           pch = 15)
  }

  else
    points(optValues)


  plot(1:n, est[,1], xlab = "Index", ylab = "Scale", type = "n",
       main = "Scale Trace Plot")
  abline(h=object$param["scale"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,1])
    points(idx, est[idx,1], col = "red", pch = 15)
  }

  else
    points(est[,1])

  plot(1:n, est[,2], xlab = "Index", ylab = "Shape", type = "n",
       main = "Shape Trace Plot")
  abline(h=object$param["shape"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,2])
    points(idx, est[idx,2], col = "red", pch = 15)
  }

  else
    points(est[,2])


  invisible(list(start.values = startValues, opt.values = optValues,
                 param = est))
    
}

