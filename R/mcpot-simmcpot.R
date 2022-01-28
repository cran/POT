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

simmcpot <- function(object, plot = TRUE, ...){

  if (all(class(object) != "mcpot"))
    stop("``object'' must be of class ``mcpot''")

  thresh <- object$threshold
  scale <- object$param["scale"]
  shape <- object$param["shape"]
  alpha <- object$param["alpha"]
  pat <- object$pat
  model <- object$model
  n <- length(object$data)

  if (model %in% c("log", "nlog"))
    param <- list(alpha = alpha)

  if (model %in% c("alog", "anlog"))
    param <- list(alpha = alpha, asCoef1 = object$param["asCoef1"],
                  asCoef2 = object$param["asCoef2"])

  if (model == "mix")
    param <- list(alpha = alpha, asCoef = 0)
  
  if (model == "amix")
    param <- list(alpha = alpha, asCoef = object$param["asCoef"])

  param <- c(param, list(n = n, model = model))

  prob <- do.call(simmc, param)

  mcpot <- rep(NA, n)
  idx <- which(prob > (1 - pat))

  prob[idx] <- (prob[idx] - (1 - pat)) / pat
  mcpot[idx] <- qgpd(prob[idx], thresh, scale, shape)

  if (plot)
    plot(mcpot, ...)
  
  return(mcpot)
}

