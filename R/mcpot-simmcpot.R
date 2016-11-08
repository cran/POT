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

simmcpot <- function(fitted, plot = TRUE, ...){

  if (all(class(fitted) != "mcpot"))
    stop("``fitted'' must be of class ``mcpot''")

  thresh <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  alpha <- fitted$param["alpha"]
  pat <- fitted$pat
  model <- fitted$model
  n <- length(fitted$data)

  if (model %in% c("log", "nlog"))
    param <- list(alpha = alpha)

  if (model %in% c("alog", "anlog"))
    param <- list(alpha = alpha, asCoef1 = fitted$param["asCoef1"],
                  asCoef2 = fitted$param["asCoef2"])

  if (model == "mix")
    param <- list(alpha = alpha, asCoef = 0)
  
  if (model == "amix")
    param <- list(alpha = alpha, asCoef = fitted$param["asCoef"])

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

