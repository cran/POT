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



dexi <- function(x, n.sim = 1000, n.mc = length(x$data),
                 plot = TRUE, ...){
  
  thresh <- x$threshold
  scale <- x$param["scale"]
  shape <- x$param["shape"]
  alpha <- x$param["alpha"]
  pat <- x$pat
  model <- x$model
  
  scale.new <- shape * thresh / (pat^(-shape) - 1)

  if (model %in% c("log", "nlog"))
    param <- list(alpha = alpha)

  if (model %in% c("alog", "anlog"))
    param <- list(alpha = alpha, asCoef1 = x$param["asCoef1"],
                  asCoef2 = x$param["asCoef2"])

  if (model == "mix")
    param <- list(alpha = alpha, asCoef = 0)
  
  if (model == "amix")
    param <- list(alpha = alpha, asCoef = x$param["asCoef"])

  param <- c(param, list(n = n.mc, model = model))

  exi <- rep(0, n.sim)
  mc <- rep(0, n.mc)
  
  for (i in 1:n.sim){
    mc <- do.call("simmc", param) #see mcpot-simmc.R
    mc <- qgpd(mc, 0, scale.new, shape)

    while(sum(mc > thresh) < 2){
      mc <- do.call("simmc", param) #see mcpot-simmc.R
      mc <- qgpd(mc, 0, scale.new, shape)
    }
    
    exi[i] <- fitexi(mc, thresh)$exi
  }

  if (plot)
    plot(density(exi, bw = sd(exi) /2), ...)
  
  invisible(exi)
}
