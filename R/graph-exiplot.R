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



exiplot <- function(data, u.range, tim.cond = 1, n.u = 50,
                    xlab, ylab, ...){

  if (all(data <= u.range[1]))
    stop("No data above u.range[1] !!! Please specify suitable value.")

  if (missing(xlab))
    xlab <- "Threshold"

  if (missing(ylab))
    ylab <- "Extremal Index"

  if (missing(u.range))
    u.range <- c(min(data[,"obs"]), quantile(data[,"obs"], 0.9))
  
  u.range <- seq(u.range[1], u.range[2], length.out = n.u)

  exi <- rep(NA, n.u)
  
  for (i in 1:n.u)
    exi[i] <- attributes(clust(data, u.range[i], tim.cond = tim.cond))$exi

  plot(u.range, exi, xlab = xlab, ylab = ylab, ...)

  invisible(cbind(thresh = u.range, exi = exi))

}
    
