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



##The return level plot for object of class ``mcpot''
retlev.mcpot <- function(fitted, opy, exi, main, xlab,
                         ylab, xlimsup, ...){
  if(!inherits(fitted, "mcpot"))
    stop("Use only with 'mcpot' objects")
  
  loc <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  data <- fitted$data
  pat <- fitted$pat

  if (missing(exi) || !is.numeric(exi) || !is.finite(exi))
    exi <- fitexi(data, loc)$exi
  
  pot.fun <- function(T){
    p <- 1 - 1 / T
    q <- (1 - p^(1/opy/exi)) / pat
    q <- loc - scale / shape * (1 - q^(-shape))
    return(q)
  }
  
  if (missing(opy)){
    warning("Argument ``opy'' is missing. Setting it to 365.25.")
    opy <- 365.25
  }
  
  if (missing(main)) main <- 'Return Level Plot'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- 500

  eps <- 10^(-3)
  plot(pot.fun, from = 1 + eps, to = xlimsup, log = "x",
       xlab = xlab, ylab = ylab, main = main, ...)

  invisible(pot.fun)
}
  
