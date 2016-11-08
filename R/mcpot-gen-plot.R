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

##The generic function for graphical diagnostic of a Markov chain
##pot object
plot.mcpot <- function(x, opy, exi, mains, which = 1:4,
                       ask = nb.fig < length(which) &&
                       dev.interactive(), acf.type = "partial",
                       ...){
  if(!inherits(x, "mcpot"))
    stop("Use only with 'mcpot' objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")

  if(missing(mains))
    mains <- c("Auto-correlation Plot", "Pickands' Dependence Function",
               "Spectral Density", "Bivariate Return Level Plot")

  else
    if (length(mains) != 4){
      warning("``mains'' must be of length two. Passing to default titles.")
      mains <- c("Auto-correlation Plot", "Pickands' Dependence Function",
                 "Spectral Density", "Bivariate Return Level Plot")
    }


  show <- rep(FALSE, 4)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask){
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  if (show[1])
    acf(x$data, main = mains[1], type = acf.type, ...)
  if (show[2])
    pickdep(x, main = mains[2], ...)
  if (show[3])
    specdens(x, main = mains[3], ...)
  if (show[4])
    retlev(x, opy = opy, exi = exi, main = mains[4], ...)
  
}

