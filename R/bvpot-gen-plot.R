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

##A function to compute spectral densities for each dependance
##models



##The generic function for graphical diagnostic of a bivariate
##pot object
plot.bvpot <- function(x, mains, which = 1:3,
                       ask = nb.fig < length(which) &&
                       dev.interactive(), ...){
  if (!inherits(x, "bvpot"))
    stop("Use only with 'bvpot' objects")

  if (!is.numeric(which) || any(which < 1) || any(which > 3)) 
        stop("`which' must be in 1:3")

  if(missing(mains))
    mains <- c("Pickands' Dependence Function",
               "Bivariate Return Level Plot", "Spectral Density")

  else
    if (length(mains) != 3){
      warning("``mains'' must be of length two. Passing to default titles.")
      mains <- c("Pickands' Dependence Function",
                 "Bivariate Return Level Plot", "Spectral Density")
    }


  show <- rep(FALSE, 3)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask){
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  if (show[1])
    pickdep(x, main = mains[1], ...)

  if (show[2])
    retlev(x, main = mains[2], ...)

  if (show[3])
    specdens(x, main = mains[3], ...)

}
