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

## This file contains several functions to plot Peaks Over a Threshold.



dens.uvpot <- function(object, main, xlab, ylab,
                       dens.adj = 1, kern.lty = 2,
                       rug = TRUE, plot.kernel = TRUE,
                       plot.hist = TRUE, hist.col = NULL,
                       ...){

  if(!inherits(object, "uvpot"))
    stop("Use only with 'uvpot' objects")
  if (object$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- object$exceed
  loc <- object$threshold[1]

  if (length(unique(loc)) != 1)
      stop("Density plot not avalaible for varying threshold...")

  scale <- object$param["scale"]
  shape <- object$param["shape"]
  n <- object$nat

  dens <- function(x) dgpd(x, loc, scale, shape)
  eps <- 10^(-5)

  if ( missing(main) ) main <- 'Density Plot'
  if ( missing(xlab) ) xlab <- 'Quantile'
  if ( missing(ylab) ) ylab <- 'Density'
  
  plot(dens, from = loc + eps, to = 1.25 * max(data), main = main,
       xlab = xlab, ylab = ylab, ..., type = "n")

  if (plot.hist)
    hist(data, add = TRUE, freq = FALSE, col = hist.col)

  if (plot.kernel){
    ##A non parametric estimate of the density from Alec Stephenson's code
    flipexceed <- c(data, 2*loc - data)
    flip.density <- density(flipexceed, adj=dens.adj, from = loc + eps,
                            to = 1.25 * max(data))
    flip.density$y <- 2 * flip.density$y
    lines(flip.density, lty = kern.lty)
  }

  if (rug) rug(data)

  plot(dens, from = loc + eps, to = 1.25 * max(data),
       add = TRUE)

}

