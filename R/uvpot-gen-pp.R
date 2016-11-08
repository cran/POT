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

## This file contains several functions to plot Peaks Over a Threshold.


pp.uvpot <- function(fitted, main, xlab,
                     ylab, ci = TRUE,...){

  if(!inherits(fitted, "uvpot"))
    stop("Use only with 'uvpot' objects")
  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceed
  loc <- fitted$threshold[1]
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  n <- fitted$nat

  p_emp <- ppoints(n)
  p_fit <- pgpd(sort(data), loc, scale, shape)
  
  if ( missing(main) ) main <- 'Probability plot'
  if ( missing(ylab) ) ylab <- 'Model'
  if ( missing(xlab) ) xlab <- 'Empirical'
  if(length(p_fit) != length(p_emp))
    stop("internal error in pp.uvpot()")
  plot(p_emp, p_fit, main = main, xlab = xlab, ylab = ylab, ...)
  abline(0,1)

  if (ci){
    p_emp <- 1:n / (n+1)
    samp <- rgpd(1000*n, loc, scale, shape)
    samp <- matrix(samp, n, 1000)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- pgpd(samp[25,], loc, scale, shape)
    ci_sup <- pgpd(samp[975,], loc, scale, shape)
    points( p_emp, ci_inf, pch = '-')
    points( p_emp, ci_sup, pch = '-')
  }
}
    

