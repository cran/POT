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


##Bivariate return level plot
retlev.bvpot <- function(fitted, p = seq(0.75,0.95,0.05), main,
                         n = 5000, only.excess = FALSE, ...){
  if (!inherits(fitted, "bvpot"))
    stop("Use only with 'bvpot' objects")
  #if (all(class(fitted) != "bvpot"))
  #  stop("``fitted'' should be an object of class ``bvpot''.")

  if (missing(main))
    main <- "Bivariate Return Level Plot"

  w <- c(0, seq(0,1, length.out = n), 1)
  ##The Pickands' dependence function
  A <- pickdep(fitted, plot = FALSE)(w)
  y1 <- y2 <- matrix(NA, ncol = length(p), nrow = n + 2)
  colnames(y1) <- colnames(y2) <- p
  rownames(y1) <- rownames(y2) <- 1:(n+2)

  for (i in 1:length(p)){
    z1 <- - A / (w * log(p[i]))
    z2 <- w / (1 - w) * z1
    
    ##We have to transform frechet observation to original
    ##scale
    y1[,i] <- ((1-exp(-1/z1)) / fitted$pat[1])^
    (-fitted$param["shape1"]) - 1
    y1[,i] <- fitted$threshold[1] + fitted$param["scale1"] /
      fitted$param["shape1"] * y1[,i]
    y2[,i] <- ((1-exp(-1/z2)) / fitted$pat[2])^
    (-fitted$param["shape2"]) - 1
    y2[,i] <- fitted$threshold[2] + fitted$param["scale2"] /
      fitted$param["shape2"] * y2[,i]
  }

  data1 <- fitted$data[,1]
  data2 <- fitted$data[,2]

  if (only.excess){
    idx <- which(data1 > fitted$threshold[1] |
                 data2 > fitted$threshold[2])
    data1 <- data1[idx]
    data2 <- data2[idx]
  }
  
  plot(data1, data2, main = main, ...)
  lines(y1, y2)
  
  ##Add the marginal threshold axis
  abline(v = fitted$threshold[1])
  abline(h = fitted$threshold[2])

  invisible(list(y1 = y1, y2 = y2))

}

