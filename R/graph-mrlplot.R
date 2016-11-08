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

## This file contains several function to select the threshold
## for which the asymptotical approximation of Peaks Over a
## Threshold by a GP distribution is quite good.


mrlplot <- function(data, u.range, main, xlab, ylab,
                    nt = max(100, length(data)),
                    lty = rep(1, 3), col = c('grey','black','grey'),
                    conf = 0.95, lwd = c(1, 1.5, 1),...){
  
  data <- sort(data[!is.na(data)])
  n <- length(data)
  
  if (n <= 5) 
    stop("`data' has too few non-missing values")
  
  if (missing(u.range)) {
    
    u.range <- c(data[1], data[n - 4])
    u.range <- u.range - .Machine$double.eps^0.5
    
  }
  
  if (all(data <= u.range[2])) 
    stop("upper limit for threshold is too high")
  
  u <- seq(u.range[1], u.range[2], length = nt)
  x <- matrix(NA, nrow = nt, ncol = 3,
              dimnames = list(NULL,c("lower", "mrl", "upper")))
  
  for (i in 1:nt) {
    
    data <- data[data > u[i]]
    x[i, 2] <- mean(data - u[i])
    sdev <- sqrt(var(data))
    sdev <- (qnorm((1 + conf)/2) * sdev)/sqrt(length(data))
    x[i, 1] <- x[i, 2] - sdev
    x[i, 3] <- x[i, 2] + sdev
    
  }

  if ( missing(main) ) main <- 'Mean Residual Life Plot'
  if ( missing(xlab) ) xlab <- 'Threshold'
  if ( missing(ylab) ) ylab <- 'Mean Excess'
  
  matplot(u, x, type = "l", lty = lty, col = col, main = main, 
          xlab = xlab, ylab = ylab, lwd =lwd, ...)
  invisible(list(x = u, y = x))

}

