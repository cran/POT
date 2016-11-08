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



tsdep.plot <- function(data, u, ..., xlab, ylab, n.boot = 100, show.lines = TRUE,
                       lag.max, ci = 0.95, block.size = 5 * lag.max, angle = 90,
                       arrow.length = 0.1){

  if (missing(xlab))
    xlab <- "Lag"

  if (missing(ylab))
    ylab <- expression(Lambda[tau])
  
  ##First we have to transform the time series to unit frechet margins
  n <- length(data)

  k <- fitexi(data, u)$tim.cond

  if (k == 1)
    stop("The time independence condition is equal to 1...\n Please check this value as exceedances may no (all) be independent")

  if (missing(lag.max))
    lag.max <- k
    
  ##First we have to transform data to unit frechet margins
  ##using empirical probabilities
  data <- rank(data) / (n + 1)
  data <- -1 / log(data)

  gamma1 <- gamma2 <- rep(NA, lag.max)
  gamma1.ci <- gamma2.ci <- matrix(NA, ncol = 2, nrow = lag.max)
  gamma1.boot <- matrix(NA, nrow = n.boot, ncol = lag.max)
  
  for (lag in 1:lag.max){
    ##Then define a new time series such as T_{i,lag} = min(X_i, X_{i+lag})
    T <- apply(cbind(data[-((n-lag+1):n)], data[-(1:lag)]), 1, min)

    Tsort <- sort(T, decreasing = TRUE)
       
    gamma1[lag] <- 2 * (mean(log(Tsort[1:(k-1)])) -
                        log(Tsort[k])) - 1

    ##Build confidence intervals by bootstrap
    n.block <- floor((n - lag) / block.size)
    
    for (i in 1:n.boot){
      idxs <- sample(1:n.block, replace = TRUE)
      Tboot <- matrix(rep(idxs-1, block.size), nrow = block.size,
                      byrow = TRUE)
      Tboot <- Tboot * block.size + 1:block.size
      Tboot <- T[as.vector(Tboot)]
      Tboot <- sort(Tboot, decreasing = TRUE)
      gamma1.boot[i, lag] <- 2 * (mean(log(Tboot[1:(k-1)])) -
                                  log(Tboot[k])) - 1
      
    }

  }

  gamma1.boot <- pmin(gamma1.boot, 1)
  gamma1.boot <- pmax(gamma1.boot, -1)
  
  gamma1.ci <- apply(gamma1.boot, 2, quantile,
                     p = c((1-ci)/2, 1 - (1 - ci)/2))

  plot(1:lag.max, gamma1, xlab = xlab, ylab = ylab,
       ylim = c(-1,1), ...)
 
  if (show.lines)
    abline(h = c(1, 0), col = "blue")

  arrows(1:lag.max, gamma1.ci[1,], 1:lag.max, gamma1.ci[2,], col = "grey",
         angle = angle, code = 3, length = arrow.length)

  
  invisible(list(gamma1 = gamma1, gamma1.ci = gamma1.ci))
}
