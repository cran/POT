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

chimeas <- function(data, u.range, n.u = 500, xlab, ylabs, ci = 0.95, boot = FALSE,
                    n.boot = 250, block.size = 50, show.bound = TRUE, which = 1:2,
                    ask = nb.fig < length(which) && dev.interactive(), ...,
                    col.ci = "grey", col.bound = "blue", lty.ci = 1, lty.bound = 1){

  if (ncol(data) != 2)
    stop("data must be a maxtrix with two columns")

  if (missing(ylabs))
    ylabs <- c(expression(chi), expression(bar(chi)))

  if (missing(xlab))
    xlab <- "u"

  show <- rep(FALSE, 2)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))

  if (ask){
    op <- par(ask=TRUE)
    on.exit(par(op))
  }

  n <- nrow(data)
  ##First we have to transform data to uniform margins
  ##using empirical probabilities
  data <- apply(data, 2, rank) / (n + 1)

  eps <- sqrt(.Machine$double.eps)

  #Marginal maximum i.e. max(X_i, Y_i)
  M.max <- apply(data, 1, max)
  #Marginal minimum i.e. min(X_i, Y_i)
  M.min <- apply(data, 1, min)
  
  if (missing(u.range))
    u.range <- c(min(M.max) + eps, max(M.min) - eps)
  

  u <- seq(u.range[1], u.range[2], length = n.u)
  probJu <- probJuBar <- rep(NA, n.u)
    
  for (i in 1:n.u){
    ##P[X <= u and Y <= u] empirical estimate
    probJu[i] <- mean(M.max <= u[i])
    ##P[X > u and Y > u] empirical estimate
    probJuBar[i] <- mean(M.min > u[i])
  }
  
  ##Compute Chi(u)
  chi <- 2 - log(probJu) / log(u)

  ##Compute Chibar(u)
  chibar <- 2 * log(1 - u) / log(probJuBar) - 1

  chiBounds <- 2 - log(pmax(2*u - 1, 0)) / log(u)

  if (!boot){
    ##Compute the variance of Chi and Chibar (Delta Method)
    ##Chi and Chibar variance obtained using the TCL
    varChi <- (log(u) * probJu)^(-2) * probJu * (1 - probJu) / n
    varChibar <- (2 * log(1 - u) / (probJuBar * log(probJuBar)^2))^2 *
      probJuBar * (1 - probJuBar) / n

    chi.ci.inf <- pmax(chi - qnorm((1+ci)/2) * sqrt(varChi),
                       chiBounds)
    chi.ci.sup <- pmin(chi + qnorm((1+ci)/2) * sqrt(varChi),
                       1)
    chibar.ci.inf <- pmax(chibar - qnorm((1+ci)/2) * sqrt(varChibar),
                          -1)
    chibar.ci.sup <- pmin(chibar + qnorm((1+ci)/2) * sqrt(varChibar),
                          1)
  }

  else{
    ##Build confidence intervals by boostraping contiguous blocks
    chi.boot <- matrix(NA, nrow = n.boot, ncol = n.u)
    chibar.boot <- matrix(NA, nrow = n.boot, ncol = n.u)
    n.block <- floor(n / block.size)

    for (j in 1:n.boot){
      idxs <- sample(1:n.block, replace = TRUE)
      idxs <- rep(idxs - 1, block.size) * block.size +
        block.size:1
      M.min.boot <- M.min[idxs]
      M.max.boot <- M.max[idxs]

      for (i in 1:n.u){
        ##P[X <= u and Y <= u] empirical estimate
        probJu[i] <- mean(M.max.boot <= u[i])
        ##P[X > u and Y > u] empirical estimate
        probJuBar[i] <- mean(M.min.boot > u[i])
      }

      chi.boot[j,] <- 2 - log(probJu) / log(u)
      chibar.boot[j,] <-  2 * log(1 - u) / log(probJuBar) - 1

    }

    chi.ci <- apply(chi.boot, 2, quantile,
                    p = c((1-ci)/2, 1 - (1 - ci)/2), na.rm = TRUE)
    chibar.ci <- apply(chibar.boot, 2, quantile,
                       p = c((1-ci)/2, 1 - (1 - ci)/2), na.rm = TRUE)

    chi.ci.inf <- pmax(chi.ci[1,], chiBounds)
    chi.ci.sup <- pmin(chi.ci[2,], 1)
    chibar.ci.inf <- pmax(chibar.ci[1,], -1)
    chibar.ci.sup <- pmin(chibar.ci[2,], 1)
  }      

  if (show[1]){
    plot(u, chi, type  ="l", xlab = xlab, ylab = ylabs[1], xlim = c(0,1),
         ylim = c(-1,1), ...)

    if (show.bound){
      lines(u, chiBounds, lty = lty.bound, col = col.bound)
      abline(h = 1, lty = lty.bound, col = col.bound)
    }
    
    lines(u, chi.ci.inf, col = col.ci)
    lines(u, chi.ci.sup, col = col.ci)

  }
  
  if (show[2]){
    plot(u, chibar, type  ="l", xlab = xlab, ylab = ylabs[2], xlim = c(0,1),
         ylim = c(-1,1), ...)

    if (show.bound)
      abline(h = c(1,-1), col = col.bound, lty = lty.bound)
   
    lines(u, chibar.ci.inf, col = col.ci)
    lines(u, chibar.ci.sup, col = col.ci)
  }

  invisible(list(chi = rbind(u = u, chi = chi, chi.ci.inf = chi.ci.inf,
                   chi.ci.sup = chi.ci.sup),
                 chibar = rbind(u = u, chibar = chibar,
                   chibar.ci.inf = chibar.ci.inf,
                   chibar.ci.sup = chibar.ci.sup)))

}

