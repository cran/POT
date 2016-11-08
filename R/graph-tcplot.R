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



tcplot <- function (data, u.range, cmax = FALSE, r = 1, 
    ulow = -Inf, rlow = 1, nt = 25, which = 1:npar, conf = 0.95, 
    lty = 1, lwd = 1, type = "b", cilty = 1, ask = nb.fig < length(which) && 
        dev.interactive(), ...){

  n <- length(data)
  data <- sort(data)
  
  if (missing(u.range)) {
    
    u.range <- c(data[1], data[n - 4])
    u.range <- u.range - .Machine$double.eps^0.5
    
  }
  
  u <- seq(u.range[1], u.range[2], length = nt)
  locs <- scls <- shps <- matrix(NA, nrow = nt, ncol = 3)
  dimnames(locs) <- list(round(u, 2), c("lower", "loc", "upper"))
  dimnames(shps) <- list(round(u, 2), c("lower", "shape", "upper"))
  
  pname <- "mscale"
  npar <- 2
  
  dimnames(scls) <- list(round(u, 2), c("lower", pname, "upper"))
  z <- gpdmle(data, u[1], corr = TRUE, ...)
  stvals <- as.list(round(fitted(z), 3))
  
  for (i in 1:nt) {
    
    z <- gpdmle(data, u[i], corr = TRUE, ...)
    stvals <- as.list(fitted(z))
    mles <- fitted(z)
    stderrs <- z$std.err
    cnst <- qnorm((1 + conf)/2)
    shp <- mles["shape"]
    scl <- mles["scale"]
    shpse <- stderrs["shape"]
    sclse <- stderrs["scale"]
    
    scl <- scl - shp * u[i]
    covar <- z$corr[1, 2] * prod(stderrs)
    sclse <- sqrt(sclse^2 - 2 * u[i] * covar + (u[i] * 
                                                shpse)^2)
    
    scls[i, ] <- c(scl - cnst * sclse, scl, scl + cnst * 
                   sclse)
    shps[i, ] <- c(shp - cnst * shpse, shp, shp + cnst * 
                   shpse)
    
  }
  
  show <- rep(FALSE, npar)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask) {
    
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  
  if (show[1]) {
    
    matplot(u, scls, type = "n", xlab = "Threshold", 
            ylab = "Modified Scale")
    lines(u, scls[, 2], lty = lty, lwd = lwd, type = type)
    segments(u, scls[, 1], u, scls[, 3], lty = cilty)
    
  }
  if (show[2]) {
    
    matplot(u, shps, type = "n", xlab = "Threshold", 
            ylab = "Shape")
    lines(u, shps[, 2], lty = lty, lwd = lwd, type = type)
      segments(u, shps[, 1], u, shps[, 3], lty = cilty)
    
  }
  
  rtlist <- list(scales = scls, shapes = shps)
  invisible(rtlist)
  
}



