## This file contains several function to select the threshold
## for which the asymptotical approximation of Peaks Over a
## Threshold by a GP distribution is quite good.

diplot <- function(data, u.range, main, xlab, ylab,
                   conf=0.95, ...){

  if ( !any(colnames(data) == "obs") )
    stop("``data'' should have a column named ``obs''...")

  if ( !any(colnames(data) == "time") )
    stop("``data'' should have a column named ``time''...")

  data <- na.omit(data)
  date <- data[,"time"]
  samp <- data[,"obs"]
                                       
  if (length(samp)<5){
    stop('Not enough data for a Dispersion Index Plot')
  }

  M <- diff(range(date))

  if (missing(u.range)) u.range <- c(min(samp),max(samp[-(1:4)]))

  thresh <- seq(u.range[1],u.range[2],
                length=max(200,length(samp)))

  DI <- NULL

  date <- floor(date)
  tim.rec <- range(date)    
 
  for (u in thresh){

    nb.occ <- NULL
    idx.excess <- samp > u
    lambda <- sum(idx.excess) / M
   
    for (year in tim.rec[1]:tim.rec[2])
      nb.occ <- c(nb.occ, sum(idx.excess &
                              (date == year)))
      
    DI <- c(DI, var(nb.occ)/lambda)

  }

  conf_sup <- qchisq(1-(1-conf)/2,M-1)/(M-1)
  conf_inf <- qchisq((1-conf)/2,M-1)/(M-1)

  if ( missing(main) ) main <- 'Dispersion Index Plot'
  if ( missing(xlab) ) xlab <- 'Threshold'
  if ( missing(ylab) ) ylab <- 'Dispersion Index'
  
  plot(c(thresh,thresh[1]), c(DI, conf_sup), xlab=xlab, ylab=ylab,
       type='n', main = main, ...)
  rect(0, conf_inf, 2*u.range[2], conf_sup, col= 'lightgrey', border = FALSE)
  lines(thresh, DI)
  return(invisible(list(thresh=thresh,DI=DI)))

}

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


lmomplot <- function(data, u.range, identify = TRUE, ...){

  data <- sort(as.numeric(data))

  n <- length(data)
  
  if ( n < 5){
    stop('Not enougth data for a L-moments Plot.\nLower limit
 for the threshold is too high')
  }
  
  if (missing(u.range)) {
    
    u.range <- c(data[1], data[n - 4])
    u.range <- u.range - .Machine$double.eps^0.5
    
  }

  data <- data[ data > u.range[1] ]
  
  if (all(data <= u.range[2])) 
    stop("upper limit for threshold is too high")
  
  thresh <- seq(u.range[1],u.range[2],
                length=max(50,length(data)))

  point <- NULL
  for ( u in thresh){

    lmoments34 <- samlmu(data[data > u ], sort.data = FALSE)[3:4]
    point <-  cbind(point,lmoments34)

  }

  courbe_theo <- function(Tau3){
    Tau3*(1+5*Tau3)/(5+Tau3)
  }

  plot(courbe_theo, main='L-Moments Plot', xlab=expression(tau[3]),
       ylab=expression(tau[4]), col='grey',lwd=2,...)
  lines(point[1,],point[2,],type='b')
  

  if (identify)  identify(point[1,],point[2,],
                          labels=round(thresh,digits=2),offset=1)
  
}
  
