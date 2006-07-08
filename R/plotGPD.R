## This file contains several functions to plot Peaks Over a Threshold.

retlev.gpd <- function(fitted, npy, main, xlab,
                   ylab, xlimsup, ci = TRUE, points = TRUE, ...){
  ## Plot the return level plot of a POT model fitted
  ## Input : ``fitted'' is a POT fitted model, result of function
  ##         ``fitgpd''
  ##         npy is the mean number of events per block -generally
  ##         per year- or equivalently the mean -intensity- of the
  ##         Poisson processus.

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceedances
  loc <- fitted$threshold[1]
  scale <- fitted$estimate[1]
  shape <- fitted$estimate[2]
  
  n <- fitted$nhigh
  
  pot.fun <- function(T){
    p <- rp2prob(T, npy)[,"prob"]
    return(qgpd(p,loc,scale,shape))
  }

  eps <- 10^(-3)

  if (missing(main)) main <- 'Return Level Plot'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- prob2rp((n - .35)/n, npy)[,"retper"] 
  
  plot(pot.fun, from= 1 / npy + eps, to = xlimsup, log='x',
       xlab = xlab, ylab = ylab, main = main, ...)

  if (points){
    p_emp <- (1:n -.35) / n
    points(1 / ( npy * (1 - p_emp) ), sort( data ), pch = 1)
  }

  if (ci){
    p_emp <- (1:n - .35 ) / n
    samp <- rgpd(999*n, loc, scale, shape)
    samp <- matrix(samp, n, 999)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- samp[30,]
    ci_sup <- samp[970,]
    lines( 1 / ( npy * (1 - p_emp) ), ci_inf, lty = 2)
    lines( 1 / ( npy * (1 - p_emp) ), ci_sup, lty = 2)
  }
}


qq.gpd <- function(fitted, main, xlab,
                   ylab, ci = TRUE,...){

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceedances
  loc <- fitted$threshold[1]
  scale <- fitted$estimate[1]
  shape <- fitted$estimate[2]
  n <- fitted$nhigh

  quant_fit <- qgpd(ppoints(n), loc, scale, shape)

  if ( missing(main) ) main <- 'QQ-plot'
  if ( missing(xlab) ) xlab <- 'Model'
  if ( missing(ylab) ) ylab <- 'Empirical'
  
  plot(quant_fit, sort(data), main = main, xlab = xlab,
       ylab = ylab, ...)
  abline(0,1)

  if (ci){
    p_emp <- 1:n / (n+1)
    samp <- rgpd(999*n, loc, scale, shape)
    samp <- matrix(samp, n, 999)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- samp[30,]
    ci_sup <- samp[970,]
    points( quant_fit, ci_inf, pch = '-')
    points( quant_fit, ci_sup, pch = '-')
  }

}

pp.gpd <- function(fitted, main, xlab,
                   ylab, ci = TRUE,...){

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceedances
  loc <- fitted$threshold[1]
  scale <- fitted$estimate[1]
  shape <- fitted$estimate[2]
  n <- fitted$nhigh

  p_emp <- ppoints(n)
  p_fit <- pgpd(sort(data), loc, scale, shape)

  if ( missing(main) ) main <- 'Probability plot'
  if ( missing(ylab) ) ylab <- 'Model'
  if ( missing(xlab) ) xlab <- 'Empirical'
  
  plot(p_emp, p_fit, main = main, xlab = xlab, ylab = ylab, ...)
  abline(0,1)

  if (ci){
    p_emp <- 1:n / (n+1)
    samp <- rgpd(999*n, loc, scale, shape)
    samp <- matrix(samp, n, 999)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- pgpd(samp[30,], loc, scale, shape)
    ci_sup <- pgpd(samp[970,], loc, scale, shape)
    points( p_emp, ci_inf, pch = '-')
    points( p_emp, ci_sup, pch = '-')
  }
}
    

dens.gpd <- function(fitted, main, xlab, ylab,
                     dens.adj = 1, kern.lty = 2,
                     rug = TRUE, ...){

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceedances
  loc <- fitted$threshold[1]

  if (length(unique(loc)) != 1)
      stop("Density plot not avalaible for varying threshold...")

  scale <- fitted$estimate[1]
  shape <- fitted$estimate[2]
  n <- fitted$nhigh

  dens <- function(x) dgpd(x, loc, scale, shape)
  eps <- 10^(-5)

  if ( missing(main) ) main <- 'Density Plot'
  if ( missing(xlab) ) xlab <- 'Quantile'
  if ( missing(ylab) ) ylab <- 'Density'
  
  plot(dens, from = loc + eps, to = 1.25 * max(data), main = main,
       xlab = xlab, ylab = ylab, ...)

  ##A non parametric estimate of the density from Alec Stephenson's code
  flipexceed <- c(data, 2*loc - data)
  flip.density <- density(flipexceed, adj=dens.adj, from = loc + eps,
                to = 1.25 * max(data))
  flip.density$y <- 2 * flip.density$y
  lines(flip.density, lty = kern.lty)

  if (rug) rug(data)
}

plotgpd <- function(fitted, npy, main, which = 1:4,
                      ask = nb.fig < length(which) && 
                     dev.interactive(),ci = TRUE, ...){
  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
  if (any(which == 4) & missing(npy))
    stop("Argument ``npy'' must be present !!!")
  
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1]) {
    pp.gpd(fitted, ci = ci, main, xlim = c(0, 1), ylim = c(0, 
                                                    1), ...)
  }
  if (show[2]) {
    qq.gpd(fitted, ci = ci, main, ...)
  }
  if (show[3]) {
    dens.gpd(fitted, main, ...)
  }
  if (show[4]) {
    retlev.gpd(fitted, npy, main, ...)
  }
}
