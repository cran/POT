## This file contains several functions to plot Peaks Over a Threshold.

retlev.uvpot <- function(fitted, npy, main, xlab,
                         ylab, xlimsup, ci = TRUE, points = TRUE,
                         ...){
  ## Plot the return level plot of a POT model fitted
  ## Input : ``fitted'' is a POT fitted model, result of function
  ##         ``fitgpd''
  ##         npy is the mean number of events per block -generally
  ##         per year- or equivalently the mean -intensity- of the
  ##         Poisson processus.

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceed
  loc <- fitted$threshold[1]
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  
  n <- fitted$nat
  
  pot.fun <- function(T){
    p <- rp2prob(T, npy)[,"prob"]
    return(qgpd(p,loc,scale,shape))
  }

  eps <- 10^(-3)

  if (missing(npy)){
    warning("Argument ``npy'' is missing. Setting it to 1.")
    npy <- 1
  }
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
    samp <- rgpd(1000*n, loc, scale, shape)
    samp <- matrix(samp, n, 1000)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- samp[25,]
    ci_sup <- samp[975,]
    lines( 1 / ( npy * (1 - p_emp) ), ci_inf, lty = 2)
    lines( 1 / ( npy * (1 - p_emp) ), ci_sup, lty = 2)
  }

  invisible(pot.fun)
}


qq.uvpot <- function(fitted, main, xlab,
                     ylab, ci = TRUE,...){

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceed
  loc <- fitted$threshold[1]
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  n <- fitted$nat

  quant_fit <- qgpd(ppoints(n), loc, scale, shape)

  if ( missing(main) ) main <- 'QQ-plot'
  if ( missing(xlab) ) xlab <- 'Model'
  if ( missing(ylab) ) ylab <- 'Empirical'
  
  plot(quant_fit, sort(data), main = main, xlab = xlab,
       ylab = ylab, ...)
  abline(0,1)

  if (ci){
    p_emp <- 1:n / (n+1)
    samp <- rgpd(1000*n, loc, scale, shape)
    samp <- matrix(samp, n, 1000)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    ci_inf <- samp[25,]
    ci_sup <- samp[975,]
    points( quant_fit, ci_inf, pch = '-')
    points( quant_fit, ci_sup, pch = '-')
  }

}

pp.uvpot <- function(fitted, main, xlab,
                     ylab, ci = TRUE,...){

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
    

dens.uvpot <- function(fitted, main, xlab, ylab,
                       dens.adj = 1, kern.lty = 2,
                       rug = TRUE, plot.kernel = TRUE,
                       plot.hist = TRUE, hist.col = NULL,
                       ...){

  if (fitted$var.thresh)
    stop("Return Level plot is available only for constant threshold !")
  
  data <- fitted$exceed
  loc <- fitted$threshold[1]

  if (length(unique(loc)) != 1)
      stop("Density plot not avalaible for varying threshold...")

  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  n <- fitted$nat

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

plot.uvpot <- function(x, npy, main, which = 1:4,
                       ask = nb.fig < length(which) && 
                       dev.interactive(),ci = TRUE, ...){
  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
  if (any(which == 4) & missing(npy)){
    warning("Argument ``npy'' should be specified !!! Setting it to 1.")
    npy <- 1
  }
  
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1]) {
    pp.uvpot(x, ci = ci, main, xlim = c(0, 1), ylim = c(0, 1),
             ...)
  }
  if (show[2]) {
    qq.uvpot(x, ci = ci, main, ...)
  }
  if (show[3]) {
    dens.uvpot(x, main, ...)
  }
  if (show[4]) {
    retlev.uvpot(x, npy, main, ...)
  }
}
