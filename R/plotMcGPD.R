##The generic function for graphical diagnostic of a Markov chain
##pot object
plot.mcpot <- function(x, opy, exi, mains, which = 1:4,
                       ask = nb.fig < length(which) &&
                       dev.interactive(), acf.type = "partial",
                       ...){

  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")

  if(missing(mains))
    mains <- c("Auto-correlation Plot", "Pickands' Dependence Function",
               "Spectral Density", "Bivariate Return Level Plot")

  else
    if (length(mains) != 4){
      warning("``mains'' must be of length two. Passing to default titles.")
      mains <- c("Auto-correlation Plot", "Pickands' Dependence Function",
                 "Spectral Density", "Bivariate Return Level Plot")
    }


  show <- rep(FALSE, 4)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask){
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  if (show[1])
    acf(x$data, main = mains[1], type = acf.type, ...)
  if (show[2])
    pickdep(x, main = mains[2], ...)
  if (show[3])
    specdens(x, main = mains[3], ...)
  if (show[4])
    retlev(x, opy = opy, exi = exi, main = mains[4], ...)
  
}

##The return level plot for object of class ``mcpot''
retlev.mcpot <- function(fitted, opy, exi, main, xlab,
                         ylab, xlimsup, ...){
  loc <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  data <- fitted$data
  pat <- fitted$pat

  if (missing(exi) || !is.numeric(exi) || !is.finite(exi))
    exi <- fitexi(data, loc)$exi
  
  pot.fun <- function(T){
    p <- 1 - 1 / T
    q <- (1 - p^(1/opy/exi)) / pat
    q <- loc - scale / shape * (1 - q^(-shape))
    return(q)
  }
  
  if (missing(opy)){
    warning("Argument ``opy'' is missing. Setting it to 365.25.")
    opy <- 365.25
  }
  
  if (missing(main)) main <- 'Return Level Plot'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- 500

  eps <- 10^(-3)
  plot(pot.fun, from = 1 + eps, to = xlimsup, log = "x",
       xlab = xlab, ylab = ylab, main = main, ...)

  invisible(pot.fun)
}
  
