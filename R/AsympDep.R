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

  eps <- .Machine$double.eps^.5

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
  
tailind.test <- function(data, c = -0.1, emp.trans = TRUE,
                         chisq.n.class = 4){

  if (ncol(data) != 2)
    stop("data must be a matrix with two columns")

  if (emp.trans){
    ##We have to transform data to reverse exponential
    ##margins using empirical transformation
    n <- nrow(data)
    data <- apply(data, 2, rank) / (n + 1)
    data <- log(data)
  }

  else
    cat("``data'' is supposed to be reverse exponential distributed\n")

  ##The margin sum
  Msum <- apply(data, 1, sum)
  idx <- which(Msum > c)
  Msum <- Msum[idx]
  n.c <- length(Msum)
  V <- Msum / c


  ######################################################
  ##                                                  ##
  ##          The Neyman-Pearson Test                 ##
  ##                                                  ##
  ######################################################
  
  NP.stat <- -log(prod(2*V))
  if (n.c <= 170){
    NP.pval <- 1
    for (j in 1:(n.c-1))
      NP.pval <- NP.pval + (-2 * sum(log(V)))^j / prod(j:2)
  
    NP.pval <- exp(2*sum(log(V))) * NP.pval
  }

  else
    NP.pval <- pnorm((2 * sum(log(V)) + n.c) / sqrt(n.c))
  
  

  ######################################################
  ##                                                  ##
  ##        Fisher's Kappa Statistic Test             ##
  ##                                                  ##
  ######################################################

  U <- (1 - (1 - Msum) * exp(Msum)) /
    (1 - (1 - c) * exp(c))
  Usort <- c(0, sort(U))
  
  Fish.stat <- (n.c + 1) * max(c(diff(Usort), 1 - max(Usort)))

  
  Fish.pval <- 1

  for (j in 1:n.c)
    Fish.pval <- Fish.pval + (-1)^j * prod(2:(n.c + 1)) /
      prod(1:j) / prod(2:(n.c + 1 - j)) *
        max(0, 1 - j * Fish.stat / (n.c  + 1))^n.c 
  
  Fish.pval <- Fish.pval + (-1)^(n.c + 1) *
    max(0, 1 - Fish.stat)^n.c 
  Fish.pval <- 1 - Fish.pval

  
  ######################################################
  ##                                                  ##
  ##            Kolmogorov-Smirnov Test               ##
  ##                                                  ##
  ######################################################


  KStest <- ks.test(U, "punif")

  KS.pval <- KStest$p.value
  KS.stat <- KStest$statistic

  
  ######################################################
  ##                                                  ##
  ##                 Chi-Square Test                  ##
  ##                                                  ##
  ######################################################

  classes <- seq(0, 1, length.out = chisq.n.class + 1)

  m <- NULL
  for (bsup in classes[-1])
    m <- c(m, sum(U < bsup))

  m[-1] <- diff(m)
  p.theo <- diff(classes)

  if (any( (n.c * p.theo) < 6))
    warning("Classes for the Chi-Square test are ill-defined")
      
  ChiSq <- chisq.test(m, p = p.theo)
  ChiSq.pval <- ChiSq$p.value
  ChiSq.stat <- ChiSq$statistic
  

  p.values <- c(NP.pval, Fish.pval,KS.pval, ChiSq.pval)
  stats <- c(NP.stat, Fish.stat, KS.stat, ChiSq.stat)
  names(p.values) <- names(stats) <- c("NP", "Fish", "KS",
                                       "ChiSq")

  return(list(stats = cbind(stats, p.values), idx = idx))
  
}

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
