chimeas <- function(data, u.range, n.u = 500, xlab, ylabs, ci = 0.95,
                    show.bound = TRUE, which = 1:2, ask = nb.fig < length(which) &&
                    dev.interactive(), ..., col.ci = "grey", col.bound = "blue",
                    lty.ci = 1, lty.bound = 1){

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
    
  ##Compute the variance of Chi and Chibar (Delta Method)
  ##Chi and Chibar variance obtained using the TCL
  varChi <- (log(u) * probJu)^(-2) * probJu * (1 - probJu) / n
  varChibar <- (2 * log(1 - u) / (probJuBar * log(probJuBar)^2))^2 *
    probJuBar * (1 - probJuBar) / n

  chiBounds <- 2 - log(pmax(2*u - 1, 0)) / log(u)

  if (show[1]){
    plot(u, chi, type  ="l", xlab = xlab, ylab = ylabs[1], xlim = c(0,1),
         ylim = c(-1,1), ...)
    lines(u, chi - qnorm((1+ci)/2) * sqrt(varChi), col = col.ci)
    lines(u, chi + qnorm((1+ci)/2) * sqrt(varChi), col = col.ci)

    if (show.bound){
      lines(u, chiBounds, lty = lty.bound, col = col.bound)
      abline(h = 1, lty = lty.bound, col = col.bound)
    }
  }
  
  if (show[2]){
    plot(u, chibar, type  ="l", xlab = xlab, ylab = ylabs[2], xlim = c(0,1),
         ylim = c(-1,1), ...)
    lines(u, chibar - qnorm((1+ci)/2) * sqrt(varChibar), col = col.ci)
    lines(u, chibar + qnorm((1+ci)/2) * sqrt(varChibar), col = col.ci)

    if (show.bound)
      abline(h = c(1,-1), col = col.bound, lty = lty.bound)
  }
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
  Msum <- Msum[Msum > c]
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

  cat("Number of observations such as x + y > c:", n.c, "\n")
  return(cbind(stats, p.values))
  
}
