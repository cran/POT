print.uvpot <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("Estimator:", x$est, "\n")

  if (x$est == "MGF")
    cat("Statistic:", x$stat, "\n")
  
  if (x$est == 'MLE'){
    cat(" Deviance:", x$deviance, "\n")
    cat("      AIC:", AIC(x), "\n")
  }

  if (x$est == 'MPLE'){
    cat("\nPenalized Deviance:", x$deviance, "\n")
    cat("     Penalized AIC:", AIC(x), "\n")
  }
  
  cat("\nVarying Threshold:", x$var.thresh, "\n")
  
  if(!x$var.thresh)
    x$threshold <- x$threshold[1]
  
  cat("\n  Threshold Call:", x$threshold.call, "\n")
  cat("    Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")
  
  cat("\nEstimates\n") 
  print.default(format(x$fitted.values, digits = digits), print.gap = 2, 
                quote = FALSE)
  if(!is.null(x$std.err)) {
    cat("\nStandard Error Type:", x$std.err.type, "\n")
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

print.bvpot <- function (x, digits = max(3, getOption("digits") -
                              3), ...) {
  cat("\nCall:", deparse(x$call), "\n")
  cat("Estimator:", x$est, "\n")

  if (x$model == "log")
    model <- "Logistic"
  if (x$model == "alog")
    model <- "Asymetric Logistic"
  if (x$model == "nlog")
    model <- "Negative Logistic"
  if (x$model == "anlog")
    model <- "Asymetric Negative Logistic"
  if (x$model == "mix")
    model <- "Mixed"
  if (x$model == "amix")
    model <- "Asymetric Mixed"
  
  cat("Dependence Model and Strenght:\n")
  cat("\tModel :", model, "\n")
  cat("\tlim_u Pr[ X_1 > u | X_2 > u] =", round(x$chi, 3),
      "\n")

  if (x$est == "MLE"){
    cat("Deviance:", x$deviance, "\n")
    cat("     AIC:", AIC(x), "\n")
  }
  cat("\nMarginal Threshold:", round(x$threshold, digits), "\n")
  cat("Marginal Number Above:", x$nat[1:2], "\n")
  cat("Marginal Proportion Above:", round(x$pat[1:2], digits), "\n")
  cat("Joint Number Above:", x$nat[3], "\n")
  cat("Joint Proportion Above:", round(x$pat[3], digits),"\n")
  cat("Number of events such as {Y1 > u1} U {Y2 > u2}:",
      x$nat[4], "\n")
  cat("\nEstimates\n")
  print.default(format(fitted(x), digits = digits), print.gap = 2, 
                quote = FALSE)
  if (!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if (!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if (!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("\tConvergence:", x$convergence, "\n")
  cat("\tFunction Evaluations:", x$counts["function"], "\n")
  if (!is.na(x$counts["gradient"])) 
    cat("\tGradient Evaluations:", x$counts["gradient"], 
        "\n")
  if (!is.null(x$message)) 
    cat("\nMessage:", x$message, "\n")
  cat("\n")
}

logLik.pot <- function(object, ...){
  llk <- object$logLik
  attr(llk, "df") <- length(fitted(object))
  class(llk) <- "logLik"
  return(llk)
}

print.mcpot <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("\nCall:", deparse(x$call), "\n")
  cat("Estimator:", x$est, "\n")

  if (x$model == "log")
    model <- "Logistic"
  if (x$model == "alog")
    model <- "Asymetric Logistic"
  if (x$model == "nlog")
    model <- "Negative Logistic"
  if (x$model == "anlog")
    model <- "Asymetric Negative Logistic"
  if (x$model == "mix")
    model <- "Mixed"
  if (x$model == "amix")
    model <- "Asymetric Mixed"
  
  cat("Dependence Model and Strenght:\n")
  cat("\tModel :", model, "\n")
  cat("\tlim_u Pr[ X_1 > u | X_2 > u] =", round(x$chi, 3),
      "\n")

  if (x$est == 'MLE'){
    cat("Deviance:", x$deviance, "\n")
    cat("     AIC:", AIC(x), "\n")
  }
  
  cat("\nThreshold Call:", x$threshold.call, "\n")
  cat("Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")
  
  cat("\nEstimates\n") 
  print.default(format(fitted(x), digits = digits), print.gap = 2, 
                quote = FALSE)
  if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

convassess.uvpot <- function(fitted, n = 50){

  if (!("uvpot" %in% class(fitted)))
    stop("``fitted'' must be of class ``uvpot''")

  if (!(fitted$est %in% c("MLE","LME","MPLE","MEDIANS","MDPD",
                          "MGF")))
    stop("\n    There's no objective function to optimize!")

  nat <- fitted$nat
  thresh <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  data <- fitted$data
  nobs <- length(data)

  if (fitted$est == "MLE")
    fun <- function()
      fitgpd(data, thresh, "mle", start = start,
             std.err.type = "none")

  if (fitted$est == "LME"){
    fun <- function()
      fitgpd(data, thresh, "lme", start = start)
    startsLME <- seq(-2, 0, length = n)
  }

  if (fitted$est == "MEDIANS")
    fun <- function()
      fitgpd(data, thresh, "med", start = start)
      
  if (fitted$est == "MPLE")
    fun <- function()
      fitgpd(data, thresh, "mple", start = start,
             std.err.type = "none")

  if (fitted$est == "MDPD")
    fun <- function() 
      fitgpd(data, thresh, "mdpd", start = start)

  if (fitted$est == "MGF")
    fun <- function()
      fitgpd(data, thresh, "mgf", stat = fitted$stat,
             start = start)

  est <- startValues <- matrix(NA, ncol = 2, nrow = n)
  colnames(est) <- colnames(startValues) <- c("scale", "shape")
  optValues <- rep(NA, n)
  
  for (i in 1:n){
    if (fitted$est == "LME")
      start <- list(x = startsLME[i])

    else{
      x <- sample(data, size = nobs, replace = TRUE)
      startValues[i,] <- fitgpd(x, thresh, "pwmu", hybrid = TRUE)$param
      start <- list(scale = startValues[i,"scale"],
                    shape = startValues[i,"shape"])
    }

    fit <- fun()
    optValues[i] <- fit$opt.value

    est[i,] <- fit$param
  }

  idx <- which(optValues == 1e6)
  optValues[idx] <- NA
  par(mfrow=c(2,2))

  if(fitted$est == "LME")
    plot(startsLME, xlab = "Index", ylab = expression(b[0]),
         main = "Starting Values")
  else{
    plot(startValues, xlab = "Scale", ylab = "Shape", main = "Starting Values",
         type = "n")
    
    
    if (length(idx) > 0){
      points(startValues[-idx,])
      points(startValues[idx,], col = "red", pch = 15)
    }
    
    else
      points(startValues)
  }
  
  plot(1:n, optValues, xlab = "Index", ylab = "Objective", type = "n",
       main = "Objective Value Trace Plot")
  abline(h=fitted$opt.value, col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], optValues[-idx])
    points(idx, rep(fitted$opt.value, length(idx)), col = "red",
           pch = 15)
  }

  else
    points(optValues)


  plot(1:n, est[,1], xlab = "Index", ylab = "Scale", type = "n",
       main = "Scale Trace Plot")
  abline(h=fitted$param["scale"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,1])
    points(idx, est[idx,1], col = "red", pch = 15)
  }

  else
    points(est[,1])

  plot(1:n, est[,2], xlab = "Index", ylab = "Shape", type = "n",
       main = "Shape Trace Plot")
  abline(h=fitted$param["shape"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,2])
    points(idx, est[idx,2], col = "red", pch = 15)
  }

  else
    points(est[,2])


  invisible(list(start.values = startValues, opt.values = optValues,
                 param = est))
    
}

convassess.bvpot <- function(fitted, n = 50){

  if (!("bvpot" %in% class(fitted)))
    stop("``fitted'' must be of class ``bvpot''")

  if (fitted$model != "log")
    return(cat("This function is only implemented for model ``log''.\nStill work to do...\n"))

  nat <- fitted$nat
  thresh <- fitted$threshold
  param <- fitted$param
  data <- fitted$data
  nobs <- nrow(data)
  model <- fitted$model
  
  fun <- function()
    fitbvgpd(data, thresh, start = start, model = model,
             std.err.type = "none")

  est <- startValues <- matrix(NA, ncol = 5, nrow = n)
  colnames(est) <- colnames(startValues) <- c("scale1", "shape1",
                                              "scale2", "shape2",
                                              "alpha")

  optValues <- rep(NA, n)
  
  for (i in 1:n){

    idx <- sample(1:nobs, size = nat, replace = TRUE)
    x <- data[idx,]
    startValues[i,1:2] <- fitgpd(x[,1], thresh[1], "pwmu",
                                 hybrid = TRUE)$param
    startValues[i,3:4] <- fitgpd(x[,2], thresh[1], "pwmu",
                                 hybrid = TRUE)$param
    startValues[i,5] <- runif(1, 0.2, 0.8)
    start <- list(scale1 = startValues[i,"scale1"],
                  shape1 = startValues[i,"shape1"],
                  scale2 = startValues[i,"scale2"],
                  shape2 = startValues[i,"shape2"],
                  alpha = startValues[i,"alpha"])
    
    fit <- fun()
    optValues[i] <- fit$deviance
    
    est[i,] <- fit$param
  }

  idx <- which(optValues == 2e6)
  optValues[idx] <- NA
  
  par(mfrow=c(3,3))
  
  ##Starting Values Marge 1
  plot(startValues[,c("scale1","shape1")], xlab = "Scale Marge 1",
       ylab = "Shape Marge 1", main = "Starting Values Marge 1",
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,c("scale1","shape1")])
    points(startValues[idx,c("scale1","shape1")], col = "red",
           pch = 15)
  }
  
  else
    points(startValues[,c("scale1","shape1")])

  ##Starting Values Marge 2
  plot(startValues[,c("scale2","shape2")], xlab = "Scale Marge 2",
       ylab = "Shape Marge 2", main = "Starting Values Marge 2",
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,c("scale2","shape2")])
    points(startValues[idx,c("scale2","shape2")], col = "red",
           pch = 15)
  }
  
  else
    points(startValues[,c("scale2","shape2")])

  ##Starting Values For Alpha
  plot(startValues[,"alpha"], xlab = "Index",
       ylab = expression(alpha),
       main = expression(paste("Starting Values for ", alpha)),
       type = "n")
    
    
  if (length(idx) > 0){
    points(startValues[-idx,"alpha"])
    points(startValues[idx,"alpha"], col = "red", pch = 15)
  }
  
  else
    points(startValues[,"alpha"])

   ##Scale Marge 1 Trace Plot
  plot(1:n, est[,"scale1"], xlab = "Index", ylab = "Scale Marge 1", type = "n",
       main = "Scale Marge 1 Trace Plot")
  abline(h=fitted$param["scale1"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"scale1"])
    points(idx, est[idx,"scale1"], col = "red", pch = 15)
  }

  else
    points(est[,"scale1"])

  ##Shape Marge 1 Trace Plot
  plot(1:n, est[,"shape1"], xlab = "Index", ylab = "Shape Marge 1", type = "n",
       main = "Shape Marge 1 Trace Plot")
  abline(h=fitted$param["shape1"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"shape1"])
    points(idx, est[idx,"shape1"], col = "red", pch = 15)
  }

  else
    points(est[,"shape1"])

  ##Alpha Trace Plot
  plot(1:n, est[,"alpha"], xlab = "Index", ylab = expression(alpha), type = "n",
       main = expression(paste(alpha, " Trace Plot")))
  abline(h=fitted$param["alpha"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"alpha"])
    points(idx, est[idx,"alpha"], col = "red", pch = 15)
  }

  else
    points(est[,"alpha"])

  ##Scale Marge 2 Trace Plot
  plot(1:n, est[,"scale2"], xlab = "Index", ylab = "Scale Marge 2", type = "n",
       main = "Scale Marge 2 Trace Plot")
  abline(h=fitted$param["scale2"], col = "blue", lty = 2)

  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"scale2"])
    points(idx, est[idx,"scale2"], col = "red", pch = 15)
  }

  else
    points(est[,"scale2"])

  ##Shape Marge 2 Trace Plot
  plot(1:n, est[,"shape2"], xlab = "Index", ylab = "Shape Marge 2", type = "n",
       main = "Shape Marge 2 Trace Plot")
  abline(h=fitted$param["shape2"], col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], est[-idx,"shape2"])
    points(idx, est[idx,"shape2"], col = "red", pch = 15)
  }

  else
    points(est[,"shape2"])

  ##Deviance Trace Plot
  plot(1:n, optValues, xlab = "Index", ylab = "Deviance", type = "n",
       main = "Deviance Trace Plot")
  abline(h=fitted$deviance, col = "blue", lty = 2)
  
  if (length(idx) > 0){
    points((1:n)[-idx], optValues[-idx])
    points(idx, rep(fitted$deviance, length(idx)), col = "red",
           pch = 15)
  }

  else
    points(optValues)

  
  invisible(list(start.values = startValues, opt.values = optValues,
                 param = est))
    
}
