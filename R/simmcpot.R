simmcpot <- function(fitted, plot = TRUE, ...){

  if (all(class(fitted) != "mcpot"))
    stop("``fitted'' must be of class ``mcpot''")

  thresh <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  alpha <- fitted$param["alpha"]
  pat <- fitted$pat
  model <- fitted$model
  n <- length(fitted$data)

  if (model %in% c("log", "nlog"))
    param <- list(alpha = alpha)

  if (model %in% c("alog", "anlog"))
    param <- list(alpha = alpha, asCoef1 = fitted$param["asCoef1"],
                  asCoef2 = fitted$param["asCoef2"])

  if (model == "mix")
    param <- list(alpha = alpha, asCoef = 0)
  
  if (model == "amix")
    param <- list(alpha = alpha, asCoef = fitted$param["asCoef"])

  param <- c(param, list(n = n, model = model))

  prob <- do.call(simmc, param)

  mcpot <- rep(NA, n)
  idx <- which(prob > (1 - pat))

  prob[idx] <- (prob[idx] - (1 - pat)) / pat
  mcpot[idx] <- qgpd(prob[idx], thresh, scale, shape)

  if (plot)
    plot(mcpot, ...)
  
  return(mcpot)
}

