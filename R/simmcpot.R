simmcpot <- function(fitted, plot = TRUE, ...){

  if (all(class(fitted) != "mcpot"))
    stop("``fitted'' must be of class ``mcpot''")

  thresh <- fitted$threshold
  scale <- fitted$param["scale"]
  shape <- fitted$param["shape"]
  alpha <- fitted$param["alpha"]
  pat <- fitted$pat


  n <- length(fitted$data)

  args <- list(n = n, alpha = alpha)
  prob <- do.call(simmc, args)

  mcpot <- rep(NA, n)
  idx <- which(prob > (1 - pat))

  prob[idx] <- (prob[idx] - (1 - pat)) / pat
  mcpot[idx] <- qgpd(prob[idx], thresh, scale, shape)

  if (plot)
    plot(mcpot, ...)
  
  return(mcpot)
}

