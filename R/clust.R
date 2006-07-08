clust <- function(data, u, tim.cond = 1, clust.max = FALSE,
                  plot = FALSE, only.excess = TRUE, ...){
  
  if ( !any(colnames(data) == "obs") )
    stop("``data'' should have a column named ``obs''...")

  if ( !any(colnames(data) == "time") )
    stop("``data'' should have a column named ``time''...")
  
  obs <- data[,"obs"]
  tim <- data[,"time"]

  if (all(obs <= u))
    stop("No data above the threshold !!!")
  
  idx.excess <- which(obs > u)
  n.excess <- length(idx.excess)
  
  diff.tim <- diff(tim[idx.excess])
  
  clust <- .C("clust", as.integer(n.excess), as.double(idx.excess),
              as.double(diff.tim), as.double(tim.cond),
              clust = double(2*n.excess), PACKAGE = "POT")$clust
  
  clust <- clust[clust != 0]
  
  clust <- matrix(clust[!is.na(clust)], ncol = 2,
                  byrow = TRUE)
  colnames(clust) <- c("start","end")
  n.clust <- length(clust[,1])
  
  if (clust.max){
    idx <- NULL
    for (i in 1:n.clust){
      temp <- which.max(obs[clust[i,1]:clust[i,2]]) - 1
      idx <- c(idx, clust[i,1] + temp)
    }
    
    events <- cbind(time = tim[idx], obs = obs[idx], idx = idx)
    rownames(events) <- 1:n.clust
  }

  else{
    events <- list()
    for (i in 1:n.clust)
      events <- c(events, list(rbind(time = tim[clust[i,1]:clust[i,2]],
                                     obs = obs[clust[i,1]:clust[i,2]])))

    names(events) <- paste("cluster", 1:n.clust)
  }

  exi <- n.clust / n.excess
  attributes(events)$exi <- exi

  if (plot) {
    plot(tim, obs, type = "n", ...)

    eps <- min(tim[clust[-1, "start"]] -
               tim[clust[-n.clust, "end"]]) / 2

    rect(tim[clust[,"start"]] - eps, rep(min(obs, na.rm = TRUE),
                                         n.clust),
         tim[clust[,"end"]] + eps, rep(max(obs, na.rm = TRUE),
                                       n.clust),
         col = "lightgrey")

    if (only.excess){
      tim <- tim[idx.excess]
      obs <- obs[idx.excess]
    }
    points(tim, obs)
  }
  
  return(events)
}



exiplot <- function(data, u.range, tim.cond = 1, n.u = 50,
                    xlab, ylab, ...){

  if (all(data <= u.range[1]))
    stop("No data above u.range[1] !!! Please specify suitable value.")

  if (missing(xlab))
    xlab <- "Threshold"

  if (missing(ylab))
    ylab <- "Extremal Index"

  if (missing(u.range))
    u.range <- c(min(data[,"obs"]), quantile(data[,"obs"], 0.9))
  
  u.range <- seq(u.range[1], u.range[2], length.out = n.u)

  exi <- rep(NA, n.u)
  
  for (i in 1:n.u)
    exi[i] <- attributes(clust(data, u.range[i], tim.cond = tim.cond))$exi

  plot(u.range, exi, xlab = xlab, ylab = ylab, ...)

  invisible(cbind(thresh = u.range, exi = exi))

}
    
  
