#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################

clust <- function (data, u, tim.cond = 1, clust.max = FALSE, plot = FALSE, 
                   only.excess = TRUE, ...){
  if (!any(colnames(data) == "obs")) 
    stop("``data'' should have a column named ``obs''...")
  if (!any(colnames(data) == "time")) 
    stop("``data'' should have a column named ``time''...")
  
  obs <- data[, "obs"]
  tim <- data[, "time"]
  
  if (any(is.na(obs))){
    warning("NA's are not allowed in object ``data''.\nReplacing them by -1e6 !!!")
    obs[is.na(obs)] <- -1e+06
  }
  
  n <- as.integer(length(obs))
  
  if (all(obs <= u)) 
    stop("No data above the threshold !!!")
  clust <- .C("clust", n, as.double(obs),
              as.double(tim), as.double(tim.cond), 
              as.double(u), clust = double(2 * n),
              PACKAGE = "POT")$clust
  clust <- clust[clust != 0]
  clust <- matrix(clust[!is.na(clust)], ncol = 2, byrow = TRUE)
  colnames(clust) <- c("start", "end")

  n.clust <- length(clust[, 1])
  n.excess <- sum(obs > u, na.rm = TRUE)
  
  ##Replace NA values in ``obs''
  obs[obs == -1e+06] <- NA
  
  if (clust.max) {
    idx <- NULL
    for (i in 1:n.clust) {
      temp <- which.max(obs[clust[i, 1]:clust[i, 2]]) - 
        1
      idx <- c(idx, clust[i, 1] + temp)
    }
    events <- cbind(time = tim[idx], obs = obs[idx], idx = idx)
    rownames(events) <- 1:n.clust
  }
  else {
    events <- list()
    for (i in 1:n.clust)
      events <- c(events,
                  list(rbind(time = tim[clust[i, 1]:clust[i, 2]],
                             obs = obs[clust[i, 1]:clust[i, 2]])))
    
    names(events) <- paste("cluster", 1:n.clust)
  }
  
  exi <- n.clust/n.excess
  attributes(events)$exi <- exi
  if (plot) {
    plot(tim, obs, type = "n", ...)
    eps <- min(tim[clust[-1, "start"]] - tim[clust[-n.clust, 
                                                   "end"]])/2
    rect(tim[clust[, "start"]] - eps, rep(min(obs, na.rm = TRUE), 
                                          n.clust),
         tim[clust[, "end"]] + eps, rep(max(obs, na.rm = TRUE),
                                        n.clust), col = "lightgrey")
    if (only.excess) {
      idx.excess <- which(obs > u)
      tim <- tim[idx.excess]
      obs <- obs[idx.excess]
    }
    points(tim, obs)
  }
  return(events)
}




  
