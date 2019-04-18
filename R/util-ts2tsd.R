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

ts2tsd <- function(ts, d, vol = FALSE, method = "linear"){
  ##ts: a time serie i.e. a matrix/data.frame with two columns ``obs''
  ##    and ``time''
  ## d: the duration for the average windows.

  if (any(is.na(ts)))
    warning("NA's are not allowed in object ``ts''.\n
Replacing them by -1e6 !!!")

  if ( !any(colnames(ts) == "obs") )
    stop("``data'' should have a column named ``obs''...")

  if ( !any(colnames(ts) == "time") )
    stop("``data'' should have a column named ``time''...")
  

  tim <- ts[,"time"]
  obs <- ts[,"obs"]

  n <- length(obs)

  ##tim.start <- apply(cbind(tim[1], tim - d/2), 1, max)
  ##tim.end <- apply(cbind(tim[n], tim + d/2), 1, min)
  tim.start <- tim - d/2
  tim.end <- tim + d/2
  
  obs.start <- approx(tim, obs, xout = tim.start,
                      method = method)$y
  obs.end <- approx(tim, obs, xout = tim.end,
                    method = method)$y

  ##Replace NAs by -1e6s as ts2tsd C code does not
  ##accept NA values
  obs.start[is.na(obs.start)] <- -1e6
  obs.end[is.na(obs.end)] <- -1e6
  obs[is.na(obs)] <- -1e6 

  obs <- .C(POT_do_ts2tsd, as.double(tim), as.double(obs),
            as.double(tim.start), as.double(tim.end),
            as.double(obs.start), as.double(obs.end),
            as.integer(n), ans = double(n))$ans

  if (!vol)
    obs <- obs / d
  
  return(cbind(time = tim, obs = obs))
}
