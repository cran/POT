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

## This file contains several function to select the threshold
## for which the asymptotical approximation of Peaks Over a
## Threshold by a GP distribution is quite good.

diplot <- function(data, u.range, main, xlab, ylab,
                   nt = max(200,nrow(data)), conf=0.95,
                   ...){

  if ( !any(colnames(data) == "obs") )
    stop("``data'' should have a column named ``obs''...")

  if ( !any(colnames(data) == "time") )
    stop("``data'' should have a column named ``time''...")

  data <- na.omit(data)
  date <- data[,"time"]
  samp <- data[,"obs"]
                                       
  if (length(samp)<5){
    stop('Not enough data for a Dispersion Index Plot')
  }

  M <- diff(range(date))

  if (missing(u.range)) u.range <- c(min(samp),max(samp[-(1:4)]))

  thresh <- seq(u.range[1],u.range[2], length = nt)

  DI <- NULL

  date <- floor(date)
  tim.rec <- range(date)    
 
  for (u in thresh){

    nb.occ <- NULL
    idx.excess <- samp > u
    lambda <- sum(idx.excess) / M
   
    for (year in tim.rec[1]:tim.rec[2])
      nb.occ <- c(nb.occ, sum(idx.excess &
                              (date == year)))
      
    DI <- c(DI, var(nb.occ)/lambda)

  }

  conf_sup <- qchisq(1-(1-conf)/2,M-1)/(M-1)
  conf_inf <- qchisq((1-conf)/2,M-1)/(M-1)

  if ( missing(main) ) main <- 'Dispersion Index Plot'
  if ( missing(xlab) ) xlab <- 'Threshold'
  if ( missing(ylab) ) ylab <- 'Dispersion Index'
  
  plot(c(thresh,thresh[1]), c(DI, conf_sup), xlab=xlab, ylab=ylab,
       type='n', main = main, ...)
  rect(0, conf_inf, 2*u.range[2], conf_sup, col= 'lightgrey', border = FALSE)
  lines(thresh, DI)
  return(invisible(list(thresh=thresh,DI=DI)))

}

