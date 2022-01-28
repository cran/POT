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


lmomplot <- function(data, u.range, nt = max(50, length(data)),
                     identify = TRUE, ...){

  data <- sort(as.numeric(data))

  n <- length(data)
  
  if ( n < 5){
    stop('Not enougth data for a L-moments Plot.\nLower limit
 for the threshold is too high')
  }
  
  if (missing(u.range)) {
    
    u.range <- c(data[1], data[n - 4])
    u.range <- u.range - sqrt(.Machine$double.eps)
    
  }

  data <- data[ data > u.range[1] ]
  
  if (all(data <= u.range[2])) 
    stop("upper limit for threshold is too high")
  
  thresh <- seq(u.range[1],u.range[2],
                length= nt)

  point <- NULL
  for ( u in thresh){

    lmoments34 <- samlmu(data[data > u ], sort.data = FALSE)[3:4] 
    point <-  cbind(point,lmoments34)

  }

  courbe_theo <- function(Tau3){
    Tau3*(1+5*Tau3)/(5+Tau3)
  }

  plot(courbe_theo, main='L-Moments Plot', xlab=expression(tau[3]),
       ylab=expression(tau[4]), col='grey',lwd=2,...)
  lines(point[1,],point[2,],type='b')
  

  if (identify)  identify(point[1,],point[2,],
                          labels=round(thresh,digits=2),offset=1)
  
}
  
