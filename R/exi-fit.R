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



fitexi <- function(data, threshold){
  
   
  if (any(is.na(data))){
    warning("NA's are not allowed in object ``data''.\nReplacing them by the threshold !!!")
    data[is.na(data)] <- threshold
  }

  idx <- which(data > threshold)
  nat <- length(idx)
  interTim <- diff(idx)

  if (max(interTim) == 1)
    exi <- 0

  else{
    if (max(interTim) <= 2){
      exi <- 2 * sum(interTim - 1)^2 / (nat - 1) /
        sum((interTim - 1) * (interTim - 2)) 
      exi <- min(1, exi)
    }
    
    else{
      exi <- 2 * sum(interTim)^2 / (nat - 1) /
        sum(interTim^2)
      exi <- min(1, exi)
    }
  }
    
  C <- floor(exi * nat) + 1
  sortInterTim <- sort(interTim, decreasing = TRUE)

  if (C <= length(interTim))
    TC <- sortInterTim[C]

  else
    TC <- max(interTim)
  
  return(list(exi = exi, tim.cond = TC))
}

  
