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



## First, we need a function which computes the samples L-moments

samlmu <- function (x, nmom = 4, sort.data = TRUE)
{
  xok <- x[!is.na(x)]
  n <- length(xok)
  if (nmom <= 0) return(numeric(0))
  if (nmom <= 2) rnames <- paste("l", 1:nmom, sep = "_")
  else rnames <- c("l_1", "l_2", paste("t", 3:nmom, sep = "_"))
  lmom <- rep(NA, nmom)
  names(lmom) <- rnames
  if (n == 0) return(lmom)
  if (sort.data == TRUE) xok <- sort(xok)
  nmom.actual <- min(nmom, n)
  
  lmom <- .C("samlmu", as.double(xok), as.integer(nmom.actual),
             as.integer(n), lmom = double(nmom.actual),
             PACKAGE = "POT")$lmom
  names(lmom) <- rnames
  return(lmom)
}



