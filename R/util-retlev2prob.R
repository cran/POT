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



##Little functions that computes return period from non exceedance
##probability and viceversa
rp2prob <- function(retper, npy){
  ##retper   : the return period
  ##npy : the mean Number of events Per Year
  
  if (any(npy <=0))
    stop("``npy'' must be positive !!!")
  if (any(retper < 1/npy))
    stop("return period incompatible with ``npy'' !!!")
  prob <- 1 - 1/(npy * retper)
  
  tab <- cbind(npy = npy, retper = retper, prob = prob)
  rownames(tab) <- 1:length(tab[,1])
  return(tab)
}
prob2rp <- function(prob, npy){
  ##prob   : the probability of non exceedance
  ##npy    : the mean (N)umber of events (P)er (Y)ear
  
  if (any(npy <=0))
    stop("``npy'' must be positive !!!")
  if (any(prob <0) | any(prob >= 1) )
    stop("``prob'' must be in [0,1) !!!")
  
  retper <- 1 / (npy * (1 - prob))
  
  tab <- cbind(npy = npy, retper = retper, prob = prob)
  rownames(tab) <- 1:length(tab[,1])
  return(tab)
} 


