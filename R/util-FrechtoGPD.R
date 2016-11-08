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

##A function to transform GPD observations to a unit Frechet scale
gpd2frech <- function(x, loc = 0, scale = 1, shape = 0, pat = 1){
  z <- pgpd(x, loc, scale, shape, lower.tail = FALSE)
  z <- -1 / log(1 - pat * z)
  return(z)
}

##A function to transform unit Frechet observations to GPD ones
frech2gpd <- function(z, loc = 0, scale = 1, shape = 0, pat = 1){
  ##First, convert to uniform(0,1)
  x <- exp(-1/z) / pat
  x <- qgpd(x, loc, scale, shape)
  return(x)  
}
