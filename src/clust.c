/*****************************************************************************
 *   Copyright (c) 2014 Mathieu Ribatet                                                                                                  
 *                                                                                                                                                                        
 *   This program is free software; you can redistribute it and/or modify                                               
 *   it under the terms of the GNU General Public License as published by                                         
 *   the Free Software Foundation; either version 2 of the License, or                                                   
 *   (at your option) any later version.                                                                                                            
 *                                                                                                                                                                         
 *   This program is distributed in the hope that it will be useful,                                                             
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
 *   GNU General Public License for more details.                                                                                    
 *                                                                                                                                                                         
 *   You should have received a copy of the GNU General Public License                                           
 *   along with this program; if not, write to the                                                                                           
 *   Free Software Foundation, Inc.,                                                                                                              
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
 *                                                                                                                                                                         
 *
 ******************************************************************************/

# include "header.h"

void do_clust(int *n, double *obs, double *tim, double *cond,
	   double *thresh, double *clust){

  int i, idx, idxMax;

  i = 0;

  while (i < *n){

    if (obs[i] > *thresh){
      clust[i] = i+1;
      idx = i + 1;
      idxMax = idx;

      while( (obs[idx] > *thresh) ||
	     ( (tim[idx]-tim[idxMax]) <= *cond)){
	if (obs[idx] >= obs[idxMax])
	  idxMax = idx;
	idx++;

	if (idx >= *n){
	  idx--;
	  break;
	}
      }
      
      clust[i+1] = idx;
      i = idx + 1;
    }

    else
      i++;    
  }
}
