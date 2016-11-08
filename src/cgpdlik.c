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

#include "header.h"

void cgpdlik(double *data, int *n, int *nn, double *lambda, double *loc,
	     double *scale, double *shape, double *dns){
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*nn, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
    *dns = -1e6;
    return;
  }

  for(i=0;i<*nn;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvec[i] = log(*lambda) - log(*scale) - data[i];
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvec[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data[i]);
    }
    //printf("dvec[%i] = %f\n", i, dvec[i]);
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns + dvec[i];

  //Now add the censored contribution to loglikelihood
  if (*n != *nn)
    *dns = *dns + (*n - *nn) * log(1 - *lambda);
   
}
