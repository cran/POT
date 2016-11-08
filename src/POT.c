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

void gpdlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns)
{
  int i;
  double *dvec;
  
  dvec = (double *)R_alloc(*n, sizeof(double));

  if(*scale <= 0) {
     *dns = -1e6;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[i]) / *scale;
    if (data[i] <= 0) {
      *dns = -1e6;
      return;
    }
    if(fabs(*shape) <= 1e-6){
      *shape = 0;
      dvec[i] = -log(*scale) - data[i];
    }
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvec[i] = -log(*scale) - (1 / *shape + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];
}

void pplik(double *data, int *n, double *loc, double *scale,
	   double *shape, double *thresh, double *noy, double *dns)
{
  int i;
  double *dvec, preg;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  if(*scale <= 0) {
     *dns = -1e6;
     return;
  }

  preg = (*thresh - *loc) / *scale;

  if (*shape == 0)
    preg = - *noy * exp(-preg);

  else {

    preg = 1 + *shape * preg;

    if ((preg <= 0) && (*shape > 0)){
      *dns = -1e6;
      return;
    }

    else {
      preg = fmax2(preg, 0);
      preg = - *noy * R_pow(preg, -1 / *shape);
    }
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    
    if(*shape == 0)
      dvec[i] = log(1 / *scale) - data[i];
      
    else {
      data[i] = 1 + *shape * data[i];
      
      if(data[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvec[i] = log(1 / *scale) - (1 / *shape + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];

  
  *dns = *dns + preg; 

}

void samlmu(double *x, int *nmom, int *n, double *lmom){

  int i, j, *temp;
  double *p1, *p, *p2;

  temp = (int *)R_alloc(*n, sizeof(int));
  p1 = (double *)R_alloc(*n, sizeof(double));
  p = (double *)R_alloc(*n, sizeof(double));
  p2 = (double *)R_alloc(*n, sizeof(double));
  
  temp[0] = 1 - *n;
  p1[0] = 1;
  p[0] = -1;

  for(i = 1 ; i < *n ; i++){
    temp[i] = temp[i-1] + 2;
    p1[i] = 1;
    p[i] = (double) temp[i] / (double) ( *n - 1);
  }

  for(i = 0 ; i < *n; i++){
    lmom[0] = lmom[0] + x[i] / (double) *n;
    lmom[1] = lmom[1] + x[i] * p[i] / (double) *n;
  }

  for (i = 2 ; i < *nmom ; i++){
    for (j = 0 ; j < *n ; j++){
      p2[j] = p1[j];
      p1[j] = p[j];
      p[j] = ( (2*(i+1)-3)*temp[j]*p1[j] - (i-1) * (*n + i - 1)  * p2[j] ) /
	(double) (i * (*n - i));
      lmom[i] = lmom[i] + x[j] * p[j] / (double) (*n * lmom[1]);
    }
  }
}
  
