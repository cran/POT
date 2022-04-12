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

void do_gpdbvlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	      double *lambda2, double *thresh, double *scale1,
	      double *shape1, double *scale2, double *shape2,
	      double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2, v12,
    censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha > 1 || *alpha < 0.05 || *scale1 < 0.01 ||
     *scale2 < 0.01){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
    
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2]) (but omitting the alpha power!!!)
    v = R_pow(z1[i], - 1 / *alpha) + R_pow(z2[i], - 1 / *alpha);
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the log negative partial derivative with
      //respect to the second component
      nv2 = -(1/ *alpha + 1) * log(z2[i]) +
	(*alpha - 1) * log(v);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = nv2 + nK2 - R_pow(v, *alpha);
      
    }
    
    if ((data1[i] > 0) && (data2[i] == 0)){ 

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the log negative partial derivative with
      //respect to the first component
      nv1 = -(1/ *alpha + 1) * log(z1[i]) +
	(*alpha - 1) * log(v);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvec[i] = nv1 + nK1 - R_pow(v, *alpha);
      
    }

    if ((data1[i] * data2[i]) > 0){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow(z1[i], -1/ *alpha - 1) * 
	R_pow(v, *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow(z2[i], -1/ *alpha - 1) * 
	R_pow(v, *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (1 - 1 / *alpha) * R_pow(z1[i]*z2[i],
				     -1/ *alpha - 1) *
	R_pow(v, *alpha  - 2);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- R_pow(v, *alpha);
      
    }
  }
      
  for (i=0;i<*nn;i++){
    *dns = *dns + dvec[i];
  }
  
  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = R_pow(*lambda1, -1 / *alpha) +
      R_pow(*lambda2, - 1/ *alpha);
    censCont = -R_pow(censCont, *alpha);
    *dns = *dns + (*n - *nn) * censCont;
  }
}
    
void do_gpdbvalog(double *data1, double *data2, int *n, int *nn,
	       double *lambda1, double *lambda2, double *thresh,
	       double *scale1, double *shape1, double *scale2,
	       double *shape2, double *alpha, double *asCoef1,
	       double *asCoef2, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2, v12,
    censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha > 1 || *alpha < 0.05 || *scale1 < 0.01 ||
     *scale2 < 0.01 || *asCoef1 < 0 || *asCoef1 > 1 ||
     *asCoef2 < 0 || *asCoef2 > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
  }

    
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = (1 - *asCoef1) / z1[i] + (1 - *asCoef2) / z2[i] +
      R_pow(R_pow(*asCoef1 / z1[i], 1 / *alpha) + 
	    R_pow(*asCoef2 / z2[i], 1 / *alpha), *alpha);
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = (1 - *asCoef2) / R_pow_di(z2[i], 2) +
	R_pow(*asCoef2, 1/ *alpha) *
	R_pow(z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = log(nv2) + nK2 - v;
     	
    }
    
    if ((data1[i] > 0) && (data2[i] == 0)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = (1 - *asCoef1) / R_pow_di(z1[i], 2) +
	R_pow(*asCoef1, 1/ *alpha) * 
	R_pow(z1[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvec[i] = log(nv1) + nK1 - v;
      
    }

    if ((data1[i] * data2[i]) > 0){
      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = (1 - *asCoef1) / R_pow_di(z1[i], 2) +
	R_pow(*asCoef1, 1/ *alpha) * 
	R_pow(z1[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = (1 - *asCoef2) / R_pow_di(z2[i], 2) +
	R_pow(*asCoef2, 1/ *alpha) *
	  R_pow(z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
	
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (1 - 1/ *alpha) * R_pow(*asCoef1 * *asCoef2, 1/ *alpha) *
	R_pow(z1[i]*z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 2);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12) - v;
      
    }
  }
        
  for (i=0;i<*nn;i++)
    *dns = *dns + dvec[i];

  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = R_pow(*lambda1 / *asCoef1, - 1 / *alpha) +
      R_pow(*lambda2 / *asCoef2, - 1 / *alpha);
    censCont = (*asCoef1 - 1) / *lambda1 + (*asCoef2 - 1) / *lambda2 -
      R_pow(censCont, *alpha);
    *dns = *dns + (*n - *nn) * censCont;
  }

}
    
void do_gpdbvnlog(double *data1, double *data2, int *n, int *nn,
	       double *lambda1, double *lambda2, double *thresh,
	       double *scale1, double *shape1, double *scale2,
	       double *shape2, double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2,
    v12, censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha < 0.01 || *alpha > 15 || *scale1 < 0.01 || *scale2 < 0.01){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
    
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = 1 / z1[i] + 1 / z2[i] -
      R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	    - 1 / *alpha);
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = log(nv2) + nK2 - v;
      
    }
    
    if ((data1[i] > 0) && (data2[i] == 0)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
	
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvec[i] = log(nv1) + nK1 - v;
      
    }

    if ((data1[i] * data2[i]) > 0){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	  1 / z2[i];
      
	//Compute the partial mixed derivative	
      v12 = -(*alpha + 1) * R_pow(z1[i] * z2[i], *alpha - 1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
		- 1 / *alpha - 2);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      
    }
  }
        
  for (i=0;i<*nn;i++)
    *dns = *dns + dvec[i];

  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = R_pow(*lambda1, *alpha) + R_pow(*lambda2, *alpha);
    censCont = - 1 / *lambda1 - 1 / *lambda2 +
      R_pow(censCont, -1 / *alpha);
    *dns = *dns + (*n - *nn) * censCont;
  }
}

void do_gpdbvanlog(double *data1, double *data2, int *n, int *nn,
		double *lambda1, double *lambda2, double *thresh,
		double *scale1, double *shape1, double *scale2,
		double *shape2, double *alpha, double *asCoef1,
		double *asCoef2, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2,
    v12, censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha < 0.2 || *alpha > 15 || *scale1 < 0.01 || *scale2 < 0.01 ||
     *asCoef1 < 0 || *asCoef1 > 1 || *asCoef2 < 0 || *asCoef2 > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
    
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = 1 / z1[i] + 1 / z2[i] -
      R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	    R_pow(z2[i] / *asCoef2, *alpha), - 1 / *alpha);
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(*asCoef2, - *alpha) *
	R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = log(nv2) + nK2 - v;
     	
    }
  
    if ((data1[i]>0) && (data2[i] == 0)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(*asCoef1, - *alpha) *
	R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvec[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] * data2[i]) > 0){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(*asCoef1, - *alpha) *
	R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(*asCoef2, - *alpha) *
	R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = -(*alpha + 1) * R_pow(*asCoef1 * *asCoef2, - *alpha) *
	R_pow(z1[i] * z2[i], *alpha - 1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha), - 1 / *alpha - 2);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      
    }
  }
      
  for (i=0;i<*nn;i++)
    *dns = *dns + dvec[i];

  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = R_pow(*lambda1 / *asCoef1, *alpha) +
      R_pow(*lambda2 / *asCoef2, *alpha);
    censCont = - 1 / *lambda1 - 1 / *lambda2 +
      R_pow(censCont, -1 / *alpha);
    *dns = *dns + (*n - *nn) * censCont;
  }
}

void do_gpdbvmix(double *data1, double *data2, int *n, int *nn,
	      double *lambda1, double *lambda2, double *thresh,
	      double *scale1, double *shape1, double *scale2,
	      double *shape2, double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2,
    v12, censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha > 1 || *alpha < 0 || *scale1 < 0.01 ||
     *scale2 < 0.01){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
    
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = R_pow_di(z1[i], - 1) + R_pow_di(z2[i], - 1) -
      *alpha / (z1[i] + z2[i]);
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = log(nv2) + nK2 - v;
      
    }
  
    if ((data1[i] > 0) && (data2[i] == 0)){
      
      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
	
      dvec[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] * data2[i]) > 0){
      
      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = -2 * *alpha * R_pow_di(z1[i] + z2[i], -3);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      
    }
  }
      
  for (i=0;i<*nn;i++)
    *dns = *dns + dvec[i];

  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = - 1 / *lambda1 - 1 / *lambda2 +
      *alpha / (*lambda1 + *lambda2);      
    *dns = *dns + (*n - *nn) * censCont;
  }
}

void do_gpdbvamix(double *data1, double *data2, int *n, int *nn,
	       double *lambda1, double *lambda2, double *thresh,
	       double *scale1, double *shape1, double *scale2,
	       double *shape2, double *alpha, double *asCoef,
	       double *dns){
  
  int i;

  double eps, *t1, *t2, *z1, *z2, *dvec, v, nv1, nK1, nv2, nK2,
    v12, c1, censCont;

  eps = R_pow(DBL_EPSILON, 0.3);
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  z1 = (double *)R_alloc(*nn, sizeof(double));
  z2 = (double *)R_alloc(*nn, sizeof(double));
  dvec = (double *)R_alloc(*nn, sizeof(double));
  
  
  if(*alpha < 0  || *scale1 < 0.01 || *scale2 < 0.01 ||
     *alpha + 2 * *asCoef > 1 || *alpha + 3 * *asCoef < 0 ||
     *alpha + *asCoef > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){

    //Margin 1
    t1[i] = (data1[i]  - thresh[0]) / *scale1;
        
    if (data1[i] <= thresh[0]){
      t1[i] = 1;
      data1[i] = 0;
    }

    else{

      if (fabs(*shape1) <= eps){
	*shape1 = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape1 * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape1);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - thresh[1]) / *scale2;

    if (data2[i] <= thresh[1]){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape2) <= eps){
	*shape2 = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape2 * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape2);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda1 * t1[i]);
    z2[i] = -1 / log(1 - *lambda2 * t2[i]);
    
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nn;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    //c1 is a common value for cases 1 - 4
    c1 = ((*alpha + *asCoef) * z1[i] + (*alpha + 2 * *asCoef) *
	  z2[i]) / R_pow_di(z1[i] + z2[i], 2);
    v = 1/z1[i] + 1/z2[i] - c1;
      
    if ((data1[i] == 0) && (data2[i] > 0)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + 2 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvec[i] = log(nv2) + nK2 - v;
      
    }
  
    if ((data1[i] > 0) && (data2[i] == 0)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvec[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] * data2[i]) > 0){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda1) - log(*scale1) + 
	(1 + *shape1) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + 2 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda2) - log(*scale2) + 
	(1 + *shape2) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (4 * *alpha + 6 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 3) - 6 * c1 /
	R_pow_di(z1[i] + z2[i], 2);
      
      dvec[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
    }
  }
  
  for (i=0;i<*nn;i++)
    *dns = *dns + dvec[i];
  
  //Now add the censored contribution to loglikelihood
  if (*nn != *n){
    *lambda1 = - 1 / log(1 - *lambda1);
    *lambda2 = - 1 / log(1 - *lambda2);
    censCont = ((*alpha + *asCoef) * *lambda1 +
		(*alpha + 2 * *asCoef) * *lambda2) /
      R_pow_di(*lambda1 + *lambda2, 2);
    censCont = censCont - 1 / *lambda1 - 1 / *lambda2;
    *dns = *dns + (*n - *nn) * censCont;
  }
}
