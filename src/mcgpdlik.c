#include "header.h"

void gpdmclog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	      double *scale, double *shape, double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
  
  
  if(*alpha > 1 || *alpha < 0.05 || *scale < 0.01 ||
     *scale < 0.01){
    *dns = -1e6;
    return;
  }
  

  //1) The joint part:

  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;
        
    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }

    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);

    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }

  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2]) (but omitting the alpha power!!!)
    v = R_pow(z1[i], - 1 / *alpha) + R_pow(z2[i], - 1 / *alpha);
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the log negative partial derivative with
      //respect to the second component
      nv2 = -(1/ *alpha + 1) * log(z2[i]) +
	(*alpha - 1) * log(v);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvecj[i] = nv2 + nK2 - R_pow(v, *alpha);
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){ 

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the log negative partial derivative with
      //respect to the first component
      nv1 = -(1/ *alpha + 1) * log(z1[i]) +
	(*alpha - 1) * log(v);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvecj[i] = nv1 + nK1 - R_pow(v, *alpha);
      
    }

    if ((data1[i] > *thresh) && (data2[i] > *thresh)){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow(z1[i], -1/ *alpha - 1) * 
	R_pow(v, *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow(z2[i], -1/ *alpha - 1) * 
	R_pow(v, *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (1 - 1 / *alpha) * R_pow(z1[i]*z2[i],
				     -1/ *alpha - 1) *
	R_pow(v, *alpha  - 2);
      
      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- R_pow(v, *alpha);
      
    }
  }
      
  for (i=0;i<*nnj;i++){
    *dns = *dns + dvecj[i];
  }
  
  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = 2 * R_pow(lambda2, -1 / *alpha);
    censCont = -R_pow(censCont, *alpha);
    *dns = *dns + (*nj - *nnj) * censCont;
  }

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}


void gpdmcalog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	       double *scale, double *shape, double *alpha, double *asCoef1,
	       double *asCoef2, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
      

  if(*alpha > 1 || *alpha < 0.05 || *scale < 0.01 ||
     *asCoef1 < 0 || *asCoef1 > 1 || *asCoef2 < 0 ||
     *asCoef2 > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);

    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }

  }

    
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = (1 - *asCoef1) / z1[i] + (1 - *asCoef2) / z2[i] +
      R_pow(R_pow(*asCoef1 / z1[i], 1 / *alpha) + 
	    R_pow(*asCoef2 / z2[i], 1 / *alpha), *alpha);
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = (1 - *asCoef2) / R_pow_di(z2[i], 2) +
	R_pow(*asCoef2, 1/ *alpha) *
	R_pow(z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      dvecj[i] = log(nv2) + nK2 - v;
     	
    }
    
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = (1 - *asCoef1) / R_pow_di(z1[i], 2) +
	R_pow(*asCoef1, 1/ *alpha) * 
	R_pow(z1[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      dvecj[i] = log(nv1) + nK1 - v;
      
    }

    if ((data1[i] > *thresh) && (data2[i] > *thresh)){
      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = (1 - *asCoef1) / R_pow_di(z1[i], 2) +
	R_pow(*asCoef1, 1/ *alpha) * 
	R_pow(z1[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = (1 - *asCoef2) / R_pow_di(z2[i], 2) +
	R_pow(*asCoef2, 1/ *alpha) *
	  R_pow(z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 1);
	
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (1 - 1/ *alpha) * R_pow(*asCoef1 * *asCoef2, 1/ *alpha) *
	R_pow(z1[i]*z2[i], -1/ *alpha - 1) *
	R_pow(R_pow(*asCoef1 / z1[i], 1/ *alpha) + 
	      R_pow(*asCoef2 / z2[i], 1/ *alpha), *alpha - 2);
      
      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12) - v;
      
    }
  }
        
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];

  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = 2 *R_pow(lambda2 / *asCoef1, - 1 / *alpha);
    censCont = (*asCoef1 - 1) / lambda2 + (*asCoef2 - 1) / lambda2 -
      R_pow(censCont, *alpha);
    *dns = *dns + (*nj - *nnj) * censCont;
  }

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}
    
void gpdmcnlog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	      double *scale, double *shape, double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
  

  if(*alpha < 0.01 || *alpha > 15 || *scale < 0.01){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);
    
    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = 1 / z1[i] + 1 / z2[i] -
      R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	    - 1 / *alpha);
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      if (!R_FINITE(log(nv2))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv2) + nK2 - v;
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
	
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];

      if (!R_FINITE(log(nv1))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv1) + nK1 - v;
      
    }

    if ((data1[i] > *thresh) && (data2[i] > *thresh)){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	  1 / z2[i];
      
	//Compute the partial mixed derivative	
      v12 = -(*alpha + 1) * R_pow(z1[i] * z2[i], *alpha - 1) *
	R_pow(R_pow(z1[i], *alpha) + R_pow(z2[i], *alpha),
		- 1 / *alpha - 2);
      
      if (!R_FINITE(log(nv1 * nv2 - v12))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      
    }
  }
        
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];

  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = 2 * R_pow(lambda2, *alpha);
    censCont = - 2/lambda2 + R_pow(censCont, -1 / *alpha);
    *dns = *dns + (*nj - *nnj) * censCont;
  }

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}


void gpdmcanlog(double *data1, double *data2, double *data3, int *nj,
		int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
		double *scale, double *shape, double *alpha, double *asCoef1,
		double *asCoef2, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
    
  if(*alpha < 0.2 || *alpha > 15 || *scale < 0.01 || 
     *asCoef1 < 0 || *asCoef1 > 1 || *asCoef2 < 0 ||
     *asCoef2 > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh){
      t2[i] = 1;
      data2[i] = 0;
    }

    else{

      if (fabs(*shape) <= eps)
	*shape = 0;

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);

    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }

  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = 1 / z1[i] + 1 / z2[i] -
      R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	    R_pow(z2[i] / *asCoef2, *alpha), - 1 / *alpha);
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(*asCoef2, - *alpha) *
	R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];

      if (!R_FINITE(log(nv2))){
      	*dns = -1e6;
      	return;
      }
            
      dvecj[i] = log(nv2) + nK2 - v;
     	
    }
  
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(*asCoef1, - *alpha) *
	R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      if (!R_FINITE(log(nv1))){
      	*dns = -1e6;
      	return;
      }
            
      dvecj[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] > *thresh)){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - R_pow(*asCoef1, - *alpha) *
	R_pow(z1[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - R_pow(*asCoef2, - *alpha) *
	R_pow(z2[i], *alpha -1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha),
	      -1 / *alpha - 1);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = -(*alpha + 1) * R_pow(*asCoef1 * *asCoef2, - *alpha) *
	R_pow(z1[i] * z2[i], *alpha - 1) *
	R_pow(R_pow(z1[i] / *asCoef1, *alpha) +
	      R_pow(z2[i] / *asCoef2, *alpha), - 1 / *alpha - 2);
      
      if (!R_FINITE(log(nv1 * nv2 - v12))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
    }
  }
      
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];

  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = R_pow(lambda2 / *asCoef1, *alpha) +
      R_pow(lambda2 / *asCoef2, *alpha);
    censCont = - 2/lambda2 + R_pow(censCont, -1 / *alpha);
    *dns = *dns + (*nj - *nnj) * censCont;
  }

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}



void gpdmcmix(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	      double *scale, double *shape, double *alpha, double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
    
  if(*alpha > 1 || *alpha < 0 || *scale < 0.01){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);
    
    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }

  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    v = R_pow_di(z1[i], - 1) + R_pow_di(z2[i], - 1) -
      *alpha / (z1[i] + z2[i]);
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      if (!R_FINITE(log(nv2))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv2) + nK2 - v;
      
    }
  
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){
      
      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
	
      if (!R_FINITE(log(nv1))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] > *thresh)){
      
      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) - *alpha * 
	R_pow_di(z1[i] + z2[i], -2);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = -2 * *alpha * R_pow_di(z1[i] + z2[i], -3);
      
      if (!R_FINITE(log(nv1 * nv2 - v12))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      
    }
  }
      
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];

  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = - 2/lambda2 + *alpha / (2 * lambda2);      
    *dns = *dns + (*nj - *nnj) * censCont;
  }

   //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}


void gpdmcamixtest(double *data1, double *data2, double *data3, int *nj,
		   int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
		   double *scale, double *shape, double *alpha, double *asCoef,
		   double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, c1, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
    
  if(*alpha < 0  || *scale < 0.01 || *alpha + 2 * *asCoef > 1 || 
     *alpha + 3 * *asCoef < 0 || *alpha + *asCoef > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);
    
    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    //c1 is a common value for cases 1 - 4
    c1 = ((*alpha + *asCoef) * z1[i] + (*alpha + 2 * *asCoef) *
	  z2[i]) / R_pow_di(z1[i] + z2[i], 2);
    v = 1/z1[i] + 1/z2[i] - c1;
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + 2 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      if (!R_FINITE(log(nv2))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv2) + nK2 - v;
      
    }
  
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      if (!R_FINITE(log(nv1))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv1) + nK1 - v;
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] > *thresh)){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + 2 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (4 * *alpha + 6 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 3) - 6 * c1 /
	R_pow_di(z1[i] + z2[i], 2);
      
      if (!R_FINITE(log(nv1 * nv2 - v12))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
    }
  }
  
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];
  
  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = ((*alpha + *asCoef) * lambda2 +
		(*alpha + 2 * *asCoef) * lambda2) /
      R_pow_di( 2 * lambda2, 2);
    censCont = censCont - 2 / lambda2;
    *dns = *dns + (*nj - *nnj) * censCont;
  }

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
}


void gpdmcamix(double *data1, double *data2, double *data3, int *nj,
	       int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	       double *scale, double *shape, double *alpha, double *asCoef,
	       double *dns){

  int i;

  double eps, *t1, *t2, *z1, *z2, *dvecj, v, nv1, nK1,
    nv2, nK2, v12, censCont, *dvecm, c1, lambda2;

  eps = R_pow(DOUBLE_EPS, 0.3);
  t1 = (double *)R_alloc(*nnj, sizeof(double));
  t2 = (double *)R_alloc(*nnj, sizeof(double));
  z1 = (double *)R_alloc(*nnj, sizeof(double));
  z2 = (double *)R_alloc(*nnj, sizeof(double));
  dvecj = (double *)R_alloc(*nnj, sizeof(double));
  dvecm = (double *)R_alloc(*nnm, sizeof(double));
    
  if(*alpha < 0  || *scale < 0.01 || *alpha + 2 * *asCoef > 1 || 
     *alpha + 3 * *asCoef < 0 || *alpha + *asCoef > 1){
    *dns = -1e6;
    return;
  }
  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //Preliminary data transformation stage//
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){

    //Margin 1
    t1[i] = (data1[i]  - *thresh) / *scale;
        
    if (data1[i] <= *thresh)
      t1[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t1[i] = exp(-t1[i]);
      }

      else {
	t1[i] = 1 + *shape * t1[i];
      
	if (t1[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t1[i] = R_pow(t1[i], -1/ *shape);
      
      }
    }

    //Margin 2
    t2[i] = (data2[i]  - *thresh) / *scale;

    if (data2[i] <= *thresh)
      t2[i] = 1;

    else{

      if (fabs(*shape) <= eps){
	*shape = 0;
	t2[i] = exp(-t2[i]);
      }

      else {
	t2[i] = 1 + *shape * t2[i];
      
	if (t2[i] <= 0){
	  *dns = -1e6;
	  return;
	}
      
	t2[i] = R_pow(t2[i], -1/ *shape);
      
      }
    }
       
    //Transform observed datas to unit frechet ones
    z1[i] = -1 / log(1 - *lambda * t1[i]);
    z2[i] = -1 / log(1 - *lambda * t2[i]);
    
    //To avoid numerical troubles
    if (!R_FINITE(z1[i]) || !R_FINITE(z2[i]) || !R_FINITE(log(z1[i])) ||
	!R_FINITE(log(z2[i]))){
      *dns = -1e6;
      return;
    }
  }

  
  //+++++++++++++++++++++++++++++++++++++//
  //                                     //
  //     Log-likelihood contribution     //
  //                                     //
  //+++++++++++++++++++++++++++++++++++++//

  for (i=0;i<*nnj;i++){ 

    //Compute the bivariate logistic distribution at point
    //(z[1], z[2])
    //c1 is a common value for cases 1 - 4
    c1 = ((2 * *alpha + *asCoef) * z1[i] + (*alpha + *asCoef) *
	  z2[i]) / R_pow_di(z1[i] + z2[i], 2);
    v = 1/z1[i] + 1/z2[i] - c1;
      
    if ((data1[i] <= *thresh) && (data2[i] > *thresh)){

      //Case 2: x1 <= threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      if (!R_FINITE(log(nv2))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = log(nv2) + nK2 - v;
      
      if (!R_FINITE(dvecj[i]))
	printf("Case 2: log(nv2) = %f, nK2 = %f, v = %f\n",
	       log(nv2), nK2, v);
      
    }
  
    if ((data1[i] > *thresh) && (data2[i] <= *thresh)){

      //Case 3: x1 > threshold1 & x2 <= threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (2 * *alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      if (!R_FINITE(log(abs(nv1)))){
	//*dns = -1e6;
	dvecj[i] = -1000 + nK1 - v;
	//nv1 = 1e-320;
	//return;
      }
      
      else
	dvecj[i] = log(abs(nv1)) + nK1 - v;

      if (!R_FINITE(dvecj[i]))
	printf("Case 3: log(nv1) = %f, nK1 = %f, v = %f\n",
	       log(nv1), nK1, v);
      
    }
    
    if ((data1[i] > *thresh) && (data2[i] > *thresh)){

      //Case 4: x1 > threshold1 & x2 > threshold2
      
      //Compute the negative partial derivative with
      //respect to the first component
      nv1 = R_pow_di(z1[i], -2) + (2 * *alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K1 constant of Ledford [1996]
      nK1 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t1[i]) + 2 * log(z1[i]) +
	1 / z1[i];
      
      //Compute the negative partial derivative with
      //respect to the second component
      nv2 = R_pow_di(z2[i], -2) + (*alpha + *asCoef) / 
	R_pow_di(z1[i] + z2[i], 2) - 2 * c1 /
	(z1[i] + z2[i]);
      
      //the log negative K2 constant of Ledford [1996]
      nK2 = log(*lambda) - log(*scale) + 
	(1 + *shape) * log(t2[i]) + 2 * log(z2[i]) +
	1 / z2[i];
      
      //Compute the partial mixed derivative	
      v12 = (4 * *alpha + 6 * *asCoef) / 
	R_pow_di(z1[i] + z2[i], 3) - 6 * c1 /
	R_pow_di(z1[i] + z2[i], 2);
      
      if (!R_FINITE(log(nv1 * nv2 - v12))){
	*dns = -1e6;
	return;
      }

      dvecj[i] = nK1 + nK2 + log(nv1 * nv2 - v12)
	- v;
      if (!R_FINITE(dvecj[i]))
	printf("Case 4: log(nv1 * nv2 - v12) = %f, nK1 = %f, nK2 = %f, v = %f\n",
	       log(nv1 * nv2 - v12), nK1, nK2, v);
    }
  }
  
  for (i=0;i<*nnj;i++)
    *dns = *dns + dvecj[i];
  
  printf("Joint part before censoring: dns = %f\n", *dns);
  
  //Now add the censored contribution to loglikelihood
  if (*nnj != *nj){
    lambda2 = - 1 / log(1 - *lambda);
    censCont = ((*alpha + *asCoef) * lambda2 +
		(*alpha + 2 * *asCoef) * lambda2) /
      R_pow_di( 2 * lambda2, 2);
    censCont = censCont - 2 / lambda2;
    *dns = *dns + (*nj - *nnj) * censCont;
  }
  printf("Joint part: dns = %f\n", *dns);

  //The marginal part:
  for(i=0;i<*nnm;i++)  {
    data3[i] = (data3[i] - *thresh) / *scale;
    
    if(fabs(*shape) <= eps) 
      dvecm[i] = log(*lambda) - log(*scale) - data3[i];
    else {
      data3[i] = 1 + *shape * data3[i];
      if(data3[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvecm[i] = log(*lambda) -log(*scale) - (1 / *shape + 1) *
	log(data3[i]);
    }
  }
  
  for(i=0;i<*nnm;i++)
    *dns = *dns - dvecm[i];
  
  printf("Marginal part before censoring: dns = %f\n", *dns);

  //Now add the censored contribution to loglikelihood
  if (*nm != *nnm)
    *dns = *dns - (*nm - *nnm) * log(1 - *lambda);
  printf("Marginal part: dns = %f\n", *dns);
}
