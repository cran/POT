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
