#include "header.h"

void nlgpd(double *data, int *n, double *loc, double *scale, 
    double *shape, double *dns)
{
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    if(fabs(*shape) <= eps) 
      dvec[i] = log(1 / *scale) - data[i];
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
        *dns = 1e6;
        return;
      }
      dvec[i] = log(1 / *scale) - (1 / *shape + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}
