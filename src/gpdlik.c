#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

void gpdlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void gpdlik2(double *data, int *n, double *loc, double *scale,
	     double *shape, double *dns);


void gpdlik(double *data, int *n, double *loc, double *scale, 
    double *shape, double *dns)
{
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
     *dns = R_NegInf;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    if(fabs(*shape) <= eps) 
      dvec[i] = -log(*scale) - data[i];
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
        *dns = R_NegInf;
        return;
      }
      dvec[i] = -log(*scale) - (1 / *shape + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];
}

void gpdlik2(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns)
{
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
     *dns = R_NegInf;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[i]) / *scale;
    if (data[i] <= 0) {
      *dns = R_NegInf;
      return;
    }
    if(fabs(*shape) <= eps) 
      dvec[i] = -log(*scale) - data[i];
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
	*dns = R_NegInf;
	return;
      }
      dvec[i] = -log(*scale) - (1 / *shape + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];
}
