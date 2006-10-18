#include "header.h"

/* Conditional copulas condition on 2nd margin. */

double ccbvlog(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,idep,u,v,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  idep = 1/dep;
  u = R_pow(tm1, idep) + R_pow(tm2, idep);
  v = R_pow(u, dep);
  fval = exp(-v) * (1 / m2) * R_pow(tm2, idep-1) * R_pow(u, dep-1) - oldm1;
  return fval;
}

double ccbvalog(double m1, double m2, double oldm1, double dep, double asy1, 
                double asy2)
{
  double tm1,tm2,idep,u,v,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  idep = 1/dep;
  u = R_pow(asy1*tm1, idep) + R_pow(asy2*tm2, idep);
  v = (1-asy1)*tm1 + (1-asy2)*tm2 + R_pow(u, dep);
  fval = exp(-v) * (1 / m2) * (1 - asy2 + R_pow(asy2, idep) * 
    R_pow(tm2, idep-1) * R_pow(u, dep-1)) - oldm1;
  return fval;
}


double ccbvnlog(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,v,idep,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = R_pow((R_pow(tm2,-dep) + R_pow(tm1,-dep)),-idep);
  fval = exp(v) * m1 * (1-R_pow(1 + R_pow(tm2/tm1,dep), -1-idep)) - oldm1;
  return fval;
}

double ccbvanlog(double m1, double m2, double oldm1, double dep, 
		 double asy1, double asy2)
{
  double tm1,tm2,v,idep,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = R_pow(asy1 * tm2, -dep) + R_pow(asy2 * tm1, -dep);
  fval = exp(R_pow(v, -idep)) * m1 * (1 - R_pow(asy1, -dep) * 
    R_pow(tm2, -dep-1) * R_pow(v, -idep-1)) - oldm1;
  return fval;
}

double ccbvamix(double m1, double m2, double oldm1, double alpha, double beta)
{
  double tm1,tm2,tm1a,v,v2,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  tm1a = tm1/(tm1 + tm2);

  v = tm1 + tm2 - tm1 * ((alpha + beta) - alpha * tm1a - beta * tm1a * tm1a);
  v2 = 1 - alpha * tm1a * tm1a - 2 * beta * tm1a * tm1a * tm1a; 
  fval = exp(-v) * (1 / m2) * v2 - oldm1;

  return fval;
}

/*
   Calculates conditional copula for any model, conditioning on
   the margin `cnd'.
*/

void ccop(double *m1, double *m2, int *cnd, double *dep, double *asy1, double *asy2, double *alpha, double *beta, int *n, int *model, double *ccop)
{
  int i;

  switch(*model) {
  case 1:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvlog(m1[i], m2[i], 0, *dep);
      else  ccop[i] = ccbvlog(m2[i], m1[i], 0, *dep);
    }
    break;
  case 2:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvalog(m1[i], m2[i], 0, *dep, *asy1, *asy2);
      else ccop[i] = ccbvalog(m2[i], m1[i], 0, *dep, *asy2, *asy1); 
    }
    break;
  case 3:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvnlog(m1[i], m2[i], 0, *dep);
      else ccop[i] = ccbvnlog(m2[i], m1[i], 0, *dep);
    }
    break;  
  case 4:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvanlog(m1[i], m2[i], 0, *dep, *asy1, *asy2);
      else ccop[i] = ccbvanlog(m2[i], m1[i], 0, *dep, *asy2, *asy1);
    }
    break;
  case 5:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvamix(m1[i], m2[i], 0, *alpha, 0);
      else ccop[i] = ccbvamix(m2[i], m1[i], 0, *alpha + 3 * *beta, - *beta);
    }
    break;
  case 6:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvamix(m1[i], m2[i], 0, *alpha, *beta);
      else ccop[i] = ccbvamix(m2[i], m1[i], 0, *alpha + 3 * *beta, - *beta);
    }
    break;
  default:
     error("no copula found for this model");
  }
}
  



