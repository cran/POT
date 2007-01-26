#include "header.h"

/* produces standard Frechet margins */
void rbvlog_shi(int *n, double *alpha, double *sim)
{
  double u,z;
  int i;
  
  RANDIN;
  for(i=0;i<*n;i++) 
  { 
    u = UNIF;
    if(UNIF < *alpha) z = EXP+EXP;
    else z = EXP;
    sim[2*i] = 1/(z * R_pow(u,*alpha));
    sim[2*i+1] = 1/(z * R_pow(1-u,*alpha));
  }
  RANDOUT;
}

/* produces standard Frechet margins */
void rbvalog_shi(int *n, double *alpha, double *asy, double *sim)
{
  double v1_1,v2_2,v1_12,v2_12,u,z;
  int i;
    
  RANDIN;

  if(*alpha == 1)
    for(i=0;i<2*(*n);i++) sim[i] = 1/EXP;
  else {
    for(i=0;i<*n;i++) 
    {
      v1_1 = (1-asy[0]) / EXP;
      v2_2 = (1-asy[1]) / EXP;
      u = UNIF;
      if(UNIF < *alpha) z = EXP+EXP;
      else z = EXP;
      v1_12 = asy[0] / (z * R_pow(u,*alpha));
      v2_12 = asy[1] / (z * R_pow(1-u,*alpha));
      sim[2*i] = fmax2(v1_1,v1_12); 
      sim[2*i+1] = fmax2(v2_2,v2_12);
    }
  }
  RANDOUT;
}

/* produces uniform margins; needed for simmc */
void rbvlog(int *n, double *dep, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvlog(llim, sim[2*i+1], sim[2*i+0], *dep);
    uval = ccbvlog(ulim, sim[2*i+1], sim[2*i+0], *dep);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvlog(midpt, sim[2*i+1], sim[2*i+0], *dep);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins; needed for simmc */
void rbvalog(int *n, double *dep, double *asy, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvalog(llim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    uval = ccbvalog(ulim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvalog(midpt, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvnlog(int *n, double *dep, double *sim){

  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvnlog(llim, sim[2*i+1], sim[2*i+0], *dep);
    uval = ccbvnlog(ulim, sim[2*i+1], sim[2*i+0], *dep);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvnlog(midpt, sim[2*i+1], sim[2*i+0], *dep);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvanlog(int *n, double *dep, double *asy, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvanlog(llim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    uval = ccbvanlog(ulim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvanlog(midpt, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}



/* produces uniform margins */
void rbvmix(int *n, double *alpha, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvmix(llim, sim[2*i+1], sim[2*i+0], *alpha);
    uval = ccbvmix(ulim, sim[2*i+1], sim[2*i+0], *alpha);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvmix(midpt, sim[2*i+1], sim[2*i+0], *alpha);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}
/* produces uniform margins */
void rbvamix(int *n, double *alpha, double *beta, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvamix(llim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    uval = ccbvamix(ulim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvamix(midpt, sim[2*i+1], sim[2*i+0], *alpha, *beta);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
    }
    sim[2*i+0] = midpt;
  }
}
