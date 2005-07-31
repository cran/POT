#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()


void nlgpd(double *data, int *n, double *loc, double *scale, 
	   double *shape, double *dns);

