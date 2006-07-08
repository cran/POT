#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

//From POT.c
void gpdlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void samlmu(double *x, int *nmom, int *n, double *lmom);

//From clust.c
void clust(int *n, double *idxExcess, double *diffTim, double *timCond,
	     double *clust);

//From ts2tsd.c
void ts2tsd(double *time, double *obs, double *start,
	    double *end, double *obsIntStart, double *obsIntEnd,
	    int *n, double *ans); 
