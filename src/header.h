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
void clust(int *n, double *obs, double *tim, double *cond,
	   double *thresh, double *clust);

//From ts2tsd.c
void ts2tsd(double *time, double *obs, double *start,
	    double *end, double *obsIntStart, double *obsIntEnd,
	    int *n, double *ans); 

//From cgpdlik.c
void cgpdlik(double *data, int *n, int *nn, double *lambda, double *loc,
	     double *scale, double *shape, double *dns);

//From bvgpdlik.c
void gpdbvlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	      double *lambda2, double *thresh, double *scale1,
	      double *shape1, double *scale2, double *shape2,
	      double *alpha, double *dns);
void gpdbvalog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *asCoef1, double *asCoef2,
	       double *dns);
void gpdbvnlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *dns);
void gpdbvanlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
		double *lambda2, double *thresh, double *scale1,
		double *shape1, double *scale2, double *shape2,
		double *alpha, double *asCoef1, double *asCoef2,
		double *dns);
void gpdbvmix(double *data1, double *data2, int *n, int *nn, double *lambda1,
	      double *lambda2, double *thresh, double *scale1,
	      double *shape1, double *scale2, double *shape2,
	      double *alpha, double *dns);
void gpdbvamix(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *asCoef, double *dns);

//From simmc.c
void rbvlog(int *n, double *dep, double *sim);
void rbvalog(int *n, double *dep, double *asy, double *sim);
void rbvnlog(int *n, double *dep, double *sim);
void rbvanlog(int *n, double *dep, double *asy, double *sim);
void rbvamix(int *n, double *alpha, double *beta, double *sim);

//From ccop.c
void ccop(double *m1, double *m2, int *cnd, double *dep,
	  double *asy1, double *asy2, double *alpha,
	  double *beta, int *n, int *model, double *ccop);
double ccbvlog(double m1, double m2, double oldm1, double dep);
double ccbvalog(double m1, double m2, double oldm1, double dep,
		double asy1, double asy2);
double ccbvnlog(double m1, double m2, double oldm1, double dep);
double ccbvanlog(double m1, double m2, double oldm1, double dep,
		 double asy1, double asy2);
double ccbvamix(double m1, double m2, double oldm1, double alpha,
		double beta);
