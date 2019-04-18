/*****************************************************************************
 *   Copyright (c) 2018 Mathieu Ribatet, Christophe Dutang                                                                                                  
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

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()



//From POT.c
void do_gpdlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void do_pplik(double *data, int *n, double *loc, double *scale,
	   double *shape, double *thresh, double *noy, double *dns);
void do_samlmu(double *x, int *nmom, int *n, double *lmom);

//From clust.c
void do_clust(int *n, double *obs, double *tim, double *cond,
	   double *thresh, double *clust);

//From ts2tsd.c
void do_ts2tsd(double *time, double *obs, double *start,
	    double *end, double *obsIntStart, double *obsIntEnd,
	    int *n, double *ans); 

//From cgpdlik.c
void do_cgpdlik(double *data, int *n, int *nn, double *lambda, double *loc,
	     double *scale, double *shape, double *dns);

//From bvgpdlik.c
void do_gpdbvlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	      double *lambda2, double *thresh, double *scale1,
	      double *shape1, double *scale2, double *shape2,
	      double *alpha, double *dns);
void do_gpdbvalog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *asCoef1, double *asCoef2,
	       double *dns);
void do_gpdbvnlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *dns);
void do_gpdbvanlog(double *data1, double *data2, int *n, int *nn, double *lambda1,
		double *lambda2, double *thresh, double *scale1,
		double *shape1, double *scale2, double *shape2,
		double *alpha, double *asCoef1, double *asCoef2,
		double *dns);
void do_gpdbvmix(double *data1, double *data2, int *n, int *nn, double *lambda1,
	      double *lambda2, double *thresh, double *scale1,
	      double *shape1, double *scale2, double *shape2,
	      double *alpha, double *dns);
void do_gpdbvamix(double *data1, double *data2, int *n, int *nn, double *lambda1,
	       double *lambda2, double *thresh, double *scale1,
	       double *shape1, double *scale2, double *shape2,
	       double *alpha, double *asCoef, double *dns);
/*never defined by M. Ribatet in C code
 * 
 void gpdbvamixtest(double *data1, double *data2, int *n, int *nn, double *lambda1,
		   double *lambda2, double *thresh, double *scale1,
		   double *shape1, double *scale2, double *shape2,
		   double *alpha, double *asCoef, double *dns);*/

//From simmc.c
void do_rbvlog_shi(int *n, double *alpha, double *sim);
void do_rbvalog_shi(int *n, double *alpha, double *asy, double *sim);
void do_rbvlog(int *n, double *dep, double *sim);
void do_rbvalog(int *n, double *dep, double *asy, double *sim);
void do_rbvnlog(int *n, double *dep, double *sim);
void do_rbvanlog(int *n, double *dep, double *asy, double *sim);
void do_rbvmix(int *n, double *alpha, double *sim);
void do_rbvamix(int *n, double *alpha, double *beta, double *sim);

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
double ccbvmix(double m1, double m2, double oldm1, double alpha);
double ccbvamix(double m1, double m2, double oldm1, double alpha,
		double beta);

//From mcgpdlik.c
void do_gpdmclog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	      double *scale, double *shape, double *alpha, double *dns);
void do_gpdmcalog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	       double *scale, double *shape, double *alpha, double *asCoef1,
	       double *asCoef2, double *dns);
void do_gpdmcnlog(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	       double *scale, double *shape, double *alpha, double *dns);
void do_gpdmcanlog(double *data1, double *data2, double *data3, int *nj,
		int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
		double *scale, double *shape, double *alpha, double *asCoef1,
		double *asCoef2, double *dns);
void do_gpdmcmix(double *data1, double *data2, double *data3, int *nj,
	      int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	      double *scale, double *shape, double *alpha, double *dns);
void do_gpdmcamix(double *data1, double *data2, double *data3, int *nj,
	       int *nnj, int *nm, int *nnm, double *lambda, double *thresh,
	       double *scale, double *shape, double *alpha, double *asCoef,
	       double *dns);
