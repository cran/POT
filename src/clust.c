# include "header.h"

void clust(int *n, double *obs, double *tim, double *cond,
	   double *thresh, double *clust){

  int i, idx, idxMax;

  i = 0;

  while (i < *n){

    if (obs[i] > *thresh){
      clust[i] = i+1;
      idx = i + 1;
      idxMax = idx;

      while( (obs[idx] > *thresh) ||
	     ( (tim[idx]-tim[idxMax]) <= *cond)){
	if (obs[idx] >= obs[idxMax])
	  idxMax = idx;
	idx++;

	if (idx >= *n){
	  idx--;
	  break;
	}
      }
      
      clust[i+1] = idx;
      i = idx + 1;
    }

    else
      i++;    
  }
}
