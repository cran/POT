# include "header.h"

void clust(int *n, double *idxExcess, double *diffTim, double *timCond,
	   double *clust){
  int i,j;
  
  j = 1;
  clust[0] = idxExcess[0];
  
  for (i=0;i<(*n-1);i++){
    if(diffTim[i] > *timCond){
      clust[j] = idxExcess[i];
      clust[j+1] = idxExcess[i+1];
      j = j + 2;
    }
  }
  
  clust[j] = idxExcess[*n-1];
}


