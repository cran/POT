#include "header.h"

void ts2tsd(double *time, double *obs, double *start,
	    double *end, double *obsIntStart, double *obsIntEnd,
	    int *n, double *ans){
  //time       : a vector for the time occurence
  //obs        : a vector for the observation
  //start      : a vector for the starting time
  //             of mobile average windows 
  //end        : the same but for the ending time
  //obsIntStart: Interpolated observation
  //             related to the starting 
  //             time ``start''
  //obsIntEnd  : the same but related to the
  //             ending time ``end''
  //n          : the number of observations
  //ans        : the response of the code

  int i, j, indexStart, indexEnd;
  
  //indexStart : index of the first observation
  //             in ``obs'' within the average
  //             window
  //indexEnd   : index of the last observation
  //             in ``obs'' within the average
  //             window


  for (i=0;i<*n;i++){
    
    indexStart = i;
    indexEnd = i;

    while ( (time[indexStart] > start[i]) && (indexStart >= 0) )
      indexStart--;
    indexStart++;

        
    while ( (time[indexEnd] <= end[i]) && (indexEnd < *n) )
      indexEnd++;
    indexEnd--;

    if ( (obsIntStart[i] == -1e6) || (obsIntEnd[i] == -1e6) ||
	 (obs[indexStart] == -1e6) || (obs[indexEnd] == -1e6))
      ans[i] = NA_REAL;

    else{
    
      for (j=indexStart;j<indexEnd;j++){
	//printf("obs[j] = %f\n", obs[j]);
	//printf("obs[j+1] = %f\n", obs[j+1]);
	//printf("time[j] = %f\n", time[j]);
	//printf("time[j+1] = %f\n", time[j+1]);
	

	if ( (obs[j] == -1e6) || (obs[j+1] == -1e6) )
	  ans[i] = NA_REAL;

	else{
	  ans[i] = ans[i] + (obs[j] + obs[j+1]) * 
	    (time[j+1] - time[j]) / 2.0;
	}
      }
      //Add the areas for interpolated points - i.e. boundaries
      //printf("obs[indexStart]=%f\n", obs[indexStart]);
      //printf("obs[indexEnd]=%f\n", obs[indexEnd]);
      //printf("obsIntStart[i]=%f\n", obsIntStart[i]);
      //printf("obsIntEnd[i]=%f\n", obsIntEnd[i]);
      
      ans[i] = ans[i] + (obs[indexStart] + obsIntStart[i]) *
	(time[indexStart] - start[i]) / 2.0;
      ans[i] = ans[i] + (obsIntEnd[i] + obs[indexEnd]) *
	(end[i] - time[indexEnd]) / 2.0;
    }
  }
}  
