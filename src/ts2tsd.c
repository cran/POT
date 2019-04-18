/*****************************************************************************
 *   Copyright (c) 2014 Mathieu Ribatet                                                                                                  
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

#include "header.h"

void do_ts2tsd(double *time, double *obs, double *start,
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
  
  //indexStart : index of the first observation
  //             in ``obs'' within the average
  //             window
  //indexEnd   : index of the last observation
  //             in ``obs'' within the average
  //             window
  
  
  for (int i=0;i<*n;i++){
    
    int indexStart = i,
      indexEnd = i;
    
    while (time[indexStart] > start[i]){
      indexStart--;
      
      if (indexStart < 0)
        break;
    }
    
    indexStart++;
    
    
    while (time[indexEnd] <= end[i]){
      indexEnd++;
      
      if (indexEnd >= *n)
        break;
    }
    
    indexEnd--;
    
    if ( (obsIntStart[i] == -1e6) || (obsIntEnd[i] == -1e6) ||
         (obs[indexStart] == -1e6) || (obs[indexEnd] == -1e6))
      ans[i] = NA_REAL;
    
    else{
      
      for (int j=indexStart;j<indexEnd;j++){
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
