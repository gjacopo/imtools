/* 
 * $Log: noise_estimation.c,v $
 * Revision 1.3  1995/11/09  13:02:33  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:20:28  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:20:28  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:56:21  ohenri
 * Initial revision
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include "tools.h"
 
#define NN 5

void NoiseEstimation(float *data,int nx,int ny,float c,int Upper,
		     float *mode,float *dmax,float *dmin,float *thresh)
{
  int i,j,N,J,JJ[NN];
  float *srt=NULL;
  float maxpx[NN],px[NN],m[NN],mode_avg,signif,sum,norm;
  
  N=nx*ny;
  J=N*0.50;
  
  srt=(float *) malloc(N*sizeof(float));
  memcpy((void *) srt,(void *) data,N*sizeof(float));
  sort(N,srt-1);
  
  for (j=0;j<NN;j++) {
    maxpx[j]=-1e38;
    /*    JJ[j]=(int) (J - j*J/(float) NN); */
   JJ[j]=(0.25 + 0.05*j)* (float) N;
  }

  for (i=0;i<N-J;i++){
    if (srt[i] > 0.0){ /* there might be a peak at 0 and we don't want it ! */
      for (j=0;j<NN;j++){
	px[j]= 1.0 / (srt[i+JJ[j]] - srt[i]);
	if (px[j] > maxpx[j]) {maxpx[j]=px[j];m[j]=0.5*( srt[i+JJ[j]] + srt[i]); }
      }
    }
  }

  for (j=0,mode_avg=sum=norm=0.0;j<NN;j++)  {
    mode_avg += m[j]; 
  }
  *mode = mode_avg / (float) NN; 
    
  signif=fsqrt(-2.0*log((100.0-c)/100.0));

  if (Upper)  *thresh=srt[N-1] - (srt[N-1] - *mode) * signif;
  else  *thresh= *mode * signif;
  
  *dmax=srt[N-1];
  *dmin=srt[0];

  free(srt);
}



