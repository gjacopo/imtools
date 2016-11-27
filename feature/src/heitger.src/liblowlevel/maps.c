/* 
 * $Log: maps.c,v $
 * Revision 1.2  1995/10/02  12:25:36  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:56:21  ohenri
 * Initial revision
 *
 * Revision 1.1  1995/06/21  12:56:21  ohenri
 * Initial revision
 * 
 *
 */

/*===========================================================================
 *
 * Institute of Communications Technology, Image Vision Group
 * Swiss Federal Institute of Technology at Zurich, Switzerland.
 *
 *==========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include "lowlevel.h"

void quality_dft(float *modulus,float *angle,int i_size,int nori,
		 float cosexp,float *quality)
{
  int *nz,i,k,*ang_int;
  float *ori_list,omega_max,mod_max,fang,tmp,*omega;
  float ang_min;
 
  /* initialize data */

  nz=(int *) malloc(nori*sizeof(int));
  ori_list=(float *) malloc(nori*sizeof(float));
  omega=(float *) malloc(91*sizeof(float));
  ang_int=(int *) malloc(nori*sizeof(int));
  
  for (k=0;k<nori;k++){
    nz[k] = i_size*k;
    ori_list[k] = (float) k * (180.0 / (float) nori);  /* degrees */
  }

  for (k=0;k<91;k++){
    tmp = (float) k * (M_PI / 180.0);  /* angular map in degrees !!! */
    omega[k]=(float) (pow (fabs (cos (tmp)),cosexp));
  }

  for (i=0;i<i_size;i++){
    
    /* calculate quality */

    mod_max=-1e38;ang_min=1e38;
    for (k=0;k<nori;k++){
      fang=ABS(ori_list[k]-angle[i]*(180.0/M_PI));
      if (fang > 90.0) fang = 180.0 - fang;
      if (fang < ang_min) ang_min = fang;
      ang_int[k] = (int) fang;
      if (modulus[i+nz[k]] > mod_max) mod_max = modulus[i+nz[k]];
    }
    
    if (mod_max > 0.0){
      tmp=1.0; 
      omega_max=omega[((int) ang_min)];
      for (k=0;k<nori;k++) 
	tmp +=SQR(modulus[i+nz[k]]/mod_max-omega[ang_int[k]]/omega_max);
      quality[i]=1.0/tmp;   /* a high value indicating good quality */
    }
    else quality[i]=0;
  }
  
}

void amplitude_phase(float *modulus,int i_size,int nori,
		     float *phase,float *amplitude,int deg)
{
  int *nz,i,k;
  float real=0.0,imag=0.0,dft_tmp,*coscos, *sinsin, *ori_list;
  float omega_max,mod_max,fang,tmp,*omega,max_q=-1e38,min_q=1e38;
  float max_mod,avg;

  /* initialize data */

  nz=(int *) malloc(nori*sizeof(int));
  coscos=(float *) malloc(nori*sizeof(float));
  sinsin=(float *) malloc(nori*sizeof(float));
  
  dft_tmp=2.0*M_PI/(float) nori;
       
  for (k=0;k<nori;k++){
    coscos[k] = fcos((float) k * dft_tmp);
    sinsin[k] = fsin((float) k * dft_tmp);
    nz[k] = i_size*k;
  }
  
  for (i=0;i<i_size;i++){
    /* calculate  dominant orientation */
    
    real=imag=max_mod=avg=0.0;
    for (k=0;k<nori;k++){
      real += modulus[i+nz[k]]*coscos[k];
      imag += modulus[i+nz[k]]*sinsin[k];
      avg  += modulus[i+nz[k]];
      if (modulus[i+nz[k]] > max_mod) max_mod=modulus[i+nz[k]];
    }
    
    if ((max_mod > 0.0)&& phase)
      if (deg) phase[i] = 90.0*(1+atan2(-imag,-real)/M_PI); 
      else phase[i] = 0.5*(M_PI+atan2(-imag,-real)); 
    
    if (amplitude) amplitude[i]=avg+fsqrt(SQR(imag)+SQR(real));
  }
}
