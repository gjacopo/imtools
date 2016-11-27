/* 
 * $Log: kpt.c,v $
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

/***********************************************************************
 *   Program Development: 
 *   
 *   Federal Institute of Technology, ETH-Zurich.
 *   Communication Technology Lab, Image Science
 *   Gloriastr. 35
 *   ETH-Zentrum
 *   CH-8092 Zurich
 *
 *   email: heitger@vision.ethz.ch
 *      
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include "lowlevel.h"

/* default:  x_ridge_gain=2.5, gamma=0.18, comp_gain=1.0, even_gain=2.0  */

int keypoint_map(float *modulus,int nori,int nx,int ny,
		 float x_ridge_gain,float gamma,float comp_gain,
		 float even_gain,float *kp_data)
{
  int *pder1_x,*pder1_y,*pder2_x,*pder2_y,*oder1_x,*oder1_y,*oder2_x,*oder2_y;
  float *pder1_val,*pder2_val,*oder1_val,*oder2_val,ang_inc,angle;
  float co,si,*kp_val,rtmp,max_kp_val;
  int i,tot2_ori,x,y,pos,max_kp_ori,*oriwrap;
  float cross_ridge,ortho_sum,oder_para,oder_ortho;
  float oder1_sum, oder2_sum,oder1_gain=0.5,oder2_gain=0.5,sms=1.0,comp, esum;
 
  /*========================================================================
   *   keypoint detection scheme ...
   *========================================================================
   */

  tot2_ori = nori / 2;
 
  pder1_x=(int *) malloc(nori*sizeof(int));
  pder1_y=(int *) malloc(nori*sizeof(int));
  pder2_x=(int *) malloc(nori*sizeof(int));
  pder2_y=(int *) malloc(nori*sizeof(int));
  oder1_x=(int *) malloc(nori*sizeof(int));
  oder1_y=(int *) malloc(nori*sizeof(int));
  oder2_x=(int *) malloc(nori*sizeof(int));
  oder2_y=(int *) malloc(nori*sizeof(int));
  oriwrap=(int *) malloc(3*nori*sizeof(int));

  ang_inc = (float) M_PI / (float) nori;
  for (i = 0; i < nori; i++) {
    angle = i*ang_inc;

    oriwrap[i] = oriwrap[i + nori] = oriwrap[i + 2*nori] = i;
    co=2.0*fcos (angle);
    si=2.0*fsin (angle);
    
    pder1_x[i] = pder2_x[i] =  IROUND (co);
    pder1_y[i] = pder2_y[i] =  -IROUND (si)*nx;
    oder1_y[i] = oder2_y[i] =  IROUND (co)*nx;
    oder1_x[i] = oder2_x[i] =  IROUND (si);
  }
    
  /*=========================================================================
   *     make scalar keypoint map
   *=========================================================================
   */
    
  pder1_val=(float *) malloc(nori*sizeof(float));
  pder2_val=(float *) malloc(nori*sizeof(float));
  oder1_val=(float *) malloc(nori*sizeof(float));
  oder2_val=(float *) malloc(nori*sizeof(float));
  kp_val   =(float *) malloc(nori*sizeof(float));
  memset ((void *) kp_data, 0, nx*ny*sizeof (float));
  
  for (y = 2; y < (ny - 2); y++) {
    for (x = 2; x < (nx - 2); x++) {
      oder1_sum = oder2_sum = 0.0;
      
      for (i = 0; i < nori; i++) {
	pos=i*nx*ny+y*nx+x;   /* center position to be evaluated */
	
	/*...................................................................
	 * first we do the 1. parallel derivative
	 */
	
	rtmp=(float) modulus[pos-pder1_x[i]-pder1_y[i]]-
	  modulus[pos+pder1_x[i]+pder1_y[i]];
	pder1_val[i] = ABS (rtmp);
	
	/*.....................................................................
	 * now we do the 2. parallel derivative
	 */
	
	rtmp=(float)  modulus[pos]-(modulus[pos+pder2_x[i]+pder2_y[i]] +
				    modulus[pos-pder2_x[i]-pder2_y[i]])*0.5;
	pder2_val[i] = CLIP (rtmp)*even_gain;
	
	/*.....................................................................
	 * now we do the orthogonal derivatives
	 */
	
	rtmp = modulus[pos+oder1_x[i]+oder1_y[i]] -
	  modulus[pos-oder1_x[i]-oder1_y[i]];
	oder1_val[i] = ABS (rtmp);
	
	rtmp =  modulus[pos]- (modulus[pos+oder2_x[i]+oder2_y[i]]+
			       modulus[pos-oder2_x[i]-oder2_y[i]])*0.5;
	oder2_val[i] = CLIP (rtmp);
      }
      
      /*.....................................................................
       * now we do the compensation map
       */
      
      cross_ridge = 0.0;
      ortho_sum = 0.0;
      for (i = 0; i < tot2_ori; i++) {
	oder_para  = oder2_val[i];
	oder_ortho = oder2_val[oriwrap[i + tot2_ori]];
	cross_ridge += fsqrt (oder_para*oder_ortho);  
	ortho_sum += (oder_para + oder_ortho);        
      }
      rtmp = cross_ridge - ortho_sum*gamma;    
      cross_ridge = CLIP (rtmp)*x_ridge_gain;         
      oder1_sum = oder2_sum = 0.0;   
      for (i = 0; i < nori; i++) {
	rtmp = (oder1_val[i] - cross_ridge);          
	oder1_sum += CLIP (rtmp);
	
	rtmp = (oder2_val[i] - cross_ridge);
	oder2_sum += CLIP (rtmp);
      }
      comp = oder1_sum*oder1_gain + oder2_sum*oder2_gain;
      
      max_kp_val = 0.0;
      for (i = 0; i < nori; i++){ 
	kp_val[i]=fsqrt(SQR(pder1_val[i])+SQR(pder2_val[i]))-comp*comp_gain;
	kp_val[i] = CLIP (kp_val[i]);
	if ( kp_val[i] > max_kp_val) max_kp_val = kp_val[i];
      }
      kp_data[y*nx+x] = max_kp_val;
    }
  }
  
  free(pder1_x);free(pder1_y);free(pder2_x);free(pder2_y);
  free(oder1_x);free(oder1_y);free(oder2_x);free(oder2_y);
  free(oriwrap);free(pder1_val);
  free(pder2_val);free(oder1_val);free(oder2_val);free(kp_val);
  

  
  if (sms > 0.0) {
    return (int) smoothing (kp_data, nx, ny, sms);
  }
  else return 2;
  
}

