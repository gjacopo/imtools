/* 
 * $Log: ll_processing.c,v $
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include "lowlevel.h"
#define RESOLUTION 512
#define LAMBDA_LOW 0.0
#define LAMBDA_UP 2.0

void Low_level_Processing(unsigned char *image,int nx,int ny,int nori,
			  float sigma,float sweep,float orisel,float supfact,
			  float enhfact,float cosweight,
			  float Xridge,float gamma,float comp_gain,
			  float even_gain,
			  float *maxmap,
			  float *mod_maxmap,
			  float *locori,
			  float *nms,
			  unsigned char *type,
			  float *kpt_map,
			  float *quality,
			  int Version)
{
  int x,y,xypos,k,f_no,support,i,ii,offset[24],maxori;

  /* new */
  float *fimage=NULL;
  int i_size=nx*ny;
  float lambda;
  float *even,*odd,*modulus;
  int border;

  /* contour detection */
  float fder,sder,tmp[3],ttmp[3];
  float sico_x[24],sico_y[24],cosa[24],sina[24],cosasq[24],sinasq[24],ang_incr;
  float coscos[24],sinsin[24];
  float realO,imagO,realE,imagE,maxvalO,maxvalE,amplO,amplE,par,sumO,sumE;
  float e,o,maxmod,a[6];
  float *ONE,*TWO;
 
  /* intern allocated maps */
  fimage=(float *) malloc(nx*ny*sizeof(float));
  even=(float *) malloc(nx*ny*nori*sizeof(float));
  odd=(float *) malloc(nx*ny*nori*sizeof(float));
  modulus=(float *) malloc(nx*ny*nori*sizeof(float));
  
  /* cleared maps */
  memset((void *) even,    0, nx*ny*nori*sizeof(float));
  memset((void *) odd,     0, nx*ny*nori*sizeof(float));
  memset((void *) modulus, 0, nx*ny*nori*sizeof(float));
  
  memset((void *) maxmap,  0, nx*ny*sizeof(float));
  if (mod_maxmap) memset((void *) mod_maxmap,  0, nx*ny*sizeof(float));
  memset((void *) locori,  0, nx*ny*sizeof(float));
  memset((void *) nms   ,  0, nx*ny*sizeof(float));
  memset((void *) type  ,  0, nx*ny*sizeof(char));

  /* convert char image to float */
  for (i=0;i<i_size;i++) fimage[i]=(float) image[i];

  /* call egabor ... */  

  lambda=egabor_params(sigma,sweep,LAMBDA_LOW,LAMBDA_UP,RESOLUTION,0); 
  egabor_filter(fimage,nx,ny,sigma,sweep,lambda,nori,orisel,even,odd,modulus);
  
  /* Defining the sampling distance for 1st and 2nd derivatives */
  ang_incr = M_PI / ((float) nori);
  for (i = 0; i < nori; i++) {
    sico_x[i] =  0.5*sigma*fcos ((float)i * ang_incr + M_PI/2.0);
    sico_y[i] = -0.5*sigma*fsin ((float)i * ang_incr + M_PI/2.0);
    cosa[i]   = (fcos ( (float)i * ang_incr + M_PI/2.));
    sina[i]   = -(fsin ( (float)i * ang_incr + M_PI/2.));
    cosasq[i] = SQR(cosa[i]);
    sinasq[i] = SQR(sina[i]);
  }
  support=(int) sigma + 1;
   
  /* Set weights for Discrete Fourier Transformation */
  for (k=0;k<nori;k++){
    coscos[k] = fcos((float) k * 2.0*M_PI/(float) nori);
    sinsin[k] = fsin((float) k * 2.0*M_PI/(float) nori);
  }
  
  fprintf(stderr,"Compute SE operator (version %d)\n",Version);
  fflush(stderr);
  for (y = support; y < (ny - support); y++) {
    for (x = support; x < (nx - support); x++) {
      xypos=x+y*nx;
      sumO=sumE=realO=imagO=maxvalO=realE=imagE=maxvalE=maxmod=0.0;
      for (k=maxori=0;k<nori;k++){
	/* Calculate first and second derivatives in each channel */

	/* old method 
	   tmp[1]=fbilin(modulus+k*nx*ny,(x+sico_x[k]),(y+sico_y[k]),nx);
	   tmp[2]=fbilin(modulus+k*nx*ny,(x-sico_x[k]),(y-sico_y[k]),nx);
	   sder = 2.0*modulus[k*nx*ny+xypos]-tmp[1]-tmp[2];
	   fder = tmp[1]-tmp[2];
	*/
	sofit(modulus+k*nx*ny,x,y,nx,a);

	fder = fabs(cosa[k]*a[1] + sina[k]*a[2]);
	sder = -2.0*(cosasq[k]*a[3] + sinasq[k]*a[4] +
				 sina[k]*cosa[k]*a[5]); 

	/* Calculate the modified sin/cosine maps */

	if (Version == 1) o = odd[k*nx*ny+xypos] = CLIP(fabs(odd[k*nx*ny+xypos]) -
							supfact*fabs(fder) + enhfact*sder);
	else if (Version == 2) o = odd[k*nx*ny+xypos] = CLIP(SQR(odd[k*nx*ny+xypos]) - 
							     supfact*SQR(fder) + 
							     enhfact*SQR(sder)*SGN(sder));
	sumO += o;
	if (o > maxvalO) maxvalO=o;
	if (Version == 1) e = even[k*nx*ny+xypos] = CLIP(fabs(even[k*nx*ny+xypos]) - 
							 supfact*fabs(fder) + enhfact*sder);
	else if (Version == 2) e=even[k*nx*ny+xypos]=CLIP(SQR(even[k*nx*ny+xypos]) - 
							  supfact*SQR(fder) + 
							  enhfact*SQR(sder)*SGN(sder));

	
	sumE += e;
	if (e > maxvalE) maxvalE=e;
	
	/* analyse first Fourier harmonic of channel distribution 
	   to estimate local orientation */

	realO += o*coscos[k];imagO += o*sinsin[k];
	realE += e*coscos[k];imagE += e*sinsin[k];
	if (modulus[k*nx*ny+xypos] > maxmod){
	  mod_maxmap[xypos]=maxmod=modulus[k*nx*ny+xypos];
	  maxori=k;
	}
	
      }
      
      locori[xypos] = ang_incr*maxori;
      if (sumO > sumE) type[xypos]=EDGETYPE;
      else type[xypos]=LINETYPE;

      if ((maxvalO != 0.0) || (maxvalE != 0.0)){
	if (Version == 1) amplO=fsqrt(SQR(imagO)+SQR(realO));
	else if (Version == 2) amplO=fsqrt(fsqrt(SQR(imagO)+SQR(realO)));
	amplE=cosweight*fsqrt(SQR(imagE)+SQR(realE));
	if (amplE > amplO){
	  if (imagE != 0.0 && realE !=0.0)
	    locori[xypos] = atan2(imagE,realE)*0.5;
	  if(locori[xypos] < 0.0) locori[xypos] += M_PI;
	  maxmap[xypos] = fabs(amplE);
	}
	else if (amplO > amplE){
	  if (imagO != 0.0 && realO !=0.0)
	    locori[xypos] = atan2(imagO,realO)*0.5;
	  if(locori[xypos] < 0.0) locori[xypos] += M_PI;
	  maxmap[xypos] = fabs(amplO);
	}
      }
    }
  }
  
  /* Non-maximum supression of odd and even maps */
  /* Sampling distance is smaller here ... */
  for (i = 0; i < nori; i++) {
    sico_x[i] = fcos ((float)i * ang_incr + M_PI/2.0);
    sico_y[i] = -fsin ((float)i * ang_incr + M_PI/2.0);
  }
  
  fprintf(stderr,"Compute non-maximum suppression\n");
  fflush(stderr);
  for (y = support; y < (ny - support); y++) {
    for (x = support; x < (nx - support); x++) {
      xypos=x+y*nx;
      i=angle2ori(locori[xypos],nori);
           
      if (type[xypos] == EDGETYPE) 
	{ONE=odd;TWO=even;par=1.0/cosweight;}
      else {ONE=even;TWO=odd;par=cosweight;}
      
      tmp[0] = fabs(ONE[xypos+nx*ny*i]);
      tmp[1]  = fabs(fbilin(ONE + nx*ny*i,(x+sico_x[i]),(y+sico_y[i]), nx));
      tmp[2]  = fabs(fbilin(ONE + nx*ny*i,(x-sico_x[i]),(y-sico_y[i]), nx));
      ttmp[0] = par*tmp[0];
      ttmp[1] = fabs(fbilin(TWO + nx*ny*i,(x+sico_x[i]),(y+sico_y[i]), nx));
      ttmp[2] = fabs(fbilin(TWO + nx*ny*i,(x-sico_x[i]),(y-sico_y[i]), nx));
      
      if (type[xypos] == EDGETYPE && maxind(tmp,3) == 0)
	nms[xypos] = maxmap[xypos];
      else if ((maxind(tmp,3) == 0)&&(maxind(ttmp,3) == 0)) 
	nms[xypos] = maxmap[xypos];
      /* old method 
	 if ((maxind(tmp,3) == 0)&&(maxind(ttmp,3) == 0)) 
	 nms[xypos] = maxmap[xypos];
	 */
    }
  }
  
  
  /* key-points */
  if (kpt_map) {
    border=keypoint_map(modulus,nori,nx,ny,Xridge,gamma,
			comp_gain,even_gain,kpt_map);
    fprintf(stderr,"Compute key-points\n");
    fflush(stderr);
  }

  /* quality */
  if (quality) {
    quality_dft(modulus,locori,nx*ny,nori,orisel,quality);
    fprintf(stderr,"Compute general edge quality\n");
    fflush(stderr);
  }

  /* if (deg) for (i=0;i<nx*ny;i++) locori[i] *= PPi;  */
  
  /* free intern maps */
  
  free(fimage);
  free(even);
  free(odd);
  free(modulus);
}

