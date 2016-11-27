/* 
 * $Log: egabor.c,v $
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
 *          ohenri@vision.ethz.ch
 * 
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include "lowlevel.h"

#define PI  3.141592654
#define PI2 6.283185308
#define SWEEPFUNC(X,K,A) ((K)*((float) erf(sqrt((double)(A))*((double)(X)))))
#define GAUSS(XX,S) (fexp(-0.5*(XX)*(XX)/((S)*(S))))
#define ORIWRAP(i,n) (oriwrap[(i) + (n)])
#define ZERO_NUM 100
#define ARR_SIZE 1024
#define FFT_SIZE 1024
#define Z(I) (ARR_SIZE + (I))

#define RESOLUTION 512
#define LAMBDA_LOW 0.0
#define LAMBDA_UP 2.0
#define NULL_MALLOC(x,n,type) { if (NULL == ((x)=(type *)malloc((n)*sizeof(type)))) {return (NULL);}}

#define COORD(x,y) ((y)*nx+(x))
#define REAL(buf,x,y,ny) (buf[2 * ( (x) * (ny) + (y) )])
#define IMAG(buf,x,y,ny) (buf[2 * ( (x) * (ny) + (y) ) + 1])

static float *wr, *wi;
static int *bitrev;

/*********************************************************************/
int embed(float *fimage,int nx,int ny,float **eimage,int *newx,int *newy)
{
  int i,j,horiz,verti,xoff,yoff,nnx,nny;
   
  for (nnx=2;nnx<nx;nnx += nnx) {}
  for (nny=2;nny<ny;nny += nny) {}
  
  if ((nnx == nx) && (nny == ny)){
    *eimage= (float *) malloc(nnx*nny*sizeof(float));
    memcpy((void *) *eimage,fimage,nnx*nny*sizeof(float));
    *newx=nx;*newy=ny;
    return 0;
  }
  else {
    *newx=nnx;*newy=nny;
    *eimage = (float *) malloc(nnx*nny*sizeof(float));
    xoff = (nnx-nx)/2;
    yoff = (nny-ny)/2;
    
    for (j=yoff;j<ny+yoff;j++){
      for (i=xoff;i<nx+xoff;i++)
	(*eimage)[i+j*nnx]=fimage[i-xoff+(j-yoff)*nx];
    }
    horiz = xoff%nx;
    verti = yoff%ny;
    
    for (i=horiz;i<nnx-horiz;i++) {
      for (j=0;j<yoff;j++)
	(*eimage)[i+j*nnx]=fimage[(i-horiz)%nx+((j+ny-verti)%ny)*nx];

      for (j=ny+yoff;j<nny;j++)
	(*eimage)[i+j*nnx]=fimage[(i-horiz)%nx+((j-ny-verti)%ny)*nx];
    }

    for (j=verti;j<nny-verti;j++) {
      for (i=0;i<xoff;i++)
	(*eimage)[i+j*nnx]=fimage[(i+nx-horiz)%nx+((j-verti)%ny)*nx];
      for (i=nx+xoff;i<nnx;i++)
	(*eimage)[i+j*nnx]=fimage[(i-nx-horiz)%nx+((j-verti)%ny)*nx];
    }
    
    for (i=0;i<horiz;i++) {
      for (j=0;j<verti;j++)
	(*eimage)[i+j*nnx]=fimage[(i+nx-horiz)%nx+((j+ny-verti)%ny)*nx];
      for (j=nny-verti;j<nny;j++)
	(*eimage)[i+j*nnx]=fimage[(i+nx-horiz)%nx+((j-ny-verti)%ny)*nx];
    }

    for (i=nnx-horiz;i<nnx;i++) {
      for (j=0;j<verti;j++)
	(*eimage)[i+j*nnx]=fimage[(i-nx-horiz)%nx+((j+ny-verti)%ny)*nx];
      for (j=nny-verti;j<nny;j++)
	(*eimage)[i+j*nnx]=fimage[(i-nx-horiz)%nx+((j-ny-verti)%ny)*nx];
    }	
    
    return 1;
  }

}

void crop(float *map,int nx,int ny,int nnx,int nny)
{
  int xoff,yoff,x,y,i=0;
  
  xoff = (nnx-nx)/2;
  yoff = (nny-ny)/2;
    
  if (xoff || yoff){
    for (y=yoff;y<yoff+ny;y++)
      for (x=xoff;x<xoff+nx;x++)
	map[i++]=map[x+y*nnx];
  }
}

void egabor_filter(float *fimage,int nx,int ny,
		   float sigma,float sweep,float lambda,int nori,
		   float ori_sel,float *even,float *odd,float *modulus)
{
  int eflag,inx,iny;
  float *fft_image,*eimage=NULL;

  /*=======================================================================
   * variables declarations egabor_avs.c
   *=======================================================================
   */
  
  float ebuffer[2*ARR_SIZE+1], obuffer[2*ARR_SIZE+1],
        gbuffer1[2*ARR_SIZE+1], gbuffer2[2*ARR_SIZE+1];
  float espect_r[ARR_SIZE+1], espect_i[ARR_SIZE+1],emax,ebawi,sp_sum;
  float ospect_r[ARR_SIZE+1], ospect_i[ARR_SIZE+1],omax,obawi,sp_max;
  float gval, swval[2*ARR_SIZE+1];
  float x_arr[2*ARR_SIZE+1], xx;
  int i, n, iemax, iomax, iehalf1, iehalf2, iohalf1, iohalf2;
  
  /*=======================================================================
   * variables declarations egabor_filter_avs.c
   *=======================================================================
   */
  
  
  float start_angle=0.0, ang;

  int nu, nv, pos, ori, sx, sy, i_size, offset,offset_real,offset_imag;
  float *out_e_data, *out_o_data, *e_data, *o_data;
  float si, co, cosx,cosy,sisx,sisy,angle, alpha, radius_scale;
  float rotco, rotsi, dumx, dumy, dumr, dumdumr, rsize;
  float ori_sel_weight, even_val_r, even_val_i, odd_val_r, odd_val_i;
  
  inx=nx;iny=ny;
  eflag=embed(fimage,inx,iny,&eimage,&nx,&ny);
  fft_image=cvrtFloatToFourier(eimage,nx,ny);
  free(eimage);

  i_size=nx*ny;
  memset ((void *) ebuffer, 0, (2*ARR_SIZE + 1)*sizeof (float));
  memset ((void *) obuffer, 0, (2*ARR_SIZE + 1)*sizeof (float));

  /*===========================================================================
   * calculate the sweep function and spatial filter functions
   *==========================================================================
   */
  
  swval[0] = 0.0;
  for (i = 1; i <= ARR_SIZE; i++) {
    xx = (float) i / (2.0*sigma);
    x_arr[Z(i)]  = (float)  i*0.5;
    x_arr[Z(-i)] = (float) -i*0.5;
    swval[Z(i)]  = SWEEPFUNC (xx, sweep, lambda);
    swval[Z(-i)] = -swval [Z(i)];
  }

  ebuffer[Z(0)]  = 1.0;
  obuffer[Z(0)]  = 0.0;
  gbuffer1[Z(0)] = 1.0;
  gbuffer2[Z(0)] = -1.0;

  for (i = 1; i <= ARR_SIZE; i++) {
    gval = GAUSS (x_arr[Z(i)], sigma);

    gbuffer1[Z(i)]  = gval;
    gbuffer1[Z(-i)] = gval;

    ebuffer[Z(i)]   = fcos(2.0*M_PI*swval[Z(i)])*gval;
    ebuffer[Z(-i)]  = ebuffer[Z(i)];

    obuffer[Z(i)]   = fsin(2.0*M_PI*swval[Z(i)])*gval;
    obuffer[Z(-i)]  = -obuffer[Z(i)];

    gbuffer2[Z(i)]  = -1*fsqrt (SQR(ebuffer[Z(i)]) + SQR(obuffer[Z(i)]));
    gbuffer2[Z(-i)] = gbuffer2[Z(i)];

    if (gval < 1.0e-5) break;
  }

  /*=======================================================================
   * calculate 1D radial spectra of filter functions
   *=======================================================================
   */

  espect_r[0] = ebuffer[Z(0)];
  espect_i[0] = 0.0;
  ospect_r[0] = 0.0;
  ospect_i[0] = 0.0;

  for (i = 1; i <= (FFT_SIZE / 2); i++) {
    espect_r[i]            = ebuffer[Z(i)];
    espect_r[FFT_SIZE - i] = ebuffer[Z(i)];
    espect_i[i]            = 0.0;
    espect_i[FFT_SIZE - i] = 0.0;
    

    ospect_r[i]            = obuffer[Z(i)];
    ospect_r[FFT_SIZE - i] = -obuffer[Z(i)];
    ospect_i[i]            = 0.0;
    ospect_i[FFT_SIZE - i] = 0.0;
  }

  fft_alloc (FFT_SIZE, 1);
  fft (espect_r, espect_i, FFT_SIZE);
  for (i = 0; i < FFT_SIZE; i++) {
    espect_r[i] /= (float) FFT_SIZE;
    espect_i[i] /= (float) FFT_SIZE;
  }
  fft_free ();

  fft_alloc (FFT_SIZE, 1);
  fft (ospect_r, ospect_i, FFT_SIZE);
  for (i = 0; i < FFT_SIZE; i++) {
    ospect_r[i] /= (float) FFT_SIZE;
    ospect_i[i] /= (float) FFT_SIZE;
  }
  fft_free ();
   
  /*
   * We have to respect that these radial filter files are prepared for
   * 512 x 512 images. The actual range of the radial frequency axis is
   * [-512,512]. If the image has a different size than 512 x 512, we have 
   * to adjust the radial frequency.
   */

  rsize = (float) (ARR_SIZE/2);
  
  /*==========================================================================
  * Main loop over all nori orientations
  *==========================================================================
  */
  
  e_data=(float *) malloc((nx/2+1)*2*ny*sizeof(float));
  o_data=(float *) malloc((nx/2+1)*2*ny*sizeof(float));
  fprintf(stderr,"\nFilter orientation (in degrees from x-axis) [0]\b\b\b");
  fflush(stderr);
  
  sx=((int) MAX(nx,ny))/nx;
  sy=((int) MAX(nx,ny))/ny;
  
  for (ori = 0; ori < nori; ori++) {
    
    alpha=180.0 * (float) ori/(float) nori;
    if (alpha < 10.0) fprintf(stderr,"[%1.0f]\b\b\b",alpha);
    else if (alpha < 100.0) fprintf(stderr,"[%2.0f]\b\b\b\b",alpha);
    else fprintf(stderr,"[%3.0f]\b\b\b\b\b",alpha);
    fflush(stderr);
    
    angle = M_PI*(float) ori / (float) nori + M_PI * start_angle / 180.0;
    
    memset ((void *) e_data, 0, (nx/2+1)*2*ny*sizeof (float));
    memset ((void *) o_data, 0, (nx/2+1)*2*ny*sizeof (float));
    
    si = fsin (angle);
    co = fcos (angle);
    sisx = si*sx; sisy = si*sy; cosx = co*sx; cosy = co*sy;
    radius_scale=rsize / (float) MAX(nx,ny);

    /*========================================================================
     * start loops over nu & nv (coordinates in Fourier plane)
     *========================================================================
     */
    
    for (nu = 0; nu < (nx / 2 + 1); nu++) {
      rotco = cosx*nu;
      rotsi = sisx*nu;
      offset= nu * ny;

      for (nv = 0; nv < (ny / 2); nv++) {
	dumx =  rotco + sisy*nv;
	dumy = -rotsi + cosy*nv;
	dumr = fsqrt (SQR (dumx) + SQR (dumy));
	if (dumr > 0.0) {
	  ang = fatan2 (dumx, dumy);
	  ori_sel_weight = fpow (ffabs (fcos (ang)), ori_sel);
	}
	else {
	  ori_sel_weight = 0.0;
	}
	
	dumdumr = dumr * radius_scale;
	if (dumdumr > rsize) {
	  even_val_r = 0.0;
	  odd_val_i  = 0.0;
	}
	else {
	  even_val_r = FBI (espect_r, dumdumr) * ori_sel_weight;
	  odd_val_i  = FBI (ospect_i, dumdumr) * ori_sel_weight;
	}
	
	if (dumy > 0.0) {
	  odd_val_i = -odd_val_i;
	}
	
	offset_real=2*(nv+offset); offset_imag=offset_real+1;
	e_data[offset_real]=fft_image[offset_real]*even_val_r;
	e_data[offset_imag]=fft_image[offset_imag]*even_val_r;
	o_data[offset_real]=fft_image[offset_imag]*(-odd_val_i);
	o_data[offset_imag]=fft_image[offset_real]*odd_val_i;

      }
      
      for (nv = (-ny / 2); nv < 0; nv++) {
	dumx =  rotco + sisy*nv;
	dumy = -rotsi + cosy*nv;
	dumr = fsqrt (SQR (dumx) + SQR (dumy));
	if (dumr > 0.0) {
	  ang = fatan2 (dumx, dumy);
	  ori_sel_weight = fpow (ffabs (fcos(ang)), ori_sel);
	}
	else {
	  ori_sel_weight = 0.0;
	}
	
	dumdumr = dumr * radius_scale;
	if (dumdumr > rsize) {
	  even_val_r = 0.0;
	  odd_val_i  = 0.0;
	}
	else {
	  even_val_r = FBI (espect_r, dumdumr) * ori_sel_weight;
	  odd_val_i  = FBI (ospect_i, dumdumr) * ori_sel_weight;	
	}
	
	if (dumy > 0.0) {
	  odd_val_i = -odd_val_i;
	}
	
	offset_real=2*(nv+ny+offset); offset_imag=offset_real+1;
	e_data[offset_real]=fft_image[offset_real]*even_val_r;
	e_data[offset_imag]=fft_image[offset_imag]*even_val_r;
	o_data[offset_real]=fft_image[offset_imag]*(-odd_val_i);
	o_data[offset_imag]=fft_image[offset_real]*odd_val_i;

      }
    }
    
    /*========================================================================
     * inverse fourier transformations and local modulus calculation
     *========================================================================
     */
    
    out_o_data=cvrtFourierToFloat(o_data,nx,ny);
    out_e_data=cvrtFourierToFloat(e_data,nx,ny);
    crop(out_o_data,inx,iny,nx,ny);
    crop(out_e_data,inx,iny,nx,ny); 
    
    pos=((nori-ori)%nori)*inx*iny;
    if (even) memcpy ((void *) (even+pos),(void *) out_e_data,
		      inx*iny*sizeof (float));
    if (odd)  memcpy ((void *) (odd+pos), (void *) out_o_data, 
		      inx*iny*sizeof (float));  
    for (i = 0; i < inx*iny; i++) {
      modulus[i+pos] = fsqrt (SQR (out_e_data[i]) + SQR (out_o_data[i]));
    }
    free(out_e_data);
    free(out_o_data); 
  }
  
  free(e_data);
  free(o_data);
  free(fft_image);

  fprintf(stderr,"\n");
}


/* ------------------------------------------------------------------
 *
 *     function to calculate lambda to integrate even-filter to 0 ...
 *
 */

float egabor_params(float sigma,float sweep,float lambda1,float lambda2,int count,int which_min)
{
  unsigned int i, k, l, min_count;
  float lambda, x;
  float ebuffer[1024], lambda_arr[1024], sumval[1024], swval[1024], gval[1024];
  float lambda_min[ZERO_NUM], sumval_min[ZERO_NUM];

  for (i = 0; i < 1024; i++) {
    gval[i]  = GAUSS ((float) i*0.5, sigma);
    if (gval[i] < 1.0e-38) break;
  }
  for (; i < 1024; i++) {
    gval[i] = 0.0;
  }

  min_count = 0;
  
  for (i = 0; i < count; i++) {
    lambda = lambda1 + (lambda2 - lambda1)*(float) i / (float) count;
    lambda_arr[i] = lambda;

    sumval[i]  = 1.0;
    ebuffer[0] = 1.0;
    for (k = 1; k < 1024; k++) {
      x = (float) k / (2.0*sigma);
      swval[k] = SWEEPFUNC (x, sweep, lambda);
      ebuffer[k] = fcos(2.0*M_PI*swval[k])*gval[k];

      sumval[i] += (2.0*ebuffer[k]);
      if (gval[k] < 1.0e-38) break;
    }

    if (i < 1) continue;
    if ((sumval[i]*sumval[i-1]) < 0.0) {
      if (min_count >= ZERO_NUM) {
	printf ("too many zero-crossings!\n");
	exit (-1);
      }
      lambda_min[min_count] = lambda_arr[i-1] - sumval[i-1]*
	(lambda_arr[i] - lambda_arr[i-1])/(sumval[i] - sumval[i-1]);
      ebuffer[0] = 1.0;
      sumval_min[min_count] = 1.0;
      for (k = 1; k < 1024; k++) {
	x = (float) k / (2.0*sigma);
	swval[k] = SWEEPFUNC (k, sweep, lambda_min[min_count]);
	ebuffer[k] = fcos(2*M_PI*swval[k])*gval[k];
	sumval_min[min_count] += (2.0*ebuffer[k]);
	if (gval[k] < 1.0e-38) break;
      }
      min_count++;
    }
  }
  
  /* for (i = 0; i < min_count; i++) {
     printf ("Minimum at: ( %e , %e )\n", lambda_min[i],sumval_min[i]);
     }
     */
  return(lambda_min[which_min]);
}


/*****************************************************************************
 * fourier-trafo
 * =============
 *
 * Fourier-transformation is based on algorithm presented in the book
 * "Numerical Recipes in C - The Art of Scientific Computing" by
 * William H. Press, Brian P. Flannery, Saul A. Teukolsky &
 * William T. Vetterling
 *----------------------------------------------------------------------------
 *
 * general complex 1-dim Fourier-trafo
 *
 *****************************************************************************
 */

/***************************************************************************
 * void four1 (float data[], int nn, int isign)
 *---------------------------------------------------------------------------
 * ATTENTION!!! data[1..n],   N O T data [0..(n-1)]
 * The data array conforms to the FORTRAN-Standard, so in C the first element
 * data[0] is void !!
 **************************************************************************
 */

void four1 (float *data, int nn, int isign)
{
  int n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi;
  
  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      SWAP (data [j], data [i]);
      SWAP (data [j+1], data [i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = 2*mmax;
    theta = 6.28318530717959 / (isign*mmax);
    wtemp = sin (0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin (theta);
    wr  = 1.0;
    wi  = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
	j = i + mmax;
	tempr = wr*data [j]   - wi*data [j+1];
	tempi = wr*data [j+1] + wi*data [j];
	data [j]   = data [i] - tempr;
	data [j+1] = data [i + 1] - tempi;
	data [i] += tempr;
	data [i + 1] += tempi;
      }
      wr = (wtemp=wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}

/*****************************************************************************
 * void realft (float data[], int n, int isign)
 *----------------------------------------------------------------------------
 * ATTENTION!!! data[1..2n],   N O T data [0..(2n-1)]
 * The data array conforms to the FORTRAN-Standard, so in C the first element
 * data[0] is void !!
 * ATTENTION!!! n is actual n/2. So this is a real transform of 2n elements !!
 ******************************************************************************
 */
void realft (float *data, int n, int isign)
{
  int i, i1, i2, i3, i4, n2p3;
  float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  double wr, wi, wpr, wpi, wtemp, theta;
  
  theta = 3.141592653589793 / (double) n;
  if (isign == 1) {
    c2 = -0.5;
    four1 (data, n, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin (0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin (theta);
  wr = 1.0 + wpr;
  wi = wpi;
  n2p3 = 2*n + 3;
  for (i = 2; i <= n / 2; i++) {
    i4 = 1 + (i3=n2p3 - (i2=1 + (i1=i + i - 1)));
    h1r = c1*(data [i1] + data [i3]);
    h1i = c1*(data [i2] - data [i4]);
    h2r = -c2*(data [i2] + data [i4]);
    h2i =  c2*(data [i1] - data [i3]);
    data [i1] =  h1r + wr*h2r - wi*h2i;
    data [i2] =  h1i + wr*h2i + wi*h2r;
    data [i3] =  h1r - wr*h2r + wi*h2i;
    data [i4] = -h1i + wr*h2i + wi*h2r;
    wr = (wtemp=wr)*wpr - wi*wpi + wr;
    wi = wi*wpr + wtemp*wpi + wi;
  }
  if ( isign == 1) {
    data [1] = (h1r = data [1]) + data [2];
    data [2] = h1r - data [2];
  } else {
    data [1] = c1*((h1r=data [1]) + data [2]);
    data [2] = c1*(h1r - data [2]);
    four1(data, n, -1);
  }
}

/*****************************************************************************
 * float *cvrtFloatToFourier (float orig[], int nx, int ny)
 *----------------------------------------------------------------------------
 * The first parameter orig is the input array (which is one-dimesnional, but
 * represents a two-dimesnional image, row after row. It starts with the
 * element orig [0], as the standard C convention!
 *****************************************************************************
 */

float *cvrtFloatToFourier (float *orig, int nx, int ny)
{
  int i, k, nx2, ny2, nn;
  float *data, *tmp, *tmpptr, *res, *resptr;
  
  nx2 = nx / 2;
  ny2 = ny / 2;
  nn = (nx2 + 1)*ny*2;
  NULL_MALLOC (data, (nx*ny + 1), float);
  NULL_MALLOC (res, nn, float);
  NULL_MALLOC (tmp, (ny*2 + 1), float);
  memcpy((void *) &data[1],(void *) orig,nx*ny*sizeof(float));
  
  for (i = 0; i < ny; i++) {
    realft (&data [COORD(0,i)], nx2, 1);
  }
  
  resptr = res;
  tmpptr = tmp + 1;
  for (k = 0; k < ny; k++) {
    *tmpptr++ = data [COORD(1,k)];
    *tmpptr++ = 0.0;
  }
  four1 (tmp, ny, 1);
  memcpy ((void *) resptr, (void *) &tmp [1], 2*ny*sizeof (float));
  resptr += 2*ny;
  
  for ( i = 3; i < nx; i += 2) {
    tmpptr = tmp + 1;
    for (k = 0; k < ny; k++) {
      *tmpptr++ = data [COORD(i,k)];
      *tmpptr++ = data [COORD(i + 1,k)];
    }
    four1 (tmp, ny, 1);
    memcpy ((void *) resptr, (void *) &tmp [1], 2*ny*sizeof (float));
    resptr += 2*ny;
  }

  tmpptr = tmp + 1;
  for (k = 0; k < ny; k++) {
    *tmpptr++ = data [COORD(2,k)];
    *tmpptr++ = 0.0;
  }
  four1 (tmp, ny, 1);
  memcpy ((void *) resptr, (void *) &tmp [1], 2*ny*sizeof (float));
  
  free (tmp);
  free (data);
  return (res);
}

/*****************************************************************************
 * float *cvrtFourierToFloat (float orig[], int nx, int ny)
 *----------------------------------------------------------------------------
 * The first parameter orig is the input array (which is one-dimensional, but
 * represents a the fourier-trafo of a two-dimesnional image. There are
 * (nx / 2 + 1) columns (one after the other), each 2*ny elements long. The
 * complex number are stored with the real part immediately before the
 * imaginary part. The first element is at orig [0] (standard C convention) !
 *****************************************************************************
 */

float *cvrtFourierToFloat (float *orig, int nx, int ny)
{
  int i, k, nx2, ny2, nn, nyggr2, i_size=nx*ny;
  float *data, *tmp, *tmpptr, *res, *resptr;

  nyggr2= ny * 2; 
  nx2 = nx / 2;
  ny2 = ny / 2;
  nn = (nx2 + 1)*nyggr2;
  NULL_MALLOC (data, (nn + 1), float);
  NULL_MALLOC (tmp, (MAX(ny,nx)*2 + 1), float);
  NULL_MALLOC (res, (i_size), float);
  memcpy ((void *) &data[1],(void *) orig, nn*sizeof (float));

  resptr = data;
  for (i = 0; i < (nx2 + 1); i++) {
    four1 (resptr, ny, -1);
    resptr += nyggr2;
  }
  
  resptr = res;
  for (i = 1; i <= nyggr2; i += 2) {
    tmpptr = tmp + 1;
    *tmpptr++ = data [i];
    *tmpptr++ = data [nx2*nyggr2 + i];
    for (k = 1; k < nx2; k++) {
      *tmpptr++ = data [k*nyggr2 + i];
      *tmpptr++ = data [k*nyggr2 + i + 1];
    }
    realft (tmp, nx2, -1);
    memcpy ((void *) resptr,(void *) &tmp [1], nx*sizeof (float));
    resptr += nx;
  }
  
  for (i = 0; i < i_size; i++) {
    res [i] = res [i] / (float) (i_size);
  }

  free (tmp);
  free (data);
  
  return (res);
}

void fft_alloc (int N,int fft_type)
{
   int i, test_bit, set_bit, max_bit;
   float r;

   r = PI2 / (float) N;
   wr = (float *) malloc(N * sizeof(float));
   wi = (float *) malloc(N * sizeof(float));

   if (fft_type == 1) {
      for (i = 0; i < N; i++) {
         wr [i] =  cos (r * (float) i);
         wi [i] = -sin (r * (float) i);
      }
   } else if (fft_type == -1) {
      for (i = 0; i < N; i++) {
         wr [i] = cos (r * (float) i);
         wi [i] = sin (r * (float) i);
      }
   }
   

   max_bit = N >> 1;
   bitrev = (int *) malloc (N * sizeof(int));
   for (i = 0; i < N; i++) {
      bitrev[i] = 0;
      set_bit = max_bit;
      for (test_bit = 1; set_bit > 0; test_bit <<= 1, set_bit >>= 1) {
         if ((test_bit & i) == test_bit) {
            bitrev[i] = bitrev[i] | set_bit;
         }
      }
   }
   return;
}

void fft_free (void)
{
  free (wr);
  free (wi);
   free (bitrev);
  return;
}


void fft (float *re, float *im, int N)
{
   int n, N2;
   int c2;
   register int i, ii, a, c;

   register float tr, ti;

   N2 = N >> 1;

   for (i = 0; i < N; i++) {
      if (bitrev [i] > i) {
         tr = re [i];
         ti = im [i];
         re [i] = re [bitrev [i]];
         im [i] = im [bitrev [i]];
         re [bitrev [i]] = tr;
         im [bitrev [i]] = ti;
      }
   }


   n = N;
   for (c = 1; c < N; c <<= 1) { /* c  = 1, 2, 4, ..., N/4 */
      c2 = c << 1;               /* c2 = 2, 4, 8, ..., N/8 */
      n >>= 1;                   /* n  = N/2, N/4, N/8, ..., 1 */

      for (a = 0; a < N; a += c2) {
         ii = 0;
         for (i = 0; i < c; i++) {

/*
            d [a + i]     = d [a + i] + d [a + c + i] * w [i*n]; 
            d [a + c + i] = d [a + i] + d [a + c + i] * w [N2 + i*n];
                                                            \___,___/
                                                                !
                                                          - w[(i-1)*n]
*/


            tr = re [a + c + i] * wr [ii] - im [a + c + i] * wi [ii];
            ti = re [a + c + i] * wi [ii] + im [a + c + i] * wr [ii];

            re [a + c + i] = re [a + i] - tr;
            im [a + c + i] = im [a + i] - ti;

            re [a + i] += tr;
            im [a + i] += ti;

            ii += n;

         } /* i */
      } /* a */
   } /* c */
   return;
}
