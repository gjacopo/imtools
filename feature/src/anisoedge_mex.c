/*
 * This file contains the source code for the Robust Anisotropic Diffusion
 * edge detector described in "Robust Anisotropic Diffusion", by M. Black,
 * G. Sapiro, D. Marimont and D. Heeger, IEEE Trans. Image Processing, vol.7,
 * no. 3, pp. 421-432, Mar 1998. It was written by Christine Kranenburg
 * (kranenbu@bigpine.csee.usf.edu).
 *
 * Compilation
 * 
 * gcc -o black_edge black_edge.c pgm.c non_max.c -lm -O3
 * usage: black_edge <image> <coefficient>
 * coefficient = float > 0.0
 * = 1.0  is default
 * outputs: <image>_smooth.pgm  - the smoothed image
 * <image>_s_<coeff>.pgm  - the binary edge image
 *
 *
 * Features: 
 *
 * The implementation was verified by comparison of the MAD when run on the
 * Canal image. The detector was then supplemented with non-maximal
 * suppression to produce single pixel wide edges.
 *
 * One other modification occurs in the function "spatial_disconts".
 * We threshold the intesity difference between the center pixel and its 4
 * neighbors. The original implementation only thresholded the top and left
 * neighbors. Using 4 neighbors instead of 2 produces thicker but more
 * continuous edges. The non-max then performs edge thinning.
 *
 * Known issues:
 *
 * Anything that causes sigma = 0 may produce unpredictable results because
 * of the division by sigma squared term in the tukey function. For this
 * reason, 0 is an invalid parameter to supply to the detector. The same
 * effect may appear on images in which a substantial part of the image is
 * the same color, causing the MAD function to return 0.
 *
 * CJK - 5/8/98
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

void black(double *smooth, int rows, int cols, double sigma, double lambda);
double psi(double x, double sigma);
int find_median(int *grad, int N);
int MAD(double *smooth, int rows, int cols);
void spatial_disconts(double *image, unsigned char *nms, unsigned char
			*edge, int rows, int cols, double sigma);

void derivative_x_y(double *smoothedim, int rows, int cols,
        double *delta_x, double *delta_y);
void magnitude_x_y(double *delta_x, double *delta_y, int rows, int cols,
        double *magnitude);
void non_max_supp(double *mag, double *gradx, double *grady, int nrows, int ncols,
    unsigned char *result);
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   double *input,*image;   // The input image
   double *smooth;	         // The *smoothed image
   unsigned char *nms;
   unsigned char *edge;	         // The output edge image
   int rows, cols;               // The dimensions of the image.
   double sigma_e, sigma, slide = 1.0;
   double lambda;	         // smoothing rate
   int i, k;			 // loop counters
   double *im_int, *dx, *dy, *mag;
   int niter, flag_max;

// inputs
    /*  create a pointer to the input matrix y */
    input = mxGetPr(prhs[0]);;
    /*  get the dimensions of the matrix input y */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);
    
// <coefficient> slide floating point # > 0.0 (default=1.0)\n"); 
       if( !mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgTxt("Input sigma must be a scalar.");
    } else
        slide = mxGetScalar(prhs[1]);
    // number of iterative smoothing: defult iter=100
    if( !mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgTxt("Input low must be a scalar.");
    } else
        niter = floor(mxGetScalar(prhs[2]));
    if( !mxIsLogicalScalar(prhs[3])) {
        mexErrMsgTxt("Input flag_max must be a logical scalar.");
    } else
        flag_max = (int)(mxGetScalar(prhs[3]));

// outputs
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericMatrix(rows, cols, mxUINT8_CLASS, mxREAL);
    /*  create a C pointer to a copy of the output matrix */
    edge = mxGetData(plhs[0]);
    if(nlhs == 2) {
         plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        mag = mxGetPr(plhs[1]);
    } else if((mag = mxMalloc(rows*cols*sizeof(double))) == NULL){
      mexErrMsgTxt("Error allocating the magnitude image.\n");
      exit(1);
   }

// Copy original image to working-copy-image
   if((smooth=mxMalloc(rows*cols*sizeof(double)))==NULL)
      {fprintf(stderr, "Error allocating the smooth image.\n");
       exit(1);
      }

   for(i = 0; i < rows*cols; i++)
      smooth[i] = input[i];

   sigma_e = 1.4826 * (double)(MAD(smooth, rows, cols));
   sigma_e = slide * sigma_e;
   sigma = sigma_e * sqrt(5.0);
   lambda = 1.0 / psi(sigma_e, sigma);

// Smooth iteratively niter times
   
   for(k = 0; k < niter; k++)
       black(smooth, rows, cols, sigma, lambda);

   if((nms=mxMalloc(rows*cols*sizeof(unsigned char)))==NULL)
      {mexErrMsgTxt("Error allocating the nms image.\n");
       exit(1);
      }

   
   if((dx = mxMalloc(rows*cols*sizeof(double))) == NULL){
      mexErrMsgTxt("Error allocating the delta_x image.\n");
      exit(1);
   }
   if((dy = mxMalloc(rows*cols*sizeof(double))) == NULL){
      mexErrMsgTxt("Error allocating the delta_x image.\n");
      exit(1);
   }

   derivative_x_y(smooth, rows, cols, dx, dy);
   magnitude_x_y(dx, dy, rows, cols, mag);
   non_max_supp(mag, dx, dy, rows, cols, nms);

   if (!flag_max)	memset(nms, 128, rows*cols);

   spatial_disconts(smooth, nms, edge, rows, cols, sigma_e);
 
    free(dx);
   free(dy);
   if(nlhs!=2) free(mag);
   free(nms);
   free(smooth);
}


void black(double *smooth, int rows, int cols, double sigma, double lambda)
{
int i, j, index;
double north, south, east, west, delta_i;


// Divide pixels into checkerboard and update white and black pixels separately.

// Iterate over white sites (odd pixels)
   for(i = 1; i < rows-1; i++)
      for(j = 1; j < cols-1; j=j+2)
	{
	index = i * cols + j;
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	east = smooth[index + 1] - smooth[index];
	west = smooth[index - 1] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(east, sigma) + psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}


// Iterate over black sites (even pixels)
   for(i = 1; i < rows-1; i++)
      for(j = 2; j < cols-1; j=j+2)
	{
	index = i * cols + j;
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	east = smooth[index + 1] - smooth[index];
	west = smooth[index - 1] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(east, sigma) + psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

// Update edge pixels (3 neighbors)
   // Top edge pixels
      i = 0;		
      for(j = 1; j < cols-1; j=j+2)
	{
	west = smooth[j - 1] - smooth[j];
	east = smooth[j + 1] - smooth[j];
	south = smooth[j + cols] - smooth[j];
	delta_i = 0.25 * lambda * (psi(south, sigma) + psi(east, sigma) +
		  psi(west, sigma));
	smooth[j] = smooth[j] + delta_i;
	}

      for(j = 2; j < cols-1; j=j+2)
	{
	west = smooth[j - 1] - smooth[j];
	east = smooth[j + 1] - smooth[j];
	south = smooth[j + cols] - smooth[j];
	delta_i = 0.25 * lambda * (psi(south, sigma) + psi(east, sigma) +
		  psi(west, sigma));
	smooth[j] = smooth[j] + delta_i;
	}

   // Bottom edge pixels
      i = rows - 1;
      for(j = 1; j < cols-1; j=j+2)
	{
	index = i * cols + j;
	west = smooth[index - 1] - smooth[index];
	east = smooth[index + 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(east, sigma) +
		  psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

      for(j = 2; j < cols-1; j=j+2)
	{
	index = i * cols + j;
	west = smooth[index - 1] - smooth[index];
	east = smooth[index + 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(east, sigma) +
		  psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

   // Left edge pixels
      j = 0;
      for(i = 1; i < rows-1; i=i+2)
	{
	index = i * cols + j;
	east = smooth[index + 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(east, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

      for(i = 2; i < rows-1; i=i+2)
	{
	index = i * cols + j;
	east = smooth[index + 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(east, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

   // Right edge pixels
      j = cols - 1;
      for(i = 1; i < rows-1; i=i+2)
	{
	index = i * cols + j;
	west = smooth[index - 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

      for(i = 2; i < rows-1; i=i+2)
	{
	index = i * cols + j;
	west = smooth[index + 1] - smooth[index];
	north = smooth[index - cols] - smooth[index];
	south = smooth[index + cols] - smooth[index];
	delta_i = 0.25 * lambda * (psi(north, sigma) + psi(south, sigma) +
		  psi(west, sigma));
	smooth[index] = smooth[index] + delta_i;
	}

// Update corner pixels
     i = 0;
     j = 0;
     index = i * cols + j;
     east = smooth[index + 1] - smooth[index];
     south = smooth[index + cols] - smooth[index];
     delta_i = 0.25 * lambda * (psi(east, sigma) + psi(south, sigma));
     smooth[index] = smooth[index] + delta_i;

     j = cols - 1;
     index = i * cols + j;
     west = smooth[index - 1] - smooth[index];
     south = smooth[index + cols] - smooth[index];
     delta_i = 0.25 * lambda * (psi(west, sigma) + psi(south, sigma));
     smooth[index] = smooth[index] + delta_i;

     i = rows - 1;
     j = 0;
     index = i * cols + j;
     east = smooth[index + 1] - smooth[index];
     north = smooth[index - cols] - smooth[index];
     delta_i = 0.25 * lambda * (psi(east, sigma) + psi(north, sigma));
     smooth[index] = smooth[index] + delta_i;

     j = cols - 1;
     index = i * cols + j;
     west = smooth[index - 1] - smooth[index];
     north = smooth[index - cols] - smooth[index];
     delta_i = 0.25 * lambda * (psi(west, sigma) + psi(north, sigma));
     smooth[index] = smooth[index] + delta_i;
}


double psi(double x, double sigma)
{
double /*perona,*/tukey, temp;

if (fabs(x) <= sigma)
   {temp = x / sigma;
    tukey = x * (1.0 - temp * temp) * (1.0 - temp * temp);
   }
else
    tukey = 0.0;

//perona = (2.0 * x * sigma * sigma) / ((2.0 * sigma * sigma) + (x * x));

return tukey;
}


void spatial_disconts(double *image, unsigned char *nms,
			unsigned char *edge, int rows, int cols, double sigma)
{
int i, j, index;

for (i = 0; i < rows*cols; i++)
    edge[i] = 255;

for (i = 1; i < rows-1; i++)
    for (j = 1; j < cols-1; j++)
	{
	index = i * cols + j;
	if(nms[index] == 128 && 
		((fabs(image[index] - image[index-cols]) >= sigma) ||
		(fabs(image[index] - image[index-1]) >= sigma) ||
		(fabs(image[index] - image[index+cols]) >= sigma) ||
		(fabs(image[index] - image[index+1]) >= sigma)))
		edge[index] = 0;
	}
}


// Median absolute deviation function - used to calculate sigma_e
int MAD(double *smooth, int rows, int cols)
{
int *north, *south, *east, *west, *gm;
int median;
int i, j, index;

north=(int*)calloc((rows*cols),sizeof(int));
south=(int*)calloc((rows*cols),sizeof(int));
east=(int*)calloc((rows*cols),sizeof(int));
west=(int*)calloc((rows*cols),sizeof(int));
gm=(int*)calloc((rows*cols*4),sizeof(int));		// Gradient magnitude

if(north == NULL || south == NULL || east == NULL || west == NULL || gm == NULL)
   {mexErrMsgTxt("Error allocating the differences matrix.\n");
    exit(1);
    }

// Calculate the neighbor differences and store in 4 matricies
for (i = 1; i < rows; i++)
    for (j = 0; j < cols; j++)
         {index = i * cols + j;
	  north[index] = (int)(smooth[index - cols] - smooth[index]);
	  }
for (i = 0; i < rows-1; i++)
    for (j = 0; j < cols; j++)
         {index = i * cols + j;
	  south[index] = (int)(smooth[index + cols] - smooth[index]);
	  }
for (i = 0; i < rows; i++)
    for (j = 0; j < cols-1; j++)
         {index = i * cols + j;
	  east[index] = (int)(smooth[index + 1] - smooth[index]);
	  }
for (i = 0; i < rows; i++)
    for (j = 1; j < cols; j++)
         {index = i * cols + j;
	  west[index] = (int)(smooth[index - 1] - smooth[index]);
	  }

// Calcuate magnitude of the gradient
for (i = 0; i < rows*cols; i++)
    {
    gm[i * 4] = abs(north[i]);
    gm[i * 4 + 1] = abs(south[i]);
    gm[i * 4 + 2] = abs(east[i]);
    gm[i * 4 + 3] = abs(west[i]);
    }

// Find the median gradient across the entire image
median = find_median(gm, rows*cols*4);

// Normalize the neighbor differences w.r.t. the median
for(i = 1; i < rows-1; i++)
    for(j = 1; j < cols-1; j++)    
        {index = i * cols + j;
        north[index] -= median;
        south[index] -= median;
        east[index]  -= median;
        west[index]  -= median;
        }

// Recompute the gradient w/ the normalized differences
for (i = 0; i < rows*cols; i++)
    {
    gm[i * 4] = abs(north[i]);
    gm[i * 4 + 1] = abs(south[i]);
    gm[i * 4 + 2] = abs(east[i]);
    gm[i * 4 + 3] = abs(west[i]);
    }

median = find_median(gm, rows*cols*4);

free(north);
free(south);
free(east);
free(west);
free(gm);

return median;
}


// Modified quicksort to find median (see Sedgewick, pp 128)
int find_median(int *grad, int N)
{
int left, right, i, j, k;
int v, temp, median;
int *list;

if((list = (int*)calloc(N+1, sizeof(double)))==NULL)
    {mexErrMsgTxt("Error allocating the smooth image.\n");
     exit(1);
     }

// Copy the gradient matrix to temp storage for partial sort
for (i = 1; i <= N; i++)
    list[i] = grad[i-1];
    
left = 1;
right = N;
k = N/2;

while(right > left)
    {v = list[right];
     i = left - 1;
     j = right;
     for(;;)
         {while (list[++i] < v);
	  while (list[--j] > v);
	  if (i >= j) break;
	  temp = list[i];
	  list[i] = list[j];
	  list[j] = temp;
	  }

     temp = list[i];
     list[i] = list[right];
     list[right] = temp;
     if (i >= k)
         right = i - 1;
     if (i <= k)
         left = i + 1;
     }
median = list[k];
free(list);
return median;
}


/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void non_max_supp(double *mag, double *gradx, double *grady, int nrows, int ncols,
    unsigned char *result) 
{
    int rowcount, colcount,count;
    double *magrowptr,*magptr;
    double *gxrowptr,*gxptr;
    double *gyrowptr,*gyptr, z1,z2;
    short m00,gx,gy;
    double mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;
    

   /****************************************************************************
   * Zero the edges of the result image.
   ****************************************************************************/
    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1); 
        count<ncols; resultptr++,resultrowptr++,count++){
        *resultrowptr = *resultptr = (unsigned char) 0;
    }

    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
        *resultptr = *resultrowptr = (unsigned char) 0;
    }

   /****************************************************************************
   * Suppress non-maximum points.
   ****************************************************************************/
   for(rowcount=1, magrowptr=mag+ncols+1, gxrowptr=gradx+ncols+1,
      gyrowptr=grady+ncols+1, resultrowptr=result+ncols+1;
      rowcount<nrows-2; 
      rowcount++, magrowptr+=ncols, gyrowptr+=ncols, gxrowptr+=ncols,
      resultrowptr+=ncols){   
      for(colcount=1, magptr=magrowptr, gxptr=gxrowptr, gyptr=gyrowptr,
         resultptr=resultrowptr; colcount<ncols-2; 
         colcount++, magptr++, gxptr++, gyptr++, resultptr++){   
         m00 = *magptr;
         if(m00 == 0){
            *resultptr = (unsigned char) NOEDGE;
         }
         else{
            xperp = -(gx = *gxptr)/((double)m00);
            yperp = (gy = *gyptr)/((double)m00);
         }

         if(gx >= 0){
            if(gy >= 0){
                    if (gx >= gy)
                    {  
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                        
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {    
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp; 
                    }
                }
                else
                {
                    if (gx >= -gy)
                    {
                        /* 101 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;
            
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {    
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp; 
                    }
                }
            }
            else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {          
                        /* 011 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (-gx > -gy)
                    {
                        /* 001 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            } 

            /* Now determine if the current point is a maximum point */

            if ((mag1 > 0.0) || (mag2 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }
            else
            {    
                if (mag2 == 0.0)
                    *resultptr = (unsigned char) NOEDGE;
                else
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
            }
        } 
    }
}


/*******************************************************************************
* PROCEDURE: magnitude_x_y
* PURPOSE: Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void magnitude_x_y(double *delta_x, double *delta_y, int rows, int cols,
        double *magnitude)
{
   int r, c, pos;
   double sq1, sq2;

   /****************************************************************************
   * Allocate an image to store the magnitude of the gradient.
   ****************************************************************************/

   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         sq1 = delta_x[pos] * delta_x[pos];
         sq2 = delta_y[pos] * delta_y[pos];
         magnitude[pos] = (double)(0.5 + sqrt(sq1 + sq2));
      }
   }

}

/*******************************************************************************
* PROCEDURE: derivative_x_y
* PURPOSE: Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void derivative_x_y(double *smoothedim, int rows, int cols,
        double *delta_x, double *delta_y)
{
   int r, c, pos;

   /****************************************************************************
   * Compute the x-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      pos = r * cols;
      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      for(c=1;c<(cols-1);c++,pos++){
         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      }
      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
   }

   /****************************************************************************
   * Compute the y-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   for(c=0;c<cols;c++){
      pos = c;
      delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
      pos += cols;
      for(r=1;r<(rows-1);r++,pos+=cols){
         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
      }
      delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
   }
}
