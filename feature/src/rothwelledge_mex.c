/*******************************************************************************
* Purpose: This program contains an implementation of the Edge Detector
* described in the technical report: "Driving Vision by Topology" by Charlie
* Rothwell, Joe Mundy, Bill Hoffman and Van-Duc Nguyen. A shorter version of
* the report was published as a paper with the same name.
* This implementation of the program was also aided by the code Charlie Rothwell
* sent us. He sent us the edge detector modulules of a larger vision package
* that was coded in C++. This program just pulled together those modules into
* a complete edge detector program coded in 'C'. I tried to use as much of his
* code as possible to both simplify the coding task and to make sure that this
* implementation is consistent with theirs.
* Name: Mike Heath
* Date: 5/2/96
*
* To Compile: gcc -O2 -o rothwell rothwell.c -lm
*
* NOTE: ALL OF THE 2-DIMENSIONAL ARRAYS IN THIS PROGRAM ARE STORED TRANSPOSED
* FROM THE USUAL WAY ARRAYS ARE STORED. I HAVE, WITH SOME PAIN, FOLLED THIS
* METHOD.
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#define GAUSS_TAIL 0.015       /* As recommended in the technical report. */
#define FAR 65535
#define DUMMYTHETA 10000.0
#define DUMMYSIGMA 0.0
#define VERBOSE 0
typedef struct{
    int x;
    int y;
    double thin;
}XYFLOAT;

void Sub_pixel_interpolation(double **_grad, double **_dx, double **_dy, int cols,
        int rows, int kwidth, double _low, double ALPHA, double ***_thresh, int ***_dist,
        double ***_theta);
void Thicken_threshold(double **_thresh, int **_dist, int x, int y, double _low,
        int kwidth);
void Compute_gradient(double **dx, double **dy, int cols, int rows, int kwidth,
        double ***grad);
void Compute_x_gradient(double **smoothedimage, int cols, int rows, int kwidth,
        double ***dx);
void Compute_y_gradient(double **smoothedimage, int cols, int rows, int kwidth,
        double ***dy);
void Smooth_image(double **image, int cols, int rows,
        double **smoothed, double sigma, int *kwidth, double gauss_tail);
void Set_kernel(double **kernel, double sigma, int width, int k_size);
int **Make_int_image(int x, int y);
void Set_int_image(int **image, int val, int cols, int rows);
void Copy_int_image(int **image1, int **image2, int cols, int rows);
void Free_int_image(int ***ptr);
double **Make_double_image(int x, int y);
void Set_double_image(double **image, double val, int cols, int rows);
void Copy_double_image(double **image1, double **image2, int cols, int rows);
void Free_double_image(double ***ptr);
void Set_thresholds(int **_dist, double **_grad, double **_thresh, double ***_thin,
        int cols, int rows, double _low);
void Forward_chamfer(int m, int n, int **dist, double **param);
void Backward_chamfer(int m, int n, int **dist, double **param);
void Alt1_chamfer(int m, int n, int **dist, double **param);
void Alt2_chamfer(int m, int n, int **dist, double **param);
int Minimum4(int a, int b, int c, int d);
int Minimum5(int a, int b, int c, int d, int e);
int compare(const void *xyf1, const void *xyf2);
void Thin_edges(double **_thin, double **_thresh, int cols, int rows, int kwidth);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int rows, cols, x, y, pos, kwidth, **dist=NULL;
    double sigma=DUMMYSIGMA, low, alpha, **image=NULL, **smoothed=NULL, **dx=NULL,
            **dy=NULL, **grad=NULL, **thresh=NULL, **theta=NULL, **thin=NULL;
    double *input, *_input=NULL, *outthe;
    unsigned char *output;
    
    /* check for the proper number of arguments */
    if(nrhs != 5)
        mexErrMsgTxt("Five inputs required.");
    if(nlhs > 3)
        mexErrMsgTxt("Too many output arguments.");
    
    // inputs
    
    if( !mxIsNumeric(prhs[0]))
        mexErrMsgTxt("Input image must be a matrix.");
    else 
        input = mxGetPr(prhs[0]);
    
    if( !mxIsEmpty(prhs[1])) 
        if( !mxIsNumeric(prhs[1]))
            mexErrMsgTxt("Input image must be a matrix.");
        else {
            _input = mxGetPr(prhs[1]);
            kwidth = 1;
        }
    
    if( mxIsEmpty(prhs[1]) ) {
        if(!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 ) 
            mexErrMsgTxt("Input sigma must be a scalar.");
        else
            sigma = mxGetScalar(prhs[2]);
    } else if(!mxIsEmpty(prhs[2])) 
        mexErrMsgTxt("Input sigma must be empty when the partial derivatives are passed.");

    if( !mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1 ) 
        mexErrMsgTxt("Input low must be a scalar.");
    else
        low = mxGetScalar(prhs[3]);
    
    if( !mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4])!=1 ) 
        mexErrMsgTxt("Input alpha must be a scalar.");
    else
        alpha = mxGetScalar(prhs[4]);
    
    /*  create a pointer to the input matrix */
    input = mxGetPr(prhs[0]);
    /*  get the dimensions of the matrix input */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);
    
    // principal output
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericMatrix(rows, cols, mxUINT8_CLASS, mxREAL);;
    /*  create a C pointer to a copy of the output matrix */
    output = mxGetData(plhs[0]);
    
    // Process
    if (mxIsEmpty(prhs[1])) {
        image = Make_double_image(cols, rows);
        for(y=0, pos=0;y<rows;y++)
            for(x=0;x<cols;x++, pos++)
                image[x][y] = input[pos];
        
        smoothed = Make_double_image(cols, rows);
        
        Smooth_image(image, cols, rows, smoothed, sigma, &kwidth, GAUSS_TAIL);
        
        Compute_x_gradient(smoothed, cols, rows, kwidth, &dx);
        Compute_y_gradient(smoothed, cols, rows, kwidth, &dy);
        Free_double_image(&smoothed);

    } else {
        dx = Make_double_image(cols, rows);
        dy = Make_double_image(cols, rows);
        for(y=0, pos=0;y<rows;y++)
            for(x=0;x<cols;x++, pos++) {
            dx[x][y] = input[pos];
            dy[x][y] = _input[pos];
            }
    }
    
    /* compute the magnitude of the gradient */
    //grad = Make_double_image(cols, rows);
    Compute_gradient(dx, dy, cols, rows, kwidth, &grad);
        
   Sub_pixel_interpolation(grad, dx, dy, cols, rows, kwidth, low, alpha,
            &thresh, &dist, &theta);
    
    Free_double_image(&dx);
    Free_double_image(&dy);
    
    Set_thresholds(dist, grad, thresh, &thin, cols, rows, low);
    
  Thin_edges(thin, thresh, cols, rows, kwidth);
    
    /****************************************************************************
     * Copy the thinned image from the 2-Dimensional array back into the original
     * image array. This will overwrite it but that is o.k.
     ****************************************************************************/
    
    /* edge map */
    for(y=0, pos=0;y<rows;y++)
        for(x=0;x<cols;x++, pos++)
            if(thin[x][y] != 0.0) output[pos] = (unsigned char)255;
            else output[pos] = (unsigned char)0;
    
    /* other outputs */
    if(nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        outthe = mxGetPr(plhs[1]);
        for(y=0, pos=0;y<rows;y++)
            for(x=0;x<cols;x++, pos++) {
            if (theta[x][y] == DUMMYTHETA)
                outthe[pos] = 0;
            else
                outthe[pos] = M_PI * theta[x][y] / 180.;
            }
    }
    
    if(nlhs==3) {
        plhs[2] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        memcpy(mxGetPr(plhs[2]), grad, rows*cols*sizeof(double));
    }
    
    
    Free_double_image(&grad);
    Free_double_image(&theta);
    Free_int_image(&dist);

}

/*******************************************************************************
* Method to thin the image using the variation of Tsai-Fu thinning used
* by Van-Duc Nguyen in Geo-Calc. This relies on computing the genus of
* an edge location, and removing it if it is not a dangling chain as has
* genus zero. We also order the edges by strength and try to remove the weaker
* ones first. This accounts for non-maximal supression, and does it in a
* topology preserving way. Note that we are creating a CoolList with a large
* number of elements, and then sorting it -- this is likely to be quite slow.
* An alternative implementation would be better.
*******************************************************************************/
void Thin_edges(double **_thin, double **_thresh, int cols, int rows, int kwidth)
{
   /* Find all of the edgels with a strength > _low */
   int x, y, a, b, c, d, e, f, g, h, genus, count;
   XYFLOAT *edgel_array=NULL;
   int edgel_array_len = 0;
   int pos = 0;

   if((edgel_array = (XYFLOAT *) calloc(cols*rows, sizeof(XYFLOAT)))==NULL){
      fprintf(stderr, "Error allocating the xydouble array in Thin_edges().\n");
      exit(1);
   }

   count = 1;      /* count set to dummy value */
   while(count){   /* Thin until no Pixels are removed */
      count = 0;
      edgel_array_len = 0;
      for(x=kwidth;x<(cols-kwidth);x++){
         for(y=kwidth;y<(rows-kwidth);y++){
            if(_thin[x][y] > _thresh[x][y]){
               edgel_array[edgel_array_len].x = x;
               edgel_array[edgel_array_len].y = y;
               edgel_array[edgel_array_len].thin = _thin[x][y];
               edgel_array_len++;
            }
         }
      }

      /*************************************************************************
      * Now sort the list; this could be slow if we have a lot of potential.
      * edges - surely we have to do number of elements (not -1)?
      *      qsort(edgel_array, edgel_array_len-1, sizeof(xydouble), &compare);
      *************************************************************************/

      qsort(edgel_array, edgel_array_len, sizeof(XYFLOAT), compare);

      /* To assist in setting the thresholds: */
      if(VERBOSE)
         printf("Edgel strengths range from %f to %f\n", edgel_array[0].thin,
            edgel_array[edgel_array_len-1].thin);
     

      /************************************************************************* 
      * Do the thinning taking the weakest edges first and works
      * up through the list strengthwise.
      **************************************************************************/
      for(pos=0;pos<edgel_array_len;pos++){
         x = edgel_array[pos].x;
         y = edgel_array[pos].y;
         
         if( _thin[x-1][y-1] > _thresh[x-1][y-1] )  a = 1; else a = 0;
         if( _thin[x][y-1]   > _thresh[x][y-1]   )  b = 1; else b = 0;
         if( _thin[x+1][y-1] > _thresh[x+1][y-1] )  c = 1; else c = 0;
         if( _thin[x+1][y]   > _thresh[x+1][y]   )  d = 1; else d = 0;
         if( _thin[x+1][y+1] > _thresh[x+1][y+1] )  e = 1; else e = 0;
         if( _thin[x][y+1]   > _thresh[x][y+1]   )  f = 1; else f = 0;
         if( _thin[x-1][y+1] > _thresh[x-1][y+1] )  g = 1; else g = 0;
         if( _thin[x-1][y]   > _thresh[x-1][y]   )  h = 1; else h = 0;
         
         genus = a+b+c+d+e+f+g+h;

         /* Continue if the pixel is not dangling. */

         if((genus!=1) && (genus!=8)){

            genus += h*a*b+b*c*d+d*e*f+f*g*h-a*b-b*c-c*d-d*e-e*f-f*g -
                     g*h-h*a-h*b-b*d-d*f-f*h-1;

            /* If the genus is zero delete the edge */

            if(genus == 0){
               count++;
               _thin[x][y] = 0.0;
            }
         }
      }
   }

   free(edgel_array);
}

/*******************************************************************************
* This is the compare function for the quicksort function provided with 'C'.
*******************************************************************************/
int compare(const void *xyf1, const void *xyf2)
{
  const XYFLOAT *f1 = (const XYFLOAT *)xyf1;
  const XYFLOAT *f2 = (const XYFLOAT *)xyf2;
  if (f1->thin < f2->thin) return(-1);
  if (f1->thin == f2->thin) return(0);
  return 1;
}

/*******************************************************************************
* Takes the _thresh image that contains threshold values near to where
* non-maximal suppression succeeded, and zero elsewhere, and extend the
* values to all areas of the image. This is done using chamfer masks so that
* the final threshold assigned at any one point (ie. a point that was initially
* zero) is functionally dependent on the the strengths of the nearest good
* edges. At present we linearly interpolate between the two (approximately)
* closest edges.
* 
* Try to do the same process using Delauney triangulation (CAR, March 1995), in
* an attempt to image the efficiency from a memeory management point of view.
* However, the triangulation becomes so complex that the computation time
* becomes incredibably long. Therefore putting up with the Chamfer method for
* the moment.
* 
* The histogram calculation was added to support
* Edgel change detection-JLM May 1995
*******************************************************************************/
void Set_thresholds(int **_dist, double **_grad, double **_thresh, double ***_thin,
   int cols, int rows, double _low)
{
   int **fdist=NULL, **bdist=NULL, **a1dist=NULL, **a2dist=NULL;
   double **fth=NULL, **bth=NULL, **a1th=NULL, **a2th=NULL;
   int x, y, option;
   double num, den;
   double max_gradient = _low;

   *_thin = Make_double_image(cols, rows);

   fdist = Make_int_image(cols, rows);
   bdist = Make_int_image(cols, rows);
   a1dist = Make_int_image(cols, rows);
   a2dist = Make_int_image(cols, rows);
   Copy_int_image(_dist, fdist, cols, rows);
   Copy_int_image(_dist, bdist, cols, rows);
   Copy_int_image(_dist, a1dist, cols, rows);
   Copy_int_image(_dist, a2dist, cols, rows);

   fth = Make_double_image(cols, rows);
   bth = Make_double_image(cols, rows);
   a1th = Make_double_image(cols, rows);
   a2th = Make_double_image(cols, rows);
   Copy_double_image(_thresh, fth, cols, rows);
   Copy_double_image(_thresh, bth, cols, rows);
   Copy_double_image(_thresh, a1th, cols, rows);
   Copy_double_image(_thresh, a2th, cols, rows);

   Forward_chamfer(cols, rows, fdist, fth);
   Backward_chamfer(cols, rows, bdist, bth);
   Alt1_chamfer(cols, rows, a1dist, a1th);
   Alt2_chamfer(cols, rows, a2dist, a2th);

   /****************************************************************************
   * The range of the effect of the smoothing kernel, including the scale
   * factor we have ignored up to now for the chamfer masks
   *   int range = 3*_width;
   ****************************************************************************/
   for(x=0;x<cols;x++){
      for(y=0;y<rows;y++){
         if(_thresh[x][y] == _low){

            /* Determine the two closest edge points. */

            option = Minimum4(fdist[x][y],bdist[x][y],a1dist[x][y],a2dist[x][y]);
            switch(option){
               case 1:
               case 2:
                  den = (fdist[x][y]+bdist[x][y]);
                  num = (bdist[x][y]*fth[x][y]+fdist[x][y]*bth[x][y]);
                  break;

               case 3:
               case 4:
                  den = (a1dist[x][y]+a2dist[x][y]);
                  num = (a2dist[x][y]*a1th[x][y]+a1dist[x][y]*a2th[x][y]);
                  break;

               default:
                  den = num = 1.0; /* Dummy values */
                  break;
            }
            if(den != 0.0)
               _thresh[x][y] = num / den;
            else if(_thresh[x][y] <= _low) _thresh[x][y] = _low;
         }

         if(_grad[x][y] > _thresh[x][y]){
            if(_grad[x][y] > max_gradient) max_gradient = _grad[x][y];
            (*_thin)[x][y] = _grad[x][y];
         }
      }
   }

   /*****************************************************************************
   * I don't see this histogram being used for anything in the program so I think
   * that I can just avoid it by commenting it out.   - Mike Heath
   *  if(_gradient_histogram){
   *     _ghist = new Histogram(_histogram_resolution, _low, max_gradient);
   *     for(x=0;x<_image->GetSizeX();x++)
   *        for(y=0;y<_image->GetSizeY();y++)
   *           _ghist->UpCount(_grad[x][y]);
   *  }
   ****************************************************************************/

   Free_int_image(&fdist);
   Free_int_image(&bdist);
   Free_int_image(&a1dist);
   Free_int_image(&a2dist);
   Free_double_image(&fth);
   Free_double_image(&bth);
   Free_double_image(&a1th);
   Free_double_image(&a2th);
}

/*******************************************************************************
* Performs a forward chamfer convolution on the dist image and associates
* a send image (param) that reports on some parameter of the nearest pixel.
* The image sizes are mxn
*******************************************************************************/
void Forward_chamfer(int m, int n, int **dist, double **param)
{
   int i, j, val;

   for(i=1;i<(m-1);i++){
      for(j=1;j<(n-1);j++){

         val = Minimum5(dist[i-1][j-1]+4, dist[i-1][j]+3, dist[i-1][j+1]+4,
                        dist[i][j-1]+3, dist[i][j]);

         switch(val){
            case 1:
               dist[i][j] = dist[i-1][j-1]+4;
               param[i][j] = param[i-1][j-1];
               break;

            case 2:
               dist[i][j] = dist[i-1][j]+3;
               param[i][j] = param[i-1][j];
               break;

            case 3:
               dist[i][j] = dist[i-1][j+1]+4;
               param[i][j] = param[i-1][j+1];
               break;

            case 4:
               dist[i][j] = dist[i][j-1]+3;
               param[i][j] = param[i][j-1];
               break;

            case 5:
               break;
         }
      }
   }
}

/*******************************************************************************
* Performs a backward chamfer convolution on the dist and param images.
*******************************************************************************/
void Backward_chamfer(int m, int n, int **dist, double **param)
{
   int i,j,val;

   for(i=m-2;i>0;i--){
      for(j=n-2;j>0;j--){

         val = Minimum5(dist[i][j], dist[i][j+1]+3, dist[i+1][j-1]+4,
                        dist[i+1][j]+3, dist[i+1][j+1]+4 );

         switch(val){
            case 1:
               break;

            case 2:
               dist[i][j] = dist[i][j+1]+3;
               param[i][j] = param[i][j+1];
               break;

            case 3:
               dist[i][j] = dist[i+1][j-1]+4;
               param[i][j] = param[i+1][j-1];
               break;

            case 4:
               dist[i][j] = dist[i+1][j]+3;
               param[i][j] = param[i+1][j];
               break;

            case 5:
               dist[i][j] = dist[i+1][j+1]+4;
               param[i][j] = param[i+1][j+1];
               break;
         }
      }
   }
}

/*******************************************************************************
* Performs a chamfer convolution starting from (minx,maxy) on the dist image
* and associates a send image (param) that reports on some parameter of the
* nearest pixel. The image sizes are mxn
********************************************************************************/
void Alt1_chamfer(int m, int n, int **dist, double **param)
{
   int i,j,val;

   for(i=1;i<m-1;i++){
      for(j=n-2;j>0;j--){

         val = Minimum5(dist[i-1][j+1]+4, dist[i-1][j]+3, dist[i-1][j-1]+4,
                        dist[i][j+1]+3, dist[i][j]);

         switch (val){

            case 1:
               dist[i][j] = dist[i-1][j+1]+4;
               param[i][j] = param[i-1][j+1];
               break;

            case 2:
               dist[i][j] = dist[i-1][j]+3;
               param[i][j] = param[i-1][j];
               break;

            case 3:
               dist[i][j] = dist[i-1][j-1]+4;
               param[i][j] = param[i-1][j-1];
               break;

            case 4:
               dist[i][j] = dist[i][j+1]+3;
               param[i][j] = param[i][j+1];
               break;

            case 5:
               break;
         }
      }
   }
}

/*******************************************************************************
* Performs a chamfer convolution starting from (maxx,miny) on the dist image
* and associates a send image (param) that reports on some parameter of the
* nearest pixel. The image sizes are mxn
*******************************************************************************/
void Alt2_chamfer(int m, int n, int **dist, double **param)
{
   int i,j,val;

   for(i=m-2;i>0;i--){
      for(j=1;j<n-1;j++){

         val = Minimum5(dist[i][j], dist[i][j+1]+3, dist[i+1][j-1]+4,
                        dist[i+1][j]+3, dist[i+1][j+1]+4);

         switch (val){

            case 1:
               break;

            case 2:
               dist[i][j] = dist[i][j+1]+3;
               param[i][j] = param[i][j+1];
               break;

            case 3:
               dist[i][j] = dist[i+1][j-1]+4;
               param[i][j] = param[i+1][j-1];
               break;

            case 4:
               dist[i][j] = dist[i+1][j]+3;
               param[i][j] = param[i+1][j];
               break;

            case 5:
               dist[i][j] = dist[i+1][j+1]+4;
               param[i][j] = param[i+1][j+1];
               break;
         }
      }
   }
}

/*******************************************************************************
* Determines the minimum of four ints.
*******************************************************************************/
int Minimum4(int a, int b, int c, int d)
{
    if((a<=b) && (a<=c) && (a<=d)) return(1);
    else if((b<=c) && (b<=d)) return(2);
    else if((c<=d)) return(3);
    else return(4);
}

/*******************************************************************************
* Determines the minimum of five ints.
*******************************************************************************/
int Minimum5(int a, int b, int c, int d, int e)
{
   if((a<=b) && (a<=c) && (a<=d) && (a<=e)) return(1);
    else if((b<=c) && (b<=d) && (b<=e)) return(2);
    else if((c<=d) && (c<=e)) return(3);
    else if(d<=e) return(4);
    else return(5);
}

/*******************************************************************************
* A procedure that performs sub-pixel interpolation for all edges greater than
* the threshold by parabolic fitting. Writes edges into the _thresh image if they
* are maxima and above _low. This gives a good indication of the local edge
* strengths. Stores sub-pixel positions in _dx and _dy, and set the orientations
* in _theta.
*******************************************************************************/
void Sub_pixel_interpolation(double **_grad, double **_dx, double **_dy, int cols,
   int rows, int kwidth, double _low, double ALPHA, double ***_thresh,
   int ***_dist, double ***_theta)
{
   double *g0=NULL, *g1=NULL, *g2=NULL, *dx=NULL, *dy=NULL;
   double h1,h2;
   double k = 180.0/M_PI;
   int x, y, orient;
   double theta, grad;
   double fraction, dnewx, dnewy;

   *_thresh = Make_double_image(cols, rows);
   Set_double_image((*_thresh), _low, cols, rows);

   *_dist = Make_int_image(cols, rows);
   Set_int_image((*_dist), FAR, cols, rows);

   *_theta = Make_double_image(cols, rows);
   Set_double_image((*_theta), DUMMYTHETA, cols, rows);

   /* Add 1 to get rid of border effects. */
   for(x=(kwidth+1);x<(cols-kwidth-1);x++){
      g0 = _grad[x-1];  g1 = _grad[x];  g2 = _grad[x+1];
      dx = _dx[x];      dy = _dy[x];

      for(y=(kwidth+1);y<(rows-kwidth-1);y++){
         /* First check that we have a potential edge. */
         if(g1[y] > _low){
            theta = k*atan2(dy[y],dx[y]);

            /* Now work out which direction wrt the eight-way */
            /* neighbours the edge normal points */
            if(theta >= 0.0) orient = (int)(theta/45.0);
            else orient = (int)(theta/45.0+4);

            /* if theta == 180.0 we will have orient = 4 */
            orient = orient%4;

            /* And now compute the interpolated heights */
            switch(orient){
               case 0:
                  grad = dy[y]/dx[y];
                  h1 = grad*g0[y-1] + (1 - grad)*g0[y];
                  h2 = grad*g2[y+1] + (1 - grad)*g2[y];
                  break;

               case 1:
                  grad = dx[y]/dy[y];
                  h1 = grad*g0[y-1] + (1 - grad)*g1[y-1];
                  h2 = grad*g2[y+1] + (1 - grad)*g1[y+1];
                  break;

               case 2:
                  grad = -dx[y]/dy[y];
                  h1 = grad*g2[y-1] + (1 - grad)*g1[y-1];
                  h2 = grad*g0[y+1] + (1 - grad)*g1[y+1];
                  break;

               case 3:
                  grad = -dy[y]/dx[y];
                  h1 = grad*g2[y-1] + (1 - grad)*g2[y];
                  h2 = grad*g0[y+1] + (1 - grad)*g0[y];
                  break;

               default:
                  h1 = h2 = 0.0;  /* Dummy value; */
                  fprintf(stderr, "*** ERROR ON SWITCH IN NMS ***\n");
            }

            /* Do subpixel interpolation by fitting a parabola */
            /* along the NMS line and finding its peak */

            fraction = (h1-h2)/(2.0*(h1-2.0*g1[y]+h2));
            switch(orient){
               case 0:
                  dnewx = fraction;
                  dnewy = dy[y]/dx[y]*fraction;
                  break;

               case 1:
                  dnewx = dx[y]/dy[y]*fraction;
                  dnewy = fraction;
                  break;

               case 2:
                  dnewx = dx[y]/dy[y]*fraction;
                  dnewy = fraction;
                  break;

               case 3:
                  dnewx = - fraction;
                  dnewy = - dy[y]/dx[y]*fraction;
                  break;

               default:
                  dnewx = dnewy = 0.0; /* Dummy values */
                  fprintf(stderr, "*** ERROR ON SWITCH IN NMS ***\n");
            }

            /*******************************************************************
            * Now store the edge data, re-use _dx[][] and _dy[][]
            * for sub-pixel locations (don't worry about the junk
            * that is already in them). Use any edgels that get
            * non-maximal suppression to bootstrap the image
            * thresholds. The >= is used rather than > for reasons
            * involving non-generic images. Should this be interpolated
            * height -- = g1[y] + frac*(h2-h1)/4 ?
            *******************************************************************/
            if((g1[y]>=h1)&&(g1[y]>=h2)&&(fabs(dnewx)<=0.5)&&(fabs(dnewy)<=0.5)){
               if(g1[y]*ALPHA > _low)
                  (*_thresh)[x][y] = ALPHA * g1[y]; /* Use a ALPHA% bound */

		  
               /* _thresh image starts off as being equal to _low */
               /* else                                            */
               /*   _thresh[x][y] = _low;                         */

               Thicken_threshold((*_thresh), (*_dist), x, y, _low, kwidth);
            }
                       
            /* + 0.5 is to account for targetjr display offset */

            if((fabs(dnewx)<=0.5) && (fabs(dnewy)<=0.5)){
               dx[y] = x + dnewx + 0.5;
	       dy[y] = y + dnewy + 0.5;
            }
            else{
               dx[y] = x + 0.5;
	       dy[y] = y + 0.5;
            }

            (*_theta)[x][y] = theta;
         }

         /* For consistency assign these values even though the */
         /* edge is below strength.                             */

         else{
            dx[y] = x + 0.5;
            dy[y] = y + 0.5;
         }
      }
   }

   /****************************************************************************
   * Clean up around the border to ensure consistency in the _dx and _dy values.
   ****************************************************************************/
   for(x=0;x<cols;x++){
      for(y=0;y<=kwidth;y++){
         _dx[x][y] = x + 0.5;
         _dy[x][y] = y + 0.5;
      }
      for(y=(rows-kwidth-1);y<rows;y++){
         _dx[x][y] = x + 0.5;
         _dy[x][y] = y + 0.5;
      }
   }

   for(y=(kwidth+1);y<(rows-kwidth-1);y++){
      for(x=0;x<=kwidth;x++){
         _dx[x][y] = x + 0.5;
         _dy[x][y] = y + 0.5;
      }
      for(x=(cols-kwidth-1);x<cols;x++){
         _dx[x][y] = x + 0.5;
         _dy[x][y] = y + 0.5;
      }
   }
}


/*******************************************************************************
* Thickens the threshold image around each good pixel to take account for the
* smoothing kernel (almost a dilation with a square structuring element).
*******************************************************************************/
void Thicken_threshold(double **_thresh, int **_dist, int x, int y, double _low,
   int kwidth)
{
   int i,j;

   /* Experimental change 13/4/95 by CAR */
   /* int width = _width; Changed back because not mentioned in the paper MH */
   int width = 0;

   for(i=(x-width);i<=(x+width);i++){
      for(j=(y-width);j<=(y+width);j++){
         _dist[i][j] = 0;
         if(_thresh[i][j] != _low){
            if(_thresh[x][y] < _thresh[i][j]) _thresh[i][j] = _thresh[x][y];
         }
         else _thresh[i][j] = _thresh[x][y];
      }
   }
}

/*******************************************************************************
* Compute the absolute intensity surface gradient, _grad[][].
*******************************************************************************/
void Compute_gradient(double **dx, double **dy, int cols, int rows, int kwidth,
   double ***grad)
{
    int x, y;

   *grad = Make_double_image(cols, rows);

   /****************************************************************************
   * JLM limits here are _width-1 because of _ksize being 2*_width + 1.
   * I don't understand why to use _width-1 but I will go along with it.
   *     - Mike Heath
   ****************************************************************************/

   for(x=kwidth;x<(cols-kwidth-1);x++){
      for(y=kwidth;y<(rows-kwidth-1);y++){
         (*grad)[x][y] = sqrt(dx[x][y]*dx[x][y] + dy[x][y]*dy[x][y]);
      }
   }
}

/*******************************************************************************
* Convolves with the kernel in the x direction, to compute the local derivative
* in that direction
*******************************************************************************/
void Compute_x_gradient(double **smoothedimage, int cols, int rows, int kwidth,
   double ***dx)
{
   int x, y;

   *dx = Make_double_image(cols, rows);

   for(y=(kwidth+1);y<(rows-kwidth-1);y++){
      for(x=(kwidth+1);x<(cols-kwidth-1);x++){
         (*dx)[x][y] = smoothedimage[x+1][y] - smoothedimage[x-1][y];
      }
   }
}

/*******************************************************************************
* Convolves the original image with the kernel in the y direction to give the
* local y derivative.
*******************************************************************************/
void Compute_y_gradient(double **smoothedimage, int cols, int rows, int kwidth,
   double ***dy)
{
    int x,y;

   *dy = Make_double_image(cols, rows);

   for(x=(kwidth+1);x<(cols-kwidth-1);x++){
      for(y=(kwidth+1);y<(rows-kwidth-1);y++){
         (*dy)[x][y] = smoothedimage[x][y+1] - smoothedimage[x][y-1];
      }
   }
}

/*******************************************************************************
* Convolves the image with the smoothing kernel.
*******************************************************************************/
void Smooth_image(double **image, int cols, int rows,
   double **smoothedimage, double sigma, int *kwidth, double gauss_tail)
{
   int width = (int)(sigma*sqrt(2*log(1/gauss_tail))+1);
   int k_size = 2*width+ 1;
   double *kernel=NULL;
   double **tmp=NULL;
   int x,y,xx,yy,i;

   Set_kernel(&kernel, sigma, width, k_size);
   *kwidth = width;

   tmp = Make_double_image(cols, rows);
  
   /****************************************************************************
   * x direction
   ****************************************************************************/
   for(y=0;y<rows;y++){
      for(x=width;x<(cols-width);x++){
         for(i=0,xx=(x-width);i<k_size;i++,xx++)
            tmp[x][y] += image[xx][y]*kernel[i];
      }
   }

   /****************************************************************************
   * y direction
   ****************************************************************************/
   for(y=width;y<(rows-width);y++){
      for(x=0;x<cols;x++){
         for(i=0,yy=(y-width);i<k_size;i++,yy++)
            smoothedimage[x][y] += tmp[x][yy]*kernel[i];
      }
   }

   Free_double_image(&tmp);
   free(kernel);
}


/*******************************************************************************
* Sets up the Gaussian convolution kernel.
*******************************************************************************/
void Set_kernel(double **kernel, double sigma, int width, int k_size)
{
   int i,x;
   double s2 = 2.0*sigma*sigma;
   double det = sigma*sqrt(2.0*M_PI);
   // double n = 0.;
   
   if(((*kernel) = (double *) calloc(k_size, sizeof(double))) == NULL){
      fprintf(stderr, "Error allocating the smoothing filter array.\n");
      exit(1);
   }

   for(i=0,x=(-width);i<k_size;i++,x++) {
       (*kernel)[i] = exp(-x*x/s2)/det;
       // n += (*kernel)[i];
   }
  // for(i=0,x=(-width);i<k_size;i++,x++) (*kernel)[i] /= n;
   
}


/*******************************************************************************
* Returns an m*n array of ints
*******************************************************************************/
int **Make_int_image(int x, int y)
{
   int **image;
   int i;

   if((image = (int **) calloc(x, sizeof(int *))) == NULL){
      fprintf(stderr, "Error allocating an array in Make_int_image().\n");
      exit(1);
   }
   if((image[0] = (int *) calloc(x*y, sizeof(int))) == NULL){
      fprintf(stderr, "Error allocating an array in Make_int_image().\n");
      exit(1);
   }
   for(i=0;i<x;i++) image[i] = image[0] + (y * i);

   return(image);
}

/*******************************************************************************
* Sets an int image to val.
*******************************************************************************/
void Set_int_image(int **image, int val, int cols, int rows)
{
   int x, y;
   int *ptr = image[0];

   /* copy first col */
   for(y=0;y<rows;y++) ptr[y] = val;

   for(x=1;x<cols;x++){
      ptr = image[x];
      memcpy((char*)ptr, (char*)image[x-1], rows*sizeof(int));
   }
}

/*******************************************************************************
* Copies int image1 to image2.
*******************************************************************************/
void Copy_int_image(int **image1, int **image2, int cols, int rows)
{
  memcpy((char*)image2[0], (char*)image1[0], cols*rows*sizeof(int));
}

/*******************************************************************************
* Frees an m*n array of ints
*******************************************************************************/
void Free_int_image(int ***ptr)
{
   free((*ptr)[0]);
   free(*ptr);
   *ptr = NULL;
}

/*******************************************************************************
* Returns an m*n array of doubles
*******************************************************************************/
double **Make_double_image(int x, int y)
{
   double **image;
   int i;

   if((image = (double **) calloc(x, sizeof(double *))) == NULL){
      fprintf(stderr, "Error allocating an array in Make_double_image().\n");
      exit(1);
   }
   if((image[0] = (double *) calloc(x*y, sizeof(double))) == NULL){
      fprintf(stderr, "Error allocating an array in Make_double_image().\n");
      exit(1);
   }
   for(i=0;i<x;i++) image[i] = image[0] + (y * i);

   return(image);
}

/*******************************************************************************
* Sets a doubleing point image to val.
*******************************************************************************/
void Set_double_image(double **image, double val, int cols, int rows)
{
   int x, y;
   double *ptr = image[0];

   /* copy first col */
   for(y=0;y<rows;y++) ptr[y] = val;

   for(x=1;x<cols;x++){
      ptr = image[x];
      memcpy((char*)ptr, (char*)image[x-1], rows*sizeof(double));
   }
}

/*******************************************************************************
* Copies double image1 to image2.
*******************************************************************************/
void Copy_double_image(double **image1, double **image2, int cols, int rows)
{
  memcpy((char*)image2[0], (char*)image1[0], cols*rows*sizeof(double));
}

/*******************************************************************************
* Frees an m*n array of doubles
*******************************************************************************/
void Free_double_image(double ***ptr)
{
   free((*ptr)[0]);
   free(*ptr);
   *ptr = NULL;
}
