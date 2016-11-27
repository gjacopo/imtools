#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <configs.h>
#include <pyramidfilters.h>
#include <pyramidtools.h>
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   double *input, *image;   // The input image
   double *output;	         // The output expanded image
   char *Filter; 				
	long Order, ng, nh;				
	double g[MAXF], h[MAXF];
   short IsCentered;
   int rows, cols;               // The dimensions of the image.
	int n;  
	
	// inputs
    /*  create a pointer to the input matrix y */
    input = mxGetPr(prhs[0]);;
    /*  get the dimensions of the matrix input y */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);
	/* Filter identification */
	n = mxGetN(prhs[1])+1;
	Filter = mxCalloc(n,sizeof(char));
	mxGetString(prhs[1],Filter,n);
	/* Order*/
	Order = (long) mxGetScalar(prhs[2]);
	
	/* Initialize the filter */
   	ng = -1L;
	nh = -1L;
	IsCentered = FALSE;

	if ( !strcmp(Filter, "spl"))	{	// "Spline"
		PyramidFilterSplinel2(g, &ng, h, &nh, Order); 
		IsCentered = FALSE;	
	} else if ( !strcmp(Filter, "spl2")) { // "Spline L2"
		PyramidFilterSplineL2(g, &ng, h, &nh, Order);
		IsCentered = FALSE;	
	} else if ( !strcmp(Filter, "cspl")) {	// "Centered Spline"
		PyramidFilterCentered(g, &ng, h, &nh, Order); 
		IsCentered = TRUE;	
	} else if ( !strcmp(Filter, "cspl2"))	{  // "Centered Spline L2"
		PyramidFilterCenteredL2(g, &ng, h, &nh, Order); 
		IsCentered = TRUE;	
	} else if ( ng == -1L && nh == -1L) {
        MessageDisplay( "This familly filters is unknown");
        exit(-1);
    }

	// outputs
    /*  set the output pointer to the output matrix */
	plhs[0] = mxCreateDoubleMatrix(rows/2, cols/2, mxREAL);
	output = mxGetPr(plhs[0]);

	Reduce_2D(input, rows, cols, output, h, nh, IsCentered);
}

